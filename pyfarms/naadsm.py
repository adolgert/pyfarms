"""farmspread

Usage:
  farmspread [options] --scenario=<scenariofile> --herd=<herdfile>

Options:
  -h --help       Show this screen.
  -v              Verbose.
  -q              Quiet.
  --runs=CNT      How many trajectories to gather. [default: 1]
  --stream=STREAM The integer number of the random number stream. [default: 1]
  --out=OUT       Name of an output file. [default: out.hdf5]
  --chunk=CHUNK   How many runs to do before saving intermediate results.
"""
import sys
import os
import os.path
import re
import xml.etree.ElementTree as etree
import xml.parsers.expat.errors
import logging
import numpy as np
import docopt
import randomstate
from pyfarms.default_parser import DefaultArgumentParser
import pyfarms.util as util
import pyfarms.farms as farms
import gspn
import pyfarms.dataformat as dataformat

logger=logging.getLogger("pyfarms.naadsm")


class Monitors(object):
    def __init__(self):
        self.infection=False
        self.exposure=False

    def from_naadsm_file(self, root, ns):
        models=root.find("models", ns)
        for monitor in ["infection", "exposure", "detection",
                "trace", "destruction", "destruction-list", "vaccination",
                "vaccination-list", "nonexistent"]:
            search_string="{0}-monitor".format(monitor)
            im=models.find(search_string, ns)
            if im is not None:
                setattr(self, monitor, True)
            else:
                setattr(self, monitor, False)
            logging.debug("Monitors {0} {1}".format(search_string,
                getattr(self, monitor)))

        outputs=root.find("output", ns)
        # Don't know how to read outputs yet. Only entry says "all".


class InitialConditionsNAADSM(object):
    def __init__(self):
        self.infections=dict()

    def from_naadsm_file(self, root, ns):
        """
        Pass in a Herd file.
        """
        for unit in root.findall("herd", ns):
            unit_name=unit.find("id").text
            status=unit.find("status").text
            if status!="Susceptible":
                self.infections[unit_name]=status

    def apply(self, scenario):
        for f in scenario.farms:
            if f.name in self.infections.keys():
                f.disease.initial_infection()
            else:
                f.disease.initial_susceptible()

class TransID(object):
    infected=0
    infectious=1
    recover=4
    wane=5
    detect=10
    quarantine=11
    restriction=12


_transition_id={
    (farms.DiseaseState.susceptible, farms.DiseaseState.latent) : TransID.infected,
    (farms.DiseaseState.latent, farms.DiseaseState.clinical) : TransID.infectious,
    (farms.DiseaseState.clinical, farms.DiseaseState.recovered) : TransID.recover,
    (farms.DiseaseState.recovered, farms.DiseaseState.susceptible) : TransID.wane,
}


class StateObserver(object):
    def __init__(self):
        self.events=list()
        self.runs=list()

    def init(self):
        if self.events:
            trace=dict()
            trace["what"]=[x[0] for x in self.events]
            trace["who"]=[x[1] for x in self.events]
            trace["whom"]=[x[2] for x in self.events]
            trace["when"]=[x[3] for x in self.events]
            self.runs.append(trace)
            self.events=list()

    def results(self):
        if self.events:
            self.init()
        return self.runs

    def __call__(self, transition, when):
        tname=transition.__class__.__name__
        if tname=="DiseaseABTransition":
            tid=_transition_id[(transition.a, transition.b)]
            whom=transition.farm[0].name
            who=whom
        elif tname=="DetectionTransition":
            tid=TransID.detect
            whom=-1
            who=whom
        elif tname=="QuarantineTransition":
            tid=TransID.quarantine
            whom=transition.farm[0].name
            who=whom
        elif tname=="InfectTransition":
            tid=TransID.infected
            who=transition.farm[0].name
            whom=transition.farm[1].name
        elif tname=="IndirectTransition":
            tid=TransID.infected
            who=transition.farm[0].name
            whom=transition.farm_models[transition.affected_idx].name
        elif tname=="RestrictionTransition":
            tid=TransID.restriction
            who=-1
            whom=who
        else:
            logger.error("Unknown transition {0}".format(tname))

        logger.info((tid, who, whom, when))
        self.events.append((tid, int(who), int(whom), when))
        return when<365


# This is the part that runs the SIR
def observe(transition, when):
    tname=transition.__class__.__name__
    if len(transition.farm)>1:
        who, whom=[x.name for x in transition.farm]
    else:
        whom=transition.farm[0].name
        who=whom
    if tname=="DiseaseABTransition":
        print("{0} {1}".format(transition, when))
    elif tname=="IndirectTransition":
        print("{0} {1}".format(transition, when))
    else:
        print("{0} {1} {2} {3}".format(tname, who, whom, when))
    return when<365


def load_naadsm_scenario(scenario_filename, herd_filename):
    ns={"naadsm" : "http://www.naadsm.org/schema",
        "xsd" : "http://www.w3.org/2001/XMLSchema",
        "xml" : "http://www.w3.org/XML/1998/namespace",
        "gml" : "http://www.opengis.net/gml",
        "xsi" : "http://www.w3.org/2001/XMLSchema-instance"}
    try:
        sxml=etree.parse(scenario_filename)
    except etree.ParseError as err:
        logger.error("Could not parse {0} at line {1} with error \"{2}\"".format(
            scenario_filename, err.position,
            xml.parsers.expat.errors.messages[err.code]))
        sys.exit(1)
    try:
        hxml=etree.parse(herd_filename)
    except etree.ParseError as err:
        logger.error("Could not parse {0} at line {1} with error {2}".format(
            herd_filename, err.position,
            xml.parsers.expat.errors.messages[err.code]))
        sys.exit(1)
    
    landscape=farms.Landscape()
    landscape.from_naadsm_file(hxml, ns)

    scenario=farms.Scenario()
    scenario.from_naadsm_file(sxml, ns)

    # Turn off all airborne spread.
    #scenario.spread_models.clear()

    net=farms.Build(scenario, landscape)

    initial=InitialConditionsNAADSM()
    initial.from_naadsm_file(hxml, ns)

    monitors=Monitors()
    monitors.from_naadsm_file(sxml, ns)

    return net, scenario, initial, monitors


def mainloop(net, scenario, initial, observer, runs, stream):
    rng=randomstate.prng.pcg64.RandomState(seed=3333334, inc=stream)
    for run_idx in range(runs):
        logger.debug("mainloop run {0}".format(run_idx))
        net.init()
        observer.init()
        initial.apply(scenario)
        sampler=gspn.NextReaction(net, rng)
        run=gspn.RunnerFSM(sampler, observer)
        run.init()
        run.run()


def multirun(scfile, hfile, runs, stream):
    net, scenario, initial, monitors=load_naadsm_scenario(scfile, hfile)
    observer=StateObserver()
    mainloop(net, scenario, initial, observer, runs, stream)
    return observer.results(), stream



def load_naadsm():
    arguments = docopt.docopt(__doc__, version="farmspread")
    if arguments["-v"]:
        logging.basicConfig(level=logging.DEBUG)
    elif arguments["-q"]:
        logging.basicConfig(level=logging.ERROR)
    else:
        logging.basicConfig(level=logging.INFO)

    metadata=dict()

    if arguments["--runs"]:
        runs=int(arguments["--runs"])
        assert(runs>0)
    else:
        runs=1
    if arguments["--stream"]:
        print(arguments)
        stream=int(arguments["--stream"])
        assert(stream>0)
    else:
        stream=1

    if arguments["--out"]:
        outfile=arguments["--out"]
        outfile=re.sub("^\~", os.environ["HOME"], outfile)
    else:
        outfile="out.hdf5"

    if arguments["--chunk"]:
        chunk_size=int(arguments["--chunk"])
    else:
        chunk_size=runs

    dataformat.clear_datafile(outfile)

    scfile=util.check_filename(arguments["--scenario"], "scenario file")
    hfile=util.check_filename(arguments["--herd"], "herd file")

    metadata.update({"run_cnt" : runs, "stream_idx" : stream,
        "datafile" : outfile, "chunk_size" : chunk_size,
        "herd" : hfile, "scenario" : scfile })

    for chunk_idx, run_cnt in util.ChunkIter(runs, chunk_size):
        results, stream=multirun(scfile, hfile, run_cnt, stream+chunk_size)
        dataformat.save_runs(outfile, results, metadata)

if __name__ == "__main__":
    load_naadsm()

