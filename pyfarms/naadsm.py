"""farmspread

Usage:
  farmspread [-v] [-q] --scenario=<scenariofile> --herd=<herdfile> [--runs=<runs>] [--stream=<stream>]

Options:
  -h --help    Show this screen.
  -v           Verbose
  -q           Quiet
"""
import sys
import os.path
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


# This is the part that runs the SIR
def observer(transition, when):
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


def mainloop(net, scenario, initial, monitors, runs, stream):

    # rng=np.random.RandomState()
    # rng.seed(33333)
    rng=randomstate.prng.pcg64.RandomState(seed=3333334, inc=stream)
    for run_idx in range(runs):
        logger.debug("mainloop run {0}".format(run_idx))
        net.init()
        initial.apply(scenario)
        sampler=gspn.NextReaction(net, rng)
        run=gspn.RunnerFSM(sampler, observer)
        run.init()
        run.run()


def load_naadsm():
    arguments = docopt.docopt(__doc__, version="farmspread")
    if arguments["-v"]:
        logging.basicConfig(level=logging.DEBUG)
    elif arguments["-q"]:
        logging.basicConfig(level=logging.ERROR)
    else:
        logging.basicConfig(level=logging.INFO)

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

    scfile=util.check_filename(arguments["--scenario"], "scenario file")
    hfile=util.check_filename(arguments["--herd"], "herd file")
    net, scenario, initial, monitors=load_naadsm_scenario(scfile, hfile)
    mainloop(net, scenario, initial, monitors, runs, stream)

if __name__ == "__main__":
    load_naadsm()

