import sys
import os.path
import xml.etree.ElementTree as etree
import xml.parsers.expat.errors
import logging
import numpy as np
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
    scenario.spread_models.clear()

    net=farms.Build(scenario, landscape)

    initial=InitialConditionsNAADSM()
    initial.from_naadsm_file(hxml, ns)
    initial.apply(scenario)

    monitors=Monitors()
    monitors.from_naadsm_file(sxml, ns)

    rng=np.random.RandomState()
    rng.seed(33333)
    sampler=gspn.NextReaction(net, rng)
    run=gspn.RunnerFSM(sampler, observer)
    run.init()
    run.run()


def load_naadsm():
    parser=DefaultArgumentParser(description="Load NAADSM XML Files")
    parser.add_function("load", "Load XML files")
    parser.add_argument("--scenario", dest="scenario_file", type=str,
        action="store", default=None, help="Scenario file")
    parser.add_argument("--herd", dest="herd_file", type=str,
        action="store", default=None, help="Herd file")

    args=parser.parse_args()
    if args.load:
        util.check_filename(args.scenario_file, "scenario file")
        util.check_filename(args.herd_file, "herd file")
        load_naadsm_scenario(args.scenario_file, args.herd_file)

if __name__ == "__main__":
    load_naadsm()

