"""discreteplot

Make plots from data for discrete distributions.

Usage:
    plot.py [-v] [-q] infect <filename1> <filename2>
    plot.py [-v] [-q] detect <filename1>
    plot.py [-v] [-q] quarantine <filename1>
    plot.py [-v] [-q] locations <filename1>

Options:
    -h --help             Print this screen.
    -v                    Verbose
    -q                    Quiet
    -t                    Testing flag
"""
import logging
import sys
import docopt
import shutil
import os
import itertools
import tempfile
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import h5py
import yaml
from docopt import docopt
import lifelines
from PyPDF2 import PdfFileWriter, PdfFileReader
import pyfarms.util as util
import pyfarms.dataformat as dataformat
import pyfarms.naadsm as naadsm

logger=logging.getLogger("pyfarms.plot")


def AddPdfKeys(pdf_file, keydict):
    """
    Add a dictionary of key-value pairs to a pdf file.
    Adds metadata to PDF files so you know what generated them.
    """
    multifiles=list()
    if isinstance(pdf_file, list):
        multifiles=pdf_file
    else:
        multifiles=[pdf_file]
    for of in multifiles:
        outstream=tempfile.NamedTemporaryFile(mode="wb", delete=False)
        logger.debug("AddPdfKeys created {0}".format(outstream.name))
        CopyPdfKeys(of, outstream, keydict)
        outstream.close()
        shutil.copyfile(outstream.name, of)
        os.remove(outstream.name)

def SaveFig(pdfname, md):
    logger.info("Writing {0}".format(pdfname))
    plt.savefig(pdfname)
    AddPdfKeys(pdfname, md)

def CopyPdfKeys(pdf_file, outstream, keydict):
    pdfDict={u"/{0}".format(k) : v for (k, v) in keydict.items()}
    infile=PdfFileReader(pdf_file)
    outfile=PdfFileWriter()
    outfile.addMetadata(pdfDict)
    for page in range(infile.getNumPages()):
        outfile.addPage(infile.getPage(page))
    outfile.write(outstream)


def normalize_name(title):
    return re.sub(r"-\\\"\(\)\{\}", "", re.sub(r"\s", "_", title))


def plot_unit_survival(kmf, ax, fired, fired_max, when_max, name):
    if len(fired)<fired_max:
        fired.extend([when_max]*(fired_max-len(fired)))
        P=[1]*len(fired) + [0]*(fired_max-len(fired))
    else:
        P=[1]*len(fired)
    kmf.fit(fired, P, label=name)
    if ax:
        ay=kmf.plot(ax=ax)
    else:
        ay=kmf.plot()
    return ay


def compare_unit_survival(infect0, infect1, unit, fired_max, when_max, md):
    kmf=lifelines.KaplanMeierFitter()
    ax=plot_unit_survival(kmf, None, infect0, fired_max, when_max, "Continuous")
    plot_unit_survival(kmf, ax, infect1, fired_max, when_max, "NAADSM")
    SaveFig("unit_survival{0}.pdf".format(unit), md)
    plt.clf()
    plt.close()


def plot_infect(filenames, md):
    """
    The goal is to make a comparison plot between two simulations,
    one plot per unit. Each plot shows the distribution of the number
    of days it takes for that unit to become infected. It's
    a survival plot, so it starts at one. If it doesn't go to zero,
    that means that, for some runs, that unit doesn't get infected.
    """
    infection_times0=dataformat.infection_time(filenames[0])
    # The NAADSM output is numbered from 0.
    in1=dataformat.infection_time(filenames[1])
    infection_times1={n+1 : v for (n,v) in in1.items()}
    logger.debug(infection_times0.keys())
    logger.debug(infection_times1.keys())
    unitsa=set([int(x) for x in infection_times0.keys()])
    unitsb=set([int(x) for x in infection_times1.keys()])
    units=unitsa&unitsb
    logger.info("Units in a, not b {0}".format(unitsa-unitsb))
    logger.info("Units in b, not a {0}".format(unitsb-unitsa))

    fired_max=-1
    when_max=-1
    for events in itertools.chain(infection_times0.values(),
            infection_times1.values()):
        fired_max=max(fired_max, len(events))
        when_max=max(when_max, max(events))

    for unit in sorted(units):
        inf0=infection_times0[unit]
        inf1=infection_times1[unit]
        compare_unit_survival(inf0, inf1, unit, fired_max, when_max, md)


def plot_detect(filename, name, event_id, md):
    """
    What is the distribution of times that infection is first detected.
    """
    detection_times, none_detected=dataformat.first_of_event(filename, event_id)
    logger.info("Detected {0} times out of {1}".format(len(detection_times),
            len(detection_times)+none_detected))
    if len(detection_times) is 0:
        logger.info("The event {0} did not happen.".format(event_id))
        sys.exit(0)
    kmf=lifelines.KaplanMeierFitter()
    last=max(detection_times)+1
    detection=np.hstack([np.array(detection_times),
            last*np.ones((none_detected,), dtype=np.double)])
    P=[1]*len(detection_times)+[0]*none_detected
    kmf.fit(detection, P, label=name)
    kmf.plot()
    SaveFig("{0}_survival.pdf".format(name), md)
    plt.clf()
    plt.close()


def plot_locations(filename, md):
    pass



if __name__ == "__main__":
    arguments = docopt(__doc__, version="pyfarms.plot 1.0")
    if arguments["-v"]:
        logging.basicConfig(level=logging.DEBUG)
    elif arguments["-q"]:
        logging.basicConfig(level=logging.ERROR)
    else:
        logging.basicConfig(level=logging.INFO)

    md={
        "script" : __file__, "arguments" : str(arguments),
        "date" : util.when()
    }

    if arguments["infect"]:
        plot_infect([arguments["<filename1>"], arguments["<filename2>"]], md)

    if arguments["detect"]:
        plot_detect(arguments["<filename1>"], "Detection", 10, md)

    if arguments["quarantine"]:
        plot_detect(arguments["<filename1>"], "Quarantine", 11, md)

    if arguments["locations"]:
        plot_locations(arguments["<filename1>"], md)