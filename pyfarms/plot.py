"""discreteplot

Make plots from data for discrete distributions.

Usage:
    plot.py [-v] [-q] infect <filename1> <filename2>
    plot.py [-v] [-q] detect <filename1>
    plot.py [-v] [-q] quarantine <filename1>
    plot.py [-v] [-q] locations <filename1> [--cutoff=CUTOFF]
    plot.py [-v] [-q] disease <filename1> <filename2> [--runcnt=RUNCNT]

runc:
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
import math
import re
import tempfile
import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import h5py
import yaml
from docopt import docopt
import lifelines
from PyPDF2 import PdfFileWriter, PdfFileReader
import pyfarms.util as util
import pyfarms.dataformat as dataformat
import pyfarms.naadsm as naadsm
import pyfarms.farms as farms

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
    pdfname=normalize_name(pdfname)
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
        P=[1]*len(fired) + [0]*(fired_max-len(fired))
        ftoplot=fired+([when_max]*(fired_max-len(fired)))
    else:
        P=[1]*len(fired)
        ftoplot=fired
    kmf.fit(ftoplot, P, label=name)
    if ax:
        ay=kmf.plot(ax=ax)
    else:
        ay=kmf.plot()
    return ay


def compare_unit_survival(infect0, infect1, unit,
        traj_cnt0, traj_cnt1, when_max, md):
    kmf=lifelines.KaplanMeierFitter()
    ax=plot_unit_survival(kmf, None, infect0, traj_cnt0, when_max, "Continuous")
    plot_unit_survival(kmf, ax, infect1, traj_cnt1, when_max, "NAADSM")
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
    infection_times0, traj_cnt0=dataformat.infection_time(filenames[0])
    # The NAADSM output is numbered from 0.
    in1, traj_cnt1=dataformat.infection_time(filenames[1])
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

    # a last event much later that ends the simulation.
    for unit in sorted(units):
        inf0=infection_times0[unit]
        inf1=infection_times1[unit]
        compare_unit_survival(inf0, inf1, unit,
                traj_cnt0, traj_cnt1, when_max, md)

    # Now make a bubble plot.
    just_sizes0=[len(x) for x in infection_times0.values()]
    largest=[max(just_sizes0)]
    smallest=[min(just_sizes0)]
    just_sizes1=[len(x) for x in infection_times1.values()]
    largest.append(max(just_sizes1))
    smallest.append(min(just_sizes1))

    print("Trajectory counts for 0 {0} 1 {1}".format(traj_cnt0, traj_cnt1))
    print("unit, times0 times1 percent0 percent1")
    for sunit in sorted(units):
        print("{0} {1} {2} {3} {4}".format(sunit, len(infection_times0[sunit]),
            len(infection_times1[sunit]), len(infection_times0[sunit])/traj_cnt0,
            len(infection_times1[sunit])/traj_cnt1))

    landscape=locations_from_filename(filenames[0])
    x=list()
    y=list()
    color=list()
    colors=["b", "r"]
    size=list()
    events=list()
    for farm_name in units:
        events.append((farm_name-1,
                len(infection_times0[farm_name])/traj_cnt0, 0))
        events.append((farm_name-1,
                len(infection_times1[farm_name])/traj_cnt1, 1))
    # sort by inverse size, so print large ones first.
    events.sort(key=lambda x: -x[1])
    logger.debug("sorted events {0}".format(events))
    for (idx, big, who) in events:
        loc=landscape.farm_locations[idx]
        x.append(loc[1])
        y.append(loc[0])
        color.append(colors[who])
        size.append(3.14*(10*big)**2)

    fig=plt.figure(1)
    ax=fig.add_subplot(111)
    ax.set_title("Infection Frequency")
    logger.debug("sizes {0}".format(size))
    logger.debug("colors {0}".format(color))
    plt.scatter(np.array(x), np.array(y), s=np.array(size),
            marker='o', c=np.array(color))
    SaveFig("scatter_compare.pdf", md)


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
    ax=kmf.plot()
    ax.set_title(name)
    ax.set_xlabel("Days")
    ax.set_ylabel("Survival")
    SaveFig("{0}_survival.pdf".format(name), md)
    plt.clf()
    plt.close()


def locations_from_filename(filename):
    f=h5py.File(filename, "r")["/trajectory"]
    logger.debug(f.keys())
    run=f[list(f.keys())[0]]
    herd=run.attrs["herd"]
    scenario=run.attrs["scenario"]
    sxml, hxml, ns, initial, monitor=naadsm.load_naadsm_scenario(scenario, herd)
    landscape=farms.Landscape()
    landscape.from_naadsm_file(hxml, ns)
    return landscape


def plot_locations(filename, color_cutoff, md):
    landscape=locations_from_filename(filename)
    for idx, farm in enumerate(landscape.premises):
        if int(farm.name) == 21:
            print("farm 21 is idx {0}".format(idx))
    print(landscape.distances[20])
    total=0.0
    for dx in landscape.distances[20]:
        if dx>1e-6:
            total+=0.05**dx
    print("Total hazard for exponential is {0}".format(total))
    total=0.0
    for dx in landscape.distances[20]:
        if dx>1e-6 and dx<20:
            total+=0.05*(20-dx)/19
    print("Total hazard for linear is {0}".format(total))

    infects=dataformat.causal_infection(filename, 0)

    fig=plt.figure(1)
    ax=fig.add_subplot(111)
    ax.set_title("Cause of Infection")
    locations=landscape.farm_locations
    names=[f.name for f in landscape.premises]
    prods=[f.production_type for f in landscape.premises]
    prod_colors={"backyard" : "r", "broilers" : "b",
        "layers" : "g", "turkeys" : "c"}
    logger.info("Unit sizes {0}".format([f.size for f in landscape.premises]))
    for idx in range(len(names)):
        x0, y0=locations[idx]
        plt.plot((y0,), (x0,), prod_colors[prods[idx]]+"o")

    idxof={int(n) : i for (i, n) in enumerate(names)}
    totals=list()
    for (i, j) in itertools.combinations(range(1, len(locations)+1), 2):
        if i < j:
            total=0
            if (i, j) in infects:
                total+=infects[(i,j)]
            if (j, i) in infects:
                total+=infects[(j,i)]
            totals.append((i, j, total))
    maximum=max([x[2] for x in totals])
    totals.sort(key=lambda x: x[2])
    for (i, j, t) in totals:
        x0, y0=locations[idxof[i]]
        x1, y1=locations[idxof[j]]
        dx=x0-x1
        dy=y0-y1
        if math.sqrt(dx*dx+dy*dy)<color_cutoff:
            plt.plot([y0, y1], [x0, x1],
                color=plt.cm.gray((1-t/maximum)**3),
                alpha=t/maximum, linewidth=3)
        else:
            if t>0:
                plt.plot([y0, y1], [x0, x1],
                    color=plt.cm.autumn(1-t/maximum), alpha=t/maximum,
                    linewidth=3)
    ax.legend([Circle((0,0), radius=3, color="r"),
               Circle((0,0), radius=3, color="b"),
               Circle((0,0), radius=3, color="g"),
               Circle((0,0), radius=3, color="c")],
               ["backyard", "broilers", "layers", "turkeys"],
               loc="upper left")
    SaveFig("connectivity.pdf", md)
    plt.clf()

    fig=plt.figure(1)
    ax=fig.add_subplot(111)
    ax.set_title("Farm Locations by Name")
    idxof={int(n) : i for (i, n) in enumerate(names)}
    for idx in range(len(names)):
        x0, y0=locations[idx]
        plt.text(y0, x0, str(names[idx]))
        plt.plot((y0,), (x0,), prod_colors[prods[idx]]+"o")
    SaveFig("locations.pdf", md)


def disease_comparison(times0, times1, name, md):
    logger.debug("times0 len {0} times1 len {1}".format(
            len(times0), len(times1)))
    plt.clf()
    fig=plt.figure(1, figsize=(4,3))
    ax=fig.add_subplot(111)
    kmf=lifelines.KaplanMeierFitter()
    logger.info("Truncating times at 50.")
    for tidx in range(len(times0)):
        if times0[tidx]>50:
            times0[tidx]=50
    P0=[1]*len(times0)
    kmf.fit(times0, P0, label="Continuous")
    ax=kmf.plot(ax=ax)
    ax.set_title(name)
    P1=[1]*len(times1)
    kmf.fit(times1, P1, label="NAADSM")
    kmf.plot(ax=ax)
    plt.tight_layout()
    SaveFig("disease_comparison{0}.pdf".format(name), md)


def plot_disease(filenames, runcnt, md):
    ds0=dataformat.disease_states(filenames[0], 21, runcnt)
    ds1=dataformat.disease_states(filenames[1], 21, runcnt)

    disease_comparison(ds0["latclin"], ds1["latclin"], "Latent to Clinical", md)
    disease_comparison(ds0["clinrec"], ds1["clinrec"], "Clinical to Recovered", md)


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
        plot_detect(arguments["<filename1>"], "First Detection", 10, md)

    if arguments["quarantine"]:
        plot_detect(arguments["<filename1>"], "First Quarantine", 11, md)

    if arguments["locations"]:
        if arguments["--cutoff"]:
            color_cutoff=float(arguments["--cutoff"])
        else:
            color_cutoff=float("inf")
        plot_locations(arguments["<filename1>"], color_cutoff, md)

    if arguments["disease"]:
        if arguments["--runcnt"]:
            runcnt=int(arguments["--runcnt"])
            logger.info("Using runcnt={0}".format(runcnt))
        else:
            runcnt=float("inf")
        plot_disease([arguments["<filename1>"], arguments["<filename2>"]],
                runcnt, md)
