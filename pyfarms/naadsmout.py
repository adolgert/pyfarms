"""
NAADSM output looks something like this:
node 0 run 5620
1 0 0 0
node 0 run 5620
3 0 0 0
node 0 run 5620
3 0 1 0
node 0 run 5620
3 0 3 0

This module helps read it and can write HDF5.
"""
import logging
import collections
import os.path
import pickle
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import h5py
from default_parser import DefaultArgumentParser

logger=logging.getLogger(__file__)


def combine_counts(starting, combine):
    """
    A dictionary tracks how many of each transition happened
    in each run. This combines those counts.
    """
    for k, v in combine.items():
        if k not in starting:
            starting[k]=v
        else:
            starting[k]+=v
    return starting


def read_multiple_naadsmsc(filename, outfile, initial):
    hdf=h5py.File(outfile, "w")
    hdf.create_group("/trajectory")
    allowed=set()
    for events in naadsm_flow(filename, initial):
        save_h5(hdf, events)
        for e in events:
            allowed.add(e[0])
    return allowed



def naadsm_flow(filename, initial):
    """
    Filename is the name of the text file listing states at each time.
    Outfile is the hdf5 file.
    Initial is the initial set of states at time 0.
    """
    vals=list()
    run=-1
    transitions_type={(0,1) : 0, (1,3) : 1, (3,4) : 3,
        (4,0) : 4, (0,3) : 5, (0,4) : 6, (1,4) : 8}
    with open(filename, "r") as f:
        line=f.readline()
        while line is not "":
            if line.startswith("node"):
                x, n, x, r=line.strip().split()
                n=int(n)
                r=int(r)

                if r!=run:
                    if len(vals)>0:
                        vals=np.vstack(vals)
                        events=events_from_states(vals,
                                transitions_type, initial)
                        yield events
                    vals=list()
                    run=r

                line=f.readline()
                states=[int(y) for y in line.strip().split()]
                vals.append(np.array(states))
            line=f.readline()
    if len(vals)>0:
        vals=np.vstack(vals)
        events=events_from_states(np.vstack(vals), transitions_type, initial)
        yield events



def show_transitions(state_array, initial):
    previous=initial
    allowed=dict()
    for i in range(1, len(state_array)):
        for j in np.where(state_array[i]!=previous)[0]:
            key=(int(previous[j]), int(state_array[i][j]))
            if key not in allowed:
                allowed[key]=0
            allowed[key]+=1
        previous=state_array[i]
    return allowed


def events_from_states(state_array, transitions_dict, initial):
    """
    The initial argument is a list of states at time 0. The first
    output is at day 1.
    """
    events=list()
    previous=initial
    for i in range(0, len(state_array)):
        for j in np.where(state_array[i]!=previous)[0]:
            key=(int(previous[j]), int(state_array[i][j]))
            event=transitions_dict[key]    
            day=i+1
            who=j
            whom=j
            events.append((event, whom, who, day))
        previous=state_array[i]
    logger.debug(events)
    return events

def next_dset(openh5):
    maxnum=-1
    for name in openh5["/trajectory"]:
        if name.startswith("dset"):
            maxnum=max(maxnum, int(name[4:]))
    return maxnum+1


def save_h5(openh5, events):
    dset_idx=next_dset(openh5)
    name="/trajectory/dset{0}".format(dset_idx)
    logger.debug("Writing dataset {0}".format(name))
    group=openh5.create_group(name)
    event=group.create_dataset("Event", (len(events),), dtype="i")
    whom=group.create_dataset("Who", (len(events),), dtype="i")
    who=group.create_dataset("Whom", (len(events),), dtype="i")
    when=group.create_dataset("When", (len(events),), dtype=np.float64)
    for eidx in range(len(events)):
        aevent, awhom, awho, aday=events[eidx]
        event[eidx]=aevent
        whom[eidx]=awhom
        who[eidx]=awho
        when[eidx]=aday


def plot_histogram(x, y, ab, name):
    fig=plt.figure(1, figsize=(8,5))
    ax=fig.add_subplot(111)
    ax.plot(x, y/np.sum(y), 'o')
    xx=np.linspace(0, x[-1], num=100)
    yy=scipy.stats.gamma.pdf(xx, a=ab[0], scale=ab[1])
    print(xx)
    print(yy)
    ax.plot(xx, yy, lw=1, color="black")
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Event Count")
    ax.set_title(name)
    plt.savefig("{0}.pdf".format(name.split()[0]), format="pdf")
    plt.clf()


def single_disease_stats(filename):
    cache_file="single_disease_stats.pickle"
    if not os.path.isfile(cache_file):
        initial=np.array([1])
        times={1 : collections.defaultdict(int),
               3 : collections.defaultdict(int)}
        for events in naadsm_flow(filename, initial):
            last_event_time=0
            for e in events:
                event_time=e[3]
                times[e[0]][event_time-last_event_time]+=1
                last_event_time=event_time
        pickle.dump(times, open(cache_file, "wb"))
    else:
        logger.info("Loading times from {0}".format(cache_file))
        times=pickle.load(open(cache_file, "rb"))
    names={1 : "Latent Period Transition Times",
            3 : "Clinical Period Transition Times"}
    params={1 : (1.34, 0.18), 3 : (13.36, 1.57)}
    for times_idx in [1, 3]:
        print("XXXXXXXXXXXXXXXXX {0}".format(times_idx))
        x=np.zeros(len(times[times_idx]))
        y=np.zeros(len(times[times_idx]))
        i=0
        for time, count in sorted(times[times_idx].items()):
            print("{0} {1}".format(time, count))
            x[i]=time
            y[i]=count
            i+=1
        plot_histogram(x, y, params[times_idx], names[times_idx])


if __name__ == "__main__":
    parser=DefaultArgumentParser(description="Produces csv of total outbreak size")
    parser.add_argument("--input", dest="infile", action="store",
        default="naadsm.out", help="Input trace from NAADSM/SC")
    parser.add_argument("--output", dest="outfile", action="store",
        default="naadsm.h5", help="HDF5 file with events")
    parser.add_function("multiple", "Copy all events to output file")
    parser.add_function("showstates", "Take a look at state transitions")
    parser.add_function("singlefarm", "One farm over and over.")
    args=parser.parse_args()

    if args.multiple:
        initial=np.array([1])
        allowed_transitions=read_multiple_naadsmsc(args.infile, args.outfile,
            initial)
        logger.info("allowed transitions are {0}.".format(allowed_transitions))
    if args.singlefarm:
        single_disease_stats(args.infile)
