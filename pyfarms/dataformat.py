import logging
import collections
import os.path
import sys
import numpy as np
import h5py
from pyfarms.gracefulinterrupthandler import GracefulInterruptHandler

logger=logging.getLogger("dataformat")


def next_dset(openh5):
    """
    Finds the next available dataset name in an HDF5 file.
    """
    maxnum=-1
    for name in openh5:
        if name.startswith("dset"):
            maxnum=max(maxnum, int(name[4:]))
    return maxnum+1


def save_single_run(traj_dir, dset_idx, run, metadata):
    grp_name="dset{0:0>4d}".format(dset_idx)
    group=traj_dir.create_group(grp_name)
    event_cnt=len(run["what"])
    event=group.create_dataset("Event", (event_cnt,),
            dtype="i", data=run["what"])
    whom=group.create_dataset("Who", (event_cnt,),
            dtype="i", data=run["who"])
    who=group.create_dataset("Whom", (event_cnt,),
            dtype="i", data=run["whom"])
    when=group.create_dataset("When", (event_cnt,),
            dtype=np.float64, data=run["when"])
    for k, v in metadata.items():
        group.attrs[k]=v

def trajectory_iter(h5stream, limit=None):
    traj_dir=h5stream.file["/trajectory"]
    idx=0
    for name in traj_dir:
        if name.startswith("dset") or name.startswith("run"):
            yield traj_dir[name]
            idx+=1
            if limit and idx>limit:
                break


def save_runs(outfilename, runs, metadata):
    """
    Takes the name of an HDF5 output file. Appends a trajectory
    to it for each entry in the list of runs. Adds metadata
    as an attribute to each of the runs.
    This can be interrupted with a Ctrl-C without corrupting
    the hdf5 file. It waits for writing to complete.
    """
    with GracefulInterruptHandler() as handler:
        out_data=h5py.File(outfilename, "a")
        base_group="/trajectory"
        if not base_group in dir(out_data):
            traj_dir=out_data.create_group(base_group)
        else:
            traj_dir=out_data["/trajectory"]
        dset_idx=next_dset(traj_dir)
        for (idx, run) in enumerate(runs):
            run_idx=dset_idx+idx
            metadata.update({ "run_idx" : run_idx})
            save_single_run(traj_dir, run_idx, run, metadata)
        out_data.close()

def clear_datafile(outfilename):
    try:
        f=open(outfilename, "wb")
        f.close()
        os.remove(outfilename)
    except Exception as e:
        raise RuntimeError("Could not write to {0} {1}".format(outfilename, e))


def infection_time(h5filename):
    infect_id=set([0, 5, 6])
    infection_times=collections.defaultdict(list)
    trajectory_cnt=0
    with h5py.File(h5filename, "r") as h5stream:
        for trajectory in trajectory_iter(h5stream):
            events=trajectory["Event"]
            who=trajectory["Who"]
            whom=trajectory["Whom"]
            when=trajectory["When"]
            seen=set()
            for idx in range(events.shape[0]):
                if events[idx] in infect_id:
                    infection_times[whom[idx]].append(when[idx])
                    if whom[idx] in seen:
                        logger.error("Found {0} twice".format(whom[idx]))
                        sys.exit(0)
                    seen.add(whom[idx])
            trajectory_cnt+=1
    return infection_times, trajectory_cnt


def first_of_event(h5filename, event_id):
    detection_times=list()
    none_detected=0
    with h5py.File(h5filename, "r") as h5stream:
        for trajectory in trajectory_iter(h5stream):
            events=trajectory["Event"]
            who=trajectory["Who"]
            whom=trajectory["Whom"]
            when=trajectory["When"]
            first=None
            for idx in range(events.shape[0]):
                if events[idx]==event_id:
                    first=when[idx]
                    break
            if first:
                detection_times.append(first)
            else:
                none_detected+=1
    return detection_times, none_detected


def causal_infection(h5filename, event_id):
    event=collections.defaultdict(int)
    none_detected=0
    with h5py.File(h5filename, "r") as h5stream:
        for trajectory in trajectory_iter(h5stream):
            events=trajectory["Event"]
            who=trajectory["Who"]
            whom=trajectory["Whom"]
            when=trajectory["When"]
            for idx in range(events.shape[0]):
                if events[idx]==event_id:
                    event[(who[idx], whom[idx])]+=1
    return event


def disease_states(h5filename, initial_unit, limit=float("inf")):
    """
    Returns a dictionary where keys are the transitions,
    latclin for latent to clinical or clinrec for clinical to recovered,
    and the value is a list of times at which the transition happened.
    """
    lat_clin=set([1, 5, 6, 8])
    clin_rec=set([3,4])
    event=collections.defaultdict(list)
    none_detected=0
    with h5py.File(h5filename, "r") as h5stream:
        for trajectory in trajectory_iter(h5stream, limit):
            events=trajectory["Event"]
            who=trajectory["Who"]
            whom=trajectory["Whom"]
            when=trajectory["When"]
            units=collections.defaultdict(int)
            units[initial_unit]=0
            for idx in range(events.shape[0]):
                if events[idx] in lat_clin:
                    if events[idx]==1 or events[idx]==8:
                        event["latclin"].append(when[idx]-units[whom[idx]])
                    else:
                        event["latclin"].append(0)
                if events[idx] in clin_rec:
                    event["clinrec"].append(when[idx]-units[whom[idx]])
                if whom[idx]>=0:
                    units[whom[idx]]=when[idx]
    return event
