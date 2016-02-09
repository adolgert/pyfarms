import logging
import os.path
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
    grp_name="run{0:0>4d}".format(dset_idx)
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
