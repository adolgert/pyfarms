import os
import os.path
import re
import glob
import logging
import numpy as np

logger=logging.getLogger("pyframs.util")

def effectively_readable(path):
    import os, stat

    uid = os.getuid()
    euid = os.geteuid()
    gid = os.getgid()
    egid = os.getegid()

    # This is probably true most of the time, so just let os.access()
    # handle it.  Avoids potential bugs in the rest of this function.
    if uid == euid and gid == egid:
        return os.access(path, os.R_OK)

    st = os.stat(path)

    # This may be wrong depending on the semantics of your OS.
    # i.e. if the file is -------r--, does the owner have access or not?
    if st.st_uid == euid:
        return st.st_mode & stat.S_IRUSR != 0

    # See comment for UID check above.
    groups = os.getgroups()
    if st.st_gid == egid or st.st_gid in groups:
        return st.st_mode & stat.S_IRGRP != 0

    return st.st_mode & stat.S_IROTH != 0



def check_filename(filename, argument):
    """
    Given a command line argument that should be a readable file,
    ensure it is a readable file and give informative messages if not.
    If it's not a file but a directory, tell me. If the path is wrong,
    tell me which part of th epath is wrong. If the file exists but
    I don't have permission to read it, tell me that. Gosh!
    """
    if filename is None or filename=="":
        raise RuntimeError("No filename given for {0}".format(argument))
    # Expand the ~ into the home directory in case shell expansion failed.
    filename=re.sub("^\~", os.environ["HOME"], filename)
    candidates=glob.glob(filename)
    if len(candidates) is 1:
        filename=candidates[0]
    elif len(candidates) > 1:
        raise RuntimeError("More than one file matches {0}".format(candidates))
    else:  #len(candidates) is 0:
        basepath=filename
        prevpath=basepath
        while basepath is not "" and not os.path.exists(basepath):
            prevpath=basepath
            basepath=os.path.dirname(basepath)
        if prevpath!=filename:
            raise RuntimeError(("The path to {0} doesn't exist so {1} "+
                "cannot be read.").format(basepath, filename))
        else:
            raise RuntimeError(("The file {0} doesn't exist in that "+
                "directory").format(filename))
    if not os.path.isfile(filename):
        raise RuntimeError("The path {0} isn't a file.".format(filename))
    if not effectively_readable(filename):
        raise RuntimeError(("The file {0} exists, "+
            "but this process cannot read it").format(filename))
    return filename



_degrees_to_radians=np.pi/180
_radians_km=180*60*1.852/np.pi

def distancekm(latlon1, latlon2):
    """Distance computed on a spherical earth.
    Taken from http://williams.best.vwh.net/avform.htm."""
    ll1=latlon1*_degrees_to_radians
    ll2=latlon2*_degrees_to_radians
    return _radians_km*(2*np.arcsin(np.sqrt(np.power(np.sin((ll1[0]-ll2[0])/2),2)+
        np.cos(ll1[0])*np.cos(ll2[0])*np.power(np.sin((ll1[1]-ll2[1])/2), 2))))
