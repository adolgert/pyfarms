import os
import os.path
import re
import subprocess
import datetime
import sys
import glob
import math
import logging
import numpy as np
import pyproj

logger=logging.getLogger("pyfarms.util")


def gitinfo():
    try:
        git_ret=subprocess.Popen(['git','log','--pretty=%H','HEAD^..HEAD'],
            stdout=subprocess.PIPE)
        git_hash = git_ret.communicate()[0]
    except Exception as e:
        print(e)
        return None
    if git_hash:
        git_hash=git_hash.strip().decode()
        try:
            url_ret=subprocess.Popen(['git','remote','show','origin'],
                stdout=subprocess.PIPE)
            remote=url_ret.communicate()[0].decode()
        except Exception as e:
            return None
        match=re.search('URL:\s*(\S+)\n',remote)
        if match:
            git_url=match.group(1)
            scmversion='{0}:{1}'.format(git_url, git_hash)
        else:
            scmversion=git_hash
        return scmversion
    else:
        return None

def makefile():
    return open("Makefile").read()


def when():
    return datetime.datetime.now().isoformat()

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

def GIS_default_projection(latlon):
    """
    Given an array of latitude and longitude, return the same projection
    NAADSM uses.
    """
    minlat=np.min(latlon[:,0])
    maxlat=np.max(latlon[:,0])
    avglat=np.mean(latlon[:,0])
    minlon=np.min(latlon[:,1])
    maxlon=np.max(latlon[:,1])
    avglon=np.mean(latlon[:,1])

    projstr="+ellps=WGS84 +units=km +lon_0={0} +proj=aea +lat_0={1} +lat_1={2} +lat_2={3}".format(
        avglon, minlat, minlat+(maxlat-minlat)/6, maxlat-(maxlat-minlat)/6)
    logger.debug("Projection string {0}".format(projstr))
    pp=pyproj.Proj(projstr)
    projected=np.zeros(latlon.shape, dtype=np.double)
    for idx in range(latlon.shape[0]):
        projected[idx,:]=np.array(pp(latlon[idx,1], latlon[idx,0]))/1000
    logger.debug("latlon {0}".format(latlon))
    logger.debug("Projected GIS {0}".format(projected))
    return projected

def GIS_distance(latlon1, latlon2):
    """
    This function, named GIS_distance in NAADSM 3.2, is exactly the same.
    Not kidding.
    """
    x=latlon1[0]-latlon2[0]
    y=latlon1[1]-latlon2[1]
    return math.sqrt(x*x + y*y)

class ChunkIter(object):
    """
    Someone asks for 2100 iterations but wants them in chunks
    of 250, so this parcels that out, from (0, 250), (1, 250),
    to (8, 100).
    """
    def __init__(self, chunk, total):
        self.total=total
        self.chunk=chunk

    def __iter__(self):
        self.idx=0
        return self

    def __next__(self):
        begin=self.idx*self.chunk
        end=(self.idx+1)*self.chunk
        if begin<self.total:
            end=min(self.total, end)
            self.idx+=1
            return (self.idx, end-begin)
        else:
            raise StopIteration()
