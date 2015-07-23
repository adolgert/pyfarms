import os
import os.path


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
    if filename is None or filename=="":
        raise RuntimeError("No filename given for {0}".format(argument))
    if not os.path.exists(filename):
        basepath=filename
        prevpath=basepath
        while basepath is not "" and not os.path.exists(basepath):
            prevpath=basepath
            basepath=os.path.dirname(basepath)
        if prevpath!=filename:
            raise RuntimeError("The path to {0} doesn't exist so {1} "+
                "cannot be read.".format(basepath, filename))
        else:
            raise RuntimeError("The file {0} doesn't exist in that "+
                "directory".format(filename))
    if not os.path.isfile(filename):
        raise RuntimeError("The path {0} isn't a file.".format(filename))
    if not effectively_readable(filename):
        raise RuntimeError("The file {0} exists, "+
            "but this process cannot read it").format(filename)
