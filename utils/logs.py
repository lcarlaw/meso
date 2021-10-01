import os
import logging
from configs import LOGDIR

def logfile(logname):
    """
    Initiate a logging instance and pass back for writing to the file system. 

    """
    if not os.path.exists(LOGDIR): os.makedirs(LOGDIR)
    logging.basicConfig(filename="%s/%s.log" % (LOGDIR, logname),
                        format='%(levelname)s %(asctime)s :: %(message)s',
                        datefmt="%Y-%m-%d %H:%M:%S")
    log = logging.getLogger()
    log.setLevel(logging.INFO)
    return log
