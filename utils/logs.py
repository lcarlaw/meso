import logging
from configs import LOGDIR

def logfile(logname):
    logging.basicConfig(filename="%s/%s.log" % (LOGDIR, logname),
                        format='%(levelname)s %(asctime)s :: %(message)s',
                        datefmt="%Y-%m-%d %H:%M:%S")
    log = logging.getLogger()
    log.setLevel(logging.INFO)
    return log
