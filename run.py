"""This controls the realtime execution of the mesoanalysis scripts.

This replaces the need to specify file and system PATHS via crontab, and instead bundles
everything within this module's directory.

Currently electing to utilize the HRRR for realtime analysis for a few reasons:
    - It's available consistently around :53-54 for all 24 cycles/day as opposed to :30
      for the RAP on its 12z and 00z cycles
    - Can quickly be "upscaled" to the same 13-km grid spacing as the RAP using WGRIB2
      to ease computations
    - The Google Cloud archive for the HRRR is more extensive as it goes back to 2014.
      The RAP currently only goes back to 2021-02-22
    - The HRRR dataset on the Google Cloud is very near real time as opposed to the RAP,
      which seems to be delayed by ~4 hours. This is (potentially) a benefit if both the
      NOMADS and FTPPRD servers are down, although I'm not sure if that also means the
      upload feed to Google will also suffer...
"""

import schedule
import time
import os
import logging
from glob import glob
from datetime import datetime, timedelta

from utils.cmd import execute
from utils.timing import timeit
from configs import PYTHON


script_path = os.path.dirname(os.path.realpath(__file__))
log_dir = "%s/logs" % (script_path)
if not os.path.exists(log_dir): os.makedirs(log_dir)
logging.basicConfig(filename='%s/logs/master.log' % (script_path),
                    format='%(levelname)s %(asctime)s :: %(message)s',
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger()
log.setLevel(logging.INFO)

@timeit
def download_data():
    """Pass arguments to the get_data.py script to download data in realtime from either
    the NOMADS, FTPPRD, or GOOGLE servers
    """
    log.info("Starting download function")
    loop_is_done = False

    start = time.time()
    delta = time.time() - start

    data_path = "%s/IO/data" % (script_path)
    while not loop_is_done and delta < 1800:
        arg = "%s %s/get_data.py -rt -m HRRR -n %s -p %s" % (PYTHON, script_path, 2,
                                                            data_path)
        execute(arg)

        # Determine download status. If we failed, wait and try again.
        with open("%s/download_status.txt" % (script_path)) as f: file_status = f.read()
        if file_status == 'True':
            loop_is_done = True
            log.info("File found and downloaded")
        else:
            log.info("File not found. Clearing directory and sleeping...")
            time_str = (datetime.utcnow()-timedelta(minutes=51)).strftime('%Y-%m-%d/%H')
            files = glob("%s/%s/*" % (data_path, time_str))
            for f in files:
                arg = "rm %s" % (f)
                log.info(arg)
                execute(arg)
            time.sleep(60)
        delta = time.time() - start
    return loop_is_done

def make_placefiles():
    """Pass arguments to the process.py script to create GR-readable placefiles"""
    arg = "%s %s/process.py -rt -meso" % (PYTHON, script_path)
    execute(arg)

# We can schedule each of the tasks at the same time since the execution of task1
# (download_data) will block the execution of task2 (make_placefiles) until the first task
# has been completed.
task1 = schedule.Scheduler()
task2 = schedule.Scheduler()
task1.every().hour.at(":53").do(download_data)
task2.every().hour.at(":53").do(make_placefiles)
while True:
    task1.run_pending()
    task2.run_pending()
    time.sleep(1)
