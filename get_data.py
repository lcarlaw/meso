from datetime import datetime, timedelta
import requests
from urllib.request import urlretrieve
import os, sys
import subprocess
import shutil
from distutils.spawn import find_executable
import argparse
from multiprocessing import Pool, freeze_support
import pandas as pd
import numpy as np

from configs import WGRIB2, WGET
from configs import sources, google_configs, thredds_configs, vars, grid_info
from utils.timing import timeit
from utils.cmd import execute

import logging

# Get the script's path and set up our log file
script_path = os.path.dirname(os.path.realpath(__file__))
log_dir = "%s/logs" % (script_path)
if not os.path.exists(log_dir): os.makedirs(log_dir)
logging.basicConfig(filename='%s/logs/download.log' % (script_path),
                    format='%(levelname)s %(asctime)s :: %(message)s',
                    datefmt="%Y-%m-%d %H:%M:%S")
log = logging.getLogger()
log.setLevel(logging.INFO)

def test_url(url):
    """Test for online file existence.

    Parameters
    ----------
    url : string
        URL we're testing for
    """
    ru = requests.head(url)
    return ru.ok

def execute_regrid(full_name):
    """
    Performs subsetting of the original data files to help speed up plotting routines.
    """

    save_name = "%s.reduced" % (full_name)
    if not os.path.exists(save_name):
        if 'ncei' in full_name:
            # WGRIB2 interpolation doesn't seem to work on some of the UGRD/VGRD fields.
            arg = "%s %s -new_grid_order - junk | %s - -new_grid %s %s" % (WGRIB2,
                                                                           full_name,
                                                                           WGRIB2,
                                                                           grid_info,
                                                                           save_name)
        else:
            arg = "%s %s -new_grid %s %s" % (WGRIB2, full_name, grid_info, save_name)

        log.info("CMD %s" % (arg))
        execute(arg)

    # Remove the original file
    execute("rm %s" % (full_name))

def execute_download(full_name, url):
    """Download the files in parallel"""

    arg2 = None
    if 'storage.googleapis.com' in url:
        arg1 = '%s -O %s %s' % (WGET, full_name, url)
        arg2 = "%s %s -s | egrep '%s' | %s -i %s -grib %s.tmp" % (WGRIB2, full_name, vars,
                                                                  WGRIB2, full_name,
                                                                  full_name)
    elif 'ncei' in url:
        arg1 = '%s -O %s %s' % (WGET, full_name, url)
    else:
        arg1 = script_path + '/IO/get_inv.pl ' + url + '.idx | egrep ' + "'" + vars +  \
              "'" + ' | ' + script_path + '/IO/get_grib.pl ' + url + ' ' + full_name

    if not os.path.exists(full_name + '.reduced'):
        execute(arg1)
    else:
        log.info("Data already exists locally")
        pass

    if arg2 is not None:
        execute(arg2)
        execute("rm %s" % (full_name))
        execute("mv %s.tmp %s" % (full_name, full_name))

def make_dir(run_time, data_path):
    download_dir = "%s/%s" % (data_path, run_time)
    if not os.path.exists(download_dir): os.makedirs(download_dir)
    return download_dir

@timeit
def download_data(dts, data_path, model, realtime=True):
    """Function called by main() to control download of model data.

    Returns
    -------
    status : boolean
        Availability of the dataset. False if unavailable.

    """
    if model is None: model = 'RAP'
    if data_path is None: data_path = '%s/data' % (script_path)

    fhr = 1
    downloads = {}
    urls = {}
    for dt in dts:
        download_dir = make_dir(dt.strftime('%Y-%m-%d/%H'), data_path)
        if model == 'HRRR':
            CONUS = 'conus'
            filename = "hrrr.t%sz.wrfnatf%s.grib2" % (str(dt.hour).zfill(2),
                                                      str(fhr).zfill(2))
        elif model == 'RAP':
            CONUS = ''
            filename = "rap.t%sz.awp130bgrbf%s.grib2" % (str(dt.hour).zfill(2),
                                                         str(fhr).zfill(2))
        else:
            log.error("Bad model type passed")
            sys.exit(1)

        full_name = "%s/%s" % (download_dir, filename)
        status = False

        for source in sources.keys():
            base_url = sources[source]

            # NOMADS or backup FTPPRD site. Priority 1 and 2
            if source in ['NOMADS', 'FTPPRD']:
                if model == 'RAP': CONUS = ''
                url = "%s/%s/prod/%s.%s/%s/%s" % (base_url, model.lower(), model.lower(),
                                                  dt.strftime('%Y%m%d'), CONUS, filename)
            # GOOGLE. Priority 3
            elif source == 'GOOGLE':
                # RAP data is ~4 hours behind realtime, while HRRR data is near-realtime.
                # Additionally, the RAP archive is limited, so we're probably ok forcing
                # use of the HRRR here and upscaling to 13 km, in the event that NOMADS &
                # the ftpprd site are down.
                #model = 'HRRR'
                model_name = google_configs[model]
                filename = "hrrr.t%sz.wrfnatf%s.grib2" % (str(dt.hour).zfill(2),
                                                          str(fhr).zfill(2))
                url = "%s/%s/%s.%s/conus/%s" % (base_url, model_name, model.lower(),
                                             dt.strftime('%Y%m%d'), filename)
                full_name = "%s/%s" % (download_dir, filename)

            # THREDDS. Priority 4--just for reanalysis/archive runs. Only RAP available.

            elif source == 'THREDDS':
                # We have two cases here: the RAP and the old RUC. The RAP took over for
                # the 2021-05-01/12z cycle.
                for case in thredds_configs.keys():
                    name = 'rap'
                    if case == 'RUC': name = 'ruc2anl'

                    base_name = thredds_configs[case]
                    filename = "%s-ncei.t%sz.awp130pgrbf%s.grib2" % (case.lower(),
                                                                    str(dt.hour).zfill(2),
                                                                    str(fhr).zfill(2))
                    url = "%s/%s/%s/%s/%s_130_%s_%s00_%s.grb2" % (base_url,
                                                                  thredds_configs[case],
                                                                  dt.strftime('%Y%m'),
                                                                  dt.strftime('%Y%m%d'),
                                                                  name,
                                                                  dt.strftime('%Y%m%d'),
                                                                  str(dt.hour).zfill(2),
                                                                  str(fhr).zfill(3))
                    if test_url(url):
                        full_name = "%s/%s" % (download_dir, filename)
                        break


            status = test_url(url)
            if status:
                log.info("Download source: %s" % (source))
                downloads[full_name] = url
                break

        log.info("Target file: %s" % (full_name))

    if status:
        my_pool = Pool(np.clip(1, len(downloads), 6))
        my_pool.starmap(execute_download, zip(downloads.keys(), downloads.values()))
        my_pool.map(execute_regrid, downloads.keys())
        my_pool.close()
        my_pool.terminate()
    else:
        log.error("*** ERROR finding files ***")

    return status

if __name__ == '__main__':
    freeze_support()    # Needed for multiprocessing.Pool

    ap = argparse.ArgumentParser()
    ap.add_argument('-rt', '--realtime', dest="realtime", action='store_true',
                    help='Realtime mode')
    ap.add_argument('-m', '--model', dest='model', help='RAP or HRRR. Leaving this flag  \
                    off will default to downloading RAP data.')
    ap.add_argument('-t', '--time-str', dest='time_str', help='YYYY-MM-DD/HH')
    ap.add_argument('-s', dest='start_time', help='Initial valid time for analysis of    \
                    multiple hours. Form is YYYY-MM-DD/HH. MUST be accompanied by the    \
                    "-e" flag. No -t flag is taken.')
    ap.add_argument('-e', dest='end_time', help='Last valid time for analysis')
    ap.add_argument('-p', '--data_path', dest='data_path', help='Directory to store data')
    args = ap.parse_args()

    if args.data_path is None: args.data_path = "%s/IO/data" % (script_path)
    timestr_fmt = '%Y-%m-%d/%H'

    if args.realtime or args.time_str is not None:
        if args.realtime:
            cycle_dt = [datetime.utcnow() - timedelta(minutes=51)]
        else:
            cycle_dt = [datetime.strptime(args.time_str, timestr_fmt)]

    elif args.start_time is not None and args.end_time is not None:
        start_dt = datetime.strptime(args.start_time, timestr_fmt)
        end_dt = datetime.strptime(args.end_time, timestr_fmt)
        if start_dt > end_dt:
            log.error("Requested start time is after the end time")
            sys.exit(1)

        cycle_dt = []
        while start_dt <= end_dt:
            cycle_dt.append(start_dt)
            start_dt += timedelta(hours=1)

    else:
        log.error("Missing time flags. Need one of -rt, -t, or -s and -e")
        sys.exit(1)

    log.info("################### New download processing ###################")
    status = download_data(list(cycle_dt), data_path=args.data_path, model=args.model,
                           realtime=args.realtime)
    with open("%s/download_status.txt" % (script_path), 'w') as f: f.write(str(status))
