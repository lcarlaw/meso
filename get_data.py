from datetime import datetime, timedelta
import requests
import os, sys
from glob import glob
import argparse
from multiprocessing import Pool, freeze_support
import pandas as pd
import numpy as np

import timeout_decorator

from configs import WGRIB2, WGET, TIMEOUT
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

def interpolate_in_time(download_dir):
    """Interpolate 1 and 2 hour forecasts in time every 15 minutes"""
    #specs = {
    #    '75 min fcst': [0.75, 0.25],
    #    '90 min fcst': [0.5, 0.5],
    #    '105 min fcst': [0.25, 0.75]
    #}
    arg = "%s %s %s %s %s %s %s %s/%s" % (
         WGRIB2, files[0], '-rpn sto_1 -import_grib', files[1],
         '-rpn sto_2 -set_grib_type same -if',
         '":1 hour fcst:" -rpn "rcl_1:0.5:*:rcl_2:0.5:*:+"',
         '-set_ftime "90 min fcst" -set_scaling same same -grib_out',
         download_dir, '0_50.grib2'
         )
    execute(arg)

    arg = "cp %s %s/%s" % (files[0], download_dir, '0_0.grib2')
    execute(arg)
    arg = "cp %s %s/%s" % (files[1], download_dir, '1_0.grib2')
    execute(arg)

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
            # Need to re-order the ugrd and vgrd entries
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
        arg1 = script_path + '/IO/get_inv.pl ' + url + '.idx | egrep ' + "'" + vars +    \
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

# The goal is to catch hung download processes here.
@timeout_decorator.timeout(TIMEOUT, timeout_exception=StopIteration)
def download_data(dts, data_path, model, num_hours=None, realtime=True):
    """Function called by main() to control download of model data.

    Parameters
    ----------
    dts : list
        List of datetime objects, each corresponding to a desired model cycle
    data_path: string
        Path to download data to. Defaults to IO/data
    model : string
        Desired model. Options are HRRR and RAP
    num_hours : int
        Last forecast hour. Default: 1
    realtime : boolean
        If this is a realtime or archived run. Default: True

    Returns
    -------
    status : boolean
        Availability of the dataset. False if unavailable.

    """

    return_status = False
    if num_hours is None: num_hours = 1
    if model is None: model = 'RAP'
    if data_path is None: data_path = '%s/data' % (script_path)

    fhrs = np.arange(1, int(num_hours)+1, 1)
    downloads = {}
    urls = {}
    expected_files = 0
    for dt in dts:
        for fhr in fhrs:
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
            for source in sources.keys():
                base_url = sources[source]

                # NOMADS or backup FTPPRD site. Priority 1 and 2
                if source in ['NOMADS', 'FTPPRD']:
                    if model == 'RAP': CONUS = ''
                    url = "%s/%s/prod/%s.%s/%s/%s" % (base_url,model.lower(),
                                                      model.lower(),dt.strftime('%Y%m%d'),
                                                      CONUS, filename)
                # GOOGLE. Priority 3
                elif source == 'GOOGLE' and model not in ['RAP']:
                    # RAP data is ~4 hours behind, while HRRR data is near-realtime.
                    # Additionally, the RAP archive is limited, so we're ok forcing
                    # use of the HRRR here and upscaling to 13 km, in the event NOMADS
                    # and the ftpprd site are down.
                    model_name = google_configs[model]
                    filename = "hrrr.t%sz.wrfnatf%s.grib2" % (str(dt.hour).zfill(2),
                                                              str(fhr).zfill(2))
                    url = "%s/%s/%s.%s/conus/%s" % (base_url, model_name, model.lower(),
                                                 dt.strftime('%Y%m%d'), filename)
                    full_name = "%s/%s" % (download_dir, filename)

                # THREDDS. Priority 4--just for reanalysis runs. Only RAP available.
                elif source == 'THREDDS':
                    # We have two cases: the RAP and the old RUC. The RAP took over for
                    # the 2021-05-01/12z cycle.
                    for case in thredds_configs.keys():
                        name = 'rap'
                        if case == 'RUC': name = 'ruc2anl'

                        base_name = thredds_configs[case]
                        filename = "%s-ncei.t%sz.%sf%s.%s" % (case.lower(),
                                                              str(dt.hour).zfill(2),
                                                              'awp130pgrb',
                                                              str(fhr).zfill(2),
                                                              'grib2')

                        url = "%s/%s/%s/%s/%s_%s_%s_%s00_%s.%s" % (base_url,
                                                                   thredds_configs[case],
                                                                   dt.strftime('%Y%m'),
                                                                   dt.strftime('%Y%m%d'),
                                                                   name, '130',
                                                                   dt.strftime('%Y%m%d'),
                                                                   str(dt.hour).zfill(2),
                                                                   str(fhr).zfill(3),
                                                                   'grb2')

                        if test_url(url):
                            full_name = "%s/%s" % (download_dir, filename)
                            break

                status = test_url(url)
                if status:
                    log.info("Download source: %s" % (source))
                    downloads[full_name] = url
                    break

            log.info("Target file: %s" % (full_name))
            expected_files += 1

    if len(downloads.keys()) >= 1:
        my_pool = Pool(np.clip(1, len(downloads), 4))
        my_pool.starmap(execute_download, zip(downloads.keys(), downloads.values()))
        my_pool.map(execute_regrid, downloads.keys())
        my_pool.close()
        my_pool.terminate()
    else:
        log.error("*** ERROR finding files ***")

    # Make sure we've got all the expected files.
    knt = 0
    for f in downloads.keys():
        filename = f + '.reduced'
        if os.path.exists(filename):
            if os.stat(filename).st_size > 1024000.:
                knt += 1
            else:
                log.error("%s is %s MB" % (filename, os.stat(filename).st_size/1024000.))
        else:
            log.error("%s doesn't exist" % (filename))

    if knt == expected_files:
        return_status = True
        log.info("Success downloading files")

    return return_status, download_dir

if __name__ == '__main__':
    freeze_support()    # Needed for multiprocessing.Pool

    ap = argparse.ArgumentParser()
    ap.add_argument('-rt', '--realtime', dest="realtime", action='store_true',
                    help='Realtime mode')
    ap.add_argument('-m', '--model', dest='model', help='RAP or HRRR. Leaving this flag  \
                    off will default to downloading RAP data.')
    ap.add_argument('-n', '--num-hours', dest='num_hours', help='Number of forecast      \
                    hours to download')
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

    # Warning if user has selected archived (non-native coordinate) RAP data
    if args.model == 'RAP' and not args.realtime:
        warn = """*WARNING*: Archived RAP data is only available on isobaric coordinates.
           Due to the sensitivity of mixed-layer and most-unstable parcel
           calculations to near-surface vertical resolution, this coarse
           dataset may result in erroneous thermodynamic calculations. Are
           you sure you want to proceed? [y|N]"""
        print(warn)
        resp = input()
        if resp not in ['y', 'Y', 'yes'] or resp in ['n', 'N']: sys.exit(1)

    log.info("################### New download processing ###################")
    with open("%s/download_status.txt" % (script_path), 'w') as f: f.write(str(False))
    status, download_dir = download_data(list(cycle_dt), data_path=args.data_path,
                                         model=args.model, num_hours=args.num_hours,
                                         realtime=args.realtime)
    with open("%s/download_status.txt" % (script_path), 'w') as f: f.write(str(status))

    # If this is realtime, interpolate the 1 and 2-hour forecasts in time
    if status and args.realtime: interpolate_in_time(download_dir)
