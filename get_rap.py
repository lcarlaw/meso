from datetime import datetime, timedelta
import requests
import os, sys
import subprocess
import shutil
from distutils.spawn import find_executable
import argparse
from multiprocessing import Pool, freeze_support
import pandas as pd
import numpy as np

from configs import WGRIB2, WGET
from utils.timing import timeit

script_path = os.path.dirname(os.path.realpath(__file__))
#radar_metadata = pd.read_pickle('%s/IO/databases/wsr88d.pkl' % (script_path))

#WGRIB2 = '/usr/local/bin/wgrib2'
#WGET = '/usr/local/bin/wget'
NOMADS = 'https://nomads.ncep.noaa.gov/pub/data/nccf/com/hrrr/prod/hrrr'
FTPPRD = 'https://ftpprd.ncep.noaa.gov/data/nccf/com/hrrr/prod/hrrr' # Haven't added this yet
GOOGLE = 'https://storage.googleapis.com/high-resolution-rapid-refresh'
vars = ':(HGT|TMP|SPFH|UGRD|VGRD|PRES):'
grid_info = 'lambert:262.5:38.5 -105:360:9000 25:310:9000'
DELTA = 3

def execute(arg):
    process = subprocess.Popen(arg, shell=True, stdout=subprocess.PIPE,
              stderr=subprocess.PIPE)
    process.communicate()

def rchop(s, suffix):
    if suffix and s.endswith(suffix):
        return s[:-len(suffix)]
    return s

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
    #lon = radar_metadata[radar]['lon']
    #lat = radar_metadata[radar]['lat']
    save_name = "%s.reduced" % (full_name)
    if not os.path.exists(save_name):
        arg = "%s %s -new_grid %s %s" % (WGRIB2, full_name, grid_info, save_name)
        #bbox_string = "%s:%s %s:%s" % (lon-DELTA, lon+DELTA, lat-DELTA, lat+DELTA)
        #arg = "%s %s -small_grib %s %s" % (WGRIB2, full_name, bbox_string, save_name)
        print(arg)
        execute(arg)

    # Remove the original file
    execute("rm %s" % (full_name))

def execute_download(full_name, url):
    """Download the files in parallel from NOMADS"""

    arg2 = None
    if 'storage.googleapis.com' in url:
        arg1 = '%s -O %s %s' % (WGET, full_name, url)
        arg2 = "%s %s -s | egrep '%s' | %s -i %s -grib %s.tmp" % (WGRIB2,
                                                                      full_name,
                                                                      vars,
                                                                      WGRIB2,
                                                                      full_name,
                                                                      full_name)
    else:
        arg1 = script_path + '/IO/get_inv.pl ' + url + '.idx | egrep ' + "'" + vars +  \
              "'" + ' | ' + script_path + '/IO/get_grib.pl ' + url + ' ' + full_name
    
    if not os.path.exists(full_name + '.reduced'):
        execute(arg1)
    else:
        print("Data already exists locally.")

    if arg2 is not None:
        execute(arg2)
        execute("rm %s" % (full_name))
        execute("mv %s.tmp %s" % (full_name, full_name))

def make_dir(run_time, data_path):
    download_dir = "%s/%s" % (data_path, run_time)
    if not os.path.exists(download_dir): os.makedirs(download_dir)
    return download_dir

@timeit
def download_data(dts, data_path):
    """Function called by main() to control download of model data.

    Surprisingly, HRRR data (both on native and isobaric coordinates), is available in
    near-real time on the Google Cloud Platform, with a latency of only a few minutes over
    NOMADs. In addition, there is an archive that seems to stretch back to 2014!
    """

    #radar = radar.upper()
    if data_path is None: data_path = '%s/data' % (script_path)

    fhr = 1
    downloads = {}
    urls = {}
    for dt in dts:
        download_dir = make_dir(dt.strftime('%Y-%m-%d/%H'), data_path)
        delta = (datetime.utcnow() - dt).total_seconds() / 3600.
        filename = "hrrr.t%sz.wrfnatf%s.grib2" % (str(dt.hour).zfill(2), str(fhr).zfill(2))
        nomads = "%s.%s/conus/%s" % (NOMADS, dt.strftime('%Y%m%d'), filename)
        google = "%s/hrrr.%s/conus/hrrr.t%sz.wrfnatf%s.grib2" % (GOOGLE,
                                                                 dt.strftime('%Y%m%d'),
                                                                 dt.strftime('%H'),
                                                                 str(fhr).zfill(2))
        full_name = "%s/%s" % (download_dir, filename)
        urls = {
            'nomads': nomads,
            'google': google
        }

        # NOMADs
        url_type = None
        if delta < 30:
            if test_url(nomads + '.idx'):
                url_type = 'nomads'
            else:
                if test_url(google):
                    url_type = 'google'

        # Google
        else:
            if test_url(google): url_type = 'google'
        if url_type is not None:
            downloads[full_name] = urls[url_type]


    my_pool = Pool(np.clip(1, len(downloads), 6))
    files = list(downloads.keys())
    urls = list(downloads.values())
    my_pool.starmap(execute_download, zip(files, urls))
    my_pool.map(execute_regrid, files)
    #my_pool.starmap(execute_regrid, zip(files, [radar]*len(downloads)))
    my_pool.close()
    my_pool.terminate()

if __name__ == '__main__':
    freeze_support()

    ap = argparse.ArgumentParser()
    ap.add_argument('-rt', '--realtime', dest="realtime", action='store_true',
                    help='Realtime mode')
    ap.add_argument('-t', '--time-str', dest='time_str', help='YYYY-MM-DD/HH')
    ap.add_argument('-s', dest='start_time', help='Initial valid time for analysis of    \
                    multiple hours. Form is YYYY-MM-DD/HH. MUST be accompanied by the    \
                    "-e" flag. No -t flag is taken.')
    ap.add_argument('-e', dest='end_time', help='Last valid time for analysis')
    #ap.add_argument('-r', '--radar', dest='radar', help='Radar site for local domain')
    ap.add_argument('-p', '--data_path', dest='data_path', help='Directory to store data')
    args = ap.parse_args()

    #if args.radar is None: args.radar = 'KLOT'
    if args.data_path is None: args.data_path = "%s/IO/data" % (script_path)
    timestr_fmt = '%Y-%m-%d/%H'

    if args.realtime or args.time_str is not None:
        if args.realtime:
            cycle_dt = [datetime.utcnow() - timedelta(minutes=55)]
        else:
            cycle_dt = [datetime.strptime(args.time_str, timestr_fmt)]

    elif args.start_time is not None and args.end_time is not None:
        start_dt = datetime.strptime(args.start_time, timestr_fmt)
        end_dt = datetime.strptime(args.end_time, timestr_fmt)
        if start_dt > end_dt:
            print("Requested start time is after the end time")
            sys.exit(1)

        cycle_dt = []
        while start_dt <= end_dt:
            cycle_dt.append(start_dt)
            start_dt += timedelta(hours=1)

    else:
        print("Missing time flags. Need one of -rt, -t, or -s and -e")
        sys.exit()

    #download_data(list(cycle_dt), data_path=args.data_path, radar=args.radar)
    download_data(list(cycle_dt), data_path=args.data_path)
