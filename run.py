import os, sys
import argparse
import numpy as np
from glob import glob
from datetime import datetime, timedelta
from time import time
import re

import sharptab.calcs as calcs
import sharptab.interp as interp
import sharptab.winds as winds
import utils.plot_hodos as plot_hodos
import utils.hodographs as hodographs
import utils.plot as plot

import IO.read as read
from configs import PYTHON, plotinfo

from utils.timing import timeit
from utils.cmd import execute

script_path = os.path.dirname(os.path.realpath(__file__))

def create_hodograph(data, point, storm_motion='right-mover', storm_relative=False):
    idx = interp.nearest_idx(point, data['lons'], data['lats'])
    u, v = winds.vec2comp(data['wdir'][:,idx[0][0], idx[0][1]],
                          data['wspd'][:,idx[0][0], idx[0][1]])
    hodo_data = {}
    heights = data['hght'][:,idx[0][0], idx[0][1]] / 1000.
    mask = np.where(heights < 12)
    hodo_data['hght'] = heights[mask]
    hodo_data['uwnd'] = u[mask]
    hodo_data['vwnd'] = v[mask]
    hodo_data['valid_time'] = data['valid_time']
    hodo_data['cycle_time'] = data['cycle_time']
    hodo_data['fhr'] = data['fhr']
    hodo_data['lon'], hodo_data['lat'] = point[0][0], point[0][1]

    params = hodographs.compute_parameters(hodo_data, storm_motion)
    plot_hodos.plot_hodograph(hodo_data, params, storm_relative=storm_relative)

@timeit
def create_placefiles(data, realtime=False):
    prof_data = {'pres':data['pres'], 'tmpc':data['tmpc'],
                 'dwpc':data['dwpc'], 'hght':data['hght'],
                 'wdir':data['wdir'], 'wspd':data['wspd']}
    arrs = calcs.sharppy_calcs(**prof_data)
    for item in ['valid_time', 'cycle_time', 'fhr', 'lons', 'lats']:
        arrs[item] = data[item]
    plot.write_placefile(arrs, plotinfo, realtime=realtime)

def find_nearest(dt, datadir):
    """Find the nearest available RAP forecast to the current time or user-requested time

    Parameters:
    -----------
    dt : datetime
        The current time or user-requested time

    """
    regex_str = 'f\d{2}.grib2.reduced'

    deltas = {}
    dirs = glob(datadir + '/*')
    min_delta = 99999999
    for dirname in dirs:
        files = glob(dirname + '/**/*.reduced', recursive=True)

        for f in files:
            shortname = f.lstrip(datadir)
            cycle = datetime.strptime(shortname[0:13], '%Y-%m-%d/%H')
            fhr = int(re.findall(regex_str, shortname)[0][1:3])
            valid_time = cycle + timedelta(hours=fhr)

            # Apply a penalty to older model cycles valid at the same time
            increment = (dt.hour - cycle.hour)
            delta = np.fabs((valid_time - dt).total_seconds()) + increment
            deltas[delta] = f
            if np.fabs(delta) < min_delta: min_delta = delta
    return deltas[min_delta]

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    #ap.add_argument('radar', help='Radar site ID to localize to')
    ap.add_argument('-rt', '--realtime', dest="realtime", action='store_true',
                    help='Realtime mode. Script will look for the most recent model run.')
    ap.add_argument('-s', dest='start_time', help='Initial valid time for analysis of    \
                    multiple hours. Form is YYYY-MM-DD/HH. MUST be accompanied by the    \
                    "-e" flag. No -t flag is taken.')
    ap.add_argument('-e', dest='end_time', help='Last valid time for analysis')
    ap.add_argument('-t', '--time-str', dest='time_str', help='Valid time for archived   \
                    model runs. Script will default to 1-hr forecasts. YYYY-MM-DD/HH')
    ap.add_argument('-p', '--data_path', dest='data_path', help='Directory to store model\
                    data. Default will be set to the IO/data directory.')
    ap.add_argument('-meso', dest='meso', action='store_true',
                    help='Flag to output mesoanalysis placefiles for GR2/Analyst.')
    ap.add_argument('-hodo', '--hodograph', dest='hodo',
                    help='Hodograph plot at a point ex: 35.03/-101.23')
    ap.add_argument('-sr', '--storm-relative', dest='storm_relative', action='store_true',
                    help='Flag to create a Storm Relative hodograph')
    ap.add_argument('-m', '--storm-motion', dest='storm_motion',
                    help='Storm motion vector. It takes one of two forms. The first is   \
                          either "BRM" for the Bunkers right mover vector, or "BLM" for  \
                          the Bunkers left mover vector. The second is the form DDD/SS,  \
                          where DDD is the direction the storm is coming from, and SS is \
                          the speed in knots (e.g. 240/25).', default='right-mover')
    args = ap.parse_args()

    timestr_fmt = '%Y-%m-%d/%H'
    #args.radar = args.radar.upper()
    if args.realtime:
        dt = datetime.utcnow()
        time_flag = "-rt"
    else:
        if args.time_str is not None:
            dt = datetime.strptime(args.time_str, timestr_fmt)
            time_flag = "-t %s" % ((dt-timedelta(hours=1)).strftime(timestr_fmt))
        elif args.start_time is not None and args.end_time is not None:
            dt_start = datetime.strptime(args.start_time, timestr_fmt)
            dt_start -= timedelta(hours=1)
            dt_end = datetime.strptime(args.end_time, timestr_fmt)
            dt_end -= timedelta(hours=1)
            time_flag = "-s %s -e %s" % (dt_start.strftime(timestr_fmt),
                                         dt_end.strftime(timestr_fmt))
        else:
            print('[ERROR] Missing one of -rt, -t, or -start/-end flags')
            sys.exit(1)

    path_flag = "-p %s/IO/data" % (script_path)
    if args.data_path is not None:
        path_flag = "-p %s" % (args.data_path)
    else:
        args.data_path = "%s/IO/data" % (script_path)

    arg = "%s %s/get_rap.py %s %s" % (PYTHON, script_path, time_flag, path_flag)
    execute(arg)
    filename = find_nearest(dt, args.data_path)

    print("Closest file in time is: %s" % (filename))
    data = read.read_data(filename)

    # Direct us to the plotting/output functions
    if args.hodo:
        point = args.hodo.split('/')
        point = [[float(point[1]), float(point[0])]]
        create_hodograph(data, point, storm_motion=args.storm_motion,
                         storm_relative=args.storm_relative)
    if args.meso: create_placefiles(data, realtime=args.realtime)
