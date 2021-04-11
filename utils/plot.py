import matplotlib.pyplot as plt
import numpy as np
import geojsoncontour
import json

import sharptab.winds as winds
from utils.plot_configs import metadata as meta

import os
script_path = os.path.dirname(os.path.realpath(__file__))
parent_path = os.path.dirname(script_path)

ALPHA = 90
outdir = "%s/output" % (parent_path)
WINDICONS = 'https://jupiter-dev.ngrok.io/windicons.png'
def contour(lon, lat, data, time_str, **kwargs):
    """Contour plot using geojsoncontour.

    Parameters:
    -----------
    lon : array
        2-D array of longitudes. Must be same shape as data
    lat : array
        2-D array of latitudes. Must be same shape as data
    data : array [N, M]
        Values over which contour is drawn.
    time_str : string
        Valid time for this plot. Included in the placefile title.

    Other Parameters:
    -----------------
    levels : list, array
        Contour levels to plot.
    linewidths: list, array
        Linewidths corresponding to each contour level
    colors : color string (hexademicals)
        Colors corresponding to each contour level
    plotinfo : string
        Brief description of the plot

    Returns:
    --------
    out : list
        List of strings, each corresponding to a new line for the placefile
    """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    levels = kwargs.get('levels')
    colors = kwargs.get('colors')
    plotinfo = kwargs.get('plotinfo', 'None')
    if levels is not None and colors is not None:
        c = ax.contour(lon, lat, data, levels, colors=colors)
    else:
        c = ax.contour(lon, lat, data)
    geojson = json.loads(geojsoncontour.contour_to_geojson(contour=c, ndigits=2))

    out = []
    out.append('Title: %s %s\n' % (plotinfo, time_str))
    out.append('RefreshSeconds: 60\n')
    for feature in geojson['features']:
        coords = feature['geometry']['coordinates']
        level = '%s' % (feature['properties']['level-value'])
        idx = feature['properties']['level-index']
        if int(levels[idx]) == int(float(level)):
            try:
                lws = list(kwargs['linewidths'])
                linewidth = lws[idx]
            except:
                linewidth = kwargs['linewidths']
        else:
            linewidth = 1
        #linewidth += 2
        rgb = hex2rgb(feature['properties']['stroke'])
        out.append('Color: %s 255\n' % (' '.join(rgb)))
        out.append('Line: %s, 0, "%s"\n' % (linewidth, level))
        for coord in coords:
            out.append(' %s, %s\n' % (coord[1], coord[0]))
        out.append('End:\n\n')
    return out

def contourf(lon, lat, data, time_str, **kwargs):
    """Contour-filled plot using geojsoncontour.

    Parameters:
    -----------
    lon : array
        2-D array of longitudes. Must be same shape as data
    lat : array
        2-D array of latitudes. Must be same shape as data
    data : array [N, M]
        Values over which contour fill is drawn
    time_str : string
        Valid time for this plot. Included in the placefile title

    Other Parameters:
    -----------------
    levels : list, array
        Number and positions of the contour lines / regions
    colors : color string (hexademicals)
        Colors corresponding to each contour-fill level
    plotinfo : string
        Brief description of the plot

    Returns:
    --------
    out : list
        List of strings, each corresponding to a new line for the placefile
    """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    levels = kwargs.get('levels')
    colors = kwargs.get('colors')
    plotinfo = kwargs.get('plotinfo', 'None')

    c = ax.contourf(lon, lat, data, levels, colors=colors)
    geojson = json.loads(geojsoncontour.contourf_to_geojson(contourf=c, ndigits=8))

    out = []
    out.append('Title: %s Filled Contour %s\n' % (plotinfo, time_str))
    out.append('RefreshSeconds: 60\n')
    for feature in geojson['features']:
        rgb = hex2rgb(feature['properties']['fill'])
        groups = feature['geometry']['coordinates']
        for group in groups:
            out.append('Polygon:\n')
            first_line = True
            for item in group[0]:
                if first_line:
                    COLOR ='%s, %s\n' % (', '.join(rgb), ALPHA)
                else:
                    COLOR = ''
                out.append(' %s, %s, %s\n' % (item[1], item[0], COLOR))
                first_line = False
            out.append('End:\n')
            out.append('\n')
    return out

def write_placefile(arr, plotinfo, realtime=False):
    """Main function controlling the plotting of GR2/Analyst-readable placefiles. Called
    by the primary run.py script.

    Parameters:
    -----------
    arr : dictionary
        Dictionary storing values necessary for plotting. This includes longitudes,
        latitudes, valid times, and any 2-d arrays.
    plotinfo : string

    Other Parameters:
    -----------------
    realtime: bool (default: False)
        If this is a realtime run (True), appends the valid time to the end of the
        filename.

    """
    parms = plotinfo.keys()
    lon, lat = arr['lons'], arr['lats']
    save_time = arr['valid_time'].strftime('%Y%m%d%H')
    valid_str = arr['valid_time'].strftime('%HZ %a %b %d %Y')
    time_str = "%sZ HRRR | F%s Valid: %s" % (str(arr['cycle_time'].hour).zfill(2),
                                            str(arr['fhr']).zfill(2), valid_str)
    for parm in parms:
        meta[parm]['plotinfo'] = plotinfo[parm]
        plot_type = meta[parm]['plot_type']
        if plot_type == 'contour':
            out = contour(lon, lat, arr[parm], time_str, **meta[parm])
        elif plot_type == 'contourf':
            out = contourf(lon, lat, arr[parm], time_str, **meta[parm])
        elif plot_type == 'barb':
            out = barbs(lon, lat, arr['vectors'], parm, time_str, **meta[parm])
        else:
            raise ValueError("%s is an invalid plot_type entry" % (plot_type))

        if not realtime:
            out_file = '%s/%s_%s.txt' % (outdir, parm, save_time)
        else:
            out_file = '%s/%s.txt' % (outdir, parm)
        with open(out_file, 'w') as f: f.write("".join(out))

def barbs(lon, lat, data, parm, time_str, **kwargs):
    plotinfo = kwargs.get('plotinfo', 'None')
    skip = kwargs.get('skip', 6)

    out = []
    out.append('Title: %s %s\n' % (plotinfo, time_str))
    out.append('RefreshSeconds: 60\n')
    out.append('Color: 255 255 255\n')
    out.append('IconFile: 1, 28, 28, 3, 28, "%s"\n' % (WINDICONS))
    out.append('Font: 1, 10, 4, "Arial"\n\n')
    for j in range(data[parm + '_u'].shape[0])[::skip]:
        for i in range(data[parm + '_u'].shape[1])[::skip]:
            wdir, wspd = winds.comp2vec(float(data[parm + '_u'][j,i]),
                                        float(data[parm + '_v'][j,i]))
            if wspd > 20:
                wspd_rounded = 5 * round(wspd/5)
                numref = str(int(wspd_rounded//5))
                out.append('Object: ' + str(lat[j,i]) + ',' + str(lon[j,i]) +'\n')
                out.append('  Threshold: 999\n')
                out.append('  Icon: 0,0,%s,1,%s,%s\n' % (wdir, numref, round(wspd,1)))
                out.append('End:\n\n')
    return out

def hex2rgb(hex):
    """Convert hexadecimal string to rgb tuple
    """
    h = hex.lstrip('#')
    return tuple(str(int(h[i:i+2], 16)) for i in (0, 2, 4))
