import matplotlib.pyplot as plt
import numpy as np
import geojsoncontour
import json
import os
from datetime import timedelta
from collections import defaultdict
import logging as log

import sharptab.winds as winds
from configs import (SCALAR_PARAMS, VECTOR_PARAMS, BUNDLES, barbconfigs, contourconfigs,
                     plotconfigs)

parent_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
PARAMS = {**SCALAR_PARAMS, **VECTOR_PARAMS}
outdir = "%s/output" % (parent_path)
if not os.path.exists(outdir): os.makedirs(outdir)

def contour(lon, lat, data, time_str, timerange_str, **kwargs):
    """
    Contour plot using geojsoncontour.

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
    timerange_str : string
        Valid time range over which to display in GR

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
    plotinfo = kwargs.get('varname', 'None')
    if levels is not None and colors is not None:
        c = ax.contour(lon, lat, data, levels, colors=colors)
    else:
        c = ax.contour(lon, lat, data)
    geojson = json.loads(geojsoncontour.contour_to_geojson(contour=c, ndigits=2))

    out = []
    out.append('Title: %s | %s\n' % (plotinfo, time_str))
    out.append('RefreshSeconds: 60\n')
    out.append('Font: 1, 14, 1, "Arial"\n')
    out.append('TimeRange: %s\n' % (timerange_str))

    clabs = defaultdict(list) # Store contour labels
    for feature in geojson['features']:
        coords = feature['geometry']['coordinates']
        level = '%s' % (round(feature['properties']['level-value'], 1))
        idx = feature['properties']['level-index']
        if int(levels[idx]) == int(float(level)):
            try:
                lws = list(kwargs['linewidths'])
                linewidth = lws[idx]
            except:
                linewidth = kwargs['linewidths']
        else:
            linewidth = 1

        rgb = hex2rgb(feature['properties']['stroke'])
        out.append('Color: %s 255\n' % (' '.join(rgb)))
        out.append('Line: %s, 0, "%s"\n' % (linewidth, level))

        KNT = 0
        for coord in coords:
            out.append(' %s, %s\n' % (coord[1], coord[0]))
            if KNT % 30 == 0: clabs[level].append([coord[1], coord[0]])
            KNT += 1
        out.append('End:\n\n')

    # Contour labels
    for lev in clabs.keys():
        for val in clabs[lev]:
            if float(lev) >= 9: lev = int(float(lev))
            out.append('Text: %s, %s, 1, "%s", ""\n' % (val[0], val[1], lev))

    plt.close(fig)
    return out

'''
# Kill this one for now. Figure out rotation order for placefile polygons
def contourf(lon, lat, data, time_str, timerange_str, **kwargs):
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
    timerange_str : string
        Valid time range over which to display in GR

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
    plotinfo = kwargs.get('varname', 'None')

    c = ax.contourf(lon, lat, data, levels, colors=colors)
    geojson = json.loads(geojsoncontour.contourf_to_geojson(contourf=c, ndigits=8))

    out = []
    out.append('Title: %s Filled Contour %s\n' % (plotinfo, time_str))
    out.append('RefreshSeconds: 60\n')
    out.append('TimeRange: %s\n' % (timerange_str))
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
    plt.close(fig)
    return out
'''

def write_placefile(arrs, realtime=False):
    """
    Main function controlling the plotting of GR2/Analyst-readable placefiles. Called
    by the primary run.py script.

    Parameters:
    -----------
    arrs : [dictionary]
        List of dictionaries storing values necessary for plotting. This includes
        longitudes, latitudes, valid times, and any 2-d arrays. Each list entry
        corresponds to a new forecast time.
    plotinfo : string

    Other Parameters:
    -----------------
    realtime: bool (default: False)
        If this is a realtime run (True), appends the valid time to the end of the
        filename.

    """
    parms = list(SCALAR_PARAMS.keys()) + list(VECTOR_PARAMS.keys())
    out_dict = defaultdict(list)

    for i in range(len(arrs)):
        arr = arrs[i]
        lon, lat = arr['lons'], arr['lats']
        valid_str = arr['valid_time'].strftime('%H:%MZ %a %b %d %Y')

        # Construct the time range string
        if realtime:
            start = arr['valid_time'] - timedelta(minutes=15)
            end = arr['valid_time'] + timedelta(minutes=15)
        else:
            start = arr['valid_time'] - timedelta(seconds=1799)
            end = arr['valid_time'] + timedelta(seconds=1800)

        timerange_str = "%s %s" % (start.strftime('%Y-%m-%dT%H:%M:%SZ'),
                                   end.strftime('%Y-%m-%dT%H:%M:%SZ'))
        time_str = "%sZ RAP" % (str(arr['cycle_time'].hour).zfill(2))
        for parm in parms:
            matches = [s for s in arr.keys() if parm in s]
            if len(matches) > 1:
                plot_type = 'barb'
                base_configs = barbconfigs
            else:
                plot_type = 'contour'
                base_configs = contourconfigs

            # Add any user-requested ploting configurations to the base. Add the variable
            # name to the configuration dictionary.
            configs = base_configs.copy()
            if parm in plotconfigs:
                config_overrides = plotconfigs[parm]
                configs.update(config_overrides)
            configs['varname'] = PARAMS[parm]

            if plot_type == 'contour':
                out = contour(lon, lat, arr[parm], time_str, timerange_str, **configs)
            #elif plot_type == 'contourf':
            #    out = contourf(lon, lat, arr[parm], time_str, timerange_str, **configs)
            elif plot_type == 'barb':
                out = barbs(lon, lat, arr[parm+'_u'], arr[parm+'_v'], parm, time_str,
                            timerange_str,
                            **configs)
            else:
                raise ValueError("%s is an invalid plot_type entry" % (plot_type))
            out_dict[parm].extend(out)

    for parm in parms:
        output = out_dict[parm]
        if not realtime:
            save_time = "%s-%s" % (arrs[0]['valid_time'].strftime('%Y%m%d%H'),
                                   arrs[-1]['valid_time'].strftime('%Y%m%d%H'))
            out_file = '%s/%s_%s.txt' % (outdir, parm, save_time)
        else:
            out_file = '%s/%s.txt' % (outdir, parm)
        with open(out_file, 'w') as f: f.write("".join(output))

    # If entries exist in the BUNDLES dictionary, output bundled placefiles
    for bundle_name, parameters in BUNDLES.items():
        log.info("Writing bundle: %s with components: %s" % (bundle_name, parameters))
        out_file = '%s/%s.txt' % (outdir, bundle_name)
        with open(out_file, 'w') as f:
            for parm in parameters:
                try:
                    infile = open('%s/%s.txt' % (outdir, parm), 'r')
                    f.write(infile.read())
                except IOError:
                    log.error("%s not found. Skipped during bundling step" % (infile))

def barbs(lon, lat, U, V, parm, time_str, timerange_str, **kwargs):
    iconfile = kwargs.get('windicons', 'Missing: Specify `WIND_ICONDS` in cofigs.py!')
    plotinfo = kwargs.get('varname', 'None')
    skip = kwargs.get('skip', 6)
    out = []
    out.append('Title: %s | %s\n' % (plotinfo, time_str))
    out.append('RefreshSeconds: 60\n')
    out.append('TimeRange: %s\n' % (timerange_str))
    out.append('Color: 255 255 255\n')
    out.append('IconFile: 1, 28, 28, 14, 14, "%s"\n' % (iconfile))
    out.append('Font: 1, 10, 4, "Arial"\n\n')
    for j in range(U.shape[0])[::skip]:
        for i in range(V.shape[1])[::skip]:
            wdir, wspd = winds.comp2vec(float(U[j,i]),
                                        float(V[j,i]))
            if wspd > 0:
                wspd_rounded = 5 * round(wspd/5)
                numref = int(wspd_rounded//5)
                out.append('Object: ' + str(lat[j,i]) + ',' + str(lon[j,i]) +'\n')
                out.append('  Threshold: 999\n')
                out.append('  Icon: 0,0,%s,1,%s,%s\n' % (wdir, np.clip(numref,1,999),
                           round(wspd,1)))
                out.append('End:\n\n')
    return out

def hex2rgb(hex):
    """Convert hexadecimal string to rgb tuple
    """
    h = hex.lstrip('#')
    return tuple(str(int(h[i:i+2], 16)) for i in (0, 2, 4))
