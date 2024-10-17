import pickle, lzma
from numba import set_num_threads
from calc import derived
from sharptab.winds import vec2comp
from sharptab.constants import KTS2MS
from calc.compute import worker 

from numba.typed import List
from plotconfigs import SCALAR_PARAMS, VECTOR_PARAMS
from plotconfigs import PLOTCONFIGS as cfg
from time import time

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.contour import GeoContourSet

import numpy as np 
from scipy.ndimage import gaussian_filter

def cache_funcs(data):
    pres = data['pres']
    tmpc = data['tmpc']
    dwpc = data['dwpc']
    wspd = data['wspd']
    wdir = data['wdir']
    hght = data['hght']
    lons = data['lons']
    lats = data['lats']

    # Vorticity calculations for NST parameter. 0th index from hybrid files is ~10m agl.
    u, v = vec2comp(wdir[0,:,:], wspd[0,:,:]*KTS2MS)
    vort = derived.vorticity(u, v, lons, lats)

    print("=====================================")
    print("In to slow compilation loop 1")
    t1 = time()
    results = worker(pres, tmpc, hght, dwpc, wspd, wdir, vort, 
                     List(SCALAR_PARAMS.keys()), List(VECTOR_PARAMS.keys()))
    
    print("In to fast loop 2")
    t2 = time()
    results = worker(pres, tmpc, hght, dwpc, wspd, wdir, vort,
                     List(SCALAR_PARAMS.keys()), List(VECTOR_PARAMS.keys()))
    t3 = time()

    print(f"Slow compilation loop 1: {t2-t1} seconds")
    print(f"Fast loop 2: {t3-t2} seconds")
    print("=====================================")

    return results

def plot_background(ax,):
    ax.set_extent([-120, -74, 25, 47.5])
    ax.coastlines(resolution='50m', color='blue', linewidths=1.2)
    ax.add_feature(cfeature.BORDERS, linewidths=1.2, edgecolor='black')
    land_mask = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                            edgecolor='face',
                                            facecolor=cfeature.COLORS['land'])
    sea_mask = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                            edgecolor='face',
                                            facecolor=cfeature.COLORS['water'])
    lake_mask = cfeature.NaturalEarthFeature('physical', 'lakes', '50m',
                                            edgecolor='face',
                                            facecolor=cfeature.COLORS['water'])
    #ax.add_feature(USCOUNTIES.with_scale('5m'), linewidth=0.5, edgecolor='#563d22')
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5)
    ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=1.2, edgecolor='k')
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(sea_mask, zorder=0)
    ax.add_feature(land_mask, zorder=0)
    ax.add_feature(lake_mask, zorder=0)
    return ax

def plot_tests(data, lons, lats):
    print("Plotting test data...")
    data_crs = ccrs.PlateCarree()
    crs = ccrs.LambertConformal(central_longitude=-99, central_latitude=35.0, 
                                standard_parallels=[35.])
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18, 12), dpi=200, 
                        subplot_kw={'projection': crs})
    plt.subplots_adjust(wspace=0, hspace=0)
    plot_background(ax)
    ax.set_adjustable('datalim')

    varlist = ['esrh', 'estp', 'mucape', 'mlcape', 'mlcin', 'cape3km', 'lr03km', 
               'srh01km', 'nst']
    for var in varlist:
        plot_data = gaussian_filter(data[var], sigma=1, mode='nearest')
        c = ax.contour(lons, lats, plot_data, transform=data_crs,
                       levels=cfg[var]['levels'], colors=cfg[var]['colors'], 
                       linewidths=cfg[var]['linewidths'])
        ax.clabel(c, c.levels, inline=True, fontsize=10)
        plt.savefig(f'./tests/{var}.png', bbox_inches='tight')

        for collection in ax.collections: 
            if isinstance(collection, GeoContourSet):
                collection.remove()
    
    plt.close()

if __name__ == '__main__':
    set_num_threads(8)
    with lzma.open('./tests/standard.xz', 'rb') as f: data = pickle.load(f)

    results = cache_funcs(data)
    plot_tests(results, data['lons'], data['lats'])

    