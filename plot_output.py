"""
Attempt to speed up plotting by pre-rendering the map background as a transparent PNG 
and then plotting this under the meteorological data. Was only able to get this to work 
using the ccrs.PlateCarree (Mercator) projection, as Cartopy and/or Matplotlib performed
some internal bbox_inches alterations for the Lambert Conformal projection -- unable to 
get the two images to align.  
"""

import pickle, lzma
from pathlib import Path 

import matplotlib.pyplot as plt
from cartopy.mpl.contour import QuadContourSet
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from metpy.plots import USCOUNTIES

import matplotlib.image as mpimg
from plotconfigs import PLOTCONFIGS, SCALAR_PARAMS

# Coordinate reference systems
data_crs = ccrs.PlateCarree()
#plot_crs = ccrs.LambertConformal(central_longitude=-99, central_latitude=35.0,
#                                 standard_parallels=[35.])
plot_crs = ccrs.PlateCarree()

# Map extent
#EXTENT = [-120, -74, 25, 47.5]     # Lambert Comformal
EXTENT = [-125, -67, 25, 50]        # Mercator

FIGSIZE = (24, 12) 
DPI = 100

# Variables to plot -- must match the dictionary keys in SCALAR_PARAMS from plotconfigs
varlist = ['esrh', 'estp', 'mucape', 'mlcape', 'mlcin', 'cape3km', 'lr03km', 'srh01km', 
           'nst', 'oprh', '925fgen', '850fgen', '700fgen', 'snsq',  'dgzomega', 'dgzdepth', 
           'sfctw', '925T', '850T', '700T']

def make_background(plot_dir):
    background_image = f'{plot_dir}/background.png'

    if not Path(background_image).is_file():
        fig = plt.figure(figsize=FIGSIZE, dpi=DPI)
        ax = fig.add_axes([0, 0, 1, 1], projection=plot_crs)
        ax.set_extent(EXTENT, crs=data_crs)
        ax.set_aspect('auto')
        ax.axis('off')
        #ax.set_autoscale_on(False)

        land_mask = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                                 edgecolor='face',
                                                 facecolor=cfeature.COLORS['land'])
        sea_mask = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                                 edgecolor='face',
                                                 facecolor=cfeature.COLORS['water'])
        lake_mask = cfeature.NaturalEarthFeature('physical', 'lakes', '50m',
                                                 edgecolor='face',
                                                 facecolor=cfeature.COLORS['water'])
        
        # Add features        
        ax.coastlines(resolution='50m', color='blue', linewidth=1.2)
        ax.add_feature(cfeature.BORDERS, linewidth=1.2, edgecolor='black')
        ax.add_feature(USCOUNTIES.with_scale('20m'), linewidth=0.5, edgecolor='#9b9b9b')
        ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=1.2, edgecolor='k')
        ax.add_feature(sea_mask, zorder=0)
        ax.add_feature(land_mask, zorder=0)
        ax.add_feature(lake_mask, zorder=0)

        # Save transparent image
        plt.savefig(f'{plot_dir}/background.png', transparent=True, dpi=DPI, 
                    bbox_inches=None) 
        plt.close(fig)


def create_map_output(data, plot_dir='./tests'):
    lons, lats = data[0]['lons'], data[0]['lats']
    make_background(plot_dir)

    fig = plt.figure(figsize=FIGSIZE, dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])

    bg_img = mpimg.imread(f'{plot_dir}/background.png')
    ax.imshow(bg_img, origin='upper', extent=EXTENT, zorder=0)
    ax.set_xlim(EXTENT[0], EXTENT[1])
    ax.set_ylim(EXTENT[2], EXTENT[3])
    ax.axis('off')
    ax.set_aspect('auto')
    for var in varlist:
        for fhr in range(len(data)):
            plot_data = data[fhr]
            valid_time = plot_data['valid_time'].strftime('%a %Y-%m-%d %HZ')
            valid_time_info = f"F{str(plot_data['fhr']).zfill(2)} Valid: {valid_time}"
            init_time_info = f"Init: {plot_data['cycle_time'].strftime('%a %Y-%m-%d %HZ')}"

            c = ax.contour(lons, lats, plot_data[var], levels=PLOTCONFIGS[var]['levels'],
                           colors=PLOTCONFIGS[var]['colors'], 
                           linewidths=PLOTCONFIGS[var]['linewidths'])
            ax.clabel(c, inline=True, fontsize=12)

            t1 = ax.text(EXTENT[0]+0.1, EXTENT[3]-0.1, f"{SCALAR_PARAMS[var]}",
                         ha='left', va='top', fontsize=16, fontweight='bold', 
                         color='black', bbox=dict(facecolor='white', edgecolor='none'))
            t2 = ax.text(EXTENT[0]+0.1, EXTENT[3]-0.8, valid_time_info,
                         ha='left', va='top', fontsize=16,
                         color='black', bbox=dict(facecolor='white', edgecolor='none'))
            t3 = ax.text(EXTENT[1]-0.1, EXTENT[3]-0.1, init_time_info,
                         ha='right', va='top', fontsize=16,
                         color='black', bbox=dict(facecolor='white', edgecolor='none'))

            vt_save = plot_data['valid_time'].strftime('%Y%m%d%HZ')
            plt.savefig(f'{plot_dir}/{var}_{vt_save}.png', dpi=DPI, bbox_inches=None) 
            for collection in ax.collections: 
                if isinstance(collection, QuadContourSet):
                    collection.remove()
            t1.remove()
            t2.remove()
            t3.remove()

    plt.close(fig)

'''
# Code from original plot_output script. Very slow on savefig step with counties on due 
# to extra vertices. 
def make_background(ax, plot_counties):
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
    if plot_counties:
        from metpy.plots import USCOUNTIES
        ax.add_feature(USCOUNTIES.with_scale('20m'), linewidth=0.5, 
                        edgecolor='#9b9b9b')
    
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.5)
    ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=1.2, edgecolor='k')
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(sea_mask, zorder=0)
    ax.add_feature(land_mask, zorder=0)
    ax.add_feature(lake_mask, zorder=0)
    return ax

def create_map_output(data, lons, lats, plot_counties=False, plot_dir='./tests',):
    """
    Plots parameters on a map.  
    """
    # Variables to plot. These are dict keys in SCALAR_PARAMS in the plotconfigs.py file.
    varlist = ['esrh', 'estp', 'mucape', 'mlcape', 'mlcin', 'cape3km', 'lr03km', 
               'srh01km', 'nst', 'oprh', '925fgen', '850fgen', '700fgen', 'snsq', 
               'dgzomega', 'dgzdepth', 'sfctw', '925T', '850T', '700T']

    data_crs = ccrs.PlateCarree()
    crs = ccrs.LambertConformal(central_longitude=-99, central_latitude=35.0, 
                                standard_parallels=[35.])
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18, 12), dpi=200, 
                         subplot_kw={'projection': crs})
    plt.subplots_adjust(wspace=0, hspace=0)
    make_background(ax, plot_counties)
    ax.set_adjustable('datalim')

    valid_time = data['valid_time'].strftime('%a %Y-%m-%d %HZ')
    valid_time_info = f"F{str(data['fhr']).zfill(2)} Valid: {valid_time}"
    init_time_info = f"Init: {data['cycle_time'].strftime('%a %Y-%m-%d %HZ')}"
    for var in varlist:
        c = ax.contour(lons, lats, data[var], transform=data_crs,
                       levels=PLOTCONFIGS[var]['levels'], colors=PLOTCONFIGS[var]['colors'], 
                       linewidths=PLOTCONFIGS[var]['linewidths'])
        ax.clabel(c, c.levels, inline=True, fontsize=10)
        
        t1 = ax.annotate(f"{SCALAR_PARAMS[var]}", xy=(0, 1.02), va='bottom', 
                         xycoords='axes fraction', fontsize=16, fontweight='bold')
        t2 = ax.annotate(valid_time_info, xy=(0, 1), va='bottom',
                         xycoords='axes fraction', fontsize=14)
        t3 = ax.annotate(init_time_info, xy=(1, 1), va='bottom', ha='right',
                         xycoords='axes fraction', fontsize=14)

        plt.savefig(f'{plot_dir}/{var}.png', bbox_inches='tight')

        for collection in ax.collections: 
            if isinstance(collection, GeoContourSet):
                collection.remove()
        t1.remove()
        t2.remove()
        t3.remove()
    
    plt.close()
'''

filename = '../tmp/tests/sharppy.xz'
with lzma.open(filename, 'rb') as f: data = pickle.load(f)
create_map_output(data, plot_dir='../tmp/tests')

