import pickle, lzma
import cartopy.crs as ccrs
from cartopy.mpl.contour import GeoContourSet
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

from plotconfigs import PLOTCONFIGS, SCALAR_PARAMS

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

filename = './tests/sharppy.xz'
with lzma.open(filename, 'rb') as f: data = pickle.load(f)
create_map_output(data[0], data[0]['lons'], data[0]['lats'], plot_counties=False, 
                  plot_dir='/Users/leecarlaw/Desktop/test')