import pickle, lzma
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.contour import GeoContourSet

from plotconfigs import PLOTCONFIGS as cfg

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
               'srh01km', 'nst', 'oprh', '925fgen', '850fgen', '700fgen', 'snsq', 
               'dgzomega', 'dgzdepth', 'sfctw']
    for var in varlist:
        #plot_data = gaussian_filter(data[var], sigma=1, mode='nearest')
        c = ax.contour(lons, lats, data[var], transform=data_crs,
                       levels=cfg[var]['levels'], colors=cfg[var]['colors'], 
                       linewidths=cfg[var]['linewidths'])
        ax.clabel(c, c.levels, inline=True, fontsize=10)
        plt.savefig(f'./tests/{var}.png', bbox_inches='tight')

        for collection in ax.collections: 
            if isinstance(collection, GeoContourSet):
                collection.remove()
    
    plt.close()

if __name__ == '__main__':
    with lzma.open('./tests/sharppy.xz', 'rb') as f: data = pickle.load(f)
    plot_tests(data[0], data[0]['lons'], data[0]['lats'])