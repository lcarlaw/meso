from collections import OrderedDict

##########################################################################################
# User configurations: Must be edited
#
# Adjust the following variables to point to Python, WGRIB2, and WGET executables on the
# filesystem. See README for notes on WGRIB2 requirements. Specify where to send logfiles.
##########################################################################################
PYTHON = '/Users/leecarlaw/anaconda3/envs/meso/bin/python'
LOG_DIR = '/Users/leecarlaw/scripts/meso/output/logs'
WGRIB2 = '/usr/local/bin/wgrib2'
WGET = '/usr/local/bin/wget'

##########################################################################################
# Plotting configurations
##########################################################################################
SCALAR_PARAMS = {
    'esrh': 'Effective Storm-Relative Helicity',
    'mucape': 'Most-Unstable CAPE (J/kg)',
    'mlcin': 'Mixed-Layer CIN (J/kg)',
    'mlcape': 'Mixed-Layer CAPE (J/kg)',
    'cape3km': '0-3 km MLCAPE (J/kg)',
    'estp': 'Effective SigTor Parameter',
    'srh500': '0-500 m Storm-Relative Helicity',
    'srh01km': '0-1 km Storm-Relative Helicity',
    'lr03km': '0-3 km Lapse Rate (C/km)',
    'mllcl': 'Mixed-Layer LCL (m)'
}

VECTOR_PARAMS = {
    'ebwd': 'Effective Bulk Shear (kt)',
    'shr1': '0-1 km shear (kt)',
    'shr3': '0-3 km shear (kt)',
    'shr6': '0-6 km shear (kt)',
    'shr8': '0-8 km shear (kt)',
    'rm5': 'Bunkers Right Vectors (kt)',
    'lm5': 'Bunkers Left Vectors (kt)',
    'devtor': 'Deviant Tornado Motion (kt)',
}

BUNDLES = {
    'QLCS Tornado': ['cape3km_cf', 'cape3km', 'estp', 'shr3'],
    'ML CAPE-CIN': ['mlcin_cf', 'mlcin', 'mlcape'],
#    'bundle_1': ['cape3km', 'shr3'],
#    'low-level-lapse-rates': ['lr03km', 'lr03km_cf']
}

##########################################################################################
# User overrides for plotting configurations. If not specified, plotting specifications
# default to those in barbconfigs and contourconfigs.
#
# Specifying either the `fill_levels` or `fill_colors` keywords will cause a separate
# contour-filled placefile to be output. These will have `_cf` appended to the filename.
##########################################################################################
# URL or local path to wind icon files
WIND_ICONS = 'https://jupiter-dev.ngrok.io/windicons.png'
SHEAR1_ICONS = 'http://jupiter-dev.ngrok.io/shr1icons.png'

plotconfigs = {
    'mllcl': {
        'levels': [250, 500, 1000, 1500, 2000, 3000]
    },

    'lr03km': {
        'colors': 'k',
        'levels': [6, 7, 8, 9, 10],
        'linewidths': 0.75,
        'fill_levels': [6.5, 999],
        'fill_colors': ['#f1a95d', '#f1a95d']
    },

    'esrh': {
        'colors': ['#81b6f7', '#81b6f7', '#3c6193', '#3c6193', '#3c6193', '#3c6193'],
        'levels': [50, 100, 200, 300, 400, 500, 600, 700, 800],
        'linewidths': [1, 1, 2, 2, 2, 3, 3, 3, 3]
    },

    'srh01km': {
        'colors': ['#81b6f7', '#81b6f7', '#3c6193', '#3c6193', '#3c6193', '#3c6193'],
        'levels': [50, 100, 200, 300, 400, 500, 600, 700, 800],
        'linewidths': [1, 1, 2, 2, 2, 3, 3, 3, 3]
    },

    'srh500': {
        'colors': ['#81b6f7', '#81b6f7', '#3c6193', '#3c6193', '#3c6193', '#3c6193'],
        'levels': [25, 50, 100, 200, 300, 400, 500, 600, 800],
        'linewidths': [1, 1, 2, 2, 2, 3, 3, 3, 3]
    },

    'shr1': {
        'windicons': SHEAR1_ICONS,
    },

    'mucape': {
        'colors': ['#ea908c', '#ea908c', '#da453a', '#da453a', '#da453a', '#c03c32',
                   '#c03c32', '#c03c32', '#c03c32'],
        'levels': [100, 250, 500, 1000, 2000, 3000, 4000, 5000, 6000],
        'linewidths': [1, 1, 1, 2, 2, 3, 3, 3, 3]
    },

    'mlcape': {
        'colors': ['#ea908c', '#ea908c', '#da453a', '#da453a', '#da453a', '#c03c32',
                   '#c03c32', '#c03c32', '#c03c32'],
        'levels': [100, 250, 500, 1000, 2000, 2500, 3000, 3500, 4000],
        'linewidths': [1, 1, 1, 2, 2, 3, 3, 3, 3]
    },

    'cape3km': {
        #'colors': ['#dd564e', '#dd564e', '#dd564e', '#bb2d1d', '#bb2d1d',
        #           '#bb2d1d', '#841f18', '#841f18'],
        'colors': 'k',
        'levels': [25, 75, 100, 150, 300, 500],
        'linewidths': [0.5, 0.5, 1, 1, 1.5, 1.5],
        'fill_levels': [50, 9999],
        'fill_colors': '#f1b1bc'
    },

    'mlcin': {
        'colors': ['#6beaea', '#4eb0e7', '#0000cf', '#0000cf', '#0000cf'],
        'levels': [25, 75, 150, 300],
        'linewidths': [1, 1, 1, 2, 2],
        'fill_levels': [25, 9999],
        'fill_colors': ['#4eb0e7', '#4eb0e7'],
    },

    'estp': {
        'levels': [0.25, 0.5, 1, 2, 4, 6, 8, 10],
        'colors': ['#ec904a', '#ec904a', '#ec904a', '#e94639', '#c23f34',
                   '#841f18', '#957cca', '#e951f5', '#e951f5'],
        'linewidths': [1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    },
}

##########################################################################################
# Base contour and wind barb configurations. These will be overriden by any corresponding
# entries in the PLOTCONFIGS dictionary
##########################################################################################
barbconfigs = {
    'skip': 6,
    'windicons': WIND_ICONS,
}

contourconfigs = {
    'colors': ['#dd564e', '#dd564e', '#dd564e', '#dd564e', '#bb2d1d', '#bb2d1d',
               '#bb2d1d', '#841f18', '#841f18'],
    'levels': [100, 250, 500, 1000, 2000, 3000, 4000, 5000, 6000],
    'linewidths': [1, 1, 1, 2, 2, 2, 3, 4, 4],
    'fill_levels': None,
    'fill_colors': None,
}

##########################################################################################
# Filtering specifications
##########################################################################################
FILTER_SPECS = {
    'mlcin': {
        'mlcape': ['>', 5],
    },

    'ebwd': {
        'ebwd': ['>', 20],
    },

    'shr1': {
        'shr1': ['>', 15],
    },

    'shr3': {
        'shr3': ['>', 20],
    },

    'shr6': {
        'shr3': ['>', 20],
    },

    'shr8': {
        'shr3': ['>', 25],
    },

    'rm5': {
        'mucape': ['>', 100],
        'ebwd': ['>', 25],
    },

    'lm5': {
        'mucape': ['>', 100],
        'ebwd': ['>', 25],
    },

    'devtor': {
        'mucape': ['>', 100],
        'ebwd': ['>', 25],
    },
}

##########################################################################################
# Download configurations
#
# You likely won't need (or want) to change these. Download priority is set by the order
# of the dictionary keys in the DATA_SOURCES variable.
##########################################################################################
TIMEOUT = 200       # Seconds after which to timeout the data download function
MAXSECONDS = 1800   # Number of seconds after which to abort data download in run.py
MINSIZE = 5         # Grib files under this size (MB) will result in a download error
SIGMA = 1.5         # For smoothing function. Larger = more smoothing, but amplitude loss
ALPHA = 50          # Alpha level for filled placefiles. 0 = transparent; 255 = opaque

DATA_SOURCES = OrderedDict({
    'NOMADS': 'https://nomads.ncep.noaa.gov/pub/data/nccf/com',
    'GOOGLE': 'https://storage.googleapis.com',
    'FTPPRD': 'https://ftpprd.ncep.noaa.gov/data/nccf/com',
    'THREDDS': 'https://www.ncei.noaa.gov/thredds/fileServer',
})

GOOGLE_CONFIGS = {
    'RAP': 'rapid-refresh',
    'HRRR': 'high-resolution-rapid-refresh'
}

THREDDS_CONFIGS = {
    'RAP-current': 'model-rap130anl',
    'RAP': 'model-rap130anl-old',
    'RUC': 'model-ruc130anl'
}

vars = ':(HGT|TMP|SPFH|UGRD|VGRD|PRES):'
grid_info = 'lambert:262.5:38.5 -105:240:13545.000000 25:210:13545.000000'
#grid_info = 'lambert:262.5:38.5 -120:360:13545.000000 23:210:13545.000000'
