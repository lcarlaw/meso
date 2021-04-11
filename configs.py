# The USER must edit these file paths
PYTHON = '/Users/leecarlaw/anaconda3/envs/meso/bin/python'
WGRIB2 = '/usr/local/bin/wgrib2'
WGET = '/usr/local/bin/wget'
TIMEOUT = 120 # Seconds after which to timeout the data download function

# ----------------------------------------------------------------------------------------
# Plotting configs
# ----------------------------------------------------------------------------------------
plotinfo = {
    'cape3km': '0-3 km MLCAPE (J/kg)',
    'estp': 'Significant Tornado (Effective)',
    'lr3km_cf': '0-3 km Lapse Rate (Filled > 6 C/km)',
    'lr3km': '0-3 km Lapse Rate (> 7 C/km)',
    'mucape': 'Most-unstable CAPE (J/kg)',
    'mlcape': 'Mixed-layer CAPE (J/kg)',
    'mlcin' : 'Mixed-layer CIN (J/kg)',
    'esrh': 'Effective Storm-Relative Helicity',
    'tts': 'Tornadic Tilting and Stretching Parameter',
    'shr3': '0-3 km shear vector (kt)',
    'ebwd': 'Effective bulk shear (kt)',
    'rm5': 'Right Bunkers motion vector',
    'lm5': 'Left Bunkers motion vector',
}

# ----------------------------------------------------------------------------------------
# Download configs. You likely won't need (or want) to change these. Priority is set by
# the order of the dictionary keys in the sources variable.
# ----------------------------------------------------------------------------------------
sources = {
    'NOMADS': 'https://nomads.ncep.noaa.gov/pub/data/nccf/com',
    'GOOGLE': 'https://storage.googleapis.com',
    'FTPPRD': 'https://ftpprd.ncep.noaa.gov/data/nccf/com',
    'THREDDS': 'https://www.ncei.noaa.gov/thredds/fileServer',
}

google_configs = {
    'RAP': 'rapid-refresh',
    'HRRR': 'high-resolution-rapid-refresh'
}

thredds_configs = {
    'RAP-current': 'model-rap130anl',
    'RAP': 'model-rap130anl-old',
    'RUC': 'model-ruc130anl'
}

vars = ':(HGT|TMP|SPFH|UGRD|VGRD|PRES):'
grid_info = 'lambert:262.5:38.5 -105:240:13545.000000 25:210:13545.000000'
