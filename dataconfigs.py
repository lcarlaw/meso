"""
Variables and references used by get_data.py.  This file is tracked, so changes will
reflect during git updates.
"""

from collections import OrderedDict

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

# Variables for WGRIB2 download. Used in get_data.py, and only applicable for realtime 
# data, or files from Google cloud storage. NCEI datafiles are downloaded in their 
# entirety since no .idx reference data exists. 
data_vars = ':(HGT|TMP|SPFH|UGRD|VGRD|PRES|VVEL):'

# 20-km CONUS - used by WGRIB2 in get_data.py to regrid data to a standard reference.
grid_info = 'lambert:262.5:38.5 -120:245:20318.000000 23:145:20318.000000'