""" 
Unit test for output data. Reference is RAP data from 2021-12-11/03 model cycle, fhr 1 
valid at 2021-12-11/04 UTC.  
"""

import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pickle, lzma
import numpy as np 
import plot.plots as plotting

# Read in data to test
with lzma.open('sharppy.xz', 'rb') as f: data = pickle.load(f)

# Reference data for comparison
with lzma.open('reference_data/sharppy.xz', 'rb') as f: reference = pickle.load(f)

plot_vars = ['esrh', 'estp', 'mucape', 'mlcape', 'mlcin', 'cape3km', 'lr03km', 
             'srh01km', 'nst', 'oprh', '925fgen', '850fgen', '700fgen', 'snsq', 
             'dgzomega', 'dgzdepth', 'sfctw']

for var in plot_vars:
    test = data[0][var]
    ref = reference[0][var]
    test = np.nan_to_num(test, nan=0.)
    ref = np.nan_to_num(ref, nan=0.)
    try:
        assert (test == ref).all()
        print(f'Test {var}: {var} - PASS')
    except AssertionError:
        print(f'Test {var}: {var} - FAIL')

#plotting.create_map_output(data[0], data[0]['lons'], data[0]['lats'], plot_dir='./')