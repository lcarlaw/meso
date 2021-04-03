"""Simple script to read in raw table of WSR88D locations and output a python-readable
dictionary of lat/lon locations"""

import pandas as pd
import pickle

def dms2dd(degrees, minutes, seconds, direction):
    """Grabbed from somewhere online. Can't remember now...
    """
    dd = float(degrees) + float(minutes)/60 + float(seconds)/(60*60);
    if direction == 'S' or direction == 'W':
        dd *= -1
    return dd;

df = pd.read_csv('wsr88d.txt', delimiter='\t')
radars = {}
for idx, row in df.iterrows():
    lat, lon = row['LOCATION'].split('/ ')
    lat = dms2dd(lat[0:2], lat[2:4], lat[4:6], 'N')

    direction = 'W'
    if len(lon) == 8: direction = lon[7]
    lon = dms2dd(lon[0:3], lon[3:5], lon[5:7], direction)
    radars[row['STATION ID']] = {'lat': lat, 'lon': lon}

with open('wsr88d.pkl', 'wb') as f: pickle.dump(radars, f)
