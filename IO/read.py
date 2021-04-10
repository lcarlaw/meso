import pygrib
import numpy as np

import sharptab.winds as winds
import sharptab.thermo as thermo
from sharptab.constants import MS2KTS, ZEROCNK

def make_data_monotonic(data):
    """For data in isobaric coordinates, need to ensure pressure monotonically decreases
    and height monotonically increases with height

    Parameters:
    -----------
    data : dictionary
        Dictionary containing the processed data cube.

    Returns:
    data : dictionary
        Monotonic dataset with the same dimensions as the input. NaNs inserted where data
        has been shifted.

    """
    original_data = data
    full_keys = list(original_data.keys())
    keys = ['tmpc', 'dwpc', 'hght', 'wdir', 'wspd', 'pres']
    arr_shape = data['tmpc'].shape
    for j in range(arr_shape[1]):
        for i in range(arr_shape[2]):
            shift = np.where(data['pres'][:,j,i] < data['pres'][0,j,i])[0][0]
            for key in keys:
                last_val = data[key][-1,j,i]
                data[key][:,j,i] = np.append(np.insert(data[key][shift:,j,i], 0,
                                             data[key][0,j,i]), [last_val]*(shift-1))

    remaining_keys = list(np.setdiff1d(full_keys, keys))
    for key in remaining_keys:
        data[key] = original_data[key]
    return data

def read_data(filename):
    grb = pygrib.open(filename)
    sfc = grb.select(name='Orography')[0]
    lats, lons = sfc.latlons()

    temperatures = grb.select(name='Temperature')[0]
    order = 1
    delta = 1

    # For realtime runs, this is easy if we use the native hybrid-sigma coordinate files.
    # By definition, we're guaranteed each vertical slice is higher than the one below.
    if temperatures.typeOfLevel == 'hybrid':
        native = True
        level_type = 'hybrid'
        temperatures = grb.select(name='Temperature', typeOfLevel=level_type)
        heights = grb.select(name='Geopotential Height', typeOfLevel=level_type)
        pressure = grb.select(name='Pressure', typeOfLevel=level_type)
        specific_humidity = grb.select(name='Specific humidity', typeOfLevel=level_type)

    # If this is isobaric data, we need to "drop in" the surface fields and mask out
    # data that's below the ground. Additionally, we need to check if the data is stored
    # top-down (changes with certain RAP cycles)
    elif temperatures.typeOfLevel == 'isobaricInhPa':
        native = False
        level_type = 'isobaricInhPa'
        heights = grb.select(name='Geopotential Height', typeOfLevel=level_type)
        rh = grb.select(name='Relative humidity', typeOfLevel=level_type)
        psfc = grb.select(name='Surface pressure')[0].values
        t2m = grb.select(name='2 metre temperature')[0].values
        td2m = grb.select(name='2 metre dewpoint temperature')[0].values
        u10m = grb.select(name='10 metre U wind component')[0].values
        v10m = grb.select(name='10 metre V wind component')[0].values

        # Determine how the data is oriented vertically
        if heights[0].values[0,0] > heights[-1].values[0,0]:
            order = -1
            delta = 0

    temperatures = grb.select(name='Temperature', typeOfLevel=level_type)
    uwind = grb.select(name='U component of wind', typeOfLevel=level_type)
    vwind = grb.select(name='V component of wind', typeOfLevel=level_type)

    arr_shape = heights[0].values.shape
    num_levels = len(heights)

    if not native: num_levels += 1
    hght = np.zeros((num_levels, arr_shape[0], arr_shape[1]))
    tmpc = np.full_like(hght, 0)
    pres = np.full_like(hght, 0)
    dwpc = np.full_like(hght, 0)
    wdir = np.full_like(hght, 0)
    wspd = np.full_like(hght, 0)

    # First build the height array. Do this in order, starting with the surface level
    # Note that array will still be flipped if order = -1, we'll take care of that later.
    knt = 0
    for lev in range(num_levels)[::order]:
        # Isobaric coordinates
        if not native:
            if knt == 0:
                pres[lev] = psfc / 100.
                hght[lev] = sfc.values
                tmpc[lev] = t2m - ZEROCNK
                dwpc[lev] = td2m - ZEROCNK
                wdir[lev], wspd[lev] = winds.wind_vecs(u10m, v10m)
            else:
                pres[lev] = np.full_like(heights[0].values, heights[lev-delta].level)
                tmpc[lev] = temperatures[lev-delta].values - 273.15
                #wdir[lev], wspd[lev] = winds.wind_vecs(uwind[lev-delta].values,
                #                                       vwind[lev-delta].values)
                hght[lev] = heights[lev-delta].values
                dwpc[lev] = thermo.dewpoint_from_rh(temperatures[lev-delta].values,
                                                    rh[lev-delta].values/100)

        # Native coordinates
        else:
            hght[lev] = heights[lev].values
            tmpc[lev] = temperatures[lev].values - 273.15
            wdir[lev], wspd[lev] = winds.wind_vecs(uwind[lev].values, vwind[lev].values)
            pres[lev] = pressure[lev].values / 100.
            dwpc[lev] = thermo.dewpoint_from_q(specific_humidity[lev].values,
                                               temperatures[lev].values,
                                               pres[lev])
        knt += 1

    tmpc = np.array(tmpc[::order], dtype='float64')
    dwpc = np.array(dwpc[::order], dtype='float64')
    hght = np.array(hght[::order], dtype='float64')
    wdir = np.array(wdir[::order], dtype='float64')
    wspd = np.array(wspd[::order], dtype='float64') * MS2KTS
    pres = np.array(pres[::order], dtype='float64')

    data = {
        'tmpc': tmpc,
        'dwpc': dwpc,
        'hght': hght,
        'wdir': wdir,
        'wspd': wspd,
        'pres': pres,
        'lons': lons,
        'lats': lats,
        'cycle_time': sfc.analDate,
        'valid_time': sfc.validDate,
        'fhr': sfc.forecastTime
    }

    if not native: data = make_data_monotonic(data)
    return data
