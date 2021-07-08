from numba import njit, prange
from numba.typed import List, Dict
from numba.core import types

from collections import defaultdict
import numpy as np
from scipy.ndimage import gaussian_filter

from multiprocessing import Pool
import multiprocessing

from configs import SCALAR_PARAMS, VECTOR_PARAMS, SIGMA
import sharptab.profile as profile
import sharptab.params as params
from sharptab import derived
from utils.timing import timeit

def filter(data):
    """
    Perform smoothing and data filtering

    Parameters:
    -----------
    data: list of dictionaries
        A list of dictionaries, with each entry corresponding to a model time step or
        forecast hour, containing the data to be filtered and smoothed

    Returns:
    --------
    data: list of dictionaries
        Filtered and smoothed data. Same form as input.
    """

    for t in range(len(data)):
        vars = list(data[t].keys())
        for v in vars:
            if v in SCALAR_PARAMS.keys():
                data[t][v] = gaussian_filter(data[t][v], sigma=SIGMA)

        # Filter the effective bulk wind difference (< 20 kts)
        idx = [i for i, s in enumerate(vars) if 'ebwd' in s]
        if len(idx) > 0:
            tmp = np.sqrt(data[t][vars[idx[0]]]**2 + data[t][vars[idx[0]]]**2)
            data[t][vars[idx[0]]] = np.where(tmp < 20, 0, data[t][vars[idx[0]]])
            data[t][vars[idx[1]]] = np.where(tmp < 20, 0, data[t][vars[idx[1]]])

        # Filter the 0-1 km shear (< 10 kts)
        idx = [i for i, s in enumerate(vars) if 'shr1' in s]
        if len(idx) > 0:
            tmp = np.sqrt(data[t][vars[idx[0]]]**2 + data[t][vars[idx[0]]]**2)
            data[t][vars[idx[0]]] = np.where(tmp < 10, 0, data[t][vars[idx[0]]])
            data[t][vars[idx[1]]] = np.where(tmp < 10, 0, data[t][vars[idx[1]]])

        # Filter the 0-3 km shear (< 20 kts)
        idx = [i for i, s in enumerate(vars) if 'shr3' in s]
        if len(idx) > 0:
            tmp = np.sqrt(data[t][vars[idx[0]]]**2 + data[t][vars[idx[0]]]**2)
            data[t][vars[idx[0]]] = np.where(tmp < 20, 0, data[t][vars[idx[0]]])
            data[t][vars[idx[1]]] = np.where(tmp < 20, 0, data[t][vars[idx[1]]])

    return data

float_array = types.float64[:,:] # No type expressions allowed in jitted functions
@timeit
@njit(parallel=True)
def worker(pres, tmpc, hght, dwpc, wspd, wdir, SCALARS, VECTORS):
    """
    While numba massively speeds up our computations, we're limited in how we store and
    process our data. A dictionary registry of parameter calculations doesn't work, so
    instead, embed specific if blocks (CASE?) to direct our workflow appropriately, while
    maintaining some degree of expandability to different parameters.

    Parameters:
    -----------
    pres: numpy array
        Array of pressure values (MB) [z,y,x]
    tmpc: numpy array
        Array of temperature values (C) [z,y,x]
    hght: numpy array
        Array of geopotential height values (m) [z,y,x]
    dwpc: numpy array
        Array of dewpoint values (C) [z,y,x]
    wspd: numpy array
        Array of wind speed values (KTS) [z,y,x]
    wdir: numpy array
        Array of wind direction values (DEG) [z,y,x]
    SCALARS: Numba typed List
        Scalar parameter keys, passed in from SCALAR_PARAMS in the config file
    VECTORS: Numba typed List
        Vector parameter keys, passed in from VECTOR_PARAMS in the config file

    Returns:
    --------
    d : Numba typed Dictionary
        Dictionary containing the derived values. [Key,j,i]

    """

    # Declare 'jit-able' dictionary and fill it with empty arrays
    d = Dict.empty(
        key_type=types.unicode_type,
        value_type=float_array,
    )
    for scalar in SCALARS:
        d[scalar] = np.zeros((tmpc.shape[1], tmpc.shape[2]), dtype='float64')
    for vector in VECTORS:
        d[vector+'_u'] = np.zeros((tmpc.shape[1], tmpc.shape[2]), dtype='float64')
        d[vector+'_v'] = np.zeros((tmpc.shape[1], tmpc.shape[2]), dtype='float64')

    for j in prange(tmpc.shape[1]):
        for i in prange(tmpc.shape[2]):
            prof = profile.create_profile(pres=pres[:,j,i], tmpc=tmpc[:,j,i],
                                          hght=hght[:,j,i], dwpc=dwpc[:,j,i],
                                          wspd=wspd[:,j,i], wdir=wdir[:,j,i])

            # Compile the big jitted methods. This will slow down the very 1st iteration.
            mlpcl = params.parcelx(prof, flag=4)
            eff_inflow = params.effective_inflow_layer(prof)
            mupcl = params.parcelx(prof, flag=3)

            # Scalars
            if 'mlcape' in SCALARS: d['mlcape'][j,i] = mlpcl.bplus
            if 'mlcin' in SCALARS: d['mlcin'][j,i] = mlpcl.bminus * -1
            if 'mucape' in SCALARS: d['mucape'][j,i] = mupcl.bplus
            if 'cape3km' in SCALARS: d['cape3km'][j,i] = mlpcl.b3km
            if 'esrh' in SCALARS: d['esrh'][j,i] = derived.esrh(prof, eff_inflow)

            # Vectors: returned as (u,v) tuples
            if 'ebwd' in VECTORS:
                d['ebwd_u'][j,i], d['ebwd_v'][j,i] = derived.ebwd(prof, mupcl, eff_inflow)
            if 'shr1' in VECTORS:
                d['shr1_u'][j,i], d['shr1_v'][j,i] = derived.bulk_shear(prof, height=1000)
            if 'shr3' in VECTORS:
                d['shr3_u'][j,i], d['shr3_v'][j,i] = derived.bulk_shear(prof, height=3000)
            if 'rm5' in VECTORS:
                d['rm5_u'][j,i], d['rm5_v'][j,i] = derived.rm5(prof)
            if 'lm5' in VECTORS:
                d['lm5_u'][j,i], d['lm5_v'][j,i] = derived.lm5(prof)
            if 'devtor' in VECTORS:
                d['devtor_u'][j,i], d['devtor_v'][j,i] = derived.devtor(prof)

            # Special parameters: prohibitive to re-compute all of the inputs...
            if 'estp' in SCALARS: d['estp'][j,i] = derived.estp(d['mlcape'][j,i],
                                                                d['mlcin'][j,i],
                                                                d['esrh'][j,i],
                                                                d['ebwd_u'][j,i],
                                                                d['ebwd_v'][j,i],
                                                                mlpcl)
    return d

def sharppy_calcs(**kwargs):
    """
    Perform the parcel lifting and calculations. Leverages numba's automatic parallel
    processing and prange.

    """

    tmpc = kwargs.get('tmpc')
    dwpc = kwargs.get('dwpc')
    hght = kwargs.get('hght')
    wdir = kwargs.get('wdir')
    wspd = kwargs.get('wspd')
    pres = kwargs.get('pres')

    # Convert list of dictionary keys to numba typed list.
    ret = worker(pres, tmpc, hght, dwpc, wspd, wdir, List(SCALAR_PARAMS.keys()),
                 List(VECTOR_PARAMS.keys()))

    # Converison back to a 'normal' Python dictionary
    output = {}
    for k, v in ret.items(): output[k] = v
    return output
