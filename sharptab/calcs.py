from scipy.ndimage import gaussian_filter
import numpy as np
from numba import njit, prange
from functools import reduce
#import ray

from sharptab.constants import *
import sharptab.profile as profile
import sharptab.params as params
import sharptab.interp as interp
import sharptab.winds as winds
import sharptab.utils as utils

from time import time
import logging as log

#ray.shutdown()
#ray.init(address='auto', _redis_password='5241590000000000')

sigma = 1.75
def componentsTo(u,v):
    return (u, v)

def add(*args):
    """From AWIPS. Perform scalar or vector addition
    """
    def scalarAddition(args):
        return reduce(np.add, args)

    def vectorAddition(args):
        uResult = np.zeros_like(args[0][0])
        vResult = np.zeros_like(args[0][0])
        for u, v in args:
            uResult += u
            vResult += v
        return componentsTo(uResult, vResult)

    if len(args)==1 and isinstance(args[0], list):
        return add(*args[0])
    elif isinstance(args[0], tuple):
        return vectorAddition(args)
    else:
        return scalarAddition(args)

def multiply(*args):
    """From AWIPS. Perform multiplication of any number of scalars or of a vector and a
    scalar.

    """
    def scalarMultiply(args):
        return reduce(multiply, args)

    def vectorMultiply(args):
        return componentsTo(scalarMultiply((args[0][0],  args[1])),
                            scalarMultiply((args[0][1],  args[1])))

    if type(args[0]) == tuple:
        return vectorMultiply(args)
    else:
        return scalarMultiply(args)

def divide(*args):
    """From AWIPS. Divide a scalar by a scalar or a vector by a scalar.
    """
    divArgs = list(args)

    for i in range(1,len(divArgs)):
        divArgs[i] = np.where(divArgs[i] == 0, np.float32(np.nan), 1/divArgs[i])
    return multiply(divArgs)

@njit
def transform(uComponent, vComponent, s1, s2, s3, s4):
    """From AWIPS. Rotate a vector

    Rotate vector by a number of degrees, or transform the vector with a matrix.
    The arguments after the vector need to be constants.
    """
    #uComponentShape = uComponent.shape
    #vComponentShape = vComponent.shape

    # in-place flatten the arrays
    #uComponent = uComponent.reshape((uComponent.size,))
    #vComponent = vComponent.reshape((vComponent.size,))

    vector = np.array([uComponent, vComponent])
    transform = np.array([[s1, s2], [s3, s4]])
    u, v = np.dot(transform, vector)

    # resize the arrays back to appropriate size
    #u.resize(uComponentShape)
    #v.resize(vComponentShape)
    return (u, v)

@njit(parallel=True)
def worker(pres, tmpc, hght, dwpc, wspd, wdir):
    mucape = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    mlcape = np.zeros_like(mucape)
    cape3km = np.zeros_like(mucape)
    lr3km = np.zeros_like(mucape)
    mlcin = np.zeros_like(mucape)
    mllcl = np.zeros_like(mucape)
    ebot = np.zeros_like(mucape)
    etop = np.zeros_like(mucape)
    eshr = np.zeros_like(mucape)
    esrh = np.zeros_like(mucape)
    estp = np.zeros_like(mucape)
    tts = np.zeros_like(mucape)

    rm5_u = np.zeros_like(mucape)
    rm5_v = np.zeros_like(mucape)
    lm5_u = np.zeros_like(mucape)
    lm5_v = np.zeros_like(mucape)
    ebwd_u = np.zeros_like(mucape)
    ebwd_v = np.zeros_like(mucape)
    shr1_u = np.zeros_like(mucape)
    shr1_v = np.zeros_like(mucape)
    shr3_u = np.zeros_like(mucape)
    shr3_v = np.zeros_like(mucape)

    for j in prange(tmpc.shape[1]):
        for i in prange(tmpc.shape[2]):
            prof = profile.create_profile(pres=pres[:,j,i], tmpc=tmpc[:,j,i],
                                          hght=hght[:,j,i], dwpc=dwpc[:,j,i],
                                          wspd=wspd[:,j,i], wdir=wdir[:,j,i])

            sfc = prof.pres[prof.sfc]
            p6km = interp.pres(prof, interp.to_msl(prof, 6000.))
            blkshr_u, blkshr_v = winds.wind_shear(prof, pbot=sfc, ptop=p6km)
            mean_u, mean_v = winds.mean_wind(prof, pbot=sfc,ptop=p6km)

            # Bunkers Right and Left motion vectors (approximations here...)
            temp = transform(blkshr_u, blkshr_v, 0., 7.5*MS2KTS, -7.5*MS2KTS, 0.)
            BlkMag = np.hypot(blkshr_u, blkshr_v)
            rm5_u[j,i] = mean_u + (temp[0]/BlkMag)
            rm5_v[j,i]  = mean_v  + (temp[1]/BlkMag)
            lm5_u[j,i]  = mean_u  - (temp[0]/BlkMag)
            lm5_v[j,i]  = mean_v  - (temp[1]/BlkMag)

            eff_inflow = params.effective_inflow_layer(prof)
            ebot[j,i] = interp.to_agl(prof, interp.hght(prof, eff_inflow[0]))
            etop[j,i] = interp.to_agl(prof, interp.hght(prof, eff_inflow[1]))

            # Parcel buoyancy calculations
            mupcl = params.parcelx(prof, flag=3)
            mlpcl = params.parcelx(prof, flag=4)

            mllcl[j,i] = mlpcl.lclhght
            mlcape[j,i] = mlpcl.bplus
            mlcin[j,i] = mlpcl.bminus
            mucape[j,i] = mupcl.bplus
            cape3km[j,i] = mlpcl.b3km

            # Effective BWD
            height_bot = interp.pres(prof, mupcl.pres)
            height_top = (mupcl.elhght + height_bot) / 2.
            ptop = interp.pres(prof, interp.to_msl(prof, height_top))
            ebwd_u[j,i], ebwd_v[j,i] = winds.wind_shear(prof, pbot=eff_inflow[0], ptop=ptop)
            eshr[j,i] = utils.mag(ebwd_u[j,i], ebwd_v[j,i])

            p1km = interp.pres(prof, interp.to_msl(prof, 1000.))
            shr1_u[j,i], shr1_v[j,i] = winds.wind_shear(prof, pbot=sfc, ptop=p1km)

            p3km = interp.pres(prof, interp.to_msl(prof, 3000.))
            shr3_u[j,i], shr3_v[j,i] = winds.wind_shear(prof, pbot=sfc, ptop=p3km)
            t3km = interp.temp(prof, p3km)
            lr3km[j,i] = (tmpc[0,j,i] - t3km) / 3.

            esrh[j,i] = winds.helicity(prof, ebot[j,i], etop[j,i], stu=rm5_u[j,i],
                                       stv=rm5_v[j,i])[0]
            estp[j,i] = params.stp_cin(mlcape[j,i], esrh[j,i], eshr[j,i],
                                       mllcl[j,i], mlcin[j,i])
            p1km = interp.pres(prof, interp.to_msl(prof, 1000.))
            srh1 = winds.helicity(prof, ebot[j,i], p1km, stu=rm5_u[j,i],
                                       stv=rm5_v[j,i])[0]


            # Tornadic Tilting and Stretching parameter (TTS)
            A = (srh1 * np.clip(cape3km[j,i], 0, 150)) / 6500.
            B = np.clip(mlcape[j,i] / 2000., 1, 1.5)
            ebwd = np.hypot(ebwd_u[j,i], ebwd_v[j,i]) * KTS2MS
            ebwd = np.where(ebwd < 12.5, 0, ebwd)
            C = np.clip(ebwd / 20., 0, 1.5)
            TTS = np.clip(A*B*C, 0, 9999)
            if mllcl[j,i] > 1700 or mlcin[j,i] < -100: TTS = 0
            tts[j,i] = TTS

    return esrh, estp, cape3km, lr3km, mlcape, mlcin, mucape, tts, ebwd_u, ebwd_v, \
           rm5_u, rm5_v, lm5_u, lm5_v, shr1_u, shr1_v, shr3_u, shr3_v

def sharppy_calcs(**kwargs):
    tmpc = kwargs.get('tmpc')
    dwpc = kwargs.get('dwpc')
    hght = kwargs.get('hght')
    wdir = kwargs.get('wdir')
    wspd = kwargs.get('wspd')
    pres = kwargs.get('pres')

    mucape = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    vectors = {
        'ebwd_u': np.zeros_like(mucape),
        'ebwd_v': np.zeros_like(mucape),
        'shr3_u': np.zeros_like(mucape),
        'shr3_v': np.zeros_like(mucape),
        'shr1_u': np.zeros_like(mucape),
        'shr1_v': np.zeros_like(mucape),
        'rm5_u': np.zeros_like(mucape),
        'rm5_v': np.zeros_like(mucape),
        'lm5_u': np.zeros_like(mucape),
        'lm5_v': np.zeros_like(mucape)
    }

    ret = worker(pres, tmpc, hght, dwpc, wspd, wdir)
    esrh, estp, cape3km, lr3km, mlcape, mlcin, mucape, tts, vectors['ebwd_u'],           \
    vectors['ebwd_v'], vectors['rm5_u'], vectors['rm5_v'], vectors['lm5_u'],             \
    vectors['lm5_v'], vectors['shr1_u'], vectors['shr1_v'], vectors['shr3_u'],           \
    vectors['shr3_v'] = ret

    # Apply some data masks
    lr3km = gaussian_filter(lr3km, sigma=sigma)
    mlcin = gaussian_filter(np.where(mlcape>20, mlcin, 0), sigma=sigma)
    mucape = gaussian_filter(mucape, sigma=sigma)
    mlcape = gaussian_filter(mlcape, sigma=sigma)
    estp = np.where(np.isnan(estp), 0, estp)
    estp = gaussian_filter(estp, sigma=sigma)
    cape3km = gaussian_filter(cape3km, sigma=sigma)
    esrh = gaussian_filter(np.where(np.isnan(esrh), 0, esrh), sigma=sigma)

    vectors['shr1_u'] = gaussian_filter(vectors['shr1_u'], sigma=sigma)
    vectors['shr1_v'] = gaussian_filter(vectors['shr1_v'], sigma=sigma)
    vectors['shr3_u'] = gaussian_filter(vectors['shr3_u'], sigma=sigma)
    vectors['shr3_v'] = gaussian_filter(vectors['shr3_v'], sigma=sigma)
    vectors['rm5_u'] = gaussian_filter(vectors['rm5_u'], sigma=sigma)
    vectors['rm5_v'] = gaussian_filter(vectors['rm5_v'], sigma=sigma)
    vectors['lm5_u'] = gaussian_filter(vectors['lm5_u'], sigma=sigma)
    vectors['lm5_v'] = gaussian_filter(vectors['lm5_v'], sigma=sigma)
    vectors['ebwd_u'] = gaussian_filter(vectors['ebwd_u'], sigma=sigma)
    vectors['ebwd_v'] = gaussian_filter(vectors['ebwd_v'], sigma=sigma)
    vectors['ebwd_u'] = np.where(np.isnan(vectors['ebwd_u']), 0, vectors['ebwd_u'])
    vectors['ebwd_v'] = np.where(np.isnan(vectors['ebwd_v']), 0, vectors['ebwd_v'])

    tts = gaussian_filter(np.where(np.isnan(tts), 0, tts), sigma=sigma)

    ret = {
        'esrh': esrh,
        'estp': estp,
        'cape3km': cape3km,
        'lr3km_cf': lr3km,
        'lr3km': lr3km,
        'mlcape': mlcape,
        'mlcin': mlcin * -1,
        'mucape': mucape,
        'vectors': vectors,
        'tts': tts
    }

    return ret
