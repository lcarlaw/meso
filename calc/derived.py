"""Functions used during the parameter calculation steps. See calc.compute.
"""

from numba import njit
import numpy as np

from sharptab.constants import MS2KTS
import sharptab.interp as interp
import sharptab.winds as winds
import sharptab.utils as utils
import sharptab.params as params
from calc.vector import transform

@njit
def srh500(prof):
    RM5 = rm5(prof)
    srh = winds.helicity(prof, 0, 500, stu=RM5[0], stv=RM5[1])[0]
    return srh

@njit
def estp(mlcape, mlcin, esrh, ebwd_u, ebwd_v, mlpcl):
    eshr = utils.mag(ebwd_u, ebwd_v)
    estp = params.stp_cin(mlcape, esrh, eshr, mlpcl.lclhght, mlcin)
    return estp

@njit
def devtor(prof):
    """Deviant tornado motion following Cameron Nixon's work:
    https://cameronnixonphotography.wordpress.com/research/anticipating-deviant-tornado-motion/
    """
    sfc = prof.pres[prof.sfc]
    A = winds.mean_wind(prof, pbot=sfc, ptop=interp.pres(prof, interp.to_msl(prof, 500)))
    B = rm5(prof)
    u = 0.5 * (A[0] + B[0])
    v = 0.5 * (A[1] + B[1])
    return (u, v)

@njit
def lm5(prof):
    sfc = prof.pres[prof.sfc]
    A = winds.mean_wind(prof, pbot=sfc, ptop=interp.pres(prof, interp.to_msl(prof, 500.)))
    B = winds.mean_wind(prof, pbot=interp.pres(prof, interp.to_msl(prof, 5500.)),
                        ptop=interp.pres(prof, interp.to_msl(prof, 6000.)))
    blkshr_u = B[0] - A[0]
    blkshr_v = B[1] - A[1]

    p6km = interp.pres(prof, interp.to_msl(prof, 6000.))
    mean_u, mean_v = winds.mean_wind(prof, pbot=sfc,ptop=p6km)

    temp = transform(blkshr_u, blkshr_v, 0., 7.5*MS2KTS, -7.5*MS2KTS, 0.)
    BlkMag = np.hypot(blkshr_u, blkshr_v)

    u = mean_u - (temp[0]/BlkMag)
    v = mean_v  - (temp[1]/BlkMag)
    return (u, v)

@njit
def rm5(prof):
    sfc = prof.pres[prof.sfc]
    A = winds.mean_wind(prof, pbot=sfc, ptop=interp.pres(prof, interp.to_msl(prof, 500.)))
    B = winds.mean_wind(prof, pbot=interp.pres(prof, interp.to_msl(prof, 5500.)),
                        ptop=interp.pres(prof, interp.to_msl(prof, 6000.)))
    blkshr_u = B[0] - A[0]
    blkshr_v = B[1] - A[1]

    p6km = interp.pres(prof, interp.to_msl(prof, 6000.))
    mean_u, mean_v = winds.mean_wind(prof, pbot=sfc,ptop=p6km)

    temp = transform(blkshr_u, blkshr_v, 0., 7.5*MS2KTS, -7.5*MS2KTS, 0.)
    BlkMag = np.hypot(blkshr_u, blkshr_v)

    u = mean_u + (temp[0]/BlkMag)
    v = mean_v  + (temp[1]/BlkMag)
    return (u, v)

@njit
def esrh(prof, eff_inflow):
    ebot = interp.to_agl(prof, interp.hght(prof, eff_inflow[0]))
    etop = interp.to_agl(prof, interp.hght(prof, eff_inflow[1]))

    RM5 = rm5(prof)
    esrh = winds.helicity(prof, ebot, etop, stu=RM5[0], stv=RM5[1])[0]
    return esrh

@njit
def bulk_shear(prof, height=1000):
    sfc = prof.pres[prof.sfc]
    hght = interp.pres(prof, interp.to_msl(prof, height))
    u, v = winds.wind_shear(prof, pbot=sfc, ptop=hght)
    return (u, v)

@njit
def ebwd(prof, mupcl, eff_inflow):
    ebot_hght = interp.to_agl(prof, interp.hght(prof, eff_inflow[0]))
    height_top = (mupcl.elhght + ebot_hght) / 2.
    ptop = interp.pres(prof, interp.to_msl(prof, height_top))
    u, v = winds.wind_shear(prof, pbot=eff_inflow[0], ptop=ptop)
    if ~np.isfinite(u): u = 0.
    if ~np.isfinite(v): v = 0.
    return (u, v)
