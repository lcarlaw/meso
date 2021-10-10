"""Functions used during the parameter calculation steps. See calc.compute and the
`worker` function for specifics.
"""

from numba import njit
import numpy as np

from sharptab.constants import MS2KTS, KTS2MS
import sharptab.interp as interp
import sharptab.winds as winds
import sharptab.utils as utils
import sharptab.params as params
from calc.vector import transform

@njit
def lapse_rate(prof, lower=0, upper=3000):
    """
    Compute the lapse rate between two layers.

    Parameters:
    -----------
    prof: SHARPpy Profile object
    lower: number (optional; default 0)
        Bottom layer (m) for lapse rate calculation.
    upper: number (optional; default 0)
        Top layer (m) for lapse rate calculation.

    Returns:
    --------
    gamma: number
        Lapse rate compute between top and bottom layers (C/km)

    """
    z1 = interp.to_msl(prof, lower)
    z2 = interp.to_msl(prof, upper)
    p1 = interp.pres(prof, z1)
    p2 = interp.pres(prof, z2)
    tv1 = interp.vtmp(prof, p1)
    tv2 = interp.vtmp(prof, p2)
    return (tv2 - tv1) / (z2 - z1) * -1000.

@njit
def srh(prof, lower=None, upper=None, effective_inflow_layer=None):
    """
    Compute storm-relative helicity.

    Parameters:
    -----------
    prof: SHARPpy Profile object
    lower: number (optional; default None)
        Bottom layer (m) for SRH calculation
    upper: number (optional; default None)
        Top layer (m) for SRH calculation
    effective_inflow_layer: array_like, tuple, or list (optional, default None)
        Bottom and Top (in mb) of the effective inflow layer

    Returns:
    --------
    srh: number
        Storm-relative helicity value (m2/s2)

    """
    RM5 = rm5(prof)
    if effective_inflow_layer is not None:
        ebot = interp.to_agl(prof, interp.hght(prof, effective_inflow_layer[0]))
        etop = interp.to_agl(prof, interp.hght(prof, effective_inflow_layer[1]))
        srh = winds.helicity(prof, ebot, etop, stu=RM5[0], stv=RM5[1])[0]
    elif lower is not None and upper is not None:
        srh = winds.helicity(prof, lower, upper, stu=RM5[0], stv=RM5[1])[0]
    return srh

@njit
def estp(mlcape, mlcin, esrh, ebwd_u, ebwd_v, mlpcl, eff_inflow_base, prof):
    eshr = utils.mag(ebwd_u, ebwd_v)
    estp = params.stp_cin(mlcape, esrh, eshr*KTS2MS, mlpcl.lclhght, mlcin)

    ebot = interp.to_agl(prof, interp.hght(prof, eff_inflow_base))
    if ebot > 10 or ~np.isfinite(ebot): estp = 0.
    return estp

@njit
def devtor(prof):
    """
    Deviant tornado motion following Cameron Nixon's work:
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
