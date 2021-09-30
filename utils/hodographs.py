"""Code taken from Tim Supinie's vad-plotter GitHub repo"""
import numpy as np
from sharptab.interp import generic_interp_hght
from sharptab.winds import comp2vec, vec2comp

def parse_vector(vec_str):
    return tuple(int(v) for v in vec_str.strip().split("/"))

def _clip_profile(prof, alt, clip_alt, intrp_prof):
    try:
        idx_clip = np.where((alt[:-1] <= clip_alt) & (alt[1:] > clip_alt))[0][0]
    except IndexError:
        return np.nan * np.ones(prof.size)

    prof_clip = prof[:(idx_clip + 1)]
    prof_clip = np.append(prof_clip, intrp_prof)

    return np.array(prof_clip)

def compute_bunkers(data):
    d = 7.5 * 1.94     # Deviation value emperically derived as 7.5 m/s
    hght = 6

    # SFC-6km Mean Wind
    u, v = data['uwnd'], data['vwnd']
    u_hght = generic_interp_hght(hght, data['hght'], u)
    v_hght = generic_interp_hght(hght, data['hght'], v)
    u_clip = _clip_profile(u, data['hght'], hght, u_hght)
    v_clip = _clip_profile(v, data['hght'], hght, v_hght)

    mnu6 = u_clip.mean()
    mnv6 = v_clip.mean()

    # SFC-6km Shear Vector
    shru = u_hght - u[0]
    shrv = v_hght - v[0]

    # Bunkers Right Motion
    tmp = d / np.hypot(shru, shrv)
    rstu = mnu6 + (tmp * shrv)
    rstv = mnv6 - (tmp * shru)
    lstu = mnu6 - (tmp * shrv)
    lstv = mnv6 + (tmp * shru)

    return comp2vec(rstu, rstv), comp2vec(lstu, lstv), comp2vec(mnu6, mnv6)

def compute_shear_mag(data, hght1, hght2=None):
    u, v = data['uwnd'], data['vwnd']
    u1 = generic_interp_hght(hght1, data['hght'], u)
    v1 = generic_interp_hght(hght1, data['hght'], v)
    V1 = u1, v1
    if hght2 is not None:
        u2 = generic_interp_hght(hght2, data['hght'], u)
        v2 = generic_interp_hght(hght2, data['hght'], v)
        V2 = u2, v2
    else:
        V2 = u[0], v[0]
    return np.hypot(V1[0] - V2[0], V1[1] - V2[1])

def compute_srh(data, storm_motion, hght):
    u, v = data['uwnd'], data['vwnd']
    if len(u) < 2 and len(v) < 2:
        return np.nan

    storm_u, storm_v = vec2comp(*storm_motion)

    sru = (u - storm_u) / 1.94
    srv = (v - storm_v) / 1.94

    sru_hght = generic_interp_hght(hght, data['hght'], sru)
    srv_hght = generic_interp_hght(hght, data['hght'], srv)
    sru_clip = _clip_profile(sru, data['hght'], hght, sru_hght)
    srv_clip = _clip_profile(srv, data['hght'], hght, srv_hght)

    layers = (sru_clip[1:] * srv_clip[:-1]) - (sru_clip[:-1] * srv_clip[1:])
    return layers.sum()

def compute_crit_angl(data, storm_motion):
    """Modified to compute a 0-3-km weighted critical angle?
    """

    u, v = data['uwnd'], data['vwnd']
    storm_u, storm_v = vec2comp(*storm_motion)

    hght_base = 0
    heights = np.arange(.25, 1.75, .25)
    angles = []
    for hght in heights:
        u_base = generic_interp_hght(hght_base, data['hght'], u)
        v_base = generic_interp_hght(hght_base, data['hght'], v)

        u_intrp = generic_interp_hght(hght, data['hght'], u)
        v_intrp = generic_interp_hght(hght, data['hght'], v)

        base_u = storm_u - u_base
        base_v = storm_v - v_base

        ang_u = u_intrp - u_base
        ang_v = v_intrp - v_base

        len_base = np.hypot(base_u, base_v)
        len_ang = np.hypot(ang_u, ang_v)

        base_dot_ang = base_u * ang_u + base_v * ang_v
        angle = np.degrees(np.arccos(base_dot_ang / (len_base * len_ang)))
        angles.append(angle)
        hght_base = hght

    mean_angle = np.nanmean(angles)
    return round(mean_angle,1)

def compute_parameters(data, storm_motion):
    params = {}

    try:
        params['bunkers_right'], params['bunkers_left'], params['mean_wind'] = compute_bunkers(data)
    except (IndexError, ValueError):
        params['bunkers_right'] = (np.nan, np.nan)
        params['bunkers_left'] = (np.nan, np.nan)
        params['mean_wind'] = (np.nan, np.nan)

    if storm_motion.lower() in ['blm', 'left-mover']:
        params['storm_motion'] = params['bunkers_left']
    elif storm_motion.lower() in ['brm', 'right-mover']:
        params['storm_motion'] = params['bunkers_right']
    elif storm_motion.lower() in ['mnw', 'mean-wind']:
        params['storm_motion'] = params['mean_wind']
    else:
        params['storm_motion'] = tuple(int(v) for v in storm_motion.split('/'))


    try:
        params['critical'] = compute_crit_angl(data, params['storm_motion'])
    except (IndexError, ValueError):
        params['critical'] = np.nan

    for hght in [1, 2, 3, 6]:
        try:
            params["shear_mag_%dkm" % hght] = compute_shear_mag(data, hght)
        except (IndexError, ValueError):
            params["shear_mag_%dkm" % hght] = np.nan

    # Mid-level shear
    params["shear_mag_16km"] = compute_shear_mag(data, 1, 6)

    for hght in [1, 2, 3]:
        params["srh_%dkm" % hght] = compute_srh(data, params['storm_motion'], hght)

    return params
