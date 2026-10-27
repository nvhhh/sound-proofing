"""Normal-incidence model of a micro-perforated panel (MPP) over a rigid back.

Maa (1998) microperforated-panel impedance, air gap D behind a rigid wall.
All quantities SI: metres, seconds, Pa.
"""
import numpy as np
from scipy.special import jv

C = 343.0          # speed of sound, m/s
RHO = 1.21         # air density, kg/m^3
MU = 1.81e-5       # dynamic viscosity of air, Pa*s
RHOC = RHO * C


def mpp_absorption(f, d, t, pitch, gap):
    """Absorption coefficient alpha(f) of one MPP layer.

    d      hole diameter (m)
    t      plate thickness (m)
    pitch  hole centre spacing on a square grid (m)
    gap    air cavity depth behind the plate (m)
    """
    f = np.atleast_1d(np.asarray(f, dtype=float))
    w = 2 * np.pi * f
    sigma = np.pi * d**2 / (4 * pitch**2)          # perforation ratio
    k = d * np.sqrt(w * RHO / (4 * MU))            # perforation shear number
    r_mpp = (32 * MU * t) / (sigma * d**2) * (np.sqrt(1 + k**2 / 32) + (np.sqrt(2) / 32) * k * d / t)
    x_mpp = (w * RHO * t) / sigma * (1 + 1 / np.sqrt(9 + k**2 / 2) + (np.sqrt(2) / 8) * k * d / t)
    k0 = w / C
    z_air = -1j * RHOC / np.tan(k0 * gap)           # rigid-backed air gap
    z = (r_mpp + 1j * x_mpp) + z_air
    refl = (z - RHOC) / (z + RHOC)
    return 1 - np.abs(refl) ** 2, refl


def band_mean(alpha, f, lo=300.0, hi=3000.0):
    m = (f >= lo) & (f <= hi)
    return float(np.mean(alpha[m]))


def _mpp_z_exact(w, d, t, pitch, t_end=0.0):
    """Exact viscous tube impedance per leaf (Bessel form), with an end correction.

    Resolves the viscous boundary layer in each hole, unlike the Maa (1998) fit.
    t_end adds an effective plate thickness for the open ends of each hole.
    """
    sigma = np.pi * d**2 / (4 * pitch**2)
    kap = (d / 2) * np.sqrt(-1j * w * RHO / MU)
    f = 1 - 2 * jv(1, kap) / (kap * jv(0, kap))
    return 1j * w * RHO * (t + t_end) / sigma / f


def _mpp_z(w, d, t, pitch):
    sigma = np.pi * d**2 / (4 * pitch**2)
    k = d * np.sqrt(w * RHO / (4 * MU))
    r = (32 * MU * t) / (sigma * d**2) * (np.sqrt(1 + k**2 / 32) + (np.sqrt(2) / 32) * k * d / t)
    x = (w * RHO * t) / sigma * (1 + 1 / np.sqrt(9 + k**2 / 2) + (np.sqrt(2) / 8) * k * d / t)
    return r + 1j * x


def stack_field(f, layers, open_frac=None, leaf='exact', t_end_factor=0.85):
    """Field and energy split through MPP leaves in series, each with its own air gap.

    layers     list of (d, t, pitch, gap) from the front face to the rigid back;
               gap is the cavity behind that leaf.
    leaf       'exact' (Bessel tube, default) or 'maa98' (Maa 1998 fit, no end correction)
    t_end_factor  end correction as a multiple of the hole diameter (0.85 is the usual value)
    open_frac  per-leaf fraction of holes that open into the cavity (None = all).
               Holes sitting on a rib open onto solid plastic, so they do not count.

    Returns a dict with, per frequency:
      alpha, refl  absorption and reflection coefficient at the front face
      share        absorbed fraction taken by each leaf (sums to alpha)
      cav          per cavity, (P, Q) at its top face, so that inside the cavity
                   p(x) = P exp(-j k x) + Q exp(+j k x), x measured down from the top
    """
    f = np.atleast_1d(np.asarray(f, dtype=float))
    w = 2 * np.pi * f
    k0 = w / C
    open_frac = [1.0] * len(layers) if open_frac is None else open_frac
    if leaf == 'exact':
        zs = [_mpp_z_exact(w, d, t, p, t_end_factor * d) / o for (d, t, p, _), o in zip(layers, open_frac)]
    else:
        zs = [_mpp_z(w, d, t, p) / o for (d, t, p, _), o in zip(layers, open_frac)]
    gaps = [g for _, _, _, g in layers]

    # Impedance from the rigid back toward the front face.
    z_load = zs[-1] + (-1j * RHOC / np.tan(k0 * gaps[-1]))
    for z, gap in zip(reversed(zs[:-1]), reversed(gaps[:-1])):
        tn = np.tan(k0 * gap)
        z_load = z + RHOC * (z_load + 1j * RHOC * tn) / (RHOC + 1j * z_load * tn)
    refl = (z_load - RHOC) / (z_load + RHOC)
    alpha = 1 - np.abs(refl) ** 2

    # Walk the pressure p and velocity u (positive into the panel) down the stack.
    p = 1 + refl
    u = (1 - refl) / RHOC
    share, cav = [], []
    for z, gap in zip(zs, gaps):
        share.append(z.real * np.abs(u) ** 2 * RHOC)
        p = p - z * u                                  # pressure drop across the leaf
        P = (p + RHOC * u) / 2
        Q = (p - RHOC * u) / 2
        cav.append((P, Q))
        e_in = np.exp(-1j * k0 * gap)
        e_out = np.exp(1j * k0 * gap)
        p, u = P * e_in + Q * e_out, (P * e_in - Q * e_out) / RHOC
    return {'alpha': alpha, 'refl': refl, 'share': np.array(share), 'cav': cav}


def double_mpp(f, layers, open_frac=None, **kw):
    """Two (or more) MPP leaves in series. Returns (alpha, reflection coefficient)."""
    out = stack_field(f, layers, open_frac, **kw)
    return out['alpha'], out['refl']
