"""Normal-incidence model of a micro-perforated panel (MPP) over a rigid back.

Maa (1998) microperforated-panel impedance, air gap D behind a rigid wall.
All quantities SI: metres, seconds, Pa.
"""
import numpy as np

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


def _mpp_z(w, d, t, pitch):
    sigma = np.pi * d**2 / (4 * pitch**2)
    k = d * np.sqrt(w * RHO / (4 * MU))
    r = (32 * MU * t) / (sigma * d**2) * (np.sqrt(1 + k**2 / 32) + (np.sqrt(2) / 32) * k * d / t)
    x = (w * RHO * t) / sigma * (1 + 1 / np.sqrt(9 + k**2 / 2) + (np.sqrt(2) / 8) * k * d / t)
    return r + 1j * x


def double_mpp(f, layers):
    """Two (or more) MPP leaves in series, each followed by its own air gap.

    layers: list of (d, t, pitch, gap) from the front face to the rigid back.
    Returns (alpha, reflection coefficient) per frequency.
    """
    f = np.atleast_1d(np.asarray(f, dtype=float))
    w = 2 * np.pi * f
    k0 = w / C
    # Start at the rigid back and work toward the front.
    d, t, p, gap = layers[-1]
    z_load = np.zeros_like(w, dtype=complex)
    z_load = _mpp_z(w, d, t, p) + (-1j * RHOC / np.tan(k0 * gap))
    for d, t, p, gap in reversed(layers[:-1]):
        tn = np.tan(k0 * gap)
        z_gap = RHOC * (z_load + 1j * RHOC * tn) / (RHOC + 1j * z_load * tn)
        z_load = _mpp_z(w, d, t, p) + z_gap
    refl = (z_load - RHOC) / (z_load + RHOC)
    return 1 - np.abs(refl) ** 2, refl
