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


def neck_z(w, r, length, cell, t_end_factor=0.68):
    """Specific impedance on one cell of a single neck (tube of radius r, length), with viscous loss.

    The neck has area Sn = pi r^2 in a cell of area Sc = cell^2; sigma = Sn / Sc. Gives the
    Helmholtz mass term and the viscous resistance of the neck, per unit cell area.
    """
    sigma = np.pi * r**2 / cell**2
    d = 2 * r
    kap = (d / 2) * np.sqrt(-1j * w * RHO / MU)
    F = 1 - 2 * jv(1, kap) / (kap * jv(0, kap))
    return 1j * w * RHO * (length + t_end_factor * d) / sigma / F


def slit_branch_z(w, h, span, length, cell, D, pr=0.71, gamma=1.4):
    """Specific impedance on one cell of a slot (gap h, span) into a closed cavity of depth D.

    Wide-slit thermoviscous duct (Stinson / Allard form, the approach used in the
    engineering report): the viscous and thermal penetration functions give the effective
    density and bulk modulus of the slot air. Edge impedance uses the isolated-slit end
    correction. Slot air area h*span is a fraction phi of the cell area.
    """
    Fv_of = lambda x: np.tanh(x) / x
    def fv_ft(x):
        q = (x / 2) * np.sqrt(1j * w * RHO / MU)
        return Fv_of(q), np.tanh(q * np.sqrt(pr)) / (q * np.sqrt(pr))
    Fvh, Fth = fv_ft(h)
    Fvb, Ftb = fv_ft(span)
    Fv = Fvh + Fvb - Fvh * Fvb
    Ft = Fth + Ftb - Fth * Ftb
    rho_eff = RHO / (1 - Fv)
    K = RHO * C**2 / (1 + (gamma - 1) * Ft)
    Zc = np.sqrt(rho_eff * K)
    keff = w * np.sqrt(rho_eff / K)
    phi = h * span / cell**2
    k0 = w / C
    Zcav = -1j * RHOC / np.tan(k0 * D)
    ZL = phi * Zcav
    t = np.tan(keff * length)
    Zin = Zc * (ZL + 1j * Zc * t) / (Zc + 1j * ZL * t)
    delta = (2 * h / np.pi) * np.log(2 * span / h + 0.5)
    z_edge = 0.5 * np.sqrt(2 * MU * RHO * w) + 1j * w * RHO * delta
    return Zin / phi + z_edge


def slit_absorption(f, slits, D, cell):
    """Absorption of a face of slot resonators, one type per equal share of the face.

    slits: list of (h, span, length, share) with h, span, length, cell in metres;
    shares are fractions of the total face area (the rest is solid, admittance zero).
    """
    f = np.atleast_1d(np.asarray(f, dtype=float))
    w = 2 * np.pi * f
    Y = sum(share / slit_branch_z(w, h, span, L, cell, D) for h, span, L, share in slits)
    Z = 1 / Y
    R = (Z - RHOC) / (Z + RHOC)
    return 1 - np.abs(R) ** 2


def _line_load(Z_T, k0, gap):
    """Impedance at the top of a cavity of depth gap, terminated by Z_T at its bottom."""
    tn = np.tan(k0 * gap)
    return RHOC * (Z_T + 1j * RHOC * tn) / (RHOC + 1j * Z_T * tn)


def stack_field(f, layers, open_frac=None, leaf='exact', t_end_factor=0.68, neck=None, plate_mass=None):
    """Field and energy split through MPP leaves in series, each with its own air gap.

    layers     list of (d, t, pitch, gap) from the front face to the rigid back;
               gap is the cavity behind that leaf.
    leaf       'exact' (Bessel tube, default) or 'maa98' (Maa 1998 fit, no end correction)
    neck       None (rigid back) or dict(r, length, cell) for a Helmholtz neck in the back wall
    t_end_factor  end correction as a multiple of the hole diameter. 0.68 is calibrated to a
               2-D resolved Stokes simulation of one hole (0.66-0.70 d); the common 0.85 d overestimates it
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
    if plate_mass is not None:       # areal mass (kg/m^2) of each unsupported leaf; adds w*m to its impedance
        zs = [z + 1j * w * m for z, m in zip(zs, plate_mass)]

    # Impedance from the back toward the front face. The rear cavity ends in a rigid wall,
    # or in a neck (Helmholtz) if neck=dict(r=, length=, cell=) is given.
    if neck is None:
        z_load = zs[-1] + (-1j * RHOC / np.tan(k0 * gaps[-1]))
    else:
        z_load = zs[-1] + _line_load(neck_z(w, neck['r'], neck['length'], neck['cell'], t_end_factor), k0, gaps[-1])
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
    # Power leaving through a neck goes into the wall space behind the tile; it is counted
    # as absorbed (the wall space is assumed to take it up, not return it to the room).
    neck_share = np.zeros_like(alpha)
    if neck is not None:
        Zn = neck_z(w, neck['r'], neck['length'], neck['cell'], t_end_factor)
        neck_share = Zn.real * np.abs(u) ** 2 * RHOC
    return {'alpha': alpha, 'refl': refl, 'share': np.array(share), 'cav': cav, 'neck_share': neck_share}


def double_mpp(f, layers, open_frac=None, **kw):
    """Two (or more) MPP leaves in series. Returns (alpha, reflection coefficient)."""
    out = stack_field(f, layers, open_frac, **kw)
    return out['alpha'], out['refl']
