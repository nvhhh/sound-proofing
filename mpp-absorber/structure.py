"""Plate-vibration check: can the leaves or rib walls resonate in 300-3000 Hz and re-radiate?

Thin-plate theory, simply supported edges (ribs and walls are the supports):
    f_mn = (pi/2) * sqrt(D / (rho h)) * ((m/a)^2 + (n/b)^2),  D = E h^3 / (12 (1 - nu^2))

Mass-law estimate of the sound a leaf re-radiates, as a fraction of the incident
intensity, for a leaf driven at pressure p with areal mass m'. Above the lowest
mode the leaf is mass-controlled and its velocity is p / (w m').

Material defaults are typical PLA: E 3.5 GPa, nu 0.36, rho 1240 kg/m^3.
"""
import numpy as np

from model import RHOC

E_PLA, NU_PLA, RHO_PLA = 3.5e9, 0.36, 1240.0

# Typical printing-grade values (datasheet ranges, rounded). Stiffness E, Poisson nu, density rho.
MATERIALS = {
    'PLA':   dict(E=3.5e9, nu=0.36, rho=1240.0),
    'PETG':  dict(E=2.1e9, nu=0.38, rho=1270.0),
    'ABS':   dict(E=2.3e9, nu=0.35, rho=1040.0),
    'Nylon': dict(E=1.5e9, nu=0.40, rho=1140.0),   # PA6/PA12, dry; stiffens less when wet
}


def plate_modes(a, b, h, E=E_PLA, nu=NU_PLA, rho=RHO_PLA, f_max=10000.0, m_max=12, n_max=12):
    D = E * h**3 / (12 * (1 - nu**2))
    c = np.sqrt(D / (rho * h))
    f = []
    for m in range(1, m_max + 1):
        for n in range(1, n_max + 1):
            fm = np.pi / 2 * c * ((m / a) ** 2 + (n / b) ** 2)
            if fm <= f_max:
                f.append((fm, m, n))
    return sorted(f)


def mass_law_radiation(f, areal_mass):
    """Radiated / incident intensity for a mass-controlled leaf (normal incidence)."""
    w = 2 * np.pi * np.asarray(f, dtype=float)
    return (RHOC / (w * areal_mass)) ** 2


def plate_velocity_ratio(f, a, h, E=E_PLA, nu=NU_PLA, rho=RHO_PLA):
    """Velocity of a clamped square leaf cell relative to the acoustic particle velocity.

    Uses the clamped-plate centre deflection w0 = 0.00126 p a^4 / D per Pa of pressure.
    A ratio well below 1 means the leaf stays still and the air moves through the holes.
    """
    D = E * h**3 / (12 * (1 - nu**2))
    w = 2 * np.pi * np.asarray(f, dtype=float)
    v_plate = w * 0.00126 * a**4 / D
    v_air = 1 / RHOC
    return v_plate / v_air


def material_table(a=0.018, h=0.8e-3):
    rows = {}
    f = np.geomspace(300, 3000, 240)
    for name, m in MATERIALS.items():
        modes = plate_modes(a, a, h, E=m['E'], nu=m['nu'], rho=m['rho'])
        ratio = plate_velocity_ratio(f, a, h, E=m['E'], nu=m['nu'], rho=m['rho'])
        rows[name] = dict(first_mode_Hz=round(modes[0][0]), max_velocity_ratio=round(float(ratio.max()), 4))
    return rows


if __name__ == '__main__':
    h = 0.8e-3
    # Leaves span one rib cell (18 mm), since the ribs support them.
    cell_modes = plate_modes(0.018, 0.018, h)
    print('leaf panel, 18 mm cell, lowest modes (Hz):', [round(x[0]) for x in cell_modes[:4]])
    # The same leaf unsupported across the whole 108 mm tile, for comparison.
    full = plate_modes(0.108, 0.108, h)
    print('leaf panel, 108 mm unsupported, lowest modes (Hz):', [round(x[0]) for x in full[:4]])
    rib = plate_modes(0.025, 0.018, 1.2e-3)
    print('rib wall 25 x 18 mm, 1.2 mm, lowest modes (Hz):', [round(x[0]) for x in rib[:4]])
    f = np.geomspace(300, 3000, 240)
    r = plate_velocity_ratio(f, 0.018, h)
    print('leaf velocity / air velocity, 300-3000 Hz: max', round(float(r.max()), 4))
    for name, row in material_table().items():
        print(name, row)
    m_leaf = RHO_PLA * h
    f = np.array([300.0, 1000.0, 3000.0])
    print('mass-law re-radiation of a leaf (fraction of incident intensity):',
          dict(zip(f.astype(int), np.round(mass_law_radiation(f, m_leaf), 4))))

