"""Hybrid 200 mm SCAD: geometry assertions re-checked in Python, then the acoustic model.

Face: 16 x 8 cells of 12 x 24 mm (192 mm active square). 40 resonator cells (4 tunings x 10)
and 88 channel cells (8 tunings x 11). Resonator: slot span 10.4 mm, gap h, neck length; a wide
cavity behind it (slit model, A). Channel: 14 closed ducts of 10.4 x 1.17 mm, length L (channel
model, A, rough). Face admittance: each cell type's share of the face, solid fraction zero.
"""
import numpy as np

from model import slit_absorption
from scad_check import channel_z, oct_mean

COLS, ROWS = 16, 8
CELL_W, CELL_H = 12.0, 24.0
SIDE, BACK, FRONT_RING = 0.8, 2.0, 8.0
TOTAL = 200.0
RES_GAP = [0.76, 1.04, 0.76, 2.72]
RES_NECK = [2.00, 20.76, 4.68, 53.85]
HF_GAP, FIN = 1.17, 0.45
HF_LEN = [164.42, 92.12, 73.34, 146.59, 108.79, 192.00, 192.00, 128.40]
SPAN = CELL_W - 2 * SIDE            # 10.4
WIDE = CELL_H - 2 * SIDE            # 22.4
RES_CELLS, HF_CELLS = 40, 88
CAV_AIR = TOTAL - 8                 # 192 mm, neck included
CELL_AREA = CELL_W * CELL_H         # 288 mm^2


def check_geometry():
    out = []
    out.append(('grid fills 192 mm', COLS * CELL_W == 192 and ROWS * CELL_H == 192))
    hf_rows = int(np.floor((WIDE + FIN) / (HF_GAP + FIN)))
    hf_stack = hf_rows * HF_GAP + (hf_rows - 1) * FIN
    hf_y0 = (CELL_H - hf_stack) / 2
    out.append(('HF slots per cell = %d' % hf_rows, hf_rows == 14))
    out.append(('HF slot stack %.2f mm fits 24 mm, margin %.2f' % (hf_stack, hf_y0), hf_y0 >= SIDE))
    out.append(('resonator cells = %d (31.25%%)' % RES_CELLS, abs(RES_CELLS / (COLS * ROWS) - 0.3125) < 1e-9))
    out.append(('HF cells = %d (4 * 11 per tuning count checks)' % HF_CELLS, HF_CELLS + RES_CELLS == 128))
    out.append(('max neck %.2f < air depth %.0f' % (max(RES_NECK), CAV_AIR), max(RES_NECK) < CAV_AIR))
    out.append(('max HF length %.0f + back %.0f <= %.0f' % (max(HF_LEN), BACK, TOTAL), max(HF_LEN) + BACK <= TOTAL))
    out.append(('smallest resonator gap %.2f mm (0.4 mm nozzle: >= 0.4 mm)' % min(RES_GAP), min(RES_GAP) >= 0.4))
    out.append(('HF slot 1.17 mm, fin 0.45 mm (line width 0.45 mm)', FIN >= 0.4))
    return out


def face_alpha(f, rng=None):
    jit = (lambda s: 0.0) if rng is None else (lambda s: rng.uniform(-s, s))
    w = 2 * np.pi * f
    Y = 0
    # resonators: each tuning has 10 of 128 cells
    for h, neck in zip(RES_GAP, RES_NECK):
        hh = h + jit(0.05); nn = neck + jit(0.5)
        Y = Y + (10 / 128) / (slit_absorption_branch(f, hh, nn, jit) )
    # channels: each tuning has 11 of 128 cells; mouth fraction 14 slots / cell
    phi = 14 * SPAN * HF_GAP / CELL_AREA
    for L in HF_LEN:
        Y = Y + (11 / 128) / (channel_z(w, (L + jit(0.5)) * 1e-3, (HF_GAP + jit(0.05)) * 1e-3, SPAN * 1e-3) / phi)
    Z = 1 / Y
    R = (Z - 412.972) / (Z + 412.972)
    return 1 - np.abs(R) ** 2


def slit_absorption_branch(f, h, neck, jit):
    """Branch impedance of one resonator (slot + wide cavity), per cell area."""
    from model import slit_branch_z
    w = 2 * np.pi * f
    return slit_branch_z(w, h * 1e-3, SPAN * 1e-3, neck * 1e-3, np.sqrt(CELL_AREA) * 1e-3, (CAV_AIR - neck) * 1e-3)


if __name__ == '__main__':
    print('Geometry (re-checked in Python, not OpenSCAD):')
    for msg, ok in check_geometry():
        print('  [%s] %s' % ('ok' if ok else 'FAIL', msg))
    f = np.geomspace(100, 6000, 500)
    a = face_alpha(f)
    print('\nAcoustic model (A, rough; no mouth end correction or coupling in channels):')
    print('  125-250 Hz linear mean %.2f' % a[(f >= 125) & (f <= 250)].mean())
    print('  500-4000 Hz linear mean %.2f, min %.2f' % (a[(f >= 500) & (f <= 4000)].mean(), a[(f >= 500) & (f <= 4000)].min()))
    for fc in (125, 250, 500, 1000, 2000, 4000):
        print('  octave %5d Hz  %.2f' % (fc, oct_mean(a, f, fc)))
    rng = np.random.default_rng(0)
    m = [face_alpha(f, rng)[(f >= 500) & (f <= 4000)].mean() for _ in range(40)]
    print('  robust 500-4000 Hz (40 draws): median %.2f, 5th pct %.2f' % (np.median(m), np.percentile(m, 5)))
