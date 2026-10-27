"""Shared hole and rib layout for the coupon and the tile.

Ribs are 1.2 mm walls on a cell grid; the first rib is at the left edge and a rib closes
the far edge, so each cell is a closed column. Holes sit midway between rib centrelines and
any hole closer than hole_r + rib_t/2 + WEB to a rib is dropped, so every kept hole has at
least WEB mm of wall to a rib. All lengths in mm.
"""
import numpy as np

WEB = 0.4          # minimum wall between a hole edge and a rib edge, mm


def rib_centres(width, cell, rib_t):
    """Rib centrelines along one axis: k*cell + rib_t/2, plus a closing rib at the far edge."""
    c = np.arange(0, width, cell) + rib_t / 2
    return np.append(c, width - rib_t / 2)


def hole_axis(width, pitch, cell, rib_t, hole_r):
    """Hole centres along one axis, midway between rib centrelines, with the web rule applied."""
    if abs(cell / pitch - round(cell / pitch)) > 1e-9:
        raise ValueError(f"cell {cell} must be a multiple of pitch {pitch}")
    phase = rib_t / 2 + pitch / 2
    a = phase + pitch * np.arange(int(np.floor((width - phase) / pitch)) + 1)
    ribs = rib_centres(width, cell, rib_t)
    gap = np.min(np.abs(a[:, None] - ribs[None, :]), axis=1)
    return a[gap >= hole_r + rib_t / 2 + WEB]


def min_web(axis, width, cell, rib_t, hole_r):
    """Smallest wall between any kept hole and any rib (mm)."""
    ribs = rib_centres(width, cell, rib_t)
    return float(np.min(np.abs(axis[:, None] - ribs[None, :])) - hole_r - rib_t / 2)
