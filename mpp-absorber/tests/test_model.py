import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from model import stack_field, band_mean, C, RHOC  # noqa: E402
from measure import alpha_from_tf  # noqa: E402

F = np.geomspace(300, 3000, 240)
LAYERS = [(0.5e-3, 0.8e-3, 3e-3, 25e-3), (0.5e-3, 0.8e-3, 6e-3, 25e-3)]


def test_leaf_shares_sum_to_absorption():
    out = stack_field(F, LAYERS)
    assert np.max(np.abs(out['share'].sum(0) - out['alpha'])) < 1e-12


def test_wall_velocity_is_zero():
    out = stack_field(F, LAYERS)
    P, Q = out['cav'][-1]
    k = 2 * np.pi * F / C
    p = P * np.exp(-1j * k * 25e-3) + Q * np.exp(1j * k * 25e-3)
    u = (P * np.exp(-1j * k * 25e-3) - Q * np.exp(1j * k * 25e-3)) / RHOC
    assert np.max(np.abs(u)) < 1e-9


def test_absorption_bounded():
    a = stack_field(F, LAYERS)['alpha']
    assert np.all((a >= 0) & (a <= 1 + 1e-12))


def test_two_microphone_recovers_known_reflection():
    s, x1 = 0.03, 0.05
    k = 2 * np.pi * F / C
    R_true = 0.6 * np.exp(0.7j)
    # Standing-wave field in front of the sample: p = e^{-jkx} + R e^{jkx}, x from the face.
    def p(x):
        return np.exp(-1j * k * x) + R_true * np.exp(1j * k * x)
    H12 = p(x1 + s) / p(x1)
    alpha, R = alpha_from_tf(H12, F, s, x1)
    assert np.allclose(np.abs(R), 0.6, atol=1e-9)
    assert np.allclose(alpha, 1 - 0.36, atol=1e-9)


def test_leaf_cell_stays_still_in_band():
    from structure import plate_velocity_ratio
    assert np.max(plate_velocity_ratio(F, 0.018, 0.8e-3)) < 0.02


def test_leaf_cell_mode_above_band():
    from structure import plate_modes
    assert plate_modes(0.018, 0.018, 0.8e-3)[0][0] > 3000


def test_wall_mode_above_band():
    from structure import plate_modes
    assert plate_modes(0.025, 0.018, 1.2e-3)[0][0] > 3000


def test_holes_clear_ribs_with_web():
    from geometry import hole_axis, min_web
    for pitch in (2.5, 5.0):
        a = hole_axis(120, pitch, 20, 1.2, 0.25)
        assert min_web(a, 120, 20, 1.2, 0.25) >= 0.4 - 1e-9




def test_neck_energy_balance():
    neck = dict(r=2e-3, length=5e-3, cell=20e-3)
    layers = [(0.5e-3, 1e-3, 2.5e-3, 25e-3), (0.5e-3, 1e-3, 5e-3, 60e-3)]
    out = stack_field(F, layers, neck=neck)
    assert np.max(np.abs(out['share'].sum(0) - out['alpha'])) < 1e-12
    assert np.all((out['alpha'] >= 0) & (out['alpha'] <= 1 + 1e-12))
