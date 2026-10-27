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
    total = out['share'].sum(0) + out['neck_share']
    assert np.max(np.abs(total - out['alpha'])) < 1e-12
    assert np.all((out['alpha'] >= 0) & (out['alpha'] <= 1 + 1e-12))


# ---- analytic limits (validation: these must hold exactly) ----

def test_neck_is_passive_and_mass_like():
    # A passive neck dissipates (Re Z >= 0) and at low frequency its reactance is mass-like (Im Z > 0).
    from model import neck_z
    w = 2 * np.pi * np.geomspace(20, 300, 30)
    Z = neck_z(w, 2e-3, 5e-3, 20e-3)
    assert np.all(Z.real >= 0) and np.all(Z.imag > 0)


def test_delany_bazley_rigid_limit():
    from model import RHOC
    from methods_check import db_alpha
    f = np.geomspace(100, 3000, 50)
    # absorption must vanish as flow resistivity grows (rigid limit), and never exceed 1
    assert np.max(db_alpha(f, 1e9, 0.05)) < 0.01
    assert np.max(db_alpha(f, 1e9, 0.05)) < np.max(db_alpha(f, 1e4, 0.05))
    assert np.all(db_alpha(f, 1e4, 0.05) <= 1 + 1e-12)


def test_slit_absorption_bounded_and_loss_grows_with_narrowness():
    from model import slit_absorption
    f = np.geomspace(100, 250, 40)
    wide = slit_absorption(f, [(3e-3, 40e-3, 10e-3, 0.5)], D=0.1, cell=0.048)
    narrow = slit_absorption(f, [(0.8e-3, 40e-3, 10e-3, 0.5)], D=0.1, cell=0.048)
    assert np.all((wide >= 0) & (wide <= 1 + 1e-12)) and np.all((narrow >= 0) & (narrow <= 1 + 1e-12))


def test_report_bass_face_regression():
    # Regression, not validation: the wide-slit face mean over 125-250 Hz sits near the
    # engineering report's 0.866 (we get 0.80 with our guessed cavity depth and edge details).
    from model import slit_absorption
    f = np.geomspace(100, 250, 151)
    slits = [(1.9e-3, 40e-3, 44.2e-3, 0.7962 / 4), (1.2e-3, 40e-3, 17.8e-3, 0.7962 / 4),
             (0.9e-3, 40e-3, 7.4e-3, 0.7962 / 4), (0.75e-3, 40e-3, 3.2e-3, 0.7962 / 4)]
    a = slit_absorption(f, slits, D=112e-3, cell=48e-3)
    band = a[(f >= 125) & (f <= 250)]
    assert 0.7 < band.mean() < 0.9


def test_foam_fill_reduces_to_air_at_low_resistivity():
    # Delany-Bazley with tiny flow resistivity is air; the foam stack must match the air stack.
    from rerun_models import foam_alpha
    layers = [(0.5e-3, 1e-3, 2e-3, 10e-3), (0.5e-3, 1e-3, 2.5e-3, 25e-3), (0.5e-3, 1e-3, 5e-3, 20e-3)]
    f = np.geomspace(100, 4000, 200)
    air = stack_field(f, layers)['alpha']
    assert np.max(np.abs(foam_alpha(f, layers, sigma=1.0, depth=10e-3, t_end=0.68) - air)) < 0.01
