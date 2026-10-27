"""Methods check: robustness of the products to parameter uncertainty, and a Delany-Bazley helper.

Uncertain inputs sampled uniformly for each product:
  end correction factor 0.60-0.85 d (the 2-D simulation range plus the common 0.85 d value)
  every cavity depth +/- 2 mm
Reports median and 5th percentile of the 125 Hz - 4 kHz mean, so a design is judged on its
worst plausible case, not its nominal value.
"""
import numpy as np

from model import stack_field, C, RHO, RHOC
from products import PRODUCTS, F, oct_mean

BANDS = (125, 250, 500, 1000, 2000, 4000)


def db_alpha(f, sigma, depth):
    """Normal-incidence absorption of a porous layer (Delany-Bazley) over a rigid back."""
    f = np.atleast_1d(np.asarray(f, dtype=float))
    X = RHO * f / sigma
    Zn = 1 + 0.0571 * X**-0.754 - 1j * 0.087 * X**-0.732
    kn = 1 + 0.0978 * X**-0.700 - 1j * 0.189 * X**-0.595
    k = 2 * np.pi * f / C * kn
    Zin = -1j * RHOC * Zn / np.tan(k * depth)
    return 1 - np.abs((Zin - RHOC) / (Zin + RHOC)) ** 2


def sample_layers(layers, rng):
    out = []
    for d, t, p, g in layers:
        out.append((d, t, p, max(g + rng.uniform(-2e-3, 2e-3), 1e-3)))
    return out


def robust(name, n=40, seed=0):
    rng = np.random.default_rng(seed)
    p = PRODUCTS[name]
    means, octs = [], []
    for _ in range(n):
        tf = rng.uniform(0.60, 0.85)
        layers = sample_layers(p['layers'], rng)
        a = stack_field(F, layers, neck=p['neck'], t_end_factor=tf)['alpha']
        octs.append([oct_mean(a, b) for b in BANDS])
        means.append(np.mean(octs[-1]))
    means = np.array(means)
    worst_band = np.min(np.array(octs), axis=1)
    return dict(median=float(np.median(means)), p05=float(np.percentile(means, 5)),
                worst_band_p05=float(np.percentile(worst_band, 5)))


if __name__ == '__main__':
    for name in PRODUCTS:
        r = robust(name)
        print('%-14s median %.2f  5th percentile %.2f  worst-band 5th pct %.2f' % (
            name, r['median'], r['p05'], r['worst_band_p05']))
