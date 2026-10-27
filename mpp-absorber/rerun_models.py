"""Rerun every proposed design under the methodology in METHODOLOGY.md.

For each design: nominal octave-band absorption (63 Hz to 4 kHz), robustness over its
uncertain inputs (seeded, 40 draws), and the evidence level of its model.
Output: rerun_results.json and a printed table.
"""
import json
import os

import numpy as np

from model import stack_field, neck_z, slit_absorption, C
from products import PRODUCTS, bed_mix, product_bands, F, oct_mean, BANDS
from methods_check import db_alpha, sample_layers

HERE = os.path.dirname(os.path.abspath(__file__))
N = 40
RNG_SEED = 0
OCT = [63, 125, 250, 500, 1000, 2000, 4000]


def stack_designs():
    """Designs that are layer stacks: name -> (layers, neck, evidence, notes)."""
    bench = dict(
        layers=[(0.5e-3, 1.0e-3, 2.5e-3, 25e-3), (0.5e-3, 1.0e-3, 5e-3, 60e-3)],
        neck=dict(r=2e-3, length=5e-3, cell=20e-3), level='V/A',
        note='bench design: 25/60 mm, 4 mm bore necks')
    out = {'bench (necked 25/60)': bench}
    for name, p in PRODUCTS.items():
        out[name] = dict(layers=p['layers'], neck=p['neck'], level='V/A', note=p['note'])
    return out


def foam_combined(sigma=10e3, depth=30e-3):
    """P2-style three-leaf stack with a 30 mm printed foam fill in the front cavity (U: sigma)."""
    return dict(layers=[(0.5e-3, 1e-3, 2e-3, depth), (0.5e-3, 1e-3, 2.5e-3, 25e-3),
                        (0.5e-3, 1e-3, 5e-3, 20e-3)], neck=None, level='U',
                note='P2 with foam fill, sigma 10 kRayl/m (assumed)', foam=(sigma, depth))


def stack_alpha(design, f, layers=None, t_end=0.68):
    layers = design['layers'] if layers is None else layers
    if 'foam' in design:
        sigma, depth = design['foam']
        return foam_alpha(f, layers, sigma, depth, t_end)
    return stack_field(f, layers, neck=design['neck'], t_end_factor=t_end)['alpha']


def foam_alpha(f, layers, sigma, depth, t_end):
    """Front cavity filled with a Delany-Bazley porous layer; rear cavity is rigid (as P2)."""
    from model import _mpp_z_exact, RHOC
    w = 2 * np.pi * f
    d, t, p1, _ = layers[0]
    z1 = _mpp_z_exact(w, d, t, p1, t_end * d)
    # rear part: the last two leaves in air, rigid back (P2 structure)
    zr = stack_field(f, layers[1:], t_end_factor=t_end)
    z_rear = RHOC * (1 + zr["refl"]) / (1 - zr["refl"])        # impedance at rear-leaf top: z = rho c (1+R)/(1-R)
    X = 1.21 * f / sigma
    Zn = 1 + 0.0571 * X**-0.754 - 1j * 0.087 * X**-0.732
    kn = 1 + 0.0978 * X**-0.700 - 1j * 0.189 * X**-0.595
    k = w / C * kn
    zt = RHOC * Zn * (z_rear + 1j * RHOC * Zn * np.tan(k * depth)) / (RHOC * Zn + 1j * z_rear * np.tan(k * depth))
    zf = z1 + zt
    R = (zf - RHOC) / (zf + RHOC)
    return 1 - np.abs(R) ** 2


def slit_design():
    """The SCAD bass face: four slot types, 79.62% of face in cells (A: wide-slit model)."""
    types = [(1.9e-3, 40e-3, 44.2e-3), (1.2e-3, 40e-3, 17.8e-3), (0.9e-3, 40e-3, 7.4e-3), (0.75e-3, 40e-3, 3.2e-3)]
    return dict(types=types, D=112e-3, cell=48e-3, level='A',
                note='SCAD bass face, wide-slit model (report method)')


def octave_table(a, f):
    return {b: float(a[(f >= b / np.sqrt(2)) & (f <= b * np.sqrt(2))].mean()) for b in OCT}


def robust_stack(design, rng):
    means, lows, mids = [], [], []
    for _ in range(N):
        t = rng.uniform(0.60, 0.85)
        layers = sample_layers(design['layers'], rng)
        a = stack_alpha(design, F, layers=layers, t_end=t)
        ob = octave_table(a, F)
        means.append(np.mean([ob[b] for b in (125, 250, 500, 1000, 2000, 4000)]))
        lows.append(ob[125] * 0 + np.mean([ob[125], ob[250]]))
        mids.append(np.mean([ob[b] for b in (250, 500, 1000, 2000, 4000)]))
    return dict(median=float(np.median(means)), p05=float(np.percentile(means, 5)))


def robust_slit(design, rng):
    f = np.geomspace(63, 4000, 300)
    means = []
    for _ in range(N):
        D = design['D'] + rng.uniform(-2e-3, 2e-3)
        types = [(h + rng.uniform(-0.1e-3, 0.1e-3), s, L) for h, s, L in design['types']]
        shares = [0.7962 / 4] * 4
        a = slit_absorption(f, [(h, s, L, sh) for (h, s, L), sh in zip(types, shares)], D=D, cell=design['cell'])
        means.append(np.mean([octave_table(a, f)[b] for b in (125, 250)]))
    return dict(median=float(np.median(means)), p05=float(np.percentile(means, 5)))


def main():
    rng = np.random.default_rng(RNG_SEED)
    results = {}
    f = np.geomspace(60, 4000, 400)
    print('%-28s lvl  %s   | 125-250 mean  | 250-4k mean | robust p05 (125-4k)' % ('design', '  '.join('%5d' % b for b in OCT)))
    designs = stack_designs()
    designs['foam combined (U)'] = foam_combined()
    for name, d in designs.items():
        a = stack_alpha(d, f)
        ob = octave_table(a, f)
        r = robust_stack(d, rng)
        m_low = np.mean([ob[125], ob[250]])
        m_mid = np.mean([ob[b] for b in (250, 500, 1000, 2000, 4000)])
        results[name] = dict(level=d['level'], note=d['note'], octaves=ob, mean_125_250=m_low,
                             mean_250_4k=m_mid, robust=r)
        print('%-28s %-4s %s   | %.2f         | %.2f         | %.2f' % (
            name, d['level'], '  '.join('%5.2f' % ob[b] for b in OCT), m_low, m_mid, r['p05']))
    sd = slit_design()
    a = slit_absorption(f, [(h, s, L, 0.7962 / 4) for h, s, L in sd['types']], D=sd['D'], cell=sd['cell'])
    ob = octave_table(a, f)
    r = robust_slit(sd, rng)
    results[sd['note']] = dict(level=sd['level'], note=sd['note'], octaves=ob,
                               mean_125_250=np.mean([ob[125], ob[250]]), mean_250_4k=float('nan'), robust=r)
    print('%-28s %-4s %s   | %.2f         | (n/a)       | %.2f' % (
        'SCAD bass face (A)', 'A', '  '.join('%5.2f' % ob[b] for b in OCT), np.mean([ob[125], ob[250]]), r['p05']))
    mix = bed_mix(product_bands())
    results['bed mix (3xP1 + 1xP2)'] = dict(level='V/A', note='240 mm bed mix', bed=mix)
    print('\nbed mix 240 mm: counts %s, combined bands %s' % (mix['counts'], mix['bands']))
    with open(os.path.join(HERE, 'rerun_results.json'), 'w') as fh:
        json.dump(results, fh, indent=1, default=float)
    print('wrote rerun_results.json')


if __name__ == '__main__':
    main()
