"""The three printable absorber products and the 240 mm print-bed mix.

Each product is one 120 x 120 mm tile in PLA/PETG/ABS. A 240 x 240 mm bed holds four tiles;
the combined absorption is the area-weighted mean of the tiles in the bed (equal tile areas).
Octave bands 63 Hz to 4 kHz.
"""
import itertools
import json
import os

import numpy as np

from model import stack_field

HERE = os.path.dirname(os.path.abspath(__file__))
F = np.geomspace(50, 8000, 500)
BANDS = [63, 125, 250, 500, 1000, 2000, 4000]
NECK_BASE = dict(r=2.5e-3, length=5e-3, cell=20e-3)   # product 1: 5 mm bore, 5 mm long

PRODUCTS = {
    'P1 bass trap': dict(
        note='<250 Hz: deep rear cell with a 5 mm neck',
        layers=[(0.5e-3, 1.0e-3, 2.5e-3, 20e-3), (0.5e-3, 1.0e-3, 5.0e-3, 60e-3)],
        neck=NECK_BASE, depth_mm=87),
    'P2 mid-high': dict(
        note='250 Hz to 4 kHz: three leaves, shallow front cells',
        layers=[(0.5e-3, 1.0e-3, 2.0e-3, 10e-3), (0.5e-3, 1.0e-3, 2.5e-3, 25e-3),
                (0.5e-3, 1.0e-3, 5.0e-3, 20e-3)],
        neck=None, depth_mm=65),
    'P3 low-mid': dict(
        note='125 to 500 Hz: rigid rear cell, 70 mm',
        layers=[(0.5e-3, 1.0e-3, 3.0e-3, 20e-3), (0.5e-3, 1.0e-3, 5.0e-3, 70e-3)],
        neck=None, depth_mm=97),
}


def oct_mean(a, fc):
    lo, hi = fc / np.sqrt(2), fc * np.sqrt(2)
    m = (F >= lo) & (F <= hi)
    return float(a[m].mean())


def product_bands():
    out = {}
    for name, p in PRODUCTS.items():
        a = stack_field(F, p['layers'], neck=p['neck'])['alpha']
        out[name] = [round(oct_mean(a, fc), 3) for fc in BANDS]
    return out


def bed_mix(bands, slots=4):
    """Best counts of each product on a 240 mm bed (4 tiles) over 125 Hz to 4 kHz.

    Scored worst band first (maximin), then the mean, so no band is left dead.
    """
    names = list(bands)
    arr = np.array([bands[n] for n in names])            # rows: products, cols: bands
    idx = slice(1, None)                                  # 125 Hz to 4 kHz
    best = None
    for counts in itertools.product(range(slots + 1), repeat=len(names)):
        if sum(counts) != slots:
            continue
        comb = np.array(counts) @ arr / slots
        score = (comb[idx].min(), comb[idx].mean())
        if best is None or score > best[0]:
            best = (score, counts, comb)
    score, counts, comb = best
    score = score[1]
    return dict(counts={n: int(c) for n, c in zip(names, counts)},
                bands=[round(float(x), 3) for x in comb], mean_125_4k=round(float(score), 3),
                min_125_4k=round(float(comb[idx].min()), 3))


if __name__ == '__main__':
    bands = product_bands()
    print('product  ' + '  '.join('%5d' % b for b in BANDS))
    for n, v in bands.items():
        print('%-12s' % n + '  '.join('%5.2f' % x for x in v))
    mix = bed_mix(bands)
    print('bed 240 mm, 4 tiles:', mix['counts'])
    print('combined  ' + '  '.join('%5.2f' % x for x in mix['bands']))
    print('mean 125-4k %.3f, min %.3f' % (mix['mean_125_4k'], mix['min_125_4k']))
    out = dict(bands=BANDS, products={n: dict(note=p['note'], depth_mm=p['depth_mm'], alpha=bands[n])
                                      for n, p in PRODUCTS.items()}, bed=mix)
    with open(os.path.join(HERE, 'products.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
