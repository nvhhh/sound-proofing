"""Compute the curves and field data, then fill template.html into the visualization page."""
import json
import os
import numpy as np
from model import stack_field, band_mean, C

HERE = os.path.dirname(os.path.abspath(__file__))

HOLE = 0.5e-3
T = 0.8e-3
GAPS = (25e-3, 25e-3)
PITCHES = (3e-3, 6e-3)                      # front leaf, rear leaf
LAYERS = [(HOLE, T, PITCHES[0], GAPS[0]), (HOLE, T, PITCHES[1], GAPS[1])]

# Rib layout from build_stl.py (mm): 1.2 mm ribs on a 20 mm grid over a 120 mm coupon.
SIZE_MM, RIB_T_MM, RIB_PITCH_MM, HOLE_R_MM = 120.0, 1.2, 20.0, 0.25


def open_fraction(pitch_mm):
    """Fraction of holes that do not sit on a rib (a hole on a rib opens onto solid plastic)."""
    n = int(SIZE_MM // pitch_mm)
    off = (SIZE_MM - (n - 1) * pitch_mm) / 2
    c = off + pitch_mm * np.arange(n)
    ribs = np.array([min(i * RIB_PITCH_MM, SIZE_MM - RIB_T_MM) for i in range(int(SIZE_MM // RIB_PITCH_MM) + 1)])
    over = np.any((ribs[None, :] - HOLE_R_MM < c[:, None] + HOLE_R_MM) &
                  (c[:, None] - HOLE_R_MM < ribs[None, :] + RIB_T_MM), axis=1)
    blocked = over[:, None] | over[None, :]
    return float(1 - blocked.mean())


f = np.geomspace(300, 3000, 240)
field = stack_field(f, LAYERS)
alpha, refl = field['alpha'], field['refl']
open_frac = [open_fraction(p * 1e3) for p in PITCHES]
built = stack_field(f, LAYERS, open_frac)
alpha_built = built['alpha']

print('open fraction front %.3f rear %.3f' % tuple(open_frac))
print('band mean alpha %.3f  min %.3f  max %.3f' % (band_mean(alpha, f), alpha.min(), alpha.max()))
print('band mean alpha, holes under ribs removed %.3f' % band_mean(alpha_built, f))
share = field['share']
print('band mean share front %.3f rear %.3f' % (band_mean(share[0], f), band_mean(share[1], f)))
for hz in [300, 500, 800, 1000, 1500, 2000, 2500, 3000]:
    print(hz, round(float(np.interp(hz, f, alpha)), 2))

def rnd(a, n=5):
    return np.round(np.asarray(a), n).tolist()

data = {
    'f': rnd(f, 2),
    'alpha': rnd(alpha, 4),
    'alpha_built': rnd(alpha_built, 4),
    'reR': rnd(refl.real),
    'reI': rnd(refl.imag),
    'share_front': rnd(share[0], 4),
    'share_rear': rnd(share[1], 4),
    'cav1': [rnd(v) for v in (field['cav'][0][0].real, field['cav'][0][0].imag,
                              field['cav'][0][1].real, field['cav'][0][1].imag)],
    'cav2': [rnd(v) for v in (field['cav'][1][0].real, field['cav'][1][0].imag,
                              field['cav'][1][1].real, field['cav'][1][1].imag)],
    'gap_mm': [g * 1e3 for g in GAPS],
    'open': [round(o, 4) for o in open_frac],
    'mean': round(band_mean(alpha, f), 3),
    'mean_built': round(band_mean(alpha_built, f), 3),
}
with open(os.path.join(HERE, 'data.json'), 'w') as fh:
    json.dump(data, fh)
print('json bytes', len(json.dumps(data)))

with open(os.path.join(HERE, 'template.html')) as fh:
    page = fh.read().replace('__DATA__', json.dumps(data))
with open(os.path.join(HERE, 'mpp_absorber_visualization.html'), 'w') as fh:
    fh.write(page)
print('wrote mpp_absorber_visualization.html')
