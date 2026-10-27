"""Compute the curves and field data, then fill template.html into the visualization page."""
import json
import os
import numpy as np
from model import stack_field, band_mean, C
from structure import plate_modes
from geometry import hole_axis, rib_centres, min_web

HERE = os.path.dirname(os.path.abspath(__file__))

# Design (see tile.py / README): best printable design found by the parameter search.
GEOM = dict(hole_mm=0.5, t_mm=1.0, p1_mm=2.5, p2_mm=5.0, g1_mm=25.0, g2_mm=30.0,
            cell_mm=20.0, back_mm=5.0, width_mm=120.0, rib_t_mm=1.2)
HOLE = GEOM['hole_mm'] * 1e-3
T = GEOM['t_mm'] * 1e-3
GAPS = (GEOM['g1_mm'] * 1e-3, GEOM['g2_mm'] * 1e-3)
PITCHES = (GEOM['p1_mm'] * 1e-3, GEOM['p2_mm'] * 1e-3)      # front leaf, rear leaf
LAYERS = [(HOLE, T, PITCHES[0], GAPS[0]), (HOLE, T, PITCHES[1], GAPS[1])]

W = GEOM['width_mm']
AX1 = hole_axis(W, GEOM['p1_mm'], GEOM['cell_mm'], GEOM['rib_t_mm'], GEOM['hole_mm'] / 2)
AX2 = hole_axis(W, GEOM['p2_mm'], GEOM['cell_mm'], GEOM['rib_t_mm'], GEOM['hole_mm'] / 2)
RIBS = rib_centres(W, GEOM['cell_mm'], GEOM['rib_t_mm'])
MIN_WEB = min(min_web(AX1, W, GEOM['cell_mm'], GEOM['rib_t_mm'], GEOM['hole_mm'] / 2),
              min_web(AX2, W, GEOM['cell_mm'], GEOM['rib_t_mm'], GEOM['hole_mm'] / 2))
open_frac = [1.0, 1.0]        # holes are placed clear of ribs, so none are blocked (see geometry.py)


f = np.geomspace(300, 3000, 240)
field = stack_field(f, LAYERS)
alpha, refl = field['alpha'], field['refl']
built = stack_field(f, LAYERS, open_frac)
alpha_built = built['alpha']

print('holes front %d rear %d, min web to rib %.2f mm' % (len(AX1) ** 2, len(AX2) ** 2, MIN_WEB))
print('plate-mode check (leaf cell, first mode Hz):', plate_modes(GEOM['cell_mm'] * 1e-3, GEOM['cell_mm'] * 1e-3, T)[0][0])
print('band mean alpha %.3f  min %.3f  max %.3f' % (band_mean(alpha, f), alpha.min(), alpha.max()))
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
    'geom': GEOM,
    'axis1': rnd(AX1, 3),
    'axis2': rnd(AX2, 3),
    'ribs': rnd(RIBS, 3),
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
