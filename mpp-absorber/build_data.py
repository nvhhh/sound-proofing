import json
import numpy as np
from model import double_mpp, band_mean, C

HOLE = 0.5e-3
T = 0.8e-3
LAYERS = [(HOLE, T, 3e-3, 25e-3), (HOLE, T, 6e-3, 25e-3)]   # front leaf, rear leaf
f = np.geomspace(300, 3000, 240)
alpha, refl = double_mpp(f, LAYERS)
print('band mean alpha %.3f  min %.3f  max %.3f' % (band_mean(alpha, f), alpha.min(), alpha.max()))
for hz in [300, 500, 800, 1000, 1500, 2000, 2500, 3000]:
    print(hz, round(float(np.interp(hz, f, alpha)), 2))
data = {
    'f': np.round(f, 2).tolist(),
    'alpha': np.round(alpha, 4).tolist(),
    'reR': np.round(refl.real, 5).tolist(),
    'reI': np.round(refl.imag, 5).tolist(),
    'mean': round(band_mean(alpha, f), 3),
}
json.dump(data, open('data.json', 'w'))
print('json bytes', len(json.dumps(data)))
