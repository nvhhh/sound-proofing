# Micro-perforated sound absorber (300–3000 Hz), proof of concept

A single-material, 3D-printable absorber: a two-leaf micro-perforated panel
(PLA or PETG). The coupon is 120 × 120 × 53.6 mm.

## Files
- `mpp_absorber.stl`: watertight printable model.
- `mpp_absorber_visualization.html`: interactive 3D view of the wave and absorption curve (open in a browser).
- `model.py`: Maa (1998) micro-perforated panel impedance model, normal incidence.
- `build_stl.py`: generates the STL (needs `numpy trimesh manifold3d shapely mapbox_earcut`).
- `build_data.py`: computes the absorption and reflection data embedded in the HTML.
- `template.html`: HTML template; `build_data.py` output (`data.json`) is inserted at `__DATA__`.

## Design
- Front leaf: 0.8 mm, 0.5 mm holes, 3 mm pitch.
- Air gap 25 mm, ribbed into 20 mm cells.
- Rear leaf: 0.8 mm, 0.5 mm holes, 6 mm pitch.
- Air gap 25 mm, ribbed.
- Back wall: 2 mm.

## Performance (model, not measured)
Band-average absorption over 300–3000 Hz ≈ 0.56. Weak near 300 Hz and above 2.5 kHz.

## Validation still needed
Print a coupon and measure it in an impedance tube (ASTM E1050 / ISO 10534-2).
0.5 mm holes are at the limit of a 0.4 mm nozzle; check hole diameter first.
