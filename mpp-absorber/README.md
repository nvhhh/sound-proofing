# Micro-perforated sound absorber (300–3000 Hz), proof of concept

A single-material, 3D-printable absorber: a two-leaf micro-perforated panel
(PLA or PETG). The coupon is 120 × 120 × 53.6 mm.

## Files
- `mpp_absorber.stl`: watertight printable model.
- `mpp_absorber_visualization.html`: interactive 3D view (explode the panel, watch the wave in the room and in each air gap) and absorption curve. Open in a browser.
- `model.py`: micro-perforated panel model at normal incidence. Default is the exact viscous tube (Bessel) impedance with an end correction of 0.85 d; the Maa (1998) fit is kept as `leaf='maa98'`. Also gives the field and energy split through the stack (`stack_field`).
- `measure.py`: two-microphone impedance-tube analysis (ISO 10534-2) from measured transfer functions.
- `tests/`: unit tests (`python3 -m pytest tests`).
- `build_stl.py`: generates the STL (needs `numpy trimesh manifold3d shapely mapbox_earcut`).
- `build_data.py`: computes the absorption, reflection and field data, writes `data.json`, and fills `template.html` into `mpp_absorber_visualization.html`. Run `python3 build_data.py` (needs numpy).
- `template.html`: page template; the data goes in at `__DATA__`.

## Design
- Front leaf: 0.8 mm, 0.5 mm holes, 3 mm pitch.
- Air gap 25 mm, ribbed into 20 mm cells.
- Rear leaf: 0.8 mm, 0.5 mm holes, 6 mm pitch.
- Air gap 25 mm, ribbed.
- Back wall: 2 mm.

Perforation ratio: front σ ≈ 2.2% (1600 holes), rear σ ≈ 0.55% (400 holes).

## Performance (model, not measured)
Band-average absorption over 300–3000 Hz ≈ 0.65 (exact tube model, end correction 0.85 d). The older Maa (1998) fit gave 0.56. Weak near 300 Hz and above 2.5 kHz.

Energy split inside the stack, band average: front leaf ≈ 0.29, rear leaf ≈ 0.36 (sum = α). The rest is reflected.

## Checks made on the calculations
- The Maa impedance, the air-gap impedance and the two-leaf cascade agree with an independent transfer-matrix calculation (max difference ~1e-15).
- Energy balance: the heat dissipated in the two leaves equals 1 − |R|² at every frequency. The wall velocity is zero, so the boundary condition holds.
- Holes on ribs: 231 of 1600 front holes and 76 of 400 rear holes sit over a rib and are blocked. The model ignores this. Removing them lowers the band average to about 0.64, but it is wasted area. `open_fraction` in `build_data.py` computes these counts; the as-built curve (dashed on the chart) uses them.
- End correction 0.85 d is a standard value, not measured. The 3-D mouth flow is not resolved.

## Validation still needed
Print a coupon and measure it in an impedance tube (ASTM E1050 / ISO 10534-2).
0.5 mm holes are at the limit of a 0.4 mm nozzle; check hole diameter first.

## Stackable parametric tile (`tile.py`)
A square tile of the absorber, sized by rib cells so it fits your printer.

    python3 tile.py                    # 6 x 6 cells = 108 mm, screws and tape
    python3 tile.py --bed 220          # largest tile that fits a 220 mm bed
    python3 tile.py --cells 3 --mount tape

- Size: `cells` x 18 mm. Holes sit between the rib lines, so none are blocked by a rib.
- Depth: 56.6 mm (two 0.8 mm leaves, two 25 mm cavities, 5 mm back).
- Stacking: tiles butt edge to edge. Two 3 mm dowels (about 6 mm into each side socket) align neighbours.
- Mounting: `screw` uses four countersunk M3 holes in the corner bosses; `tape` leaves the back flat for tape; `both` does both.
- Corner bosses (8 mm) carry the screws and remove the holes under them.
- The generator checks watertightness, but the design is not yet printed or measured.
