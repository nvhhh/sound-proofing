# Micro-perforated sound absorber (300–3000 Hz)

A single-material, 3D-printable absorber: two perforated leaves, each with an air cavity behind it.
Best printable design found so far: 120 × 120 × 92 mm coupon, PLA, PETG or ABS.

## Files
- `model.py`: micro-perforated panel model at normal incidence. Default leaf impedance is the exact viscous tube (Bessel) with an end correction of 0.68 d. `leaf='maa98'` gives the older Maa (1998) fit. `stack_field` returns absorption, reflection, the energy split per leaf and the field inside each cavity.
- `geometry.py`: shared hole and rib layout for the coupon and the tile. Holes sit midway between rib centrelines with at least 0.4 mm of wall to each rib.
- `build_data.py`: computes the curves and field data for the best design, writes `data.json`, and fills `template.html` into `mpp_absorber_visualization.html`. Run `python3 build_data.py` (numpy, scipy).
- `template.html`: bench page template. The geometry and data are read from `data.json`.
- `tile.py`: parametric, stackable tile with screw and tape mounting (`python3 tile.py --bed 220`, needs trimesh, manifold3d).
- `structure.py`: plate-vibration and printing-material checks.
- `measure.py`: two-microphone impedance-tube analysis (ISO 10534-2).
- `tests/`: unit tests (`python3 -m pytest tests`).
- `build_stl.py`, `mpp_absorber.stl`: earlier single-coupon STL (previous design, 0.8 mm leaves, 3/6 mm pitch). Regenerate from `tile.py` for the current design.

## Best design
| | Front leaf | Air gap 1 | Rear leaf | Air gap 2 | Back |
|---|---|---|---|---|---|
| Thickness / size | 1.0 mm | 25 mm | 1.0 mm | 60 mm (see neck) | 5 mm |
| Holes | 0.5 mm on 2.5 mm pitch, 2209 | ribs 1.2 mm on 20 mm | 0.5 mm on 5 mm pitch, 576 | ribs 1.2 mm on 20 mm | solid |
| Perforation ratio σ | 3.1% | | 0.8% | | |

Rear cells are necked: each 20 mm cell ends in a 4 mm bore, 5 mm long neck through the back wall (Helmholtz resonator). The coupon is 120 × 120 × 92 mm.

Rule for printability: holes ≥ 0.5 mm (0.4 mm nozzle), leaves 1.0 mm (three perimeters), pitch a divisor of the 20 mm rib cell, at least 0.4 mm of wall between every hole and every rib.

## Performance (model, not measured)
| Design | Band mean α (300–3000 Hz) | Worst α |
|---|---|---|
| Previous (0.8 mm leaves, 3/6 mm pitch, 25/25 mm gaps) | 0.65 | 0.20 |
| **Best (table above)** | **0.66** | **0.24** |
| 0.5 mm, 2.5/5 mm, 1.2 mm leaves, 25/30 mm | 0.66 | 0.21 |

Energy split inside the stack, band average: front leaf ≈ 0.27, rear leaf ≈ 0.39 (sum = α).
Weakest near 300 Hz (α ≈ 0.28) and at 2.5–3 kHz (0.24–0.37).

## Checks made on the calculations
- The leaf impedance and the two-leaf cascade agree with an independent transfer-matrix calculation (max difference ~1e-15).
- Energy balance: heat in the two leaves equals 1 − |R|² at every frequency; the wall velocity is zero.
- End correction 0.68 d comes from a 2-D resolved Stokes simulation of one hole (0.66–0.70 d). The usual 0.85 d overestimates the mass term by about 20%. The simulation used an equal-area cylinder, not the square array, and is not full 3-D.
- Plate vibration (`structure.py`): the 20 mm cell leaf's first mode is about 4.1 kHz and the rib walls' is about 4.6 kHz, both above the band. A clamped leaf moves at under 1% of the air velocity, so its mass does not enter the impedance. Mass-law re-radiation from the leaf is about 5% of incident intensity at 300 Hz and under 0.5% above 1 kHz.
- Holes on ribs: none. Earlier designs placed holes under ribs (231 of 1600 front holes); the current layout clears every rib.

## Printing material (`structure.py`)
| Material | Leaf first mode (18 mm cell, 0.8 mm) | Leaf velocity / air velocity (max) |
|---|---|---|
| PLA | 4.0 kHz | 0.006 |
| PETG | 3.1 kHz | 0.010 |
| ABS | 3.6 kHz | 0.009 |
| Nylon (PA) | 2.8 kHz | 0.014 |

The acoustic result does not depend on material; the air in the holes sets the impedance. Material matters for stiffness, print fidelity (warp, hole size) and moisture (nylon). Nylon's first mode falls inside the upper edge of the band. Re-check the leaf modes for the 1.0 mm, 20 mm cell design before printing with any material.

## Limits
- The model is 1-D normal incidence. It does not include oblique sound, room-mode coupling or the diffuse field.
- The 0.68 d end correction is from one 2-D simulation. The 3-D hole-mouth flow is not resolved.
- Structure-borne transmission through the mounts is not modelled.
- Nothing has been measured. Print a coupon and measure it in an impedance tube (ASTM E1050 / ISO 10534-2) with `measure.py`.
- 0.5 mm holes at a 0.4 mm nozzle need a test print to confirm diameter; the 0.4 mm web is at the printer's minimum.

## Room estimates (`room_sim.html`)
Sabine estimate, tiles on the walls, 0.8 diffuse correction, untreated absorption assumed typical.

| Room | Tiled area | Tiles | RT60 500 Hz before → after | RT60 1 kHz before → after |
|---|---|---|---|---|
| Living room 5 × 4 × 2.5 m | 6 m² | 417 | 0.32 → 0.29 s | 0.27 → 0.26 s |
| Studio 6 × 5 × 3 m | 12 m² | 833 | 0.96 → 0.68 s | 0.77 → 0.62 s |

A studio target of 0.3–0.4 s at 500 Hz to 2 kHz needs about 30 m² of tile on this estimate. Low frequencies (125–250 Hz) need bass traps; 62 mm tiles do little there.

### 240 × 240 mm tiles
- Geometry: `python3 tile.py --cells 12` builds a 240 × 240 × 62 mm tile, watertight. The rib grid (20 mm cells) and hole layout are the same, so the leaf and wall modes are unchanged.
- Holes per tile: 9025 front, 2304 rear (vs. 2209 and 576 for 120 mm).
- Solid PLA weight: about 1.0 kg per 240 mm tile (0.27 kg for 120 mm, solid print volume before infill). Use screws, not tape, for the larger tile.
- Room estimate: the model is per unit area, so the Sabine times do not change with tile size. Living room 6 m² (104 tiles): 0.29 s at 500 Hz. Studio 12 m² (208 tiles): 0.68 s; 30 m² (521 tiles): 0.47 s.
- Not captured by the 1-D model: edge diffraction and tile-to-tile seams. These matter more for the larger tile, and need an impedance-tube test of a 240 mm sample (or a 100 mm tube with a 120 mm sample) to check.

## Necked rear cells (low end, 60 mm rear cavity)
| Rear design | Depth | 120–250 Hz α | 300–3000 Hz mean α |
|---|---|---|---|
| Rigid back, 30 mm | 62 mm | 0.08 | 0.66 |
| Rigid back, 60 mm | 92 mm | 0.27 | 0.64 |
| Necked, 60 mm, 4 mm bore × 5 mm (**bench design**) | 92 mm | 0.51 | 0.53 |
| Necked, 60 mm, 5 mm bore × 5 mm | 92 mm | 0.59 | 0.51 |

Necks trade mid-band absorption for low-end absorption at the same depth. At 100 mm rigid depth, the rigid cavity (0.57 at 120–250 Hz) beats the neck; the necks are for when the depth is limited to about 60 mm. Modelled with the exact tube impedance for the neck (`neck_z` in `model.py`), a resonance-broadened Helmholtz check; not measured.

Energy: the neck's power goes into the wall space behind the tile, and the model counts it as absorbed (about 3% of the band energy). That assumes the wall space takes up the sound and does not return it to the room. A sealed or reflective back would change this.
