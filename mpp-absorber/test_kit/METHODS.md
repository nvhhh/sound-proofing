# Testing methodology

## 1. Simulation (before printing)
| Check | Tool | Pass criterion |
|---|---|---|
| Leaf impedance vs. full 3-D viscous hole | scikit-fem or COMSOL/Elmer, linearized Navier-Stokes, periodic unit cell | Re and Im of Z within 10% of `model.py` at 300, 1000, 3000 Hz |
| Mouth end correction in the square array | same, full 3-D unit cell with 3x3 holes | fitted t_end within 0.66–0.70 d (current model value 0.68 d) |
| Stack absorption | `model.py` (exact tube), transfer-matrix cross-check in `tests/` | energy balance and cascade agree to 1e-12 |
| Plate and wall modes | `structure.py` (thin-plate theory) and a modal FEM (e.g. CalculiX, scikit-fem shell) | first mode above 3 kHz for all leaves and rib walls |
| Print geometry | `tile.py`, trimesh | watertight, min wall ≥ 0.4 mm, holes clear of ribs |

Software: Python 3.11+, numpy, scipy, pytest, trimesh, manifold3d. Optional: scikit-fem for the 3-D hole solve, CalculiX for plate modes, PrusaSlicer or Cura to check slicing (wall count, hole size after slicing).

## 2. Printing and dimensional checks
1. Print the hole test strip (`holes`, in `mpp_test_kit.scad`). Measure each hole with pin gauges (0.40–0.60 mm) or a microscope. Pick the smallest diameter that passes; that is the real hole size.
2. Measure leaf thickness with a micrometer at 5 points: target 1.0 mm.
3. Check warp after cooling: flatness with a granite plate and feeler gauge, target < 0.2 mm across 120 mm.
4. Print the sample (`sample`), measure its outer diameter (target 28.8 mm to fit a 29 mm tube) and cavity depths.

## 3. Impedance tube (absorption, the main test)
Standard: ISO 10534-2 (two-microphone transfer function) or ASTM E1050. Use the 29 mm tube for the sample (in `mpp_test_kit.scad`) and the 100 mm tube for larger tile samples.
Setup:
- Two ¼-inch or ½-inch microphones in the tube, spacing s set for the band (small tube: s = 10 mm and 50 mm).
- Broadband random or sine-sweep source at the closed end.
- Calibrate the two microphones in swapped positions (transfer-function calibration) before each series.
- Measure with the rigid reference first (`rigid` part): absorption must be < 0.05.
- Measure the sample in the same tube, sample flush to the tube end, sealed with clay or tape on the rim.

Analysis: `measure.py` takes a CSV of frequency and complex transfer function H12 and returns α(f). Compare α with `model.py` at the same frequencies. Report the band mean, the worst frequency and the frequencies of the peaks.

Repeat each sample three times (remove and reseat) and report the spread. Expect ±0.03 in band mean from seating alone.

Tolerances to accept:
- Band mean within ±0.08 of the model (0.66) is a pass; outside that, check hole size first (most likely cause), then the end correction.
- Peak frequency within ±10% of the model peak.

## 4. Real-life (room) checks
- Reverberation chamber or a small room: mount four tiles on a wall, measure RT60 in the 1/3-octave bands 315 Hz–2.5 kHz (ISO 3382-2) with and without the tiles. Expect a noticeable RT60 drop mid-band.
- Structure-borne check: accelerometer (or phone-based vibration app as a rough guide) on the back wall while driving the front with a 500 Hz, 1 kHz and 2 kHz tone at 80 dB SPL. Confirm the back wall does not ring.
- Mounting: screw pull-out on your wall type and tape peel test (time to failure at 20 °C, 50% RH).
- Temperature and humidity: repeat one sample after 48 h at 40 °C, 80% RH. Nylon especially: record any change in dimensions and absorption.

## 5. Reporting
For each sample record: material and batch, print settings (nozzle, layer height, walls, infill, orientation), measured hole diameter, leaf thickness, cavity depths, tube and microphone set, calibration date, α(f) table, and the model curve on the same plot. Keep raw transfer functions as CSV.
