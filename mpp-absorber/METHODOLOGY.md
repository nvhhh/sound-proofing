# Methodology

Every number in this project is a model estimate until it has been measured. Each claim carries one of three evidence levels.

| Level | Meaning | Examples |
|---|---|---|
| **V** Validated | Checked against an exact limit, an independent calculation, or a resolved simulation | cascade vs transfer matrix (1e-15); energy balance; wall velocity zero; neck passivity; Delany–Bazley vanishes as resistivity grows; end correction 0.60–0.70 d from a 2-D Stokes simulation |
| **A** Approximate | Engineering model with known simplifications; reproduces a reference only partly | wide-slit slot model (bass face mean 0.80 vs report 0.866); Helmholtz neck; Kozeny–Carman; plate-mass membrane |
| **U** Unmeasured | Assumed input or unchecked output | foam flow resistivity (eSUN TPU-LW, Filaflex Foamy); room absorption baselines; diffuse factor 0.8; 125 Hz and 4 kHz tile values in the room model; the 3-D hole mouth; structure-borne path |

No product or design should be described as having a measured performance until it has passed the gates below.

## 1. Models (`model.py`)
- Perforated leaves: exact viscous tube (Bessel) with end correction 0.68 d (V). Maa (1998) fit kept for comparison.
- Cavities: rigid-backed transmission line; optional Helmholtz neck termination (A for the mass and loss, V for passivity).
- Slots: wide-slit thermoviscous duct with isolated-slit end correction (A).
- Fill: Delany–Bazley porous layer (A for the model form, U for the resistivity).
- Energy: the leaf shares, and the neck share, sum to the absorption (V). The neck share assumes the wall space behind the tile absorbs the power (U).

## 2. Validation tests (`tests/`)
Run `python3 -m pytest tests`. The suite includes exact-limit checks (energy balance, zero wall velocity, neck passivity, Delany–Bazley rigid limit), an independent cascade check, a two-microphone recovery test, geometry clearance, and plate-mode checks. The bass-face test is a regression with a documented tolerance, not a validation.

## 3. Robustness (`methods_check.py`)
Each product is sampled over end correction 0.60–0.85 d and every cavity depth ±2 mm (40 draws, seeded). Reported: median and 5th percentile of the 125 Hz–4 kHz mean, and the 5th percentile of the worst band.
- Result: parameter uncertainty moves the mean by under 0.01. The worst-band value is more sensitive, and P2's worst band falls to 0.03 at the 5th percentile.
- This underestimates real uncertainty, because model-form errors (slot model, neck wall space, foam resistivity) are not sampled. Model-form uncertainty is the larger risk, and it is reported separately.

## 3b. Band definitions
Targets are scored on the linear band they name (for example 125–250 Hz, 500–4000 Hz), not on octave bins. An octave bin for 125 Hz spans about 88–177 Hz, which mixes in frequencies outside the target. This changed the bass face from 0.80 to 0.64 in the same model.

## 4. Design optimisation (`products.py`, `geometry.py`)
- Objective: maximise the worst octave band over the target range first, then the mean. A mean-only objective hides dead bands and picked four bass traps for the 240 mm bed.
- Constraints: holes ≥ 0.5 mm, leaves 0.8–1.2 mm, pitch divides the 20 mm rib cell, ≥ 0.4 mm web, depth limits.
- Search is grid or random over these constraints. Each candidate is scored on the model only, so the optimum is an optimum of the model, not of the physical product.
- Before a design is used, its robustness score (section 3) should be checked and its worst band reported.

## 5. Gates before any performance claim
1. **Coupon measurement** (ISO 10534-2, `measure.py`): one perforated coupon and one bass cell measured, with the model curve plotted on the same axes. Accept a band mean within ±0.08 of the model, and a peak frequency within ±10%.
2. **Foam flow resistivity** measured on a printed sample before any foam product is rated.
3. **Slot model** checked against one printed bass cell before the SCAD bass face is rated.
4. **Room measurement** (ISO 3382-2) of at least one installed panel before any room-level claim.

Until these gates pass, report model results as "model, not measured" with their evidence level.
