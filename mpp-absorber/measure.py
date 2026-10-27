"""Impedance-tube analysis: absorption from a two-microphone transfer function.

Standard two-microphone method (ISO 10534-2 / ASTM E1050), normal incidence.
Inputs are complex transfer functions H12 = p2/p1 measured at the two microphones.

    s       microphone spacing (m)
    x1      distance from the sample face to the nearer microphone (m)
    C       speed of sound (m/s)

Usage from a CSV with columns f, H12_re, H12_im:
    python3 measure.py data.csv --s 0.03 --x1 0.05
"""
import argparse
import csv

import numpy as np

from model import C


def reflection_from_tf(H12, f, s, x1, c=C):
    k = 2 * np.pi * np.asarray(f, dtype=float) / c
    H12 = np.asarray(H12, dtype=complex)
    r = (H12 - np.exp(-1j * k * s)) / (np.exp(1j * k * s) - H12)
    return r * np.exp(2j * k * x1)


def alpha_from_tf(H12, f, s, x1, c=C):
    R = reflection_from_tf(H12, f, s, x1, c)
    return 1 - np.abs(R) ** 2, R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--s", type=float, required=True, help="mic spacing, m")
    ap.add_argument("--x1", type=float, required=True, help="sample to nearer mic, m")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.csv)))
    f = np.array([float(r["f"]) for r in rows])
    H = np.array([float(r["H12_re"]) + 1j * float(r["H12_im"]) for r in rows])
    alpha, _ = alpha_from_tf(H, f, a.s, a.x1)
    for fi, al in zip(f, alpha):
        print(f"{fi:8.1f} Hz  alpha {al:.3f}")


if __name__ == "__main__":
    main()
