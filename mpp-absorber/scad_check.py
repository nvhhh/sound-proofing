"""Run the absorbers.scad design (bass slot face + HF channel banks) through the methodology.

Bass face: wide-slit model (A), as in rerun_models.py.
HF channels: each channel is a closed, lossy duct (wide-slit thermoviscous) of width 1.0 mm
and span 10.4 mm, length L, mouth open to the front. The bank's specific impedance is the
duct impedance divided by the bank's mouth fraction; banks combine as admittances, with the
solid carrier area contributing zero admittance (A, rough: no mouth end correction, no
inter-channel coupling, no retainer lip).
"""
import numpy as np

from model import slit_absorption, RHOC, RHO, MU, C

HF_L = [148.0, 148.0, 126.8, 111.6, 98.0, 81.8, 70.2, 57.6]   # channel lengths, mm
HF_H, HF_SPAN = 1.0, 10.4                                    # channel gap and span, mm
BANK_AREA = 24 * 192                                          # mm^2 per bank
FACE = 208 * 208
MOUTHS_PER_BANK = 264
GAMMA, PR = 1.4, 0.71
BASS_H = [1.9, 1.2, 0.9, 0.75]
BASS_L = [44.2, 17.8, 7.4, 3.2]


def channel_z(w, L, h, span):
    """Closed lossy duct, specific impedance at its mouth (Pa s/m)."""
    def fvft(x):
        q = (x / 2) * np.sqrt(1j * w * RHO / MU)
        return np.tanh(q) / q, np.tanh(q * np.sqrt(PR)) / (q * np.sqrt(PR))
    Fvh, Fth = fvft(h)
    Fvb, Ftb = fvft(span)
    Fv = Fvh + Fvb - Fvh * Fvb
    Ft = Fth + Ftb - Fth * Ftb
    rho_eff = RHO / (1 - Fv)
    K = RHO * C**2 / (1 + (GAMMA - 1) * Ft)
    Zc = np.sqrt(rho_eff * K)
    keff = w * np.sqrt(rho_eff / K)
    return -1j * Zc / np.tan(keff * L)


def hf_bank_alpha(f, L_mm, rng=None):
    """HF bank absorption; rng adds the robustness draws (channel length +/-0.5 mm, gap +/-0.05 mm)."""
    L = L_mm if rng is None else L_mm + rng.uniform(-0.5, 0.5)
    h = HF_H if rng is None else HF_H + rng.uniform(-0.05, 0.05)
    w = 2 * np.pi * f
    phi = MOUTHS_PER_BANK * h * HF_SPAN / BANK_AREA      # mouth fraction within the bank
    Zb = channel_z(w, L * 1e-3, h * 1e-3, HF_SPAN * 1e-3) / phi
    return Zb


def hf_alpha(f, rng=None):
    share = BANK_AREA / FACE
    Y = sum(share / hf_bank_alpha(f, L, rng) for L in HF_L)
    Z = 1 / Y
    R = (Z - RHOC) / (Z + RHOC)
    return 1 - np.abs(R) ** 2


def bass_alpha(f, rng=None):
    types = [(h * 1e-3, 40e-3, L * 1e-3) for h, L in zip(BASS_H, BASS_L)]
    D = 112e-3 if rng is None else 112e-3 + rng.uniform(-2e-3, 2e-3)
    return slit_absorption(f, [(h, s, L, 0.7962 / 4) for h, s, L in types], D=D, cell=48e-3)


def oct_mean(a, f, fc):
    m = (f >= fc / np.sqrt(2)) & (f <= fc * np.sqrt(2))
    return float(a[m].mean())


if __name__ == '__main__':
    rng = np.random.default_rng(0)
    f = np.geomspace(60, 4000, 400)
    print('HF channel bank (A, rough) octave means:')
    hf = hf_alpha(f)
    for fc in (500, 1000, 2000, 4000):
        print('  %5d Hz  %.2f' % (fc, oct_mean(hf, f, fc)))
    print('Bass slot face (A) octave means:')
    bs = bass_alpha(f)
    for fc in (63, 125, 250):
        print('  %5d Hz  %.2f' % (fc, oct_mean(bs, f, fc)))
    print('Robustness (40 draws):')
    hf_m, bs_m = [], []
    for _ in range(40):
        hf_m.append(np.mean([oct_mean(hf_alpha(f, rng), f, fc) for fc in (500, 1000, 2000, 4000)]))
        bs_m.append(np.mean([oct_mean(bass_alpha(f, rng), f, fc) for fc in (125, 250)]))
    print('  HF 500-4k mean: median %.2f, 5th pct %.2f' % (np.median(hf_m), np.percentile(hf_m, 5)))
    print('  bass 125-250 mean: median %.2f, 5th pct %.2f' % (np.median(bs_m), np.percentile(bs_m, 5)))
