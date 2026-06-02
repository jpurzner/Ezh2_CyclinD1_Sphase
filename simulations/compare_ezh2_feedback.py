"""Test the EZH2 -> CycD1 feedback under the new S-phase-controlled EZH2 expression.

The loop:  S-phase (Ma) -> EZH2 (Me+Ma gate) -> represses CycD1 transcription
           -> CycD1/CDK4-6 (Md) -> pRB/E2F -> cell-cycle entry -> S-phase.

EZH2i toggles ONLY the EZH2->CycD1 repression arm (EZH2*(1-EZH2i) in CycD1_transcription):
    EZH2i=0  -> feedback PRESENT   (EZH2 represses CycD1)
    EZH2i=1  -> feedback ABSENT    (repression factor = 1; EZH2 protein still accumulates)

Compared under DMSO (HU=0) and HU=1.0 (longer S-phase -> more EZH2 -> stronger
repression, if the feedback is intact). Readouts: CycD1 protein [Cd], active
CycD1/CDK4-6 [Md], EZH2, cycle period, division count, G1 time-fraction.

Run: ./venv/bin/python simulations/compare_ezh2_feedback.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from calibrate_ezh2_rewire import run, period, mean_last, mean_ezh2_by_phase

REW = dict(K_Ma_EZH2=0.35, K_Me_EZH2=0.35, w_Me=0.65,
           k_EZH2_deg=0.20, k_EZH2_translation=0.80, k_EZH2_stress=0.16)
BASAL, E2F_AMP = 0.02, 3.5
T_END, NPTS = 300, 30000


def go(hu, ezh2i):
    return run(shh=0.5, hu=hu, ezh2i=ezh2i, ezh2_rewire_params=REW,
               basal=BASAL, e2f_amp=E2F_AMP, t_end=T_END, n_pts=NPTS)


def n_div(r, t_start=50):
    m = r['time'] >= t_start; cb = r['[Cb]'][m]
    if cb.max() <= 0: return 0
    pk, _ = find_peaks(cb, prominence=0.05*cb.max(), distance=int(0.04*len(cb)))
    return len(pk)


def phase_fracs(r, t_start=50):
    """Fraction of time in G1 / S / G2M using fixed fraction-of-max thresholds."""
    m = r['time'] >= t_start
    ma = r['[Ma]'][m]; cb = r['[Cb]'][m]
    thr_ma = 0.20*ma.max(); thr_cb = 0.25*cb.max()
    g2 = cb >= thr_cb
    s = (~g2) & (ma >= thr_ma)
    g1 = (~g2) & (ma < thr_ma)
    n = len(ma)
    return g1.sum()/n, s.sum()/n, g2.sum()/n


def mean_win(r, sp, t_start=50):
    m = r['time'] >= t_start
    return float(np.mean(r[f'[{sp}]'][m]))


def report(hu):
    on = go(hu, 0.0)    # feedback present
    off = go(hu, 1.0)   # feedback absent
    print(f"\n--- HU={hu} ---")
    print(f"{'metric':<22}{'feedback ON':>14}{'feedback OFF':>14}{'OFF/ON':>10}")
    rows = [
        ('CycD1 protein [Cd]', mean_win(on,'Cd'),  mean_win(off,'Cd')),
        ('active CycD [Md]',   mean_win(on,'Md'),  mean_win(off,'Md')),
        ('CycD1 mRNA',         mean_win(on,'Cd_mRNA'), mean_win(off,'Cd_mRNA')),
        ('EZH2',               mean_win(on,'EZH2'),mean_win(off,'EZH2')),
        ('E2F',                mean_win(on,'E2F'), mean_win(off,'E2F')),
    ]
    for name, a, b in rows:
        ratio = b/a if a else float('nan')
        print(f"{name:<22}{a:>14.3f}{b:>14.3f}{ratio:>10.2f}")
    per_on, per_off = period(on), period(off)
    d_on, d_off = n_div(on), n_div(off)
    g1_on, s_on, g2_on = phase_fracs(on)
    g1_off, s_off, g2_off = phase_fracs(off)
    print(f"{'period (h)':<22}{per_on:>14.1f}{per_off:>14.1f}")
    print(f"{'divisions/250h':<22}{d_on:>14d}{d_off:>14d}")
    print(f"{'G1 frac':<22}{g1_on:>14.2f}{g1_off:>14.2f}")
    print(f"{'S frac':<22}{s_on:>14.2f}{s_off:>14.2f}")
    print(f"{'G2M frac':<22}{g2_on:>14.2f}{g2_off:>14.2f}")


if __name__ == "__main__":
    print("="*64)
    print("EZH2 -> CycD1 FEEDBACK: presence vs absence (S-phase-controlled EZH2)")
    print("="*64)
    print(f"rewire: {REW}")
    for hu in (0.0, 1.0):
        report(hu)
