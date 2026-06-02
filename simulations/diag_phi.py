"""Diagnose why phi(HU) on Ma_activation does not slow the cycle.
Compare HU=0 vs HU=1.0: period, Ma/Mb/E2F/EZH2 stats, and time Ma spends above a
FIXED (DMSO) threshold per cycle = true S-dwell.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from calibrate_ezh2_rewire import run, period

REW = dict(K_Ma_EZH2=0.35, K_Me_EZH2=0.35, w_Me=0.65,
           k_EZH2_deg=0.20, k_EZH2_translation=0.80)
KW = dict(ezh2_rewire_params=REW, basal=0.02, e2f_amp=3.5, t_end=300, n_pts=30000)


def go(hu):
    return run(shh=0.5, hu=hu, checkpoint=None, **KW)


def s_dwell(r, thr_ma, t_start=50):
    """Mean fraction of time per cycle that Ma is above a FIXED threshold."""
    m = r['time'] >= t_start
    ma = r['[Ma]'][m]
    return (ma > thr_ma).mean()


base = go(0.0)
thr_ma = 0.20 * base['[Ma]'][base['time'] >= 50].max()
print(f"fixed thr_ma (DMSO) = {thr_ma:.3f}\n")
print(f"{'HU':>5}{'period':>8}{'phi':>6}{'Ma_mean':>9}{'Ma_max':>8}"
      f"{'S_dwell':>9}{'E2F_mean':>9}{'EZH2_mean':>10}")
for hu in (0.0, 0.5, 1.0, 2.0):
    r = go(hu)
    m = r['time'] >= 50
    phi = 0.45 + 0.55/(1+(hu/1.0)**3)
    ma = r['[Ma]'][m]; e2f = r['[E2F]'][m]; ez = r['[EZH2]'][m]
    print(f"{hu:>5.1f}{period(r):>8.1f}{phi:>6.2f}{ma.mean():>9.3f}{ma.max():>8.3f}"
          f"{s_dwell(r,thr_ma):>9.3f}{e2f.mean():>9.2f}{ez.mean():>10.3f}")
