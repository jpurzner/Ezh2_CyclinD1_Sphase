"""HU dose-response under the algebraic phi(HU) S-progression brake (de-hacked v43).

phi(HU) scales the CycA/CDK2 (Ma) activation flux -> slows S-progression -> cells dwell
longer in the high-CycA S window. EZH2 rise is EMERGENT (no HU term in EZH2 eqn).

Validation targets (from the spec):
  - graded slowing, NOT arrest (period up, divisions still happen at all doses)
  - elevated S-phase fraction with dose
  - EZH2 in S/G2 rises with HU; target ~1.31x at HU=1.0 (emergent)

Run: ./venv/bin/python simulations/hu_phi_response.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from calibrate_ezh2_rewire import run, period, mean_last, mean_ezh2_by_phase

REW = dict(K_Ma_EZH2=0.35, K_Me_EZH2=0.35, w_Me=0.65,
           k_EZH2_deg=0.20, k_EZH2_translation=0.80)
BASAL, E2F_AMP = 0.02, 3.5
T_END, NPTS = 300, 30000


def go(hu, phi=None):
    cp = phi  # phi overrides dict or None
    return run(shh=0.5, hu=hu, checkpoint=cp, ezh2_rewire_params=REW,
               basal=BASAL, e2f_amp=E2F_AMP, t_end=T_END, n_pts=NPTS)


def n_div(r, t_start=50):
    m = r['time'] >= t_start; cb = r['[Cb]'][m]
    if cb.max() <= 0: return 0
    pk, _ = find_peaks(cb, prominence=0.05*cb.max(), distance=int(0.04*len(cb)))
    return len(pk)


def s_frac(r, t_start=50):
    m = r['time'] >= t_start
    ma = r['[Ma]'][m]; cb = r['[Cb]'][m]
    thr_ma = 0.20*ma.max(); thr_cb = 0.25*cb.max()
    g2 = cb >= thr_cb
    s = (~g2) & (ma >= thr_ma)
    return s.mean(), (~g2 & (ma < thr_ma)).mean(), g2.mean()


def sg2_ezh2(r, thr_ma, t_start=50):
    m = r['time'] >= t_start
    ma = r['[Ma]'][m]; ez = r['[EZH2]'][m]
    sel = ma > thr_ma
    return float(np.mean(ez[sel])) if sel.any() else np.nan


if __name__ == "__main__":
    print("="*78)
    print("HU DOSE-RESPONSE via phi(HU) S-progression brake (emergent EZH2)")
    print("="*78)
    print(f"phi defaults | REW={REW}\n")
    # fix the S/G2 threshold from DMSO and reuse for all HU (one gate, all conditions)
    base = go(0.0)
    _, thr = mean_ezh2_by_phase(base)
    thr_ma = thr['thr_ma']
    sg2_0 = sg2_ezh2(base, thr_ma)
    print(f"{'HU':>5}{'period':>9}{'div':>5}{'G1f':>7}{'Sf':>7}{'G2f':>7}"
          f"{'meanEZH2':>10}{'S/G2 EZH2':>11}{'boost':>8}")
    for hu in (0.0, 0.5, 1.0, 1.5, 2.0):
        r = go(hu)
        per = period(r); d = n_div(r)
        sf, g1f, g2f = s_frac(r)
        mez = mean_last(r, 'EZH2')
        sg2 = sg2_ezh2(r, thr_ma)
        boost = sg2/sg2_0
        arr = '  <-- ARREST' if d < 2 else ''
        print(f"{hu:>5.1f}{per:>9.1f}{d:>5d}{g1f:>7.2f}{sf:>7.2f}{g2f:>7.2f}"
              f"{mez:>10.3f}{sg2:>11.3f}{boost:>8.2f}{arr}")
