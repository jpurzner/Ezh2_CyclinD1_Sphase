"""
Re-calibrate the rewired EZH2 module after adding the switch-like Hill gate
(Ma^n/(K^n+Ma^n)). The sharper gate zeroes arrested transcription (restoring the GDC
reduction and sharpening the within-cycle gradient) but lowers integrated cycling
transcription -- so amplitude (E2F_amp) and the stress coupling (k_EZH2_stress) must
be retuned. This sweep reports ALL targets at once for each (n, E2F_amp, basal):

  invariants (HU=0): mean cycling EZH2 ~0.49 | GDC reduction 25-47% | G0/cyc < 1
  gradient   (HU=0): G1/S/G2 vs G0 -> target 1.16 / 1.32 / 1.48, monotone
  HU=1.0 boost     : EZH2 in S/G2 (high-CycA) -> target ~1.31, no arrest, period up

Run:  ./venv/bin/python simulations/sim_recalibrate_hill.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from calibrate_ezh2_rewire import run, period, mean_last, mean_ezh2_by_phase

BASE = dict(K_Ma_EZH2=0.35, k_EZH2_deg=0.20, k_EZH2_translation=0.80, k_EZH2_stress=0.16)
T_END, NPTS = 240, 24000


def go(n, e2f_amp, basal, hu=0.0, shh=0.5, gdc=0.0, k_stress=None):
    rew = dict(BASE, n_Ma_EZH2=n)
    if k_stress is not None:
        rew['k_EZH2_stress'] = k_stress
    return run(shh=shh, hu=hu, gdc=gdc, ezh2_rewire_params=rew,
               basal=basal, e2f_amp=e2f_amp, t_end=T_END, n_pts=NPTS)


def n_div(r, t_start=40):
    t = r['time']; m = t >= t_start; cb = r['[Cb]'][m]
    if cb.max() <= 0:
        return 0
    pk, _ = find_peaks(cb, prominence=0.05 * cb.max(), distance=int(0.05 * len(cb)))
    return len(pk)


def sg2_ezh2(r, thr_ma, t_start=40):
    t = r['time']; m = t >= t_start
    ma = r['[Ma]'][m]; ez = r['[EZH2]'][m]
    sel = ma > thr_ma
    return float(np.mean(ez[sel])) if sel.any() else np.nan


def report(n, e2f_amp, basal):
    cyc = go(n, e2f_amp, basal)
    g0 = go(n, e2f_amp, basal, shh=0.0)
    gdc = go(n, e2f_amp, basal, gdc=1.0)
    ez_cyc = mean_last(cyc, 'EZH2'); ez_g0 = mean_last(g0, 'EZH2')
    ez_gdc = mean_last(gdc, 'EZH2')
    by, thr = mean_ezh2_by_phase(cyc)
    rG1, rS, rG2 = by['G1']/ez_g0, by['S']/ez_g0, by['G2M']/ez_g0
    gdc_red = 100 * (1 - ez_gdc/ez_cyc)
    mono = ez_g0 < by['G1'] < by['S'] < by['G2M']

    hu1 = go(n, e2f_amp, basal, hu=1.0)
    b = sg2_ezh2(hu1, thr['thr_ma']) / sg2_ezh2(cyc, thr['thr_ma'])
    per0, per1, div1 = period(cyc), period(hu1), n_div(hu1)

    print(f"n={n} E2F_amp={e2f_amp} basal={basal}:")
    print(f"   mean_cyc={ez_cyc:.3f}(~.49)  GDC_red={gdc_red:+.0f}%(25-47)  "
          f"G0/cyc={ez_g0/ez_cyc:.2f}(<1)")
    print(f"   gradient G1={rG1:.2f} S={rS:.2f} G2={rG2:.2f} (1.16/1.32/1.48) "
          f"mono={'Y' if mono else 'N'}")
    print(f"   HU=1 boost={b:.2f}(~1.31) per {per0:.1f}->{per1:.1f}h div={div1} "
          f"{'ARREST' if div1<2 else ''}")


if __name__ == "__main__":
    print("=" * 78)
    print("RE-CALIBRATE EZH2 REWIRE WITH HILL GATE")
    print("=" * 78)
    print(f"base: {BASE}\n")
    for n in (2, 4):
        for e2f_amp in (3.5, 5.0, 7.0):
            for basal in (0.02,):
                report(n, e2f_amp, basal)
                print()
