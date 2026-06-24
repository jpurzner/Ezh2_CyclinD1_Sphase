"""
Diagnose the two failing v43 invariants (HU=0):
  A3  HHi (HHi) should LOWER EZH2 (v42: 25-47% reduction); currently flat.
  B1  within-cycle gradient should be monotone G0<G1<S<G2; currently G1,S dip below G0.

Both hinge on what GDC and cell-cycle phase do to the EZH2 transcription DRIVERS
(E2F and Ma) and to EZH2 protein. Print the drivers side-by-side so we can see which
coupling is too weak.

Run:  ./venv/bin/python simulations/diag_gdc_gradient.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from calibrate_ezh2_rewire import run, period, mean_last, phase_call

REW = dict(K_Ma_EZH2=0.35, k_EZH2_deg=0.20, k_EZH2_translation=0.80, k_EZH2_stress=0.16)
BASAL, E2F_AMP = 0.02, 3.5
KW = dict(ezh2_rewire_params=REW, basal=BASAL, e2f_amp=E2F_AMP, t_end=240, n_pts=24000)


def stats(r, t_start=40):
    t = r['time']; m = t >= t_start
    out = {}
    for sp in ('E2F', 'Ma', 'EZH2', 'EZH2_mRNA', 'pRBpp', 'CycD1'):
        key = f'[{sp}]'
        if key in r.colnames:
            v = r[key][m]
            out[sp] = (float(v.mean()), float(v.max()), float(v.min()))
    return out


def show(label, r):
    s = stats(r)
    print(f"\n{label}: period={period(r):.1f}h  EZH2_last={mean_last(r,'EZH2'):.3f}")
    for sp in ('E2F', 'Ma', 'EZH2_mRNA', 'EZH2', 'CycD1'):
        if sp in s:
            mean, mx, mn = s[sp]
            print(f"    {sp:10s} mean={mean:.3f}  max={mx:.3f}  min={mn:.3f}")
    # the E2F Hill factor that gates EZH2 transcription (K_E2F_EZH2 = 5.0 in v42)
    e2f_mean = s['E2F'][0] if 'E2F' in s else float('nan')
    print(f"    -> E2F/(5+E2F) at mean E2F = {e2f_mean/(5.0+e2f_mean):.3f} "
          f"(saturates near 1 -> GDC's E2F drop won't propagate)")


if __name__ == "__main__":
    print("=" * 78)
    print("DIAGNOSE GDC + within-cycle gradient (HU=0)")
    print("=" * 78)
    print(f"rewire: {REW}  basal={BASAL} e2f_amp={E2F_AMP}")

    cyc = run(shh=0.5, **KW)
    hhi = run(shh=0.5, hhi=1.0, **KW)
    g0 = run(shh=0.0, **KW)

    show("CYCLING (SHH)", cyc)
    show("HHi (HHi)", hhi)
    show("G0 (no SHH)", g0)

    # within-cycle gradient
    by, thr = (lambda x: x)(__import__('calibrate_ezh2_rewire').mean_ezh2_by_phase(cyc))
    ez_g0 = mean_last(g0, 'EZH2')
    print("\nWITHIN-CYCLE EZH2 (means by phase, ratio vs G0):")
    print(f"    G0={ez_g0:.3f} (1.00)  G1={by['G1']:.3f} ({by['G1']/ez_g0:.2f})  "
          f"S={by['S']:.3f} ({by['S']/ez_g0:.2f})  G2={by['G2M']:.3f} ({by['G2M']/ez_g0:.2f})")
    print(f"    target ratios: G1=1.16 S=1.32 G2=1.48")
    print(f"    thr_ma={thr['thr_ma']:.3f} thr_cb={thr['thr_cb']:.3f}")
