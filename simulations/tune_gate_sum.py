"""Tune the Me+Ma SUM gate (candidate C) for EZH2 transcription.

Gate = w_Me*Me/(K_Me+Me) + (1-w_Me)*Ma/(K_Ma+Ma), K_Ma fixed 0.35.
Raising w_Me shifts weight to the early Me (S-onset) spike -> lifts G1/S, softens
the G2 overshoot. Lowering basal drops G0 toward the 0.6 target (but inflates the
G0-normalized gradient ratios). Map the (w_Me, K_Me, basal) tradeoff against all
targets at once.

Run: ./venv/bin/python simulations/tune_gate_sum.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from src.build_model_v43_checkpoint import build_model_v43
from calibrate_ezh2_rewire import period, mean_last, mean_ezh2_by_phase
import tellurium as te

T_END, NPTS = 240, 24000
GATE_OLD = "Ma^n_Ma_EZH2 / (K_Ma_EZH2^n_Ma_EZH2 + Ma^n_Ma_EZH2)"
GATE_NEW = "(w_Me*Me/(K_Me_EZH2+Me) + (1-w_Me)*Ma/(K_Ma_EZH2+Ma))"


def build(K_Me, w_Me):
    m = build_model_v43()
    m = m.replace(GATE_OLD, GATE_NEW)
    m = m.replace("K_Ma_EZH2 = 0.35;",
                  f"K_Ma_EZH2 = 0.35;\nK_Me_EZH2 = {K_Me};\nw_Me = {w_Me};\n")
    return m


def simulate(m, basal, shh=0.5, hhi=0.0, hu=0.0, e2f=3.5):
    rr = te.loada(m)
    rr['SHH'] = shh; rr['HHi'] = hhi; rr['HU'] = hu
    rr['k_EZH2_mRNA_synth_basal'] = basal
    rr['k_EZH2_mRNA_synth_E2F'] = e2f
    if shh == 0.0:
        rr['SHH_Ptch'] = 0.0; rr['Ptch1_free'] = 1.0; rr['Smo_active'] = 0.01
    return rr.simulate(0, T_END, NPTS)


def sg2(r, thr_ma, t_start=40):
    mask = r['time'] >= t_start
    ma = r['[Ma]'][mask]; ez = r['[EZH2]'][mask]
    sel = ma > thr_ma
    return float(np.mean(ez[sel])) if sel.any() else np.nan


def ndiv(r, t_start=40):
    mask = r['time'] >= t_start; cb = r['[Cb]'][mask]
    if cb.max() <= 0: return 0
    pk, _ = find_peaks(cb, prominence=0.05*cb.max(), distance=int(0.05*len(cb)))
    return len(pk)


def row(K_Me, w_Me, basal, e2f=3.5):
    m = build(K_Me, w_Me)
    cyc = simulate(m, basal, e2f=e2f); g0 = simulate(m, basal, shh=0.0, e2f=e2f)
    hhi = simulate(m, basal, hhi=1.0, e2f=e2f); hu1 = simulate(m, basal, hu=1.0, e2f=e2f)
    ez_cyc = mean_last(cyc,'EZH2'); ez_g0 = mean_last(g0,'EZH2'); ez_gdc = mean_last(hhi,'EZH2')
    by, thr = mean_ezh2_by_phase(cyc)
    rG1, rS, rG2 = by['G1']/ez_g0, by['S']/ez_g0, by['G2M']/ez_g0
    gdc_red = 100*(1-ez_gdc/ez_cyc)
    mono = ez_g0 < by['G1'] < by['S'] < by['G2M']
    b = sg2(hu1, thr['thr_ma'])/sg2(cyc, thr['thr_ma'])
    print(f"K_Me={K_Me} w={w_Me} basal={basal} e2f={e2f}: "
          f"cyc={ez_cyc:.3f} G0/cyc={ez_g0/ez_cyc:.2f} GDC={gdc_red:+.0f}% "
          f"| G1={rG1:.2f} S={rS:.2f} G2={rG2:.2f} mono={'Y' if mono else 'N'} "
          f"| HU={b:.2f} div={ndiv(hu1)}")


if __name__ == "__main__":
    print("="*90)
    print("TUNE Me+Ma SUM gate  (targets: cyc~.49 G0/cyc~.6 GDC 25-47 | "
          "G1/S/G2 1.16/1.32/1.48 mono | HU~1.31)")
    print("="*90)
    for w_Me in (0.5, 0.65, 0.8):
        for K_Me in (0.20, 0.35):
            for basal in (0.02, 0.013):
                row(K_Me, w_Me, basal)
        print()
