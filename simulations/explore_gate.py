"""Explore EZH2 transcription gate drivers (Me vs Ma) to fix the within-cycle
gradient (S must exceed G1) AND the GDC reduction simultaneously.

Builds v43 (current defaults), then string-replaces the EZH2_mRNA_synthesis gate
factor with candidate forms, runs CYCLING / GDC / G0 / HU=1, and reports:
  mean_cyc(~.49)  G0/cyc(~.6)  GDC_red(25-47%)  gradient G1/S/G2(1.16/1.32/1.48,mono)
  HU=1 S/G2 boost(~1.31)

Run: ./venv/bin/python simulations/explore_gate.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from src.build_model_v43_checkpoint import build_model_v43
from calibrate_ezh2_rewire import (period, mean_last, mean_ezh2_by_phase)
import tellurium as te

BASAL, E2F_AMP = 0.02, 3.5
T_END, NPTS = 240, 24000

# the v43 default gate factor that we will replace
GATE_OLD = "Ma^n_Ma_EZH2 / (K_Ma_EZH2^n_Ma_EZH2 + Ma^n_Ma_EZH2)"

# extra constants we may need (Me half-sat); declare next to K_Ma_EZH2
EXTRA_CONSTS = "K_Me_EZH2 = {K_Me};\nw_Me = {w_Me};\n"


def build(gate, K_Me=0.30, w_Me=0.5):
    m = build_model_v43()
    if GATE_OLD not in m:
        raise RuntimeError("gate factor not found")
    m = m.replace(GATE_OLD, gate)
    m = m.replace("K_Ma_EZH2 = 0.35;",
                  "K_Ma_EZH2 = 0.35;\n" + EXTRA_CONSTS.format(K_Me=K_Me, w_Me=w_Me))
    return m


def simulate(m, shh=0.5, hhi=0.0, hu=0.0):
    rr = te.loada(m)
    rr['SHH'] = shh; rr['HHi'] = hhi; rr['HU'] = hu
    rr['k_EZH2_mRNA_synth_basal'] = BASAL
    rr['k_EZH2_mRNA_synth_E2F'] = E2F_AMP
    if shh == 0.0:
        rr['SHH_Ptch'] = 0.0; rr['Ptch1_free'] = 1.0; rr['Smo_active'] = 0.01
    return rr.simulate(0, T_END, NPTS)


def sg2_ezh2(r, thr_ma, t_start=40):
    t = r['time']; mask = t >= t_start
    ma = r['[Ma]'][mask]; ez = r['[EZH2]'][mask]
    sel = ma > thr_ma
    return float(np.mean(ez[sel])) if sel.any() else np.nan


def n_div(r, t_start=40):
    t = r['time']; mask = t >= t_start; cb = r['[Cb]'][mask]
    if cb.max() <= 0:
        return 0
    pk, _ = find_peaks(cb, prominence=0.05*cb.max(), distance=int(0.05*len(cb)))
    return len(pk)


def report(label, gate, K_Me=0.30, w_Me=0.5):
    m = build(gate, K_Me=K_Me, w_Me=w_Me)
    cyc = simulate(m); g0 = simulate(m, shh=0.0); hhi = simulate(m, hhi=1.0)
    hu1 = simulate(m, hu=1.0)
    ez_cyc = mean_last(cyc, 'EZH2'); ez_g0 = mean_last(g0, 'EZH2')
    ez_gdc = mean_last(hhi, 'EZH2')
    by, thr = mean_ezh2_by_phase(cyc)
    rG1, rS, rG2 = by['G1']/ez_g0, by['S']/ez_g0, by['G2M']/ez_g0
    gdc_red = 100*(1 - ez_gdc/ez_cyc)
    mono = ez_g0 < by['G1'] < by['S'] < by['G2M']
    b = sg2_ezh2(hu1, thr['thr_ma'])/sg2_ezh2(cyc, thr['thr_ma'])
    print(f"[{label}]")
    print(f"   mean_cyc={ez_cyc:.3f}(~.49) G0/cyc={ez_g0/ez_cyc:.2f}(~.6) "
          f"GDC_red={gdc_red:+.0f}%(25-47)")
    print(f"   grad G1={rG1:.2f} S={rS:.2f} G2={rG2:.2f} (1.16/1.32/1.48) "
          f"mono={'Y' if mono else 'N'}  HUboost={b:.2f}(~1.31) div={n_div(hu1)}")
    print()


if __name__ == "__main__":
    print("="*70)
    print("EXPLORE EZH2 GATE DRIVERS")
    print("="*70)
    # first: Me vs Ma magnitude/timing
    m = build("Ma^n_Ma_EZH2 / (K_Ma_EZH2^n_Ma_EZH2 + Ma^n_Ma_EZH2)")
    cyc = simulate(m)
    for sp in ('Me', 'Ma'):
        v = cyc[f'[{sp}]']
        print(f"{sp}: max={v.max():.3f} mean={v.mean():.3f}")
    print()

    # A. baseline Ma n=1
    report("A Ma n1", "Ma^n_Ma_EZH2 / (K_Ma_EZH2^n_Ma_EZH2 + Ma^n_Ma_EZH2)")
    # B. Ma n=2 (Hill)
    report("B Ma n2", "Ma^2 / (K_Ma_EZH2^2 + Ma^2)")
    # C. Me OR Ma (sum of two n=1 gates, scaled by w_Me/(1-w_Me))
    report("C Me+Ma sum", "(w_Me*Me/(K_Me_EZH2+Me) + (1-w_Me)*Ma/(K_Ma_EZH2+Ma))")
    # D. Me OR Ma, both Hill n=2
    report("D Me+Ma n2", "(w_Me*Me^2/(K_Me_EZH2^2+Me^2) + (1-w_Me)*Ma^2/(K_Ma_EZH2^2+Ma^2))")
