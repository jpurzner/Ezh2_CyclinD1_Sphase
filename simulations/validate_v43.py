"""
Validation suite for model v43 (intra-S checkpoint + EZH2->Ma transcription rewire).

v43 = published v42  +  (1) EZH2 transcription rewired from the pRBpp spike to a
weighted Me+Ma (CycE/CycA-CDK2) gate with faster protein turnover, and (2) an
algebraic, phase-gated HU effect: phi(HU) slows the cell-cycle CLOCK (eps)
SPECIFICALLY within the replication window (sphase = Ma-high), prolonging S/G2.
There is NO HU term in the EZH2 equation -- the HU -> EZH2 response is EMERGENT
(longer dwell in the high-synthesis S window). At HU=0, phi=1 -> eps=eps0 in every
phase -> the model is numerically identical to v42 (+ rewire), so the published
validation is preserved exactly.

What this checks (each prints PASS/FAIL):

  A. v42-PRESERVING INVARIANTS at HU=0 (the rewire must not regress published results)
       A1 mean cycling EZH2     ~0.49   (sets EZH2i->CycD1 repression strength)
       A2 EZH2 G0/cycling ratio  <1     (quiescent EZH2 below cycling)
       A3 GDC0449 (HHi) lowers EZH2     (Hedgehog dependence retained)
  B. NEW within-cycle gradient unlocked by the rewire (v42 could not show this)
       G0 < G1 < S < G2  with G2/G0 approaching the measured 1.48
  C. HU DOSE RESPONSE (the v43 hypothesis: longer/stressed S -> more EZH2, EMERGENT)
       C1 HU=0 is inert to the phi params (checkpoint silent at baseline)
       C2 10 uM HU (HU=1.0) -> ~1.31x EZH2 in the S/G2 (high-CycA) compartment
       C3 cells still divide (NO arrest) at every HU dose
       C4 period lengthens with HU (clock slowed in a prolonged S)
       C5 S-phase dwell lengthens with HU (the emergent driver of the EZH2 rise)
  D. STRUCTURAL: builder leaves the v42 core intact (string-surgery assertions pass)

Run:  ./venv/bin/python simulations/validate_v43.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from calibrate_ezh2_rewire import (
    run, period, mean_last, mean_ezh2_by_phase, s_phase_duration)
from src.build_model_v43_checkpoint import (
    build_model_v43, EZH2_REWIRE_DEFAULTS, HU_PHI_DEFAULTS)

# Canonical Phase-3 calibrated configuration ---------------------------------
REW = dict(K_Ma_EZH2=0.35, K_Me_EZH2=0.35, w_Me=0.65,
           k_EZH2_deg=0.20, k_EZH2_translation=0.80)
BASAL, E2F_AMP = 0.02, 3.5
T_END = 300
NPTS = T_END * 100

_n_pass = 0
_n_fail = 0


def check(name, ok, detail=""):
    global _n_pass, _n_fail
    tag = "PASS" if ok else "FAIL"
    if ok:
        _n_pass += 1
    else:
        _n_fail += 1
    print(f"  [{tag}] {name}{('  -- ' + detail) if detail else ''}")


def go(hu, shh=0.5, gdc=0.0, checkpoint=None):
    return run(shh=shh, hu=hu, gdc=gdc, ezh2_rewire_params=REW, checkpoint=checkpoint,
               basal=BASAL, e2f_amp=E2F_AMP, t_end=T_END, n_pts=NPTS)


def n_div(r, t_start=40):
    t = r['time']; m = t >= t_start; cb = r['[Cb]'][m]
    if cb.max() <= 0:
        return 0
    pk, _ = find_peaks(cb, prominence=0.05 * cb.max(), distance=int(0.05 * len(cb)))
    return len(pk)


def s_g2_ezh2(r, thr, t_start=40):
    """Mean EZH2 over the high-CycA (Ma>thr) S/G2 compartment = 'S-phase cells'."""
    t = r['time']; m = t >= t_start
    ma = r['[Ma]'][m]; ez = r['[EZH2]'][m]
    sel = ma > thr
    return float(np.mean(ez[sel])) if sel.any() else np.nan


def s_dwell(r, thr, t_start=40):
    """Fraction of time Ma is above the FIXED DMSO threshold (S/G2 occupancy)."""
    t = r['time']; m = t >= t_start
    return float((r['[Ma]'][m] > thr).mean())


if __name__ == "__main__":
    print("=" * 78)
    print("v43 VALIDATION SUITE")
    print("=" * 78)
    print(f"calibrated rewire: {REW}")
    print(f"phi(HU) phase-gated clock: {HU_PHI_DEFAULTS}")
    print(f"runtime tx: basal={BASAL} E2F_amp={E2F_AMP}\n")

    # ---- reference runs (HU=0) ----
    cyc = go(0.0)                      # cycling (GNP + SHH)
    g0 = go(0.0, shh=0.0)              # quiescent reference
    gdc = go(0.0, gdc=1.0)            # Hedgehog inhibitor
    ez_cyc = mean_last(cyc, 'EZH2')
    ez_g0 = mean_last(g0, 'EZH2')
    ez_gdc = mean_last(gdc, 'EZH2')
    by, thr = mean_ezh2_by_phase(cyc)
    thr_ma = thr['thr_ma']

    print("A. v42-PRESERVING INVARIANTS (HU=0)")
    check("A1 mean cycling EZH2 ~0.49", 0.40 <= ez_cyc <= 0.58, f"{ez_cyc:.3f}")
    check("A2 G0/cycling < 1", ez_g0 < ez_cyc, f"G0={ez_g0:.3f} cyc={ez_cyc:.3f} "
          f"ratio={ez_g0/ez_cyc:.2f}")
    check("A3 GDC0449 lowers EZH2", ez_gdc < ez_cyc,
          f"GDC={ez_gdc:.3f} ({100*(1-ez_gdc/ez_cyc):.0f}% reduction)")

    print("\nB. WITHIN-CYCLE GRADIENT (HU=0)")
    rG1, rS, rG2 = by['G1']/ez_g0, by['S']/ez_g0, by['G2M']/ez_g0
    check("B1 monotone G0<G1<S<G2", ez_g0 < by['G1'] < by['S'] < by['G2M'],
          f"G1={rG1:.2f} S={rS:.2f} G2={rG2:.2f} (tgt 1.16/1.32/1.48)")
    check("B2 G2/G0 rises toward 1.48", rG2 >= 1.30, f"G2/G0={rG2:.2f}")

    print("\nC. HU DOSE RESPONSE (emergent: phi(HU) slows the clock only in S)")
    # C1: at HU=0, phi=1 regardless of phi params -> changing them must not move EZH2.
    cyc_altphi = go(0.0, checkpoint=dict(phi_min=0.10, K_HU_phi=0.3))
    d = abs(mean_last(cyc_altphi, 'EZH2') - ez_cyc)
    check("C1 HU=0 inert to phi params (checkpoint silent at baseline)", d < 1e-6,
          f"|delta|={d:.2e}")

    sg2_dmso = s_g2_ezh2(cyc, thr_ma)
    per_dmso, div_dmso, dwell_dmso = period(cyc), n_div(cyc), s_dwell(cyc, thr_ma)
    print(f"     DMSO: per={per_dmso:.1f}h div={div_dmso} S_dwell={dwell_dmso:.3f} "
          f"S/G2_EZH2={sg2_dmso:.3f}")
    boosts = {}
    for hu in (0.5, 1.0, 2.0):
        r = go(hu)
        b = s_g2_ezh2(r, thr_ma) / sg2_dmso
        boosts[hu] = (b, period(r), n_div(r), s_dwell(r, thr_ma))
        print(f"     HU={hu}: per={boosts[hu][1]:.1f}h div={boosts[hu][2]} "
              f"S_dwell={boosts[hu][3]:.3f} S/G2_EZH2={s_g2_ezh2(r, thr_ma):.3f}  "
              f"boost={b:.2f}")

    b1, per1, div1, dwell1 = boosts[1.0]
    check("C2 HU=1.0 -> ~1.31x EZH2 in S/G2", 1.25 <= b1 <= 1.40, f"boost={b1:.2f}")
    check("C3 no arrest at any HU dose (div>=2)",
          all(v[2] >= 2 for v in boosts.values()),
          "div=" + ",".join(f"{k}:{v[2]}" for k, v in boosts.items()))
    check("C4 period lengthens with HU", per1 > per_dmso,
          f"DMSO={per_dmso:.1f}h HU1={per1:.1f}h")
    check("C5 S-dwell lengthens with HU", dwell1 > dwell_dmso,
          f"DMSO={dwell_dmso:.3f} HU1={dwell1:.3f}")

    print("\nD. STRUCTURAL (v42 core intact)")
    m43 = build_model_v43(ezh2_rewire_params=REW)
    core_ok = ("model ezh2_cyclind1_v43_checkpoint" in m43
               and "phi_HU := phi_min" in m43
               and "sphase := Ma^n_S" in m43
               and "eps := eps0 * (1 - sphase * (1 - phi_HU))" in m43
               and "(w_Me * Me / (K_Me_EZH2 + Me) + " in m43
               and "k_EZH2_stress" not in m43
               and "Pa * phi_HU * eps" not in m43
               and "eps = 150;" not in m43)
    check("D1 v43 builds with phase-gated clock + Me+Ma rewire (no stress hack)", core_ok)
    # rewire OFF must reproduce exact v42 EZH2 wiring (pRBpp)
    m43_norw = build_model_v43(rewire_ezh2=False)
    check("D2 rewire_ezh2=False recovers v42 pRBpp wiring",
          "pRBpp / (K_pRBpp_EZH2 + pRBpp)" in m43_norw
          and "k_EZH2_stress" not in m43_norw)

    print("\n" + "=" * 78)
    print(f"RESULT: {_n_pass} passed, {_n_fail} failed")
    print("=" * 78)
    sys.exit(1 if _n_fail else 0)
