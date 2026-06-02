"""
Test the Cdc20/CycA-stability checkpoint lever (v43): does CHK1 -> slowed CycA
destruction genuinely PROLONG the high-Ma window and raise EZH2 (real, not artifact)?

The Cdc25 levers fail in the Gerard-Goldbeter core (Pe arrests; Pb only inflates free
CycB and fooled a max-scaled phase detector). Here CHK1 slows Cdc20-mediated CycA
destruction, holding Ma high longer. We measure the boost with ARTIFACT-PROOF metrics:

  * phase thresholds are derived ONCE from DMSO and reused for HU (one gate, both
    conditions -- the experimental analog), so CycB inflation cannot relabel cells.
  * whole-cycle mean EZH2 and EZH2 over a FIXED absolute Ma window are reported as
    ground truth: if these rise under HU, the boost is real.

Targets: 10 uM HU -> 1.31x EZH2 in S; S-phase prolonged but not abolished; period up.

Run:  ./venv/bin/python simulations/sim_cdc20_lever.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from calibrate_ezh2_rewire import run, period, mean_last

REW = dict(k_EZH2_deg=0.20, K_Ma_EZH2=0.35, k_EZH2_translation=0.80)
BASAL, E2F_AMP = 0.02, 3.5
OFF = 1.0e6


def go(hu, Ki_Cdc20, t_end=240):
    chk = dict(Ki_CHK1_Pe=OFF, Ki_CHK1_Pb=OFF, Ki_CHK1_Cdc20=Ki_Cdc20)
    return run(shh=0.5, hu=hu, ezh2_rewire_params=REW, basal=BASAL, e2f_amp=E2F_AMP,
               checkpoint=chk, t_end=t_end, n_pts=int(t_end * 100))


def fixed_thresholds(dmso, t_start=40):
    """Derive S/G2 gates from the DMSO run; reuse for every condition."""
    t = dmso['time']; m = t >= t_start
    ma = dmso['[Ma]'][m]; cb = dmso['[Cb]'][m]
    return dict(thr_ma=0.20 * ma.max(), thr_cb=0.25 * cb.max())


def metrics(r, thr, t_start=40):
    t = r['time']; m = t >= t_start
    ma = r['[Ma]'][m]; cb = r['[Cb]'][m]; ez = r['[EZH2]'][m]
    s_sel = (ma > thr['thr_ma']) & (cb <= thr['thr_cb'])      # S: replicating, pre-mitotic
    sg2_sel = (ma > thr['thr_ma'])                             # broad high-CycA window
    # contiguous high-Ma duration per cycle
    durs, on, start = [], False, None
    tt = t[m]
    for i in range(len(tt)):
        if sg2_sel[i] and not on:
            on, start = True, tt[i]
        elif not sg2_sel[i] and on:
            on = False; durs.append(tt[i] - start)
    return dict(
        ezh2_S=float(np.mean(ez[s_sel])) if s_sel.any() else np.nan,
        ezh2_maHigh=float(np.mean(ez[sg2_sel])) if sg2_sel.any() else np.nan,
        ezh2_cycle=float(np.mean(ez)),
        maHigh_dur=float(np.mean(durs)) if durs else 0.0,
        maHigh_frac=float(sg2_sel.mean()),
    )


if __name__ == "__main__":
    print("=" * 78)
    print("Cdc20/CycA-STABILITY LEVER  (Pe & Pb OFF; fixed DMSO thresholds)")
    print("=" * 78)
    print(f"baseline rewire: {REW}  basal={BASAL} E2F_amp={E2F_AMP}\n")

    for Ki in (OFF, 0.50, 0.25, 0.12):
        tag = "OFF" if Ki >= 1e5 else f"{Ki}"
        print(f"\n### Ki_CHK1_Cdc20 = {tag} ###")
        dmso = go(0.0, Ki)
        thr = fixed_thresholds(dmso)
        md = metrics(dmso, thr)
        print(f"  thresholds (from DMSO): thr_ma={thr['thr_ma']:.3f} thr_cb={thr['thr_cb']:.3f}")
        print(f"  DMSO : per={period(dmso):.1f}h  maHi_dur={md['maHigh_dur']:.2f}h "
              f"maHi_frac={md['maHigh_frac']:.3f}  EZH2: S={md['ezh2_S']:.3f} "
              f"maHi={md['ezh2_maHigh']:.3f} cycle={md['ezh2_cycle']:.3f}")
        for hu in (0.5, 1.0, 2.0):
            r = go(hu, Ki)
            mh = metrics(r, thr)
            bS = mh['ezh2_S'] / md['ezh2_S'] if md['ezh2_S'] else float('nan')
            bC = mh['ezh2_cycle'] / md['ezh2_cycle'] if md['ezh2_cycle'] else float('nan')
            print(f"  HU={hu:<4}: per={period(r):.1f}h  maHi_dur={mh['maHigh_dur']:.2f}h "
                  f"maHi_frac={mh['maHigh_frac']:.3f}  EZH2: S={mh['ezh2_S']:.3f} "
                  f"maHi={mh['ezh2_maHigh']:.3f} cycle={mh['ezh2_cycle']:.3f}  "
                  f"| boost_S={bS:.2f} boost_cycle={bC:.2f} (tgt 1.31)")
