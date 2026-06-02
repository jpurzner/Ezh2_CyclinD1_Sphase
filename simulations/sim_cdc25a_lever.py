"""
Test the Cdc25A (Pa) checkpoint lever (v43) -- the canonical CHK1 -> CDC25A -> CDK2
axis acting on S-PROGRESSION, not S-entry or M-exit.

Cdc25A activates Ma = active CycA/CDK2, the kinase that drives DNA replication
(S-phase progression). Throttling it should STRETCH the high-Ma window (time-in-S)
rather than block S-entry (Pe lever -> arrest) or merely inflate free CycB (Pb/Cdc20
levers -> M-exit artifact). This is the lever the user flagged as correct:
"the Cdc25E/Cdc25B/Cdc20 is too strong ... ezh2 transcription [should be] dependent
on the time spent in s-phase."

Artifact-proof metrics (same discipline as sim_cdc20_lever.py):
  * phase thresholds derived ONCE from DMSO, reused for every HU dose (one gate)
  * whole-cycle mean EZH2 and EZH2 over a FIXED absolute Ma window are ground truth
  * DIVISION COUNT per condition confirms cells still cycle (no full arrest)

Targets: 10 uM HU -> 1.31x EZH2 in S; S-phase prolonged but NOT abolished; period up.

Run:  ./venv/bin/python simulations/sim_cdc25a_lever.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from calibrate_ezh2_rewire import run, period

REW = dict(k_EZH2_deg=0.20, K_Ma_EZH2=0.35, k_EZH2_translation=0.80)
BASAL, E2F_AMP = 0.02, 3.5
OFF = 1.0e6


def go(hu, Ki_Pa, t_end=240):
    chk = dict(Ki_CHK1_Pa=Ki_Pa, Ki_CHK1_Pe=OFF, Ki_CHK1_Pb=OFF, Ki_CHK1_Cdc20=OFF)
    return run(shh=0.5, hu=hu, ezh2_rewire_params=REW, basal=BASAL, e2f_amp=E2F_AMP,
               checkpoint=chk, t_end=t_end, n_pts=int(t_end * 100))


def n_divisions(r, t_start=40):
    """Count cell divisions = CycB (Cb) peaks after settling."""
    t = r['time']; m = t >= t_start
    cb = r['[Cb]'][m]
    if cb.max() <= 0:
        return 0
    pk, _ = find_peaks(cb, prominence=0.05 * cb.max(),
                       distance=int(0.05 * len(cb)))
    return len(pk)


def fixed_thresholds(dmso, t_start=40):
    t = dmso['time']; m = t >= t_start
    ma = dmso['[Ma]'][m]; cb = dmso['[Cb]'][m]
    return dict(thr_ma=0.20 * ma.max(), thr_cb=0.25 * cb.max())


def metrics(r, thr, t_start=40):
    t = r['time']; m = t >= t_start
    ma = r['[Ma]'][m]; cb = r['[Cb]'][m]; ez = r['[EZH2]'][m]
    s_sel = (ma > thr['thr_ma']) & (cb <= thr['thr_cb'])   # S: replicating, pre-mitotic
    sg2_sel = (ma > thr['thr_ma'])                          # broad high-CycA window
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
        ma_peak=float(ma.max()),
        maHigh_dur=float(np.mean(durs)) if durs else 0.0,
        maHigh_frac=float(sg2_sel.mean()),
    )


if __name__ == "__main__":
    print("=" * 78)
    print("Cdc25A (Pa) S-PROGRESSION LEVER  (Pe/Pb/Cdc20 OFF; fixed DMSO thresholds)")
    print("=" * 78)
    print(f"baseline rewire: {REW}  basal={BASAL} E2F_amp={E2F_AMP}\n")

    for Ki in (OFF, 5.75, 5.5, 5.25, 5.0, 4.75, 4.5):
        tag = "OFF" if Ki >= 1e5 else f"{Ki}"
        print(f"\n### Ki_CHK1_Pa = {tag} ###")
        dmso = go(0.0, Ki)
        thr = fixed_thresholds(dmso)
        md = metrics(dmso, thr)
        nd0 = n_divisions(dmso)
        print(f"  thresholds (from DMSO): thr_ma={thr['thr_ma']:.3f} thr_cb={thr['thr_cb']:.3f}")
        print(f"  DMSO : per={period(dmso):.1f}h  div={nd0}  maHi_dur={md['maHigh_dur']:.2f}h "
              f"maHi_frac={md['maHigh_frac']:.3f} ma_pk={md['ma_peak']:.2f}  "
              f"EZH2: S={md['ezh2_S']:.3f} maHi={md['ezh2_maHigh']:.3f} cyc={md['ezh2_cycle']:.3f}")
        for hu in (0.5, 1.0):
            r = go(hu, Ki)
            mh = metrics(r, thr)
            nd = n_divisions(r)
            bS = mh['ezh2_S'] / md['ezh2_S'] if md['ezh2_S'] else float('nan')
            bC = mh['ezh2_cycle'] / md['ezh2_cycle'] if md['ezh2_cycle'] else float('nan')
            arr = "  <-- ARREST" if nd <= 1 else ""
            print(f"  HU={hu:<4}: per={period(r):.1f}h  div={nd}  maHi_dur={mh['maHigh_dur']:.2f}h "
                  f"maHi_frac={mh['maHigh_frac']:.3f} ma_pk={mh['ma_peak']:.2f}  "
                  f"EZH2: S={mh['ezh2_S']:.3f} maHi={mh['ezh2_maHigh']:.3f} cyc={mh['ezh2_cycle']:.3f} "
                  f"| bS={bS:.2f} bC={bC:.2f} (tgt 1.31){arr}")
