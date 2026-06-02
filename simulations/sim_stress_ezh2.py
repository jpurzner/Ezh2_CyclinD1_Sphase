"""
Calibrate the replication-stress -> EZH2-transcription coupling (v43, k_EZH2_stress).

Design (per user directive: "ezh2 transcription [should be] dependent on the time
spent in s-phase"; "Cdc25E/Cdc25B/Cdc20 is too strong"):

  * ALL CHK1->cell-cycle levers OFF -> the validated Gerard-Goldbeter core is never
    perturbed, so cells NEVER arrest and the v42 validation is preserved at HU=0.
  * HU -> replication stress D_stress (HU-dosed AND Ma-gated = S-phase-specific) ->
    multiplies the Ma-driven EZH2 transcription term by (1 + k_EZH2_stress * D_stress).
    EZH2 transcription therefore depends on the replication stress accumulated while
    the cell is IN S-phase: HU-graded, S-specific, no arrest.

We sweep k_EZH2_stress to find the value giving 10 uM HU (HU~1) -> 1.31x EZH2 in S,
using ARTIFACT-PROOF metrics: phase thresholds fixed from DMSO and reused for HU;
whole-cycle and fixed-Ma-window EZH2 reported as ground truth; division count to
confirm cells still cycle.

Run:  ./venv/bin/python simulations/sim_stress_ezh2.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from calibrate_ezh2_rewire import run, period

BASE_REW = dict(k_EZH2_deg=0.20, K_Ma_EZH2=0.35, k_EZH2_translation=0.80)
BASAL, E2F_AMP = 0.02, 3.5


def go(hu, k_stress, t_end=240):
    rew = dict(BASE_REW, k_EZH2_stress=k_stress)
    # checkpoint=None -> CHECKPOINT_DEFAULTS (all cell-cycle levers OFF); D_stress still active
    return run(shh=0.5, hu=hu, ezh2_rewire_params=rew, basal=BASAL, e2f_amp=E2F_AMP,
               t_end=t_end, n_pts=int(t_end * 100))


def n_divisions(r, t_start=40):
    t = r['time']; m = t >= t_start
    cb = r['[Cb]'][m]
    if cb.max() <= 0:
        return 0
    pk, _ = find_peaks(cb, prominence=0.05 * cb.max(), distance=int(0.05 * len(cb)))
    return len(pk)


def fixed_thresholds(dmso, t_start=40):
    t = dmso['time']; m = t >= t_start
    ma = dmso['[Ma]'][m]; cb = dmso['[Cb]'][m]
    return dict(thr_ma=0.20 * ma.max(), thr_cb=0.25 * cb.max())


def metrics(r, thr, t_start=40):
    t = r['time']; m = t >= t_start
    ma = r['[Ma]'][m]; cb = r['[Cb]'][m]; ez = r['[EZH2]'][m]
    s_sel = (ma > thr['thr_ma']) & (cb <= thr['thr_cb'])
    sg2_sel = (ma > thr['thr_ma'])
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
        Dpk=float(r['[D_stress]'][m].max()),
        maHigh_dur=float(np.mean(durs)) if durs else 0.0,
        maHigh_frac=float(sg2_sel.mean()),
    )


if __name__ == "__main__":
    print("=" * 80)
    print("REPLICATION-STRESS -> EZH2-TRANSCRIPTION COUPLING  (all cell-cycle levers OFF)")
    print("=" * 80)
    print(f"baseline rewire: {BASE_REW}  basal={BASAL} E2F_amp={E2F_AMP}\n")

    for ks in (0.0, 0.14, 0.16, 0.18, 0.20):
        print(f"\n### k_EZH2_stress = {ks} ###")
        dmso = go(0.0, ks)
        thr = fixed_thresholds(dmso)
        md = metrics(dmso, thr)
        nd0 = n_divisions(dmso)
        print(f"  thr_ma={thr['thr_ma']:.3f} thr_cb={thr['thr_cb']:.3f}")
        print(f"  DMSO : per={period(dmso):.1f}h div={nd0}  Dpk={md['Dpk']:.2f} "
              f"maHi_frac={md['maHigh_frac']:.3f}  EZH2: S={md['ezh2_S']:.3f} "
              f"maHi={md['ezh2_maHigh']:.3f} cyc={md['ezh2_cycle']:.3f}")
        for hu in (0.5, 1.0, 2.0):
            r = go(hu, ks)
            mh = metrics(r, thr)
            nd = n_divisions(r)
            bS = mh['ezh2_S'] / md['ezh2_S'] if md['ezh2_S'] else float('nan')
            bC = mh['ezh2_cycle'] / md['ezh2_cycle'] if md['ezh2_cycle'] else float('nan')
            arr = "  <-- ARREST" if nd <= 1 else ""
            print(f"  HU={hu:<4}: per={period(r):.1f}h div={nd}  Dpk={mh['Dpk']:.2f} "
                  f"maHi_frac={mh['maHigh_frac']:.3f}  EZH2: S={mh['ezh2_S']:.3f} "
                  f"maHi={mh['ezh2_maHigh']:.3f} cyc={mh['ezh2_cycle']:.3f} "
                  f"| bS={bS:.2f} bC={bC:.2f} (tgt 1.31){arr}")
