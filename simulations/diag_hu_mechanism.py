"""
Mechanism check for the HU->EZH2-in-S boost (v43, EZH2 rewired to Ma).

Questions:
  1. Is the boost caused by the checkpoint? -> control with both Cdc25 levers OFF
     should give boost ~1.0.
  2. Which lever causes it (and with which sign)?  Pe-only vs Pb-only vs both.
  3. What does HU actually do to the trajectory (prolong S? arrest in G2 with low CycB?
     raise Ma amplitude?) -> dump Ma, Cb, Me, EZH2 over a few cycles, DMSO vs HU.
  4. Is the boost a phase-detector artifact? -> also report EZH2 over the broader
     Ma-high (S+G2M) window and the whole-cycle mean.

Run:  ./venv/bin/python simulations/diag_hu_mechanism.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from calibrate_ezh2_rewire import (
    run, phase_call, mean_ezh2_by_phase, period, mean_last)

REW = dict(k_EZH2_deg=0.20, K_Ma_EZH2=0.35, k_EZH2_translation=0.80)
BASAL, E2F_AMP = 0.02, 3.5
OFF = 1.0e6


def ezh2_windows(result, t_start=40):
    """EZH2 in S only, in S+G2M (broad Ma-high), and whole-cycle mean."""
    phase, _ = phase_call(result, t_start)
    ez = result['[EZH2]']
    s_sel = np.array([p == 'S' for p in phase])
    sg2_sel = np.array([p in ('S', 'G2M') for p in phase])
    cyc_sel = np.array([p is not None for p in phase])
    return (float(np.mean(ez[s_sel])) if s_sel.sum() else np.nan,
            float(np.mean(ez[sg2_sel])) if sg2_sel.sum() else np.nan,
            float(np.mean(ez[cyc_sel])) if cyc_sel.sum() else np.nan)


def go(hu, Ki_Pe, Ki_Pb, t_end=240):
    chk = dict(Ki_CHK1_Pe=Ki_Pe, Ki_CHK1_Pb=Ki_Pb)
    return run(shh=0.5, hu=hu, ezh2_rewire_params=REW, basal=BASAL, e2f_amp=E2F_AMP,
               checkpoint=chk, t_end=t_end, n_pts=int(t_end * 100))


if __name__ == "__main__":
    print("=" * 78)
    print("HU MECHANISM CHECK")
    print("=" * 78)

    levers = [
        ("both OFF (control)", OFF, OFF),
        ("Pe only",           0.40, OFF),
        ("Pb only",           0.40, OFF if False else OFF),  # placeholder, set below
    ]
    # explicit configs
    configs = [
        ("both OFF (control)", OFF,  OFF),
        ("Pe only (Pb off)",   0.20, OFF),
        ("Pb only (Pe off)",   OFF,  0.20),
        ("both on",            0.20, 0.20),
    ]
    print("\n[1-2] Boost decomposition (EZH2 in S):  HU=1 vs DMSO")
    print(f"{'config':22s} {'DMSO_S':>8} {'HU_S':>8} {'boost':>7} {'per_DMSO':>9} {'per_HU':>8}")
    for name, kpe, kpb in configs:
        d = go(0.0, kpe, kpb); h = go(1.0, kpe, kpb)
        ds, _, _ = ezh2_windows(d); hs, _, _ = ezh2_windows(h)
        boost = hs / ds if ds else float('nan')
        print(f"{name:22s} {ds:8.3f} {hs:8.3f} {boost:7.2f} "
              f"{period(d):9.1f} {period(h):8.1f}")

    print("\n[4] Window robustness (Pb only, Ki=0.2):  EZH2 in S / in S+G2M / whole-cycle")
    for hu in (0.0, 1.0):
        r = go(hu, OFF, 0.20)
        s, sg2, cyc = ezh2_windows(r)
        print(f"  HU={hu}: EZH2_S={s:.3f}  EZH2_S+G2M={sg2:.3f}  EZH2_cycle={cyc:.3f}")

    # [3] trajectory dump: DMSO vs HU, Pb only
    print("\n[3] Trajectory over ~3 cycles (Pb only, Ki=0.2). t, Ma, Cb, Me, EZH2")
    for hu in (0.0, 1.0):
        r = go(hu, OFF, 0.20, t_end=240)
        t = r['time']
        cb = r['[Cb]']
        pk, _ = find_peaks(cb, prominence=0.05, distance=200)
        # show window spanning last ~3 divisions
        if len(pk) >= 4:
            i0, i1 = pk[-4], pk[-1]
        else:
            i0, i1 = len(t) - 6000, len(t) - 1
        print(f"\n  --- HU={hu} ---  (divisions at t={np.round(t[pk][-4:],1)})")
        idxs = np.linspace(i0, i1, 22).astype(int)
        for i in idxs:
            print(f"   t={t[i]:6.1f}  Ma={r['[Ma]'][i]:.3f}  Cb={r['[Cb]'][i]:.3f}  "
                  f"Me={r['[Me]'][i]:.3f}  CHK1={r['[CHK1_active]'][i]:.2f}  "
                  f"EZH2={r['[EZH2]'][i]:.3f}")
