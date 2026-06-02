"""
Phase-3 linchpin test: does checkpoint-driven S-phase prolongation boost EZH2-in-S?

With EZH2 transcription rewired to track Ma (CycA/CDK2), the plan's central hypothesis
("longer S -> more EZH2") becomes mechanistically testable. The intra-S checkpoint
should prolong the high-Ma S/G2 window; if EZH2 tracks Ma, that window's mean EZH2 rises.

Experimental target (Fig. 4G): 10 uM HU -> 1.31-fold EZH2 in S-phase cells (p=0.005).

Two checkpoint levers (CHK1 -> Cdc25):
  Pe (Cdc25E, CycE/CDK2): throttles S-ENTRY. Expected to LOWER Ma -> wrong sign.
  Pb (Cdc25B, CycB/CDK1): delays MITOSIS -> cell held in high-Ma G2 -> EZH2 up. Right sign.

This script sweeps HU dose x Ki_CHK1_Pb (Pe disabled) and reports, vs DMSO:
  - cell-cycle period
  - duration of the Ma-high (S/G2) window
  - mean EZH2 in S-phase  -> HU/DMSO ratio (target 1.31)

Run:  ./venv/bin/python simulations/sim_hu_dose_response.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from calibrate_ezh2_rewire import (
    run, phase_call, mean_ezh2_by_phase, s_phase_duration, period, mean_last)

# Calibrated baseline EZH2 rewire (HU=0 shape; see calibrate_ezh2_rewire.py).
REW = dict(k_EZH2_deg=0.20, K_Ma_EZH2=0.35, k_EZH2_translation=0.80)
BASAL, E2F_AMP = 0.02, 3.5

# Disable the Pe lever (wrong sign); drive the checkpoint through Pb (mitotic delay).
DISABLE_PE = 1.0e6


def ma_high_duration(result, t_start=40):
    """Mean contiguous duration (h) of the Ma-high window (S+G2), per cycle.

    This is the window over which EZH2 accumulates; the checkpoint should lengthen it.
    Uses S or G2M phase (i.e. NOT G1) as 'replicating/post-replicative high-CycA'.
    """
    t = result['time']
    phase, _ = phase_call(result, t_start)
    durs, on, start = [], False, None
    for i in range(len(t)):
        hi = phase[i] in ('S', 'G2M')
        if hi and not on:
            on, start = True, t[i]
        elif not hi and on:
            on = False
            durs.append(t[i] - start)
    return (float(np.mean(durs)) if durs else 0.0)


def ezh2_in_S(result, t_start=40):
    """Mean EZH2 over the replication window (S phase: Ma-high, Cb-low)."""
    by, _ = mean_ezh2_by_phase(result, t_start)
    return by['S']


def assess(hu, Ki_Pb, t_end=240):
    chk = dict(Ki_CHK1_Pe=DISABLE_PE, Ki_CHK1_Pb=Ki_Pb)
    r = run(shh=0.5, hu=hu, ezh2_rewire_params=REW, basal=BASAL, e2f_amp=E2F_AMP,
            checkpoint=chk, t_end=t_end, n_pts=int(t_end * 100))
    return dict(
        period=period(r),
        ma_dur=ma_high_duration(r),
        ezh2_S=ezh2_in_S(r),
        ezh2_mean=mean_last(r, 'EZH2'),
        chk1=mean_last(r, 'CHK1_active'),
        dstress=mean_last(r, 'D_stress'),
    )


if __name__ == "__main__":
    print("=" * 78)
    print("HU DOSE RESPONSE  (EZH2 rewired to Ma; Pe lever OFF, Pb lever active)")
    print("=" * 78)
    print(f"baseline rewire: {REW}  basal={BASAL} E2F_amp={E2F_AMP}\n")

    for Ki_Pb in (0.40, 0.20, 0.10, 0.05):
        print(f"\n### Ki_CHK1_Pb = {Ki_Pb}  (smaller = stronger mitotic delay) ###")
        base = assess(hu=0.0, Ki_Pb=Ki_Pb)
        print(f"  DMSO : period={base['period']:.1f}h  Ma-win={base['ma_dur']:.2f}h  "
              f"EZH2_S={base['ezh2_S']:.3f}  EZH2_mean={base['ezh2_mean']:.3f}")
        for hu in (0.5, 1.0, 2.0):
            a = assess(hu=hu, Ki_Pb=Ki_Pb)
            boost = a['ezh2_S'] / base['ezh2_S'] if base['ezh2_S'] else float('nan')
            print(f"  HU={hu:<4}: period={a['period']:.1f}h  Ma-win={a['ma_dur']:.2f}h  "
                  f"EZH2_S={a['ezh2_S']:.3f}  boost={boost:.2f} (target 1.31)  "
                  f"CHK1={a['chk1']:.2f} Dstr={a['dstress']:.2f}")
