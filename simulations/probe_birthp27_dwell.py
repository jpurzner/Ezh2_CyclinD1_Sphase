"""Calibration probe: map BIRTH p27 (P21_div) -> G0 dwell at the baked threshold.

Tells us the P21_div range that spans immediate re-entry (dwell < 2h) <-> transient-G0
<-> arrest, so the bifurcation ensemble's birth-p27 distribution can be set to straddle it.
Sweeps P21_div at the GNP CyclinD1 setpoint (cd_scale=1), and again with EZH2i.
Fast: ~7 values x 2 = 14 sims.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

T_END, N_PTS, SETTLE = 11000, 22000, 5500
KTL0 = 0.801
_RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
_RR.integrator.setValue("absolute_tolerance", 1e-9)
_RR.integrator.setValue("relative_tolerance", 1e-6)
SEL = ["time", "MPF", "P21", "aRc", "Dna", "Cd"]


def dwell_of(p21d, shh=0.5, ezh2i=0, cd_scale=1.0):
    _RR.reset()
    _RR['SHH'] = shh; _RR['EZH2i'] = ezh2i; _RR['MYCN_amplification'] = 1.0
    _RR['Ptch1_copy_number'] = 1.0; _RR['P21_div'] = p21d
    _RR['k_Cd_translation'] = KTL0 * cd_scale
    _RR['p16'] = 0.0; _RR['p18'] = 0.464; _RR['kSyP21'] = 0.002
    try:
        r = _RR.simulate(0, T_END, N_PTS, selections=SEL)
    except Exception:
        return None, None
    t = r['time']; m = t >= SETTLE; tt = t[m]; dt = tt[1] - tt[0]
    pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=int(200 / dt))
    if len(pk) < 2:
        return 'arrest', np.nan
    per = float(np.mean(np.diff(tt[pk]))) / 60.0
    P21 = r['P21'][m]; aRc = r['aRc'][m]; Dna = r['Dna'][m]
    inS = (aRc > 0.05) & (Dna < 0.98); inG2 = Dna >= 0.98; preS = ~inS & ~inG2
    dwell = per * (preS & (P21 > 0.1)).mean() / ((preS | inS | inG2).mean() or 1)
    return ('immediate' if dwell < 2.0 else 'transient_G0'), dwell


for tag, kw in [("GNP (SHH0.5)", dict()), ("GNP+EZH2i", dict(ezh2i=1))]:
    print(f"\n{tag}:  P21_div -> (class, dwell h)")
    for p in [0.15, 0.3, 0.45, 0.6, 0.9, 1.2, 1.6, 2.2, 3.0]:
        cls, dw = dwell_of(p, **kw)
        print(f"   P21_div={p:4.2f}  ->  {cls:13s} {dw if dw is not None else float('nan'):5.2f}")
