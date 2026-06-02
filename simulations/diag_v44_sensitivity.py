"""v44 decisive diagnostic: WHICH rate constants actually set S-phase duration?

Mitotic-module concentration levers (Cdc25B, Cdk1/Mb, Cdc20) all failed to lengthen S
(diag_v44_*_test.py). Hypothesis: in this parameterization the PACEMAKER is the
Rb-E2F-CyclinA loop, and CyclinB/Cdk1/APC is a SLAVED output -- so the S/CyclinA window
is terminated by the E2F<->CyclinA negative feedback (CyclinA/Ma phosphorylates E2F ->
E2Fp -> CyclinA synthesis stops -> CyclinA decays), NOT by mitosis. If true, S-duration
is controlled by the E2F/CyclinA rate constants (esp. V1e2f = rate CyclinA inactivates
E2F), and the biological checkpoint target (Cdc25/Cdk1) has no leverage here.

This scan multiplies each candidate rate constant by 0.5x and 2x and reports the effect
on period and S_dwell (fixed DMSO Ma threshold). Whichever knobs move S_dwell are where
S-duration is actually set -> the only places a concentration-dependent checkpoint could
act.

Run:  ./venv/bin/python simulations/diag_v44_sensitivity.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from src.build_model_v43_checkpoint import build_model_v43
import tellurium as te
import re

REW = dict(K_Ma_EZH2=0.35, K_Me_EZH2=0.35, w_Me=0.65,
           k_EZH2_deg=0.20, k_EZH2_translation=0.80)
_BASE = build_model_v43(rewire_ezh2=True, ezh2_rewire_params=REW).replace(
    "eps := eps0 * (1 - sphase * (1 - phi_HU));", "eps := eps0;")

# (param, baseline value, module)
PARAMS = [
    ("V1e2f", "E2F<-Ma inactivation (S terminator)"),
    ("vse2f", "E2F synthesis"),
    ("V2e2f", "E2Fp->E2F dephos"),
    ("kde2fp", "E2Fp decay"),
    ("kca",   "CyclinA synthesis"),
    ("Vda",   "CyclinA degradation"),
    ("Vm1a",  "Ma activation (Cdc25A)"),
    ("kce",   "CyclinE synthesis"),
    ("V3",    "pRBp->pRBpp by Me (G1/S)"),
    ("V1",    "pRB->pRBp by Md (G1)"),
    ("vcb",   "CyclinB synthesis"),
    ("Vm1b",  "Mb activation (Cdc25B/Cdk1)"),
]


def get_val(name):
    m = re.search(rf"\b{name}\s*=\s*([0-9.]+)\s*;", _BASE)
    return float(m.group(1)) if m else None


def set_val(name, val):
    return re.sub(rf"\b{name}\s*=\s*[0-9.]+\s*;", f"{name} = {val};", _BASE, count=1)


def sim(model, t_end=400, n_pts=40000):
    rr = te.loada(model)
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 1.0; rr['HU'] = 0.0
    rr['k_EZH2_mRNA_synth_basal'] = 0.02; rr['k_EZH2_mRNA_synth_E2F'] = 3.5
    return rr.simulate(0, t_end, n_pts)


def period(r, t0=80):
    t = r['time']; cb = r['[Cb]']; m = t >= t0
    pk, _ = find_peaks(cb[m], prominence=0.05, distance=100); pt = t[m][pk]
    return float(np.mean(np.diff(pt))) if len(pt) > 1 else float('nan')


base = sim(_BASE)
thr_ma = 0.20 * base['[Ma]'][base['time'] >= 80].max()
per0 = period(base)
sd0 = (base['[Ma]'][base['time'] >= 80] > thr_ma).mean()
print("=" * 88)
print(f"WHERE IS S-DURATION SET?  baseline period={per0:.1f}h  S_dwell={sd0:.3f}  thr_ma={thr_ma:.3f}")
print("=" * 88)
print(f"{'param':>8} {'module':<34}{'val':>8} | "
      f"{'per(.5x)':>9}{'Sdw(.5x)':>9} | {'per(2x)':>9}{'Sdw(2x)':>9}")


def metrics(model):
    r = sim(model); m = r['time'] >= 80
    return period(r), (r['[Ma]'][m] > thr_ma).mean()


for name, desc in PARAMS:
    v = get_val(name)
    if v is None:
        print(f"{name:>8} {desc:<34}{'NOT FOUND':>8}")
        continue
    pl, sl = metrics(set_val(name, v * 0.5))
    ph, sh = metrics(set_val(name, v * 2.0))
    print(f"{name:>8} {desc:<34}{v:>8.4g} | {pl:>9.1f}{sl:>9.3f} | {ph:>9.1f}{sh:>9.3f}")
