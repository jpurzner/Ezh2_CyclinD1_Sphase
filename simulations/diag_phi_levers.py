"""Which GG flux is the PERIOD actually sensitive to?

Both Ma and Mb ACTIVATION-amplitude scaling leave the period flat (the bistable
Cdc25/Wee1 switch fires at the same time regardless of linear activation rate).
To make HU genuinely slow the cycle / prolong S we must couple phi to a flux the
period is sensitive to. Scan candidate levers at a FIXED throttle factor f=0.6
(mimics phi at HU~1.3) and see which one lengthens the period and/or prolongs the
Ma-high (S/G2) window WITHOUT inflating free CycB.

Each lever multiplies (or divides) one v42 reaction rate by f via string surgery
on the built model. Reports period, S_dwell (Ma>fixed DMSO thr), Ma_max, Cb_max.

Run:  ./venv/bin/python simulations/diag_phi_levers.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from src.build_model_v43_checkpoint import build_model_v43
import tellurium as te

REW = dict(K_Ma_EZH2=0.35, K_Me_EZH2=0.35, w_Me=0.65,
           k_EZH2_deg=0.20, k_EZH2_translation=0.80)

# Build base, strip the builder's phi-on-Ma coupling -> clean GG (phi_HU unused).
MA_ON = "Ma_activation: Mai -> Ma; Vm1a * (Mai / (K1a + Mai)) * Pa * phi_HU * eps;"
MA_OFF = "Ma_activation: Mai -> Ma; Vm1a * (Mai / (K1a + Mai)) * Pa * eps;"
BASE = build_model_v43(rewire_ezh2=True, ezh2_rewire_params=REW).replace(MA_ON, MA_OFF)

# Candidate levers: (label, old_substring, new_substring) with {f} formatted in.
# f<1 = throttle. We pick the SIGN biologically: intra-S checkpoint SLOWS progression
# toward mitosis and SLOWS cyclin destruction is NOT it -- test both directions.
LEVERS = {
    # slow CycA/CDK2 activation (current; control)
    "Ma_act*f":  ("Pa * eps;",          "Pa * {f} * eps;"),
    # slow CycB/CDK1 activation (mitotic entry)
    "Mb_act*f":  ("Pb * eps;",          "Pb * {f} * eps;"),
    # slow CycA degradation -> prolong high-CycA window
    "CaDeg*f":   ("(Cdc20a / (Kacdc20 + Cdc20a)) * eps;",
                  "(Cdc20a / (Kacdc20 + Cdc20a)) * {f} * eps;"),
    # slow Cdc20 activation (the APC/C exit trigger)
    "Cdc20act*f":("(Cdc20i / (K3b + Cdc20i)) * Mb * eps;",
                  "(Cdc20i / (K3b + Cdc20i)) * Mb * {f} * eps;"),
    # slow CycB synthesis
    "CbSyn*f":   ("Cb_synthesis: -> Cb; vcb * eps;",
                  "Cb_synthesis: -> Cb; vcb * {f} * eps;"),
}


def period(r, t_start=50):
    t = r['time']; cb = r['[Cb]']; m = t >= t_start
    pk, _ = find_peaks(cb[m], prominence=0.05, distance=10)
    pt = t[m][pk]
    return float(np.mean(np.diff(pt))) if len(pt) > 1 else float('nan')


def sim(model, hu=0.0, t_end=300, n_pts=30000):
    rr = te.loada(model)
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 1.0; rr['HU'] = hu
    rr['k_EZH2_mRNA_synth_basal'] = 0.02; rr['k_EZH2_mRNA_synth_E2F'] = 3.5
    return rr.simulate(0, t_end, n_pts)


# DMSO baseline + fixed Ma threshold
base0 = sim(BASE)
thr_ma = 0.20 * base0['[Ma]'][base0['time'] >= 50].max()
per0 = period(base0)
print(f"DMSO baseline: period={per0:.1f}h  thr_ma={thr_ma:.3f}")
print(f"{'lever':>14}{'period':>8}{'dPer%':>7}{'S_dwell':>9}{'Ma_max':>8}{'Cb_max':>8}")

def show(label, model):
    r = sim(model)
    m = r['time'] >= 50
    per = period(r)
    sdw = (r['[Ma]'][m] > thr_ma).mean()
    print(f"{label:>14}{per:>8.1f}{100*(per/per0-1):>7.1f}{sdw:>9.3f}"
          f"{r['[Ma]'][m].max():>8.3f}{r['[Cb]'][m].max():>8.3f}")

show("baseline", BASE)
f = 0.6
for lbl, (old, new) in LEVERS.items():
    if old is None:
        continue
    if old not in BASE:
        print(f"{lbl:>14}  [SKIP: substring not found]")
        continue
    show(lbl, BASE.replace(old, new.format(f=f)))

# global eps sanity: replace the eps assignment value
import re
m_eps = re.search(r"eps\s*=\s*([0-9.]+)\s*;", BASE)
if m_eps:
    eps_val = float(m_eps.group(1))
    show("eps*f_all", BASE.replace(m_eps.group(0), f"eps = {eps_val*f};"))
