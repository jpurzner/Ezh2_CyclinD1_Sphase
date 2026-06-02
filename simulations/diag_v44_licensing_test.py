"""v44 follow-up: can ANY concentration lever hold mitosis / lengthen S?

diag_v44_checkpoint_test.py showed inhibiting Cdc25B (Pb) -- the switch ACTIVATOR --
even dynamically and near-completely, does NOT lengthen S (relaxation robustness;
free CycB just inflates). This script tests the last concentration option: a HARD
replication-LICENSING gate placed DIRECTLY on the Cdk1 (Mb) activation flux, i.e.
mitotic entry literally cannot proceed until replication (Rep) completes -- the
Lemmens-2018 "brake" idea. Rep integrates during S at fork speed v_fork(HU); HU
slows forks so the gate stays shut longer.

Arms:
  CONSTANT Mb gate   : Mb_activation * f  (does a direct constant throttle hold it?)
  DYNAMIC licensing  : Mb_activation * Rep^n/(K_lic^n+Rep^n), swept over HU
Also try gating Cdc20 activation (the APC exit trigger) by the same license, since
Cdc20 drives the cyclin destruction that ends the cycle.

If even a hard direct gate fails -> concentration approach is not viable in this model
topology and eps-clock-slowing is the only S-lever. If it works -> viable path found.

Run:  ./venv/bin/python simulations/diag_v44_licensing_test.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from src.build_model_v43_checkpoint import build_model_v43
import tellurium as te

REW = dict(K_Ma_EZH2=0.35, K_Me_EZH2=0.35, w_Me=0.65,
           k_EZH2_deg=0.20, k_EZH2_translation=0.80)
MB_RXN = "Mb_activation: Mbi -> Mb; Vm1b * (Mbi / (K1b + Mbi)) * Pb * eps;"
CDC20_RXN = "Cdc20_activation: Cdc20i -> Cdc20a; Vm3b * (Cdc20i / (K3b + Cdc20i)) * Mb * eps;"
EPS_GATED = "eps := eps0 * (1 - sphase * (1 - phi_HU));"
EPS_PLAIN = "eps := eps0;"
_BASE = build_model_v43(rewire_ezh2=True, ezh2_rewire_params=REW).replace(EPS_GATED, EPS_PLAIN)
assert MB_RXN in _BASE and CDC20_RXN in _BASE

REP_BLOCK = """
// === replication licensing variable Rep ===
species Rep = 0;
v0_fork = 1.0;
K_HU_fork = 0.6;
h_fork = 3;
k_rep = {k_rep};
k_reset = 3.0;
K_lic = 0.7;
n_lic = 8;
v_fork := v0_fork * K_HU_fork^h_fork / (K_HU_fork^h_fork + HU^h_fork);
license := Rep^n_lic / (K_lic^n_lic + Rep^n_lic);
Rep_progress: -> Rep; k_rep * v_fork * sphase * (1 - Rep) * eps;
Rep_clearance: Rep -> ; k_reset * (1 - sphase) * Rep * eps;

"""


def make_lic(k_rep, target):
    m = _BASE.replace("// INITIAL CONDITIONS",
                      REP_BLOCK.format(k_rep=k_rep) + "// INITIAL CONDITIONS")
    if target == "Mb":
        m = m.replace(MB_RXN, "Mb_activation: Mbi -> Mb; Vm1b * (Mbi / (K1b + Mbi)) "
                              "* Pb * license * eps;")
    elif target == "Cdc20":
        m = m.replace(CDC20_RXN, "Cdc20_activation: Cdc20i -> Cdc20a; Vm3b * "
                                 "(Cdc20i / (K3b + Cdc20i)) * Mb * license * eps;")
    return m


def make_const_mb(f):
    return _BASE.replace(MB_RXN, f"Mb_activation: Mbi -> Mb; Vm1b * (Mbi / (K1b + Mbi)) * Pb * {f} * eps;")


def sim(model, hu, t_end=500, n_pts=50000):
    rr = te.loada(model)
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 1.0; rr['HU'] = hu
    rr['k_EZH2_mRNA_synth_basal'] = 0.02; rr['k_EZH2_mRNA_synth_E2F'] = 3.5
    return rr.simulate(0, t_end, n_pts)


def period(r, t0=80):
    t = r['time']; cb = r['[Cb]']; m = t >= t0
    pk, _ = find_peaks(cb[m], prominence=0.05, distance=100); pt = t[m][pk]
    return float(np.mean(np.diff(pt))) if len(pt) > 1 else float('nan')


def ndiv(r, t0=80):
    t = r['time']; cb = r['[Cb]']; m = t >= t0
    pk, _ = find_peaks(cb[m], prominence=0.05, distance=100); return len(pk)


base0 = sim(_BASE, 0.0)
thr_ma = 0.20 * base0['[Ma]'][base0['time'] >= 80].max()


def row(label, r):
    m = r['time'] >= 80
    sdw = (r['[Ma]'][m] > thr_ma).mean(); p = period(r)
    print(f"  {label:>16}{p:>8.1f}{ndiv(r):>5d}{sdw:>9.3f}{sdw*p:>7.1f}"
          f"{r['[Cb]'][m].max():>8.3f}{r['[Ma]'][m].max():>8.3f}")


print("=" * 80)
print(f"v44 LICENSING TEST  (thr_ma={thr_ma:.3f})  -- can a HARD gate hold mitosis?")
print("=" * 80)
print(f"  {'arm':>16}{'period':>8}{'div':>5}{'S_dwell':>9}{'S_hrs':>7}{'Cb_max':>8}{'Ma_max':>8}")
row("BASELINE", base0)

print("\nCONSTANT Mb-activation throttle (does direct constant gate hold?):")
for f in (0.5, 0.1, 0.02):
    row(f"Mb*{f}", sim(make_const_mb(f), 0.0))

print("\nDYNAMIC licensing gate on Mb-activation (Rep^8/(0.7^8+Rep^8)) vs HU:")
for k_rep in (0.3, 1.0):
    print(f"  -- k_rep={k_rep} --")
    for hu in (0.0, 0.5, 1.0, 2.0):
        row(f"Mb-lic HU={hu}", sim(make_lic(k_rep, "Mb"), hu))

print("\nDYNAMIC licensing gate on Cdc20-activation vs HU (k_rep=1.0):")
for hu in (0.0, 1.0, 2.0):
    row(f"Cdc20-lic HU={hu}", sim(make_lic(1.0, "Cdc20"), hu))
