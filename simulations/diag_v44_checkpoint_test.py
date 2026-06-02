"""v44 VERIFICATION: does a DYNAMIC, replication-gated CHK1 -> Cdc25B inhibition
lengthen S, where a CONSTANT Cdc25B inhibition does NOT?

This is the gating experiment for the whole concentration-dependent-checkpoint plan.
We KNOW (diag_phi_levers.py) that scaling Cdc25/CDK activation by a CONSTANT factor
leaves the period flat (relaxation-oscillator robustness to amplitude). The claim to
test: a checkpoint is a STATE-DEPENDENT gate, not a constant -- CHK1 holds the system
at the mitotic threshold until replication completes, so the dwell (S/G2) lengthens in
a CONCENTRATION-DEPENDENT way and ESCAPES the robustness.

Design (string surgery on the v42/v43 core; eps phase-gating disabled so HU drives ONLY
the new checkpoint):
  Rep      : replication progress, integrates during S (sphase = CycA-high) at fork speed
             v_fork(HU); HU slows forks (dNTP depletion) so Rep lags.
  CHK1     := sphase * K_Rep^m/(K_Rep^m + Rep^m)   (active = replicating AND incomplete)
  Cdc25B   : Cdc25B_activation * K_inh/(K_inh + CHK1)  -> CHK1 holds Cdc25B down -> Cdk1
             (Mb) activation delayed -> mitosis delayed -> CycA(Ma) persists -> S prolonged.
  When Rep completes -> CHK1 clears -> Pb activates -> mitosis fires -> reset. No permanent
  arrest (Rep always completes eventually).

Three arms:
  BASELINE  : no checkpoint
  CONSTANT  : Cdc25B_activation * f  (f fixed)  -> expect period ~flat (robustness)
  DYNAMIC   : the CHK1 checkpoint above, swept over HU -> expect graded S prolongation

Run:  ./venv/bin/python simulations/diag_v44_checkpoint_test.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from src.build_model_v43_checkpoint import build_model_v43
import tellurium as te

REW = dict(K_Ma_EZH2=0.35, K_Me_EZH2=0.35, w_Me=0.65,
           k_EZH2_deg=0.20, k_EZH2_translation=0.80)
PB_RXN = "Cdc25B_activation: Pbi -> Pb; Vm5b * (Mb + ab) * (Pbi / (K5b + Pbi)) * eps;"
EPS_GATED = "eps := eps0 * (1 - sphase * (1 - phi_HU));"
EPS_PLAIN = "eps := eps0;"

# Base model with EZH2 rewire; disable the v43 eps phase-gating so HU is free for the
# new checkpoint. sphase (Ma^n_S/(K_S_phi^n_S+Ma^n_S)) stays defined and is reused.
_BASE = build_model_v43(rewire_ezh2=True, ezh2_rewire_params=REW).replace(EPS_GATED, EPS_PLAIN)
assert PB_RXN in _BASE and EPS_PLAIN in _BASE

CHK_BLOCK = """
// === v44 test: replication-gated CHK1 -> Cdc25B checkpoint ===
species Rep = 0;
v0_fork = 1.0;
K_HU_fork = 0.6;
h_fork = 3;
k_rep = {k_rep};
k_reset = 3.0;
K_Rep_chk = 0.8;
m_chk = 6;
K_inh = {K_inh};
v_fork := v0_fork * K_HU_fork^h_fork / (K_HU_fork^h_fork + HU^h_fork);
CHK1 := sphase * (K_Rep_chk^m_chk / (K_Rep_chk^m_chk + Rep^m_chk));
Rep_progress: -> Rep; k_rep * v_fork * sphase * (1 - Rep) * eps;
Rep_clearance: Rep -> ; k_reset * (1 - sphase) * Rep * eps;

"""


def make_dynamic(k_rep, K_inh):
    m = _BASE.replace("// INITIAL CONDITIONS",
                      CHK_BLOCK.format(k_rep=k_rep, K_inh=K_inh) + "// INITIAL CONDITIONS")
    pb_new = ("Cdc25B_activation: Pbi -> Pb; Vm5b * (Mb + ab) * (Pbi / (K5b + Pbi)) "
              "* (K_inh / (K_inh + CHK1)) * eps;")
    return m.replace(PB_RXN, pb_new)


def make_constant(f):
    pb_new = ("Cdc25B_activation: Pbi -> Pb; Vm5b * (Mb + ab) * (Pbi / (K5b + Pbi)) "
              f"* {f} * eps;")
    return _BASE.replace(PB_RXN, pb_new)


def sim(model, hu, t_end=400, n_pts=40000):
    rr = te.loada(model)
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 1.0; rr['HU'] = hu
    rr['k_EZH2_mRNA_synth_basal'] = 0.02; rr['k_EZH2_mRNA_synth_E2F'] = 3.5
    return rr.simulate(0, t_end, n_pts)


def period(r, t0=80):
    t = r['time']; cb = r['[Cb]']; m = t >= t0
    pk, _ = find_peaks(cb[m], prominence=0.05, distance=100)
    pt = t[m][pk]
    return float(np.mean(np.diff(pt))) if len(pt) > 1 else float('nan')


def ndiv(r, t0=80):
    t = r['time']; cb = r['[Cb]']; m = t >= t0
    pk, _ = find_peaks(cb[m], prominence=0.05, distance=100)
    return len(pk)


# fixed DMSO Ma threshold from baseline
base0 = sim(_BASE, 0.0)
thr_ma = 0.20 * base0['[Ma]'][base0['time'] >= 80].max()


def stats(r):
    m = r['time'] >= 80
    sdw = (r['[Ma]'][m] > thr_ma).mean()
    return period(r), ndiv(r), sdw, r['[Cb]'][m].max(), r['[Ma]'][m].max()


print("=" * 78)
print("v44 CHECKPOINT VERIFICATION  (thr_ma=%.3f)" % thr_ma)
print("=" * 78)

p, d, s, cbm, mam = stats(base0)
print(f"\nBASELINE (no checkpoint):  period={p:.1f}h  div={d}  S_dwell={s:.3f}  "
      f"Cb_max={cbm:.3f}  Ma_max={mam:.3f}")

print("\nCONSTANT Cdc25B inhibition (expect period ~FLAT = robustness):")
print(f"{'f':>6}{'period':>8}{'div':>5}{'S_dwell':>9}{'Cb_max':>8}{'Ma_max':>8}")
for f in (1.0, 0.5, 0.2, 0.1):
    p, d, s, cbm, mam = stats(sim(make_constant(f), 0.0))
    print(f"{f:>6.2f}{p:>8.1f}{d:>5d}{s:>9.3f}{cbm:>8.3f}{mam:>8.3f}")

print("\nDYNAMIC replication-gated CHK1 checkpoint (expect S GROWS with HU):")
for k_rep in (0.3, 0.7, 1.5):
    for K_inh in (0.2,):
        print(f"\n  [k_rep={k_rep} K_inh={K_inh}]")
        print(f"  {'HU':>5}{'v_fork':>8}{'period':>8}{'div':>5}{'S_dwell':>9}"
              f"{'S_hrs':>7}{'Cb_max':>8}{'Ma_max':>8}")
        for hu in (0.0, 0.5, 1.0, 2.0):
            r = sim(make_dynamic(k_rep, K_inh), hu)
            p, d, s, cbm, mam = stats(r)
            vf = 1.0 * 0.6**3 / (0.6**3 + hu**3)
            print(f"  {hu:>5.1f}{vf:>8.3f}{p:>8.1f}{d:>5d}{s:>9.3f}"
                  f"{s*p:>7.1f}{cbm:>8.3f}{mam:>8.3f}")
