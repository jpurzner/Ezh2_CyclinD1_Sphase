"""v44 crux experiment: explicit replication timer that SUSTAINS S until DNA is replicated.

Verification showed the CyclinA/S window self-terminates via E2F autoinhibition
(CyclinA -> phosphorylates E2F -> E2F off -> CyclinA synthesis stops), independent of
mitosis. So to make S concentration-dependent we GATE THE S-EXIT on replication state:
while replication is incomplete (Rep < 1), block CyclinA's inactivation of E2F, so E2F
stays on -> CyclinA sustained -> the cell stays in S. Replication progresses at fork
speed v_fork(HU); HU (dNTP depletion) slows forks -> Rep lags -> S genuinely lengthens.

This is real biology (S-phase length set by replication kinetics; HU prolongs S), and
acts on a node the pacemaker actually feels (V1e2f moved S in the sensitivity scan).

Rep dynamics:
  Rep' = k_rep * v_fork(HU) * sphase * (1 - Rep) * eps     (replicate while in S, to 100%)
  Rep' -= k_reset * (1 - sphase) * Rep * eps               (reset in G1/M)
  R_done := Rep^n/(K_done^n + Rep^n)                        (gate: ~1 only when replicated)
Gate S-exit:
  E2F_phosphorylation (E2F off by CyclinA)  *= R_done       (A: sustain S until replicated)
  optionally Mb_activation *= R_done                        (B: also block mitosis until done)

Arms vs HU; report period, S_dwell, S_hours, divisions, Cb_max (artifact watch).
If S grows with HU and the cycle stays healthy -> surgical restructure is viable.

Run:  ./venv/bin/python simulations/diag_v44_replication_timer.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from src.build_model_v43_checkpoint import build_model_v43
import tellurium as te

REW = dict(K_Ma_EZH2=0.35, K_Me_EZH2=0.35, w_Me=0.65,
           k_EZH2_deg=0.20, k_EZH2_translation=0.80)
E2F_RXN = "E2F_phosphorylation: E2F -> E2Fp; V1e2f * (E2F / (K1e2f + E2F)) * Ma * eps;"
MB_RXN = "Mb_activation: Mbi -> Mb; Vm1b * (Mbi / (K1b + Mbi)) * Pb * eps;"
_BASE = build_model_v43(rewire_ezh2=True, ezh2_rewire_params=REW).replace(
    "eps := eps0 * (1 - sphase * (1 - phi_HU));", "eps := eps0;")
assert E2F_RXN in _BASE and MB_RXN in _BASE

REP_BLOCK = """
// === explicit DNA-replication timer (sustains S until replicated) ===
species Rep = 0;
v0_fork = 1.0;
K_HU_fork = 0.6;
h_fork = 3;
k_rep = {k_rep};
k_reset = 4.0;
K_done = 0.9;
n_done = 8;
v_fork := v0_fork * K_HU_fork^h_fork / (K_HU_fork^h_fork + HU^h_fork);
R_done := Rep^n_done / (K_done^n_done + Rep^n_done);
Rep_progress: -> Rep; k_rep * v_fork * sphase * (1 - Rep) * eps;
Rep_clearance: Rep -> ; k_reset * (1 - sphase) * Rep * eps;

"""


def make(k_rep, gate_mb=False):
    m = _BASE.replace("// INITIAL CONDITIONS",
                      REP_BLOCK.format(k_rep=k_rep) + "// INITIAL CONDITIONS")
    m = m.replace(E2F_RXN, "E2F_phosphorylation: E2F -> E2Fp; V1e2f * "
                           "(E2F / (K1e2f + E2F)) * Ma * R_done * eps;")
    if gate_mb:
        m = m.replace(MB_RXN, "Mb_activation: Mbi -> Mb; Vm1b * (Mbi / (K1b + Mbi)) "
                              "* Pb * R_done * eps;")
    return m


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


def row(lbl, r):
    m = r['time'] >= 80
    sdw = (r['[Ma]'][m] > thr_ma).mean(); p = period(r)
    ez_s = r['[EZH2]'][m & (r['[Ma]'] > thr_ma)]
    print(f"  {lbl:>14}{p:>8.1f}{ndiv(r):>5d}{sdw:>9.3f}{sdw*p:>7.1f}"
          f"{r['[Cb]'][m].max():>8.3f}{r['[Ma]'][m].max():>8.3f}"
          f"{(ez_s.mean() if ez_s.size else float('nan')):>9.3f}")


print("=" * 86)
print(f"v44 REPLICATION TIMER -> concentration-dependent S  (thr_ma={thr_ma:.3f})")
print("=" * 86)
print(f"  {'arm':>14}{'period':>8}{'div':>5}{'S_dwell':>9}{'S_hrs':>7}{'Cb_max':>8}"
      f"{'Ma_max':>8}{'EZH2_S':>9}")
row("BASELINE", base0)

for k_rep in (0.3, 0.7, 1.5):
    print(f"\n[A] gate S-exit only (E2F-shutoff * R_done), k_rep={k_rep}:")
    for hu in (0.0, 0.5, 1.0, 2.0):
        row(f"HU={hu}", sim(make(k_rep, gate_mb=False), hu))

print(f"\n[B] gate S-exit + mitosis (E2F-shutoff & Mb * R_done), k_rep=0.7:")
for hu in (0.0, 0.5, 1.0, 2.0):
    row(f"HU={hu}", sim(make(0.7, gate_mb=True), hu))
