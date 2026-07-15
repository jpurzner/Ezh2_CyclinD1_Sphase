"""STRUCTURE ANALYSIS — step B: timescale separation & the reduced core.

Per-species relaxation timescale tau_i = 1/|J_ii| (how fast species i forgets a perturbation), plus the full
Jacobian eigenvalue spectrum, at settled operating points. FAST (small tau) species are quasi-steady-state
followers; SLOW (large tau) species are the real degrees of freedom that must be tracked. The spectral gap
locates the fast/slow cut. Output: a ranked timescale table + histogram, and the proposed slow core.

Run: ./venv/bin/python simulations/struct/B_timescales.py
"""
import os, sys, json
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1')
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import tellurium as te
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.build_model_v44_heldt import build_model_v44
from A_loops import MODULES, MOD_OF  # reuse the module map

OUT = os.path.dirname(os.path.abspath(__file__))
rr = te.loada(build_model_v44())
sp = rr.getFloatingSpeciesIds()
rr.reset(); rr.simulate(0, 3000, 3000)

# collect diagonal timescales across a few phases (min); take the SLOWEST (max tau) each species reaches
taus = {s: [] for s in sp}
eig_all = []
t = 3000.0
for k in range(8):
    rr.simulate(t, t + 175, 20); t += 175
    J = np.array(rr.getFullJacobian())
    for i, s in enumerate(sp):
        d = abs(J[i, i])
        taus[s].append(1.0 / d if d > 1e-9 else np.inf)
    ev = np.linalg.eigvals(J)
    eig_all.append(ev)

# per-species representative tau = MEDIAN over phases (min). Convert to hours.
tau_h = {s: float(np.median([x for x in taus[s] if np.isfinite(x)]) / 60.0)
         if any(np.isfinite(x) for x in taus[s]) else np.inf for s in sp}
ranked = sorted(sp, key=lambda s: tau_h[s], reverse=True)

print("=== SPECIES RELAXATION TIMESCALES (tau = 1/|J_ii|, hours; slowest first) ===")
print(f"{'species':>12} {'tau(h)':>9}  module")
for s in ranked:
    th = tau_h[s]
    tag = f"{th:8.2f}" if np.isfinite(th) and th < 1e4 else "   slow*"
    print(f"{s:>12} {tag}  {MOD_OF.get(s,'?')}")

# fast/slow split: find the biggest log-gap in the RELEVANT band (0.1-40 h) -- the minutes/hours boundary,
# not the spurious gap among the sub-second rapid-equilibrium species.
band = sorted([(s, tau_h[s]) for s in sp if np.isfinite(tau_h[s]) and 0.1 < tau_h[s] < 40], key=lambda x: x[1])
logs = np.log10([v for _, v in band])
gaps = np.diff(logs)
cut_i = int(np.argmax(gaps)) if len(gaps) else 0
cut_tau = float(np.sqrt(band[cut_i][1] * band[cut_i + 1][1])) if len(band) > cut_i + 1 else 1.0
slow = [s for s in sp if np.isfinite(tau_h[s]) and tau_h[s] >= cut_tau] + \
       [s for s in sp if not np.isfinite(tau_h[s])]   # Dna (inf) is slow
print(f"\nbiggest spectral gap at tau ~ {cut_tau:.2f} h  ->  {len(slow)} SLOW species")
print("SLOW (real dof):", sorted(slow, key=lambda s: -tau_h[s]))
print("slow modules:", sorted(set(MOD_OF.get(s,'?') for s in slow)))

# full eigenvalue timescales (mixed modes) -> the true dynamical spectrum
ev = eig_all[3]
re = np.real(ev)
mode_tau = np.sort(1.0 / (np.abs(re[np.abs(re) > 1e-9]) + 1e-12))[::-1] / 60.0
print(f"\nslowest 8 eigen-mode timescales (h): {np.round(mode_tau[:8],2)}")
print(f"fastest 5 eigen-mode timescales (h): {np.round(mode_tau[-5:],4)}")

# ---- figure: timescale bar (log), colored fast/slow ----
fig, ax = plt.subplots(figsize=(10, 11))
ys = [s for s in ranked if np.isfinite(tau_h[s])]
vals = [tau_h[s] for s in ys]
cols = ['#b34a43' if tau_h[s] > cut_tau else '#12756c' for s in ys]
ax.barh(range(len(ys)), vals, color=cols, edgecolor='none')
ax.set_yticks(range(len(ys))); ax.set_yticklabels([f"{s}  ·{MOD_OF.get(s,'?')}" for s in ys], fontsize=7)
ax.set_xscale('log'); ax.invert_yaxis()
ax.axvline(cut_tau, color='k', ls='--', lw=1, alpha=.6)
ax.set_xlabel('relaxation timescale  tau = 1/|J_ii|   (hours, log)')
ax.set_title('v44 species timescales — SLOW drivers (red) vs FAST slaved (teal)\nfast/slow cut at the largest spectral gap', fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(OUT, 'B_timescales.png'), dpi=140, bbox_inches='tight')
print("wrote B_timescales.png")

json.dump(dict(tau_hours={s: (tau_h[s] if np.isfinite(tau_h[s]) else None) for s in sp},
               ranked=ranked, cut_tau_h=float(cut_tau), slow=sorted(slow, key=lambda s: -tau_h[s]),
               slow_modules=sorted(set(MOD_OF.get(s,'?') for s in slow)),
               slow_eigen_tau_h=[float(x) for x in mode_tau[:8]]),
          open(os.path.join(OUT, 'B_timescales.json'), 'w'), indent=2)
print("wrote B_timescales.json")
