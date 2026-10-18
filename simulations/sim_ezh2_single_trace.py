"""Single-cell EZH2 mRNA vs protein over time (JP 2026-10-09). Illustrates the bursting mechanism:
the EZH2 gene fires in telegraph bursts (on/off), so mRNA (t1/2 ~1.2 h) tracks the bursts closely and
is spiky, while the STABLE protein (t1/2 ~42 h) integrates them into slow drift. One representative DMSO
cell (mult=1, fExt=1), fine time resolution. Records EZH2m, EZH2, telegraph g, division sawtooth.
Run:  venv/bin/python simulations/sim_ezh2_single_trace.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

DT = 15.0                 # min per chunk (fine, to resolve the ~1.2 h mRNA decay)
WARM = 5400.0             # 90 h warmup (start near stationary)
REC = 9000.0             # 150 h of recorded trace
TAU_BURST = 200.0; DUTY = 0.62
K_ON = DUTY / TAU_BURST; K_OFF = (1 - DUTY) / TAU_BURST
SEED = int(os.environ.get('TRACE_SEED', '7'))
SEL = ['EZH2m', 'EZH2', 'Dna', 'pRb']

rr = te.loada(build_model_v44())
rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
rr.integrator.setValue('maximum_num_steps', 12000000)
B0 = {p: rr[p] for p in ['kEZbas', 'kEZbas_Cd', 'kEZE2f']}
rng = np.random.default_rng(SEED)

rr.reset()
rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
except Exception: pass
rr['EZH2'] = 6.0; rr['EZH2m'] = 0.3          # representative cell: mult=1, fExt=1
g = 1 if rng.random() < DUTY else 0

t = 0.0; T_END = WARM + REC
T, MRNA, PROT, DNA, GON = [], [], [], [], []
while t < T_END:
    g = 0 if (g and rng.random() < K_OFF * DT) else (1 if (not g and rng.random() < K_ON * DT) else g)
    eta = g / DUTY                            # telegraph only (fExt=1); mean(eta)=1
    for p in B0: rr[p] = B0[p] * eta
    r = rr.simulate(t, t + DT, 2, selections=SEL)
    if t >= WARM:
        T.append((t - WARM) / 60.0)           # hours since recording start
        MRNA.append(float(r['EZH2m'][-1])); PROT.append(float(r['EZH2'][-1]))
        DNA.append(float(r['Dna'][-1])); GON.append(g)
    t += DT

T = np.array(T); MRNA = np.array(MRNA); PROT = np.array(PROT); GON = np.array(GON)
cv_m = MRNA.std() / MRNA.mean(); cv_p = PROT.std() / PROT.mean()

# telegraph ON spans for shading
spans = []; on0 = None
for k in range(len(T)):
    if GON[k] and on0 is None: on0 = T[k]
    if (not GON[k]) and on0 is not None: spans.append((on0, T[k])); on0 = None
if on0 is not None: spans.append((on0, T[-1]))

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 6.6), sharex=True)
fig.suptitle("Single cell: EZH2 mRNA (fast, bursty) vs protein (stable, slow drift)\n"
             "transcription fires in telegraph bursts (shaded = gene ON); the stable protein low-pass-filters them",
             fontsize=12, fontweight='bold')
for (a, b) in spans:
    ax1.axvspan(a, b, color='#f0c45a', alpha=0.22, lw=0)
    ax2.axvspan(a, b, color='#f0c45a', alpha=0.22, lw=0)
ax1.plot(T, MRNA, color='#c0392b', lw=1.3)
ax1.set_ylabel('EZH2 mRNA (a.u.)')
ax1.set_title(f"mRNA  —  t1/2 ~1.2 h  →  tracks the bursts  (temporal CV {cv_m:.2f})", fontsize=10, loc='left')
ax2.plot(T, PROT, color='#1f6f54', lw=1.8)
ax2.set_ylabel('EZH2 protein (a.u.)'); ax2.set_xlabel('time (h)')
ax2.set_title(f"protein  —  t1/2 ~42 h  →  integrates mRNA into slow drift + division sawtooth  "
              f"(temporal CV {cv_p:.2f})", fontsize=10, loc='left')
plt.tight_layout(rect=[0, 0, 1, 0.93])
out = "/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/79019eb0-7c15-4650-b7d3-81347ccd6a13/scratchpad/fig_ezh2_single_trace.png"
plt.savefig(out, dpi=130); print('saved', out)
print(f"seed={SEED}  mRNA: mean {MRNA.mean():.3f} CV {cv_m:.2f} | protein: mean {PROT.mean():.3f} CV {cv_p:.2f} | n_bursts {len(spans)}")
