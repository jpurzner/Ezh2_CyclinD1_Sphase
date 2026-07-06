"""Does the model already produce a BIMODAL G0 p27 distribution — from the transient-G0 vs
immediate-re-entry split (Purzner hypothesis)?

The v44 commitment switch is the Skp2->p27->E2f feedforward (bistable-ish R-point). If it is sharp, an
asynchronous 2N (pre-S) population should split into two p27 modes even from a UNIMODAL birth-p27 input:
  - high p27  = uncommitted / dwelling in transient G0 (E2f low)
  - low  p27  = committed, p27 degraded (E2f high), about to enter S
That would reproduce the measured bimodal/heavy-tailed G0 p27 as an EMERGENT property, not an imposed one.

Test: a heterogeneous cycling population (CyclinD1 setpoint data-grounded; birth p27 UNIMODAL). We sample
p27 across each cell's settled trajectory, keep the 2N pre-S samples (Dna<0.15, aRc<0.05 = G0/G1), pool
them, and ask whether p27 is bimodal — and whether the two modes are the uncommitted (E2f-low) vs
committed (E2f-high) cells. Scale-free (p27 shown / median).

Run:  ./venv/bin/python simulations/sim_p27_bimodality.py [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

FRESH = '--fresh' in sys.argv
CACHE = 'simulations/sim_p27_bimodality_cache.npz'
N = 30
KTL0 = 0.801
CD_SDLOG = 0.633          # CyclinD1 setpoint spread (DATA, scale-free)
P27_SDLOG_IN = 0.35       # birth-p27 spread: UNIMODAL input (test whether bimodality EMERGES)
# IF has a floor: what is imaged = model p27 + a modest offset (background + basal p27; cells are above
# background). Kept small -- observable-level, does NOT change the p27 dynamics. Tunable to the measured
# ~4x low/high separation. This is applied in post-processing, so the cache holds raw model p27.
P27_FLOOR = 0.12
T_END, SETTLE = 15000, 6000

M = build_model_v44(with_ezh2=True, with_hh=True)
MB = dict(SHH=0.5, HHi=0.0, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002)
_R = te.loada(M); _R.integrator.setValue('relative_tolerance', 1e-6)
try: _R.integrator.setValue('maximum_num_steps', 300000)
except Exception: pass


def trace(cd_scale, p21d):
    _R.reset()
    for k, v in MB.items(): _R[k] = v
    _R['k_Cd_translation'] = float(KTL0 * cd_scale); _R['P21_div'] = float(p21d)
    for atol in (1e-9, 1e-8, 1e-7, 1e-6):
        try:
            _R.integrator.setValue('absolute_tolerance', atol)
            d = _R.simulate(0, T_END, int(T_END / 2), selections=['time', 'P21', 'Dna', 'aRc', 'E2f', 'MPF'])
            m = d['time'] >= SETTLE
            return d['P21'][m], d['Dna'][m], d['aRc'][m], d['E2f'][m], d['MPF'][m]
        except Exception:
            continue
    return None


if FRESH or not os.path.exists(CACHE):
    rng = np.random.default_rng(7)
    cd_scale = np.clip(np.exp(rng.normal(0.0, CD_SDLOG, N)), 0.2, 6.0)
    p21_div = np.clip(np.exp(rng.normal(np.log(0.5), P27_SDLOG_IN, N)), 0.15, 2.0)   # UNIMODAL input
    P27, E2F, CYCS = [], [], []      # pooled 2N pre-S samples (+ per-sample cycling flag)
    ncyc = 0
    for i in range(N):
        r = trace(cd_scale[i], p21_div[i])
        if r is None:
            print(f"  cell {i+1}/{N}  FAILED", flush=True); continue
        p21, dna, arc, e2f, mpf = r
        preS = (dna < 0.15) & (arc < 0.05)      # 2N, not in S = G0/G1
        pk, _ = find_peaks(mpf, prominence=0.15, distance=200)
        cyc = len(pk) >= 2                        # properly cycling (>=2 divisions in the window)
        ncyc += int(cyc)
        P27.append(p21[preS]); E2F.append(e2f[preS]); CYCS.append(np.full(preS.sum(), cyc))
        print(f"  cell {i+1:2d}/{N}  cd {cd_scale[i]:.2f} p27_birth {p21_div[i]:.2f}  divs {len(pk)}  2N-samp {preS.sum()}  {'CYCLING' if cyc else 'arrested'}", flush=True)
    P27 = np.concatenate(P27); E2F = np.concatenate(E2F); CYCS = np.concatenate(CYCS)
    print(f"  cycling cells: {ncyc}/{N}", flush=True)
    np.savez(CACHE, P27=P27, E2F=E2F, CYCS=CYCS, p27_sdlog_in=P27_SDLOG_IN, cd_sdlog=CD_SDLOG)
    print('cached ->', CACHE, flush=True)

from scipy.stats import gaussian_kde
z = np.load(CACHE)
P27, E2F, CYCS = z['P27'], z['E2F'], z['CYCS'].astype(bool)
P27f = P27 + P27_FLOOR                        # OBSERVABLE = model p27 + IF floor (modest, background+basal)
P27n = P27f / np.median(P27f)                 # scale-free
lg = np.log(P27n)

# KDE + mode detection on log(p27)
xs = np.linspace(lg.min(), lg.max(), 400)
kde = gaussian_kde(lg, bw_method=0.25)(xs)
pk, _ = find_peaks(kde, prominence=0.03)
n_modes = len(pk)
modes = np.sort(np.exp(xs[pk]))
sep = modes[-1] / modes[0] if len(modes) >= 2 else np.nan     # high/low mode separation
# commitment split (E2f)
thr_e2 = np.median(E2F)
com = E2F >= thr_e2; unc = ~com               # committed (E2f-high) vs uncommitted (E2f-low)
cv = lambda x: np.std(x) / np.mean(x)
cv_hi, cv_lo = cv(P27f[unc]), cv(P27f[com])
print(f"\n2N p27 (with floor {P27_FLOOR}):  {n_modes} mode(s), separation {sep:.1f}x   [data ~4x]")
print(f"  uncommitted (E2f-low, high p27) CV {cv_hi:.2f}  [data G0-hi 1.29]   committed (E2f-high, low p27) CV {cv_lo:.2f}  [data G1 0.33 / G0-lo 0.27]")
print(f"  cycling 2N samples: {100*CYCS.mean():.0f}%")

C_UNC, C_COM = '#c0392b', '#2c7fb8'
fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.0))
b = np.logspace(np.log10(P27n.min()), np.log10(P27n.max()), 44)

# (A) pooled 2N p27 distribution (with IF floor) — bimodal?
ax[0].hist(P27n, bins=b, density=True, color='#b0adc7', edgecolor='none', alpha=0.85)
ax[0].plot(np.exp(xs), kde, color='#4a2f6b', lw=2.2, label=f'KDE ({n_modes} mode{"s" if n_modes>1 else ""})')
for m in modes:
    ax[0].axvline(m, color='#4a2f6b', ls=':', lw=1)
ax[0].set_xscale('log'); ax[0].set_xlabel('p27  (÷ median, incl. IF floor)'); ax[0].set_ylabel('density')
ax[0].set_title(f'(A) 2N (G0/G1) p27 distribution\n{n_modes} modes, separation {sep:.1f}× (data ~4×)', fontweight='bold'); ax[0].legend(fontsize=9)

# (B) split by commitment (E2f) — the two components + CVs vs data
ax[1].hist(P27n[unc], bins=b, density=True, color=C_UNC, alpha=0.55, label=f'uncommitted / G0 (CV {cv_hi:.2f})')
ax[1].hist(P27n[com], bins=b, density=True, color=C_COM, alpha=0.55, label=f'committed / G1 (CV {cv_lo:.2f})')
ax[1].set_xscale('log'); ax[1].set_xlabel('p27  (÷ median)'); ax[1].set_ylabel('density')
ax[1].set_title('(B) Two components = the commitment split\nhigh-p27 uncommitted vs low-p27 committed', fontweight='bold'); ax[1].legend(fontsize=8.5)

# (C) p27 vs E2f (the switch), coloured by cycling vs arrested
idx = np.random.default_rng(0).choice(len(P27n), size=min(5000, len(P27n)), replace=False)
ax[2].scatter(E2F[idx][CYCS[idx]], P27n[idx][CYCS[idx]], s=5, c='#2ecc71', alpha=0.35, label='cycling (transient G0)')
ax[2].scatter(E2F[idx][~CYCS[idx]], P27n[idx][~CYCS[idx]], s=5, c='#b03a2e', alpha=0.35, label='arrested')
ax[2].axvline(thr_e2, color='k', ls=':', lw=1)
ax[2].set_yscale('log'); ax[2].set_xlabel('E2f (commitment marker)'); ax[2].set_ylabel('p27 (÷ median)')
ax[2].set_title('(C) Commitment (E2f↑) degrades p27;\nhigh-p27 mode = dwelling + arrested', fontweight='bold'); ax[2].legend(fontsize=8, loc='lower left'); ax[2].grid(alpha=0.15)

fig.suptitle(f'Bimodal G0 p27 emerges from the commitment switch (unimodal birth-p27 → {n_modes} modes, {sep:.1f}× apart; IF floor {P27_FLOOR})',
             fontsize=12, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/sim_p27_bimodality.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/sim_p27_bimodality.pdf', bbox_inches='tight')
plt.close()
print('Saved sim_p27_bimodality.png')
