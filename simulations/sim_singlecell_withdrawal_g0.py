"""SINGLE-CELL ensemble through Hh withdrawal — heterogeneous TEMPORAL properties.

The cycle-averaged ramp figures smooth over two things real cells show: (1) the inter-division interval
LENGTHENING as Hh is withdrawn, and (2) transient G0 arrests (p27-high pre-S dwells) that recur before a
cell finally arrests. Here we simulate actual individual cells, each with different temporal parameters,
and read those out — WITH vs WITHOUT the H3K27me3 mark on CyclinD1.

Per-cell heterogeneity (drawn once, PAIRED across WITH/WITHOUT — same cells, only f0_mk differs):
  mu        cell-cycle period / growth rate   (THE temporal axis)      lognormal
  del_mk    H3K27me3 mark turnover                                     lognormal
  kDeEZ     EZH2 protein half-life                                     lognormal
  P21_div   inherited (birth) p27  — G0 propensity                     lognormal (Spencer/Fan-Meyer axis)
  k_Cd_tl   CyclinD1 setpoint                                          lognormal

Protocol: settle cycling at SHH=1.0 → withdraw to SHH=0.25 over 100 h → hold 150 h.
Readouts (probe_transient_g0 conventions): divisions = MPF peaks; G0 = pre-S ∧ p27>0.1; a G0 run ≥2 h that
is later followed by a division = TRANSIENT G0, one that never resolves = ARREST.

Run:  ./venv/bin/python simulations/sim_singlecell_withdrawal_g0.py [--pilot] [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.signal import find_peaks
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

PILOT = '--pilot' in sys.argv
FRESH = '--fresh' in sys.argv
N = 6 if PILOT else 28
CACHE = 'simulations/sim_singlecell_withdrawal_g0_cache.npz'
SHH_HI, SHH_LO = 1.0, 0.35                                        # partial withdrawal (effect maximal here)
SETTLE, D_WD, HOLD, CHUNK, NP = 6000.0, 6000.0, 12000.0, 120.0, 11
DWELL_CUT, P27_HI = 2.0, 0.1
SEL = ['time', 'SHH', 'MPF', 'P21', 'aRc', 'Dna', 'Cd']

rng = np.random.default_rng(11)
mu_i   = np.clip(0.0005 * np.exp(rng.normal(0, 0.22, N)), 0.00033, 0.00075)
delmk_i = np.clip(0.0007 * np.exp(rng.normal(0, 0.30, N)), 0.0003, 0.0016)
kdeez_i = np.clip(0.00015 * np.exp(rng.normal(0, 0.30, N)), 0.00007, 0.00033)
p21d_i = np.clip(np.exp(rng.normal(np.log(0.45), 0.55, N)), 0.12, 2.0)
ktl_i  = np.clip(0.801 * np.exp(rng.normal(0, 0.33, N)), 0.35, 1.6)

_RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
_RR.integrator.setValue('absolute_tolerance', 1e-9); _RR.integrator.setValue('relative_tolerance', 1e-6)
try: _RR.integrator.setValue('maximum_num_steps', 200000)
except Exception: pass


def shh_of(t):                                                    # t in post-settle frame
    if t <= 0: return SHH_HI
    return SHH_HI + (SHH_LO - SHH_HI) * min(t / D_WD, 1.0)


def run_cell(i, f0):
    _RR.reset()
    _RR['SHH'] = SHH_HI; _RR['MYCN_amplification'] = 1.0; _RR['Ptch1_copy_number'] = 1.0
    _RR['p16'] = 0.0; _RR['p18'] = 0.464; _RR['kSyP21'] = 0.002; _RR['EZH2i'] = 0
    _RR['mu'] = mu_i[i]; _RR['del_mk'] = delmk_i[i]; _RR['kDeEZ'] = kdeez_i[i]
    _RR['P21_div'] = p21d_i[i]; _RR['k_Cd_translation'] = ktl_i[i]; _RR['f0_mk'] = f0
    try:
        _RR.simulate(0, SETTLE, 50, selections=['time'])         # settle cycling at high Hh
    except Exception:
        pass
    t = 0.0; rows = {s: [] for s in SEL}
    while t < D_WD + HOLD:
        _RR['SHH'] = float(shh_of(t)); r = None
        for atol in (1e-9, 1e-8, 1e-7):
            try:
                _RR.integrator.setValue('absolute_tolerance', atol)
                r = _RR.simulate(SETTLE + t, SETTLE + t + CHUNK, NP, selections=SEL); break
            except Exception:
                r = None
        if r is None: break
        for s in SEL: rows[s].append(r[s][1:])
        t += CHUNK
    if not rows['time']: return None
    for s in SEL: rows[s] = np.concatenate(rows[s])
    rows['time'] = rows['time'] - SETTLE                          # withdrawal frame (0 = withdrawal start)
    return rows


def analyze(T, MPF, P21, aRc, Dna):
    dt = np.median(np.diff(T))
    pk, _ = find_peaks(MPF, prominence=0.15, distance=max(1, int(200 / dt)))
    pt = T[pk]
    inS = (aRc > 0.05) & (Dna < 0.98); inG2 = Dna >= 0.98; preS = ~inS & ~inG2
    g0 = (preS & (P21 > P27_HI)).astype(np.int8)
    state = np.zeros(len(T), np.int8)                            # 0 cycling, 1 transient-G0, 2 arrested
    d = np.diff(np.concatenate([[0], g0, [0]]))
    for s, e in zip(np.where(d == 1)[0], np.where(d == -1)[0]):
        ee = min(e, len(T) - 1)
        if T[ee] - T[s] < DWELL_CUT * 60: continue               # normal short G1, not a G0 excursion
        if np.any(pt > T[ee]):
            state[s:e] = 1
        else:
            state[s:] = 2; break
    inter = np.diff(pt) / 60.0; mids = (pt[:-1] + pt[1:]) / 2 if len(pt) > 1 else np.array([])
    arrest_t = float(T[np.argmax(state == 2)]) if (state == 2).any() else np.nan
    return pt, inter, mids, state, arrest_t


T_REF = np.arange(0.0, D_WD + HOLD + 1e-6, 12.0)                 # fixed 12-min withdrawal-frame grid
L = len(T_REF); EXSET = {0, N // 2, N - 1}
cat = lambda xs: np.concatenate(xs) if xs else np.array([])

if FRESH or not os.path.exists(CACHE):
    store = {}
    for cond, f0 in [('W', 0.233), ('N', 1.0)]:
        state_rows = np.full((N, L), 2, np.int8)                 # default = arrested (covers cells that die at settle)
        cd_rows = np.full((N, L), np.nan); p21_rows = np.full((N, L), np.nan)
        inters = []; mids = []; arrest = np.full(N, np.nan); ndiv = np.zeros(N, int)
        dft = []; dfc = []
        for i in range(N):
            r = run_cell(i, f0)
            if r is None:
                arrest[i] = 0.0; print(f'  {cond} cell {i+1}/{N}  FAILED-at-settle -> arrested', flush=True); continue
            pt, inter, mid, state, at = analyze(r['time'], r['MPF'], r['P21'], r['aRc'], r['Dna'])
            state_rows[i] = np.clip(np.round(np.interp(T_REF, r['time'], state)), 0, 2).astype(np.int8)  # interp holds endpoints
            cd_rows[i] = np.interp(T_REF, r['time'], r['Cd']); p21_rows[i] = np.interp(T_REF, r['time'], r['P21'])
            inters.append(inter); mids.append(mid); arrest[i] = at; ndiv[i] = len(pt)
            dft.append(pt); dfc.append(np.full(len(pt), i))
            print(f'  {cond} cell {i+1}/{N}  arrest_t={(at/60 if np.isfinite(at) else float("nan")):.0f}h  ndiv={len(pt)}', flush=True)
        store[f'{cond}_state'] = state_rows; store[f'{cond}_cd'] = cd_rows; store[f'{cond}_p21'] = p21_rows
        store[f'{cond}_inter'] = cat(inters); store[f'{cond}_mids'] = cat(mids)
        store[f'{cond}_arrest'] = arrest; store[f'{cond}_ndiv'] = ndiv
        store[f'{cond}_divflat_t'] = cat(dft); store[f'{cond}_divflat_cell'] = cat(dfc)
    np.savez(CACHE, T_ref=T_REF, mu=mu_i, delmk=delmk_i, kdeez=kdeez_i, p21d=p21d_i, ktl=ktl_i, N=N,
             D_WD=D_WD, SHH_HI=SHH_HI, SHH_LO=SHH_LO, **store)
    print('cached ->', CACHE, flush=True)

z = np.load(CACHE, allow_pickle=True)
T = z['T_ref'] / 60.0                                            # h (withdrawal frame)
Nc = int(z['N']); D_WD = float(z['D_WD'])
shh_t = np.where(T * 60 <= 0, SHH_HI, SHH_HI + (SHH_LO - SHH_HI) * np.minimum(T * 60 / D_WD, 1.0))

fig = plt.figure(figsize=(17, 11))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1])
axR = fig.add_subplot(gs[0, 0]); axE = fig.add_subplot(gs[0, 1])
axI = fig.add_subplot(gs[1, 0]); axF = fig.add_subplot(gs[1, 1])

# (A) raster WITH mark — order cells by arrest time (latest at top)
St = z['W_state']; order = np.argsort(np.nan_to_num(z['W_arrest'], nan=1e9))
cmap = ListedColormap(['#2e8b57', '#e08a1e', '#8b1a1a'])         # cycling / transient-G0 / arrested
axR.imshow(St[order], aspect='auto', cmap=cmap, vmin=0, vmax=2,
           extent=[T.min(), T.max(), 0, Nc], interpolation='nearest')
dt_, dc_ = z['W_divflat_t'] / 60.0, z['W_divflat_cell']
pos = {c: r for r, c in enumerate(order)}
axR.plot(dt_, [Nc - pos[int(c)] - 0.5 for c in dc_], '|', color='k', ms=4, mew=0.7, alpha=0.5)
axR.axvline(D_WD / 60, color='k', ls=':', lw=1); axR.axvline(0, color='k', lw=0.6)
ax2 = axR.twinx(); ax2.plot(T, shh_t, color='#2b6cb0', lw=2); ax2.set_ylabel('SHH', color='#2b6cb0'); ax2.set_ylim(0, 1.15)
axR.set_ylabel('cells (sorted by arrest)'); axR.set_xlabel('time since withdrawal start (h)')
axR.set_title('(A) WITH mark — single-cell states through Hh withdrawal\ngreen cycling · orange TRANSIENT G0 · dark ARREST · ticks=divisions', fontweight='bold', fontsize=11)

# (B) example single-cell traces (Cd + p27) — period lengthening + G0
ex = list(np.argsort(z['W_ndiv'])[::-1][:3]); colors = ['#1f77b4', '#9467bd', '#2ca02c']
cdM, p21M = z['W_cd'], z['W_p21']
for k, i in enumerate(ex):
    axE.plot(T, cdM[i], color=colors[k], lw=1.4, label=f'cell {i} CyclinD1 (μ={z["mu"][i]*1e4:.1f}e-4, {int(z["W_ndiv"][i])} div)')
axE.set_ylabel('CyclinD1'); axE.set_xlabel('time since withdrawal start (h)')
axp = axE.twinx()
for k, i in enumerate(ex):
    axp.plot(T, p21M[i], color=colors[k], lw=1.0, ls=':', alpha=0.85)
axp.axhline(P27_HI, color='#999', ls='--', lw=0.8); axp.set_ylabel('p27 (dotted; dashed=G0 cut)', color='#666')
axE.axvline(D_WD / 60, color='k', ls=':', lw=1); axE.set_xlim(0, 200)
axE.set_title('(B) Example cells (WITH mark): CyclinD1 oscillations slow & shrink;\np27 (dotted) rises into G0 as Hh withdraws', fontweight='bold', fontsize=11)
axE.legend(fontsize=8, loc='upper right')

# (C) inter-division interval vs time — WITH vs WITHOUT
def binned(mids, inter, edges):
    idx = np.digitize(mids, edges); m = np.array([np.median(inter[idx == k]) if (idx == k).any() else np.nan for k in range(1, len(edges))])
    return 0.5 * (edges[:-1] + edges[1:]), m
edges = np.linspace(0, T.max(), 11)
for cond, col, lab in [('W', '#8b1a1a', 'WITH mark'), ('N', '#666', 'WITHOUT mark')]:
    mid = z[f'{cond}_mids'] / 60.0; iv = z[f'{cond}_inter']
    axI.plot(mid, iv, '.', color=col, ms=3, alpha=0.28)
    bx, bm = binned(mid, iv, edges); axI.plot(bx, bm, '-o', color=col, lw=2.4, ms=5, label=lab)
axI.axvline(D_WD / 60, color='k', ls=':', lw=1)
axI.set_xlabel('time since withdrawal start (h)'); axI.set_ylabel('inter-division interval (h)')
axI.set_title('(C) Cell-cycle period LENGTHENS through withdrawal\n(median per time-bin; points = individual cycles)', fontweight='bold', fontsize=11)
axI.legend(fontsize=9); axI.grid(alpha=0.15)

# (D) fraction of cells transient-G0 and arrested vs time — WITH vs WITHOUT
for cond, col, lab in [('W', '#8b1a1a', 'WITH mark'), ('N', '#666', 'WITHOUT mark')]:
    St = z[f'{cond}_state']
    frac_g0 = (St == 1).mean(0); frac_ar = (St == 2).mean(0)
    axF.plot(T, frac_ar, color=col, lw=2.6, label=f'{lab} — arrested')
    axF.plot(T, frac_g0, color=col, lw=1.6, ls='--', label=f'{lab} — transient G0')
axF.axvline(D_WD / 60, color='k', ls=':', lw=1)
axF.set_xlabel('time since withdrawal start (h)'); axF.set_ylabel('fraction of cells')
axF.set_title('(D) G0 & arrest through withdrawal — WITH vs WITHOUT mark', fontweight='bold', fontsize=11)
axF.legend(fontsize=8.5, loc='center right'); axF.grid(alpha=0.15); axF.set_ylim(0, 1)

fig.suptitle(f'Single-cell ensemble (N={Nc}, heterogeneous temporal properties) through Hh withdrawal — transient G0 & period lengthening',
             fontsize=13.5, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/sim_singlecell_withdrawal_g0.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/sim_singlecell_withdrawal_g0.pdf', bbox_inches='tight')
plt.close()
print('WITH   median ndiv=%.1f  arrested@end=%.0f%%' % (np.median(z['W_ndiv']), 100 * (z['W_state'][:, -1] == 2).mean()))
print('WITHOUT median ndiv=%.1f  arrested@end=%.0f%%' % (np.median(z['N_ndiv']), 100 * (z['N_state'][:, -1] == 2).mean()))
print('Saved sim_singlecell_withdrawal_g0.png')
