"""SINGLE-CELL ensemble through Hh withdrawal — heterogeneous TEMPORAL properties (PARALLEL, hundreds–thousands).

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

Protocol: settle cycling at SHH=1.0 → withdraw to SHH_LO over 100 h → hold 200 h.
Readouts (probe_transient_g0 conventions): divisions = MPF peaks; G0 = pre-S ∧ p27>0.1; a G0 run ≥2 h that
is later followed by a division = TRANSIENT G0, one that never resolves = ARREST.

Run:  ./venv/bin/python simulations/sim_singlecell_withdrawal_g0.py [--n=800] [--shallow] [--workers=9] [--fresh]
      (--pilot => n=6, serial-ish; --shallow => withdraw to 0.50 not 0.35)
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
import multiprocessing as mp
from src.build_model_v44_heldt import build_model_v44

SEL = ['time', 'SHH', 'MPF', 'P21', 'aRc', 'Dna', 'Cd']
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, EZH2i=0)
MAXK = 10

# --- worker globals (populated by _init_worker in each process) ---
_RR = None
SHH_HI = SHH_LO = SETTLE = D_WD = HOLD = CHUNK = NP = DWELL_CUT = P27_HI = None
T_REF = None


def shh_of(t):                                                    # t in post-settle (withdrawal) frame
    if t <= 0: return SHH_HI
    return SHH_HI + (SHH_LO - SHH_HI) * min(t / D_WD, 1.0)


def run_cell(pars, f0):
    mu, delmk, kdeez, p21d, ktl = pars
    _RR.reset()
    _RR['SHH'] = SHH_HI
    for k, v in GNP.items(): _RR[k] = v
    _RR['mu'] = mu; _RR['del_mk'] = delmk; _RR['kDeEZ'] = kdeez
    _RR['P21_div'] = p21d; _RR['k_Cd_translation'] = ktl; _RR['f0_mk'] = f0
    try:
        _RR.simulate(0, SETTLE, 50, selections=['time'])          # settle cycling at high Hh
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
    rows['time'] = rows['time'] - SETTLE
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
        if T[ee] - T[s] < DWELL_CUT * 60: continue
        if np.any(pt > T[ee]):
            state[s:e] = 1
        else:
            state[s:] = 2; break
    inter = np.diff(pt) / 60.0
    at = float(T[np.argmax(state == 2)]) if (state == 2).any() else np.nan
    return pt, inter, state, at


def _init_worker(cfg):
    global _RR, SHH_HI, SHH_LO, SETTLE, D_WD, HOLD, CHUNK, NP, DWELL_CUT, P27_HI, T_REF
    SHH_HI = cfg['SHH_HI']; SHH_LO = cfg['SHH_LO']; SETTLE = cfg['SETTLE']; D_WD = cfg['D_WD']
    HOLD = cfg['HOLD']; CHUNK = cfg['CHUNK']; NP = cfg['NP']; DWELL_CUT = cfg['DWELL_CUT']
    P27_HI = cfg['P27_HI']; T_REF = cfg['T_REF']
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('absolute_tolerance', 1e-9); _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 200000)
    except Exception: pass


def _work(task):
    i, f0, mu, delmk, kdeez, p21d, ktl = task
    r = run_cell((mu, delmk, kdeez, p21d, ktl), f0)
    z32 = np.array([], np.float32)
    if r is None:
        return (i, f0, None, None, None, z32, 0.0, 0, z32)
    pt, inter, state, at = analyze(r['time'], r['MPF'], r['P21'], r['aRc'], r['Dna'])
    srow = np.clip(np.round(np.interp(T_REF, r['time'], state)), 0, 2).astype(np.int8)
    cd = np.interp(T_REF, r['time'], r['Cd']).astype(np.float32)
    p21 = np.interp(T_REF, r['time'], r['P21']).astype(np.float32)
    return (i, f0, srow, cd, p21, inter.astype(np.float32), at, len(pt), pt.astype(np.float32))


def _arg(flag, default):
    for a in sys.argv:
        if a.startswith(flag + '='): return a.split('=', 1)[1]
    return default


if __name__ == '__main__':
    PILOT = '--pilot' in sys.argv
    FRESH = '--fresh' in sys.argv
    SHALLOW = '--shallow' in sys.argv
    N = int(_arg('--n', 6 if PILOT else 800))
    WORKERS = int(_arg('--workers', max(1, mp.cpu_count() - 1)))
    SUFFIX = '_shallow' if SHALLOW else ''
    CACHE = f'simulations/sim_singlecell_withdrawal_g0_cache{SUFFIX}.npz'
    OUT = f'simulations/sim_singlecell_withdrawal_g0{SUFFIX}'
    cfg = dict(SHH_HI=1.0, SHH_LO=(0.50 if SHALLOW else 0.35), SETTLE=6000.0, D_WD=6000.0, HOLD=12000.0,
               CHUNK=120.0, NP=11, DWELL_CUT=2.0, P27_HI=0.1)
    cfg['T_REF'] = np.arange(0.0, cfg['D_WD'] + cfg['HOLD'] + 1e-6, 12.0)
    T_REF_M = cfg['T_REF']; L = len(T_REF_M)

    if FRESH or not os.path.exists(CACHE):
        rng = np.random.default_rng(11)
        mu_i = np.clip(0.0005 * np.exp(rng.normal(0, 0.22, N)), 0.00033, 0.00075)
        delmk_i = np.clip(0.0007 * np.exp(rng.normal(0, 0.30, N)), 0.0003, 0.0016)
        kdeez_i = np.clip(0.00015 * np.exp(rng.normal(0, 0.30, N)), 0.00007, 0.00033)
        p21d_i = np.clip(np.exp(rng.normal(np.log(0.45), 0.55, N)), 0.12, 2.0)
        ktl_i = np.clip(0.801 * np.exp(rng.normal(0, 0.33, N)), 0.35, 1.6)
        tasks = ([(i, 0.233, mu_i[i], delmk_i[i], kdeez_i[i], p21d_i[i], ktl_i[i]) for i in range(N)] +
                 [(i, 1.0, mu_i[i], delmk_i[i], kdeez_i[i], p21d_i[i], ktl_i[i]) for i in range(N)])
        print(f'running {2*N} cells on {WORKERS} workers ...', flush=True)
        ctx = mp.get_context('spawn')
        with ctx.Pool(WORKERS, initializer=_init_worker, initargs=(cfg,)) as pool:
            results = []
            for k, res in enumerate(pool.imap_unordered(_work, tasks, chunksize=4)):
                results.append(res)
                if (k + 1) % 100 == 0: print(f'  {k+1}/{2*N} done', flush=True)

        store = {}
        for f0, tag in [(0.233, 'W'), (1.0, 'N')]:
            rows = sorted([r for r in results if r[1] == f0], key=lambda r: r[0])
            state = np.full((N, L), 2, np.int8); cd = np.full((N, L), np.nan, np.float32); p21 = np.full((N, L), np.nan, np.float32)
            ivk = np.full((N, MAXK), np.nan); arrest = np.full(N, np.nan); ndiv = np.zeros(N, int)
            inters = []; dft = []; dfc = []
            for k, r in enumerate(rows):
                i, _f0, srow, cdr, p21r, inter, at, nd, pt = r
                if srow is None:
                    arrest[k] = 0.0; continue
                state[k] = srow; cd[k] = cdr; p21[k] = p21r; arrest[k] = at; ndiv[k] = nd
                ivk[k, :min(len(inter), MAXK)] = inter[:MAXK]
                inters.append(inter); dft.append(pt); dfc.append(np.full(len(pt), k))
            store[f'{tag}_state'] = state; store[f'{tag}_cd'] = cd; store[f'{tag}_p21'] = p21
            store[f'{tag}_ivk'] = ivk; store[f'{tag}_arrest'] = arrest; store[f'{tag}_ndiv'] = ndiv
            store[f'{tag}_divflat_t'] = np.concatenate(dft) if dft else np.array([])
            store[f'{tag}_divflat_cell'] = np.concatenate(dfc) if dfc else np.array([])
        np.savez(CACHE, T_ref=T_REF_M, N=N, D_WD=cfg['D_WD'], SHH_HI=cfg['SHH_HI'], SHH_LO=cfg['SHH_LO'], **store)
        print('cached ->', CACHE, flush=True)

    z = np.load(CACHE, allow_pickle=True)
    T = z['T_ref'] / 60.0; Nc = int(z['N']); D_WD = float(z['D_WD']); SHH_HI = float(z['SHH_HI']); SHH_LO = float(z['SHH_LO'])
    shh_t = np.where(T * 60 <= 0, SHH_HI, SHH_HI + (SHH_LO - SHH_HI) * np.minimum(T * 60 / D_WD, 1.0))

    fig = plt.figure(figsize=(17, 11))
    gs = fig.add_gridspec(2, 2)
    axR = fig.add_subplot(gs[0, 0]); axE = fig.add_subplot(gs[0, 1]); axI = fig.add_subplot(gs[1, 0]); axF = fig.add_subplot(gs[1, 1])

    # (A) raster WITH mark — cells ordered by arrest time
    St = z['W_state']; order = np.argsort(np.nan_to_num(z['W_arrest'], nan=1e9))
    cmap = ListedColormap(['#2e8b57', '#e08a1e', '#8b1a1a'])
    axR.imshow(St[order], aspect='auto', cmap=cmap, vmin=0, vmax=2, extent=[T.min(), T.max(), 0, Nc], interpolation='nearest')
    if Nc <= 60:
        dt_, dc_ = z['W_divflat_t'] / 60.0, z['W_divflat_cell']; pos = {c: r for r, c in enumerate(order)}
        axR.plot(dt_, [Nc - pos[int(c)] - 0.5 for c in dc_], '|', color='k', ms=4, mew=0.7, alpha=0.5)
    axR.axvline(D_WD / 60, color='k', ls=':', lw=1)
    ax2 = axR.twinx(); ax2.plot(T, shh_t, color='#2b6cb0', lw=2); ax2.set_ylabel('SHH', color='#2b6cb0'); ax2.set_ylim(0, 1.15)
    axR.set_ylabel('cells (sorted by arrest)'); axR.set_xlabel('time since withdrawal start (h)')
    axR.set_title(f'(A) WITH mark — {Nc} single-cell states through Hh withdrawal\ngreen cycling · orange TRANSIENT G0 · dark ARREST', fontweight='bold', fontsize=11)

    # (B) example cells (most divisions) — CyclinD1 + p27
    ex = list(np.argsort(z['W_ndiv'])[::-1][:3]); colors = ['#1f77b4', '#9467bd', '#2ca02c']
    cdM, p21M = z['W_cd'], z['W_p21']
    for k, i in enumerate(ex):
        axE.plot(T, cdM[i], color=colors[k], lw=1.4, label=f'cell {i} ({int(z["W_ndiv"][i])} div)')
    axp = axE.twinx()
    for k, i in enumerate(ex):
        axp.plot(T, p21M[i], color=colors[k], lw=1.0, ls=':', alpha=0.85)
    axp.axhline(float(z['SHH_LO']) * 0 + 0.1, color='#999', ls='--', lw=0.8); axp.set_ylabel('p27 (dotted; dashed=G0 cut)', color='#666')
    axE.axvline(D_WD / 60, color='k', ls=':', lw=1); axE.set_xlim(0, 200)
    axE.set_ylabel('CyclinD1'); axE.set_xlabel('time since withdrawal start (h)')
    axE.set_title('(B) Example cells (WITH mark): CyclinD1 oscillations slow & shrink;\np27 (dotted) rises into G0 as Hh withdraws', fontweight='bold', fontsize=11)
    axE.legend(fontsize=8, loc='upper right')

    # (C) interval ALIGNED TO EACH CELL'S ARREST — the slowdown into G0 (survivorship-free: only arrested cells,
    #     each aligned to its own last completed cycle). By-division-number-from-START is survivorship-confounded
    #     (only intrinsically-fast cells reach high k), so it spuriously DECREASES; aligning to arrest removes that.
    MK = z['W_ivk'].shape[1]; NB = 6                             # show up to 6 cycles before arrest
    for cond, col, lab in [('W', '#8b1a1a', 'WITH mark'), ('N', '#4a4a4a', 'WITHOUT mark')]:
        ivk = z[f'{cond}_ivk']; arr = z[f'{cond}_state'][:, -1] == 2
        M = np.full((int(arr.sum()), NB), np.nan)               # rows: arrested cells; col 0 = last cycle before arrest
        for j, i in enumerate(np.where(arr)[0]):
            v = ivk[i][~np.isnan(ivk[i])]
            if len(v): M[j, :min(len(v), NB)] = v[::-1][:NB]
        pos = np.arange(NB); nk = np.sum(~np.isnan(M), 0); valid = nk >= 10
        med = np.nanmedian(M, 0); q1 = np.nanpercentile(M, 25, 0); q3 = np.nanpercentile(M, 75, 0)
        axI.plot(-pos[valid], med[valid], '-o', color=col, lw=2.6, ms=6, label=lab)
        axI.fill_between(-pos[valid], q1[valid], q3[valid], color=col, alpha=0.14)
        for k in pos[valid]:
            axI.annotate(f'n={nk[k]}', (-k, med[k]), textcoords='offset points', xytext=(0, 9), fontsize=6.5, color=col, ha='center')
    axI.set_xlabel('completed cycle relative to arrest (0 = last cycle before G0)'); axI.set_ylabel('inter-division interval (h)')
    axI.set_title('(C) Cells SLOW approaching arrest (each aligned to its own G0)\narrested cells only; median ± IQR — the last cycle is the longest', fontweight='bold', fontsize=11)
    axI.set_xticks(-pos); axI.set_xticklabels([f'{-k}' if k else 'last' for k in pos]); axI.legend(fontsize=9); axI.grid(alpha=0.15)

    # (D) fraction transient-G0 & arrested vs time — WITH vs WITHOUT
    for cond, col, lab in [('W', '#8b1a1a', 'WITH mark'), ('N', '#4a4a4a', 'WITHOUT mark')]:
        St = z[f'{cond}_state']
        axF.plot(T, (St == 2).mean(0), color=col, lw=2.6, label=f'{lab} — arrested')
        axF.plot(T, (St == 1).mean(0), color=col, lw=1.5, ls='--', label=f'{lab} — transient G0')
    axF.axvline(D_WD / 60, color='k', ls=':', lw=1)
    axF.set_xlabel('time since withdrawal start (h)'); axF.set_ylabel('fraction of cells')
    axF.set_title('(D) G0 & arrest through withdrawal — WITH vs WITHOUT mark', fontweight='bold', fontsize=11)
    axF.legend(fontsize=8.5, loc='center right'); axF.grid(alpha=0.15); axF.set_ylim(0, 1)

    fig.suptitle(f'Single-cell ensemble (N={Nc}, heterogeneous temporal properties) through Hh withdrawal to SHH={SHH_LO:.2f} — transient G0 & period lengthening',
                 fontsize=13.5, fontweight='bold', y=1.0)
    plt.tight_layout()
    plt.savefig(f'{OUT}.png', dpi=150, bbox_inches='tight'); plt.savefig(f'{OUT}.pdf', bbox_inches='tight'); plt.close()
    print('WITH   median ndiv=%.1f  arrested@end=%.0f%%' % (np.median(z['W_ndiv']), 100 * (z['W_state'][:, -1] == 2).mean()))
    print('WITHOUT median ndiv=%.1f  arrested@end=%.0f%%' % (np.median(z['N_ndiv']), 100 * (z['N_state'][:, -1] == 2).mean()))
    print(f'Saved {OUT}.png')
