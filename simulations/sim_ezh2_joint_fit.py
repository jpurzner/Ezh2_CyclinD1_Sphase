"""Joint 4-param EZH2 rate optimization: kEZbas x kEZE2f x kDeEZ x kDeEZm vs DMSO+HU protein + mRNA
(MB55, MB context). Parallel sample-and-refine (broad log-uniform sweep -> tighten around best), which
also gives per-parameter identifiability. Kez_cd/wCe fixed (minor). Model file untouched (runtime)."""
import sys, os
sys.path.insert(0, '.')
import numpy as np, tellurium as te, multiprocessing as mp
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.build_model_v44_heldt import build_model_v44

MB = dict(SHH=0.5, HHi=0.0, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002)
DRUGS = {'DMSO': {}, 'HU': {'HU': 1.0}}
SETTLE, TREAT, CHUNK, NP = 4800.0, 1440.0, 60.0, 9
KTL0, KTLEZ0, P21_MED, MU0 = 0.801, 0.004, 0.42, 0.0005
CD_SDLOG, P27_SDLOG, MU_SDLOG, PART = 0.633, 0.32, 0.22, 0.30
FIXED = dict(Kez_cd=2.94, wCe=0.582, kTlEZ=KTLEZ0)
FREE = ['kEZbas', 'kEZE2f', 'kDeEZ', 'kDeEZm']
BOX = {'kEZbas': (0.0005, 0.020), 'kEZE2f': (0.008, 0.12), 'kDeEZ': (0.00008, 0.003), 'kDeEZm': (0.0015, 0.03)}
DEFAULTS = {'kEZbas': 0.00027, 'kEZE2f': 0.022, 'kDeEZ': 0.00015, 'kDeEZm': 0.02}
TARGETS = [('DMSO', 'p', 'G1', 1.159, 1.0), ('DMSO', 'p', 'S', 1.31, 1.0), ('DMSO', 'p', 'G2', 1.50, 1.0),
           ('HU', 'p', 'G0', 1.052, 0.5), ('HU', 'p', 'G1', 1.501, 1.0), ('HU', 'p', 'S', 1.709, 1.0),
           ('DMSO', 'm', 'S', 1.87, 0.7), ('DMSO', 'm', 'G2', 2.17, 0.7)]
PHASES = ['G0', 'G1', 'S', 'G2']
SEL = ['time', 'EZH2', 'EZH2m', 'Dna', 'aRc', 'E2f']
_RR = None; _CELLS = None


def _init(cells):
    global _RR, _CELLS
    _CELLS = cells
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('absolute_tolerance', 1e-9); _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 300000)
    except Exception: pass


def _work(task):
    pid, pars, drug, i = task
    ktl, p21d, mu = _CELLS[i]
    rng = np.random.default_rng((hash((pid, drug, i)) & 0xFFFFFFFF))
    _RR.reset(); _RR['SHH'] = MB['SHH']
    for k, v in MB.items(): _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['mu'] = float(mu); _RR['P21_div'] = float(p21d)
    for k, v in FIXED.items(): _RR[k] = float(v)
    for k, v in pars.items(): _RR[k] = float(v)
    draw = lambda: float(p21d * np.exp(rng.normal(-PART * PART / 2.0, PART)))
    tm = 0.0; settle_end = SETTLE + int(rng.integers(0, 34)) * CHUNK
    while tm < settle_end:
        _RR['P21_div'] = draw()
        try: _RR.simulate(tm, tm + CHUNK, 2, selections=['time'])
        except Exception: return (pid, drug, None)
        tm += CHUNK
    for k, v in DRUGS[drug].items(): _RR[k] = float(v)
    EZ, EM, DN, AR, E2 = [], [], [], [], []
    t = 0.0
    while t < TREAT:
        _RR['P21_div'] = draw(); r = None
        for atol in (1e-9, 1e-8, 1e-7):
            try:
                _RR.integrator.setValue('absolute_tolerance', atol)
                r = _RR.simulate(tm, tm + CHUNK, NP, selections=SEL); break
            except Exception: r = None
        if r is None: break
        EZ.append(r['EZH2'][1:]); EM.append(r['EZH2m'][1:]); DN.append(r['Dna'][1:]); AR.append(r['aRc'][1:]); E2.append(r['E2f'][1:])
        t += CHUNK; tm += CHUNK
    if not EZ: return (pid, drug, None)
    return (pid, drug, (np.concatenate(EZ), np.concatenate(EM), np.concatenate(DN), np.concatenate(AR), np.concatenate(E2)))


def classify(dna, arc, e2f, thr):
    ph = np.full(len(dna), '', dtype='<U2')
    inS = arc > 0.05; inG2 = (~inS) & (dna > 0.9); in2N = (~inS) & (dna <= 0.9)
    ph[inS] = 'S'; ph[inG2] = 'G2'; ph[in2N & (e2f >= thr)] = 'G1'; ph[in2N & (e2f < thr)] = 'G0'
    return ph


def predict_from(byd):
    dm = byd.get('DMSO'); hu = byd.get('HU')
    if dm is None or hu is None: return None
    EZ, EM, DN, AR, E2 = dm
    thr = float(np.median(E2[(AR <= 0.05) & (DN <= 0.9)]))
    phd = classify(DN, AR, E2, thr)
    if not np.any(phd == 'G0'): return None
    g0p = np.median(EZ[phd == 'G0']); g0m = np.median(EM[phd == 'G0'])
    out = {}
    for drug, d in byd.items():
        e, m, dn, ar, e2 = d; ph = classify(dn, ar, e2, thr)
        for p in PHASES:
            msk = ph == p
            out[(drug, 'p', p)] = (np.median(e[msk]) / g0p) if np.any(msk) else np.nan
            out[(drug, 'm', p)] = (np.median(m[msk]) / g0m) if np.any(msk) else np.nan
    return out


def sse(pred):
    if pred is None: return np.inf
    s = w = 0.0
    for cond, obs, ph, val, wt in TARGETS:
        m = pred.get((cond, obs, ph), np.nan)
        if not np.isfinite(m): return np.inf
        s += wt * (m - val)**2; w += wt
    return s / w


def evaluate(param_sets, cells, workers=9):
    """Run all param_sets x 2 drugs x N cells in one pool; return list of (pred, sse)."""
    N = len(cells)
    tasks = [(pid, ps, drug, i) for pid, ps in enumerate(param_sets) for drug in DRUGS for i in range(N)]
    agg = {pid: {} for pid in range(len(param_sets))}
    raw = {pid: {d: [] for d in DRUGS} for pid in range(len(param_sets))}
    with mp.get_context('spawn').Pool(workers, initializer=_init, initargs=(cells,)) as pool:
        for pid, drug, dat in pool.imap_unordered(_work, tasks, chunksize=6):
            if dat is not None: raw[pid][drug].append(dat)
    out = []
    for pid, ps in enumerate(param_sets):
        byd = {}
        for drug in DRUGS:
            lst = raw[pid][drug]
            if not lst: byd = None; break
            byd[drug] = tuple(np.concatenate([x[k] for x in lst]) for k in range(5))
        pred = predict_from(byd) if byd else None
        out.append((pred, sse(pred)))
    return out


def sample_box(box, n, rng, center=None, frac=1.0):
    """log-uniform samples; if center given, sample within +-frac dex around center (clipped to box)."""
    S = []
    for _ in range(n):
        ps = {}
        for k in FREE:
            lo, hi = box[k]
            if center is None:
                ps[k] = float(np.exp(rng.uniform(np.log(lo), np.log(hi))))
            else:
                c = np.log(center[k]); ps[k] = float(np.clip(np.exp(rng.uniform(c - frac*np.log(10), c + frac*np.log(10))), lo, hi))
        S.append(ps)
    return S


if __name__ == '__main__':
    N = 22
    rng = np.random.default_rng(7)
    ktl = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)
    p21 = np.clip(P21_MED * np.exp(rng.normal(0, P27_SDLOG, N)), 0.15, 2.5)
    mu = np.clip(MU0 * np.exp(rng.normal(0, MU_SDLOG, N)), 0.00028, 0.00085)
    cells = list(zip(ktl, p21, mu))

    all_ps, all_sse, all_pred = [], [], []
    # stage 1: broad
    S1 = sample_box(BOX, 240, rng)
    print(f'stage 1: {len(S1)} param sets ...', flush=True)
    r1 = evaluate(S1, cells)
    for ps, (pred, s) in zip(S1, r1):
        all_ps.append(ps); all_sse.append(s); all_pred.append(pred)
    best_i = int(np.argmin(all_sse)); best = all_ps[best_i]
    print(f'  stage1 best SSE={all_sse[best_i]:.4f} at {[f"{k}={best[k]:g}" for k in FREE]}', flush=True)
    # stage 2: refine around top-5 centroids
    order = np.argsort(all_sse)[:5]
    cen = {k: float(np.exp(np.mean([np.log(all_ps[j][k]) for j in order]))) for k in FREE}
    S2 = sample_box(BOX, 140, rng, center=cen, frac=0.35)
    print(f'stage 2: {len(S2)} refine around {[f"{k}={cen[k]:g}" for k in FREE]} ...', flush=True)
    r2 = evaluate(S2, cells)
    for ps, (pred, s) in zip(S2, r2):
        all_ps.append(ps); all_sse.append(s); all_pred.append(pred)

    all_sse = np.array(all_sse)
    bi = int(np.argmin(all_sse)); best = all_ps[bi]; bpred = all_pred[bi]
    print(f'\n=== JOINT BEST  SSE={all_sse[bi]:.4f} ===')
    for k in FREE:
        th = f'  (t½ {np.log(2)/best[k]/60:.0f}h)' if k in ('kDeEZ', 'kDeEZm') else ''
        print(f'  {k:8s} = {best[k]:.5g}   (default {DEFAULTS[k]:g}){th}')
    print('  fit vs data:')
    for cond, obs, ph, val, wt in TARGETS:
        print(f'    {cond:4s} {obs} {ph}: model {bpred[(cond,obs,ph)]:.2f}  target {val:.2f}')

    import pickle
    with open('simulations/sim_ezh2_joint_fit_cache.pkl', 'wb') as f:
        pickle.dump(dict(all_ps=all_ps, all_sse=all_sse, best=best, bpred=bpred, FREE=FREE, BOX=BOX, DEFAULTS=DEFAULTS, TARGETS=TARGETS), f)

    # ---- figure ----
    fig = plt.figure(figsize=(16, 8.5))
    gs = fig.add_gridspec(2, 4)
    finite = np.isfinite(all_sse)
    for j, k in enumerate(FREE):
        a = fig.add_subplot(gs[0, j])
        vals = np.array([p[k] for p in all_ps])[finite]; ss = all_sse[finite]
        sc = a.scatter(vals, ss, c=np.log10(ss), cmap='viridis_r', s=14, alpha=0.6)
        a.axvline(best[k], color='#ff3b3b', lw=2, label='best')
        a.axvline(DEFAULTS[k], color='k', ls=':', lw=1.5, label='default')
        a.set_xscale('log'); a.set_yscale('log'); a.set_xlabel(k); a.set_ylabel('fit SSE' if j == 0 else '')
        th = f'\nt½ {np.log(2)/best[k]/60:.0f}h' if k in ('kDeEZ', 'kDeEZm') else ''
        a.set_title(f'{k} = {best[k]:.4g}{th}', fontweight='bold', fontsize=10)
        a.legend(fontsize=7); a.grid(alpha=0.15, which='both')
    # best-fit vs data
    a = fig.add_subplot(gs[1, :2]); xx = np.arange(4)
    a.plot(xx, [1, 1.159, 1.31, 1.50], 'k--s', lw=2.2, ms=8, label='DMSO protein — data', zorder=9)
    a.plot(xx, [bpred[('DMSO', 'p', p)] for p in PHASES], '-o', color='#8b1a1a', lw=2.2, ms=6, label='DMSO protein — model')
    a.plot(xx, [1.052, 1.501, 1.709, 1.56], 'k--^', lw=1.8, ms=8, alpha=0.55, label='HU protein — data', zorder=9)
    a.plot(xx, [bpred[('HU', 'p', p)] for p in PHASES], '-^', color='#d1495b', lw=1.8, ms=6, label='HU protein — model')
    a.plot(xx, [1, 1.87, 1.87, 2.17], 'k--D', lw=1.6, ms=7, alpha=0.4, label='mRNA — scRNAseq', zorder=9)
    a.plot(xx, [bpred[('DMSO', 'm', p)] for p in PHASES], '-D', color='#3b7bbf', lw=1.6, ms=5, label='mRNA — model')
    a.set_xticks(xx); a.set_xticklabels(PHASES); a.set_xlabel('cell-cycle phase'); a.set_ylabel('÷ DMSO-G0')
    a.set_title(f'Joint 4-param best fit  (SSE {all_sse[bi]:.4f} vs 2-param 0.064)', fontweight='bold', fontsize=11)
    a.legend(fontsize=8, ncol=2); a.grid(alpha=0.15)
    # residual table as text
    a = fig.add_subplot(gs[1, 2:]); a.axis('off')
    lines = ['JOINT BEST parameters (÷ default):', '']
    for k in FREE:
        th = f'   t½ {np.log(2)/best[k]/60:.0f}h' if k in ('kDeEZ', 'kDeEZm') else ''
        lines.append(f'  {k:8s} = {best[k]:.4g}   ({best[k]/DEFAULTS[k]:.1f}× default){th}')
    lines += ['', 'residuals (model − target):']
    for cond, obs, ph, val, wt in TARGETS:
        m = bpred[(cond, obs, ph)]
        lines.append(f'  {cond:4s} {obs} {ph:2s}: {m:.2f}  (tgt {val:.2f}, Δ {m-val:+.2f})')
    a.text(0.02, 0.98, '\n'.join(lines), va='top', ha='left', fontsize=10, family='monospace', transform=a.transAxes)

    fig.suptitle('EZH2 joint 4-parameter fit (kEZbas × kEZE2f × kDeEZ × kDeEZm) vs DMSO+HU+mRNA (MB55) — Pal/RO pending arrest fix',
                 fontsize=13, fontweight='bold', y=1.0)
    plt.tight_layout()
    plt.savefig('simulations/sim_ezh2_joint_fit.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_ezh2_joint_fit.pdf', bbox_inches='tight')
    print('Saved sim_ezh2_joint_fit.png')
