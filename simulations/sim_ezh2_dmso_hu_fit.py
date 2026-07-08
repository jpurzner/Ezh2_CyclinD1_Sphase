"""Finer EZH2 rate fit vs DMSO + HU + scRNA-seq (MB55, MB context). Pal/RO excluded (model cell cycle
collapses to 2N under those arrests). Fine (kEZbas x kDeEZ) grid + sensitivity of the OTHER rate vars
(kEZE2f, Kez_cd, kDeEZm, wCe). Model file untouched (runtime overrides).

Targets (fold vs DMSO-G0):
  protein DMSO  G0 1.00 / G1 1.159 / S 1.31 / G2 1.50   (IF)
  protein HU    G1 1.501 / S 1.709                       (IF; longer S -> more EZH2)
  mRNA DMSO     S 1.87 / G2 2.17                          (scRNA-seq)
"""
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
# EZH2 rate defaults (for sensitivity)
DEF = dict(kEZbas=0.00027, kEZE2f=0.022, Kez_cd=2.94, kDeEZm=0.02, wCe=0.582, kDeEZ=0.00015)
KEZBAS = list(np.round(np.geomspace(0.0008, 0.012, 14), 5))
KDEEZ = list(np.round(np.geomspace(0.00015, 0.004, 11), 6))
# fit targets: (condition, obs, phase, value, weight)
TARGETS = [('DMSO', 'p', 'G1', 1.159, 1), ('DMSO', 'p', 'S', 1.31, 1), ('DMSO', 'p', 'G2', 1.50, 1),
           ('HU', 'p', 'G1', 1.501, 1), ('HU', 'p', 'S', 1.709, 1),
           ('DMSO', 'm', 'S', 1.87, 0.7), ('DMSO', 'm', 'G2', 2.17, 0.7)]
PHASES = ['G0', 'G1', 'S', 'G2']
SEL = ['time', 'EZH2', 'EZH2m', 'Dna', 'aRc', 'E2f']
CACHE = 'simulations/sim_ezh2_dmso_hu_fit_cache.npz'
_RR = None


def _init(_):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('absolute_tolerance', 1e-9); _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 300000)
    except Exception: pass


def _work(task):
    i, pars, drug, ktl, p21d, mu = task
    key = tuple(sorted(pars.items()))
    rng = np.random.default_rng((hash((i, key, drug)) & 0xFFFFFFFF))
    _RR.reset(); _RR['SHH'] = MB['SHH']
    for k, v in MB.items(): _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['kTlEZ'] = KTLEZ0; _RR['mu'] = float(mu); _RR['P21_div'] = float(p21d)
    for k, v in pars.items(): _RR[k] = float(v)
    draw = lambda: float(p21d * np.exp(rng.normal(-PART * PART / 2.0, PART)))
    tm = 0.0; settle_end = SETTLE + int(rng.integers(0, 34)) * CHUNK
    while tm < settle_end:
        _RR['P21_div'] = draw()
        try: _RR.simulate(tm, tm + CHUNK, 2, selections=['time'])
        except Exception: return None
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
    if not EZ: return None
    return (key, drug, np.concatenate(EZ), np.concatenate(EM), np.concatenate(DN), np.concatenate(AR), np.concatenate(E2))


def classify(dna, arc, e2f, thr):
    ph = np.full(len(dna), '', dtype='<U2')
    inS = arc > 0.05; inG2 = (~inS) & (dna > 0.9); in2N = (~inS) & (dna <= 0.9)
    ph[inS] = 'S'; ph[inG2] = 'G2'; ph[in2N & (e2f >= thr)] = 'G1'; ph[in2N & (e2f < thr)] = 'G0'
    return ph


def predict(res, key):
    """Return protein & mRNA folds (÷ DMSO-G0) by (condition, phase)."""
    byd = {}
    for drug in DRUGS:
        sub = [r for r in res if r[0] == key and r[1] == drug]
        if not sub: return None
        byd[drug] = dict(EZ=np.concatenate([r[2] for r in sub]), EM=np.concatenate([r[3] for r in sub]),
                         DN=np.concatenate([r[4] for r in sub]), AR=np.concatenate([r[5] for r in sub]), E2=np.concatenate([r[6] for r in sub]))
    dm = byd['DMSO']; thr = float(np.median(dm['E2'][(dm['AR'] <= 0.05) & (dm['DN'] <= 0.9)]))
    phd = classify(dm['DN'], dm['AR'], dm['E2'], thr)
    g0p = np.median(dm['EZ'][phd == 'G0']); g0m = np.median(dm['EM'][phd == 'G0'])
    out = {}
    for drug, d in byd.items():
        ph = classify(d['DN'], d['AR'], d['E2'], thr)
        for p in PHASES:
            msk = ph == p
            out[(drug, 'p', p)] = (np.median(d['EZ'][msk]) / g0p) if np.any(msk) else np.nan
            out[(drug, 'm', p)] = (np.median(d['EM'][msk]) / g0m) if np.any(msk) else np.nan
    return out


def sse(pred):
    s = 0.0; w = 0.0
    for cond, obs, ph, val, wt in TARGETS:
        m = pred.get((cond, obs, ph), np.nan)
        if np.isfinite(m): s += wt * (m - val)**2; w += wt
    return s / w if w else np.nan


def _arg(f, d):
    for a in sys.argv:
        if a.startswith(f + '='): return a.split('=', 1)[1]
    return d


if __name__ == '__main__':
    N = int(_arg('--n', 24))
    if '--fresh' in sys.argv or not os.path.exists(CACHE):
        rng = np.random.default_rng(11)
        ktl = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)
        p21 = np.clip(P21_MED * np.exp(rng.normal(0, P27_SDLOG, N)), 0.15, 2.5)
        mu = np.clip(MU0 * np.exp(rng.normal(0, MU_SDLOG, N)), 0.00028, 0.00085)
        # grid over (kEZbas, kDeEZ) at default other params
        combos = [dict(DEF, kEZbas=kb, kDeEZ=kd) for kb in KEZBAS for kd in KDEEZ]
        # sensitivity: vary each secondary param around a provisional best (kEZbas=0.004, kDeEZ=0.0005)
        PROV = dict(DEF, kEZbas=0.004, kDeEZ=0.0005)
        sens = {'kEZE2f': [0.011, 0.022, 0.033], 'Kez_cd': [1.5, 2.94, 5.0], 'kDeEZm': [0.01, 0.02, 0.04], 'wCe': [0.35, 0.582, 0.8]}
        for pname, vals in sens.items():
            for v in vals: combos.append(dict(PROV, **{pname: v}))
        seen = set(); ucombos = []
        for c in combos:
            k = tuple(sorted(c.items()))
            if k not in seen: seen.add(k); ucombos.append(c)
        tasks = [(i, c, drug, ktl[i], p21[i], mu[i]) for c in ucombos for drug in DRUGS for i in range(N)]
        print(f'{len(ucombos)} param sets x2 drugs, {len(tasks)} cells (N={N})', flush=True)
        with mp.get_context('spawn').Pool(9, initializer=_init, initargs=(None,)) as pool:
            res = [r for r in pool.imap_unordered(_work, tasks, chunksize=4) if r is not None]
        # store predictions per param set
        store = {}
        for c in ucombos:
            key = tuple(sorted(c.items())); pred = predict(res, key)
            if pred is None: continue
            store[key] = (pred, sse(pred))
        import pickle
        with open(CACHE.replace('.npz', '.pkl'), 'wb') as f:
            pickle.dump(dict(store=store, KEZBAS=KEZBAS, KDEEZ=KDEEZ, DEF=DEF, PROV=PROV, sens=sens), f)
        print('cached', flush=True)

    import pickle
    with open(CACHE.replace('.npz', '.pkl'), 'rb') as f:
        D = pickle.load(f)
    store = D['store']; KEZBAS = D['KEZBAS']; KDEEZ = D['KDEEZ']; DEF = D['DEF']; PROV = D['PROV']; sens = D['sens']

    def getkey(**kw):
        return tuple(sorted(dict(DEF, **kw).items()))
    # landscape
    L = np.full((len(KEZBAS), len(KDEEZ)), np.nan)
    for a, kb in enumerate(KEZBAS):
        for b, kd in enumerate(KDEEZ):
            k = getkey(kEZbas=kb, kDeEZ=kd)
            if k in store: L[a, b] = store[k][1]
    ba = np.unravel_index(np.nanargmin(L), L.shape)
    bkb, bkd = KEZBAS[ba[0]], KDEEZ[ba[1]]
    bpred = store[getkey(kEZbas=bkb, kDeEZ=bkd)][0]
    print(f'BEST (kEZbas x kDeEZ): kEZbas={bkb:g}  kDeEZ={bkd:g} (t½ {np.log(2)/bkd/60:.0f}h)  SSE={L[ba]:.4f}')
    for cond, obs, ph, val, wt in TARGETS:
        print(f'  {cond:4s} {obs} {ph}: model {bpred[(cond,obs,ph)]:.2f}  target {val:.2f}')

    fig = plt.figure(figsize=(16, 5.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.15, 1])
    # (A) landscape
    a = fig.add_subplot(gs[0, 0])
    im = a.imshow(np.log10(L), aspect='auto', origin='lower', cmap='viridis_r')
    a.set_xticks(range(len(KDEEZ))); a.set_xticklabels([f'{k:g}\n{np.log(2)/k/60:.0f}h' for k in KDEEZ], fontsize=7, rotation=0)
    a.set_yticks(range(len(KEZBAS))); a.set_yticklabels([f'{k:g}' for k in KEZBAS], fontsize=7)
    a.set_xlabel('kDeEZ (protein t½)'); a.set_ylabel('kEZbas (baseline transcription)')
    a.plot(ba[1], ba[0], '*', ms=22, color='#ff3b3b', mec='white', mew=1)
    a.axhline(KEZBAS.index(min(KEZBAS, key=lambda x: abs(x-0.00027))) if 0.00027 in KEZBAS else 0, color='w', ls=':', lw=0.8, alpha=0.4)
    fig.colorbar(im, ax=a, label='log10 fit SSE (DMSO+HU+mRNA)')
    a.set_title(f'(A) Fine fit landscape — ★ best\nkEZbas={bkb:g}, kDeEZ={bkd:g} (t½ {np.log(2)/bkd/60:.0f}h)', fontweight='bold', fontsize=10.5)
    # (B) best fit vs data
    a = fig.add_subplot(gs[0, 1])
    xx = np.arange(4)
    a.plot(xx, [1, 1.159, 1.31, 1.50], 'k--s', lw=2.2, ms=8, label='DMSO protein — data', zorder=9)
    a.plot(xx, [bpred[('DMSO','p',p)] for p in PHASES], '-o', color='#8b1a1a', lw=2.2, ms=6, label='DMSO protein — model')
    a.plot(xx, [1.052, 1.501, 1.709, 1.56], 'k--^', lw=1.8, ms=8, alpha=0.55, label='HU protein — data', zorder=9)
    a.plot(xx, [bpred[('HU','p',p)] for p in PHASES], '-^', color='#d1495b', lw=1.8, ms=6, label='HU protein — model')
    a.plot(xx, [1, 1.87, 1.87, 2.17], 'k--D', lw=1.6, ms=7, alpha=0.4, label='mRNA — scRNAseq', zorder=9)
    a.plot(xx, [bpred[('DMSO','m',p)] for p in PHASES], '-D', color='#3b7bbf', lw=1.6, ms=5, label='mRNA — model')
    a.set_xticks(xx); a.set_xticklabels(PHASES); a.set_xlabel('cell-cycle phase'); a.set_ylabel('÷ DMSO-G0')
    a.set_title('(B) Best fit vs data — DMSO & HU protein + mRNA\n(HU longer-S boost captured; Pal/RO excluded — model 2N-collapse)', fontweight='bold', fontsize=10.5)
    a.legend(fontsize=7.5, ncol=1); a.grid(alpha=0.15)
    # (C) secondary-param sensitivity (ΔSSE around provisional best)
    a = fig.add_subplot(gs[0, 2])
    base_sse = store[tuple(sorted(PROV.items()))][1]
    rows = []
    for pname, vals in sens.items():
        ss = [store[tuple(sorted(dict(PROV, **{pname: v}).items()))][1] for v in vals]
        rows.append((pname, vals, ss))
    ylab = []
    for j, (pname, vals, ss) in enumerate(rows):
        a.plot([min(ss), max(ss)], [j, j], '-', color='#999', lw=1, zorder=1)
        a.scatter(ss, [j]*len(ss), c=range(len(ss)), cmap='coolwarm', s=70, zorder=3, edgecolor='k', lw=0.5)
        ylab.append(f'{pname}\n{vals}')
    a.axvline(base_sse, color='#ff3b3b', ls='--', lw=1.5, label=f'provisional best SSE={base_sse:.3f}')
    a.set_yticks(range(len(rows))); a.set_yticklabels(ylab, fontsize=8)
    a.set_xlabel('fit SSE (lower = better)')
    a.set_title('(C) Secondary rate-var sensitivity\n(around kEZbas 0.004 / kDeEZ 0.0005; blue→red = low→high value)', fontweight='bold', fontsize=10.5)
    a.legend(fontsize=8); a.grid(alpha=0.15, axis='x')

    fig.suptitle('EZH2 rate fit vs DMSO + HU + scRNA-seq (MB55, MB context, N=%d) — Pal/RO pending arrest-analog fix' % N, fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig('simulations/sim_ezh2_dmso_hu_fit.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_ezh2_dmso_hu_fit.pdf', bbox_inches='tight')
    print('Saved sim_ezh2_dmso_hu_fit.png')
