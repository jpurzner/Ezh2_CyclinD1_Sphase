"""Fit EZH2 protein-by-phase (baseline transcription kEZbas x protein stability kDeEZ) and make the
figure: (A) fit landscape, (B) best-fit protein & mRNA by phase vs data, (C) single-cell EZH2 across
the cycle (S-rise/G2-peak/division-drop on ongoing G0 baseline). Model file untouched (runtime params)."""
import sys, os
sys.path.insert(0, '.')
import numpy as np, tellurium as te, multiprocessing as mp
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.build_model_v44_heldt import build_model_v44

GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, EZH2i=0)
HH, F0 = 0.75, 0.233
SETTLE, RECORD, CHUNK, NP = 6000.0, 6000.0, 60.0, 9
KTL0, KTLEZ0, P21_MED, MU0 = 0.801, 0.004, 0.42, 0.0005
CD_SDLOG, P27_SDLOG, MU_SDLOG, PART = 0.633, 0.32, 0.22, 0.30
KEZBAS = [0.00027, 0.001, 0.002, 0.003, 0.004, 0.007, 0.011]
KDEEZ = [0.00015, 0.0005, 0.001, 0.002, 0.004]
PT = {'G1': 1.159, 'S': 1.31, 'G2': 1.50}                 # protein targets (IF)
RT = {'S': 1.87, 'G2': 2.17}                              # mRNA targets (scRNA-seq)
BEST_KB, BEST_KD = 0.004, 0.0005                          # joint pick for panels B/C (protein+RNA+some G0 drop)
SEL = ['time', 'EZH2', 'EZH2m', 'Dna', 'aRc', 'E2f']
CACHE = 'simulations/sim_ezh2_phase_fit_cache.npz'
_RR = None


def _init(_):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('absolute_tolerance', 1e-9); _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 300000)
    except Exception: pass


def _setup(kb, kd, ktl, p21d, mu):
    _RR.reset(); _RR['SHH'] = HH
    for k, v in GNP.items(): _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['kTlEZ'] = KTLEZ0; _RR['mu'] = float(mu)
    _RR['f0_mk'] = F0; _RR['P21_div'] = float(p21d); _RR['kEZbas'] = float(kb); _RR['kDeEZ'] = float(kd)


def _work(task):
    i, kb, kd, ktl, p21d, mu = task
    rng = np.random.default_rng((hash((i, round(kb, 6), round(kd, 6))) & 0xFFFFFFFF))
    _setup(kb, kd, ktl, p21d, mu)
    draw = lambda: float(p21d * np.exp(rng.normal(-PART * PART / 2.0, PART)))
    tm = 0.0; settle_end = SETTLE + int(rng.integers(0, 34)) * CHUNK
    while tm < settle_end:
        _RR['P21_div'] = draw()
        try: _RR.simulate(tm, tm + CHUNK, 2, selections=['time'])
        except Exception: return None
        tm += CHUNK
    EZ, EM, DN, AR, E2 = [], [], [], [], []
    t = 0.0
    while t < RECORD:
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
    return (kb, kd, np.concatenate(EZ), np.concatenate(EM), np.concatenate(DN), np.concatenate(AR), np.concatenate(E2))


def classify(dna, arc, e2f, thr):
    ph = np.full(len(dna), '', dtype='<U2')
    inS = arc > 0.05; inG2 = (~inS) & (dna > 0.9); in2N = (~inS) & (dna <= 0.9)
    ph[inS] = 'S'; ph[inG2] = 'G2'; ph[in2N & (e2f >= thr)] = 'G1'; ph[in2N & (e2f < thr)] = 'G0'
    return ph


def ratios(x, ph):
    g0 = np.median(x[ph == 'G0'])
    return {p: (np.median(x[ph == p]) / g0 if np.any(ph == p) else np.nan) for p in ['G0', 'G1', 'S', 'G2']}


def trajectory(kb, kd, ncyc_min=17000):
    """Single deterministic-ish cell: EZH2 protein/mRNA/Dna over several cycles at best params."""
    _setup(kb, kd, KTL0, P21_MED, MU0)
    tm = 0.0
    while tm < SETTLE:
        try: _RR.simulate(tm, tm + CHUNK, 2, selections=['time'])
        except Exception: break
        tm += CHUNK
    for atol in (1e-9, 1e-8, 1e-7):
        try:
            _RR.integrator.setValue('absolute_tolerance', atol)
            d = _RR.simulate(tm, tm + ncyc_min, 3000, selections=['time', 'EZH2', 'EZH2m', 'Dna', 'aRc'])
            return (d['time'] - tm) / 60.0, d['EZH2'], d['EZH2m'], d['Dna'], d['aRc']
        except Exception: continue
    return None


if __name__ == '__main__':
    FRESH = '--fresh' in sys.argv
    if FRESH or not os.path.exists(CACHE):
        N = 40
        rng = np.random.default_rng(11)
        ktl = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)
        p21 = np.clip(P21_MED * np.exp(rng.normal(0, P27_SDLOG, N)), 0.15, 2.5)
        mu = np.clip(MU0 * np.exp(rng.normal(0, MU_SDLOG, N)), 0.00028, 0.00085)
        tasks = [(i, kb, kd, ktl[i], p21[i], mu[i]) for kb in KEZBAS for kd in KDEEZ for i in range(N)]
        print(f'grid {len(KEZBAS)}x{len(KDEEZ)}, {len(tasks)} cells', flush=True)
        with mp.get_context('spawn').Pool(9, initializer=_init, initargs=(None,)) as pool:
            res = [r for r in pool.imap_unordered(_work, tasks, chunksize=4) if r is not None]
        pSSE = np.full((len(KEZBAS), len(KDEEZ)), np.nan); rSSE = np.full_like(pSSE, np.nan)
        pr_grid = {}; mr_grid = {}
        for a, kb in enumerate(KEZBAS):
            for b, kd in enumerate(KDEEZ):
                sub = [r for r in res if abs(r[0]-kb) < 1e-9 and abs(r[1]-kd) < 1e-9]
                if not sub: continue
                EZ = np.concatenate([r[2] for r in sub]); EM = np.concatenate([r[3] for r in sub])
                DN = np.concatenate([r[4] for r in sub]); AR = np.concatenate([r[5] for r in sub]); E2 = np.concatenate([r[6] for r in sub])
                thr = float(np.median(E2[(AR <= 0.05) & (DN <= 0.9)]))
                ph = classify(DN, AR, E2, thr); pr = ratios(EZ, ph); mr = ratios(EM, ph)
                pSSE[a, b] = sum((pr[p]-PT[p])**2 for p in PT); rSSE[a, b] = sum((mr[p]-RT[p])**2 for p in RT)
                pr_grid[(a, b)] = [pr[p] for p in ['G0', 'G1', 'S', 'G2']]; mr_grid[(a, b)] = [mr[p] for p in ['G0', 'G1', 'S', 'G2']]
        # best-param single-cell trajectories (two kDeEZ at best baseline)
        with mp.get_context('spawn').Pool(1, initializer=_init, initargs=(None,)) as pool:
            trS = pool.apply(trajectory, (BEST_KB, 0.00015))     # slow (77h)
            trF = pool.apply(trajectory, (BEST_KB, 0.0005))      # 23h
        ai = KEZBAS.index(BEST_KB)
        np.savez(CACHE, KEZBAS=KEZBAS, KDEEZ=KDEEZ, pSSE=pSSE, rSSE=rSSE,
                 pr_best=np.array(pr_grid[(ai, KDEEZ.index(0.0005))]), mr_best=np.array(mr_grid[(ai, KDEEZ.index(0.0005))]),
                 trS=np.array(trS, dtype=object), trF=np.array(trF, dtype=object))
        print('cached ->', CACHE, flush=True)

    z = np.load(CACHE, allow_pickle=True)
    KB = list(z['KEZBAS']); KD = list(z['KDEEZ']); pSSE = z['pSSE']; rSSE = z['rSSE']
    thh = [np.log(2)/k/60 for k in KD]

    fig = plt.figure(figsize=(16, 5.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 1.3])

    # (A) protein-fit SSE heatmap with RNA~2x band
    a = fig.add_subplot(gs[0, 0])
    im = a.imshow(np.log10(pSSE), aspect='auto', origin='lower', cmap='viridis_r')
    a.set_xticks(range(len(KD))); a.set_xticklabels([f'{k:g}\n{t:.0f}h' for k, t in zip(KD, thh)], fontsize=8)
    a.set_yticks(range(len(KB))); a.set_yticklabels([f'{k:g}' for k in KB], fontsize=8)
    a.set_xlabel('kDeEZ (protein t½)'); a.set_ylabel('kEZbas (baseline transcription)')
    # mark RNA-consistent cells (rSSE small) with a dot; protein-best with a star
    for aa in range(len(KB)):
        for bb in range(len(KD)):
            if np.isfinite(rSSE[aa, bb]) and rSSE[aa, bb] < 0.15:
                a.plot(bb, aa, 'o', ms=5, mfc='none', mec='white', mew=1.5)
    ba = np.unravel_index(np.nanargmin(pSSE + 2*rSSE), pSSE.shape)   # joint best
    a.plot(ba[1], ba[0], '*', ms=20, color='#ff3b3b', mec='white', mew=1.0)
    a.axhline(KB.index(0.00027), color='w', ls=':', lw=1, alpha=0.6)
    a.text(len(KD)-0.5, KB.index(0.00027), ' current', color='w', va='center', ha='right', fontsize=8)
    fig.colorbar(im, ax=a, label='log10 protein-fit SSE')
    a.set_title('(A) Fit landscape — protein SSE\n○ = mRNA also ~2× (scRNA-seq);  ★ = joint best', fontweight='bold', fontsize=10.5)

    # (B) best-fit protein & mRNA by phase vs data
    a = fig.add_subplot(gs[0, 1])
    phs = ['G0', 'G1', 'S', 'G2']; xx = np.arange(4)
    a.plot(xx, [1, PT['G1'], PT['S'], PT['G2']], 'k--s', lw=2.4, ms=9, label='protein — MEASURED (IF)', zorder=10)
    a.plot(xx, z['pr_best'], '-o', color='#8b1a1a', lw=2.4, ms=7, label='protein — model (fit)')
    a.plot(xx, [1, np.nan, RT['S'], RT['G2']], 'k--^', lw=1.8, ms=8, alpha=0.6, label='mRNA — scRNA-seq')
    a.plot(xx, z['mr_best'], '-^', color='#3b7bbf', lw=1.8, ms=6, alpha=0.85, label='mRNA — model (fit)')
    a.set_xticks(xx); a.set_xticklabels(phs); a.set_xlabel('cell-cycle phase'); a.set_ylabel('÷ G0')
    a.set_title(f'(B) Best fit: kEZbas={BEST_KB:g}, kDeEZ=0.0005 (t½ 23h)\nprotein 1.5× & mRNA ~2× both matched', fontweight='bold', fontsize=10.5)
    a.legend(fontsize=8); a.grid(alpha=0.15)

    # (C) single-cell EZH2 across the cycle
    a = fig.add_subplot(gs[0, 2])
    for tr, lab, col in [(z['trS'], 'kDeEZ 0.00015 (t½ 77h) — over-stable', '#c0392b'),
                         (z['trF'], 'kDeEZ 0.0005 (t½ 23h) — fit', '#8b1a1a')]:
        th, ez, em, dna, arc = tr
        th = np.asarray(th, float); ez = np.asarray(ez, float)
        m = th < 90
        a.plot(th[m], ez[m], '-', color=col, lw=2.2, label=lab)
    th, ez, em, dna, arc = z['trF']; th = np.asarray(th, float)
    m = th < 90
    ab = a.twinx()
    ab.fill_between(th[m], 0, np.asarray(arc, float)[m] > 0.05, color='#2e8b57', alpha=0.10, step='mid')
    ab.plot(th[m], np.asarray(em, float)[m], ':', color='#3b7bbf', lw=1.5, label='EZH2 mRNA')
    ab.set_ylabel('EZH2 mRNA (dotted) / S-phase (green)', color='grey', fontsize=9)
    a.set_xlabel('time (h) — successive cycles'); a.set_ylabel('EZH2 protein')
    a.set_title('(C) EZH2 across the cell cycle (fitted baseline)\nrise in S → peak G2 → drop at division, on ongoing G0 baseline', fontweight='bold', fontsize=10.5)
    a.legend(fontsize=8, loc='lower right')

    fig.suptitle('EZH2 protein-by-phase fit — raise baseline transcription (kEZbas), then kDeEZ set by the G0/withdrawal drop  (GNP, Hh 0.75, mark ON)',
                 fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig('simulations/sim_ezh2_phase_fit.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_ezh2_phase_fit.pdf', bbox_inches='tight')
    print('Saved sim_ezh2_phase_fit.png')
    print('joint-best (★): kEZbas=%g  kDeEZ=%g (t½ %.0fh)' % (KB[ba[0]], KD[ba[1]], np.log(2)/KD[ba[1]]/60))
