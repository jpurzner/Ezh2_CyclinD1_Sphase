"""Arrest model STAGE 1 (no core-model edit): EdU-analog S-gate (dDna/dt, not aRc), p27-based G0/G1
split, robust solver for strong Palbociclib blocks. Calibrate EdU + p27 thresholds on DMSO; score HU/Pal
vs measured MB55 proportions. RO3306 dropped."""
import sys, os
sys.path.insert(0, '.')
import numpy as np, tellurium as te, multiprocessing as mp
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.build_model_v44_heldt import build_model_v44

MB = dict(SHH=0.5, HHi=0.0, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002)
KPHRBCD0 = 0.35
SETTLE, TREAT, CHUNK, NP = 4800.0, 1440.0, 60.0, 13   # NP=13 -> 5-min sampling for dDna/dt
DT = CHUNK / (NP - 1)
KTL0, P21_MED, MU0 = 0.801, 0.42, 0.0005
CD_SDLOG, P27_SDLOG, MU_SDLOG, PART = 0.633, 0.32, 0.22, 0.30
TGT = {'DMSO': dict(G0=.2165, G1=.3788, S=.137, G2=.141),
       'HU':   dict(G0=.2571, G1=.4406, S=.1862, G2=.0328),
       'Pal':  dict(G0=.4239, G1=.253,  S=.0404, G2=.0956)}
CONDS = [('DMSO', {})]
CONDS += [(f'HU|d1|v{v}', {'HU': 1.0, 'vmin_fork': v}) for v in (0.06, 0.03, 0.012, 0.005)]
CONDS += [(f'Pal|{f}', {'kPhRbCd': round(KPHRBCD0*f, 4)}) for f in (0.5, 0.3, 0.15, 0.08)]
PHASES = ['G0', 'G1', 'S', 'G2']
SEL = ['time', 'Dna', 'P21']
_RR = None


def _init(_):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 800000)
    except Exception: pass


def robust_sim(tm, dur, npts):
    """Try progressively coarser tol, then sub-step, to survive stiff strong-block arrests."""
    for atol in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
        try:
            _RR.integrator.setValue('absolute_tolerance', atol)
            return _RR.simulate(tm, tm + dur, npts, selections=SEL)
        except Exception:
            continue
    # sub-step fallback: split into 4
    try:
        _RR.integrator.setValue('absolute_tolerance', 1e-6)
        sub = dur / 4.0; out = None
        for j in range(4):
            _RR.simulate(tm + j*sub, tm + (j+1)*sub, 2, selections=['time'])
        return _RR.simulate(tm, tm + dur, npts, selections=SEL)   # final short pass for output
    except Exception:
        return None


def _work(task):
    i, label, pars, ktl, p21d, mu = task
    rng = np.random.default_rng((hash((i, label)) & 0xFFFFFFFF))
    _RR.reset(); _RR['SHH'] = MB['SHH']
    for k, v in MB.items(): _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['mu'] = float(mu); _RR['P21_div'] = float(p21d)
    draw = lambda: float(p21d * np.exp(rng.normal(-PART * PART / 2.0, PART)))
    tm = 0.0; settle_end = SETTLE + int(rng.integers(0, 34)) * CHUNK
    while tm < settle_end:
        _RR['P21_div'] = draw()
        try: _RR.simulate(tm, tm + CHUNK, 2, selections=['time'])
        except Exception: return None
        tm += CHUNK
    for k, v in pars.items(): _RR[k] = float(v)
    T, DN, P2 = [], [], []
    t = 0.0
    while t < TREAT:
        _RR['P21_div'] = draw(); r = robust_sim(tm, CHUNK, NP)
        if r is None: break
        T.append(r['time'][1:] - tm + t); DN.append(r['Dna'][1:]); P2.append(r['P21'][1:])
        t += CHUNK; tm += CHUNK
    if not T: return None
    T = np.concatenate(T); DN = np.concatenate(DN); P2 = np.concatenate(P2)
    # synthesis rate (EdU analog): positive dDna/dt; mask division resets (negative jumps)
    dd = np.zeros_like(DN); dd[1:] = (DN[1:] - DN[:-1]) / DT
    dd[dd < 0] = 0.0
    return (label, DN, dd, P2)


def classify(DN, dd, P2, edu_thr, p27_thr):
    """EdU/DNA gating: S=EdU+; else 2N(Dna<0.5)->G0/G1 by p27; 4N(Dna>=0.5)->G2."""
    inS = dd > edu_thr
    twoN = (~inS) & (DN < 0.5); fourN = (~inS) & (DN >= 0.5)
    ph = np.full(len(DN), '', dtype='<U2')
    ph[inS] = 'S'; ph[fourN] = 'G2'
    ph[twoN & (P2 >= p27_thr)] = 'G0'; ph[twoN & (P2 < p27_thr)] = 'G1'
    return ph


if __name__ == '__main__':
    N = 40
    rng = np.random.default_rng(11)
    ktl = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)
    p21 = np.clip(P21_MED * np.exp(rng.normal(0, P27_SDLOG, N)), 0.15, 2.5)
    mu = np.clip(MU0 * np.exp(rng.normal(0, MU_SDLOG, N)), 0.00028, 0.00085)
    tasks = [(i, lab, pars, ktl[i], p21[i], mu[i]) for lab, pars in CONDS for i in range(N)]
    print(f'{len(CONDS)} conditions x N={N} = {len(tasks)} cells (robust solver, EdU+p27 classifier)', flush=True)
    with mp.get_context('spawn').Pool(9, initializer=_init, initargs=(None,)) as pool:
        res = [r for r in pool.imap_unordered(_work, tasks, chunksize=4) if r is not None]
    pool_by = {}
    for lab, _ in CONDS:
        sub = [r for r in res if r[0] == lab]
        if sub: pool_by[lab] = (np.concatenate([r[1] for r in sub]), np.concatenate([r[2] for r in sub]), np.concatenate([r[3] for r in sub]))

    # calibrate EdU threshold -> DMSO S=0.137, then p27 threshold -> DMSO G0=0.2165
    DN, dd, P2 = pool_by['DMSO']
    edu_thr = float(np.quantile(dd[dd > 0], 1 - TGT['DMSO']['S'] / np.mean(dd > 0))) if np.mean(dd > 0) > TGT['DMSO']['S'] else float(np.quantile(dd, 1 - TGT['DMSO']['S']))
    # simpler: pick edu_thr so fraction(dd>thr)=0.137
    edu_thr = float(np.quantile(dd, 1 - TGT['DMSO']['S']))
    inS = dd > edu_thr; twoN = (~inS) & (DN < 0.5)
    # p27 threshold so G0 fraction (of all) = 0.2165
    p2_2n = np.sort(P2[twoN])[::-1]   # high p27 first
    k = int(round(TGT['DMSO']['G0'] * len(DN)))
    p27_thr = float(p2_2n[min(k, len(p2_2n)-1)]) if len(p2_2n) else 0.5
    print(f'calibrated: edu_thr={edu_thr:.5f} (DMSO S->{TGT["DMSO"]["S"]}), p27_thr={p27_thr:.3f} (DMSO G0->{TGT["DMSO"]["G0"]})')

    def props(t):
        DN, dd, P2 = t; ph = classify(DN, dd, P2, edu_thr, p27_thr)
        return {p: np.mean(ph == p) for p in PHASES}

    def sse(p, base): return sum((p[ph] - TGT[base][ph])**2 for ph in PHASES)

    results = {lab: (props(pool_by[lab]) if lab in pool_by else None) for lab, _ in CONDS}
    print(f'\n{"condition":>14} | {"G0":>5} {"G1":>5} {"S":>5} {"G2":>5}  SSE   target')
    best = {'HU': (None, 9), 'Pal': (None, 9)}
    for lab, _ in CONDS:
        p = results[lab]
        if p is None: print(f'{lab:>14} | FAILED'); continue
        base = lab.split('|')[0]
        if base not in TGT: base = 'DMSO'
        s = sse(p, base); t = TGT[base]
        print(f'{lab:>14} | {p["G0"]:.2f}  {p["G1"]:.2f}  {p["S"]:.2f}  {p["G2"]:.2f}  {s:.3f}  [{t["G0"]:.2f}/{t["G1"]:.2f}/{t["S"]:.2f}/{t["G2"]:.2f}]')
        if base in ('HU', 'Pal') and s < best[base][1]: best[base] = (lab, s)
    print(f'\nbest HU  = {best["HU"][0]}  (SSE {best["HU"][1]:.3f})')
    print(f'best Pal = {best["Pal"][0]} (SSE {best["Pal"][1]:.3f})')

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.5))
    show = [('DMSO', 'DMSO'), (best['HU'][0], 'HU'), (best['Pal'][0], 'Pal')]
    x = np.arange(4); w = 0.38
    for a, (lab, base) in zip(ax, show):
        if lab is None or results.get(lab) is None: a.set_title(f'{base}: no fit'); continue
        p = results[lab]
        a.bar(x - w/2, [p[ph] for ph in PHASES], w, color='#8b1a1a', label='model (Stage 1)', alpha=0.85)
        a.bar(x + w/2, [TGT[base][ph] for ph in PHASES], w, color='#3b7bbf', label='data', alpha=0.85)
        a.set_xticks(x); a.set_xticklabels(PHASES); a.set_ylim(0, 0.55)
        a.set_title(f'{base}\n({lab})', fontweight='bold', fontsize=10); a.set_ylabel('proportion' if base == 'DMSO' else '')
        a.legend(fontsize=8); a.grid(alpha=0.15, axis='y')
    fig.suptitle(f'Arrest model STAGE 1 (EdU-gate + p27 split + robust solver, NO model edit) vs MB55 — edu_thr={edu_thr:.4f}, p27_thr={p27_thr:.2f}', fontweight='bold', fontsize=11.5, y=1.02)
    plt.tight_layout()
    plt.savefig('simulations/sim_ezh2_arrest_stage1.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_ezh2_arrest_stage1.pdf', bbox_inches='tight')
    print('Saved sim_ezh2_arrest_stage1.png')
