"""Arrest-model calibration: tune HU (dose + vmin_fork) and Palbociclib (kPhRbCd fraction) to the
measured MB55 phase proportions. G0/G1 threshold calibrated on DMSO (so DMSO G0=0.217). Runtime only.
RO3306 dropped. Reports best drug strengths + achievable match + residuals; figure model vs data."""
import sys, os
sys.path.insert(0, '.')
import numpy as np, tellurium as te, multiprocessing as mp
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.build_model_v44_heldt import build_model_v44

MB = dict(SHH=0.5, HHi=0.0, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002)
KPHRBCD0 = 0.35
SETTLE, TREAT, CHUNK, NP = 4800.0, 1440.0, 60.0, 9
KTL0, P21_MED, MU0 = 0.801, 0.42, 0.0005
CD_SDLOG, P27_SDLOG, MU_SDLOG, PART = 0.633, 0.32, 0.22, 0.30
TGT = {'DMSO': dict(G0=.2165, G1=.3788, S=.137, G2=.141),
       'HU':   dict(G0=.2571, G1=.4406, S=.1862, G2=.0328),
       'Pal':  dict(G0=.4239, G1=.253,  S=.0404, G2=.0956)}
# candidate interventions
CONDS = [('DMSO', {})]
CONDS += [(f'HU|d{d}|v{v}', {'HU': d, 'vmin_fork': v}) for d in (0.6, 1.0) for v in (0.1, 0.04, 0.015)]
CONDS += [(f'Pal|{f}', {'kPhRbCd': round(KPHRBCD0*f, 4)}) for f in (0.6, 0.45, 0.35, 0.25)]
PHASES = ['G0', 'G1', 'S', 'G2']
SEL = ['time', 'Dna', 'aRc', 'E2f']
_RR = None


def _init(_):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('absolute_tolerance', 1e-9); _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 300000)
    except Exception: pass


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
    DN, AR, E2 = [], [], []
    t = 0.0
    while t < TREAT:
        _RR['P21_div'] = draw(); r = None
        for atol in (1e-9, 1e-8, 1e-7):
            try:
                _RR.integrator.setValue('absolute_tolerance', atol)
                r = _RR.simulate(tm, tm + CHUNK, NP, selections=SEL); break
            except Exception: r = None
        if r is None: break
        DN.append(r['Dna'][1:]); AR.append(r['aRc'][1:]); E2.append(r['E2f'][1:])
        t += CHUNK; tm += CHUNK
    if not DN: return None
    return (label, np.concatenate(DN), np.concatenate(AR), np.concatenate(E2))


def props(DN, AR, E2, thr):
    inS = AR > 0.05; inG2 = (~inS) & (DN > 0.9); in2N = (~inS) & (DN <= 0.9)
    n = len(DN)
    return dict(S=inS.mean(), G2=inG2.mean(),
                G0=(in2N & (E2 < thr)).mean(), G1=(in2N & (E2 >= thr)).mean())


if __name__ == '__main__':
    N = 40
    rng = np.random.default_rng(11)
    ktl = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)
    p21 = np.clip(P21_MED * np.exp(rng.normal(0, P27_SDLOG, N)), 0.15, 2.5)
    mu = np.clip(MU0 * np.exp(rng.normal(0, MU_SDLOG, N)), 0.00028, 0.00085)
    tasks = [(i, lab, pars, ktl[i], p21[i], mu[i]) for lab, pars in CONDS for i in range(N)]
    print(f'{len(CONDS)} conditions x N={N} = {len(tasks)} cells', flush=True)
    with mp.get_context('spawn').Pool(9, initializer=_init, initargs=(None,)) as pool:
        res = [r for r in pool.imap_unordered(_work, tasks, chunksize=4) if r is not None]
    pool_by = {}
    for lab, _ in CONDS:
        sub = [r for r in res if r[0] == lab]
        if sub:
            pool_by[lab] = (np.concatenate([r[1] for r in sub]), np.concatenate([r[2] for r in sub]), np.concatenate([r[3] for r in sub]))
    # calibrate E2f threshold on DMSO so G0 fraction = target 0.2165
    DN, AR, E2 = pool_by['DMSO']; in2N = (AR <= 0.05) & (DN <= 0.9)
    e2_2n = np.sort(E2[in2N])
    # want (#2N with E2<thr)/total = 0.2165  -> thr = quantile of E2[2N] at (0.2165*total)/n2N
    ntot = len(DN); k = int(round(TGT['DMSO']['G0'] * ntot))
    thr = float(e2_2n[min(k, len(e2_2n)-1)]) if k < len(e2_2n) else float(e2_2n[-1])
    print(f'calibrated e2f_thr = {thr:.3f} (DMSO G0 -> {TGT["DMSO"]["G0"]:.3f})')

    def sse(p, base):
        return sum((p[ph] - TGT[base][ph])**2 for ph in PHASES)

    # pick best HU and best Pal
    results = {}
    for lab, _ in CONDS:
        if lab not in pool_by: results[lab] = None; continue
        results[lab] = props(*pool_by[lab], thr)
    print(f'\n{"condition":>14} | {"G0":>5} {"G1":>5} {"S":>5} {"G2":>5}  SSE   target(G0/G1/S/G2)')
    best = {'HU': (None, 9), 'Pal': (None, 9)}
    for lab, _ in CONDS:
        p = results[lab]
        if p is None: print(f'{lab:>14} | FAILED'); continue
        base = lab.split('|')[0].split('_')[0] if '|' in lab else lab
        if base not in TGT: base = 'DMSO'
        s = sse(p, base)
        t = TGT[base]
        print(f'{lab:>14} | {p["G0"]:.2f}  {p["G1"]:.2f}  {p["S"]:.2f}  {p["G2"]:.2f}  {s:.3f}  [{t["G0"]:.2f}/{t["G1"]:.2f}/{t["S"]:.2f}/{t["G2"]:.2f}]')
        if base in ('HU', 'Pal') and s < best[base][1]: best[base] = (lab, s)
    print(f'\nbest HU  = {best["HU"][0]}  (SSE {best["HU"][1]:.3f})')
    print(f'best Pal = {best["Pal"][0]} (SSE {best["Pal"][1]:.3f})')

    # figure: model vs data proportions for DMSO + best HU + best Pal
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.5))
    show = [('DMSO', 'DMSO'), (best['HU'][0], 'HU'), (best['Pal'][0], 'Pal')]
    x = np.arange(4); w = 0.38
    for a, (lab, base) in zip(ax, show):
        p = results[lab]
        a.bar(x - w/2, [p[ph] for ph in PHASES], w, color='#8b1a1a', label='model', alpha=0.85)
        a.bar(x + w/2, [TGT[base][ph] for ph in PHASES], w, color='#3b7bbf', label='data', alpha=0.85)
        a.set_xticks(x); a.set_xticklabels(PHASES); a.set_ylim(0, 0.55)
        a.set_title(f'{base}\n({lab})', fontweight='bold', fontsize=10); a.set_ylabel('proportion' if base == 'DMSO' else '')
        a.legend(fontsize=8); a.grid(alpha=0.15, axis='y')
    fig.suptitle(f'Arrest model (runtime-only) vs MB55 phase proportions — DMSO/HU/Pal (e2f_thr={thr:.2f} calibrated on DMSO-G0)', fontweight='bold', fontsize=12, y=1.02)
    plt.tight_layout()
    plt.savefig('simulations/sim_ezh2_arrest_calib.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_ezh2_arrest_calib.pdf', bbox_inches='tight')
    print('Saved sim_ezh2_arrest_calib.png')
