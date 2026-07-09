"""Full EZH2-stability re-fit: kEZbas x kEZE2f x kDeEZ x kTlEZ vs ALL 9 EZH2-coupled validation targets
(6 EZH2 + 3 CyclinD1), scored with the exact validate_v44 tolerances. Keeps the HU S-entry edit
(KmHU_fire=0.25). Parallel sample-and-refine. Goal: fast kDeEZ (HU boost + withdrawal) WITHOUT breaking
EZH2i-fold / G0-cycling / transcript / MB-CyclinD1."""
import sys, os
sys.path.insert(0, '.')
import numpy as np, tellurium as te, multiprocessing as mp
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

P27 = 0.1
GNP = dict(SHH=0.5, Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0)
MB = dict(SHH=0.5, Ptch1_copy_number=0.1, MYCN_amplification=2.8, p16=0.306, p18=1.553, kSyP21=0.002)
CONDS = {'GNP': (GNP, {}), 'GNP_HHi': (GNP, {'HHi': 1.0}), 'GNP_EZH2i': (GNP, {'EZH2i': 1.0}),
         'GNP_SS': (GNP, {'k_Cd_translation': 0.0}), 'MB': (MB, {}), 'MB_HHi': (MB, {'HHi': 1.0}), 'MB_HU': (MB, {'HU': 1.0})}
SEL = ['time', 'Cd_mRNA', 'EZH2', 'EZH2m', 'MPF', 'P21', 'aRc', 'Dna', 'vfork']
FREE = ['kEZbas', 'kEZE2f', 'kDeEZ', 'kTlEZ', 'kDeEZm']
BOX = {'kEZbas': (0.0003, 0.004), 'kEZE2f': (0.01, 0.06), 'kDeEZ': (0.0002, 0.0012), 'kTlEZ': (0.002, 0.02), 'kDeEZm': (0.002, 0.02)}
# target: (value, tol); some depend on cross-condition ratios computed below
TGT = {'Cd_GNPHHi': (0.157, .30), 'Cd_MBHHi': (0.144, .30), 'Cd_MBGNP': (5.07, .35),
       'EZ_MBGNP': (2.05, .25), 'EZi_fold': (2.2, .30), 'EZ_G0cyc': (0.6, .30),
       'EZ_protG2G0': (1.48, .35), 'EZ_transSG0': (2.0, .30), 'HU_boost': (1.31, .25)}
_RR = None
_DEF = {}


def _init(_):
    global _RR, _DEF
    _RR = te.loada(build_model_v44()); _RR.integrator.setValue('relative_tolerance', 1e-6)
    # capture builder defaults for every perturbable param (reset() does NOT restore params -> avoid leak)
    _DEF = {k: float(_RR[k]) for k in ('SHH', 'Ptch1_copy_number', 'MYCN_amplification', 'p16', 'p18',
                                       'kSyP21', 'HHi', 'EZH2i', 'HU', 'k_Cd_translation', 'kPhRbCd')}


def _sim(ctx, drug, pars):
    for ms in (1e9, 20.0, 5.0):
        for atol in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
            _RR.reset()
            for k, v in _DEF.items(): _RR[k] = v          # full baseline every call (no param leak)
            for k, v in ctx.items(): _RR[k] = v
            for k, v in drug.items(): _RR[k] = v
            for k, v in pars.items(): _RR[k] = v
            try:
                _RR.integrator.setValue('maximum_num_steps', 1000000); _RR.integrator.setValue('maximum_time_step', ms)
                _RR.integrator.setValue('absolute_tolerance', atol)
                return _RR.simulate(0, 12000, 48000, selections=SEL)
            except Exception: pass
    return None


def settled(res, sp, s=4000):
    t = res['time']; return float(np.mean(res[sp][t >= s]))


def phase_ez(res, s=4000):
    t = res['time']; m = t >= s; tt = t[m]; dt = tt[1] - tt[0]
    P = res['P21'][m]; aRc = res['aRc'][m]; D = res['Dna'][m]; vf = res['vfork'][m]; ez = res['EZH2'][m]; em = res['EZH2m'][m]
    inS = (vf * aRc > 0.02) & (D < 0.98); inG2 = D >= 0.98; pre = ~inS & ~inG2
    G0 = pre & (P > P27); G2 = inG2
    def rat(v, num, den): return float(v[num].mean()) / float(v[den].mean()) if num.any() and den.any() and v[den].mean() else np.nan
    protG2G0 = rat(ez, G2, G0)
    boostS = float(ez[inS].mean()) if inS.any() else np.nan
    # transcript S/G0 uses aRc-gate (DMSO gradient, matches validate_v44.classify_grad)
    inS_a = (aRc > 0.05) & (D < 0.98); G0a = (~inS_a & ~inG2) & (P > P27)
    transSG0 = float(em[inS_a].mean()) / float(em[G0a].mean()) if inS_a.any() and G0a.any() and em[G0a].mean() else np.nan
    return protG2G0, boostS, transSG0


def _eval(pars):
    r = {}
    for nm, (ctx, drug) in CONDS.items():
        r[nm] = _sim(ctx, drug, pars)
        if r[nm] is None: return None
    cd = lambda c: settled(r[c], 'Cd_mRNA'); ez = lambda c: settled(r[c], 'EZH2')
    pg2, boost0, trans = phase_ez(r['MB']); _, boost1, _ = phase_ez(r['MB_HU'])
    out = {'Cd_GNPHHi': cd('GNP_HHi') / cd('GNP'), 'Cd_MBHHi': cd('MB_HHi') / cd('MB'), 'Cd_MBGNP': cd('MB') / cd('GNP'),
           'EZ_MBGNP': ez('MB') / ez('GNP'), 'EZi_fold': cd('GNP_EZH2i') / cd('GNP'),
           'EZ_G0cyc': ez('GNP_SS') / ez('GNP'), 'EZ_protG2G0': pg2, 'EZ_transSG0': trans,
           'HU_boost': boost1 / boost0 if boost0 else np.nan}
    return out


def score(o):
    if o is None: return (-1, -1e9)
    np_ = sum(1 for k, (t, tol) in TGT.items() if o[k] == o[k] and abs(o[k] - t) <= tol * t)
    re = sum(abs((o[k] if o[k] == o[k] else 9 * t) - t) / t for k, (t, tol) in TGT.items())
    return (np_, -re)


def _work(task):
    idx, pars = task
    return (idx, _eval(pars), pars)


def sample(box, n, rng, center=None, frac=1.0):
    S = []
    for _ in range(n):
        p = {}
        for k in FREE:
            lo, hi = box[k]
            if center is None: p[k] = float(np.exp(rng.uniform(np.log(lo), np.log(hi))))
            else: p[k] = float(np.clip(np.exp(rng.uniform(np.log(center[k]) - frac*np.log(10), np.log(center[k]) + frac*np.log(10))), lo, hi))
        S.append(p)
    return S


if __name__ == '__main__':
    rng = np.random.default_rng(3)
    S1 = sample(BOX, 120, rng)
    allp, allo = [], []
    print(f'stage 1: {len(S1)} param sets x 7 conditions ...', flush=True)
    with mp.get_context('spawn').Pool(9, initializer=_init, initargs=(None,)) as pool:
        for idx, o, pars in pool.imap_unordered(_work, list(enumerate(S1)), chunksize=2):
            allp.append(pars); allo.append(o)
    sc = [score(o) for o in allo]; bi = int(np.argmax([s[0]*100 + s[1] for s in sc]))
    print(f'  stage1 best {sc[bi][0]}/9 at {[f"{k}={allp[bi][k]:.4g}" for k in FREE]}', flush=True)
    order = sorted(range(len(allo)), key=lambda i: sc[i][0]*100 + sc[i][1], reverse=True)[:5]
    cen = {k: float(np.exp(np.mean([np.log(allp[j][k]) for j in order]))) for k in FREE}
    S2 = sample(BOX, 90, rng, center=cen, frac=0.35)
    print(f'stage 2: refine around {[f"{k}={cen[k]:.4g}" for k in FREE]} ...', flush=True)
    with mp.get_context('spawn').Pool(9, initializer=_init, initargs=(None,)) as pool:
        for idx, o, pars in pool.imap_unordered(_work, list(enumerate(S2)), chunksize=2):
            allp.append(pars); allo.append(o)
    sc = [score(o) for o in allo]; bi = int(np.argmax([s[0]*100 + s[1] for s in sc]))
    best, bo = allp[bi], allo[bi]
    print(f'\n=== BEST {sc[bi][0]}/9 EZH2-coupled targets ===')
    for k in FREE:
        th = f'  (t1/2 {np.log(2)/best[k]/60:.0f}h)' if k == 'kDeEZ' else ''
        print(f'  {k:8s} = {best[k]:.5g}{th}')
    print(f'{"target":>14} {"model":>8} {"want":>8}  ok')
    for k, (t, tol) in TGT.items():
        v = bo[k]; okc = (v == v and abs(v - t) <= tol * t)
        print(f'{k:>14} {v:8.3f} {t:8.3f}  {"OK" if okc else "--"}')
