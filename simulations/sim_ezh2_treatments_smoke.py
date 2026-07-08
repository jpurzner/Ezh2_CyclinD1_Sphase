"""EZH2 protein by phase across 4 arrest drugs (MB55, 24h) vs model — fine (kEZbas x kDeEZ) fit.
Drugs: DMSO / HU(S) / Palbociclib(G1) / RO3306(G2). Protocol: settle asynchronous cycling (no drug,
random phase) -> apply drug -> 24h -> pool samples over the window, phase-bin EZH2. Normalize to
DMSO-G0 (the data anchor). Model file untouched (runtime params).

--smoke : tiny grid/N to sanity-check phase shifts + DMSO baseline.
"""
import sys, os
sys.path.insert(0, '.')
import numpy as np, tellurium as te, multiprocessing as mp
from src.build_model_v44_heldt import build_model_v44

SMOKE = '--smoke' in sys.argv
MB = dict(SHH=0.5, HHi=0.0, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002)
DRUGS = {'DMSO': {}, 'HU': {'HU': 1.0}, 'Pal': {'kPhRbCd': 0.0}, 'RO': {'k25': 0.0}}
SETTLE, TREAT, CHUNK, NP = 4800.0, 1440.0, 60.0, 9
KTL0, KTLEZ0, P21_MED, MU0 = 0.801, 0.004, 0.42, 0.0005
CD_SDLOG, P27_SDLOG, MU_SDLOG, PART = 0.633, 0.32, 0.22, 0.30
# targets: fold vs DMSO-G0, (value, n_cells)
TGT = {
 'DMSO': {'G0': (1.00, 3376), 'G1': (1.159, 5855), 'S': (1.31, 2099), 'G2': (1.50, 2178)},
 'HU':   {'G0': (1.052, 2427), 'G1': (1.501, 4189), 'S': (1.709, 1738), 'G2': (1.56, 286)},
 'Pal':  {'G0': (0.964, 3417), 'G1': (0.967, 1996), 'S': (1.17, 317), 'G2': (1.148, 725)},
 'RO':   {'G0': (0.973, 1906), 'G1': (1.264, 2161), 'S': (1.543, 1629), 'G2': (1.737, 1422)},
}
PHASES = ['G0', 'G1', 'S', 'G2']
SEL = ['time', 'EZH2', 'EZH2m', 'Dna', 'aRc', 'E2f']
_RR = None


def _init(_):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('absolute_tolerance', 1e-9); _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 300000)
    except Exception: pass


def _work(task):
    i, kb, kd, drug, ktl, p21d, mu = task
    rng = np.random.default_rng((hash((i, round(kb, 6), round(kd, 6), drug)) & 0xFFFFFFFF))
    _RR.reset(); _RR['SHH'] = MB['SHH']
    for k, v in MB.items(): _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['kTlEZ'] = KTLEZ0; _RR['mu'] = float(mu)
    _RR['P21_div'] = float(p21d); _RR['kEZbas'] = float(kb); _RR['kDeEZ'] = float(kd)
    draw = lambda: float(p21d * np.exp(rng.normal(-PART * PART / 2.0, PART)))
    tm = 0.0; settle_end = SETTLE + int(rng.integers(0, 34)) * CHUNK
    while tm < settle_end:                                   # settle asynchronous, NO drug
        _RR['P21_div'] = draw()
        try: _RR.simulate(tm, tm + CHUNK, 2, selections=['time'])
        except Exception: return None
        tm += CHUNK
    for k, v in DRUGS[drug].items(): _RR[k] = float(v)       # apply drug
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
    return (kb, kd, drug, np.concatenate(EZ), np.concatenate(EM), np.concatenate(DN), np.concatenate(AR), np.concatenate(E2))


def classify(dna, arc, e2f, thr):
    ph = np.full(len(dna), '', dtype='<U2')
    inS = arc > 0.05; inG2 = (~inS) & (dna > 0.9); in2N = (~inS) & (dna <= 0.9)
    ph[inS] = 'S'; ph[inG2] = 'G2'; ph[in2N & (e2f >= thr)] = 'G1'; ph[in2N & (e2f < thr)] = 'G0'
    return ph


def phase_medians(res, kb, kd):
    """Return {drug: {phase: median EZH2}}, {drug:{phase:frac}}, using DMSO-2N E2f threshold."""
    byd = {}
    for drug in DRUGS:
        sub = [r for r in res if abs(r[0]-kb) < 1e-9 and abs(r[1]-kd) < 1e-9 and r[2] == drug]
        if not sub: return None, None, None
        byd[drug] = dict(EZ=np.concatenate([r[3] for r in sub]), DN=np.concatenate([r[5] for r in sub]),
                         AR=np.concatenate([r[6] for r in sub]), E2=np.concatenate([r[7] for r in sub]))
    dm = byd['DMSO']; thr = float(np.median(dm['E2'][(dm['AR'] <= 0.05) & (dm['DN'] <= 0.9)]))
    med = {}; frac = {}
    for drug, d in byd.items():
        ph = classify(d['DN'], d['AR'], d['E2'], thr)
        med[drug] = {p: (np.median(d['EZ'][ph == p]) if np.any(ph == p) else np.nan) for p in PHASES}
        frac[drug] = {p: float(np.mean(ph == p)) for p in PHASES}
    return med, frac, thr


if __name__ == '__main__':
    N = 6 if SMOKE else 26
    KEZBAS = [0.004] if SMOKE else list(np.round(np.geomspace(0.0008, 0.012, 12), 5))
    KDEEZ = [0.0005] if SMOKE else list(np.round(np.geomspace(0.00015, 0.004, 10), 6))
    rng = np.random.default_rng(11)
    ktl = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)
    p21 = np.clip(P21_MED * np.exp(rng.normal(0, P27_SDLOG, N)), 0.15, 2.5)
    mu = np.clip(MU0 * np.exp(rng.normal(0, MU_SDLOG, N)), 0.00028, 0.00085)
    tasks = [(i, kb, kd, drug, ktl[i], p21[i], mu[i]) for kb in KEZBAS for kd in KDEEZ for drug in DRUGS for i in range(N)]
    print(f'{"SMOKE " if SMOKE else ""}{len(KEZBAS)}x{len(KDEEZ)} grid x4 drugs, {len(tasks)} cells', flush=True)
    with mp.get_context('spawn').Pool(9, initializer=_init, initargs=(None,)) as pool:
        res = [r for r in pool.imap_unordered(_work, tasks, chunksize=4) if r is not None]

    if SMOKE:
        med, frac, thr = phase_medians(res, KEZBAS[0], KDEEZ[0])
        g0 = med['DMSO']['G0']
        print(f'\nDMSO-G0 EZH2 (denominator) = {g0:.3f}   e2f_thr={thr:.3f}')
        print(f'{"drug":>5} | ' + ' '.join(f'{p:>16}' for p in PHASES))
        print('       | ' + ' '.join(f'{"fold(frac)":>16}' for p in PHASES) + '   [target fold]')
        for drug in DRUGS:
            row = ' '.join(f'{med[drug][p]/g0:6.2f}({frac[drug][p]:4.2f})    ' for p in PHASES)
            tg = ' '.join(f'{TGT[drug][p][0]:.2f}' for p in PHASES)
            print(f'{drug:>5} | {row}  [{tg}]')
        sys.exit(0)
