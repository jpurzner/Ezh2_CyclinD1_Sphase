"""Sweep EZH2 protein degradation rate kDeEZ -> impact on the Hh-withdrawal ramp, WITH the
within-cycle EZH2 phase ratio (measured G0:1 / G1:1.16 / S:1.31 / G2:1.50) as the honesty check.

Protocol per cell (GNP, mark ON, Hh 0.75): settle cycling -> cycling window (phase-bin) -> ramp
Hh 0.75->0.05 over 48h -> floor 150h. Record EZH2/Dna/aRc/E2f. Sweep kDeEZ across half-lives.
Panel A: EZH2(t) normalized to pre-withdrawal cycling level (decay kinetics) per kDeEZ.
Panel B: within-cycle EZH2 phase ratio (/G0) per kDeEZ vs the measured IF data.
"""
import sys, os
sys.path.insert(0, '.')
import numpy as np, tellurium as te, multiprocessing as mp
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.build_model_v44_heldt import build_model_v44

GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, EZH2i=0)
HH_HI, HH_LO, F0 = 0.75, 0.05, 0.233
SETTLE, T_CYC, T_DOWN, FLOOR, CHUNK, NP = 6000.0, 4800.0, 2880.0, 9000.0, 60.0, 9
KTL0, KTLEZ0, P21_MED, MU0 = 0.801, 0.004, 0.42, 0.0005
CD_SDLOG, P27_SDLOG, MU_SDLOG, PART = 0.633, 0.32, 0.22, 0.30
KDEEZ = [0.00015, 0.0005, 0.001, 0.002, 0.004, 0.008]          # per-min; half-life = ln2/k/60 h
SEL = ['time', 'EZH2', 'Dna', 'aRc', 'E2f']
T_TRACE = np.arange(0.0, T_CYC + T_DOWN + FLOOR + 1e-6, 30.0)   # 0.5h grid, post-settle
DATA_RATIO = {'G0': 1.0, 'G1': 1.159, 'S': 1.31, 'G2': 1.50}
_RR = None


def _init(_):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('absolute_tolerance', 1e-9); _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 300000)
    except Exception: pass


def shh_of(t):
    if t < T_CYC: return HH_HI
    if t < T_CYC + T_DOWN: return HH_HI + (HH_LO - HH_HI) * (t - T_CYC) / T_DOWN
    return HH_LO


def _work(task):
    i, kdeez, ktl, p21d, mu = task
    rng = np.random.default_rng((hash((i, round(kdeez, 6))) & 0xFFFFFFFF))
    _RR.reset(); _RR['SHH'] = HH_HI
    for k, v in GNP.items(): _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['kTlEZ'] = KTLEZ0; _RR['mu'] = float(mu)
    _RR['f0_mk'] = F0; _RR['P21_div'] = float(p21d); _RR['kDeEZ'] = float(kdeez)
    draw = lambda: float(p21d * np.exp(rng.normal(-PART * PART / 2.0, PART)))
    tm = 0.0; settle_end = SETTLE + int(rng.integers(0, 34)) * CHUNK
    while tm < settle_end:
        _RR['P21_div'] = draw()
        try: _RR.simulate(tm, tm + CHUNK, 2, selections=['time'])
        except Exception: return None
        tm += CHUNK
    T, EZ, DN, AR, E2 = [], [], [], [], []
    t = 0.0; TEND = T_CYC + T_DOWN + FLOOR
    while t < TEND:
        _RR['P21_div'] = draw(); _RR['SHH'] = float(shh_of(t)); r = None
        for atol in (1e-9, 1e-8, 1e-7):
            try:
                _RR.integrator.setValue('absolute_tolerance', atol)
                r = _RR.simulate(tm, tm + CHUNK, NP, selections=SEL); break
            except Exception: r = None
        if r is None: break
        T.append(r['time'][1:] - tm + t); EZ.append(r['EZH2'][1:]); DN.append(r['Dna'][1:]); AR.append(r['aRc'][1:]); E2.append(r['E2f'][1:])
        t += CHUNK; tm += CHUNK
    if not T: return None
    T = np.concatenate(T); EZ = np.concatenate(EZ); DN = np.concatenate(DN); AR = np.concatenate(AR); E2 = np.concatenate(E2)
    ez_t = np.interp(T_TRACE, T, EZ).astype(np.float32)               # withdrawal timecourse
    cyc = T < T_CYC                                                   # cycling-window samples for phase bins
    return (kdeez, ez_t, EZ[cyc].astype(np.float32), DN[cyc].astype(np.float32), AR[cyc].astype(np.float32), E2[cyc].astype(np.float32))


def classify(dna, arc, e2f, thr):
    ph = np.full(len(dna), '', dtype='<U2')
    inS = arc > 0.05; inG2 = (~inS) & (dna > 0.9); in2N = (~inS) & (dna <= 0.9)
    ph[inS] = 'S'; ph[inG2] = 'G2'; ph[in2N & (e2f >= thr)] = 'G1'; ph[in2N & (e2f < thr)] = 'G0'
    return ph


if __name__ == '__main__':
    N = 60
    rng = np.random.default_rng(11)
    ktl = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)
    p21 = np.clip(P21_MED * np.exp(rng.normal(0, P27_SDLOG, N)), 0.15, 2.5)
    mu = np.clip(MU0 * np.exp(rng.normal(0, MU_SDLOG, N)), 0.00028, 0.00085)
    tasks = [(i, k, ktl[i], p21[i], mu[i]) for k in KDEEZ for i in range(N)]
    print(f'sweep: {len(tasks)} cells, {len(KDEEZ)} kDeEZ values', flush=True)
    with mp.get_context('spawn').Pool(9, initializer=_init, initargs=(None,)) as pool:
        res = [r for r in pool.imap_unordered(_work, tasks, chunksize=4) if r is not None]

    fig, ax = plt.subplots(1, 2, figsize=(15, 5.8))
    cmap = plt.cm.viridis(np.linspace(0, 0.9, len(KDEEZ)))
    trace_h = T_TRACE / 60.0
    ramp0, ramp1 = T_CYC / 60.0, (T_CYC + T_DOWN) / 60.0
    print(f'{"kDeEZ":>8} {"t1/2(h)":>8} {"G1/G0":>6} {"S/G0":>6} {"G2/G0":>6}  cyc-EZH2')
    ratios = {}
    for j, k in enumerate(KDEEZ):
        sub = [r for r in res if abs(r[0] - k) < 1e-9]
        ez = np.nanmean(np.array([r[1] for r in sub]), 0)
        # phase ratios from pooled cycling samples
        EZc = np.concatenate([r[2] for r in sub]); DNc = np.concatenate([r[3] for r in sub])
        ARc = np.concatenate([r[4] for r in sub]); E2c = np.concatenate([r[5] for r in sub])
        thr = float(np.median(E2c[(ARc <= 0.05) & (DNc <= 0.9)]))
        ph = classify(DNc, ARc, E2c, thr)
        med = {p: (np.median(EZc[ph == p]) if np.any(ph == p) else np.nan) for p in ['G0', 'G1', 'S', 'G2']}
        rr = {p: med[p] / med['G0'] for p in ['G0', 'G1', 'S', 'G2']}
        ratios[k] = rr
        th = np.log(2) / k / 60.0
        # pre-withdrawal cycling level = mean EZH2 over last 40h of cycling
        base = ez[(trace_h > ramp0 - 40) & (trace_h < ramp0)].mean()
        ax[0].plot(trace_h, ez / base, color=cmap[j], lw=2.3, label=f'kDeEZ {k:g}  (t½ {th:.0f} h)')
        ax[1].plot(range(4), [rr[p] for p in ['G0', 'G1', 'S', 'G2']], '-o', color=cmap[j], lw=2.0, ms=5, label=f't½ {th:.0f} h')
        print(f'{k:8g} {th:8.1f} {rr["G1"]:6.2f} {rr["S"]:6.2f} {rr["G2"]:6.2f}  {base:.3f}')

    a = ax[0]
    a.axvspan(ramp0, ramp1, color='#c0392b', alpha=0.08); a.axvspan(ramp1, trace_h[-1], color='grey', alpha=0.06)
    a.axhline(0.5, color='k', ls=':', lw=0.8, alpha=0.5)
    a.set_xlabel('time (h) — withdrawal ramp starts at 80 h'); a.set_ylabel('EZH2 ÷ pre-withdrawal cycling level')
    a.set_title('(A) Impact of kDeEZ on the Hh-withdrawal EZH2 decay\n(normalized; red = ramp 0.75→0.05, grey = floor)', fontweight='bold', fontsize=11)
    a.legend(fontsize=8, loc='upper right'); a.grid(alpha=0.15); a.set_ylim(0, 1.05)

    a = ax[1]
    a.plot(range(4), [DATA_RATIO[p] for p in ['G0', 'G1', 'S', 'G2']], 'k--s', lw=2.6, ms=8, label='MEASURED (IF)', zorder=10)
    a.set_xticks(range(4)); a.set_xticklabels(['G0', 'G1', 'S', 'G2'])
    a.set_xlabel('cell-cycle phase'); a.set_ylabel('EZH2 ÷ G0 (within-cycle)')
    a.set_title('(B) Within-cycle EZH2 phase ratio vs kDeEZ\nthe measured G2/G0=1.50× is the constraint (which kDeEZ still fits?)', fontweight='bold', fontsize=11)
    a.legend(fontsize=8, loc='upper left'); a.grid(alpha=0.15)

    fig.suptitle(f'EZH2 protein degradation sweep (kDeEZ) — withdrawal-ramp decay vs the measured within-cycle phase ratio  (GNP, Hh {HH_HI}, mark ON, N={N}/value)',
                 fontsize=12, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig('simulations/sim_ezh2_kdeez_sweep.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_ezh2_kdeez_sweep.pdf', bbox_inches='tight')
    print('Saved sim_ezh2_kdeez_sweep.png')
