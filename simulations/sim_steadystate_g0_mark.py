"""STEADY-STATE transient-G0 fraction vs constant Hh — REPRESSED (f0_mk=0.233) vs NOT repressed (f0_mk=1.0).

Companion to sim_singlecell_withdrawal_g0: during WITHDRAWAL the transient-G0 count is confounded by arrest
(a repressed cell that would dwell in transient G0 may instead drop to terminal arrest). Holding Hh CONSTANT
removes that — each heterogeneous cell settles into a limit cycle (or arrests) and we classify it once.

Per-cell heterogeneity (paired across conditions; same draws as the withdrawal sim):
  mu, del_mk, kDeEZ, P21_div (birth p27), k_Cd_translation (CyclinD1 setpoint).
Classification at steady state (probe_transient_g0 / sim_g0_bifurcation conventions):
  divisions = MPF peaks; dwell = period × fraction-of-cycle in (pre-S ∧ p27>0.1).
  <2 divisions               -> ARREST (deep quiescence)
  cycling, dwell ≥ 2 h       -> TRANSIENT G0 (CDK2-low branch)
  cycling, dwell < 2 h       -> IMMEDIATE re-entry (CDK2-inc branch)

Mechanistic question: the mark lowers CyclinD1 (lower CyclinD1/p27 ratio → p27-dominant). Does that give
MORE transient G0 at steady state (opposite of the arrest-confounded withdrawal appearance)?

Run:  ./venv/bin/python simulations/sim_steadystate_g0_mark.py [--n=400] [--workers=9] [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import tellurium as te
import multiprocessing as mp
from src.build_model_v44_heldt import build_model_v44

SEL = ['time', 'MPF', 'P21', 'aRc', 'Dna', 'Cd']
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, EZH2i=0)
SHH_LEVELS = [0.30, 0.40, 0.50, 0.65, 0.80, 1.00]
T_END, N_PTS, SETTLE, DWELL_CUT, P27_HI = 12000, 20000, 6000, 2.0, 0.1

_RR = None


def _init_worker(_):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('absolute_tolerance', 1e-9); _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 300000)
    except Exception: pass


def _work(task):
    i, shh, f0, mu, delmk, kdeez, p21d, ktl = task
    _RR.reset()
    _RR['SHH'] = shh
    for k, v in GNP.items(): _RR[k] = v
    _RR['mu'] = mu; _RR['del_mk'] = delmk; _RR['kDeEZ'] = kdeez
    _RR['P21_div'] = p21d; _RR['k_Cd_translation'] = ktl; _RR['f0_mk'] = f0
    try:
        r = _RR.simulate(0, T_END, N_PTS, selections=SEL)
    except Exception:
        return (i, shh, f0, 'arrest', np.nan, np.nan, np.nan)
    t = r['time']; m = t >= SETTLE; tt = t[m]; dt = tt[1] - tt[0]
    cd = float(np.mean(r['Cd'][m]))
    pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=max(1, int(200 / dt)))
    if len(pk) < 2:
        return (i, shh, f0, 'arrest', np.nan, cd, np.nan)
    per = float(np.mean(np.diff(tt[pk]))) / 60.0
    P21 = r['P21'][m]; aRc = r['aRc'][m]; Dna = r['Dna'][m]
    inS = (aRc > 0.05) & (Dna < 0.98); inG2 = Dna >= 0.98; preS = ~inS & ~inG2
    dwell = per * (preS & (P21 > P27_HI)).mean() / ((preS | inS | inG2).mean() or 1.0)
    cls = 'transient' if dwell >= DWELL_CUT else 'immediate'
    return (i, shh, f0, cls, dwell, cd, per)


def _arg(flag, default):
    for a in sys.argv:
        if a.startswith(flag + '='): return a.split('=', 1)[1]
    return default


if __name__ == '__main__':
    FRESH = '--fresh' in sys.argv
    N = int(_arg('--n', 400))
    WORKERS = int(_arg('--workers', max(1, mp.cpu_count() - 1)))
    CACHE = 'simulations/sim_steadystate_g0_mark_cache.npz'

    if FRESH or not os.path.exists(CACHE):
        rng = np.random.default_rng(11)
        mu_i = np.clip(0.0005 * np.exp(rng.normal(0, 0.22, N)), 0.00033, 0.00075)
        delmk_i = np.clip(0.0007 * np.exp(rng.normal(0, 0.30, N)), 0.0003, 0.0016)
        kdeez_i = np.clip(0.00015 * np.exp(rng.normal(0, 0.30, N)), 0.00007, 0.00033)
        p21d_i = np.clip(np.exp(rng.normal(np.log(0.45), 0.55, N)), 0.12, 2.0)
        ktl_i = np.clip(0.801 * np.exp(rng.normal(0, 0.33, N)), 0.35, 1.6)
        tasks = []
        for shh in SHH_LEVELS:
            for f0 in (0.233, 1.0):
                for i in range(N):
                    tasks.append((i, shh, f0, mu_i[i], delmk_i[i], kdeez_i[i], p21d_i[i], ktl_i[i]))
        print(f'running {len(tasks)} cells on {WORKERS} workers ...', flush=True)
        with mp.get_context('spawn').Pool(WORKERS, initializer=_init_worker, initargs=(None,)) as pool:
            results = []
            for k, res in enumerate(pool.imap_unordered(_work, tasks, chunksize=8)):
                results.append(res)
                if (k + 1) % 500 == 0: print(f'  {k+1}/{len(tasks)} done', flush=True)
        # aggregate fractions per (SHH, cond)
        cls_codes = {'immediate': 0, 'transient': 1, 'arrest': 2}
        agg = {}
        for f0, tag in [(0.233, 'W'), (1.0, 'N')]:
            frac = np.zeros((len(SHH_LEVELS), 3)); dwell_med = np.full(len(SHH_LEVELS), np.nan); cd_med = np.full(len(SHH_LEVELS), np.nan)
            for si, shh in enumerate(SHH_LEVELS):
                sub = [r for r in results if r[1] == shh and r[2] == f0]
                codes = np.array([cls_codes[r[3]] for r in sub])
                for c in range(3): frac[si, c] = np.mean(codes == c)
                dw = np.array([r[4] for r in sub if r[3] == 'transient'])
                dwell_med[si] = np.nanmedian(dw) if len(dw) else np.nan
                cd_med[si] = np.nanmedian([r[5] for r in sub])
            agg[f'{tag}_frac'] = frac; agg[f'{tag}_dwell'] = dwell_med; agg[f'{tag}_cd'] = cd_med
        np.savez(CACHE, shh=np.array(SHH_LEVELS), N=N, **agg)
        print('cached ->', CACHE, flush=True)

    z = np.load(CACHE)
    shh = z['shh']; Nc = int(z['N'])
    Wf, Nf = z['W_frac'], z['N_frac']

    fig, ax = plt.subplots(2, 2, figsize=(14, 10))
    # (A) transient-G0 fraction vs SHH
    a = ax[0, 0]
    a.plot(shh, Wf[:, 1], '-o', color='#8b1a1a', lw=2.6, ms=7, label='WITH mark (repressed)')
    a.plot(shh, Nf[:, 1], '-s', color='#4a4a4a', lw=2.6, ms=7, label='WITHOUT mark')
    a.set_xlabel('steady-state Hh (SHH)'); a.set_ylabel('fraction TRANSIENT G0')
    a.set_title('(A) Transient-G0 fraction at steady state', fontweight='bold', fontsize=12)
    a.legend(fontsize=10); a.grid(alpha=0.15)
    # (B) arrest fraction vs SHH
    a = ax[0, 1]
    a.plot(shh, Wf[:, 2], '-o', color='#8b1a1a', lw=2.6, ms=7, label='WITH mark')
    a.plot(shh, Nf[:, 2], '-s', color='#4a4a4a', lw=2.6, ms=7, label='WITHOUT mark')
    a.set_xlabel('steady-state Hh (SHH)'); a.set_ylabel('fraction ARRESTED')
    a.set_title('(B) Deep-arrest fraction at steady state', fontweight='bold', fontsize=12)
    a.legend(fontsize=10); a.grid(alpha=0.15)
    # (C) stacked composition — WITH (left offset) vs WITHOUT (right offset)
    a = ax[1, 0]
    labels = ['immediate re-entry', 'transient G0', 'arrest']; cols = ['#2e8b57', '#e08a1e', '#8b1a1a']
    w = 0.018; x = shh
    for grp, frac, off, hatch in [('WITH', Wf, -w, None), ('WITHOUT', Nf, +w, '//')]:
        bottom = np.zeros(len(shh))
        for c in range(3):
            a.bar(x + off, frac[:, c], width=w * 1.8, bottom=bottom, color=cols[c], hatch=hatch,
                  edgecolor='white', lw=0.4, label=(labels[c] if grp == 'WITH' else None))
            bottom += frac[:, c]
    a.set_xlabel('steady-state Hh (SHH)  (left bar=WITH, right hatched=WITHOUT)'); a.set_ylabel('fraction of cells')
    a.set_title('(C) Full composition: immediate / transient-G0 / arrest', fontweight='bold', fontsize=12)
    a.legend(fontsize=9, loc='lower right'); a.set_ylim(0, 1)
    # (D) mechanism: median CyclinD1 (cycling) vs SHH + median transient-G0 dwell
    a = ax[1, 1]
    a.plot(shh, z['W_cd'], '-o', color='#8b1a1a', lw=2.4, ms=6, label='CyclinD1 WITH')
    a.plot(shh, z['N_cd'], '-s', color='#4a4a4a', lw=2.4, ms=6, label='CyclinD1 WITHOUT')
    a.set_xlabel('steady-state Hh (SHH)'); a.set_ylabel('median CyclinD1', color='#333')
    a2 = a.twinx()
    a2.plot(shh, z['W_dwell'], '--^', color='#e08a1e', lw=2.0, ms=6, label='transient-G0 dwell WITH')
    a2.plot(shh, z['N_dwell'], ':v', color='#c9a227', lw=2.0, ms=6, label='dwell WITHOUT')
    a2.set_ylabel('median transient-G0 dwell (h)', color='#c47f17')
    a.set_title('(D) Mechanism: mark lowers CyclinD1 → deeper/longer p27 dwell', fontweight='bold', fontsize=12)
    h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
    a.legend(h1 + h2, l1 + l2, fontsize=8.5, loc='center right'); a.grid(alpha=0.15)

    fig.suptitle(f'Steady-state (constant Hh) transient G0 & arrest — REPRESSED vs NOT repressed (N={Nc}/level/condition)',
                 fontsize=13.5, fontweight='bold', y=1.0)
    plt.tight_layout()
    plt.savefig('simulations/sim_steadystate_g0_mark.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_steadystate_g0_mark.pdf', bbox_inches='tight')
    plt.close()
    for si, s in enumerate(shh):
        print(f'SHH={s:.2f}  WITH[imm/tG0/arr]={Wf[si].round(2)}  WITHOUT={Nf[si].round(2)}')
    print('Saved sim_steadystate_g0_mark.png')
