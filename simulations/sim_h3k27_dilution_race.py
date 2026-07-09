"""THE RACE: H3K27me3 on CyclinD1 as a competition between EZH2/PRC2 activity (methylation, fills the mark)
and REPLICATIVE DILUTION (mark halved every S; frequency = division rate). JP concept: with MORE EZH2
activity methylation OVERWHELMS dilution -> H3K27me3 accumulates -> CyclinD1 dampened (MB-like); with
INSUFFICIENT EZH2 dilution wins -> mark stays low -> CyclinD1 high (GNP-like). Slower division (MB's high
p16/p18/p21 brake) ALSO tips the race toward accumulation by removing dilution events. This is a GRADED
mechanism (not the self-reinforcing switch that broke the MB/GNP fold): the mark ends up HIGHER where EZH2
is higher and division slower, which is exactly the dose-graded feedback that yields CyclinD1 MB/GNP=5.07.

Phase map over EZH2 activity (methylation-rate x, relative to GNP baseline) x division period (h, set by
growth rate mu). Color = steady-state H3K27me3 (Mk) on CyclinD1 and the resulting CyclinD1 level. GNP and MB
operating points overlaid. Weak-self-reinforcement (baseline k_w_mk) regime = JP's rate-race, not a switch.

Run:  ./venv/bin/python simulations/sim_h3k27_dilution_race.py [--na=9] [--nmu=9] [--workers=9] [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tellurium as te
import multiprocessing as mp
from src.build_model_v44_heldt import build_model_v44

SEL = ['time', 'Mk', 'Cd_mRNA', 'EZH2', 'Dna']
GNP = dict(SHH=0.5, Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)
KW0, K00 = 0.00449887, 0.000540513          # baseline methylation-rate constants (scaled together by activity A)
SETTLE, REC, NP = 6000.0, 14000.0, 4000
CACHE = 'simulations/sim_h3k27_dilution_race_cache.npz'
_RR = None


def _init(_):
    global _RR
    _RR = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    _RR.integrator.setValue('relative_tolerance', 1e-6)
    try:
        _RR.integrator.setValue('maximum_num_steps', 1000000); _RR.integrator.setValue('maximum_time_step', 20.0)
    except Exception:
        pass


def _work(task):
    A, muv = task
    _RR.reset()
    for k, v in GNP.items():
        _RR[k] = v
    _RR['k_w_mk'] = A * KW0; _RR['k0_mk'] = A * K00; _RR['mu'] = float(muv)
    try:
        _RR.simulate(0, SETTLE, 2)
        d = _RR.simulate(SETTLE, SETTLE + REC, NP, selections=SEL)
    except Exception:
        return (A, muv, np.nan, np.nan, np.nan, 0)
    Dna = d['Dna']; div = np.where((Dna[:-1] > 0.9) & (Dna[1:] < 0.1))[0]
    t = d['time']; per = float(np.mean(np.diff(t[div])) / 60) if len(div) > 1 else np.nan
    return (A, muv, float(np.mean(d['Mk'])), float(np.mean(d['Cd_mRNA'])), per, int(len(div)))


def _arg(f, d):
    for a in sys.argv:
        if a.startswith(f + '='):
            return type(d)(a.split('=', 1)[1])
    return d


if __name__ == '__main__':
    NA = _arg('--na', 9); NMU = _arg('--nmu', 9); WORKERS = _arg('--workers', max(1, mp.cpu_count() - 1))
    A_VALS = np.round(np.exp(np.linspace(np.log(0.3), np.log(3.5), NA)), 3)   # EZH2/PRC2 activity x GNP baseline
    MU_VALS = np.round(np.exp(np.linspace(np.log(0.00030), np.log(0.00110), NMU)), 6)  # growth -> division rate
    if '--fresh' in sys.argv or not os.path.exists(CACHE):
        tasks = [(float(A), float(mu)) for A in A_VALS for mu in MU_VALS]
        print(f'{len(tasks)} cells: EZH2 activity {A_VALS[0]}-{A_VALS[-1]}x  x  mu {MU_VALS[0]}-{MU_VALS[-1]}', flush=True)
        with mp.get_context('spawn').Pool(WORKERS, initializer=_init, initargs=(None,)) as pool:
            res = list(pool.imap_unordered(_work, tasks, chunksize=4))
        MK = np.full((NA, NMU), np.nan); CD = np.full((NA, NMU), np.nan); PER = np.full((NA, NMU), np.nan); NDIV = np.zeros((NA, NMU))
        ai = {round(a, 3): i for i, a in enumerate(A_VALS)}; mj = {round(m, 6): j for j, m in enumerate(MU_VALS)}
        for A, mu, mk, cd, per, nd in res:
            i = ai[round(A, 3)]; j = mj[round(mu, 6)]; MK[i, j] = mk; CD[i, j] = cd; PER[i, j] = per; NDIV[i, j] = nd
        np.savez(CACHE, A_VALS=A_VALS, MU_VALS=MU_VALS, MK=MK, CD=CD, PER=PER, NDIV=NDIV)
        print('cached ->', CACHE, flush=True)

    z = np.load(CACHE)
    A_VALS, MU_VALS, MK, CD, PER, NDIV = z['A_VALS'], z['MU_VALS'], z['MK'], z['CD'], z['PER'], z['NDIV']
    # median period per mu column (for the y-axis in hours); use nan-robust
    per_axis = np.nanmedian(PER, axis=0)
    # operating points (from probe): GNP A=1.0 period~21.8h; MB EZH2 2.28/1.15=1.98, period~23.2h
    GNP_PT = (1.0, 21.8); MB_PT = (1.98, 23.2)

    fig, ax = plt.subplots(1, 2, figsize=(15.5, 6))
    Amesh, Pmesh = np.meshgrid(A_VALS, per_axis, indexing='ij')
    # (A) steady-state H3K27me3
    a = ax[0]
    cf = a.contourf(Amesh, Pmesh, MK, levels=np.linspace(0, 1, 21), cmap='viridis', extend='max')
    cs = a.contour(Amesh, Pmesh, MK, levels=[0.5], colors='white', linewidths=2.2, linestyles='--')
    a.clabel(cs, fmt='race balance (Mk=0.5)', fontsize=8)
    for (x, y), nm, c in [(GNP_PT, 'GNP', '#00e5ff'), (MB_PT, 'MB', '#ff3b3b')]:
        a.plot(x, y, 'o', ms=13, mfc='none', mec=c, mew=2.6); a.annotate(nm, (x, y), textcoords='offset points', xytext=(8, 6), color=c, fontweight='bold', fontsize=11)
    a.annotate('', MB_PT, GNP_PT, arrowprops=dict(arrowstyle='->', color='white', lw=1.4, alpha=0.7))
    a.set_xscale('log'); a.set_xlabel('EZH2 / PRC2 activity  (methylation rate, × GNP)'); a.set_ylabel('division period (h)  ← slower ·  faster →')
    a.set_title('(A) H3K27me3 on CyclinD1 = methylation vs dilution race\nhigh EZH2 + slow division → accumulation; low EZH2 + fast division → dilution wins', fontweight='bold', fontsize=10.5)
    fig.colorbar(cf, ax=a, label='steady-state mark Mk')
    # (B) resulting CyclinD1 (log) -- silencing corner
    a = ax[1]
    cf = a.contourf(Amesh, Pmesh, np.log10(CD), levels=20, cmap='magma')
    for (x, y), nm, c in [(GNP_PT, 'GNP', '#00e5ff'), (MB_PT, 'MB', '#39ff14')]:
        a.plot(x, y, 'o', ms=13, mfc='none', mec=c, mew=2.6); a.annotate(nm, (x, y), textcoords='offset points', xytext=(8, 6), color=c, fontweight='bold', fontsize=11)
    a.set_xscale('log'); a.set_xlabel('EZH2 / PRC2 activity  (× GNP)'); a.set_ylabel('division period (h)')
    a.set_title('(B) Resulting CyclinD1 -- damped in the high-EZH2 / slow-division corner\n(H3K27me3 accumulation), high where dilution wins', fontweight='bold', fontsize=10.5)
    cb = fig.colorbar(cf, ax=a); cb.set_label('log10 CyclinD1 (Cd_mRNA)')
    fig.suptitle('The H3K27me3-on-CyclinD1 RACE: EZH2 activity vs replicative dilution (GNP baseline, weak-self-reinforcement / graded regime)', fontsize=12.5, fontweight='bold', y=1.0)
    plt.tight_layout()
    plt.savefig('simulations/sim_h3k27_dilution_race.png', dpi=150, bbox_inches='tight')
    plt.savefig('simulations/sim_h3k27_dilution_race.pdf', bbox_inches='tight')
    print('Saved sim_h3k27_dilution_race.png')
    print('period axis (h):', np.round(per_axis, 1))
    print('Mk at GNP-ish (A=1):', np.round(MK[np.argmin(np.abs(A_VALS - 1.0))], 3))
    print('Mk at MB-ish (A~2):', np.round(MK[np.argmin(np.abs(A_VALS - 1.98))], 3))
