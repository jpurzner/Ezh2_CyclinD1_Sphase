"""Poised bivalent differentiation gene -> premature activation -> cell-cycle EXIT (JP 2026-07-21).
Diff now feeds back: it represses CyclinD1 (proneural-TF-like) AND actively erases its own H3K27me3
(bivalent resolution / KDM6B) -> once fired it LATCHES and the cell exits (post-mitotic). Generates
time courses (normal / fast-division / EZH2i) with cumulative divisions + CyclinD1, and the
division-rate dose-response (final Diff + divisions-after). Run: ./venv/bin/python simulations/sim_diff_gene.py"""
import os, sys, pickle, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1'); os.environ.setdefault('EZH2_CONC', '0')
import tellurium as te
from src.build_model_v44_heldt import build_model_v44
HERE = os.path.dirname(os.path.abspath(__file__)); UPH = 62.8
GNP = dict(SHH=0.5, HHi=0, MYCN_amplification=1.0, Ptch1_copy_number=1.0)
MODEL = build_model_v44(with_diff_gene=True)
SEL = ['time', 'me1_diff', 'me2_diff', 'Mk_diff', 'Diff', 'Dna', 'Cd']

def div_idx(dna): dna = np.asarray(dna); return np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]

def protocol(segs):
    rr = te.loada(MODEL); rr.reset()
    for k, v in GNP.items(): rr[k] = v
    rr.integrator.setValue('absolute_tolerance', 1e-7)
    T, M1, M2, MKD, DF, DNA, CD = [], [], [], [], [], [], []
    t0 = 0.0
    for ov, dur in segs:
        for k, v in ov.items(): rr[k] = v
        r = rr.simulate(t0, t0 + dur, int(dur / 5), selections=SEL)
        T += list(r['time']); M1 += list(r['me1_diff']); M2 += list(r['me2_diff'])
        MKD += list(r['Mk_diff']); DF += list(r['Diff']); DNA += list(r['Dna']); CD += list(r['Cd'])
        t0 += dur
    T = np.array(T) / UPH
    dv = div_idx(DNA)
    cumdiv = np.zeros(len(T))
    for i in dv: cumdiv[i:] += 1
    return dict(t=T, me1=np.array(M1), me2=np.array(M2), Mk_diff=np.array(MKD), Diff=np.array(DF), Cd=np.array(CD), cumdiv=cumdiv, div_t=T[dv])

def settle(ov, T=30000):
    rr = te.loada(MODEL); rr.reset()
    for k, v in GNP.items(): rr[k] = v
    for k, v in ov.items(): rr[k] = v
    for atol in (1e-8, 1e-7, 1e-6):
        try:
            rr.integrator.setValue('absolute_tolerance', atol)
            r = rr.simulate(0, T, 5000, selections=SEL); m = np.asarray(r['time']) > T * 0.55
            dv = div_idx(r['Dna']); t = np.asarray(r['time'])
            # divisions in the LAST third (rate after settling)
            late = t[dv] > T * 0.66
            Tc = float(np.mean(np.diff(t[dv])) / UPH) if len(dv) >= 2 else np.inf
            return dict(Mk_diff=float(np.asarray(r['Mk_diff'])[m].mean()), Diff=float(np.asarray(r['Diff'])[m].mean()),
                        Tc=Tc, div_late=int(np.sum(late)))
        except Exception:
            rr.reset()
            for k, v in GNP.items(): rr[k] = v
            for k, v in ov.items(): rr[k] = v
    return dict(Mk_diff=np.nan, Diff=np.nan, Tc=np.nan, div_late=np.nan)

if __name__ == '__main__':
    res = {}
    ONSET = 10000
    res['normal'] = protocol([(dict(), 40000)])
    res['fast'] = protocol([(dict(), ONSET), (dict(mu=0.0013), 34000)])
    res['ezh2i'] = protocol([(dict(), ONSET), (dict(EZH2i=1), 34000)])
    res['onset_h'] = ONSET / UPH
    for k in ('normal', 'fast', 'ezh2i'):
        d = res[k]; print(f"  {k:7s}: total div={int(d['cumdiv'][-1])}  final Diff={d['Diff'][-1]:.2f} Mk_diff={d['Mk_diff'][-1]:.2f} Cd={d['Cd'][-1]:.2f}")

    print("division-rate dose-response...")
    mus = [0.0004, 0.0005, 0.0006, 0.0007, 0.00078, 0.00086, 0.00095, 0.0011, 0.0013, 0.0016]
    sw = [settle(dict(mu=mu)) for mu in mus]
    res['sweep'] = dict(mu=np.array(mus), Tc=np.array([s['Tc'] for s in sw]), Mk_diff=np.array([s['Mk_diff'] for s in sw]),
                        Diff=np.array([s['Diff'] for s in sw]), div_late=np.array([s['div_late'] for s in sw]))
    pickle.dump(res, open(os.path.join(HERE, 'sim_diff_gene_cache.pkl'), 'wb'))
    print("-> cached. figure: fig_diff_gene.py")
    print("  sweep Tc/Diff/div_late:", [f"{tc:.0f}h:D{d:.2f}:n{n}" for tc, d, n in zip(res['sweep']['Tc'], res['sweep']['Diff'], res['sweep']['div_late'])])
