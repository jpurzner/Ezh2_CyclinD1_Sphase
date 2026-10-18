"""Probe the poised bivalent differentiation gene (JP 2026-07-21). Normal GNP: Mk_diff high, Diff OFF
(repressed/poised). Push division faster (raise mu -> shorter Tc) or block the writer (EZH2i): Mk_diff
dilutes out -> Diff PREMATURELY ACTIVATES. Run: ./venv/bin/python simulations/probe_diff_gene.py"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1'); os.environ.setdefault('EZH2_CONC', '0')
import tellurium as te
from src.build_model_v44_heldt import build_model_v44
UPH = 62.8
GNP = dict(SHH=0.5, HHi=0, MYCN_amplification=1.0, Ptch1_copy_number=1.0)

# tuning candidate: SLOW writing (restoration timescale ~ Tc, the me3 de-novo lag) -> rate-sensitive mark;
# sharp poised switch (n_diff, k_on) so the gene fires decisively once H3K27me3 crosses below H3K4me3
DIFF = dict(a0_diff=0.00045, a_rw_diff=0.0026, del_diff=0.0003, k_dil_diff=0.693,
            K_diff=0.62, n_diff=10, k_on_diff=0.05, k_off_diff=0.015)
MODEL = build_model_v44(with_diff_gene=True)  # baked diff-gene defaults

def run(over, T=22000):
    rr = te.loada(MODEL); rr.reset()
    for k, v in GNP.items(): rr[k] = v
    for k, v in over.items(): rr[k] = v
    for atol in (1e-8, 1e-7, 1e-6):
        try:
            rr.integrator.setValue('absolute_tolerance', atol)
            r = rr.simulate(0, T, 4000, selections=['time', 'Mk', 'Mk_diff', 'Diff', 'Dna', 'EZH2'])
            m = np.asarray(r['time']) > T * 0.55
            dna = np.asarray(r['Dna']); dv = np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]
            t = np.asarray(r['time'])
            Tc = float(np.mean(np.diff(t[dv])) / UPH) if len(dv) >= 2 else np.inf
            return dict(Mk=float(np.asarray(r['Mk'])[m].mean()), Mk_diff=float(np.asarray(r['Mk_diff'])[m].mean()),
                        Diff=float(np.asarray(r['Diff'])[m].mean()), EZ=float(np.asarray(r['EZH2'])[m].mean()), Tc=Tc,
                        ndiv=len(dv))
        except Exception:
            rr.reset()
            for k, v in GNP.items(): rr[k] = v
            for k, v in over.items(): rr[k] = v
    return None

if __name__ == '__main__':
    print("Poised bivalent differentiation gene in the full v44 model (GNP+SHH).\n")
    print(f"  {'condition':28s} {'Tc(h)':>6} {'Mk_cd':>6} {'Mk_diff':>8} {'Diff':>6}  interpretation")
    cases = [('normal GNP (mu=0.0005)', dict()),
             ('faster mu=0.0008', dict(mu=0.0008)),
             ('faster mu=0.0011', dict(mu=0.0011)),
             ('faster mu=0.0015', dict(mu=0.0015)),
             ('very fast mu=0.0022', dict(mu=0.0022)),
             ('EZH2i (writer blocked)', dict(EZH2i=1)),
             ('weak writer kTlEZ x0.4', dict(kTlEZ=0.006)),]
    for name, ov in cases:
        r = run(ov)
        if r is None:
            print(f"  {name:28s}  FAILED"); continue
        interp = 'PREMATURE DIFFERENTIATION' if r['Diff'] > 0.5 else ('poised (repressed)' if r['Mk_diff'] > 0.6 else 'mark falling')
        print(f"  {name:28s} {r['Tc']:>6.1f} {r['Mk']:>6.3f} {r['Mk_diff']:>8.3f} {r['Diff']:>6.3f}  {interp}")
