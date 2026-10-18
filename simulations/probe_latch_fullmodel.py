"""What happens to the FULL v44 model in the LATCH regime? (JP 2026-07-21)
Latch = low nucleation a0_prc2 + cooperative read-write n_rw>=2. Definitive bistability test:
run the SAME condition from a HIGH-mark and a LOW-mark initial state; if they settle to DIFFERENT
marks the mark is BISTABLE (a latch), if they converge it is monostable.
Run: ./venv/bin/python simulations/probe_latch_fullmodel.py"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1'); os.environ.setdefault('EZH2_CONC', '0')
import tellurium as te
from src.build_model_v44_heldt import build_model_v44
UPH = 62.8

COND = {
    'GNP+SHH': dict(SHH=0.5, HHi=0, MYCN_amplification=1.0, Ptch1_copy_number=1.0),
    'MB':      dict(SHH=0.5, HHi=0, MYCN_amplification=2.8, Ptch1_copy_number=0.1,
                   p16=0.306, p18=1.553, kSyP21=0.002, P21_div=0.6),
    'GNP+HHi': dict(SHH=0.5, HHi=1, MYCN_amplification=1.0, Ptch1_copy_number=1.0),
}
HIGH = dict(Mk=0.90, m2_me=0.06, m1_me=0.03)   # start heavily methylated
LOW = dict(Mk=0.02, m2_me=0.05, m1_me=0.05)    # start de-methylated

def run(model, cond, init, T=20000):
    rr = te.loada(model); rr.reset()
    for k, v in cond.items(): rr[k] = v
    for k, v in init.items():
        try: rr[k] = v
        except Exception: pass
    for atol in (1e-8, 1e-7, 1e-6):
        try:
            rr.integrator.setValue('absolute_tolerance', atol)
            r = rr.simulate(0, T, 4000, selections=['time', 'Mk', 'Cd', 'Dna', 'EZH2'])
            m = np.asarray(r['time']) > T * 0.55
            dna = np.asarray(r['Dna'])
            divs = int(np.sum((dna[:-1] > 0.9) & (dna[1:] < 0.1)))
            return dict(Mk=float(np.asarray(r['Mk'])[m].mean()), Cd=float(np.asarray(r['Cd'])[m].mean()),
                        EZH2=float(np.asarray(r['EZH2'])[m].mean()), divs=divs)
        except Exception:
            rr.reset()
            for k, v in cond.items(): rr[k] = v
            for k, v in init.items():
                try: rr[k] = v
                except Exception: pass
    return None

def probe(name, over):
    print(f"\n=== {name}  {over} ===")
    model = build_model_v44(params=over)
    print(f"  {'cond':9s} | {'Mk(hi)':>7} {'Mk(lo)':>7} {'Δ':>6} bistable? | {'Cd(hi)':>7} {'Cd(lo)':>7} | {'div(hi)':>7} {'div(lo)':>7}")
    for c, cd in COND.items():
        h = run(model, cd, HIGH); l = run(model, cd, LOW)
        if h is None or l is None:
            print(f"  {c:9s} | INTEGRATION FAILED"); continue
        d = abs(h['Mk'] - l['Mk']); bi = 'YES <<<' if d > 0.10 else 'no'
        print(f"  {c:9s} | {h['Mk']:7.3f} {l['Mk']:7.3f} {d:6.3f} {bi:8s} | {h['Cd']:7.3f} {l['Cd']:7.3f} | {h['divs']:7d} {l['divs']:7d}")

if __name__ == '__main__':
    probe('BAKED (monostable reference)', {})
    for a0 in (1e-4, 3e-5, 1e-5):
        probe(f'LATCH a0={a0}, n_rw=3', dict(a0_prc2=a0, n_rw=3))
    probe('LATCH a0=1e-5, n_rw=3, a_rw up', dict(a0_prc2=1e-5, n_rw=3, a_rw_prc2=0.02))
