"""Fast probe of the four 2026-07-21 mark-mechanism additions (JP): nerfed eraser, cooperative
read-write (n_rw), continuous S-dilution (k_dil_S), H3.3 turnover (k_h33). Reports the mark level,
MB/GNP CyclinD1 fold, ChIP direction (Mk MB vs GNP), and mark half-life for a given override set.
Run: ./venv/bin/python simulations/probe_four_mechanisms.py [name]"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1'); os.environ.setdefault('EZH2_CONC', '0')
import tellurium as te
from src.build_model_v44_heldt import build_model_v44
UPH = 62.8

COND = {  # match validate_v44.py exactly
    'GNP+SHH': dict(SHH=0.5, HHi=0, MYCN_amplification=1.0, Ptch1_copy_number=1.0),
    'MB':      dict(SHH=0.5, HHi=0, MYCN_amplification=2.8, Ptch1_copy_number=0.1,
                   p16=0.306, p18=1.553, kSyP21=0.002, P21_div=0.6),
    'GNP+HHi': dict(SHH=0.5, HHi=1, MYCN_amplification=1.0, Ptch1_copy_number=1.0),
}

def settled(rr, sp, t0=8000, T=16000):
    r = rr.simulate(0, T, 3000, selections=['time'] + sp)
    m = np.asarray(r['time']) > t0
    return {s: float(np.asarray(r[s])[m].mean()) for s in sp}

def run(over, cond):
    m = build_model_v44(params=over)
    rr = te.loada(m); rr.reset()
    for k, v in cond.items(): rr[k] = v
    for atol in (1e-8, 1e-7, 1e-6):
        try:
            rr.integrator.setValue('absolute_tolerance', atol)
            return settled(rr, ['Mk', 'Cd', 'EZH2', 'Cd_mRNA', 'Gli1'])
        except Exception:
            rr.reset()
            for k, v in cond.items(): rr[k] = v
    return None

def probe(name, over):
    print(f"\n=== {name}  overrides={over} ===")
    res = {c: run(over, cond) for c, cond in COND.items()}
    if any(v is None for v in res.values()):
        print("  INTEGRATION FAILED"); return
    fold = res['MB']['Cd'] / res['GNP+SHH']['Cd']
    chip = res['MB']['Mk'] / res['GNP+SHH']['Mk']
    print(f"  {'cond':10s} {'Mk':>7} {'Cd':>8} {'EZH2':>7} {'Cd_mRNA':>8} {'Gli1':>7}")
    for c, v in res.items():
        print(f"  {c:10s} {v['Mk']:7.3f} {v['Cd']:8.3f} {v['EZH2']:7.3f} {v['Cd_mRNA']:8.4f} {v['Gli1']:7.3f}")
    print(f"  --> CyclinD1 MB/GNP fold = {fold:.3f} (target 5.07)   Mk MB/GNP (ChIP, want <1) = {chip:.3f}")
    return dict(fold=fold, chip=chip, res=res)

# current baked = all four ENGAGED. OFF = revert to pre-2026-07-21 config.
OFF = dict(k_jmjd3_gli=0.0, n_rw=1, dil_frac=0.5, k_dil_S=0.0, k_h33=0.0)

if __name__ == '__main__':
    probe('ENGAGED (baked: eraser 0.005, n_rw 2, cont-dil, h33 0.003)', {})
    probe('OFF (pre-2026-07-21 config)', OFF)
    # eraser-strength sweep to show ChIP MB<GNP sensitivity
    print("\n--- eraser strength sweep (ChIP MB/GNP direction) ---")
    for kj in (0.0, 0.005, 0.02, 0.05, 0.1):
        r = probe(f'k_jmjd3_gli={kj}', dict(k_jmjd3_gli=kj))
