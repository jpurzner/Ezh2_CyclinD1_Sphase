"""Recon (JP): how the model behaves when EZH2 generation is SUFFICIENT to compensate replicative dilution
of H3K27me3, IN EXCESS, or FAILING to keep up. Knob = kTlEZ (EZH2 translation = amount of EZH2 generated),
scaled by ez_scale (baked default = 1.0). Baked default.

Regime read-out (steady cycling): mean Mk near 1 = EXCESS (writing >> dilution, mark saturated); mid =
SUFFICIENT (balanced); low = FAILING (dilution wins, mark eroded). Division-RATE sensitivity = Mk(slow
cycle, 4x CDKI) - Mk(fast cycle, 1x CDKI): large gap => the mark is a division-rate sensor (dilution bites);
~0 => saturated/insensitive. Also track EZH2, CyclinD1 (Cd), proliferation.

Run: PYTHONPATH=. ./venv/bin/python simulations/recon_ezh2_dilution_balance.py
"""
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

KTLEZ0 = 0.015098862630302894  # baked default kTlEZ


def build(ez_scale):
    rr = te.loada(build_model_v44(params={'kTlEZ': KTLEZ0 * ez_scale}))
    rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
    rr.integrator.setValue('maximum_num_steps', 12000000)
    return rr


def measure(rr, cell, cdki=1.0, T=45000):
    if cell == 'GNP':
        rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 1.0; p18, p19 = 0.464, 0.36
    else:
        rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3; p18, p19 = 1.73, 0.58
        try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
        except Exception: pass
    ks = rr['kSyP21']; rr['p18'] = p18 * cdki; rr['p19'] = p19 * cdki; rr['kSyP21'] = ks * cdki
    rr.reset()
    r = rr.simulate(0, T, T // 10, selections=['time', 'Dna', 'Mk', 'EZH2', 'Cd', 'cdk6'])
    rr['kSyP21'] = ks
    t = r['time']; mm = t >= T * 0.6; dd = t[np.where((r['Dna'][:-1] > 0.9) & (r['Dna'][1:] < 0.1))[0]]
    per = float(np.mean(np.diff(dd)) / 60) if len(dd) >= 2 else None
    nlast = int(np.sum(dd > T * 0.66))
    return dict(Mk=float(r['Mk'][mm].mean()), EZH2=float(r['EZH2'][mm].mean()),
                Cd=float(r['Cd'][mm].mean()), cdk6=float(r['cdk6'][mm].mean()),
                per=per, ndiv=len(dd), cycling=(nlast > 0))


print("EZH2-generation (kTlEZ) sweep — writer vs dilution regime (MB, then GNP).")
print("regime: Mk~1 EXCESS | mid SUFFICIENT | low FAILING;  rate-sens = Mk(4x)-Mk(1x)")
for cell in ('MB', 'GNP'):
    print(f"\n--- {cell} ---")
    print(f"{'ez_scale':>8} {'EZH2':>6} {'Mk(1x)':>7} {'Mk(4x)':>7} {'rate-sens':>9} {'period':>7} {'Cd':>6} {'cdk6':>6} {'cyc?':>5}")
    for ez in (0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0):
        rr = build(ez)
        o1 = measure(rr, cell, cdki=1.0)
        rr = build(ez)
        o4 = measure(rr, cell, cdki=4.0)
        rs = o4['Mk'] - o1['Mk']
        reg = 'EXCESS' if o1['Mk'] > 0.85 else ('SUFFIC' if o1['Mk'] > 0.45 else 'FAILING')
        print(f"{ez:>8} {o1['EZH2']:>6.2f} {o1['Mk']:>7.3f} {o4['Mk']:>7.3f} {rs:>9.3f} "
              f"{str(round(o1['per'],1) if o1['per'] else None):>7} {o1['Cd']:>6.2f} {o1['cdk6']:>6.2f} {str(o1['cycling']):>5}  {reg}")
