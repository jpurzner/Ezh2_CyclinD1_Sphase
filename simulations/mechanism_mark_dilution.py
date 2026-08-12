"""Mechanism (JP overnight): does H3K27me3 accumulation/dilution MODULATE the permanent-vs-transient boundary,
or is it set purely by the CyclinD1/CDK6 drive? Turn the mark->CyclinD1/CDK6 repression ON/OFF and vary the
replicative dilution rate (k_dil_S), for GNP and MB across the CDKI drive-reduction axis.

mark OFF = f0_prc2=1 + f0_cdk6=1 (the H3K27me3/PRC2 IFFL cannot repress -> drive at un-repressed max).
dil x2 / x0.5 = faster/slower replicative dilution -> lower/higher mark -> less/more repression.

Baked default. Classify (60000 min): permanent / transient-G0 (period>40h) / cycling.
Run: PYTHONPATH=. ./venv/bin/python simulations/mechanism_mark_dilution.py
"""
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44


def build(extra):
    m = build_model_v44(params=extra) if extra else build_model_v44()
    rr = te.loada(m); rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
    rr.integrator.setValue('maximum_num_steps', 12000000)
    return rr


def classify(rr, cell, cdki, T=60000):
    if cell == 'GNP':
        rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 1.0; p18, p19 = 0.464, 0.36
    else:
        rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3; p18, p19 = 1.73, 0.58
        try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
        except Exception: pass
    ks = rr['kSyP21']
    rr['p18'] = p18 * cdki; rr['p19'] = p19 * cdki; rr['kSyP21'] = ks * cdki
    rr.reset()
    r = rr.simulate(0, T, T // 10, selections=['time', 'Dna'])
    rr['kSyP21'] = ks
    t = r['time']; dna = r['Dna']; dd = t[np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]]
    nlast = int(np.sum(dd > T * 0.66)); per = np.mean(np.diff(dd)) / 60 if len(dd) >= 2 else None
    return 'P' if nlast == 0 else ('T' if (per and per > 40) else 'C')   # Permanent / Transient / Cycling


CONFIGS = [
    ('baked (mark on)', {}),
    ('mark OFF (f0=1)', {'f0_prc2': 1.0, 'f0_cdk6': 1.0}),
    ('dilution x2', {'k_dil_S': 1.386}),
    ('dilution x0.5', {'k_dil_S': 0.3465}),
    ('dilution OFF (0)', {'k_dil_S': 0.0}),
]
cdkis = [1, 2, 3, 4, 5, 6, 7, 8, 10]
print("State vs CDKI (C=cycling, T=transient-G0, P=permanent). Boundary = where it turns P.")
for cell in ('GNP', 'MB'):
    print(f"\n--- {cell} ---   CDKI: " + " ".join(f"{c:>2}" for c in cdkis))
    for name, extra in CONFIGS:
        rr = build(extra)
        row = " ".join(f"{classify(rr, cell, c):>2}" for c in cdkis)
        print(f"  {name:20} {row}")
