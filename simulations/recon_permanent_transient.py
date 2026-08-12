"""Recon (JP 2026-08-09 overnight): permanent cell-cycle exit vs transient G0 in GNP vs MB, as a function of
drive-reduction (CDKI or Hedgehog), tracking H3K27me3 (Mk), CyclinD1 (Cd), CDK6 (cdk6), and cycle period.

Baked default = integrated + CDK4/6-alone escape. Hypothesis: the mark accumulates as a cell slows (less
replicative dilution) and represses CyclinD1/CDK6; whether the cell PERMANENTLY exits or TRANSIENTLY arrests
depends on whether the (mark-repressed) drive stays above the escape threshold — GNP (low drive) tips to
permanent exit, MB (high drive) stays transient.

Classify (60000 min): permanent (no division in last third) / transient-G0 (still cycling, period>40h) /
cycling (period<=40h). Run: PYTHONPATH=. ./venv/bin/python simulations/recon_permanent_transient.py
"""
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

M = build_model_v44()   # baked default (integrated + escape)
RR = te.loada(M); RR.integrator.setValue('relative_tolerance', 1e-7); RR.integrator.setValue('absolute_tolerance', 1e-9)
RR.integrator.setValue('maximum_num_steps', 12000000)
try: CD2_DEF = RR['Cd2_expr']
except Exception: CD2_DEF = None


def run(celltype, cdki=1.0, shh=0.5, T=60000):
    for k in ():
        pass
    if celltype == 'GNP':
        RR['SHH'] = shh; RR['Ptch1_copy_number'] = 1.0; p18, p19 = 0.464, 0.36
        if CD2_DEF is not None: RR['Cd2_expr'] = CD2_DEF
    else:
        RR['SHH'] = shh; RR['Ptch1_copy_number'] = 0.3; p18, p19 = 1.73, 0.58
        try: RR['Cd2_expr'] = RR['CD2_EXPR_MB']
        except Exception: pass
    ks = RR['kSyP21']
    RR['p18'] = p18 * cdki; RR['p19'] = p19 * cdki; RR['kSyP21'] = ks * cdki
    RR.reset()
    r = RR.simulate(0, T, T // 10, selections=['time', 'Dna', 'Mk', 'Cd', 'cdk6', 'EZH2'])
    RR['kSyP21'] = ks
    t = r['time']; dna = r['Dna']
    dd = t[np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]]
    ntot = len(dd); nlast = int(np.sum(dd > T * 0.66))
    per = float(np.mean(np.diff(dd)) / 60) if ntot >= 2 else None
    mm = t >= T * 0.7
    state = 'permanent' if nlast == 0 else ('transient-G0' if (per and per > 40) else 'cycling')
    return dict(state=state, ndiv=ntot, period_h=(round(per, 1) if per else None),
                Mk=round(float(r['Mk'][mm].mean()), 3), Cd=round(float(r['Cd'][mm].mean()), 2),
                cdk6=round(float(r['cdk6'][mm].mean()), 2), EZH2=round(float(r['EZH2'][mm].mean()), 2))

print("=== CDKI sweep (SHH=0.5): GNP vs MB — permanent exit vs transient G0 ===")
print(f"{'cell':4} {'CDKI':>5} | {'state':13} {'ndiv':>4} {'period':>7} {'Mk':>6} {'Cd':>6} {'cdk6':>6} {'EZH2':>5}")
for cell in ('GNP', 'MB'):
    for cdki in (1, 2, 3, 4, 5, 6, 8):
        o = run(cell, cdki=cdki)
        print(f"{cell:4} {cdki:>5} | {o['state']:13} {o['ndiv']:>4} {str(o['period_h']):>7} {o['Mk']:>6} {o['Cd']:>6} {o['cdk6']:>6} {o['EZH2']:>5}")
    print()

print("=== Hedgehog sweep (CDKI=1x): GNP vs MB — permanent exit as Hh wanes ===")
print(f"{'cell':4} {'SHH':>5} | {'state':13} {'ndiv':>4} {'period':>7} {'Mk':>6} {'Cd':>6} {'cdk6':>6} {'EZH2':>5}")
for cell in ('GNP', 'MB'):
    for shh in (0.5, 0.35, 0.25, 0.15, 0.08, 0.0):
        o = run(cell, shh=shh)
        print(f"{cell:4} {shh:>5} | {o['state']:13} {o['ndiv']:>4} {str(o['period_h']):>7} {o['Mk']:>6} {o['Cd']:>6} {o['cdk6']:>6} {o['EZH2']:>5}")
    print()
