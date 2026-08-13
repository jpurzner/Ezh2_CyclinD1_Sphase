"""EZH2 writer vs replicative-dilution regimes for H3K27me3 (JP 2026-08-09): does EZH2 generation SUFFICE to
compensate dilution, run in EXCESS, or FAIL to keep up? Knob = kTlEZ (EZH2 amount), ez_scale (baked=1.0).
Baked default.

  (A) Mark level vs EZH2 generation -> the three regimes (FAILING: mark eroded; SUFFICIENT: balanced;
      EXCESS: saturated near 1). Baked model (ez=1) sits in EXCESS for both GNP and MB.
  (B) Division-RATE sensitivity of the mark [Mk(slow 4x) - Mk(fast 1x)] vs EZH2 generation: the mark is a
      division-rate SENSOR (speed-limiter) only when the writer is NOT in excess -- strong in failing/
      sufficient, ~0 in excess. This is WHY dilution doesn't move the fate at the baked (excess) setpoint.
  (C) Consequence for the governor: CyclinD1 vs EZH2 generation -- de-repressed when failing, floored in
      excess (diminishing returns; the baked model is past the useful-repression knee).
  (D) The rate-sensor curves: Mk vs cycle period (set by CDKI) for failing / sufficient / excess -- steep
      (dilution bites) vs flat (saturated).

Run:  ./venv/bin/python simulations/fig_v44_ezh2_dilution_regimes.py   (~20 min)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from src.build_model_v44_heldt import build_model_v44

KTLEZ0 = 0.015098862630302894
C_MB = '#6a2c8f'; C_GNP = '#117a65'; C_MK = '#7a5aa6'
plt.rcParams.update({'font.size': 10, 'axes.titlesize': 10.5, 'axes.labelsize': 10, 'axes.titleweight': 'bold'})


def build(ez):
    rr = te.loada(build_model_v44(params={'kTlEZ': KTLEZ0 * ez}))
    rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
    rr.integrator.setValue('maximum_num_steps', 12000000)
    return rr


def measure(rr, cell, cdki=1.0, T=42000):
    if cell == 'GNP':
        rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 1.0; p18, p19 = 0.464, 0.36
    else:
        rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3; p18, p19 = 1.73, 0.58
        try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
        except Exception: pass
    ks = rr['kSyP21']; rr['p18'] = p18 * cdki; rr['p19'] = p19 * cdki; rr['kSyP21'] = ks * cdki
    rr.reset(); r = rr.simulate(0, T, T // 10, selections=['time', 'Dna', 'Mk', 'Cd'])
    rr['kSyP21'] = ks; mm = r['time'] >= T * 0.6
    dd = r['time'][np.where((r['Dna'][:-1] > 0.9) & (r['Dna'][1:] < 0.1))[0]]
    per = float(np.mean(np.diff(dd)) / 60) if len(dd) >= 2 else np.nan
    return float(r['Mk'][mm].mean()), float(r['Cd'][mm].mean()), per


ez_grid = [0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]
data = {c: {'Mk1': [], 'Mk4': [], 'Cd': []} for c in ('MB', 'GNP')}
for c in ('MB', 'GNP'):
    for ez in ez_grid:
        m1, cd1, _ = measure(build(ez), c, 1.0)
        m4, _, _ = measure(build(ez), c, 4.0)
        data[c]['Mk1'].append(m1); data[c]['Mk4'].append(m4); data[c]['Cd'].append(cd1)

fig = plt.figure(figsize=(15, 9)); gs = GridSpec(2, 2, figure=fig, hspace=0.40, wspace=0.28)

# (A) mark level vs EZH2 generation -> regimes
axA = fig.add_subplot(gs[0, 0])
axA.axhspan(0, 0.45, color='#d9534f', alpha=0.08); axA.axhspan(0.45, 0.85, color='#e0a800', alpha=0.08)
axA.axhspan(0.85, 1.0, color='#117a65', alpha=0.08)
axA.text(3.05, 0.30, 'FAILING', color='#a02c24', fontsize=8, ha='right', style='italic')
axA.text(3.05, 0.62, 'SUFFICIENT', color='#8a6d00', fontsize=8, ha='right', style='italic')
axA.text(3.05, 0.92, 'EXCESS', color='#0b5', fontsize=8, ha='right', style='italic')
axA.plot(ez_grid, data['MB']['Mk1'], 'o-', color=C_MB, lw=1.8, label='MB')
axA.plot(ez_grid, data['GNP']['Mk1'], 's-', color=C_GNP, lw=1.8, label='GNP')
axA.axvline(1.0, color='k', ls='--', lw=1); axA.text(1.03, 0.15, 'baked', fontsize=8)
axA.set_xscale('log'); axA.set_xlabel('EZH2 generation  (kTlEZ scale)'); axA.set_ylabel('mean H3K27me3 (Mk, cycling)')
axA.set_ylim(0, 1.02); axA.set_title('(A) EZH2 generation → mark regime (baked = EXCESS)'); axA.legend(fontsize=9)

# (B) division-rate sensitivity (speed-limiter) vs EZH2 generation
axB = fig.add_subplot(gs[0, 1])
rs_mb = np.array(data['MB']['Mk4']) - np.array(data['MB']['Mk1'])
rs_gnp = np.array(data['GNP']['Mk4']) - np.array(data['GNP']['Mk1'])
axB.plot(ez_grid, rs_mb, 'o-', color=C_MB, lw=1.8, label='MB')
axB.plot(ez_grid, rs_gnp, 's-', color=C_GNP, lw=1.8, label='GNP')
axB.axvline(1.0, color='k', ls='--', lw=1); axB.text(1.03, max(rs_gnp.max(), rs_mb.max())*0.9, 'baked', fontsize=8)
axB.set_xscale('log'); axB.set_xlabel('EZH2 generation  (kTlEZ scale)')
axB.set_ylabel('rate-sensitivity  Mk(slow) − Mk(fast)')
axB.set_title('(B) Mark = division-rate sensor only when writer NOT in excess')
axB.legend(fontsize=9)
axB.annotate('strong sensor\n(dilution bites)', (0.1, rs_gnp[0]), (0.14, rs_gnp[0]*0.7), fontsize=8,
             color='#a02c24', arrowprops=dict(arrowstyle='->', color='#a02c24'))
axB.annotate('saturated:\nweak sensor', (1.0, rs_mb[5]), (1.3, rs_mb[5]+0.03), fontsize=8, color='#0b5',
             arrowprops=dict(arrowstyle='->', color='#0b5'))

# (C) governor strength: CyclinD1 vs EZH2 generation
axC = fig.add_subplot(gs[1, 0])
axC.plot(ez_grid, data['MB']['Cd'], 'o-', color=C_MB, lw=1.8, label='MB CyclinD1')
axC.plot(ez_grid, data['GNP']['Cd'], 's-', color=C_GNP, lw=1.8, label='GNP CyclinD1')
axC.axvline(1.0, color='k', ls='--', lw=1); axC.text(1.03, max(data['MB']['Cd'])*0.6, 'baked', fontsize=8)
axC.set_xscale('log'); axC.set_xlabel('EZH2 generation  (kTlEZ scale)'); axC.set_ylabel('CyclinD1 (Cd, cycling)')
axC.set_title('(C) Governor: CyclinD1 de-repressed (failing) → floored (excess)')
axC.legend(fontsize=9)
axC.annotate('diminishing returns\n(past the knee)', (1.5, data['MB']['Cd'][6]), (0.35, data['MB']['Cd'][6]+12),
             fontsize=8, color='#666', arrowprops=dict(arrowstyle='->', color='#666'))

# (D) rate-sensor curves: Mk vs period for failing/sufficient/excess (MB, vary CDKI)
axD = fig.add_subplot(gs[1, 1])
for ez, col, lab in [(0.1, '#d9534f', 'ez 0.1 (failing)'), (0.35, '#e0a800', 'ez 0.35 (sufficient)'),
                     (1.0, '#117a65', 'ez 1.0 (baked, excess)')]:
    pers, mks = [], []
    for ck in (1, 2, 3, 4, 6):
        mk, _, per = measure(build(ez), 'MB', ck)
        if per == per: pers.append(per); mks.append(mk)
    axD.plot(pers, mks, 'o-', color=col, lw=1.8, ms=6, label=lab)
axD.set_xlabel('cycle period (h)  [slower → less dilution]'); axD.set_ylabel('mean H3K27me3 (Mk)')
axD.set_ylim(0, 1.02); axD.set_title('(D) Rate-sensor: mark rises with period only when not saturated')
axD.legend(fontsize=8.5, loc='lower right')

fig.suptitle('EZH2 generation vs replicative dilution of H3K27me3: failing / sufficient / excess — '
             'the baked model is in EXCESS (mark saturated, weak rate-sensor)', fontsize=11.5, fontweight='bold', y=0.995)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_v44_ezh2_dilution_regimes')
fig.savefig(out + '.png', dpi=140, bbox_inches='tight'); fig.savefig(out + '.pdf', bbox_inches='tight')
print('wrote', out + '.png/.pdf')
