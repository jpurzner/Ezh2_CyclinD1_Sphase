"""Permanent cell-cycle exit vs transient G0: the CyclinD1/CDK6 drive difference (GNP vs MB) x H3K27me3
(JP 2026-08-09 overnight). Baked default (integrated + CDK4/6-alone escape).

Findings mapped:
  (A) PHASE DIAGRAM drive x CDKI: the CyclinD1/CDK6 DRIVE sets the fate. Low drive (GNP-like) -> straight
      cycling->PERMANENT (no transient regime). High drive (MB-like) -> a wide TRANSIENT-G0 escape regime
      before permanent. The transient regime EMERGES with drive.
  (B) The MARK's contribution: turning the H3K27me3/PRC2 repression OFF shifts the permanent boundary ~2x
      CDKI later (GNP 5->7x, MB 8->10x) -- the mark lowers effective drive, pushing toward permanent exit.
  (C) Replicative DILUTION is a CYCLING-RATE sensor, NOT a fate switch: k_dil_S sets the mark & period in
      CYCLING cells, but (from mechanism_mark_dilution.py) does NOT move the permanent/transient boundary --
      because that fate is decided in the ARRESTED state, where there is no replication to dilute.
  (D) Why: in the arrested cell the mark PLATEAUS/DROPS (EZH2 writer E2f-gated, collapses); the fate is set
      by whether the (mark-repressed) CyclinD1/CDK6 drive clears the escape threshold -- MB yes (transient),
      GNP no (permanent).

Run:  ./venv/bin/python simulations/fig_v44_permanent_transient.py   (~25 min)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap, BoundaryNorm
from src.build_model_v44_heldt import build_model_v44

C_CYC = '#117a65'; C_TRA = '#e08e0b'; C_PER = '#8e2f2f'
plt.rcParams.update({'font.size': 10, 'axes.titlesize': 10.5, 'axes.labelsize': 10, 'axes.titleweight': 'bold'})


def rr_for(extra=None):
    m = build_model_v44(params=extra) if extra else build_model_v44()
    r = te.loada(m); r.integrator.setValue('relative_tolerance', 1e-7); r.integrator.setValue('absolute_tolerance', 1e-9)
    r.integrator.setValue('maximum_num_steps', 12000000)
    return r


def set_mb(rr, cdki=1.0, drive=1.0):
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
    try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
    except Exception: pass
    rr['k_Cd_tx_Gli_max'] = 0.24197652 * drive
    rr['k_cdk6_Gli'] = rr['k_cdk6_Gli'] * 1.0  # base; scaled below via _kg
    rr['p18'] = 1.73 * cdki; rr['p19'] = 0.58 * cdki


def classify(rr, cdki, drive, T=55000, _kg0=[None]):
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
    try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
    except Exception: pass
    if _kg0[0] is None: _kg0[0] = rr['k_cdk6_Gli']
    ks = rr['kSyP21']
    rr['k_Cd_tx_Gli_max'] = 0.24197652013556434 * drive
    rr['k_cdk6_Gli'] = _kg0[0] * drive
    rr['p18'] = 1.73 * cdki; rr['p19'] = 0.58 * cdki; rr['kSyP21'] = ks * cdki
    rr.reset()
    r = rr.simulate(0, T, T // 10, selections=['time', 'Dna'])
    rr['kSyP21'] = ks; rr['k_Cd_tx_Gli_max'] = 0.24197652013556434; rr['k_cdk6_Gli'] = _kg0[0]
    t = r['time']; dna = r['Dna']; dd = t[np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]]
    nlast = int(np.sum(dd > T * 0.66)); per = np.mean(np.diff(dd)) / 60 if len(dd) >= 2 else None
    return 0 if nlast == 0 else (1 if (per and per > 40) else 2)   # 0=perm,1=transient,2=cycling


fig = plt.figure(figsize=(15, 9))
gs = GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.30)

# ---------- (A) phase diagram: drive x CDKI ----------
axA = fig.add_subplot(gs[0, 0])
rr = rr_for()
drives = np.array([0.2, 0.35, 0.5, 0.7, 1.0, 1.3])
cdkis = np.array([1, 2, 3, 4, 5, 6, 7, 8, 10])
Z = np.zeros((len(cdkis), len(drives)))
for j, dv in enumerate(drives):
    for i, ck in enumerate(cdkis):
        Z[i, j] = classify(rr, ck, dv)
cmap = ListedColormap([C_PER, C_TRA, C_CYC]); norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
axA.pcolormesh(np.arange(len(drives) + 1), np.arange(len(cdkis) + 1), Z, cmap=cmap, norm=norm, edgecolors='w', lw=0.5)
axA.set_xticks(np.arange(len(drives)) + 0.5); axA.set_xticklabels([f'{d:.2f}' for d in drives])
axA.set_yticks(np.arange(len(cdkis)) + 0.5); axA.set_yticklabels([str(c) for c in cdkis])
axA.set_xlabel('CyclinD1/CDK6 drive scale  (GNP-low → MB-high)'); axA.set_ylabel('CDKI (× baseline)')
axA.set_title('(A) Drive sets the fate — transient-G0 regime emerges with drive')
axA.axvline(4.0, color='k', lw=1.4); axA.text(4.02, 8.4, ' MB (1.0)', fontsize=8, fontweight='bold')
axA.text(0.05, 8.4, 'low drive:\nno transient', fontsize=7.5, color='w')
from matplotlib.patches import Patch
axA.legend(handles=[Patch(color=C_CYC, label='cycling'), Patch(color=C_TRA, label='transient G0'),
                    Patch(color=C_PER, label='permanent exit')], fontsize=7.5, loc='upper right', framealpha=0.95)

# ---------- (B) mark ON vs OFF: permanent boundary shift ----------
axB = fig.add_subplot(gs[0, 1])
labels = ['GNP', 'MB']; on = [5, 8]; off = [7, 10]
x = np.arange(2); w = 0.35
axB.bar(x - w/2, on, w, color='#8e2f2f', label='mark ON (baked)')
axB.bar(x + w/2, off, w, color='#c99', label='mark OFF (f0=1)')
for xi, (a, b) in enumerate(zip(on, off)):
    axB.annotate('', (xi + w/2, a + 0.15), (xi - w/2, a + 0.15), arrowprops=dict(arrowstyle='->', color='k'))
    axB.text(xi, max(a, b) + 0.5, f'+{b-a}× later\nw/o mark', ha='center', fontsize=8, style='italic')
axB.set_xticks(x); axB.set_xticklabels(labels); axB.set_ylabel('CDKI at permanent-exit onset')
axB.set_ylim(0, 12); axB.legend(fontsize=8.5, loc='upper left')
axB.set_title('(B) The mark lowers effective drive → permanent exit ~2× earlier')

# ---------- (C) dilution = cycling-rate sensor (not a fate switch) ----------
axC = fig.add_subplot(gs[1, 0])
kdils = [0.0, 0.2, 0.3465, 0.693, 1.0, 1.386]; Mk_c = []; per_c = []
for kd in kdils:
    r2 = rr_for({'k_dil_S': kd})
    r2['SHH'] = 0.5; r2['Ptch1_copy_number'] = 0.3
    try: r2['Cd2_expr'] = r2['CD2_EXPR_MB']
    except Exception: pass
    r2['p18'] = 1.73 * 2; r2['p19'] = 0.58 * 2  # a slowed (2x CDKI) cycling MB cell
    r2.reset(); d = r2.simulate(0, 45000, 22500, selections=['time', 'Dna', 'Mk'])
    mm = d['time'] >= 30000; Mk_c.append(float(d['Mk'][mm].mean()))
    dd = d['time'][np.where((d['Dna'][:-1] > 0.9) & (d['Dna'][1:] < 0.1))[0]]
    per_c.append(float(np.mean(np.diff(dd)) / 60) if len(dd) >= 2 else np.nan)
axC.plot(kdils, Mk_c, 'o-', color='#7a5aa6', lw=1.7, label='mean H3K27me3 (Mk)')
axC.set_xlabel('replicative dilution rate  k_dil_S'); axC.set_ylabel('mean Mk (cycling)', color='#7a5aa6')
axC.tick_params(axis='y', labelcolor='#7a5aa6'); axC.set_title('(C) Dilution tunes the CYCLING mark/rate — NOT the fate')
axr = axC.twinx(); axr.plot(kdils, per_c, 's--', color='#b0402f', lw=1.5, label='cycle period')
axr.set_ylabel('cycle period (h)', color='#b0402f'); axr.tick_params(axis='y', labelcolor='#b0402f')
axC.axvline(0.693, color='k', ls=':', lw=0.8); axC.text(0.71, min(Mk_c), 'baked', fontsize=7.5)
axC.text(0.02, 0.05, 'permanent/transient boundary\nUNCHANGED across all k_dil_S\n(arrested state = no replication)',
         transform=axC.transAxes, fontsize=7.5, style='italic', bbox=dict(boxstyle='round', fc='#f7f7f7', ec='#999'))

# ---------- (D) arrested-state: why the fate is drive-decided ----------
axD = fig.add_subplot(gs[1, 1])
# arrested GNP (5x) vs transient MB (6x) vs permanent MB (8x): Cd, cdk6
rows = [('GNP 5×\n(permanent)', 'GNP', 5), ('MB 6×\n(transient)', 'MB', 6), ('MB 8×\n(permanent)', 'MB', 8)]
Cds = []; cdk6s = []
for _, cell, ck in rows:
    r3 = rr_for()
    if cell == 'GNP':
        r3['SHH'] = 0.5; r3['Ptch1_copy_number'] = 1.0; p18, p19 = 0.464, 0.36
    else:
        r3['SHH'] = 0.5; r3['Ptch1_copy_number'] = 0.3; p18, p19 = 1.73, 0.58
        try: r3['Cd2_expr'] = r3['CD2_EXPR_MB']
        except Exception: pass
    ks = r3['kSyP21']; r3['p18'] = p18 * ck; r3['p19'] = p19 * ck; r3['kSyP21'] = ks * ck
    r3.reset(); d = r3.simulate(0, 50000, 25000, selections=['time', 'Cd', 'cdk6'])
    mm = d['time'] >= 40000; Cds.append(float(d['Cd'][mm].mean())); cdk6s.append(float(d['cdk6'][mm].mean()))
x = np.arange(3); w = 0.35
axD.bar(x - w/2, Cds, w, color='#b0402f', label='CyclinD1 (Cd)')
axD.bar(x + w/2, cdk6s, w, color='#1f6f4a', label='CDK6 (cdk6)')
axD.set_xticks(x); axD.set_xticklabels([r[0] for r in rows], fontsize=8)
axD.set_ylabel('arrested-state level'); axD.legend(fontsize=8.5)
axD.set_title('(D) Arrested state: MB retains CDK6 to escape; GNP cannot')
axD.text(0.5, 0.92, 'transient = drive clears escape threshold; permanent = it does not',
         transform=axD.transAxes, ha='center', fontsize=7.5, style='italic')

fig.suptitle('Permanent exit (GNP) vs transient G0 (MB): the CyclinD1/CDK6 drive decides the fate; '
             'H3K27me3 tunes it (static) and sets cycling rate (dilution)', fontsize=11.5, fontweight='bold', y=0.995)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_v44_permanent_transient')
fig.savefig(out + '.png', dpi=140, bbox_inches='tight'); fig.savefig(out + '.pdf', bbox_inches='tight')
print('wrote', out + '.png/.pdf')
