"""Soluble EZH2 is a CONCENTRATION, not an amount. At a symmetric division a daughter inherits half the
protein in half the volume -> the CONCENTRATION is preserved (with_ezh2_conc=True, DEFAULT). The legacy
convention (with_ezh2_conc=False) discretely HALVES EZH2 at every mitosis (EZH2 = EZH2/2), treating it as an
amount while every downstream Hill reads it as a concentration -- an abrupt ~2x sawtooth drop with no physical
basis. Because kTlEZ is re-fit (x0.78 -> 0.011771) so the two conventions sit at the SAME mean level, the swap
is validation-NEUTRAL: it only replaces the discontinuous jump with a smooth, turnover-limited (~10 h) decay of
the writer carried from mother to daughter.

Panels:
  (A) EZH2(t) across several divisions, legacy halving (sawtooth) vs concentration-preserved default (smooth),
      overlaid; divisions marked with vertical lines.
  (B) Zoom on ONE division: the discontinuous /2 jump vs the smooth mother->daughter EZH2 carryover (traces
      aligned to each convention's own division at relative t = 0).
  (C) Validation-neutrality: EZH2 levels / phase ratios (GNP level, MB level, MB/GNP fold, within-cycle G2:G1)
      are preserved across the two conventions.

Species: EZH2 (protein), EZH2m (mRNA). The halving lives in the E_div division event.

Run:  ./venv/bin/python simulations/fig_v44_ezh2_concentration.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt, tellurium as te
from src.build_model_v44_heldt import build_model_v44

SEL = ["time", "EZH2", "EZH2m", "Cd", "P21", "aRc", "Dna", "MPF", "mass"]
T_END, N_PTS = 15000, 60000          # dt = 0.25 min -> resolves the discontinuous halving

# ---- two model builds: concentration-preserved (default) vs legacy amount-halving ----
M_CONC = build_model_v44(with_ezh2_conc=True)    # DEFAULT: no halving, kTlEZ re-fit x0.78
M_LEG  = build_model_v44(with_ezh2_conc=False)   # legacy: EZH2 = EZH2/2 at each mitosis

# ---- conditions ----
GNP = dict(shh=0.5, ptch1=1.0, mycn=1.0, p16=0.0, p18=None, ksyp21=None, p21div=0.6)
MB  = dict(shh=0.5, ptch1=0.1, mycn=2.8, p16=0.306, p18=1.553, ksyp21=0.002, p21div=1.8)


def run(model, cond):
    """Fresh loada per condition (reset() does NOT reset parameters)."""
    rr = te.loada(model)
    rr.integrator.setValue("absolute_tolerance", 1e-8)
    rr.integrator.setValue("relative_tolerance", 1e-6)
    rr.integrator.setValue("maximum_time_step", 5.0)
    rr.integrator.setValue("maximum_num_steps", 2000000)
    rr['SHH'] = cond['shh']; rr['Ptch1_copy_number'] = cond['ptch1']
    rr['MYCN_amplification'] = cond['mycn']; rr['p16'] = cond['p16']
    if cond['p18'] is not None:    rr['p18'] = cond['p18']
    if cond['ksyp21'] is not None: rr['kSyP21'] = cond['ksyp21']
    rr['P21_div'] = cond['p21div']
    return rr.simulate(0, T_END, N_PTS, selections=SEL)


def division_times(res):
    """Divisions via the Dna 1->0 reset (robust to coarse sampling; MPF peaks are too sharp)."""
    t = res['time']; dna = res['Dna']
    idx = np.where((dna[:-1] > 0.9) & (dna[1:] < 0.5))[0]
    return t[idx]


def phase_masks(res, settle=4000):
    """Cell-cycle phase masks on the settled window."""
    t = res['time']; P21 = res['P21']; aRc = res['aRc']; Dna = res['Dna']; MPF = res['MPF']
    s = t >= settle
    in_S  = (aRc > 0.05) & (Dna < 0.98)
    in_G2 = Dna >= 0.98
    in_M  = MPF > 0.3
    preS  = (~in_S) & (~in_G2) & (~in_M)
    G0 = preS & (P21 > 0.1)
    G1 = preS & (P21 <= 0.1)
    return dict(s=s, G0=G0 & s, G1=G1 & s, S=in_S & s, G2=in_G2 & s, preS=preS & s)


def ezh2_stats(res):
    """Mean settled EZH2 protein + mRNA + within-cycle G2:G1 ratio."""
    m = phase_masks(res); EZ = res['EZH2']
    lvl = float(np.mean(EZ[m['s']]))
    mrna = float(np.mean(res['EZH2m'][m['s']]))
    ez_g2 = float(np.mean(EZ[m['G2']])) if m['G2'].sum() else np.nan
    # denominator: interphase pre-S (G1 if present, else G0) -- robust across GNP/MB
    denom_mask = m['G1'] if m['G1'].sum() > 20 else m['preS']
    ez_g1 = float(np.mean(EZ[denom_mask])) if denom_mask.sum() else np.nan
    return dict(lvl=lvl, mrna=mrna, g2g1=(ez_g2 / ez_g1 if ez_g1 and not np.isnan(ez_g1) else np.nan))


# ============================ run the four cases ============================
print("=" * 78)
print("EZH2 CONCENTRATION vs legacy AMOUNT-HALVING at division")
print("=" * 78)
res = {}
for tag, model in [('conc', M_CONC), ('leg', M_LEG)]:
    for cond_name, cond in [('GNP', GNP), ('MB', MB)]:
        res[(tag, cond_name)] = run(model, cond)

stats = {}
for k, r in res.items():
    stats[k] = ezh2_stats(r)
    ndiv = len(division_times(r))
    print(f"  {k[0]:4s} {k[1]:3s}:  EZH2 level={stats[k]['lvl']:.3f}   mRNA={stats[k]['mrna']:.3f}   "
          f"G2:G1={stats[k]['g2g1']:.2f}   divisions={ndiv}")

# MB/GNP folds
fold_conc = stats[('conc', 'MB')]['lvl'] / stats[('conc', 'GNP')]['lvl']
fold_leg  = stats[('leg',  'MB')]['lvl'] / stats[('leg',  'GNP')]['lvl']
print(f"\n  MB/GNP EZH2 fold:   conc={fold_conc:.3f}   legacy={fold_leg:.3f}   (Delta={abs(fold_conc-fold_leg):.3f})")

# ============================ figure ============================
fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2))

# ---- (A) EZH2(t) over several divisions, GNP, both conventions ----
axA = axes[0]
rC = res[('conc', 'GNP')]; rL = res[('leg', 'GNP')]
tC = rC['time'] / 60.0; tL = rL['time'] / 60.0
divC = division_times(rC) / 60.0
# show a window with several divisions after settling
w0, w1 = 100, 200   # hours
for d in divC[(divC >= w0) & (divC <= w1)]:
    axA.axvline(d, color='0.55', ls=':', lw=0.9, zorder=0)
axA.plot(tL, rL['EZH2'], color='#c0392b', lw=1.5, label='legacy amount-halving (÷2 sawtooth)')
axA.plot(tC, rC['EZH2'], color='#1f6f8b', lw=1.7, label='concentration-preserved (default)')
axA.set_xlim(w0, w1)
axA.set_xlabel('time (h)'); axA.set_ylabel('EZH2 protein (a.u.)')
axA.set_title('(a) EZH2 across divisions (GNP): sawtooth vs smooth', fontweight='bold', fontsize=11)
axA.legend(fontsize=8, loc='upper right'); axA.grid(alpha=0.2)
axA.text(0.015, 0.03, 'dotted = divisions', transform=axA.transAxes, fontsize=7.5, color='0.4')

# ---- (B) zoom on ONE division, aligned to each convention's own division ----
axB = axes[1]
def aligned_trace(r, target_h=140.0, half=3.5):
    t = r['time'] / 60.0
    dts = division_times(r) / 60.0
    d = dts[np.argmin(np.abs(dts - target_h))]
    m = (t >= d - half) & (t <= d + half)
    return t[m] - d, r['EZH2'][m], d
relC, ezC, dC = aligned_trace(rC)
relL, ezL, dL = aligned_trace(rL)
axB.axvline(0, color='k', ls='--', lw=1.0, alpha=0.6)
axB.plot(relL, ezL, color='#c0392b', lw=2.0, label=f'legacy: discontinuous ÷2 jump')
axB.plot(relC, ezC, color='#1f6f8b', lw=2.0, label='concentration: smooth carryover')
# annotate the halving jump
jump_hi = ezL[np.argmin(np.abs(relL + 0.05))]
axB.annotate('÷2 at mitosis\n(no physical basis)', xy=(0, jump_hi * 0.72), xytext=(0.6, jump_hi * 0.5),
             fontsize=8, color='#c0392b',
             arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.2))
_yb = axB.get_ylim()[0]
axB.text(relL.min() + 0.2, _yb, 'mother', fontsize=8, color='0.35', va='bottom')
axB.text(relL.max() - 0.2, _yb, 'daughter', fontsize=8, color='0.35', va='bottom', ha='right')
axB.set_xlabel('time relative to division (h)'); axB.set_ylabel('EZH2 protein (a.u.)')
axB.set_title('(b) One division: halving jump vs smooth ~10 h decay', fontweight='bold', fontsize=11)
axB.legend(fontsize=8, loc='upper right'); axB.grid(alpha=0.2)

# ---- (C) validation-neutrality bars ----
axC = axes[2]
metrics = ['GNP EZH2\nprotein', 'MB EZH2\nprotein', 'MB/GNP\nfold', 'GNP EZH2\nmRNA']
vals_conc = [stats[('conc', 'GNP')]['lvl'], stats[('conc', 'MB')]['lvl'], fold_conc,
             stats[('conc', 'GNP')]['mrna']]
vals_leg  = [stats[('leg',  'GNP')]['lvl'], stats[('leg',  'MB')]['lvl'], fold_leg,
             stats[('leg',  'GNP')]['mrna']]
x = np.arange(len(metrics)); bw = 0.38
axC.bar(x - bw/2, vals_conc, bw, color='#1f6f8b', label='concentration (default)')
axC.bar(x + bw/2, vals_leg,  bw, color='#c0392b', label='legacy halving')
for xi, (vc, vl) in enumerate(zip(vals_conc, vals_leg)):
    axC.text(xi - bw/2, vc, f'{vc:.2f}', ha='center', va='bottom', fontsize=7.5)
    axC.text(xi + bw/2, vl, f'{vl:.2f}', ha='center', va='bottom', fontsize=7.5)
axC.set_xticks(x); axC.set_xticklabels(metrics, fontsize=8.5)
axC.set_ylabel('EZH2 level / ratio')
axC.set_title('(c) Validation-neutral: levels & phase ratios preserved', fontweight='bold', fontsize=11)
axC.legend(fontsize=8, loc='upper right'); axC.grid(alpha=0.2, axis='y')
axC.set_ylim(0, max(vals_conc + vals_leg) * 1.22)

fig.suptitle('Soluble EZH2 is a concentration: preserved at symmetric division, reset by ~10 h turnover — not amount-halved',
             fontsize=13, fontweight='bold', y=1.005)
fig.text(0.5, 0.945,
         'The legacy ÷2 sawtooth is replaced by a smooth, turnover-limited writer carryover from mother to daughter; a compensating kTlEZ re-fit keeps every EZH2 level/ratio unchanged (validation-neutral).',
         ha='center', fontsize=9.5, style='italic', color='0.3')

plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig('simulations/fig_v44_ezh2_concentration.png', dpi=170)
plt.savefig('simulations/fig_v44_ezh2_concentration.pdf')
plt.close()
print("\nSaved: simulations/fig_v44_ezh2_concentration.png / .pdf")
print(f"  panel B aligned divisions:  conc @ {dC:.1f} h,  legacy @ {dL:.1f} h")
