"""Figures from the overnight run (JP 2026-10-18). Reads simulations/overnight/results.json.
fig_overnight_lock.png     — locked ~34% MB reference (two-step, model-tuned protein) vs GNP.
fig_overnight_attribution.png — 2x2: {model-tuned, measured-protein} x {burst p21+p18, p21-only}."""
import json, os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, 'overnight')
R = json.load(open(os.path.join(OUT, 'results.json')))
mb = R['phaseA_lock_MB']; gnp = R['phaseA_lock_GNP']; at = R['phaseB_attribution']
DT = 30.0; THR = 1.5

# ---------- Figure 1: locked reference ----------
fig = plt.figure(figsize=(14, 8.5))
fig.suptitle(f"CDKI mRNA bursting → MB transient-G0: snapshot-G0 {mb['occ']:.0%} (IF ~34%), {mb['transient']:.0%} of cells transiently arrest; GNPs stay cycling ({gnp['occ']:.0%}, N={mb['N']})",
             fontsize=12, fontweight='bold')
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.33, wspace=0.22)

axA = fig.add_subplot(gs[0, :])
c = mb['examples'][0]; pr = np.array(c['pRb']); t = np.arange(len(pr)) * DT * 2 / 60.0
axA.plot(t, pr, color='#1f2d3d', lw=1.8, label='pRb (commitment)', zorder=3)
axA.fill_between(t, 0, pr, where=pr < THR, color='#95a5a6', alpha=0.25, step='mid', label='transient G0')
axA.axhline(THR, color='#888', ls='--', lw=1.1); axA.set_ylabel('pRb'); axA.set_xlabel('time (h)')
ax2 = axA.twinx()
ax2.plot(t, np.array(c['p21a_m']), color='#c0392b', lw=1.3, label='p21 mRNA (t½ 5.2h)')
ax2.plot(t, np.array(c['p21a']), color='#e8a33d', lw=1.3, label='p21 protein (t½ 0.6h)')
ax2.set_ylabel('p21 mRNA / protein (a.u.)')
l1, b1 = axA.get_legend_handles_labels(); l2, b2 = ax2.get_legend_handles_labels()
axA.legend(l1 + l2, b1 + b2, fontsize=8.5, loc='upper right', ncol=4)
axA.set_title("(A) One MB cell: p21 transcription bursts → mRNA → protein → pRb dips into G0, then recovers (two-stage filter)", fontsize=10.5)

axB = fig.add_subplot(gs[1, 0])
cats = ['pure\ncycling', 'transient\nG0', 'permanent\nG0']; x = np.arange(3); w = 0.36
axB.bar(x - w/2, [mb['cyc'], mb['transient'], mb['perm']], w, color='#c0392b', label=f"MB (n={mb['n']}/{mb['N']})")
axB.bar(x + w/2, [gnp['cyc'], gnp['transient'], gnp['perm']], w, color='#2c6fa6', label=f"GNP (n={gnp['n']}/{gnp['N']})")
for xi, v in zip(x - w/2, [mb['cyc'], mb['transient'], mb['perm']]): axB.text(xi, v + 0.01, f'{v:.0%}', ha='center', fontsize=8, color='#c0392b')
for xi, v in zip(x + w/2, [gnp['cyc'], gnp['transient'], gnp['perm']]): axB.text(xi, v + 0.01, f'{v:.0%}', ha='center', fontsize=8, color='#2c6fa6')
axB.set_xticks(x); axB.set_xticklabels(cats, fontsize=9); axB.set_ylim(0, 1.08); axB.set_ylabel('fraction of cells'); axB.legend(fontsize=8.5)
axB.set_title(f"(B) Cell fate over ~115h\nsnapshot-G0 occupancy: MB {mb['occ']:.0%} (IF ~34%) vs GNP {gnp['occ']:.0%}", fontsize=10.5)

axC = fig.add_subplot(gs[1, 1])
labels = ['p21\nmRNA', 'p21\nprotein', 'p18\nprotein', 'p27\n(P21)']
vals = [mb['cv_p21a_m'], mb['cv_p21a'], mb['cv_p18_prot'], mb['cv_P21']]
axC.bar(range(4), vals, color=['#c0392b', '#e8a33d', '#2c6fa6', '#7f8c8d'])
for i, v in enumerate(vals): axC.text(i, v + 0.02, f'{v:.2f}', ha='center', fontsize=8.5)
axC.set_xticks(range(4)); axC.set_xticklabels(labels, fontsize=8.5); axC.set_ylabel('temporal CV (within-cell)')
axC.set_title(f"(C) Variance (MB): p21 carries it\nEZH2 quiescent {mb['ezh2_q']:.2f} vs cycling {mb['ezh2_c']:.2f} (fold {mb['ezh2_q']/mb['ezh2_c']:.2f})", fontsize=10.5)
fig.text(0.5, 0.012,
         f"Snapshot-G0 {mb['occ']:.0%} is a FLOOR: integrator attrition ({mb['mortality']:.0%} of MB cells) is arrest-biased, so permanent-G0 ({mb['perm']:.0%}) is under-counted; "
         f"a data-grounded measured-protein variant gives ~50–57% (both pass 30/32) → report MB snapshot-G0 as ~36–55% across INK4/CIP protein-stability. "
         f"GNP 0% reflects GNPs' ~10× lower p21 (measured CDKN1A MB>GNP), not the SHH drive.",
         ha='center', va='bottom', fontsize=7.6, color='#555', wrap=True)
plt.tight_layout(rect=[0, 0.045, 1, 0.95])
f1 = os.path.join(OUT, 'fig_overnight_lock.png'); plt.savefig(f1, dpi=135); print('saved', f1)

# ---------- Figure 2: 2x2 attribution ----------
fig2, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5.2))
fig2.suptitle("Measured-protein attribution: p21 carries the transient-G0; stable p18 (24h) filters its own bursts out",
              fontsize=12, fontweight='bold')
groups = ['model-tuned\nprotein', 'measured\nprotein (p18 24h)']
occ_full = [at['MB_tuned_p21+p18']['occ'], at['MB_measured_p21+p18']['occ']]
occ_p21 = [at['MB_tuned_p21only']['occ'], at['MB_measured_p21only']['occ']]
perm_full = [at['MB_tuned_p21+p18']['perm'], at['MB_measured_p21+p18']['perm']]
perm_p21 = [at['MB_tuned_p21only']['perm'], at['MB_measured_p21only']['perm']]
x = np.arange(2); w = 0.36
b1 = axL.bar(x - w/2, occ_full, w, color='#c0392b', label='burst p21 + p18')
b2 = axL.bar(x + w/2, occ_p21, w, color='#e8a33d', label='burst p21 only')
for bars, ocs, pms in [(b1, occ_full, perm_full), (b2, occ_p21, perm_p21)]:
    for bar, o, p in zip(bars, ocs, pms):
        axL.text(bar.get_x() + bar.get_width()/2, o + 0.01, f'{o:.0%}\n({p:.0%} perm)', ha='center', fontsize=8)
axL.axhline(0.34, color='#2c6fa6', ls='--', lw=1.2); axL.text(1.4, 0.345, 'IF ~34%', color='#2c6fa6', fontsize=8.5)
axL.set_xticks(x); axL.set_xticklabels(groups, fontsize=9.5); axL.set_ylim(0, 0.72); axL.set_ylabel('MB snapshot-G0'); axL.legend(fontsize=9)
axL.set_title("(A) p21-only ≈ p21+p18 in both → p21 is the driver\n(measured overshoots: stable p21 0.98h dwells longer)", fontsize=10)

# p18 variance collapse
p18cv_full = [at['MB_tuned_p21+p18']['cv_p18_prot'], at['MB_measured_p21+p18']['cv_p18_prot']]
axR.bar(x, p18cv_full, 0.5, color='#2c6fa6')
for xi, v in zip(x, p18cv_full): axR.text(xi, v + 0.005, f'{v:.2f}', ha='center', fontsize=9)
axR.set_xticks(x); axR.set_xticklabels(groups, fontsize=9.5); axR.set_ylabel('p18 protein temporal CV (when p18 is bursting)')
axR.set_title("(B) Stable p18 filters its own bursts:\nCV 0.34 → 0.16 when protein t½ = measured 24h", fontsize=10)
plt.tight_layout(rect=[0, 0, 1, 0.92])
f2 = os.path.join(OUT, 'fig_overnight_attribution.png'); plt.savefig(f2, dpi=135); print('saved', f2)
