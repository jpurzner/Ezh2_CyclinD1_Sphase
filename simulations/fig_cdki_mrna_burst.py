"""Two-step CDKI mRNA bursting -> MB transient-G0 (JP 2026-10-18). Reads cdki_mrna_burst_{MB,GNP}.json.
(A) one MB cell: a p21 transcription burst raises p21a mRNA, the (short) protein tracks it, and pRb dips
below the committed threshold (transient G0), then recovers -- the two-stage filter in action.
(B) MB vs GNP fate + snapshot G0 (GNP = genuinely low CDKI + high SHH -> stays cycling).
(C) two-stage temporal CV (mRNA -> protein per CDKI) + EZH2 quiescent vs cycling."""
import json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
SC = "/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/79019eb0-7c15-4650-b7d3-81347ccd6a13/scratchpad"
mb = json.load(open(f"{SC}/cdki_mrna_burst_MB.json")); gnp = json.load(open(f"{SC}/cdki_mrna_burst_GNP.json"))
DT = mb['dt']; THR = mb['prb_thr']

fig = plt.figure(figsize=(14, 9))
fig.suptitle("Two-step CDKI expression (transcription -> mRNA -> protein): p21 mRNA bursts drive MB transient-G0; GNPs stay cycling",
             fontsize=13, fontweight='bold')
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.34, wspace=0.24)

# (A) one MB cell: pRb + p21a mRNA + p21a protein (two-stage filter)
axA = fig.add_subplot(gs[0, :])
c = mb['examples'][0]
pr = np.array(c['pRb']); t = np.arange(len(pr)) * DT * 2 / 60.0
axA.plot(t, pr, color='#1f2d3d', lw=1.8, label='pRb (commitment)', zorder=3)
axA.fill_between(t, 0, pr, where=pr < THR, color='#95a5a6', alpha=0.25, step='mid', label='transient G0')
axA.axhline(THR, color='#888', ls='--', lw=1.1); axA.set_ylabel('pRb', fontsize=9); axA.set_xlabel('time (h)', fontsize=9)
ax2 = axA.twinx()
ax2.plot(t, np.array(c['p21a_m']), color='#c0392b', lw=1.3, alpha=0.9, label='p21 mRNA (t½ 5.2h)')
ax2.plot(t, np.array(c['p21a']), color='#e8a33d', lw=1.3, alpha=0.9, label='p21 protein (t½ 0.6h)')
ax2.set_ylabel('p21 mRNA / protein (a.u.)', fontsize=9)
l1, lab1 = axA.get_legend_handles_labels(); l2, lab2 = ax2.get_legend_handles_labels()
axA.legend(l1 + l2, lab1 + lab2, fontsize=8.5, loc='upper right', ncol=4)
axA.set_title("(A) One MB cell: a transcription burst raises p21 mRNA, the short-lived protein tracks it, pRb dips into G0, then recovers",
              fontsize=10.5)

# (B) fate breakdown MB vs GNP
axB = fig.add_subplot(gs[1, 0])
cats = ['pure\ncycling', 'transient\nG0', 'permanent\nG0']
mbv = [mb['cyc'], mb['transient'], mb['perm']]; gnv = [gnp['cyc'], gnp['transient'], gnp['perm']]
x = np.arange(3); w = 0.36
axB.bar(x - w/2, mbv, w, color='#c0392b', label=f"MB (n={mb['n']}/{mb['N']})")
axB.bar(x + w/2, gnv, w, color='#2c6fa6', label=f"GNP (n={gnp['n']}/{gnp['N']})")
for xi, v in zip(x - w/2, mbv): axB.text(xi, v + 0.01, f'{v:.0%}', ha='center', fontsize=8, color='#c0392b')
for xi, v in zip(x + w/2, gnv): axB.text(xi, v + 0.01, f'{v:.0%}', ha='center', fontsize=8, color='#2c6fa6')
axB.set_xticks(x); axB.set_xticklabels(cats, fontsize=9); axB.set_ylim(0, 1.08); axB.set_ylabel('fraction of cells')
axB.legend(fontsize=8.5)
axB.set_title(f"(B) Cell fate over ~100h\nsnapshot G0: MB {mb['occ']:.0%} (IF target 34%) vs GNP {gnp['occ']:.0%}", fontsize=10.5)

# (C) two-stage CV + EZH2
axC = fig.add_subplot(gs[1, 1])
cv = mb['cv']
labels = ['p21', 'p18', 'p27']
mrna_cv = [cv['p21a_m'], np.nan, np.nan]  # only p21 mRNA recorded in detail
prot_cv = [cv['p21a'], cv['p18_prot'], cv['p27']]
xc = np.arange(3); w = 0.36
axC.bar(xc - w/2, mrna_cv, w, color='#c0392b', alpha=0.85, label='mRNA CV')
axC.bar(xc + w/2, prot_cv, w, color='#e8a33d', alpha=0.9, label='protein CV')
axC.set_xticks(xc); axC.set_xticklabels(labels, fontsize=9); axC.set_ylabel('temporal CV (within-cell)')
axC.legend(fontsize=8.5, loc='upper right')
axC.set_title(f"(C) Variance per CDKI (MB)\nEZH2 quiescent {mb['ezh2_q']:.2f} vs cycling {mb['ezh2_c']:.2f} (fold {mb['ezh2_q']/mb['ezh2_c']:.2f})",
              fontsize=10.5)

axD = fig.add_axes([0.0, 0.0, 1.0, 0.05]); axD.axis('off')
txt = (f"Burst lands on TRANSCRIPTION (p21a_tx/p18_tx) -> mRNA (measured t½) -> protein (measured t½): two low-passes, "
       f"no growth leak. MB (high CDKI, SHH 0.5) sits near the committed-Rb threshold so p21 mRNA bursts tip ~{mb['occ']:.0%} of cells "
       f"transiently into G0 ({mb['perm']:.0%} permanent); proliferating GNPs (low p21, SHH 1.0) stay cycling ({gnp['occ']:.0%} G0, n={gnp['n']}/{gnp['N']}). "
       f"The mRNA low-pass also removes the parameter-step stiffness (division-free integration, GNP 0% mortality).")
axD.text(0.5, 0.5, txt, ha='center', va='center', fontsize=8.5, wrap=True, color='#333')

out = f"{SC}/fig_cdki_mrna_burst.png"
plt.savefig(out, dpi=135, bbox_inches='tight'); print('saved', out)
