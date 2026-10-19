"""CDKI mRNA bursting -> MB transient-G0 (JP 2026-10-18). Reads cdki_burst_{MB,GNP}.json from
sim_cdki_bursting_g0.py. Shows: (A) example MB cells drifting transiently into G0 as p21/p18 bursts spike
and recovering as they fall; (B) MB vs GNP cell-fate breakdown + snapshot G0 occupancy; (C) EZH2 protein
distribution, quiescent vs cycling (EZH2-low = the Ki67-low quiescent mode); (D) mechanism + the ~22% cap."""
import json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
SC = "/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/79019eb0-7c15-4650-b7d3-81347ccd6a13/scratchpad"
mb = json.load(open(f"{SC}/cdki_burst_MB.json")); gnp = json.load(open(f"{SC}/cdki_burst_GNP.json"))
DT = mb['dt']; THR = mb['prb_thr']

fig = plt.figure(figsize=(14, 9))
fig.suptitle("CDKI mRNA bursting (p21 + p18) drives MB transient-G0; proliferating GNPs stay cycling",
             fontsize=13.5, fontweight='bold')
gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.33, wspace=0.22)

# (A) example MB trajectories: pRb with CDKI bursts, G0 shaded
axA = fig.add_subplot(gs[0, :])
ex = mb['examples'][:3]
colors = ['#1f6f54', '#c0392b', '#2c6fa6']
for j, c in enumerate(ex):
    pr = np.array(c['pRb']); t = np.arange(len(pr)) * DT * 2 / 60.0   # hours (sub=::2)
    off = j * 0.0
    axA.plot(t, pr + off, color=colors[j], lw=1.6, label=f'cell {j+1}')
    g0 = pr < THR
    axA.fill_between(t, 0, pr, where=g0, color=colors[j], alpha=0.12, step='mid')
axA.axhline(THR, color='#888', ls='--', lw=1.2); axA.text(t[-1]*0.995, THR+0.05, 'committed-Rb threshold', ha='right', fontsize=8.5, color='#555')
axA.set_xlabel('time (h)', fontsize=9); axA.set_ylabel('pRb (commitment)', fontsize=9)
axA.set_title("(A) Example MB cells: p21/p18 bursts push pRb below threshold (shaded = transient G0), then cells re-enter as bursts fall",
              fontsize=10.5); axA.legend(fontsize=8.5, loc='upper right', ncol=3)

# (B) cell-fate breakdown MB vs GNP
axB = fig.add_subplot(gs[1, 0])
cats = ['pure\ncycling', 'transient\nG0', 'permanent\nG0']
mbv = [mb['cyc'], mb['transient'], mb['perm']]; gnv = [gnp['cyc'], gnp['transient'], gnp['perm']]
x = np.arange(3); w = 0.36
axB.bar(x - w/2, mbv, w, color='#c0392b', label=f"MB (n={mb['n']})")
axB.bar(x + w/2, gnv, w, color='#2c6fa6', label=f"GNP (n={gnp['n']})")
for xi, v in zip(x - w/2, mbv): axB.text(xi, v + 0.01, f'{v:.0%}', ha='center', fontsize=8, color='#c0392b')
for xi, v in zip(x + w/2, gnv): axB.text(xi, v + 0.01, f'{v:.0%}', ha='center', fontsize=8, color='#2c6fa6')
axB.set_xticks(x); axB.set_xticklabels(cats, fontsize=9); axB.set_ylim(0, 1.05); axB.set_ylabel('fraction of cells')
axB.legend(fontsize=9)
axB.set_title(f"(B) Cell fate over ~100h\nsnapshot G0: MB {mb['occ']:.0%} vs GNP {gnp['occ']:.0%}  ({mb['occ']/max(gnp['occ'],1e-6):.0f}x)",
              fontsize=10.5)

# (C) EZH2 distribution quiescent vs cycling (MB, pooled cell-time)
axC = fig.add_subplot(gs[1, 1])
ez = np.array(mb['ezh2_pool']); q = np.array(mb['q_pool'])
lez = np.log10(np.clip(ez, 1e-3, None))
bins = np.linspace(lez.min(), lez.max(), 40)
axC.hist(lez[~q], bins, density=True, color='#2c6fa6', alpha=0.55, label='cycling')
axC.hist(lez[q], bins, density=True, color='#95a5a6', alpha=0.7, label='quiescent (G0)')
axC.set_xlabel('log10 EZH2 protein', fontsize=9); axC.set_yticks([])
axC.legend(fontsize=9)
axC.set_title(f"(C) MB EZH2: quiescent below cycling\nmedian {mb['ezh2_q']:.2f} vs {mb['ezh2_c']:.2f} (fold {mb['ezh2_q']/mb['ezh2_c']:.2f})",
              fontsize=10.5)

# (D) mechanism text
axD = fig.add_axes([0.0, 0.0, 1.0, 0.045]); axD.axis('off')
cv = mb['cv']
txt = (f"Mechanism: INDEPENDENT transcriptional bursting of p21 (strongest driver) + p18 (temporal CV "
       f"p21 {cv['p21']:.1f}, p18 {cv['p18']:.1f}); MB's high CDKI + low mitogen drive sit near the committed-Rb "
       f"threshold so bursts tip ~{mb['occ']:.0%} of cells transiently into G0 (0% permanent), while high-SHH GNPs stay clear ({gnp['occ']:.0%}). "
       f"p27/Skp2 bursting do NOT cross (buffered / anti-arrest). Bursting caps ~25-30%; the last gap to IF 34% = deeper p57/OLIG2 dormancy.")
axD.text(0.5, 0.5, txt, ha='center', va='center', fontsize=8.6, wrap=True, color='#333')

out = f"{SC}/fig_cdki_bursting_g0.png"
plt.savefig(out, dpi=135, bbox_inches='tight'); print('saved', out)
