"""EZH2 variance decomposition: cell-to-cell vs within-cell-over-time, and why the split is a guess.
Reads the 4 population JSONs from sim_ezh2_variance_decomp.py. Three panels:
 A  many-cell EZH2 over time (baked) -> emergent slow protein drift riding on cell-to-cell offsets
 B  baked variance split 3 ways: cell-to-cell (frozen) / cell-cycle (deterministic) / bursting (stochastic);
    the cell-cycle part is isolated by the no-bursting control
 C  THE GUESS: baked vs mostly-frozen reach the SAME snapshot CV with opposite temporal fractions -> the IF
    snapshot (which sees only the marginal) cannot resolve the split.
"""
import json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
SC = "/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/79019eb0-7c15-4650-b7d3-81347ccd6a13/scratchpad"
D = {k: json.load(open(f"{SC}/vd_{k}.json")) for k in ['baked', 'noburst', 'frozen', 'bursty']}

fig = plt.figure(figsize=(16, 5.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1.5, 1, 1.1], wspace=0.30)
fig.suptitle("EZH2 variance: cell-to-cell vs within-cell-over-time — the split is a modeling guess the IF snapshot cannot resolve",
             fontsize=12.5, fontweight='bold')

# ---- A: many-cell traces over time (baked) ----
axA = fig.add_subplot(gs[0, 0]); b = D['baked']; tv = np.array(b['tvec']); tr = np.array(b['traces'])
for row in tr: axA.plot(tv, row, lw=0.8, alpha=0.6)
axA.plot(tv, tr.mean(axis=0), color='k', lw=2.2, label='population mean')
axA.set_xlabel('time (h)'); axA.set_ylabel('EZH2 protein (a.u.)')
axA.set_title(f"(A) {b['n']} cells over time (baked)\nemergent slow drift + cell-to-cell offsets", fontsize=10.5)
axA.legend(fontsize=8.5, loc='upper right')

# ---- B: baked variance split 3 ways (cell-to-cell / cell-cycle / bursting) ----
axB = fig.add_subplot(gs[0, 1])
Vtot = b['cv_total'] ** 2
f_c2c = b['frac_cell2cell']                                   # cell-to-cell (frozen)
cc_cv = D['noburst']['mean_within_cell_temporal_cv']         # cell-cycle-only within-cell CV (from no-burst control)
f_cyc = min(cc_cv ** 2 / Vtot, 1 - f_c2c)                    # cell-cycle share of baked total
f_burst = max(0.0, 1 - f_c2c - f_cyc)                        # bursting = remainder
parts = [('cell-to-cell\n(frozen)', f_c2c, '#3b7fb8'),
         ('cell cycle\n(deterministic)', f_cyc, '#9b59b6'),
         ('bursting\n(stochastic)', f_burst, '#c0392b')]
bot = 0
for name, frac, col in parts:
    axB.bar(0, frac, bottom=bot, color=col, width=0.6, label=f"{name}  {frac:.0%}")
    axB.text(0, bot + frac / 2, f"{frac:.0%}", ha='center', va='center', color='w', fontweight='bold', fontsize=10)
    bot += frac
axB.set_xticks([]); axB.set_ylabel('fraction of snapshot variance'); axB.set_ylim(0, 1.0); axB.set_xlim(-0.6, 0.9)
axB.set_title(f"(B) Baked split (CV {b['cv_total']:.2f})\ntemporal = cell cycle + bursting", fontsize=10.5)
axB.legend(fontsize=8.2, loc='center left', bbox_to_anchor=(0.62, 0.5), frameon=False)

# ---- C: the guess — same snapshot CV, different split ----
axC = fig.add_subplot(gs[0, 2])
cfgs = [('baked', 'baked'), ('frozen', 'mostly\nfrozen'), ('bursty', 'mostly\nbursting')]
x = np.arange(len(cfgs)); w = 0.6
for i, (k, lab) in enumerate(cfgs):
    d = D[k]; c2c = d['frac_cell2cell']; tmp = 1 - c2c
    axC.bar(i, c2c, w, color='#3b7fb8')
    axC.bar(i, tmp, w, bottom=c2c, color='#c0392b')
    axC.text(i, c2c / 2, f"{c2c:.0%}", ha='center', va='center', color='w', fontsize=8.5, fontweight='bold')
    axC.text(i, c2c + tmp / 2, f"{tmp:.0%}", ha='center', va='center', color='w', fontsize=8.5, fontweight='bold')
    axC.text(i, 1.04, f"CV {d['cv_total']:.2f}", ha='center', fontsize=9, fontweight='bold')
axC.set_xticks(x); axC.set_xticklabels([l for _, l in cfgs], fontsize=9)
axC.set_ylim(0, 1.14); axC.set_ylabel('fraction of variance')
axC.set_title("(C) THE GUESS\nbaked & mostly-frozen: SAME CV, opposite split", fontsize=10.5)
from matplotlib.patches import Patch
axC.legend(handles=[Patch(color='#3b7fb8', label='cell-to-cell'), Patch(color='#c0392b', label='within-cell (temporal)')],
           fontsize=8.2, loc='lower center', ncol=1, frameon=False)
axC.annotate('', xy=(1, 1.10), xytext=(0, 1.10), arrowprops=dict(arrowstyle='<->', color='#555'))
axC.text(0.5, 1.115, 'IF snapshot identical', ha='center', fontsize=8, color='#555', style='italic')

plt.tight_layout(rect=[0, 0, 1, 0.92])
out = f"{SC}/fig_ezh2_variance_decomp.png"
plt.savefig(out, dpi=135); print('saved', out)
print(f"baked split: cell-to-cell {f_c2c:.0%} / cell-cycle {f_cyc:.0%} / bursting {f_burst:.0%}  (total CV {b['cv_total']:.3f})")
print(f"guess: baked temporal {1-D['baked']['frac_cell2cell']:.0%} vs frozen temporal {1-D['frozen']['frac_cell2cell']:.0%} at CV {D['baked']['cv_total']:.2f}/{D['frozen']['cv_total']:.2f}")
