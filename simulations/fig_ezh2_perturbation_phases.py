"""EZH2 distribution test across perturbations: model vs IF (JP 2026-10-18). Verdict: the model gets every
drug DIRECTION right (CDK4/6i->arrest, HU->hold in G1, RO-3306->G2 block + EZH2 accumulation) but its
responses are far too COMPLETE -- it collapses ~all cells into one phase, while the IF shows modest, graded
phase shifts. Reads the model sim + the committed IF targets."""
import json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
SC = "/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/79019eb0-7c15-4650-b7d3-81347ccd6a13/scratchpad"
m = json.load(open(f"{SC}/model_ezh2_phase_dist.json"))
# IF full-population phase fractions (from extract_if_ezh2_targets)
IF = {'DMSO': [.277, .471, .125, .08], 'CDK4/6i': [.419, .367, .112, .057],
      'HU_10': [.285, .589, .084, .027], 'HU_100': [.265, .501, .123, .073], 'RO3306': [.298, .424, .123, .107]}
CONDS = ['DMSO', 'CDK4/6i', 'HU_10', 'HU_100', 'RO3306']
P = ['G0', 'G1', 'S', 'G2']; COL = {'G0': '#95a5a6', 'G1': '#3b7fb8', 'S': '#e8a33d', 'G2': '#c0392b'}

fig, (axA, axB) = plt.subplots(1, 2, figsize=(15.5, 5.2), gridspec_kw={'width_ratios': [1.7, 1]})
fig.suptitle("EZH2 distributions across perturbations: per-cell drug sensitivity -> GRADED responses (all 4 phases populated under every drug, like the IF)",
             fontsize=12, fontweight='bold')

# A: stacked phase fractions, model vs IF, per perturbation
x = np.arange(len(CONDS)); w = 0.38
for i, (src, getf) in enumerate([('model', lambda c: [m['treatments'][c]['phase_frac'][p] for p in P]),
                                  ('IF', lambda c: IF[c])]):
    off = (-w/2 if src == 'model' else w/2)
    for ci, c in enumerate(CONDS):
        bot = 0
        for p in P:
            v = getf(c)[P.index(p)]
            axA.bar(ci + off, v, w, bottom=bot, color=COL[p], edgecolor='w', lw=0.4,
                    label=p if (ci == 0 and i == 0) else None)
            bot += v
    for ci in range(len(CONDS)):
        axA.text(ci + off, 1.02, src, ha='center', fontsize=7, rotation=90, va='bottom', color='#555')
axA.set_xticks(x); axA.set_xticklabels(CONDS, fontsize=9); axA.set_ylim(0, 1.15); axA.set_ylabel('phase fraction')
axA.legend(ncol=4, fontsize=8.5, loc='upper center', bbox_to_anchor=(0.5, -0.08), frameon=False)
axA.set_title("(A) Phase fractions (left bar = model, right = IF)\nboth graded now; CDK4/6i still over-arrests a bit", fontsize=10.5)

# B: EZH2 per-phase median trajectory, model vs IF (DMSO), + drug EZH2 direction arrows
mp = m['treatments']['DMSO']['protein']
def med(d, p): return d[p]['median_log10'] if 'median_log10' in d.get(p, {}) else None
mg1 = med(mp, 'G1')
# IF DMSO per-phase medians (log10) from targets
ift = json.load(open("/Users/jpurzner/Dropbox/Q_research/py_projects/Ezh2_CyclinD1_Sphase/data/if_ezh2_phase_targets.json"))
ip = ift['protein']['DMSO']['phases']; ig1 = ip['G1']['median_log10']
for lab, g1, getm, col in [('model', mg1, lambda p: med(mp, p), '#1f6f54'),
                           ('IF', ig1, lambda p: ip[p].get('median_log10'), '#3b7fb8')]:
    ys = [10 ** (getm(p) - g1) if getm(p) is not None else np.nan for p in P]
    axB.plot(range(4), ys, '-o', color=col, lw=2, ms=6, label=lab)
axB.set_xticks(range(4)); axB.set_xticklabels(P); axB.set_yscale('log'); axB.axhline(1, color='#ccc', lw=0.8)
axB.set_ylabel('EZH2 median / G1 (log)'); axB.legend(fontsize=9)
axB.set_title("(B) EZH2 over the cycle (DMSO)\nboost shared; within-cycle model rises vs IF flat/declines", fontsize=10.5)

plt.tight_layout(rect=[0, 0, 1, 0.9])
out = f"{SC}/fig_ezh2_perturbation_phases.png"
plt.savefig(out, dpi=135); print('saved', out)
for c in CONDS:
    mfr = m['treatments'][c]['phase_frac']
    print(f"{c:9s} model {[round(mfr[p],2) for p in P]}  IF {IF[c]}")
