"""EZH2 per-phase distribution test: model vs IF (JP 2026-10-17). The test's verdict: the model reproduces
the BRIGHT cycling population (G1 width matches) but is missing the ~28% G0-low (EZH2-cleared) quiescent
fraction the IF has -> it cannot show the G0->cycling boost without more reversible-G0 cells."""
import json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
SC = "/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/79019eb0-7c15-4650-b7d3-81347ccd6a13/scratchpad"
m = json.load(open(f"{SC}/model_ezh2_phase_dist.json"))
t = json.load(open("/Users/jpurzner/Dropbox/Q_research/py_projects/Ezh2_CyclinD1_Sphase/data/if_ezh2_phase_targets.json"))
IFFR = {'G0': 0.277, 'G1': 0.471, 'S': 0.125, 'G2': 0.08}
mf = m['treatments']['DMSO']['phase_frac']

fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
fig.suptitle("EZH2 per-phase distribution test (DMSO): model reproduces cycling cells but lacks the ~28% G0-low (EZH2-cleared) fraction",
             fontsize=12, fontweight='bold')

# A: phase fractions
P = ['G0', 'G1', 'S', 'G2']; x = np.arange(4); w = 0.38
ax[0].bar(x - w/2, [mf[p] for p in P], w, color='#1f6f54', label='model')
ax[0].bar(x + w/2, [IFFR[p] for p in P], w, color='#3b7fb8', label='IF')
ax[0].set_xticks(x); ax[0].set_xticklabels(P); ax[0].set_ylabel('fraction'); ax[0].legend(fontsize=9)
ax[0].set_title("(A) Phase fractions\nmodel G0 2.5% vs IF 28% — the gap", fontsize=10.5)
ax[0].annotate('G0 gap', xy=(0, 0.28), xytext=(0.5, 0.45), fontsize=9, color='#c0392b',
               arrowprops=dict(arrowstyle='->', color='#c0392b'))

# B: pooled EZH2 protein distribution (reconstruct from per-phase hists, n-weighted), median-aligned on G1
bm = np.array(m['bins_log10']); cm = 0.5 * (bm[1:] + bm[:-1])
bi = np.array(t['bins_log10']); ci = 0.5 * (bi[1:] + bi[:-1])
def pooled(d, bins_c, frac):
    tot = np.zeros(len(bins_c))
    for p in P:
        ph = d[p] if 'median_log10' in d.get(p, {}) else None
        if ph: tot += frac[p] * np.array(ph['hist_density'])
    return tot / max(tot.sum(), 1e-9)
mp = m['treatments']['DMSO']['protein']; ip = t['protein']['DMSO']['phases']
pm = pooled(mp, cm, mf); pi = pooled(ip, ci, IFFR)
mg1 = mp['G1']['median_log10']; ig1 = ip['G1']['median_log10']
ax[1].plot(ci - ig1, pi, color='#3b7fb8', lw=2.2, label='IF (bimodal: G0-low + cycling)')
ax[1].fill_between(ci - ig1, pi, color='#3b7fb8', alpha=0.2)
ax[1].plot(cm - mg1, pm, color='#1f6f54', lw=2.2, label='model (cycling only)')
ax[1].axvline(0, color='#888', ls=':', lw=1)
ax[1].set_xlabel('log10 EZH2 (aligned to G1 peak)'); ax[1].set_ylabel('density'); ax[1].legend(fontsize=8.5)
ax[1].set_title("(B) Protein distribution\nIF has a big G0-low mode the model lacks", fontsize=10.5)
ax[1].set_xlim(-3, 1.2)

# C: per-phase median trajectory (normalized to G1), model vs IF, + G1 width match
def med(d, p): return d[p]['median_log10'] if 'median_log10' in d.get(p, {}) else None
for lab, d, col, g1 in [('model', mp, '#1f6f54', mg1), ('IF', ip, '#3b7fb8', ig1)]:
    ys = [10 ** (med(d, p) - g1) if med(d, p) is not None else np.nan for p in P]
    ax[2].plot(range(4), ys, '-o', color=col, lw=2, ms=6, label=lab)
ax[2].set_xticks(range(4)); ax[2].set_xticklabels(P); ax[2].set_yscale('log')
ax[2].axhline(1, color='#ccc', lw=0.8); ax[2].set_ylabel('EZH2 median / G1 (log)'); ax[2].legend(fontsize=9)
ax[2].set_title(f"(C) Over the cycle (median/G1)\nmodel G1 sd {mp['G1']['log10_sd']:.2f} ~ IF-bright 0.17; IF G0 cleared", fontsize=10.5)

plt.tight_layout(rect=[0, 0, 1, 0.92])
out = f"{SC}/fig_ezh2_phase_dist_vs_if.png"
plt.savefig(out, dpi=135); print('saved', out)
print(f"model phase-frac G0 {mf['G0']:.3f} vs IF 0.277 | model G1 log10_sd {mp['G1']['log10_sd']:.3f} vs IF-bright ~0.17")
print(f"within-cycle protein G2/G1: model {10**(med(mp,'G2')-mg1):.2f} vs IF {10**(med(ip,'G2')-ig1):.2f}")
