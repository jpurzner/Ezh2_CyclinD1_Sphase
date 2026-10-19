"""EZH2 protein abundance distributions per cell-cycle class x treatment: model vs IF (JP 2026-10-18).
Grid: rows = treatment, cols = phase (G0/G1/S/G2). Model (green fill) + IF (blue step), each aligned to its
OWN DMSO-G1 (cycling) median so they share a 'log10 EZH2 / cycling' axis (the arbitrary a.u. differ)."""
import json, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
SC = "/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/79019eb0-7c15-4650-b7d3-81347ccd6a13/scratchpad"
m = json.load(open(f"{SC}/model_ezh2_phase_dist.json"))
t = json.load(open("/Users/jpurzner/Dropbox/Q_research/py_projects/Ezh2_CyclinD1_Sphase/data/if_ezh2_phase_targets.json"))
mc = 0.5 * (np.array(m['bins_log10'])[1:] + np.array(m['bins_log10'])[:-1])
ic = 0.5 * (np.array(t['bins_log10'])[1:] + np.array(t['bins_log10'])[:-1])
mref = m['treatments']['DMSO']['protein']['G1']['median_log10']          # model cycling ref
iref = t['protein']['DMSO']['phases']['G1']['median_log10']              # IF cycling ref
PHASES = ['G0', 'G1', 'S', 'G2']
ROWS = [('DMSO', 'DMSO'), ('CDK4/6i', 'CDK4/6i'), ('HU_10', 'HU_10'), ('HU_100', 'HU_100'), ('RO3306', 'RO3306_G2block')]

fig, ax = plt.subplots(len(ROWS), 4, figsize=(13.5, 13), sharex=True)
fig.suptitle("EZH2 protein abundance distributions per cell-cycle class x treatment  (model = green fill, IF = blue;  aligned to DMSO-G1 cycling peak at 0)",
             fontsize=12.5, fontweight='bold')
for ri, (mt, it) in enumerate(ROWS):
    mp = m['treatments'][mt]['protein']; ip = t['protein'][it]['phases']
    for ci, p in enumerate(PHASES):
        a = ax[ri, ci]
        mh = mp.get(p, {}).get('hist_density'); ih = ip.get(p, {}).get('hist_density')
        if ih is not None:
            a.fill_between(ic - iref, np.array(ih), color='#3b7fb8', alpha=0.30, lw=0, step='mid')
            a.plot(ic - iref, np.array(ih), color='#3b7fb8', lw=1.6, drawstyle='steps-mid')
        if mh is not None:
            a.fill_between(mc - mref, np.array(mh), color='#1f6f54', alpha=0.45, lw=0, step='mid')
        a.axvline(0, color='#aaa', ls=':', lw=0.8)
        a.set_yticks([]); a.set_xlim(-2.2, 1.0)
        nm = mp.get(p, {}).get('n', 0); ni = ip.get(p, {}).get('n', 0)
        a.text(0.03, 0.82, f"m{nm}/IF{ni}", transform=a.transAxes, fontsize=6.5, color='#555')
        if ri == 0: a.set_title(p, fontsize=11, fontweight='bold')
        if ci == 0: a.set_ylabel(mt, fontsize=10, fontweight='bold')
for ci in range(4): ax[-1, ci].set_xlabel('log10 EZH2 / cycling', fontsize=8.5)
from matplotlib.patches import Patch
fig.legend(handles=[Patch(color='#1f6f54', alpha=0.5, label='model'), Patch(color='#3b7fb8', alpha=0.5, label='IF (measured)')],
           loc='upper right', fontsize=10, bbox_to_anchor=(0.995, 0.965))
plt.tight_layout(rect=[0, 0, 1, 0.95])
out = f"{SC}/fig_ezh2_phase_dist_grid.png"
plt.savefig(out, dpi=130); print('saved', out)
