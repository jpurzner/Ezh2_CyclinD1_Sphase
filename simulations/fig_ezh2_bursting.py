"""Figure for the transcriptional-bursting EZH2 model (JP 2026-10-09). Three panels:
  (A) single-cell EZH2 traces over time  -> slow drift from bursty transcription
  (B) DMSO model vs IF abundance distribution (model median-aligned to the IF bright/G1 peak)
  (C) drug response (median EZH2 / DMSO): model vs IF, IF ratios computed from per-phase G1 medians
Reads the scratchpad JSONs; writes simulations/fig_ezh2_bursting.png.
Run:  venv/bin/python simulations/fig_ezh2_bursting.py
"""
import json, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SC = "/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/79019eb0-7c15-4650-b7d3-81347ccd6a13/scratchpad"
m = json.load(open(f"{SC}/model_ezh2_bursting.json"))
iftgt = json.load(open(f"{SC}/if_EZH2_target.json"))
hist = json.load(open(f"{SC}/ezh2_hist.json"))

DMSO = m['treatments']['DMSO']
cv = DMSO['cv']; g1cv = DMSO['g1_cv']; duty = m['params']['DUTY']; tau = m['params']['TAU_BURST']
mdl_log10 = np.array(DMSO['log10_ezh2'])
mdl_med = np.median(mdl_log10)

# IF reference: bright/G1 population median (log10) from the per-phase histograms
if_g1 = np.array(hist['data']['DMSO_0_24h']['G1'])
if_g1_med = np.median(if_g1)
if_raw = np.array(iftgt['raw_log10_sample'])        # pooled 8000-cell sample (bimodal, incl. dim lump)
# IF bright-pop CV from the AMOS mixture fit (log10_sd 0.171 -> CV ~0.41); G1 CV ~0.45
IF_BRIGHT_CV = 0.45

# align model bright bulk onto IF bright peak
shift = if_g1_med - mdl_med
mdl_aligned = mdl_log10 + shift

# ---- panel C: drug ratios ----
def g1_lin_median(cond):
    return 10 ** np.median(np.array(hist['data'][cond]['G1']))
if_dmso = g1_lin_median('DMSO_0_24h')
if_ratio = {'CDK4/6i': g1_lin_median('Pal_1_24h') / if_dmso,
            'HU': g1_lin_median('HU_10_24h') / if_dmso}
mdl_ratio = {'CDK4/6i': m['treatments']['CDK4/6i']['median'] / DMSO['median'],
             'HU': m['treatments']['HU']['median'] / DMSO['median']}

fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(15, 4.6))
fig.suptitle(f"EZH2 transcriptional bursting (variable production rate) -> slow single-cell drift; "
             f"CV {cv:.2f} (target ~0.45); drug response preserved", fontsize=12, fontweight='bold')

# (A) traces
for tr in m['traces']:
    axA.plot(tr['t'], tr['ez'], lw=1.3, alpha=0.85)
axA.set_xlabel('time (h)'); axA.set_ylabel('EZH2 (a.u.)')
axA.set_title(f"(A) Single-cell EZH2 over time\n(bursty transcription, tau~{tau:.0f}min -> slow drift)")

# (B) distribution
bins = np.linspace(2.3, 5.0, 46)
axB.hist(if_raw, bins=bins, density=True, color='#8fb8d8', alpha=0.65,
         label=f"IF measured (G1 CV ~{IF_BRIGHT_CV:.2f})")
axB.hist(mdl_aligned, bins=bins, density=True, histtype='step', color='#1f6f54', lw=2.0,
         label=f"model (CV {cv:.2f}, G1 {g1cv:.2f})")
axB.set_xlabel('log10 EZH2 (aligned)'); axB.set_ylabel('density')
axB.set_title("(B) DMSO: model vs IF"); axB.legend(fontsize=9)

# (C) drug response
labels = ['DMSO', 'CDK4/6i', 'HU']
mdl_bars = [1.0, mdl_ratio['CDK4/6i'], mdl_ratio['HU']]
if_bars = [1.0, if_ratio['CDK4/6i'], if_ratio['HU']]
x = np.arange(3); w = 0.38
axC.bar(x - w/2, mdl_bars, w, color='#1f6f54', label='model')
axC.bar(x + w/2, if_bars, w, color='#3b7fb8', label='IF (G1)')
axC.set_xticks(x); axC.set_xticklabels(labels)
axC.set_ylabel('EZH2 median / DMSO'); axC.axhline(1.0, color='k', lw=0.6, alpha=0.4)
axC.set_title("(C) Drug response: model vs IF"); axC.legend(fontsize=9)

plt.tight_layout(rect=[0, 0, 1, 0.94])
out = "/Users/jpurzner/Dropbox/Q_research/py_projects/Ezh2_CyclinD1_Sphase/simulations/fig_ezh2_bursting.png"
plt.savefig(out, dpi=130); print('saved', out)
print(f"model: CV {cv:.3f} G1_CV {g1cv:.3f} | model ratios CDK4/6i {mdl_ratio['CDK4/6i']:.2f} HU {mdl_ratio['HU']:.2f}")
print(f"IF(G1) ratios: CDK4/6i {if_ratio['CDK4/6i']:.2f} HU {if_ratio['HU']:.2f} | IF G1 median log10 {if_g1_med:.3f}")
