"""Figure: resolving the PRC2/transcription paradox.
A. The paradox + first resolution: mitogen drives GLOBAL EZH2 UP (E2F-target, peaks at division) yet the LOCAL
   H3K27me3 at CyclinD1 DOWN (evicted by the higher transcription) -> high PRC2 + high expression coexist; ChIP MB<GNP.
B. The deep resolution: the transcription<->PRC2 antagonism is a bistable SWITCH (no graded intermediate in mean-
   field). But a single ALLELE is bistable with a STABLE mark (ON=high / OFF=low); PUNCTATE transcription flips
   alleles to OFF more often at higher mitogen; the POPULATION AVERAGE is then a smooth, graded MB<GNP + dose-graded
   CyclinD1 -- with NO fast eraser (each allele's mark is stable, latched). Self-contained. Run: python simulations/fig_prc2_transcription_paradox.py"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- A: mean-field locus (global EZH2 up, local mark down) ---
def steady(mit, EZH2, g=0.95, he=3, M0=0.7, repmax=0.8, a_rw=0.9, Krep=0.36, Kev=2.4, kw=0.021, a0=0.018,
           hr=4, delta=0.0007, mu=0.0006):
    M = M0
    for _ in range(120000):
        rep = repmax * M**hr / (Krep**hr + M**hr); T = mit * (1 - rep); ev = g * T**he / (Kev**he + T**he)
        M = min(1.0, max(0.0, M + (kw * EZH2 * (a0 + a_rw * M) * (1 - M) * (1 - ev) - delta * M - mu * M) * 0.5))
    return M, mit * (1 - repmax * M**hr / (Krep**hr + M**hr))
mit = np.array([2, 3, 4, 6, 8, 10, 12, 15, 18, 22]); ez = 2.0 + 0.11 * mit
mf = np.array([steady(m, e) for m, e in zip(mit, ez)])

# --- B: allele-population model. Each allele bistable (ON mark 0.9 / OFF mark 0.08, STABLE). Punctate transcription
#     flips ON->OFF at a mitogen-tuned rate; PRC2 re-methylation flips OFF->ON at a constant rate. Population fraction
#     OFF (and thus average mark, CyclinD1) grades SMOOTHLY with mitogen. ---
M_ON, M_OFF, K_SW, H_SW, LEAK = 0.90, 0.08, 9.0, 3.0, 0.15
f_off = (mit / K_SW)**H_SW / (1 + (mit / K_SW)**H_SW)          # population fraction of OFF (active) alleles
avg_mark = (1 - f_off) * M_ON + f_off * M_OFF                   # population-average H3K27me3 (what ChIP sees)
cyc = mit * ((1 - f_off) * LEAK + f_off)                        # CyclinD1: OFF alleles full, ON alleles leaky

fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.8, 5.6))
fig.subplots_adjust(left=0.07, right=0.93, top=0.82, bottom=0.14, wspace=0.44)

c_ez, c_mk, c_cd = '#6c3483', '#1e5638', '#c0392b'
axA.plot(mit, ez / ez.max(), '-^', color=c_ez, lw=2.3, ms=6, label='GLOBAL EZH2 (E2F-driven)')
axA.plot(mit, mf[:, 0], '-o', color=c_mk, lw=2.5, ms=6, label='LOCAL H3K27me3 at CyclinD1')
axA.plot(mit, mf[:, 1] / mf[:, 1].max(), '-s', color=c_cd, lw=2.3, ms=6, label='CyclinD1 expression')
axA.set_xlabel('mitogen (Hh / CyclinD1 transcription →)', fontsize=11); axA.set_ylabel('normalized level', fontsize=11)
axA.set_title('A. Global EZH2 UP, local mark DOWN — the paradox dissolves', fontsize=11.5, fontweight='bold')
axA.legend(fontsize=9, frameon=False, loc='center left'); axA.set_ylim(0, 1.08)
for s in ('top', 'right'): axA.spines[s].set_visible(False)
axA.text(8.5, 0.12, '"high PRC2 + high expression" = global EZH2 is\nE2F-driven (peaks at division); at the CyclinD1\nlocus the transcription evicts it → low mark.',
         fontsize=8.2, color='#333', va='bottom')

axB.plot(mit, avg_mark, '-o', color=c_mk, lw=2.6, ms=6, label='population-avg H3K27me3 (ChIP)')
axB.axhline(M_ON, color=c_mk, ls=':', lw=1, alpha=0.5); axB.axhline(M_OFF, color=c_mk, ls=':', lw=1, alpha=0.5)
axB.text(21.5, M_ON - 0.05, 'allele ON (latched, stable)', fontsize=7.5, color=c_mk, ha='right')
axB.text(21.5, M_OFF + 0.02, 'allele OFF (latched, stable)', fontsize=7.5, color=c_mk, ha='right')
axB2 = axB.twinx()
axB2.plot(mit, cyc, '-s', color=c_cd, lw=2.4, ms=6, label='CyclinD1 (population)')
axB2.set_ylabel('CyclinD1 (population)', color=c_cd, fontsize=11); axB2.tick_params(axis='y', labelcolor=c_cd)
axB.set_ylabel('population-avg H3K27me3', color=c_mk, fontsize=11); axB.tick_params(axis='y', labelcolor=c_mk)
axB.axvline(8, color='#888', ls=':', lw=1); axB.axvline(15, color='#888', ls=':', lw=1)
axB.text(8, 0.62, 'GNP', ha='center', fontsize=9, color='#555', fontweight='bold'); axB.text(15, 0.30, 'MB', ha='center', fontsize=9, color='#555', fontweight='bold')
axB.set_xlabel('mitogen (transcription →)', fontsize=11); axB.set_ylim(0, 1.0)
axB.set_title('B. Resolution: bistable alleles + punctate flips → smooth population', fontsize=11.5, fontweight='bold')
for s in ('top',): axB.spines[s].set_visible(False)
axB.text(2.3, 0.45, 'each allele is bistable with a STABLE mark;\npunctate transcription flips alleles to OFF more\noften at higher mitogen → the POPULATION average\ngrades smoothly (MB<GNP) + CyclinD1 dose-graded.\nNO fast eraser needed.',
         fontsize=8.2, color='#333', va='top')
print(f"  population avg mark: GNP(mit8) {avg_mark[4]:.2f}  MB(mit15) {avg_mark[7]:.2f}  MB/GNP {avg_mark[7]/avg_mark[4]:.2f}; CyclinD1 fold {cyc[7]/cyc[4]:.2f}")

fig.suptitle('Resolving the PRC2 / transcription paradox — global EZH2 (E2F) vs local eviction; a STABLE mark, no fast eraser',
             fontsize=12.6, fontweight='bold', y=0.95)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_prc2_transcription_paradox.png')
fig.savefig(out, dpi=150); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
