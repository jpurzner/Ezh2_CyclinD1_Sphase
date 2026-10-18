"""Enhanced redo figure: GNP mitogen ramp + distribution-population WITH the differentiation tripwire (JP 2026-07-21).
Reads sim_ramp_pop_diffgene_cache.pkl. Run: ./venv/bin/python simulations/fig_ramp_pop_diffgene.py"""
import os, pickle, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R = pickle.load(open(os.path.join(HERE, 'sim_ramp_pop_diffgene_cache.pkl'), 'rb'))
C_DR, C_CD, C_RATE, C_MK, C_MKD, C_DIFF = '#7f8c8d', '#2e86c1', '#e67e22', '#5b2c6f', '#8e44ad', '#c0392b'

fig, ax = plt.subplots(2, 2, figsize=(15.4, 10.6))
fig.subplots_adjust(left=0.07, right=0.955, top=0.88, bottom=0.07, hspace=0.32, wspace=0.28)


def ramp_panel(a, d, drive_label, title):
    sp = d['split']
    a.axvline(d['t'][sp], color='#bbb', ls='-', lw=1)
    a.text(d['t'][sp] * 0.5, 1.04, 'ramp UP', ha='center', fontsize=8.5, color='#888')
    a.text(d['t'][sp] + (d['t'][-1] - d['t'][sp]) * 0.5, 1.04, 'ramp DOWN', ha='center', fontsize=8.5, color='#888')
    dn = d['drive'] / d['drive'].max()
    a.plot(d['t'], dn, color=C_DR, lw=1.5, label=drive_label + ' (norm)')
    a.plot(d['t'], d['Cd'] / max(d['Cd'].max(), 1e-9), color=C_CD, lw=1.8, label='CyclinD1 (norm)')
    a.plot(d['t'], d['divrate'] / max(d['divrate'].max(), 1e-9), color=C_RATE, lw=1.5, alpha=0.8, label='division rate (norm)')
    a.plot(d['t'], d['Mk'], color=C_MK, lw=1.8, label='H3K27me3 @CyclinD1')
    a.plot(d['t'], d['Mk_diff'], color=C_MKD, lw=2.0, label='H3K27me3 @diff-gene (me3)')
    a.plot(d['t'], d['Diff'], color=C_DIFF, lw=2.4, label='diff-gene activation (Diff)')
    a.set_xlabel('time (h)'); a.set_ylabel('normalized / level'); a.set_ylim(-0.03, 1.1)
    a.set_title(title, fontsize=10.6, fontweight='bold'); a.grid(alpha=0.2); a.legend(fontsize=7.2, loc='center left', ncol=1)


ramp_panel(ax[0, 0], R['ramp_shh'], 'SHH',
           'A · Physiological SHH ramp up/down\nCyclinD1 governed by the mark; diff-gene tripwire')
ramp_panel(ax[0, 1], R['ramp_mu'], 'growth μ',
           'B · OVERDRIVE ramp (growth-rate μ up/down)\npush division fast → me3@diff collapses → PREMATURE activation')

# C: population premature-activation fraction vs growth drive
c = ax[1, 0]; P = R['pop']; mu = P['mu'] * 1e3
A = P['A']; frac = (A[:, :, 2] > 0.5).mean(axis=0) * 100
mediandiv = np.nanmedian(A[:, :, 0], axis=0)
c.plot(mu, frac, '-o', color=C_DIFF, lw=2.6, ms=7, label='% cells premature-activated')
c.set_xlabel('growth drive  μ ×10³  (→ faster division)'); c.set_ylabel('% cells with diff-gene ACTIVATED', color=C_DIFF)
c.tick_params(axis='y', labelcolor=C_DIFF); c.set_ylim(-3, 100)
c.set_title('C · Population: premature activation rises with how hard\ncells are pushed (a tail trips even at normal drive)', fontsize=10.6, fontweight='bold')
c.grid(alpha=0.2)
c2 = c.twinx(); c2.plot(mu, mediandiv, '-s', color=C_RATE, lw=1.8, ms=5, alpha=0.8)
c2.set_ylabel('median division rate (÷/day)', color=C_RATE, fontsize=9); c2.tick_params(axis='y', labelcolor=C_RATE)
c.legend(fontsize=8, loc='upper left')

# D: per-cell me3_diff vs division rate, colored by Diff, at the highest-drive column
d = ax[1, 1]; col = -1  # highest mu
div = A[:, col, 0]; me3d = A[:, col, 1]; Dcol = A[:, col, 2]
sc = d.scatter(div, me3d, c=Dcol, cmap='coolwarm', s=60, vmin=0, vmax=1, edgecolor='#444', linewidth=0.4)
cb = fig.colorbar(sc, ax=d, fraction=0.046); cb.set_label('diff-gene activation (Diff)', fontsize=9)
d.axhline(0.28, color='#888', ls='--', lw=1.3); d.text(div.max() * 0.98, 0.30, 'activation threshold', ha='right', fontsize=8, color='#666')
d.set_xlabel('per-cell division rate (÷/day)'); d.set_ylabel('H3K27me3 at diff-gene (me3)')
d.set_title('D · Per-cell (at high drive): the fastest cells dilute me3\nbelow threshold → the gene fires (blue→red)', fontsize=10.6, fontweight='bold')
d.grid(alpha=0.2)

fig.suptitle('Redo — GNP mitogen ramp + distribution-population, now through the speed-limiter with the differentiation tripwire',
             fontsize=12.5, fontweight='bold', y=0.955)
fig.text(0.5, 0.915, 'CyclinD1 is governed by its H3K27me3 mark; the poised differentiation gene stays repressed unless division is pushed past the sustainable rate (overdrive B, high-drive cells C/D), then it prematurely activates. Default OFF; validation 28/29.',
         ha='center', fontsize=8.8, color='#444')
out = os.path.join(HERE, 'fig_ramp_pop_diffgene.png')
fig.savefig(out, dpi=140); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
