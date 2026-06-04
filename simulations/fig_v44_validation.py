"""Supp Fig 7 B-M: v44 model validation grid (model vs experimental data).

Panels (matching the manuscript caption):
  B  GNP Cyclin B traces (+SHH, -SHH, +HHi, +EZH2i, +HHi+EZH2i)
  C  GNP CyclinD1 mRNA fold-change across conditions  (model vs data: +HHi=0.14)
  D  GNP MYCN model vs data (GNP, GNP+HHi=0.78)
  E  GNP EZH2 (cycling / serum-starved G0 / +HHi)
  F  MB Cyclin B traces (MB, +HHi, +EZH2i, +HHi+EZH2i)
  G  MB CyclinD1 mRNA fold-change across conditions  (model vs data: +HHi=0.40)
  H  MB MYCN model vs data (MB, MB+HHi=0.86)
  I  MB Gli1 model vs data (>99% reduction with HHi)
  J  division counts, all conditions
  K  MB CyclinD1 mRNA time courses
  L  MB EZH2 protein time courses
  M  MYCN protein (GNP+SHH, GNP+HHi, MB, MB+HHi)

Run:  ./venv/bin/python simulations/fig_v44_validation.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from simulations.validate_v44 import run, count_divisions, mean_settled, CONDITIONS

sims = {name: run(**kw) for name, kw in CONDITIONS.items()}
TH = lambda r: r['time'] / 60.0

GNP = ['GNP + SHH', 'GNP - SHH', 'GNP + HHi', 'GNP + EZH2i', 'GNP + HHi + EZH2i']
GNP_C = ['#1b9e77', '#999999', '#e7298a', '#66a61e', '#d95f02']
MB = ['MB', 'MB + HHi', 'MB + EZH2i', 'MB + HHi + EZH2i']
MB_C = ['#762A83', '#E08214', '#D6604D', '#4A1486']

cd = lambda c: mean_settled(sims[c], 'Cd_mRNA')
my = lambda c: mean_settled(sims[c], 'MYCN')
gl = lambda c: mean_settled(sims[c], 'Gli1')
ez = lambda c: mean_settled(sims[c], 'EZH2')

fig = plt.figure(figsize=(19, 15))
gs = GridSpec(4, 4, figure=fig, hspace=0.5, wspace=0.42)


def lab(ax, L, title):
    ax.set_title(f"({L}) {title}", fontsize=10, fontweight='bold')


# B: GNP CycB traces
ax = fig.add_subplot(gs[0, 0])
for n, c in zip(GNP, GNP_C):
    ax.plot(TH(sims[n]), sims[n]['Cb'], color=c, lw=0.8,
            label=n.replace('GNP + ', '').replace('GNP ', ''))
lab(ax, 'B', 'GNP: Cyclin B'); ax.set_xlabel('Time (h)'); ax.set_ylabel('Cyclin B')
ax.legend(fontsize=6); ax.grid(alpha=0.2)

# C: GNP CyclinD1 mRNA fold change
ax = fig.add_subplot(gs[0, 1])
base = cd('GNP + SHH')
vals = [1.0, cd('GNP + HHi')/base, cd('GNP + EZH2i')/base, cd('GNP + HHi + EZH2i')/base]
data = [1.0, 0.14, 2.2, None]
x = np.arange(4)
ax.bar(x - 0.17, vals, 0.34, color=GNP_C[:1]+GNP_C[2:5], edgecolor='k', alpha=0.8, label='Model')
for i, d in enumerate(data):
    if d is not None: ax.scatter(i + 0.17, d, color='k', s=55, marker='D', zorder=5)
lab(ax, 'C', 'GNP: CyclinD1 mRNA'); ax.set_xticks(x)
ax.set_xticklabels(['SHH', '+HHi', '+EZH2i', '+Both'], fontsize=8)
ax.set_ylabel('fold of GNP+SHH'); ax.legend(['Data', 'Model'], fontsize=7); ax.grid(alpha=0.2, axis='y')

# D: GNP MYCN
ax = fig.add_subplot(gs[0, 2])
vals = [1.0, my('GNP + HHi')/my('GNP + SHH')]; data = [1.0, 0.78]
x = np.arange(2)
ax.bar(x - 0.17, vals, 0.34, color=['#1b9e77', '#e7298a'], edgecolor='k', alpha=0.8)
for i, d in enumerate(data):
    ax.scatter(i + 0.17, d, color='k', s=55, marker='D', zorder=5)
    ax.text(i + 0.27, d, f'{d:.2f}', fontsize=7, va='center')
lab(ax, 'D', 'GNP: MYCN'); ax.set_xticks(x); ax.set_xticklabels(['GNP', '+HHi'], fontsize=8)
ax.set_ylabel('fold of GNP'); ax.grid(alpha=0.2, axis='y')

# E: GNP EZH2 cycling / G0 / HHi
ax = fig.add_subplot(gs[0, 3])
base = ez('GNP + SHH')
vals = [1.0, ez('GNP Serum-starved')/base, ez('GNP + HHi')/base]
ax.bar(range(3), vals, color=['#1b9e77', '#7570b3', '#e7298a'], edgecolor='k', alpha=0.8)
ax.axhline(0.6, color='red', ls='--', alpha=0.5, label='target G0 ~0.6')
lab(ax, 'E', 'GNP: EZH2'); ax.set_xticks(range(3))
ax.set_xticklabels(['Cycling', 'G0', '+HHi'], fontsize=8); ax.set_ylabel('fold of cycling')
ax.legend(fontsize=7); ax.grid(alpha=0.2, axis='y')

# F: MB CycB traces
ax = fig.add_subplot(gs[1, 0])
for n, c in zip(MB, MB_C):
    ax.plot(TH(sims[n]), sims[n]['Cb'], color=c, lw=0.8, label=n.replace('MB + ', '+').replace('MB', 'Untx'))
lab(ax, 'F', 'MB: Cyclin B'); ax.set_xlabel('Time (h)'); ax.set_ylabel('Cyclin B')
ax.legend(fontsize=7); ax.grid(alpha=0.2)

# G: MB CyclinD1 mRNA fold change
ax = fig.add_subplot(gs[1, 1])
base = cd('MB')
vals = [1.0, cd('MB + HHi')/base, cd('MB + EZH2i')/base, cd('MB + HHi + EZH2i')/base]
data = [1.0, 0.40, None, None]
x = np.arange(4)
ax.bar(x - 0.17, vals, 0.34, color=MB_C, edgecolor='k', alpha=0.8)
for i, d in enumerate(data):
    if d is not None:
        ax.scatter(i + 0.17, d, color='k', s=55, marker='D', zorder=5)
        ax.text(i + 0.27, d, f'{d:.2f}', fontsize=7, va='center')
lab(ax, 'G', 'MB: CyclinD1 mRNA'); ax.set_xticks(x)
ax.set_xticklabels(['MB', '+HHi', '+EZH2i', '+Both'], fontsize=8)
ax.set_ylabel('fold of MB'); ax.grid(alpha=0.2, axis='y')

# H: MB MYCN
ax = fig.add_subplot(gs[1, 2])
g_my = my('GNP + SHH')
vals = [my('MB')/g_my, my('MB + HHi')/g_my]; data = [2.8, 2.8 * 0.86]
x = np.arange(2)
ax.bar(x - 0.17, vals, 0.34, color=['#762A83', '#E08214'], edgecolor='k', alpha=0.8)
for i, d in enumerate(data):
    ax.scatter(i + 0.17, d, color='k', s=55, marker='D', zorder=5)
    ax.text(i + 0.27, d, f'{d:.2f}', fontsize=7, va='center')
lab(ax, 'H', 'MB: MYCN'); ax.set_xticks(x); ax.set_xticklabels(['MB', '+HHi'], fontsize=8)
ax.set_ylabel('fold of GNP'); ax.grid(alpha=0.2, axis='y')

# I: MB Gli1
ax = fig.add_subplot(gs[1, 3])
base = gl('MB')
vals = [1.0, gl('MB + HHi')/base if base else 0]
data = [1.0, 0.01]
x = np.arange(2)
ax.bar(x - 0.17, vals, 0.34, color=['#762A83', '#E08214'], edgecolor='k', alpha=0.8)
for i, d in enumerate(data):
    ax.scatter(i + 0.17, d, color='k', s=55, marker='D', zorder=5)
lab(ax, 'I', 'MB: Gli1'); ax.set_xticks(x); ax.set_xticklabels(['MB', '+HHi'], fontsize=8)
ax.set_ylabel('fold of MB'); ax.grid(alpha=0.2, axis='y')

# J: division counts
ax = fig.add_subplot(gs[2, 0])
allc = ['GNP + SHH', 'GNP - SHH', 'GNP + HHi', 'GNP + EZH2i', 'MB', 'MB + HHi', 'MB + EZH2i', 'MB + HHi + EZH2i']
allcol = ['#1b9e77', '#999999', '#e7298a', '#66a61e', '#762A83', '#E08214', '#D6604D', '#4A1486']
divs = [count_divisions(sims[c])[0] for c in allc]
bars = ax.bar(range(len(allc)), divs, color=allcol, edgecolor='k')
for b, v in zip(bars, divs):
    ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.2, str(v), ha='center', fontsize=8, fontweight='bold')
lab(ax, 'J', 'Cell divisions (200 h)'); ax.set_xticks(range(len(allc)))
ax.set_xticklabels(['GNP', '-SHH', '+HHi', '+EZH2i', 'MB', 'MB\n+HHi', 'MB\n+EZH2i', 'MB\n+Both'], fontsize=6)
ax.set_ylabel('divisions'); ax.grid(alpha=0.2, axis='y')

# K: MB CyclinD1 mRNA time courses
ax = fig.add_subplot(gs[2, 1])
for n, c in zip(MB, MB_C):
    ax.plot(TH(sims[n]), sims[n]['Cd_mRNA'], color=c, lw=0.8, label=n.replace('MB + ', '+').replace('MB', 'Untx'))
lab(ax, 'K', 'MB: CyclinD1 mRNA'); ax.set_xlabel('Time (h)'); ax.set_ylabel('CyclinD1 mRNA')
ax.legend(fontsize=7); ax.grid(alpha=0.2)

# L: MB EZH2 protein time courses
ax = fig.add_subplot(gs[2, 2])
for n, c in zip(MB, MB_C):
    ax.plot(TH(sims[n]), sims[n]['EZH2'], color=c, lw=0.8, label=n.replace('MB + ', '+').replace('MB', 'Untx'))
lab(ax, 'L', 'MB: EZH2 protein'); ax.set_xlabel('Time (h)'); ax.set_ylabel('EZH2 protein')
ax.legend(fontsize=7); ax.grid(alpha=0.2)

# M: MYCN protein dynamics
ax = fig.add_subplot(gs[2, 3])
mycn_c = ['GNP + SHH', 'GNP + HHi', 'MB', 'MB + HHi']
mycn_col = ['#1b9e77', '#e7298a', '#762A83', '#E08214']
for n, c in zip(mycn_c, mycn_col):
    ax.plot(TH(sims[n]), sims[n]['MYCN'], color=c, lw=0.8, label=n)
lab(ax, 'M', 'MYCN protein'); ax.set_xlabel('Time (h)'); ax.set_ylabel('MYCN')
ax.legend(fontsize=6); ax.grid(alpha=0.2)

fig.suptitle('Supp Fig 7 B-M: v44 model validation (model vs experimental data)',
             fontsize=14, fontweight='bold', y=1.005)
plt.savefig('simulations/fig_v44_validation.png', dpi=200, bbox_inches='tight')
plt.savefig('simulations/fig_v44_validation.pdf', bbox_inches='tight')
plt.close()
print("Saved: fig_v44_validation.png/pdf")
