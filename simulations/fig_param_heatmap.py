"""Heatmap of ALL 202 v44 model parameters, grouped by module, colored by log10(value).
The 6 cell-type-varying parameters (p16, p18, P21_div, Ptch1, MYCN, mu_eff) show BOTH GNP and MB values;
everything else is a single shared engine value. JP 2026-07-27. Writes fig_param_heatmap.pdf + .png.

Reads scratchpad allparams.json (values) + param_cats.json (module/role from the categorization workflow).
"""
import os, json, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SC = '/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/81291f88-7d19-4149-be27-53b1162ee52d/scratchpad'
INK, MUT, FAINT = '#2d2d2d', '#666666', '#9a9a9a'
C_GNP, C_MB = '#2c7fb8', '#c0392b'
plt.rcParams.update({'font.size': 8, 'font.family': 'DejaVu Sans', 'pdf.fonttype': 42})

P = json.load(open(os.path.join(SC, 'allparams.json')))['params']
CATS = json.load(open(os.path.join(SC, 'param_cats.json')))
MOD_ORDER = ['Hedgehog', 'MYCN', 'CyclinD1', 'EZH2', 'H3K27me3', 'Rb-E2F', 'Skp2-p27/CKI',
             'CDK-cyclins', 'DNA-replication', 'Mitosis-APC', 'HU', 'Growth-size', 'Drug/other']

cmap = plt.cm.viridis
norm = Normalize(vmin=-4, vmax=2)   # log10(value)
def color(v):
    if v == 0 or not np.isfinite(v): return '#e6e6e6', INK
    c = cmap(norm(math.log10(abs(v))))
    lum = 0.299 * c[0] + 0.587 * c[1] + 0.114 * c[2]
    return c, ('white' if lum < 0.55 else '#111')
def vfmt(v):
    if v == 0: return '0'
    a = abs(v)
    if a >= 100 or a < 0.01: return f"{v:.0e}".replace('e-0', 'e-').replace('e+0', 'e')
    if a >= 1: return f"{v:.2f}".rstrip('0').rstrip('.')
    return f"{v:.3f}".rstrip('0').rstrip('.')

# ordered items: header rows + param rows
items = []
for mod in MOD_ORDER:
    ps = sorted([p for p in P if CATS.get(p, {}).get('module') == mod], key=str.lower)
    if not ps: continue
    items.append(('H', mod, len(ps)))
    for p in ps:
        items.append(('P', p, P[p]))
# also catch any uncategorized
uncat = [p for p in P if p not in CATS]
if uncat:
    items.append(('H', 'other', len(uncat)))
    for p in uncat: items.append(('P', p, P[p]))

NCOL = 4
rows_per = math.ceil(len(items) / NCOL)
# balance: break only between modules where possible -> simple sequential fill
cols = [[] for _ in range(NCOL)]
ci = 0
for it in items:
    if len(cols[ci]) >= rows_per and ci < NCOL - 1 and it[0] == 'H':
        ci += 1
    elif len(cols[ci]) >= rows_per + 4 and ci < NCOL - 1:
        # forced break mid-module: repeat header as (cont.)
        cur_mod = None
        for prev in reversed(cols[ci]):
            if prev[0] == 'H': cur_mod = prev[1].replace(' (cont.)', ''); break
        ci += 1
        if cur_mod: cols[ci].append(('H', cur_mod + ' (cont.)', 0))
    cols[ci].append(it)

maxrows = max(len(c) for c in cols)
fig_h = max(9.0, maxrows * 0.205 + 1.6)
fig = plt.figure(figsize=(19, fig_h))
fig.suptitle('v44 parameter map — all 202 parameters by module', x=0.03, ha='left',
             fontsize=15, fontweight='bold', color=INK, y=0.985)
fig.text(0.03, 0.958, 'one shared engine (196 params) + 6 cell-type knobs · color = log10(value) · '
         'the 6 that differ show GNP | MB', fontsize=9.5, color=MUT, ha='left')

MODCOL = {'Hedgehog': '#1b7837', 'MYCN': '#762a83', 'CyclinD1': '#2166ac', 'EZH2': '#7d5ba6',
          'H3K27me3': '#8c6bb1', 'Rb-E2F': '#b2182b', 'Skp2-p27/CKI': '#d6604d', 'CDK-cyclins': '#e08214',
          'DNA-replication': '#35978f', 'Mitosis-APC': '#4393c3', 'HU': '#878787', 'Growth-size': '#5aae61',
          'Drug/other': '#999999', 'other': '#999999'}

axw = 0.235
for ci, col in enumerate(cols):
    ax = fig.add_axes([0.03 + ci * (axw + 0.007), 0.045, axw, 0.895])
    ax.set_xlim(0, 1); ax.set_ylim(0, maxrows); ax.axis('off'); ax.invert_yaxis()
    y = 0.0
    for it in col:
        if it[0] == 'H':
            mod = it[1]; mc = MODCOL.get(mod.replace(' (cont.)', ''), '#888')
            ax.add_patch(plt.Rectangle((0, y + 0.12), 1.0, 0.76, facecolor=mc, alpha=0.16, edgecolor='none'))
            ax.text(0.012, y + 0.5, mod + (f"  ({it[2]})" if it[2] else ''), va='center', ha='left',
                    fontsize=8.4, fontweight='bold', color=mc)
            y += 1; continue
        _, name, v = it
        diff = P[name]['diff'] if isinstance(P[name], dict) else False
        # 'v' here is the dict {gnp, mb, diff, dist}
        d = P[name]
        ax.text(0.012, y + 0.5, name, va='center', ha='left', fontsize=6.9,
                color=INK, family='DejaVu Sans')
        if d['dist']:
            ax.text(0.60, y + 0.5, '~', va='center', ha='center', fontsize=9, color='#555', fontweight='bold')
        if d['diff']:
            for j, (val, ec) in enumerate([(d['gnp'], C_GNP), (d['mb'], C_MB)]):
                fc, tc = color(val)
                x0 = 0.635 + j * 0.185
                ax.add_patch(plt.Rectangle((x0, y + 0.13), 0.175, 0.74, facecolor=fc, edgecolor=ec, lw=1.4))
                ax.text(x0 + 0.088, y + 0.5, vfmt(val), va='center', ha='center', fontsize=6.2, color=tc)
        else:
            fc, tc = color(d['gnp'])
            ax.add_patch(plt.Rectangle((0.635, y + 0.13), 0.355, 0.74, facecolor=fc, edgecolor='white', lw=0.5))
            ax.text(0.8125, y + 0.5, vfmt(d['gnp']), va='center', ha='center', fontsize=6.4, color=tc)
        y += 1

# colorbar + legend
cax = fig.add_axes([0.35, 0.012, 0.3, 0.013])
cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal')
cb.set_ticks([-4, -2, 0, 2]); cb.set_ticklabels(['1e-4', '0.01', '1', '100'])
cb.ax.tick_params(labelsize=7.5, length=2); cb.set_label('parameter value (log scale)', fontsize=8, color=MUT)
fig.text(0.68, 0.018, '~ = drawn from a distribution', fontsize=8, color='#555')
fig.text(0.03, 0.018, 'GNP', color=C_GNP, fontsize=8.5, fontweight='bold')
fig.text(0.062, 0.018, '| MB borders mark the 6 cell-type-specific knobs', color=C_MB, fontsize=8.5)

out = os.path.join(ROOT, 'simulations', 'fig_param_heatmap')
fig.savefig(out + '.pdf', bbox_inches='tight'); fig.savefig(out + '.png', dpi=140, bbox_inches='tight')
print('wrote', out + '.pdf / .png  |', len([1 for it in items if it[0] == 'P']), 'params')
