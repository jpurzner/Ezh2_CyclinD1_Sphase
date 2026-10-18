"""SVD biplot of the parameter->target sensitivity matrix (JP 2026-07-27). The right singular vectors are the
model's major axes of parameter-driven variation; parameters that point the same way are aligned, opposite =
antagonistic. Reads scratchpad sensitivity.json + param_meta.json. Writes fig_param_biplot.pdf + .png."""
import os, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SC = '/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/81291f88-7d19-4149-be27-53b1162ee52d/scratchpad'
INK, MUT, FAINT = '#2d2d2d', '#666', '#a8a8a8'
MODCOL = {'Hedgehog': '#1b7837', 'MYCN': '#762a83', 'CyclinD1': '#2166ac', 'EZH2': '#7d5ba6',
          'H3K27me3': '#8c6bb1', 'Rb-E2F': '#b2182b', 'Skp2-p27/CKI': '#d6604d',
          'DNA-replication': '#35978f', 'HU': '#e08214', 'Growth-size': '#5aae61'}
plt.rcParams.update({'font.size': 9, 'font.family': 'DejaVu Sans', 'pdf.fonttype': 42})

S = json.load(open(SC + '/sensitivity.json')); META = json.load(open(SC + '/param_meta.json'))
E, TG = S['elasticity'], S['targets']
P = [p for p in E]
M = np.array([[E[p][t] for t in TG] for p in P])           # params x targets
U, s, Vt = np.linalg.svd(M, full_matrices=False)
var = s**2 / np.sum(s**2)
infl = M.__abs__().sum(1)

def short(t):
    t = t.split('(')[0].strip()
    for k, v in {'CyclinD1':'CycD1','reduction':'red','count%':'%','duration':'dur','transcript':'mRNA',
                 'protein':'prot','GNP+HHi':'G+HHi','MB+HHi':'M+HHi',' boost':''}.items(): t = t.replace(k, v)
    return t

def biplot(ax, a, b):
    # orient each axis so its largest target loading is positive (stable, readable)
    Ua, Va = U[:, a].copy(), Vt[a].copy()
    Ub, Vb = U[:, b].copy(), Vt[b].copy()
    if Va[np.argmax(np.abs(Va))] < 0: Ua, Va = -Ua, -Va
    if Vb[np.argmax(np.abs(Vb))] < 0: Ub, Vb = -Ub, -Vb
    px, py = Ua * s[a], Ub * s[b]
    sc = 0.9 * max(np.abs(np.r_[px, py]).max(), 1e-9)
    # target reference arrows (top-loading targets in this plane)
    tload = np.hypot(Va, Vb); ti = np.argsort(-tload)[:8]
    tsc = sc / max(tload[ti].max(), 1e-9) * 0.95
    for i in ti:
        ax.annotate('', (Va[i]*tsc, Vb[i]*tsc), (0, 0),
                    arrowprops=dict(arrowstyle='-|>', color='#c9c9c9', lw=1.3, shrinkA=0, shrinkB=0), zorder=1)
        ax.text(Va[i]*tsc*1.08, Vb[i]*tsc*1.08, short(TG[i]), fontsize=6.6, color='#999',
                ha='center', va='center', style='italic')
    ax.axhline(0, color='#eee', lw=1, zorder=0); ax.axvline(0, color='#eee', lw=1, zorder=0)
    for i, p in enumerate(P):
        c = MODCOL.get(META.get(p, {}).get('module', ''), '#888')
        ms = 4 + 22 * (infl[i] / infl.max())
        ax.scatter(px[i], py[i], s=ms, color=c, edgecolor='white', lw=0.6, zorder=3)
        if infl[i] > 0.35 * infl.max():   # label only the influential params
            ax.text(px[i], py[i] + sc*0.03, p, fontsize=6.6, color=INK, ha='center', va='bottom', zorder=4)
    ax.set_xlim(-sc*1.25, sc*1.25); ax.set_ylim(-sc*1.25, sc*1.25)
    ax.set_aspect('equal'); ax.tick_params(labelsize=7, colors=MUT)
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
    for sp in ('left', 'bottom'): ax.spines[sp].set_color('#ccc')

AXNAME = {0: 'replication & cycle timing', 1: 'EZH2 governor ↔ CyclinD1 level', 2: 'S-phase duration (fork vs firing)'}
fig, ax = plt.subplots(1, 2, figsize=(16, 8.2))
fig.subplots_adjust(left=0.05, right=0.98, top=0.86, bottom=0.09, wspace=0.18)
fig.suptitle('v44 parameter biplot — the model’s major axes of variation (SVD of the sensitivity matrix)',
             x=0.05, ha='left', fontsize=14.5, fontweight='bold', color=INK, y=0.975)
fig.text(0.05, 0.92, 'grey arrows = validation targets · dots = parameters (size = total influence, colour = module) · '
         'same direction = aligned, opposite = antagonistic', fontsize=9.3, color=MUT, ha='left')

biplot(ax[0], 0, 1)
ax[0].set_xlabel(f'axis 1  ({var[0]*100:.0f}%)  —  {AXNAME[0]}', fontsize=9.5, color=INK)
ax[0].set_ylabel(f'axis 2  ({var[1]*100:.0f}%)  —  {AXNAME[1]}', fontsize=9.5, color=INK)
biplot(ax[1], 1, 2)
ax[1].set_xlabel(f'axis 2  ({var[1]*100:.0f}%)  —  {AXNAME[1]}', fontsize=9.5, color=INK)
ax[1].set_ylabel(f'axis 3  ({var[2]*100:.0f}%)  —  {AXNAME[2]}', fontsize=9.5, color=INK)

# module legend
handles = [plt.Line2D([0],[0], marker='o', ls='', mec='white', mfc=c, ms=8, label=m) for m, c in MODCOL.items()]
fig.legend(handles=handles, loc='lower center', ncol=5, fontsize=7.8, frameon=False, bbox_to_anchor=(0.5, -0.01))
out = os.path.join(ROOT, 'simulations', 'fig_param_biplot')
fig.savefig(out + '.pdf', bbox_inches='tight'); fig.savefig(out + '.png', dpi=150, bbox_inches='tight')
print('wrote', out, '| var:', ' '.join(f'{v*100:.0f}%' for v in var[:4]))
