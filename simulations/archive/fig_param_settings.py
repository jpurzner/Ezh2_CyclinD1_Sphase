"""Parameter-settings figure (JP 2026-07-27): for each load-bearing parameter show (1) PROVENANCE (measured/
inherited/soft/free), (2) the tested RANGE and the LANDED value, (3) GNP|MB for the cell-type knobs, and (4) a
SENSITIVITY strip = elasticity of each of the 29 validation targets to that parameter. Reads scratchpad
sensitivity.json + param_meta.json. Writes fig_param_settings.pdf + .png."""
import os, json, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SC = '/private/tmp/claude-501/-Users-jpurzner-Dropbox-Q-research-py-projects-ezh2-cyclind1-scheckpoint/81291f88-7d19-4149-be27-53b1162ee52d/scratchpad'
INK, MUT, FAINT, GRID = '#2d2d2d', '#666', '#9a9a9a', '#e2e2e2'
C_GNP, C_MB = '#2c7fb8', '#c0392b'
PROVC = {'MEASURED': '#2a9d5c', 'INHERITED': '#8a8a8a', 'SOFT': '#e0952b', 'FREE': '#3f7fd0'}
MODCOL = {'Hedgehog': '#1b7837', 'MYCN': '#762a83', 'CyclinD1': '#2166ac', 'EZH2': '#7d5ba6',
          'H3K27me3': '#8c6bb1', 'Rb-E2F': '#b2182b', 'Skp2-p27/CKI': '#d6604d', 'CDK-cyclins': '#e08214',
          'DNA-replication': '#35978f', 'Mitosis-APC': '#4393c3', 'HU': '#878787', 'Growth-size': '#5aae61'}
plt.rcParams.update({'font.size': 8, 'font.family': 'DejaVu Sans', 'pdf.fonttype': 42})

S = json.load(open(SC + '/sensitivity.json'))
META = json.load(open(SC + '/param_meta.json'))
E, LANDED, TGTS = S['elasticity'], S['landed'], S['targets']

# short target labels (29)
def abbr(t):
    t = t.split('(')[0].strip()
    r = {'CyclinD1':'CycD1','reduction':'red','arrest':'arr','cycles':'cyc','count%':'%','duration':'dur',
         'transcript':'mRNA','protein':'prot','GNP+HHi':'G+HHi','MB+HHi':'M+HHi','Serum':'ss'}
    for k,v in r.items(): t=t.replace(k,v)
    return t[:15]
TL = [abbr(t) for t in TGTS]

MOD_ORDER = ['Hedgehog','MYCN','CyclinD1','EZH2','H3K27me3','Rb-E2F','Skp2-p27/CKI','DNA-replication','HU','Growth-size']
KNOBS = ['Ptch1_copy_number','MYCN_expr','p16','p18','P21_div']  # cell-type identities (no global sensitivity row)
# shared load-bearing params (have sensitivity), grouped by module, sorted by influence
shared = [p for p in E if p in META and p not in KNOBS]
infl = {p: sum(abs(v) for v in E[p].values()) for p in shared}
groups = []
for mod in MOD_ORDER:
    ps = sorted([p for p in shared if META[p]['module']==mod], key=lambda p:-infl[p])
    if ps: groups.append((mod, ps))
misc = [p for p in shared if META[p]['module'] not in MOD_ORDER]
if misc: groups.append(('other', sorted(misc, key=lambda p:-infl[p])))

nrows = len(KNOBS) + 2 + sum(1+len(ps) for _,ps in groups)   # knobs + header + module(header+rows)
fig_h = nrows*0.30 + 2.2
fig = plt.figure(figsize=(19, fig_h))
axL = fig.add_axes([0.005, 0.02, 0.42, 0.855]); axL.axis('off'); axL.set_xlim(0,1); axL.set_ylim(0,nrows); axL.invert_yaxis()
axH = fig.add_axes([0.435, 0.02, 0.56, 0.855]); axH.set_xlim(0,len(TGTS)); axH.set_ylim(0,nrows); axH.invert_yaxis()

fig.suptitle('v44 parameter settings — provenance · tested range & landed value · sensitivity',
             x=0.005, ha='left', fontsize=15, fontweight='bold', color=INK, y=0.995)
fig.text(0.005, 0.975, '42 load-bearing parameters (the inherited Heldt-engine + chromatin constants are summarized, not shown) · '
         'sensitivity = elasticity of each target to a +25% change', fontsize=9, color=MUT, ha='left')

norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=2)
cmap = plt.cm.RdBu_r
# target column headers on the heatmap
for j,t in enumerate(TL):
    axH.text(j+0.5, -0.3, t, rotation=90, ha='center', va='bottom', fontsize=6.0, color=MUT)

def prov_chip(y, prov):
    axL.add_patch(plt.Rectangle((0.30, y+0.16), 0.085, 0.66, facecolor=PROVC[prov], edgecolor='none'))
    axL.text(0.3425, y+0.5, prov[:4], ha='center', va='center', fontsize=6.0, color='white', fontweight='bold')

def range_landed(y, p):
    m = META[p]; landed = LANDED[p]; rng = m['range']
    x0, x1 = 0.42, 0.995
    if rng:
        lo, hi = rng
        # log position of landed within [lo,hi]
        f = (math.log10(landed)-math.log10(lo))/(math.log10(hi)-math.log10(lo)+1e-12)
        f = min(max(f,0),1)
        axL.plot([x0,x1],[y+0.5,y+0.5], color=GRID, lw=3.2, solid_capstyle='round')
        axL.plot([x0,x0],[y+0.34,y+0.66],color=FAINT,lw=1); axL.plot([x1,x1],[y+0.34,y+0.66],color=FAINT,lw=1)
        xl = x0+f*(x1-x0)
        axL.plot(xl, y+0.5, 'o', ms=7, color=INK, mec='white', mew=1.1, zorder=3)
        axL.text(x0, y+0.5, f'{lo:g}', ha='right', va='center', fontsize=5.6, color=FAINT)
        axL.text(x1+0.002, y+0.5, f'{hi:g}', ha='left', va='center', fontsize=5.6, color=FAINT)
        axL.text(xl, y-0.02, f'{landed:.4g}', ha='center', va='bottom', fontsize=6.2, color=INK, fontweight='bold')
    else:
        axL.text(0.42, y+0.5, f'{landed:.4g}', ha='left', va='center', fontsize=7, color=INK)
        axL.text(0.55, y+0.5, 'prior opt', ha='left', va='center', fontsize=5.8, color=FAINT, style='italic')

y = 0
# ---- cell-type identities section ----
axL.add_patch(plt.Rectangle((0,y+0.1),1.0,0.8,facecolor='#efe9f2',edgecolor='none'))
axL.text(0.005,y+0.5,'CELL-TYPE IDENTITIES',fontsize=8.5,fontweight='bold',color='#5a2d82',va='center')
axL.text(0.62,y+0.5,'GNP',color=C_GNP,fontsize=8,fontweight='bold',ha='center',va='center')
axL.text(0.80,y+0.5,'MB',color=C_MB,fontsize=8,fontweight='bold',ha='center',va='center')
axH.text(len(TGTS)/2, y+0.5, '(set per condition — not in the global sensitivity sweep)', ha='center', va='center', fontsize=7, color=FAINT, style='italic')
y+=1
for p in KNOBS:
    m=META[p]; g,mb=m['celltype']
    axL.text(0.012,y+0.42,p,fontsize=7.4,color=INK,va='center'); axL.text(0.012,y+0.74,m['role'],fontsize=5.6,color=FAINT,va='center')
    prov_chip(y,m['prov'])
    axL.text(0.62,y+0.5,f'{g:g}',color=C_GNP,fontsize=8,ha='center',va='center')
    axL.text(0.80,y+0.5,f'{mb:g}',color=C_MB,fontsize=8,ha='center',va='center',fontweight='bold')
    axH.axhspan(y,y+1,color='#efe9f2',alpha=0.5)
    y+=1

# ---- shared load-bearing params with sensitivity ----
axL.text(0.005,y+0.5,'SHARED ENGINE — LOAD-BEARING (by module, sorted by influence)',fontsize=8.5,fontweight='bold',color=MUT,va='center'); y+=1
for mod, ps in groups:
    mc=MODCOL.get(mod,'#888')
    axL.add_patch(plt.Rectangle((0,y+0.15),1.0,0.7,facecolor=mc,alpha=0.13,edgecolor='none'))
    axL.text(0.008,y+0.5,mod,fontsize=8,fontweight='bold',color=mc,va='center')
    y+=1
    for p in ps:
        m=META[p]
        axL.text(0.012,y+0.42,p,fontsize=7.2,color=INK,va='center')
        axL.text(0.012,y+0.74,m['role'],fontsize=5.6,color=FAINT,va='center')
        prov_chip(y,m['prov'])
        range_landed(y,p)
        for j,t in enumerate(TGTS):
            e=E[p][t]
            if abs(e)<0.05: continue
            axH.add_patch(plt.Rectangle((j,y+0.12),1,0.76,facecolor=cmap(norm(e)),edgecolor='none'))
        y+=1
axH.set_xticks([]); axH.set_yticks([])
for s in axH.spines.values(): s.set_visible(False)
# gridlines between target columns (faint)
for j in range(len(TGTS)+1): axH.plot([j,j],[len(KNOBS)+2,nrows],color='white',lw=0.4)

# legends
lx=0.005
for i,(pv,c) in enumerate(PROVC.items()):
    fig.text(lx+i*0.055, 0.005, '  ', bbox=dict(boxstyle='square',fc=c,ec='none'), fontsize=6)
    fig.text(lx+i*0.055+0.014, 0.004, pv.lower(), fontsize=7.5, color=INK)
from matplotlib.cm import ScalarMappable
cax=fig.add_axes([0.46,0.006,0.16,0.010]); cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='horizontal')
cb.set_ticks([-2,0,2]); cb.set_ticklabels(['−2','0','+2']); cb.ax.tick_params(labelsize=6.5,length=2)
fig.text(0.63,0.004,'sensitivity: elasticity of target to +25% param (red=up, blue=down)',fontsize=7.5,color=MUT)

out=os.path.join(ROOT,'simulations','fig_param_settings')
fig.savefig(out+'.pdf',bbox_inches='tight'); fig.savefig(out+'.png',dpi=145,bbox_inches='tight')
print('wrote',out+'.pdf / .png')
