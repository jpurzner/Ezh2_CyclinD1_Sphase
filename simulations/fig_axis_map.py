"""Axis map (JP 2026-07-27): the v44 model's ~4 emergent behavioral axes, with each parameter placed on the
axis it moves and the side it pushes toward, plus the key couplings (CKI hits two axes; decouple_commit severed
the speed->threshold leak; the EZH2 governor sets the gain of A). Editorial schematic. Writes .pdf + .png."""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INK, MUT, FAINT = '#2d2d2d', '#5f5f5f', '#9a9a9a'
C_GNP, C_MB = '#2c7fb8', '#c0392b'
MODCOL = {'Hedgehog':'#1b7837','MYCN':'#762a83','CyclinD1':'#2166ac','EZH2':'#7d5ba6','H3K27me3':'#8c6bb1',
          'Rb':'#b2182b','CKI':'#d6604d','DNA-rep':'#35978f','HU':'#e08214','Growth':'#5aae61'}
plt.rcParams.update({'font.size': 9, 'font.family': 'DejaVu Sans', 'pdf.fonttype': 42})

AXES = [
 dict(y=3.4, name='A  ·  PROLIFERATION  ↔  QUIESCENCE',
      sub='the go / no-go decision  =  net (mitogen drive − commitment threshold)  =  your "G0 probability" and "mitogen-impact"',
      lo='QUIESCENT · G0 · arrest', hi='PROLIFERATIVE · divides',
      left=[('p18','CKI'),('p16','CKI'),('P21_div','CKI'),('K_CdRb','Rb'),('M_commit','Growth')],
      right=[('k_Cd_tx_Gli_max','CyclinD1'),('MYCN_expr','MYCN'),('kPhRbCd','Rb'),('Ptch1-loss→Gli','Hedgehog'),('k_Cd_tx_basal','CyclinD1')],
      gnp=0.72, mb=0.56),
 dict(y=2.4, name='B  ·  CYCLE SPEED  (Tc)',
      sub='how fast, when proliferating · Tc = ln2 / mu_eff',
      lo='SLOW · long Tc (MB ~23.5 h)', hi='FAST · short Tc (GNP ~16 h)',
      left=[('k_mu_cki','Growth'),('p18','CKI'),('p16','CKI'),('mu_min_frac↓','Growth')],
      right=[('mu','Growth')], gnp=0.82, mb=0.34),
 dict(y=1.4, name='C  ·  EZH2 GOVERNOR  (gain / buffering)',
      sub='how mitogen-dependent A is · EZH2 ⊣ CyclinD1, incoherent feed-forward + slow reservoir',
      lo='WEAK · mitogen-slaved', hi='STRONG · buffered / memory',
      left=[('K_prc2','EZH2'),('kDeEZ','EZH2'),('f0_prc2','EZH2')],
      right=[('kEZbas','EZH2'),('kEZE2f','EZH2'),('a_rw_prc2','H3K27me3')], gnp=0.44, mb=0.66),
 dict(y=0.4, name='D  ·  REPLICATION / S-PHASE',
      sub='S-phase length & fork stress · largely orthogonal to A–C (machinery)',
      lo='LONG S · HU-sensitive', hi='SHORT S · fast fork',
      left=[('KmHU_fire','HU'),('KmHU_fork','HU'),('vmin_fork','HU')],
      right=[('kSyDna','DNA-rep')], gnp=0.5, mb=0.5),
]
X0, X1 = 2.55, 9.35           # track span
def xf(f): return X0 + f*(X1-X0)

fig, ax = plt.subplots(figsize=(15.5, 9.6)); ax.set_xlim(0, 10); ax.set_ylim(-0.3, 4.35); ax.axis('off')
ax.text(0.1, 4.18, 'v44 model — the major axes of variation and how the parameters push them',
        fontsize=15, fontweight='bold', color=INK)
ax.text(0.1, 3.98, 'each parameter is placed on the axis it moves, on the side it pushes toward · '
        'a parameter on two tracks moves both · GNP▲ / MB▲ mark where each cell type sits',
        fontsize=9.3, color=MUT)

def stagger(n):
    off=[0.0,0.20,-0.20,0.38,-0.38,0.56,-0.56]; return off[:n]
for A in AXES:
    y=A['y']
    ax.annotate('', (X1+0.15,y),(X0-0.15,y), arrowprops=dict(arrowstyle='<|-|>',color='#7a7a7a',lw=1.8))
    ax.text(0.1, y+0.13, A['name'], fontsize=11.5, fontweight='bold', color=INK)
    ax.text(0.1, y-0.05, A['sub'], fontsize=7.4, color=FAINT)
    ax.text(X0-0.2, y-0.24, A['lo'], fontsize=8.2, color=MUT, ha='left', style='italic')
    ax.text(X1+0.2, y-0.24, A['hi'], fontsize=8.2, color=MUT, ha='right', style='italic')
    for side,items,base in (('L',A['left'],[0.34,0.20,0.06]),('R',A['right'],[0.66,0.80,0.94])):
        offs=stagger(len(items))
        for i,(p,mod) in enumerate(items):
            f = 0.34 - 0.10*(i//2) if side=='L' else 0.66 + 0.10*(i//2)
            f = max(0.06,min(0.94,f)); x=xf(f); dy=offs[i]
            c=MODCOL.get(mod,'#888')
            ax.plot([x,x],[y, y+dy*0.62],color=c,lw=0.8,alpha=0.5,zorder=1)
            ax.scatter(x,y+dy*0.62,s=26,color=c,edgecolor='white',lw=0.7,zorder=3)
            ax.text(x, y+dy*0.62+(0.055 if dy>=0 else -0.055), p, fontsize=6.9, color=INK,
                    ha='center', va='bottom' if dy>=0 else 'top', zorder=4)
    # GNP / MB position markers
    ax.scatter(xf(A['gnp']), y-0.16, marker='^', s=42, color=C_GNP, edgecolor='white', lw=0.6, zorder=5)
    ax.scatter(xf(A['mb']),  y-0.16, marker='^', s=42, color=C_MB,  edgecolor='white', lw=0.6, zorder=5)

# ---- couplings ----
# CKI hits A and B (p16/p18 on both, left side)
axA,axB=AXES[0]['y'],AXES[1]['y']
ax.add_patch(FancyArrowPatch((xf(0.24),axA-0.28),(xf(0.24),axB+0.28),connectionstyle='arc3,rad=0.25',
             arrowstyle='-',color=MODCOL['CKI'],lw=1.6,ls=(0,(4,2)),zorder=2))
ax.text(xf(0.05), (axA+axB)/2, 'CKI tone (p16/p18):\nthe master MB dial —\n↑ quiescent  &  ↑ slow', fontsize=7.6,
        color=MODCOL['CKI'], ha='left', va='center', fontweight='bold')
# decouple_commit severed the speed->threshold leak (between B and A)
ax.add_patch(FancyArrowPatch((xf(0.60),axB+0.30),(xf(0.60),axA-0.30),connectionstyle='arc3,rad=-0.2',
             arrowstyle='-',color='#bbb',lw=1.4,ls=(0,(2,2)),zorder=2))
ax.text(xf(0.605), (axA+axB)/2, '✂  decouple_commit\nsevered the speed→threshold leak\n(fast cycle no longer lowers the R-point)',
        fontsize=7.4, color=MUT, ha='left', va='center')
# governor C sets the gain of A
axC=AXES[2]['y']
ax.add_patch(FancyArrowPatch((xf(0.80),axC+0.30),(xf(0.80),axA-0.30),connectionstyle='arc3,rad=0.28',
             arrowstyle='-|>',color=MODCOL['EZH2'],lw=1.6,zorder=2))
ax.text(xf(0.86), (axA+axC)/2, 'the governor sets\nthe GAIN of A\n(how mitogen-\ndependent division is)', fontsize=7.4,
        color=MODCOL['EZH2'], ha='left', va='center')

# legend
ax.scatter(0.5,-0.15,marker='^',s=42,color=C_GNP,edgecolor='white',lw=0.6); ax.text(0.75,-0.17,'GNP',fontsize=8,color=C_GNP,va='center',fontweight='bold')
ax.scatter(1.7,-0.15,marker='^',s=42,color=C_MB,edgecolor='white',lw=0.6); ax.text(1.95,-0.17,'MB',fontsize=8,color=C_MB,va='center',fontweight='bold')
ax.text(3.0,-0.17,'dot colour = module',fontsize=7.6,color=MUT,va='center')

out=os.path.join(ROOT,'simulations','fig_axis_map')
fig.savefig(out+'.pdf',bbox_inches='tight'); fig.savefig(out+'.png',dpi=150,bbox_inches='tight')
print('wrote',out)
