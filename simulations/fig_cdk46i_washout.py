"""Figure: CDK4/6i WASHOUT kinetics in MB. 72h palbociclib -> washout, tracking the cycling fraction over time.
Reads cdk46i_washout_results.json (sim_cdk46i_washout.py). Reproduces JP's reversibility: 1 µM (sub-saturating)
holds a ~16-18% cycling residual and rebounds fast; 5 µM (saturating) arrests fully and re-enters more slowly --
both fully reversible. Run: ./venv/bin/python simulations/fig_cdk46i_washout.py"""
import os, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, 'cdk46i_washout_results.json')))
t = np.array(R['time_h']); palbo_h = R['palbo_h']
COL = {'1 µM (sub-saturating)': '#e0a13c', '5 µM (saturating)': '#c0392b'}

fig, ax = plt.subplots(figsize=(10.4, 5.8))
fig.subplots_adjust(left=0.09, right=0.97, top=0.88, bottom=0.20)

ax.axvspan(0, palbo_h, color='#9aa4ad', alpha=0.16, zorder=0)
ax.text(palbo_h / 2, 96, 'palbociclib on (72h)', ha='center', fontsize=10, color='#555')
ax.text(palbo_h + (t[-1] - palbo_h) / 2, 96, 'washout', ha='center', fontsize=10, color='#555')
ax.axvline(0, color='#555', ls='--', lw=1, alpha=0.5); ax.axvline(palbo_h, color='#555', ls='--', lw=1, alpha=0.5)

for label, frac in R['doses'].items():
    ax.plot(t, frac, '-', color=COL.get(label, '#333'), lw=2.6, label=label, zorder=3)
ax.axhspan(16, 18, color='#7f8c8d', alpha=0.12, zorder=0)
ax.text(t[0] + 2, 17, "JP: 1 µM residual 16–18%", fontsize=8.5, color='#555', va='center')

ax.set_xlabel('time (hours)   —   0 = palbociclib start, cell cycle ~22h', fontsize=12)
ax.set_ylabel('% of cells cycling  (divided within last cycle)', fontsize=12)
ax.set_xlim(-30, t[-1]); ax.set_ylim(-3, 103)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax.set_title('CDK4/6i washout kinetics: reversible arrest, dose-dependent re-entry',
             fontsize=13, fontweight='bold', pad=10)
ax.legend(loc='center right', fontsize=10.5, frameon=False, title='dose')

foot = ('MB population (CyclinD1 CV 0.70 + CKI CV 0.42 + 20% drug-independent reserve); palbo = graded residual '
        'CDK4/6 activity. 1 µM holds a partial cycling residual and rebounds within ~24h;\n5 µM arrests fully and '
        're-enters more slowly (deeper arrest → more p27 to clear) — both fully reversible, no senescence.')
fig.text(0.09, 0.02, foot, fontsize=7.8, color='#666', va='bottom', linespacing=1.5)
out = os.path.join(HERE, 'fig_cdk46i_washout.png')
fig.savefig(out, dpi=150)
fig.savefig(out.replace('.png', '.pdf'))
for label, frac in R['doses'].items():
    f = np.array(frac)
    print(f"  {label}: baseline {np.mean(f[t<0]):.0f}% -> palbo plateau {np.mean(f[(t>palbo_h*0.5)&(t<palbo_h)]):.0f}% -> +24h {np.interp(palbo_h+24,t,f):.0f}%")
print('->', out)
