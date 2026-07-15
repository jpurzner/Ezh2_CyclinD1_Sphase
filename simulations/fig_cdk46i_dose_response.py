"""Figure: CDK4/6i DOSE-RESPONSE in the v44 MB population -- the pRb-positive residual is dose-limited, not a
resistant clone. Data from probe_dose_response.py (N=100 MB cells, CyclinD1 CV0.70 + CKI CV0.42). Every dividing
cell under sub-saturating inhibition is pRb-positive (Ser807/811+). Run: ./venv/bin/python simulations/fig_cdk46i_dose_response.py"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# (residual CDK4/6 activity, % dividing == % dividing pRb+); dividing fraction IS the pRb+ fraction at every dose
resid = np.array([1.00,0.90,0.80,0.75,0.72,0.70,0.68,0.66,0.64,0.62,0.60,0.50,0.30])
div   = np.array([98.0,95.0,93.0,96.0,96.0,96.0,92.0,85.0,74.0,45.0, 1.0, 0.0, 0.0])
inhib = 100*(1-resid)   # % CDK4/6 inhibition = dose axis (0 = no drug, 100 = full inhibition)

fig, ax = plt.subplots(figsize=(10.4,6.0))
fig.subplots_adjust(left=0.10, right=0.97, top=0.87, bottom=0.19)
ax.plot(inhib, div, '-o', color='#c0392b', lw=2.4, ms=6, zorder=3, label='% cells dividing  (ALL pRb-Ser807/811 positive)')
ax.fill_between(inhib, 0, div, color='#c0392b', alpha=0.08, zorder=1)

# mark JP's ~16-18% residual band and the dose that produces it
ax.axhspan(16, 18, color='#2f7d54', alpha=0.18, zorder=0)
ax.text(2, 17, "JP's residual  16–18%", color='#1e5638', fontsize=10, va='center', fontweight='bold')
# the sub-saturating dose window where the curve passes through 16-18%
ax.axvspan(100*(1-0.615), 100*(1-0.605), color='#7f8c8d', alpha=0.15, zorder=0)
ax.annotate('sub-saturating dose\n(~39% CDK4/6 inhibition)\ngives ~16–18% dividing, pRb+',
            xy=(38.5, 17), xytext=(52, 40), fontsize=9.5, color='#333', ha='left',
            arrowprops=dict(arrowstyle='->', color='#555', lw=1.2))

ax.set_xlabel('CDK4/6 inhibition  (dose →)', fontsize=12)
ax.set_ylabel('% of MB cells dividing', fontsize=12)
ax.set_xlim(-2, 72); ax.set_ylim(-3, 103)
for s in ('top','right'): ax.spines[s].set_visible(False)
ax.set_title('CDK4/6i dose-response: the pRb⁺ residual is dose-limited, not a resistant clone',
             fontsize=12.5, fontweight='bold', pad=12)
ax.legend(loc='upper right', fontsize=10, frameon=False)
foot = ('v44 MB population (CyclinD1 CV 0.70 + CKI CV 0.42). Every dividing cell is pRb-Ser807/811⁺ — the residual '
        'is the high-CyclinD tail beating partial inhibition\n(gone at saturating dose, back on washout). NB the '
        'modeled transition is steeper than a real CDK4/6i curve; wider heterogeneity would spread it realistically.')
fig.text(0.10, 0.015, foot, fontsize=7.6, color='#666', va='bottom', linespacing=1.5)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_cdk46i_dose_response.png')
fig.savefig(out, dpi=150); print('->', out)
