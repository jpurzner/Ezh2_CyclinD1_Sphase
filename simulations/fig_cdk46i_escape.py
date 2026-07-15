"""Figure: CDK4/6i therapy-escape -- GNP arrests cleanly, MB has TWO escape routes.

Reads cdk46i_escape_results.json (sim_cdk46i_escape.py) and renders the composition of each cohort under baseline
vs CDK4/6i as stacked bars:
    cycling (proliferating)  |  reversible G0 reserve (high-CKI, SOX2/OLIG2-like)  |
    CDK4/6i-arrested G0      |  resistant-dividing (low-Rb / RB1-loss)
The MB 'therapy-escape reservoir' = the reversible reserve (survives, re-enters) + the resistant fraction (never
stops dividing). GNP has neither -> it arrests cleanly.

Run:  ./venv/bin/python simulations/fig_cdk46i_escape.py
"""
import os, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, 'cdk46i_escape_results.json')))

# segment colours (semantic, colour-blind-safe-ish)
C_CYCLE   = '#2f7d54'    # green  -- proliferating
C_RESERVE = '#e0a13c'    # amber  -- reversible G0 reserve (survives, re-enters)
C_ARREST  = '#9aa4ad'    # grey   -- CDK4/6i-induced arrest (drug working)
C_RESIST  = '#c0392b'    # red    -- resistant, still dividing (drug failing)

def segments(d):
    """(baseline segments, cdk4/6i segments) as ordered (label, value, colour) from the bottom up."""
    g0b, g0i, divi = d['g0_base'], d['g0_cdki'], d['div_cdki']
    reserve = g0b                       # baseline G0 = the reversible reserve
    arrested = max(0.0, g0i - g0b)      # extra G0 induced by CDK4/6i = normal cells arresting
    base = [('Cycling', 100.0 - reserve, C_CYCLE), ('Reversible G0 reserve', reserve, C_RESERVE)]
    cdki = [('Resistant (dividing)', divi, C_RESIST),
            ('CDK4/6i-arrested', arrested, C_ARREST),
            ('Reversible G0 reserve', reserve, C_RESERVE)]
    return base, cdki

fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.9), sharey=True)
fig.subplots_adjust(left=0.075, right=0.70, top=0.80, bottom=0.16, wspace=0.10)

for ax, key in zip(axes, ('GNP', 'MB')):
    d = R[key]
    base, cdki = segments(d)
    for x, segs in [(0, base), (1, cdki)]:
        bottom = 0.0
        for label, val, col in segs:
            if val <= 0.05:
                bottom += val; continue
            ax.bar(x, val, bottom=bottom, width=0.62, color=col, edgecolor='white', linewidth=1.2)
            if val >= 6:
                ax.text(x, bottom + val / 2, f'{val:.0f}%', ha='center', va='center',
                        color='white', fontsize=11, fontweight='bold')
            bottom += val
    ax.set_xticks([0, 1]); ax.set_xticklabels(['Baseline', 'CDK4/6i'], fontsize=12)
    ax.set_xlim(-0.6, 1.6); ax.set_ylim(0, 100)
    ax.set_title(key, fontsize=15, fontweight='bold', pad=8)
    for s in ('top', 'right'): ax.spines[s].set_visible(False)

mb = R['MB']; esc = mb['div_cdki'] + mb['g0_base']
axes[0].set_ylabel('% of proliferation-competent pool', fontsize=12)
fig.suptitle('CDK4/6 inhibition: GNP arrests cleanly, medulloblastoma escapes two ways',
             fontsize=15.5, fontweight='bold', x=0.39, y=0.965)
fig.text(0.39, 0.895, f'MB keeps a ~{esc:.0f}% escape reservoir (20% reversible reserve + 15% resistant); GNP ~2%',
         fontsize=10.5, color='#444', ha='center')

legend = [Patch(fc=C_CYCLE, label='Cycling (proliferating)'),
          Patch(fc=C_RESERVE, label='Reversible G0 reserve\n(high-CKI, SOX2/OLIG2-like)'),
          Patch(fc=C_ARREST, label='CDK4/6i-arrested G0\n(drug working)'),
          Patch(fc=C_RESIST, label='Resistant, still dividing\n(low-Rb / RB1-loss)')]
fig.legend(handles=legend, loc='center left', bbox_to_anchor=(0.715, 0.52), fontsize=9.5,
           frameon=False, title='cell state', title_fontsize=11, alignment='left', labelspacing=1.1, handleheight=1.6)

foot = ('Reversible G0 reserve = a high-CDK-inhibitor subpopulation (baseline G0, survives and re-enters). '
        'Resistant fraction = low total Rb (proliferates independent of cyclinD-CDK4/6). GNP carries neither.')
fig.text(0.075, 0.04, foot, fontsize=8.2, color='#555')

out = os.path.join(HERE, 'fig_cdk46i_escape.png')
fig.savefig(out, dpi=150)
print(f"  GNP: baseline G0 {R['GNP']['g0_base']:.0f}%  CDK4/6i still-dividing {R['GNP']['div_cdki']:.0f}%")
print(f"  MB : baseline G0 {R['MB']['g0_base']:.0f}%  CDK4/6i still-dividing {R['MB']['div_cdki']:.0f}%")
print(f"  -> {out}")
