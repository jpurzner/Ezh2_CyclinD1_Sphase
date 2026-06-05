"""Figure 5 (population reframing): vismodegib reduces, and EZH2i restores, the CYCLING FRACTION of a
heterogeneous Ptch+/- MB tumour. With the data-matched (de-saturated) HH, MB+HHi retains ~3x a cycling
GNP's CyclinD1, so single cells still cycle; the vismodegib effect is a population/fractional one (a
fraction of cells drop below the commitment threshold) -- matching the fractional Fig 5D/E bar charts.

Ensemble of MB cells with mother-set cyclin D1 / p27 heterogeneity (as sim_g0_bifurcation). % cycling:
  MB  ->  MB+HHi (reduced)  ->  MB+HHi+EZH2i (restored)        [the rescue]
  MB+CDK4/6i (reduced)  ->  MB+CDK4/6i+EZH2i (NOT restored)    [downstream block, no rescue]

Run (background):  ./venv/bin/python simulations/fig_v44_fig5_population.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te, matplotlib.pyplot as plt
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

N = 120
T_END, N_PTS, SETTLE = 11000, 22000, 5500
rng = np.random.default_rng(7)
cd_scale = np.clip(np.exp(rng.normal(0.0, 0.32, N)), 0.5, 2.0)        # CyclinD1 setpoint heterogeneity
p21_div = np.clip(np.exp(rng.normal(np.log(0.85), 0.33, N)), 0.5, 2.0)  # inherited p27

_RR = te.loada(build_model_v44())                                     # de-saturated defaults
_RR.integrator.setValue("absolute_tolerance", 1e-9); _RR.integrator.setValue("relative_tolerance", 1e-6)
KTL0 = 0.8   # matches builder default k_Cd_translation (recalibrated to MB_GDC0449 vismo drop)

CONDS = {  # MB context (Ptch1=0.1, MYCN_amp=2.8) + treatment
    'MB':                 dict(gdc=0.0, ezh2i=0, cdk46i=False),
    'MB+HHi':             dict(gdc=1.0, ezh2i=0, cdk46i=False),
    'MB+HHi+EZH2i':       dict(gdc=1.0, ezh2i=1, cdk46i=False),
    'MB+CDK4/6i':         dict(gdc=0.0, ezh2i=0, cdk46i=True),
    'MB+CDK4/6i+EZH2i':   dict(gdc=0.0, ezh2i=1, cdk46i=True),
}


def cycles(ktl, p21d, gdc, ezh2i, cdk46i):
    for atol in (1e-9, 1e-8):
        _RR.reset()
        _RR['SHH'] = 0.5; _RR['Ptch1_copy_number'] = 0.1; _RR['MYCN_amplification'] = 2.8
        _RR['GDC0449'] = gdc; _RR['EZH2i'] = ezh2i
        _RR['P21_div'] = p21d; _RR['k_Cd_translation'] = ktl
        if cdk46i:
            _RR['kPhRbCd'] = 0.0
        _RR.integrator.setValue("absolute_tolerance", atol)
        try:
            r = _RR.simulate(0, T_END, N_PTS, selections=["time", "MPF"])
        except Exception:
            continue
        t = r['time']; m = t >= SETTLE; tt = t[m]; dt = tt[1] - tt[0]
        pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=int(200 / dt))
        return len(pk) >= 2
    return False


frac = {}
for name, c in CONDS.items():
    nc = sum(cycles(KTL0 * cd_scale[i], p21_div[i], **c) for i in range(N))
    frac[name] = 100 * nc / N
    print(f"{name:18s}: {frac[name]:.0f}% cycling")

fig, ax = plt.subplots(figsize=(8.5, 5))
order = list(CONDS); vals = [frac[n] for n in order]
colors = ['#762A83', '#E08214', '#4A1486', '#C2185B', '#D81B60']
b = ax.bar(range(len(order)), vals, color=colors, edgecolor='k')
for bar, v in zip(b, vals):
    ax.text(bar.get_x() + bar.get_width()/2, v + 1.5, f'{v:.0f}%', ha='center', fontweight='bold', fontsize=10)
ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=20, ha='right', fontsize=9)
ax.set_ylabel('% of tumour cells cycling'); ax.set_ylim(0, 105)
ax.set_title('Figure 5 (population): EZH2i rescues the cycling fraction of HHi-treated MB,\n'
             'but cannot rescue CDK4/6i (downstream block) — matches the fractional Fig 5D/E',
             fontweight='bold', fontsize=10.5)
# bracket the rescue
ax.annotate('', xy=(2, frac['MB+HHi+EZH2i']+6), xytext=(1, frac['MB+HHi']+6),
            arrowprops=dict(arrowstyle='->', color='#4A1486', lw=2))
ax.text(1.5, max(frac['MB+HHi'], frac['MB+HHi+EZH2i'])+9, 'rescue', ha='center', color='#4A1486', fontsize=9, fontweight='bold')
plt.tight_layout()
plt.savefig('simulations/fig_v44_fig5_population.png', dpi=170, bbox_inches='tight')
plt.savefig('simulations/fig_v44_fig5_population.pdf', bbox_inches='tight')
plt.close()
print("\nSaved: fig_v44_fig5_population.png")
