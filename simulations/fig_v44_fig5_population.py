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
T_END, N_PTS, SETTLE = 12000, 48000, 4000
P16_MB = 0.15     # MB CDK-inhibitor tone: INK4 (CDK4/6) = p16 0.88 + p18 1.2 (GNP p18 baseline 0.4) ...
P18_MB = 1.5
KSYP21_MB = 0.004 # ... + CIP/KIP (CDK2) = p21/p27 (kSyP21 2x the GNP baseline 0.002). All elevated in MB.
rng = np.random.default_rng(7)
# CyclinD1 heterogeneity (sigma 0.68, search-tuned): vismo lowers MB CyclinD1 onto the p16-raised
# commitment threshold, so the tumour straddles it -> fractional arrest; EZH2i shifts most cells back.
cd_scale = np.clip(np.exp(rng.normal(0.0, 0.68, N)), 0.25, 4.0)        # CyclinD1 setpoint heterogeneity
p21_div = np.clip(np.exp(rng.normal(np.log(0.85), 0.15, N)), 0.5, 2.0)  # inherited p27 (tight; the
# proliferation-quiescence split is carried by the CyclinD1/p16 axis above, not inherited p27)

_RR = te.loada(build_model_v44())                                     # de-saturated defaults
_RR.integrator.setValue("absolute_tolerance", 1e-9); _RR.integrator.setValue("relative_tolerance", 1e-6)
KTL0 = 0.75   # matches builder default k_Cd_translation (search-tuned)

CONDS = {  # MB context (Ptch1=0.1, MYCN_amp=2.8, p16 high) + treatment
    'MB':                 dict(hhi=0.0, ezh2i=0, cdk46i=False),
    'MB+HHi':             dict(hhi=1.0, ezh2i=0, cdk46i=False),
    'MB+HHi+EZH2i':       dict(hhi=1.0, ezh2i=1, cdk46i=False),
    'MB+CDK4/6i':         dict(hhi=0.0, ezh2i=0, cdk46i=True),
    'MB+CDK4/6i+EZH2i':   dict(hhi=0.0, ezh2i=1, cdk46i=True),
}


def cycles(ktl, p21d, hhi, ezh2i, cdk46i):
    """True/False cycling, or None if integration fails at every tolerance (excluded from denominator)."""
    for atol in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
        _RR.reset()
        _RR['SHH'] = 0.5; _RR['Ptch1_copy_number'] = 0.1; _RR['MYCN_amplification'] = 2.8
        _RR['HHi'] = hhi; _RR['EZH2i'] = ezh2i; _RR['p16'] = P16_MB; _RR['p18'] = P18_MB; _RR['kSyP21'] = KSYP21_MB
        _RR['P21_div'] = p21d; _RR['k_Cd_translation'] = ktl
        if cdk46i:
            _RR['kPhRbCd'] = 0.0
        _RR.integrator.setValue("absolute_tolerance", atol)
        try:
            _RR.integrator.setValue("maximum_num_steps", 200000)   # high-CyclinD1 EZH2i cells are stiff
        except Exception:
            pass
        try:
            r = _RR.simulate(0, T_END, N_PTS, selections=["time", "MPF"])
        except Exception:
            continue
        t = r['time']; m = t >= SETTLE; tt = t[m]; dt = tt[1] - tt[0]
        pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=int(200 / dt))
        return len(pk) >= 2
    return None     # integration crash -> indeterminate (NOT counted as arrested)


frac = {}
for name, c in CONDS.items():
    res = [cycles(KTL0 * cd_scale[i], p21_div[i], **c) for i in range(N)]
    ok = [x for x in res if x is not None]      # exclude integration crashes from the denominator
    frac[name] = 100 * sum(ok) / len(ok) if ok else 0.0
    print(f"{name:18s}: {frac[name]:.0f}% cycling  (n={len(ok)}/{N}, {N-len(ok)} crashed)")

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
