"""Population EZH2i dose-response: what FRACTION of a heterogeneous tumour is pushed back into the
cycle at each EZH2i dose. Combines the cyclin D1/p27 heterogeneity (sim_g0_bifurcation) with a graded
EZH2i sweep in the MB+HHi (vismodegib-arrested) context and the MB+CDK4/6i control.

This is the clinically relevant curve behind the paper's caution that a sub-maximal EZH2i dose may
rescue a FRACTION of arrested cells (increasing proliferation) without the full benefit.
Run (background, ~15 min):  ./venv/bin/python simulations/sim_ezh2i_population_dose.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te, matplotlib.pyplot as plt
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

N = 100
T_END, N_PTS, SETTLE = 11000, 22000, 5500
DOSES = np.round(np.linspace(0, 1, 11), 2)
rng = np.random.default_rng(7)
cd_scale = np.clip(np.exp(rng.normal(0.0, 0.30, N)), 0.55, 1.9)
p21_div = np.clip(np.exp(rng.normal(np.log(0.85), 0.32, N)), 0.5, 2.0)

_RR = te.loada(build_model_v44())
_RR.integrator.setValue("absolute_tolerance", 1e-9); _RR.integrator.setValue("relative_tolerance", 1e-6)


def cycles(ktl, p21d, ezh2i, cdk46i):
    _RR.reset()
    _RR['SHH'] = 0.5; _RR['Ptch1_copy_number'] = 0.1; _RR['MYCN_amplification'] = 2.8
    _RR['GDC0449'] = 0.0 if cdk46i else 1.0           # CDK4/6i context has no HHi; HHi context has GDC=1
    _RR['EZH2i'] = ezh2i; _RR['P21_div'] = p21d; _RR['k_Cd_translation'] = ktl
    if cdk46i:
        _RR['kPhRbCd'] = 0.0
    try:
        r = _RR.simulate(0, T_END, N_PTS, selections=["time", "MPF"])
    except Exception:
        return False
    t = r['time']; m = t >= SETTLE; tt = t[m]; dt = tt[1] - tt[0]
    pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=int(200 / dt))
    return len(pk) >= 2


frac = {'MB + HHi': [], 'MB + CDK4/6i': []}
for ctx, cdk in [('MB + HHi', False), ('MB + CDK4/6i', True)]:
    for d in DOSES:
        nc = sum(cycles(0.26 * cd_scale[i], p21_div[i], float(d), cdk) for i in range(N))
        frac[ctx].append(100 * nc / N)
        print(f"{ctx:13s} EZH2i={d:.2f}: {100*nc/N:.0f}% cycling")

fig, ax = plt.subplots(figsize=(7.5, 5))
for ctx, c in [('MB + HHi', '#E08214'), ('MB + CDK4/6i', '#C2185B')]:
    ax.plot(DOSES, frac[ctx], color=c, lw=2.5, marker='o', label=ctx)
ax.set_xlabel('EZH2i dose (0 = none, 1 = full)'); ax.set_ylabel('% of tumour cells cycling (rescued)')
ax.set_title('Population EZH2i dose-response in arrested MB\n(fraction of a heterogeneous tumour pushed back into the cycle)',
             fontweight='bold', fontsize=11)
ax.legend(fontsize=9); ax.grid(alpha=0.2); ax.set_ylim(-3, 103)
plt.tight_layout()
plt.savefig('simulations/fig_v44_ezh2i_population_dose.png', dpi=170, bbox_inches='tight')
plt.savefig('simulations/fig_v44_ezh2i_population_dose.pdf', bbox_inches='tight')
plt.close()
print("\nSaved: fig_v44_ezh2i_population_dose.png")
