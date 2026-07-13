"""Does the current model support a 'collapse' to differentiation, and what role does EZH2 play?

Tests two scenarios in a cycling GNP:
  (A) Mitogen withdrawal (SHH 0.5 -> 0 at t0): the physiological exit trigger. Watch CyclinD1,
      divisions (MPF), EZH2, p27 — does the cell collapse to G0, and does EZH2 drive or follow it?
  (B) EZH2 knockdown alone (kTlEZ -> low at t0, mitogen intact): tests the user's proposed
      'reduced EZH2 -> CyclinD1 drop -> collapse'. In the current model the EZH2->CyclinD1 arm is
      REPRESSION, so this should do the OPPOSITE (CyclinD1 UP, no collapse).

Key question: is EZH2's fall a CAUSE of exit (positive-feedback collapse) or a CONSEQUENCE of it
(EZH2 is E2F-driven), and does the EZH2->CyclinD1 negative feedback resist the collapse?

Run:  ./venv/bin/python simulations/fig_v44_differentiation_collapse.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt, tellurium as te
from src.build_model_v44_heldt import build_model_v44
from simulations.validate_v44 import SEL

_M = build_model_v44()
t0 = 6000


def two_phase(perturb):
    rr = te.loada(_M); rr.integrator.setValue("absolute_tolerance", 1e-9); rr.integrator.setValue("relative_tolerance", 1e-6)
    try:
        rr.integrator.setValue("maximum_num_steps", 2000000); rr.integrator.setValue("maximum_time_step", 5.0)
    except Exception:
        pass
    rr['SHH'] = 0.5
    def _sim(a, b, n):   # small maxstep from a CLEAN state -- the two-step is stiff at the abrupt
        for atol in (1e-8, 1e-7, 1e-6):        # mid-cycle withdrawal discontinuity (maxstep 20 corrupts the CVODE state)
            for ms in (5.0, 2.0, 1.0):
                try:
                    rr.integrator.setValue("absolute_tolerance", atol)
                    rr.integrator.setValue("maximum_time_step", ms)
                    return rr.simulate(a, b, n, selections=SEL)
                except Exception:
                    continue
        return rr.simulate(a, b, n, selections=SEL)
    r1 = _sim(0, t0, 4000)
    perturb(rr)                                   # apply the perturbation at t0
    r2 = _sim(t0, 18000, 8000)                    # coarser output: dense sampling near the two-step
                                                  # withdrawal discontinuity trips CVODE convergence
    out = {}
    for k in ('time', 'Cd_mRNA', 'Cd', 'EZH2', 'MPF', 'P21', 'E2f'):
        out[k] = np.concatenate([r1[k], r2[k]])
    return out


A = two_phase(lambda rr: rr.__setitem__('SHH', 0.0))           # (A) mitogen withdrawal
B = two_phase(lambda rr: rr.__setitem__('kTlEZ', 0.004 * 0.1)) # (B) EZH2 knockdown alone

def summ(d, name):
    t = d['time']; pre = (t >= 3000) & (t < t0); post = t >= (t0 + 4000)
    div_pre = int(np.sum((d['MPF'][:-1] > 0.3) & (d['MPF'][1:] < 0.1) & (t[:-1] >= 3000) & (t[:-1] < t0)))
    div_post = int(np.sum((d['MPF'][:-1] > 0.3) & (d['MPF'][1:] < 0.1) & (t[:-1] >= t0)))
    print(f"{name}: CyclinD1 {d['Cd_mRNA'][pre].mean():.2f}->{d['Cd_mRNA'][post].mean():.2f}  "
          f"EZH2 {d['EZH2'][pre].mean():.2f}->{d['EZH2'][post].mean():.2f}  "
          f"p27 {d['P21'][pre].mean():.2f}->{d['P21'][post].mean():.2f}  divisions {div_pre}(pre)->{div_post}(post)")

print("=" * 78)
summ(A, "(A) mitogen withdrawal ")
summ(B, "(B) EZH2 knockdown only")

fig, ax = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
for axi, d, title in [(ax[0], A, '(A) Mitogen withdrawal (SHH->0): cell COLLAPSES to G0'),
                      (ax[1], B, '(B) EZH2 knockdown only (mitogen intact): CyclinD1 UP, NO collapse')]:
    th = d['time'] / 60.0
    axi.axvline(t0/60.0, color='k', ls='--', alpha=0.6)
    axi.plot(th, d['Cd_mRNA'], color='#117a65', lw=1.6, label='CyclinD1 transcript')
    axi.plot(th, d['EZH2'], color='#8e44ad', lw=1.6, label='EZH2 protein')
    axi.plot(th, d['P21'], color='#b9770e', lw=1.4, label='p27')
    axi.plot(th, d['MPF'] * 5, color='#e67e22', lw=0.8, alpha=0.8, label='divisions (MPF x5)')
    axi.set_title(title, fontweight='bold', fontsize=10.5); axi.set_xlabel('time (h)')
    axi.legend(fontsize=8, loc='upper right'); axi.grid(alpha=0.2)
ax[0].set_ylabel('level (a.u.)')
fig.suptitle('Differentiation as a collapse? EZH2 FOLLOWS exit (it does not drive it); the EZH2->CyclinD1 arm is negative feedback',
             fontsize=12, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/fig_v44_differentiation_collapse.png', dpi=170, bbox_inches='tight')
plt.savefig('simulations/fig_v44_differentiation_collapse.pdf', bbox_inches='tight')
plt.close()
print("\nSaved: fig_v44_differentiation_collapse.png")
