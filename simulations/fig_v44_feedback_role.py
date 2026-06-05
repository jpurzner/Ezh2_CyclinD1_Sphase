"""What does the EZH2 -| CyclinD1 arm DO? (present in GNP/muscle, absent in most cells)

The CyclinD1 -> Rb-E2F -> EZH2 arm is ubiquitous; the EZH2 -| CyclinD1 arm is context-specific and
CLOSES a negative feedback loop. We compare the cell WITH this arm (default) vs WITHOUT it
(K_EZH2_repression -> 1e6, i.e. EZH2 present but cannot repress CyclinD1, as in 'most cells') across
mitogen (SHH) levels, to see the functional consequences:

  - homeostasis: does the feedback BUFFER CyclinD1 (and proliferation) against mitogen level?
  - what the loop does to CyclinD1 level, division rate, and the transient G0.

Run:  ./venv/bin/python simulations/fig_v44_feedback_role.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44
from simulations.validate_v44 import count_divisions, mean_settled, classify, P27_THR, SEL

T_END, N_PTS = 10080, 20160
SHH = np.round(np.linspace(0.2, 1.0, 9), 2)

M_fb = build_model_v44()                                        # feedback ON (default)
M_no = build_model_v44(params={"K_EZH2_repression": 1e6})       # feedback OFF (EZH2 cannot repress Cd)


def run(model, shh):
    for atol in (1e-9, 1e-8, 1e-7):
        rr = te.loada(model); rr.integrator.setValue("absolute_tolerance", atol)
        rr.integrator.setValue("relative_tolerance", 1e-6)
        rr['SHH'] = shh; rr['MYCN_amplification'] = 1.0; rr['Ptch1_copy_number'] = 1.0
        try:
            return rr.simulate(0, T_END, N_PTS, selections=SEL)
        except Exception:
            continue
    return None


res = {'fb': {k: [] for k in ('cdm', 'cd', 'div', 'g0', 'ezh2')},
       'no': {k: [] for k in ('cdm', 'cd', 'div', 'g0', 'ezh2')}}
for tag, model in (('fb', M_fb), ('no', M_no)):
    for s in SHH:
        r = run(model, s)
        if r is None:
            for k in res[tag]: res[tag][k].append(np.nan)
            continue
        n, _, _ = count_divisions(r); f, _ = classify(r, P27_THR)
        res[tag]['cdm'].append(mean_settled(r, 'Cd_mRNA'))
        res[tag]['cd'].append(mean_settled(r, 'Cd'))
        res[tag]['div'].append(n); res[tag]['g0'].append(f['G0'])
        res[tag]['ezh2'].append(mean_settled(r, 'EZH2'))

# buffering: CyclinD1 fold change from lowest to highest SHH, with vs without feedback
def fold(tag): a = res[tag]['cdm']; return a[-1] / a[0]
print("=" * 70)
print("EZH2 -| CyclinD1 feedback: GNP, mitogen (SHH) sweep")
print(f"{'SHH':>5} | {'CyclinD1 (fb)':>13} {'CyclinD1 (no-fb)':>16} | {'div fb':>6} {'div no':>6} | {'G0 fb':>6} {'G0 no':>6}")
for i, s in enumerate(SHH):
    print(f"{s:>5.2f} | {res['fb']['cdm'][i]:>13.2f} {res['no']['cdm'][i]:>16.2f} | "
          f"{res['fb']['div'][i]:>6d} {res['no']['div'][i]:>6d} | {res['fb']['g0'][i]:>6.1f} {res['no']['g0'][i]:>6.1f}")
print(f"\nCyclinD1 fold over the SHH range:  WITH feedback {fold('fb'):.2f}x   "
      f"WITHOUT feedback {fold('no'):.2f}x   (feedback buffers the mitogen response if <)")

fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))
CF, CN = '#2c3e50', '#e74c3c'

ax[0, 0].plot(SHH, res['fb']['cdm'], CF, lw=2.2, marker='o', label='WITH EZH2-|CyclinD1 (GNP/muscle)')
ax[0, 0].plot(SHH, res['no']['cdm'], CN, lw=2.2, marker='s', ls='--', label='WITHOUT (most cells)')
ax[0, 0].set_xlabel('mitogen (SHH)'); ax[0, 0].set_ylabel('CyclinD1 transcript (mean)')
ax[0, 0].set_title('(a) The feedback holds CyclinD1 ~3x lower\n(under negative feedback vs unopposed mitogen drive)', fontweight='bold', fontsize=11)
ax[0, 0].legend(fontsize=8); ax[0, 0].grid(alpha=0.2)

ax[0, 1].plot(SHH, res['fb']['cd'], CF, lw=2.2, marker='o')
ax[0, 1].plot(SHH, res['no']['cd'], CN, lw=2.2, marker='s', ls='--')
ax[0, 1].axhline(0.4, color='gray', ls=':', alpha=0.6, label='~cycling threshold')
ax[0, 1].set_xlabel('mitogen (SHH)'); ax[0, 1].set_ylabel('CyclinD1 protein (mean)')
ax[0, 1].set_title('(b) Without feedback CyclinD1 stays above the cycling\nthreshold at all mitogen levels', fontweight='bold', fontsize=11)
ax[0, 1].legend(fontsize=8); ax[0, 1].grid(alpha=0.2)

ax[1, 0].plot(SHH, res['fb']['div'], CF, lw=2.2, marker='o', label='WITH feedback')
ax[1, 0].plot(SHH, res['no']['div'], CN, lw=2.2, marker='s', ls='--', label='WITHOUT')
ax[1, 0].set_xlabel('mitogen (SHH)'); ax[1, 0].set_ylabel('divisions / 168 h')
ax[1, 0].set_title('(c) WITHOUT feedback the cell cycles even at low mitogen;\nthe feedback installs a quiescence THRESHOLD', fontweight='bold', fontsize=11)
ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=0.2)

ax[1, 1].plot(SHH, res['fb']['g0'], CF, lw=2.2, marker='o', label='WITH feedback')
ax[1, 1].plot(SHH, res['no']['g0'], CN, lw=2.2, marker='s', ls='--', label='WITHOUT')
ax[1, 1].set_xlabel('mitogen (SHH)'); ax[1, 1].set_ylabel('transient G0 (p27) %')
ax[1, 1].set_title('(d) Feedback enables QUIESCENCE below a mitogen threshold\n(= cell-cycle EXIT capability)', fontweight='bold', fontsize=11)
ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=0.2)

fig.suptitle('EZH2 -| CyclinD1 feedback installs a mitogen-gated proliferation/quiescence switch (enables cell-cycle exit)',
             fontsize=13, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/fig_v44_feedback_role.png', dpi=170, bbox_inches='tight')
plt.savefig('simulations/fig_v44_feedback_role.pdf', bbox_inches='tight')
plt.close()
print("\nSaved: fig_v44_feedback_role.png")
