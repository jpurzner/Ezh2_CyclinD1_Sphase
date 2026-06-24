"""Ptch1 negative-feedback loop (Gli -> Ptch1 -| Smo -| Gli) and its disruption in MB.

Ptch1 is a Gli target gene: Gli activity induces Ptch1 transcription, and functional Ptch1 protein
re-inhibits Smo -- the canonical Hedgehog negative feedback. In GNPs the loop is intact, so a step of
SHH drives a transient Gli1 OVERSHOOT that then adapts (Ptch1 catches up). In MB the induced Ptch1 is
non-functional (low functional fraction = Ptch1_copy_number), so the loop is BROKEN: Gli1 rises
monotonically to a high plateau, Ptch1 mRNA is still induced (a Gli/SHH-MB marker) but cannot brake Smo.

  ./venv/bin/python simulations/fig_v44_ptch1_feedback.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te, matplotlib.pyplot as plt
from src.build_model_v44_heldt import build_model_v44

_RR = te.loada(build_model_v44())
_RR.integrator.setValue("absolute_tolerance", 1e-9)
_RR.integrator.setValue("relative_tolerance", 1e-6)
SEL = ["time", "Gli1", "Ptch1_mRNA", "Smo_active", "Gli_act"]


def step(func_frac, mycn, t_pre=1500, t_post=300):
    """Settle at SHH=0 (quiescent), then step SHH on; return the post-step trace."""
    _RR.reset()
    _RR['SHH'] = 0.0; _RR['Ptch1_copy_number'] = func_frac
    _RR['MYCN_amplification'] = mycn; _RR['HHi'] = 0
    _RR.simulate(0, t_pre, 1500)
    _RR['SHH'] = 0.5
    return _RR.simulate(0, t_post, 6000, selections=SEL)


gnp = step(1.0, 1.0)      # intact feedback
mb = step(0.1, 2.86)      # broken feedback (functional Ptch1 ~10%)


def norm(y):
    return y / y[-1] if y[-1] > 1e-6 else y

fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))

# A: Gli1 step response, normalized to final -> shows overshoot/adaptation
ax[0].axhline(1.0, color='gray', lw=0.8, ls=':')
ax[0].plot(gnp['time'], norm(gnp['Gli1']), color='#1B7837', lw=2.2,
           label=f"GNP (intact loop)  overshoot {gnp['Gli1'].max()/gnp['Gli1'][-1]:.2f}x")
ax[0].plot(mb['time'], norm(mb['Gli1']), color='#762A83', lw=2.2,
           label=f"MB (broken loop)  overshoot {mb['Gli1'].max()/mb['Gli1'][-1]:.2f}x")
ax[0].set_xlabel('time after SHH step (min)'); ax[0].set_ylabel('Gli1 / steady state')
ax[0].set_title('A  Gli1 step response\nGNP adapts (overshoot), MB does not', fontsize=10.5, fontweight='bold')
ax[0].legend(fontsize=8, loc='upper right')

# B: Ptch1 mRNA induction (absolute) -- the feedback signal, induced in both
ax[1].plot(gnp['time'], gnp['Ptch1_mRNA'], color='#1B7837', lw=2.2, label='GNP')
ax[1].plot(mb['time'], mb['Ptch1_mRNA'], color='#762A83', lw=2.2, label='MB')
ax[1].set_xlabel('time after SHH step (min)'); ax[1].set_ylabel('Ptch1 mRNA (a.u.)')
ax[1].set_title('B  Ptch1 mRNA (Gli target) is induced\nin both — but functional only in GNP',
                fontsize=10.5, fontweight='bold')
ax[1].legend(fontsize=8)

# C: steady-state Smo + Gli_act (GNP vs MB) -- the broken loop leaves Smo/Gli high in MB
labels = ['Smo active', 'Gli activator']
gnp_ss = [gnp['Smo_active'][-1], gnp['Gli_act'][-1]]
mb_ss = [mb['Smo_active'][-1], mb['Gli_act'][-1]]
x = np.arange(len(labels)); w = 0.36
ax[2].bar(x - w/2, gnp_ss, w, color='#1B7837', label='GNP (Ptch1 functional)')
ax[2].bar(x + w/2, mb_ss, w, color='#762A83', label='MB (Ptch1 non-functional)')
ax[2].set_xticks(x); ax[2].set_xticklabels(labels)
ax[2].set_ylabel('steady-state level (a.u.)')
ax[2].set_title('C  Broken loop in MB: Smo/Gli stay high\ndespite Gli-induced Ptch1',
                fontsize=10.5, fontweight='bold')
ax[2].legend(fontsize=8)

plt.tight_layout()
plt.savefig('simulations/fig_v44_ptch1_feedback.png', dpi=170, bbox_inches='tight')
plt.savefig('simulations/fig_v44_ptch1_feedback.pdf', bbox_inches='tight')
plt.close()
print("Saved: fig_v44_ptch1_feedback.png/pdf")
print(f"GNP overshoot {gnp['Gli1'].max()/gnp['Gli1'][-1]:.2f}x  |  MB overshoot {mb['Gli1'].max()/mb['Gli1'][-1]:.2f}x")
print(f"steady Smo: GNP {gnp['Smo_active'][-1]:.2f}  MB {mb['Smo_active'][-1]:.2f}")
