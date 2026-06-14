"""v45 Figure -- single-cell anatomy: a committed cycler vs a mass-capped arrested cell.

Two columns of stacked trajectories (one cell each, from a noisy birth):
  LEFT  = committed GNP cycler (low-p27 birth): grows, dilutes Rb, commits (CDK2 toggle flips,
          p27 cleared), fires S, builds G2, divides. EZH2 integrates over the cycle; CyclinD1 dips
          as EZH2 rises.
  RIGHT = arrested cell (low CyclinD1 drive + high-p27 birth): no mitogen bootstrap -> never commits;
          mass grows to the cap and stalls; RbC floors; CDK2 stays low -> permanent G0.

Rows: mass & RbC | CDK2 activity & p27 (the toggle) | Dna & G2 clock | EZH2 & effective CyclinD1.

Run:  ./venv/bin/python simulations/fig_v45_single_cell.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
import tellurium as te
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v45_stochastic_commitment as v45

SEL = v45.SEL


def sim_raw(Cd_drive, EZH2i, p27_0, CDK2_0, EZH2_0, t_max):
    """Direct kernel run returning the full trajectory (works for arrested cells too)."""
    rr = te.loada(v45.KERNEL)
    rr['Cd_drive'] = Cd_drive; rr['EZH2i'] = EZH2i
    rr['p27'] = p27_0; rr['CDK2'] = CDK2_0; rr['EZH2'] = EZH2_0
    rr['Cdh1'] = 1.0; rr['G2p'] = 0.0; rr['Dna'] = 0.0; rr['mass'] = 1.0
    rr['Rb'] = rr['ksRb'] / rr['kdRb']
    rr.integrator.setValue("absolute_tolerance", 1e-9); rr.integrator.setValue("relative_tolerance", 1e-7)
    r = rr.simulate(0, t_max, 3000, selections=SEL)
    return {k: np.array(r[k]) for k in SEL}


# committed cycler: GNP drive, low-p27 birth -> one cycle to division
gnp_drive, _ = v45.condition_drive('GNP')
E0 = v45.equilibrate_ezh2(gnp_drive, 0)
traj_c, _, _, _ = v45.run_cell(p27_0=0.18, CDK2_0=0.6, Cd_drive=gnp_drive, EZH2_0=E0, EZH2i=0)

# arrested: GNP+HHi drive (low CyclinD1), high-p27 birth -> never commits
hhi_drive, _ = v45.condition_drive('GNP+HHi')
traj_a = sim_raw(hhi_drive, 0, p27_0=0.6, CDK2_0=0.1, EZH2_0=E0, t_max=4000)


def cdk2act(tr):
    return tr['CDK2'] / (1 + tr['p27'] / v45._rr['Ki'])


fig, axes = plt.subplots(4, 2, figsize=(12, 9.5), sharex='col')
cols = [('Committed cycler  (GNP, low-p27 birth)', traj_c, '#1B7837'),
        ('Arrested  (GNP+HHi drive, high-p27 birth)', traj_a, '#B2182B')]

for j, (title, tr, accent) in enumerate(cols):
    t = tr['time'] / 60.0
    # shade S phase (committed col only, where it exists)
    inS = (tr['Dna'] >= 0.05) & (tr['Dna'] < 0.95)
    if inS.any():
        axes[0, j].axvspan(t[inS][0], t[inS][-1], color='#FFF3BF', alpha=0.7, zorder=0)
        for i in range(4):
            axes[i, j].axvspan(t[inS][0], t[inS][-1], color='#FFF3BF', alpha=0.6, zorder=0)
    # row 0: mass & RbC
    ax = axes[0, j]; ax.plot(t, tr['mass'], color='#333', lw=2, label='mass')
    ax.axhline(v45._rr['Mmax'], color='#B2182B', ls=':', lw=1)
    ax.text(t[-1], v45._rr['Mmax'], ' Mmax', va='bottom', ha='right', fontsize=7, color='#B2182B')
    ax.set_ylabel('mass'); axt = ax.twinx()
    axt.plot(t, tr['Rb'] / tr['mass'], color='#2166AC', lw=1.6, label='RbC')
    axt.axhline(v45._rr['KfireRb'], color='#2166AC', ls='--', lw=1, alpha=0.7)
    axt.set_ylabel('RbC', color='#2166AC'); axt.tick_params(axis='y', colors='#2166AC')
    ax.set_title(title, fontweight='bold', fontsize=10.5, color=accent)
    # row 1: CDK2act & p27 (the toggle)
    ax = axes[1, j]; ax.plot(t, cdk2act(tr), color='#1B7837', lw=2, label='CDK2 activity')
    ax.axhline(v45.COMMIT_CDK2, color='#1B7837', ls=':', lw=1, alpha=0.7)
    ax.text(t[0], v45.COMMIT_CDK2, 'commit', fontsize=7, color='#1B7837', va='bottom')
    ax.set_ylabel('CDK2 activity'); ax.set_ylim(0, max(0.9, cdk2act(traj_c).max() * 1.1))
    axt = ax.twinx(); axt.plot(t, tr['p27'], color='#E08214', lw=1.6); axt.set_ylabel('p27', color='#E08214')
    axt.tick_params(axis='y', colors='#E08214')
    # row 2: Dna & G2 clock
    ax = axes[2, j]; ax.plot(t, tr['Dna'], color='#542788', lw=2, label='DNA (replication)')
    ax.plot(t, tr['G2p'], color='#8073AC', lw=1.6, ls='--', label='G2 clock')
    ax.axhline(1.0, color='#542788', ls=':', lw=0.8, alpha=0.6)
    ax.set_ylabel('DNA / G2'); ax.set_ylim(-0.05, 1.15); ax.legend(fontsize=7, loc='center left')
    # row 3: EZH2 & effective CyclinD1
    ax = axes[3, j]; ax.plot(t, tr['EZH2'], color='#E08214', lw=2, label='EZH2')
    ax.set_ylabel('EZH2', color='#E08214'); ax.tick_params(axis='y', colors='#E08214')
    ax.set_ylim(0, max(traj_c['EZH2'].max(), traj_a['EZH2'].max()) * 1.1)
    axt = ax.twinx(); axt.plot(t, tr['Cd'], color='#1B7837', lw=1.6, label='CyclinD1 (eff.)')
    axt.set_ylabel('CyclinD1 (eff.)', color='#1B7837'); axt.tick_params(axis='y', colors='#1B7837')
    axes[3, j].set_xlabel('time (h)')

axes[0, 0].annotate('', xy=(traj_c['time'][-1] / 60, 1.0), xytext=(traj_c['time'][-1] / 60, 3.5),
                    arrowprops=dict(arrowstyle='->', color='#1B7837', lw=1.5))
axes[0, 0].text(traj_c['time'][-1] / 60 - 0.3, 3.6, 'division', color='#1B7837', fontsize=8, ha='right')
fig.suptitle('v45 single-cell anatomy: commit-and-divide vs mitogen-starved arrest (mass-capped, low CDK2)',
             fontsize=12, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.98])
fig.savefig('simulations/fig_v45_single_cell.png', dpi=165, bbox_inches='tight')
fig.savefig('simulations/fig_v45_single_cell.pdf', bbox_inches='tight')
print('wrote simulations/fig_v45_single_cell.{png,pdf}')
print(f'  committed period {traj_c["time"][-1]/60:.1f}h; arrested mass plateau {traj_a["mass"][-1]:.2f} (cap {v45._rr["Mmax"]})')
