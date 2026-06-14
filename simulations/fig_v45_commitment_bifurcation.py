"""v45 Figure -- the bistable CDK2-p27 commitment toggle and the commit-or-arrest bifurcation.

Panels:
  (A) Nullclines of the CDK2-p27 toggle at birth RbC=9 for GNP (Cd=1) and MB (Cd=7): the CDK2-from-E2F
      curve and the CDK2-from-p27 curve cross at 1 or 3 points (bistable). Marks the stable low/high
      fixed points and the separatrix.
  (B) The fixed points vs Rb concentration (the G1 timer axis): as RbC dilutes, the toggle goes from
      bistable -> monostable-high (commitment). The mass cap floors RbC so low-Cd cells stay low.
  (C) Commit-or-arrest: permanent-arrest fraction vs the CyclinD1 drive, with the conditions marked --
      a sharp threshold ~Cd 0.2-0.3. GNP+HHi falls below it (arrest); MB+HHi is MYCN-floored above it.

Run:  ./venv/bin/python simulations/fig_v45_commitment_bifurcation.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v45_stochastic_commitment as v45

# kernel params for the nullcline algebra (same as fixed_points())
import tellurium as te
_rr = te.loada(v45.KERNEL)
P = {k: _rr[k] for k in ['wE', 'dmax', 'Kd0', 'nd0', 'K_CdRb', 'Km', 'nE',
                         'ksE', 'ksE0', 'kdE', 'Ki', 'ksp0', 'Kp', 'kdp', 'kdp0']}


def nullclines(RbC, Cd):
    x = np.linspace(1e-4, 3.0, 2000)                 # x = CDK2act
    ksp = P['ksp0'] / (1 + Cd / P['Kp'])
    p27 = ksp / (P['kdp'] * x + P['kdp0'])
    cdk2_from_p27 = x * (1 + p27 / P['Ki'])          # CDK2 needed to give this CDK2act at this p27
    d0 = P['dmax'] * Cd ** P['nd0'] / (P['Kd0'] ** P['nd0'] + Cd ** P['nd0'])
    drive = d0 + P['wE'] * x
    RbP = drive / (P['K_CdRb'] * RbC + drive)
    E2F = RbP ** P['nE'] / (P['Km'] ** P['nE'] + RbP ** P['nE'])
    cdk2_from_E2F = (P['ksE'] * E2F + P['ksE0']) / P['kdE']
    return x, cdk2_from_p27, cdk2_from_E2F


fig, axes = plt.subplots(1, 3, figsize=(14, 4.7))


def cdk2_total(xf, Cd):
    return xf * (1 + (P['ksp0'] / (1 + Cd / P['Kp']) / (P['kdp'] * xf + P['kdp0'])) / P['Ki'])


# ---- (A) nullclines at birth RbC=9 (GNP, Cd=1 -- the clearest bistable example) ----
axA = axes[0]
Cd = 1.0
x, fp27, fE2F = nullclines(9.0, Cd)
axA.plot(x, fE2F, color='#1B7837', lw=2.2, label='E2F $\\to$ CDK2 (activation)')
axA.plot(x, fp27, color='#E08214', lw=2.2, ls='--', label='p27 $\\to$ CDK2 (brake)')
fps = v45.fixed_points(9.0, Cd)
names = ['OFF\n(quiescent)', 'separatrix', 'ON\n(committed)']
for j, xf in enumerate(fps):
    stable = (j != 1)
    axA.plot(xf, cdk2_total(xf, Cd), 'o' if stable else 'X', color='#333',
             ms=12 if stable else 10, mfc='#333' if stable else 'white', mew=2, zorder=5)
    axA.annotate(names[j] if len(fps) == 3 else 'ON', (xf, cdk2_total(xf, Cd)),
                 textcoords='offset points', xytext=(6, 8 if j != 0 else -2), fontsize=7.5,
                 va='top' if j == 0 else 'bottom', fontweight='bold')
axA.set_xlabel('CDK2 activity'); axA.set_ylabel('CDK2 (total)')
axA.set_xlim(0, 1.2); axA.set_ylim(0, 2.6)
axA.set_title('A  Bistable CDK2-p27 toggle\n    (GNP at birth, RbC=9)', loc='left', fontweight='bold', fontsize=10.5)
axA.legend(fontsize=8, loc='upper center')
axA.text(0.97, 0.04, 'crossings = fixed points\n$\\bullet$ stable   $\\times$ separatrix',
         transform=axA.transAxes, ha='right', fontsize=7.3, color='#444')

# ---- (B) fixed points vs RbC (the G1 timer axis) ----
axB = axes[1]
RbCs = np.linspace(1.0, 11, 60)
for Cd, c, lab in [(1.0, '#2166AC', 'GNP'), (7.0, '#B2182B', 'MB'), (0.3, '#999999', 'low Cd (HHi)')]:
    lo, hi, sep = [], [], []
    for R in RbCs:
        fps = v45.fixed_points(R, Cd)
        if len(fps) >= 3:
            lo.append((R, fps[0])); sep.append((R, fps[1])); hi.append((R, fps[-1]))
        elif len(fps) == 1:
            (hi if fps[0] > 1.0 else lo).append((R, fps[0]))
    for arr, ls in [(hi, '-'), (lo, '-'), (sep, ':')]:
        if arr:
            a = np.array(arr); axB.plot(a[:, 0], a[:, 1], ls, color=c, lw=1.8,
                                        label=lab if ls == '-' and arr is hi else None)
axB.axvline(9, color='k', ls=':', lw=1, alpha=0.5); axB.text(9, 2.55, 'birth', fontsize=7, ha='center')
axB.axvspan(1.0, 9 / v45._rr['Mmax'], color='#FDE0CE', alpha=0.5)
axB.text(9 / v45._rr['Mmax'] - 0.1, 0.15, 'mass-cap\nfloor', fontsize=6.8, ha='right', color='#B2182B')
axB.set_xlabel('Rb concentration  RbC = Rb/mass  (dilutes as cell grows $\\to$)')
axB.set_ylabel('CDK2 activity (fixed points)'); axB.invert_xaxis(); axB.set_ylim(0, 2.7)
axB.set_title('B  Rb-dilution G1 timer drives commitment', loc='left', fontweight='bold', fontsize=10.5)
axB.legend(fontsize=8, loc='center left')

# ---- (C) commit-or-arrest threshold vs Cd ----
axC = axes[2]
Cds = [0.0, 0.05, 0.1, 0.16, 0.2, 0.25, 0.3, 0.4, 0.5, 0.78, 1.0, 2.0]
arrests = []
print('scanning arrest vs Cd ...')
for cd in Cds:
    _, st = v45.ensemble(cd, N=150, return_stats=True)
    arrests.append(100 * st['arrest_frac'])
axC.plot(Cds, arrests, '-o', color='#762A83', lw=2, ms=5)
axC.set_xscale('symlog', linthresh=0.1)
# mark conditions
marks = {'GNP+HHi': 0.16, 'MB+HHi': 0.78, 'GNP': 1.0}
for name, cd in marks.items():
    a = np.interp(cd, Cds, arrests)
    axC.annotate(name, xy=(cd, a), xytext=(cd, a + 12), fontsize=7.5, ha='center',
                 arrowprops=dict(arrowstyle='->', lw=0.9), color='#B2182B' if 'MB' in name else '#2166AC')
axC.set_xlabel('CyclinD1 drive (mitogen)'); axC.set_ylabel('permanent arrest (%)')
axC.set_ylim(-3, 100)
axC.set_title('C  Commit-or-arrest bifurcation', loc='left', fontweight='bold', fontsize=10.5)

fig.suptitle('v45: the mechanistic commitment switch -- bistable CDK2-p27 toggle, Rb-dilution G1 timer, mitogen-gated arrest',
             fontsize=11, y=1.02)
fig.tight_layout()
fig.savefig('simulations/fig_v45_commitment_bifurcation.png', dpi=170, bbox_inches='tight')
fig.savefig('simulations/fig_v45_commitment_bifurcation.pdf', bbox_inches='tight')
print('wrote simulations/fig_v45_commitment_bifurcation.{png,pdf}')
