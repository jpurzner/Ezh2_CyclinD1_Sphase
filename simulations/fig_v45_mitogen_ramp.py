"""v45 Figure -- gradual mitogen INCREASE then DECREASE, with vs without the EZH2->CyclinD1 feedback.

The commitment toggle is bistable, but in v45 the ONLY state inherited across division is EZH2
(p27/CDK2 births are drawn fresh) -- so EZH2 is the cell cycle's slow MEMORY. Two experiments:

  QUASI-STATIC sweep (A, B): ramp the mitogen (CyclinD1 drive) slowly UP then DOWN through a cycling
  population, carrying the slow EZH2 between steps. WITHOUT the feedback (EZH2i=1) EZH2 cannot repress
  CyclinD1 (effective CyclinD1 = drive). WITH the feedback (EZH2i=0) EZH2 builds as cells cycle and
  REPRESSES CyclinD1, so it takes MORE mitogen to enter cycle -- the feedback RAISES the commitment
  threshold (the paper: "increases the mitogen threshold to enter cell cycle"). At a slow ramp the up
  and down sweeps coincide (EZH2 equilibrates faster than the ramp -> no static hysteresis loop).

  DYNAMIC withdrawal (C): start cycling at high mitogen, then withdraw it over ~5 cell cycles, and count
  divisions vs time. EZH2's slow fall partially de-represses CyclinD1 -> tests the "couple extra
  divisions before arrest" reading (a transient effect on the EZH2-relaxation timescale).

Run:  ./venv/bin/python simulations/fig_v45_mitogen_ramp.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v45_stochastic_commitment as v45

P21_div, sig_p, phi, sig_c, sig_ez = 0.42, 0.55, 0.55, 0.5, 0.25
K_EZ = v45._rr['K_EZH2_Cd']
kDeEZ = v45._rr['kDeEZ']


def _birth(rng):
    return P21_div * rng.lognormal(0, sig_p), phi * rng.lognormal(0, sig_c)


# ---------- (A,B) quasi-static up/down sweep, carrying the slow birth-EZH2 ----------
NSTEP, N, RELAX = 24, 26, 0.45
CD = np.concatenate([np.linspace(0.1, 3.2, NSTEP), np.linspace(3.2, 0.1, NSTEP)])
PHASE = ['up'] * NSTEP + ['down'] * NSTEP


def sweep(EZH2i, seed=1):
    rng = np.random.default_rng(seed); E = 0.6; rows = []
    for cd, ph in zip(CD, PHASE):
        divs, ez_divs, cdeff = 0, [], []
        for _ in range(N):
            p27_0, CDK2_0 = _birth(rng)
            traj, _, ez_div, status = v45.run_cell(p27_0, CDK2_0, cd, EZH2_0=E * rng.lognormal(0, sig_ez), EZH2i=EZH2i)
            if status == 'divided':
                divs += 1; ez_divs.append(ez_div); cdeff.append(float(traj['Cd'].mean()))
        E = E + RELAX * (np.mean(ez_divs) / 2 - E) if ez_divs else E * 0.6
        rows.append(dict(cd=cd, phase=ph, frac=100 * divs / N, E=E,
                         cd_eff=np.mean(cdeff) if cdeff else cd * K_EZ / (K_EZ + E * (1 - EZH2i))))
    return rows


def split(rows, key):
    up = np.array([(r['cd'], r[key]) for r in rows if r['phase'] == 'up'])
    dn = np.array([(r['cd'], r[key]) for r in rows if r['phase'] == 'down'])
    return up, dn


# ---------- (C) dynamic withdrawal: count divisions vs time as mitogen ramps down ----------
CD_HI, T_SETTLE, T_RAMP, T_MAX = 3.0, 60.0, 120.0, 320.0     # hours


def cd_of_t(t):
    if t < T_SETTLE:
        return CD_HI
    f = min(1.0, (t - T_SETTLE) / T_RAMP)
    return CD_HI * (1 - f) + 0.05


def lineage_divisions(EZH2i, M=60, seed=7):
    rng = np.random.default_rng(seed)
    div_times = []
    for _ in range(M):
        t, E = 0.0, v45.equilibrate_ezh2(CD_HI, EZH2i)
        while t < T_MAX:
            cd = cd_of_t(t)
            p27_0, CDK2_0 = _birth(rng)
            traj, _, ez_div, status = v45.run_cell(p27_0, CDK2_0, cd, EZH2_0=E * rng.lognormal(0, sig_ez), EZH2i=EZH2i)
            if status == 'divided':
                t += traj['time'][-1] / 60.0
                if t <= T_MAX:
                    div_times.append(t)
                E = ez_div / 2.0
            else:
                break                                       # arrested -> lineage stops
    return np.array(div_times), M


print('quasi-static sweeps ...')
FB, NF = sweep(0), sweep(1)
print('dynamic withdrawal lineages ...')
dtFB, M = lineage_divisions(0)
dtNF, _ = lineage_divisions(1)

CFB, CNF = '#762A83', '#1B7837'
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

# (A) threshold shift
axA = axes[0]
for rows, c, lab in [(NF, CNF, 'no feedback (EZH2i)'), (FB, CFB, 'EZH2 feedback')]:
    up, dn = split(rows, 'frac')
    axA.plot(up[:, 0], up[:, 1], '-o', color=c, lw=2, ms=4, label=f'{lab}: increasing')
    axA.plot(dn[:, 0], dn[:, 1], '--', color=c, lw=1.6, label=f'{lab}: decreasing')
for rows, c in [(NF, CNF), (FB, CFB)]:
    up, _ = split(rows, 'frac'); thr = np.interp(50, up[:, 1], up[:, 0])
    axA.axvline(thr, color=c, ls=':', lw=1, alpha=0.6)
axA.set_xlabel('mitogen (CyclinD1 drive)'); axA.set_ylabel('dividing fraction (%)'); axA.set_ylim(-3, 103)
axA.set_title('A  EZH2 raises the entry threshold', loc='left', fontweight='bold', fontsize=11)
axA.legend(fontsize=7, loc='center right')
axA.annotate('feedback shifts\nthreshold ~3x', xy=(0.40, 50), xytext=(0.9, 62),
             fontsize=7.5, color=CFB, arrowprops=dict(arrowstyle='->', color=CFB, lw=1))

# (B) effective CyclinD1 buffering
axB = axes[1]
axB.plot([0.1, 3.2], [0.1, 3.2], ':', color='#999', lw=1.2, label='no repression (y=x)')
for rows, c, lab in [(NF, CNF, 'no feedback'), (FB, CFB, 'EZH2 feedback')]:
    up, dn = split(rows, 'cd_eff')
    axB.plot(up[:, 0], up[:, 1], '-', color=c, lw=2, label=lab)
axB.set_xlabel('mitogen (CyclinD1 drive)'); axB.set_ylabel('effective CyclinD1 (cycle-mean)')
axB.set_title('B  EZH2 represses / buffers CyclinD1', loc='left', fontweight='bold', fontsize=11)
axB.legend(fontsize=8, loc='upper left')

# (C) dynamic withdrawal: cumulative divisions per cell
axC = axes[2]
tgrid = np.linspace(0, T_MAX, 200)
for dt, c, lab in [(dtNF, CNF, 'no feedback'), (dtFB, CFB, 'EZH2 feedback')]:
    cum = np.array([np.sum(dt <= tt) for tt in tgrid]) / M
    axC.plot(tgrid, cum, color=c, lw=2.2, label=lab)
axC.axvspan(T_SETTLE, T_SETTLE + T_RAMP, color='#EEE8F4', alpha=0.7)
axC.text(T_SETTLE + T_RAMP / 2, 0.3, 'mitogen\nwithdrawal', fontsize=7.5, ha='center', color='#555')
axC.set_xlabel('time (h)'); axC.set_ylabel('cumulative divisions per cell')
axC.set_title('C  Divisions during withdrawal', loc='left', fontweight='bold', fontsize=11)
axC.legend(fontsize=8, loc='upper left')
axC.annotate('feedback cells arrest SOONER\n(higher threshold dominates;\nthe de-repression buffer is weak)',
             xy=(T_SETTLE + T_RAMP, len(dtFB) / M), xytext=(120, 3.0), fontsize=6.8, color=CFB,
             arrowprops=dict(arrowstyle='->', color=CFB, lw=0.9))

fig.suptitle('v45: gradual mitogen increase $\\to$ decrease, with vs without the EZH2 $\\dashv$ CyclinD1 feedback',
             fontsize=11.5, y=1.02)
fig.tight_layout()
fig.savefig('simulations/fig_v45_mitogen_ramp.png', dpi=170, bbox_inches='tight')
fig.savefig('simulations/fig_v45_mitogen_ramp.pdf', bbox_inches='tight')
print('wrote simulations/fig_v45_mitogen_ramp.{png,pdf}')
for rows, lab in [(FB, 'EZH2 feedback'), (NF, 'no feedback')]:
    up, dn = split(rows, 'frac')
    ent = np.interp(50, up[:, 1], up[:, 0]); ext = np.interp(50, dn[::-1, 1], dn[::-1, 0])
    print(f'  {lab:16s}: entry~{ent:.2f}  exit~{ext:.2f}')
print(f'  divisions/cell by t={T_MAX:.0f}h: feedback {len(dtFB)/M:.2f}  no-feedback {len(dtNF)/M:.2f}')
