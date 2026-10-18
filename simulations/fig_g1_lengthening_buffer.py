"""Figure: graded G1/cell-cycle LENGTHENING (no transient G0) as the EZH2/CyclinD1 buffering mechanism (JP 2026-07-20).
Reads sim_g1_lengthening_buffer_cache.npz. Panels:
  A  Tc vs mitogen: the CURRENT v44 model holds Tc ~constant then abruptly ARRESTS (measured); the proposed graded
     G1 lengthening slows the cycle smoothly (GNP-grounded, Tc_min 16.25 h). B  withdrawal: CyclinD1 falls FAST,
     EZH2 (baseline + S-phase) falls SLOW via the cycle-length response -> the mark is buffered. C  a transient
     subthreshold dip: graded lengthening retains EZH2 better than an abrupt G0.
Run: ./venv/bin/python simulations/fig_g1_lengthening_buffer.py"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(HERE, 'sim_g1_lengthening_buffer_cache.npz'), allow_pickle=True)
c_cd, c_ez, c_mk, c_cur = '#c0392b', '#6c3483', '#1e8449', '#7f8c8d'

fig, ax = plt.subplots(1, 3, figsize=(16, 4.7))
fig.subplots_adjust(left=0.05, right=0.99, top=0.82, bottom=0.15, wspace=0.29)

# ---- A: Tc vs mitogen -- current model (abrupt) vs proposed graded lengthening ----
a = ax[0]; dr = Z['dr']; M, D, E, Mk, Tc, fS = dr[:, 0], dr[:, 1], dr[:, 2], dr[:, 3], dr[:, 4], dr[:, 5]
a.plot(M, Tc, '-o', color=c_ez, lw=2.4, ms=6, label='proposed: graded G1 lengthening')
# current v44 model (measured in this session): flat ~21 h above SHH~0.08, ARREST below
cur_M = np.array([0.06, 0.08, 0.1, 0.15, 0.25, 0.4, 0.6, 0.9]); cur_Tc = np.array([np.nan, np.nan, 19.7, 21.7, 21.9, 21.9, 22.0, 22.0])
a.plot(cur_M, cur_Tc, '-s', color=c_cur, lw=2.2, ms=6, label='current v44 (flat, then arrest)')
a.axvspan(0.0, 0.085, color=c_cur, alpha=0.13); a.text(0.02, 55, 'current model:\nabrupt ARREST\n(transient-G0)', fontsize=8, color='#555')
a.axhline(16.25, color='k', ls=':', lw=1, alpha=0.5); a.text(0.9, 13.5, 'GNP Tc 16.25 h', fontsize=8, ha='right')
a.set_xlabel('mitogen (Shh)'); a.set_ylabel('cell-cycle length Tc (h)'); a.set_ylim(0, 72)
a.set_title('A  Graded G1 lengthening vs the current\nabrupt arrest (no transient G0)', fontsize=11, fontweight='bold')
a.legend(fontsize=8.5, loc='upper right'); a.grid(alpha=0.2)

# ---- B: withdrawal -- CyclinD1 fast, EZH2/mark slow (buffering) ----
b = ax[1]; st = float(Z['wd_settle'][0]); t, D, E, Mk, fS, Tc = Z['wd']
b.plot(t - st, D / D[(t > st - 20) & (t < st)].mean(), color=c_cd, lw=2.6, label='CyclinD1 (fast: ½ in 1.4 h)')
b.plot(t - st, E / E[(t > st - 20) & (t < st)].mean(), color=c_ez, lw=2.6, label='EZH2 (slow: ½ in 14.8 h)')
b.axvline(0, color='k', ls='--', lw=1, alpha=0.6); b.text(3, 0.2, 'mitogen withdrawn\n(→ subthreshold)', fontsize=8)
b.annotate('', xy=(60, 0.56), xytext=(60, 0.98), arrowprops=dict(arrowstyle='<->', color='#555'))
b.text(66, 0.75, 'buffering\ngap', fontsize=8, color='#555')
b.set_xlim(-30, 300); b.set_ylim(0, 1.08); b.set_xlabel('time after withdrawal (h)'); b.set_ylabel('level (norm. to pre-withdrawal)')
b.set_title('B  Withdrawal: CyclinD1 fast, EZH2 slow\n(low-pass filter = the cycle-length response)', fontsize=11, fontweight='bold')
b.legend(fontsize=8.5, loc='center right'); b.grid(alpha=0.2)

# ---- C: transient dip -- EZH2 response by mode ----
c = ax[2]; st2 = float(Z['dip_settle'][0])
for key, col, lab in [('dip_graded', c_ez, 'graded G1 lengthening'), ('dip_abrupt', c_cd, 'abrupt G0 (current model)'), ('dip_no', c_cur, 'no lengthening (fixed Tc)')]:
    t2, D2, E2, Mk2, fS2, Tc2 = Z[key]; c.plot(t2 - st2, E2, color=col, lw=2.2, label=lab)
c.axvspan(0, 60, color='#999', alpha=0.12); c.text(4, 1.75, 'subthreshold\ndip (60 h)', fontsize=8, color='#555')
c.axhline(1.5, color='k', ls=':', lw=1, alpha=0.5); c.text(295, 1.58, 'EZH2 baseline', fontsize=7.5, ha='right')
c.set_xlim(-30, 300); c.set_xlabel('time (h)'); c.set_ylabel('EZH2')
c.set_title('C  Transient dip: graded lengthening retains\nEZH2 better than an abrupt arrest', fontsize=11, fontweight='bold')
c.legend(fontsize=8.5, loc='upper right'); c.grid(alpha=0.2)

fig.suptitle('Graded G1 lengthening (no transient G0) as the buffering mechanism: CyclinD1 tracks mitogen fast; EZH2 (baseline + S-phase) falls only as the cycle lengthens '
             '→ EZH2/H3K27me3 low-pass-filter mitogen', fontsize=11.5, fontweight='bold', y=0.975)
out = os.path.join(HERE, 'fig_g1_lengthening_buffer.png')
fig.savefig(out, dpi=150); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
