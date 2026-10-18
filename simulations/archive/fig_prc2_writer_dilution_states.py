"""Figure: two proliferation STATES from an EZH2-writing vs division-dilution balance on H3K27me3 (JP 2026-07-19).
Reads sim_prc2_writer_dilution_states_cache.npz. Panels:
  A  dose-response on the dividing (ON) branch: EZH2 tracks division; the mark stays HIGH because writing outpaces
     dilution (kw*kE > mu). B  mitogen up/down ramp -> HYSTERESIS = a bistable memory. C  withdrawal from fast
     division: EZH2 collapses, the mark persists, proliferation latches (memory) -- and division-coupled dilution
     (mu) makes it LESS sticky than constant turnover. D  the two coexisting states at one mitogen level.
Run: ./venv/bin/python simulations/fig_prc2_writer_dilution_states.py"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(HERE, 'sim_prc2_writer_dilution_states_cache.npz'), allow_pickle=True)
c_ez, c_mk, c_p = '#6c3483', '#1e8449', '#c0392b'

fig, ax = plt.subplots(1, 4, figsize=(19, 4.6))
fig.subplots_adjust(left=0.045, right=0.985, top=0.84, bottom=0.15, wspace=0.34)

# ---- A: dose-response on the ON branch ----
a = ax[0]; M = Z['dr_M']; dr = Z['dr']; E, Mk, p = dr[:, 0], dr[:, 1], dr[:, 2]
a.plot(M, E, '-o', color=c_ez, lw=2.2, ms=5, label='EZH2 (writer)')
a.plot(M, Mk, '-s', color=c_mk, lw=2.4, ms=5, label='H3K27me3 (mark)')
a.plot(M, p, '-^', color=c_p, lw=2.0, ms=5, label='proliferation p')
a.set_xlabel('mitogen M'); a.set_ylabel('level'); a.set_ylim(0, 1.05)
a.set_title('A  Dividing (ON) branch\nEZH2 tracks division; mark stays high\n(writing outpaces dilution)', fontsize=10.5, fontweight='bold')
a.legend(fontsize=8.5, loc='center right'); a.grid(alpha=0.2)

# ---- B: hysteresis (bistability = memory) ----
b = ax[1]
for key, col, lab in [('hys_9', c_p, 'division-coupled dilution (μ=0.9)'), ('hys_0', '#e67e22', 'constant turnover (μ=0)')]:
    up, pu, dn, pd, mku, mkd = Z[key]
    b.plot(up, pu, '-', color=col, lw=2.2, label=lab + ' ↑')
    b.plot(dn, pd, '--', color=col, lw=2.0, alpha=0.8, label=lab + ' ↓')
b.set_xlabel('mitogen M (ramp)'); b.set_ylabel('proliferation p'); b.set_ylim(-0.02, 1.02)
b.set_title('B  Mitogen up/down ramp → HYSTERESIS\n(bistable ON/OFF states = a memory)', fontsize=10.5, fontweight='bold')
b.legend(fontsize=7.5, loc='center right'); b.grid(alpha=0.2)
b.annotate('memory\n(hysteresis)', xy=(0.16, 0.5), xytext=(0.35, 0.35), fontsize=8,
           arrowprops=dict(arrowstyle='->', color='k'))

# ---- C: withdrawal from fast division ----
c = ax[2]; st = float(Z['wd_settle'][0])
t, E, Mk, p = Z['wd_9']
c.plot(t - st, E, color=c_ez, lw=2, label='EZH2 (collapses)')
c.plot(t - st, Mk, color=c_mk, lw=2.4, label='H3K27me3 (persists)')
c.plot(t - st, p, color=c_p, lw=2.2, label='proliferation (latch)')
t0, E0, Mk0, p0 = Z['wd_0']
c.plot(t0 - st, p0, ':', color=c_p, lw=1.8, alpha=0.7, label='p, μ=0 (stickier)')
c.axvline(0, color='k', ls='--', lw=1, alpha=0.6); c.text(6, 0.04, 'mitogen withdrawn', fontsize=8)
c.set_xlim(-40, 400); c.set_ylim(0, 1.05); c.set_xlabel('time after withdrawal'); c.set_ylabel('level')
c.set_title('C  Withdrawal from fast division\nEZH2 collapses, mark persists → division\nlatches (memory)', fontsize=10.5, fontweight='bold')
c.legend(fontsize=7.5, loc='center right'); c.grid(alpha=0.2)

# ---- D: the two states side by side ----
d = ax[3]
labels = ['EZH2', 'H3K27me3', 'prolif.', 'dilution']
on = [0.48, 0.60, 0.48, 0.43]; off = [0.00, 0.19, 0.00, 0.00]
x = np.arange(len(labels)); w = 0.38
d.bar(x - w / 2, on, w, color=c_p, label='ON: fast-dividing\n(high replacement)')
d.bar(x + w / 2, off, w, color='#7f8c8d', label='OFF: arrested\n(no replacement → diff. genes ON)')
d.set_xticks(x); d.set_xticklabels(labels, fontsize=9); d.set_ylabel('level'); d.set_ylim(0, 0.7)
d.set_title('D  Two coexisting states at M=0.15\n(same mitogen, different fate)', fontsize=10.5, fontweight='bold')
d.legend(fontsize=8, loc='upper right'); d.grid(alpha=0.2, axis='y')

fig.suptitle('EZH2-writing vs division-dilution sets H3K27me3 → two proliferation states: fast-dividing (high EZH2, mark maintained despite dilution) '
             'vs arrested (EZH2 collapses, mark lost, differentiation genes fire) — a bistable memory',
             fontsize=12, fontweight='bold', y=0.975)
out = os.path.join(HERE, 'fig_prc2_writer_dilution_states.png')
fig.savefig(out, dpi=150); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
