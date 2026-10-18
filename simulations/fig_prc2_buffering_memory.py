"""Figure: under what conditions does PRC2/H3K27me3 give a BUFFERING MEMORY for proliferation? (JP 2026-07-19)
Reads sim_prc2_buffering_memory_cache.npz. Two topologies for "PRC2 on the CyclinD1 axis":
  A = mark REPRESSES CyclinD1 (negative feedback, JP's literal description)
  C = mark SUSTAINS proliferation (double-negative: PRC2 silences anti-proliferative/differentiation genes)
Answer: A buffers the steady-state DOSE-RESPONSE only (no fluctuation buffering, no withdrawal memory -- a
repressor can't sustain CyclinD1 at zero mitogen). The MEMORY that persists through fluctuations + withdrawal
needs C: a supercritical, slow, proliferation-written mark LATCH. Its decline-with-differentiation is an ACTIVE
erasure, distinct from mitogen withdrawal -> division persists through withdrawal but stops with differentiation.
Run: ./venv/bin/python simulations/fig_prc2_buffering_memory.py"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(HERE, 'sim_prc2_buffering_memory_cache.npz'), allow_pickle=True)
st = float(Z['settle'][0])
c_off, c_A, c_C, c_mit = '#95a5a6', '#c0392b', '#1f6f8b', '#7f8c8d'

fig, ax = plt.subplots(2, 3, figsize=(16.5, 9.2))
fig.subplots_adjust(left=0.055, right=0.975, top=0.9, bottom=0.08, hspace=0.42, wspace=0.32)

# ---- A: dose-response buffering (topology A) + mark-timing (maximal at high mitogen) ----
a = ax[0, 0]; Ms, D0, Dg, Mk0, Mkg = Z['A_dr']
a.plot(Ms, D0 / D0.max(), '-o', color=c_off, lw=2, ms=5, label='CyclinD1, feedback OFF')
a.plot(Ms, Dg / D0.max(), '-o', color=c_A, lw=2.3, ms=5, label='CyclinD1, mark⊣CyclinD1 ON')
a.set_xscale('log'); a.set_xlabel('mitogen M (log)'); a.set_ylabel('CyclinD1 (norm.)')
a.set_title('A  Topology A buffers the DOSE-RESPONSE\n(flatter: gain 0.85→0.53)', fontsize=11, fontweight='bold')
a2 = a.twinx(); a2.plot(Ms, Mkg, ':', color='#6c3483', lw=2, label='H3K27me3 mark')
a2.set_ylabel('H3K27me3 (mark)', color='#6c3483'); a2.tick_params(axis='y', labelcolor='#6c3483'); a2.set_ylim(0, 1)
a.legend(fontsize=8, loc='upper left'); a.grid(alpha=0.2)
a.text(0.03, 0.02, 'mark rises with mitogen →\nmaximal at max division ✓', transform=a.transAxes, fontsize=8, color='#6c3483', va='bottom')

# ---- B: fluctuation buffering FAILS for topology A (division-rate CV rises with feedback) ----
b = ax[0, 1]; F = Z['F']; periods = F[:, 0]
for j, g in enumerate((0.0, 0.6, 0.85)):
    b.plot(periods, F[:, j + 1], '-o', lw=2, ms=5,
           color=[c_off, '#e67e22', c_A][j], label=f'feedback g={g}')
b.set_xlabel('mitogen fluctuation period'); b.set_ylabel('division-rate CV')
b.set_title('B  Topology A does NOT buffer fluctuations\n(CV RISES with feedback — negative-feedback lag)', fontsize=11, fontweight='bold')
b.legend(fontsize=8); b.grid(alpha=0.2)

# ---- C: withdrawal — no memory (A) vs latch memory (C) ----
c = ax[0, 2]
for key, col, lab in [('trace_off', c_off, 'feedback OFF'), ('trace_A', c_A, 'A: mark⊣CyclinD1'), ('trace_C', c_C, 'C: mark sustains (latch)')]:
    t, D, Mk, pi = Z[key]; c.plot(t - st, pi, color=col, lw=2.3, label=lab)
c.axvline(0, color='k', ls='--', lw=1, alpha=0.6); c.text(2, 0.05, 'mitogen withdrawn', fontsize=8, color='k')
c.set_xlim(-40, 220); c.set_xlabel('time after withdrawal'); c.set_ylabel('division rate π')
c.set_title('C  Withdrawal: only topology C PERSISTS\n(A = repressor can\'t sustain D at M=0)', fontsize=11, fontweight='bold')
c.legend(fontsize=8, loc='center right'); c.grid(alpha=0.2)

# ---- D: topology C is a BISTABLE latch (sharp bifurcation), not graded ----
d = ax[1, 0]; C = Z['C']
mono = C[(C[:, 0] < 0.5) & (np.isclose(C[:, 2], 0.012))]
d.plot(mono[:, 1], mono[:, 4], '-o', color=c_C, lw=2.3, ms=6)
d.axhline(2.1, color=c_off, ls='--', lw=1.5, label='feedback-off (no memory)')
d.set_xlabel('sustain strength  s'); d.set_ylabel('persistence (proliferation half-life)')
d.set_title('D  Topology C = BISTABLE latch\n(OFF below sₙ → permanent ON above; not graded)', fontsize=11, fontweight='bold')
d.legend(fontsize=8, loc='center right'); d.grid(alpha=0.2)
d.annotate('bifurcation sₙ', xy=(0.2, 300), xytext=(0.6, 200), fontsize=8,
           arrowprops=dict(arrowstyle='->', color='k'))

# ---- E: latch carries division THROUGH a transient dip (memory) vs collapse without ----
e = ax[1, 1]; dm = Z['dip_mem']; dn = Z['dip_nomem']
t, M, D, pi, Mk = dm
e.plot(t - st, pi, color=c_C, lw=2.3, label='with sustain-memory (latch)')
tn, Mn, Dn, pin, Mkn = dn
e.plot(tn - st, pin, color=c_off, lw=2.0, label='no memory (s=0)')
e.fill_between(t - st, 0, 1.05, where=((t - st >= 0) & (t - st < 40)), color=c_mit, alpha=0.13)
e.text(5, 0.1, 'mitogen\ndip', fontsize=8, color=c_mit)
e.set_xlim(-30, 180); e.set_ylim(0, 1.08); e.set_xlabel('time'); e.set_ylabel('division rate π')
e.set_title('E  Latch rides through a transient dip\n(buffering memory) — collapses without it', fontsize=11, fontweight='bold')
e.legend(fontsize=8, loc='lower right'); e.grid(alpha=0.2)

# ---- F: two-input logic — withdrawal keeps the latch (memory), differentiation ERASURE releases it ----
f = ax[1, 2]; e0 = Z['erase_off']; e1 = Z['erase_on']
t0, D0f, M0f, P0 = e0; t1, D1f, M1f, P1 = e1
f.plot(t0 - st, P0, color=c_C, lw=2.3, label='withdrawal only → latch holds (memory)')
f.plot(t1 - st, P1, color=c_A, lw=2.3, label='withdrawal + differentiation erasure → released')
f.plot(t1 - st, M1f, ':', color='#6c3483', lw=1.8, label='H3K27me3 (erasure case)')
f.axvline(0, color='k', ls='--', lw=1, alpha=0.5)
f.fill_between(t1 - st, 0, 1.05, where=((t1 - st >= 150) & (t1 - st < 180)), color='#8e44ad', alpha=0.13)
f.text(152, 0.1, 'diff.\nerasure', fontsize=8, color='#6c3483')
f.set_xlim(-30, 250); f.set_ylim(0, 1.08); f.set_xlabel('time after withdrawal'); f.set_ylabel('division rate π  /  mark')
f.set_title('F  Two-input logic: mitogen ≠ differentiation\n(persists through withdrawal, stops with differentiation)', fontsize=11, fontweight='bold')
f.legend(fontsize=7.5, loc='center right'); f.grid(alpha=0.2)

fig.suptitle('PRC2/H3K27me3 as a buffering MEMORY for proliferation: a CyclinD1 repressor only buffers the dose-response; '
             'the memory that persists through fluctuations + withdrawal needs a proliferation-written SUSTAIN latch',
             fontsize=12.5, fontweight='bold', y=0.975)
out = os.path.join(HERE, 'fig_prc2_buffering_memory.png')
fig.savefig(out, dpi=150); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
