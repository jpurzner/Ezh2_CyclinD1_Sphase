"""Phase-plane portraits of the triangle-ramp trajectory, WITH vs WITHOUT H3K27me3 repression on CyclinD1.
Reuses the two cached runs from sim_ramp_hysteresis_trace_control (run that first).

A phase plot draws the trajectory in STATE space (state vs state), not vs time — so the feedback and its
history-dependence show up as the shape of the orbit:
  (A) CyclinD1 vs H3K27me3 mark, WITH repression  -> a NEGATIVELY-sloped orbit (mark down => CyclinD1 up)
       that encloses area = the feedback + its lag.
  (B) CyclinD1 vs H3K27me3 mark, WITHOUT repression -> CyclinD1 is ~independent of the mark (flat cloud):
       the mark still moves, but the arrow from mark to CyclinD1 is cut, so there is no relationship.
  (C) CyclinD1 vs SHH (input-output transfer) for both -> WITH opens a loop, WITHOUT retraces itself.

Run:  ./venv/bin/python simulations/sim_ramp_hysteresis_phaseplot.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CACHE = 'simulations/sim_ramp_hysteresis_trace_control_cache.npz'
if not os.path.exists(CACHE):
    sys.exit('run sim_ramp_hysteresis_trace_control.py first (need the cached runs)')
z = np.load(CACHE)
SEL = ['time', 'SHH', 'Cd', 'Mk', 'EZH2', 'mass']
D = {'W': {s: z[f'W_{s}'] for s in SEL}, 'N': {s: z[f'N_{s}'] for s in SEL}}
tt = float(z['D_LEG']) / 60.0
PREEQ = float(z['PREEQ'])


def prep(d):
    T = (d['time'] - PREEQ) / 60.0
    dt = np.median(np.diff(d['time'])); w = max(3, int(1400.0 / dt))
    cav = lambda x: np.convolve(x, np.ones(w) / w, mode='same')
    return T, cav(d['Cd']), cav(d['Mk']), d['Cd'], d['Mk'], d['SHH']


TW, cdW, mkW, rawCdW, rawMkW, shW = prep(D['W'])
TN, cdN, mkN, rawCdN, rawMkN, shN = prep(D['N'])
upW = TW <= tt; upN = TN <= tt

fig, ax = plt.subplots(1, 3, figsize=(19, 6.2))

# (A) Cd vs Mk, WITH repression — the feedback plane
a = ax[0]
a.plot(rawMkW, rawCdW, color='#bbb', lw=0.5, alpha=0.5, zorder=1)                      # raw (fast cell-cycle limit cycles)
a.plot(mkW[upW], cdW[upW], color='#1f8a4c', lw=3, label='ramp UP', zorder=3)
a.plot(mkW[~upW], cdW[~upW], color='#c0392b', lw=3, label='ramp DOWN', zorder=3)
a.annotate('arrest:\nmark high,\nCyclinD1 low', xy=(0.255, 1.2), xytext=(0.20, 3.2),
           fontsize=8.5, color='#555', ha='center', arrowprops=dict(arrowstyle='-|>', color='#555', lw=1.4))
a.annotate('fast cycling:\nmark diluted low,\nCyclinD1 high', xy=(0.14, 6.8), xytext=(0.185, 6.3),
           fontsize=8.5, color='#555', ha='center', arrowprops=dict(arrowstyle='-|>', color='#555', lw=1.4))
a.set_xlabel('H3K27me3 mark'); a.set_ylabel('CyclinD1 (cycle-avg)')
a.set_title('(A) WITH repression: CyclinD1 vs mark\nNEGATIVE feedback (mark ↑ ⟹ CyclinD1 ↓)', fontweight='bold', fontsize=10.5)
a.legend(fontsize=9, loc='center right'); a.grid(alpha=0.15)

# (B) Cd vs Mk, WITHOUT repression — arrow cut
a = ax[1]
a.plot(rawMkN, rawCdN, color='#ccc', lw=0.5, alpha=0.5, zorder=1)
a.plot(mkN[upN], cdN[upN], color='#4a6fa5', lw=3, label='ramp UP', zorder=3)
a.plot(mkN[~upN], cdN[~upN], color='#8a8a8a', lw=3, ls='--', label='ramp DOWN', zorder=3)
a.set_xlabel('H3K27me3 mark'); a.set_ylabel('CyclinD1 (cycle-avg)')
a.set_title('(B) WITHOUT repression: CyclinD1 vs mark\nno feedback — both co-drift with SHH (slope flips POSITIVE)', fontweight='bold', fontsize=10.5)
a.legend(fontsize=9, loc='lower right'); a.grid(alpha=0.15)

# (C) Cd vs SHH — input-output transfer, both conditions
a = ax[2]
a.plot(shW[upW], cdW[upW], color='#1f8a4c', lw=2.6, label='WITH — up')
a.plot(shW[~upW], cdW[~upW], color='#c0392b', lw=2.6, label='WITH — down')
a.plot(shN[upN], cdN[upN], color='#555', lw=2.0, ls='--', label='WITHOUT — up')
a.plot(shN[~upN], cdN[~upN], color='#aaa', lw=2.0, ls=':', label='WITHOUT — down')
a.set_xlabel('mitogen SHH'); a.set_ylabel('CyclinD1 (cycle-avg)')
a.set_title('(C) CyclinD1 vs SHH (input–output)\nWITH opens a loop; WITHOUT nearly retraces itself', fontweight='bold', fontsize=10.5)
a.legend(fontsize=8.5, loc='upper left'); a.grid(alpha=0.15)

fig.suptitle('Phase-plane portraits of the triangle-ramp trajectory — WITH vs WITHOUT H3K27me3 repression on CyclinD1',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('simulations/sim_ramp_hysteresis_phaseplot.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/sim_ramp_hysteresis_phaseplot.pdf', bbox_inches='tight')
plt.close()
print('Saved sim_ramp_hysteresis_phaseplot.png')
