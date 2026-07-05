"""Linear Hedgehog WITHDRAWAL: when does a GNP's cell cycle END as Hh declines, and what happens
epigenetically at exit — with vs without EZH2 repression of CyclinD1?

Mirror of fig_v44_cellcycle_entry.py. SHH is ramped DOWN 1.0 -> 0.05 over ~150 h from a cycling start.
EZH2 raises the CyclinD1 (Hh) threshold, so with the brake intact the cell should exit EARLIER (at
higher Hh); with EZH2 inhibited it keeps cycling to lower Hh. And once divisions stop, the H3K27me3
mark is no longer diluted -> it is maintained/consolidated, epigenetically locking the arrested state
(unless EZH2 is inhibited, in which case the mark clears).

Panels:
  A  CyclinD1 (cycle-avg) + falling SHH; division ticks; the cell-cycle-EXIT time (last division) marked.
  B  Cumulative divisions vs time — the plateau marks exit.
  C  H3K27me3 mark + EZH2 vs time — after exit the mark is maintained (default; no dilution → lock-in)
     vs cleared (EZH2 inhibited).

Run:  ./venv/bin/python simulations/fig_v44_cellcycle_exit.py [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

FRESH = '--fresh' in sys.argv
CACHE = 'simulations/fig_v44_cellcycle_exit_cache.npz'
M = build_model_v44(with_ezh2=True, with_hh=True)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0)
T0, DUR, T_END = 1500.0, 9000.0, 13500.0
shh_fn = lambda t: 1.0 if t < T0 else (0.05 if t > T0 + DUR else 1.0 - 0.95 * (t - T0) / DUR)


def ramp(ezh2i, chunk=60, fine=3.0):
    rr = te.loada(M); rr.integrator.setValue('relative_tolerance', 1e-6); rr.reset()
    for k, v in GNP.items():
        rr[k] = v
    rr['EZH2i'] = ezh2i; rr['SHH'] = 1.0
    out = {k: [] for k in ['time', 'Cd', 'MPF', 'SHH', 'EZH2', 'Mk']}
    t = 0.0
    while t < T_END:
        rr['SHH'] = float(shh_fn(t)); rr.saveState('s'); ok = False
        for atol in (1e-8, 1e-7, 1e-6):
            try:
                rr.integrator.setValue('absolute_tolerance', atol)
                r = rr.simulate(t, t + chunk, max(2, int(chunk / fine)), selections=['time', 'Cd', 'MPF', 'SHH', 'EZH2', 'Mk']); ok = True; break
            except Exception:
                rr.loadState('s')
        if not ok:
            break
        for k in out:
            out[k].append(r[k][1:])
        t += chunk
    return {k: np.concatenate(v) for k, v in out.items()}


if FRESH or not os.path.exists(CACHE):
    print('withdrawing Hh 1.0 -> 0.05; default vs EZH2i ...')
    A = ramp(0.0); B = ramp(1.0)
    np.savez(CACHE, **{f'A_{k}': v for k, v in A.items()}, **{f'B_{k}': v for k, v in B.items()})
    print('  cached ->', CACHE)
z = np.load(CACHE)


def divs(t, mpf):
    pk, _ = find_peaks(mpf, prominence=0.15, distance=int(200 / (t[1] - t[0])))
    return t[pk]


def mov(t, x, win=1320):
    return np.array([x[(t >= T - win / 2) & (t <= T + win / 2)].mean() for T in t])


THR = 2.0   # nominal CyclinD1 commitment threshold (a.u.)


def exit_time(t, cd):
    """last time the cycle-averaged CyclinD1 is above the commitment threshold = cell-cycle exit."""
    cav = mov(t, cd); above = np.where(cav >= THR)[0]
    return t[above[-1]] if len(above) else np.nan


exitA, exitB = exit_time(z['A_time'], z['A_Cd']), exit_time(z['B_time'], z['B_Cd'])
shhA = float(np.interp(exitA, z['A_time'], z['A_SHH'])) if exitA == exitA else np.nan
shhB = float(np.interp(exitB, z['B_time'], z['B_SHH'])) if exitB == exitB else np.nan
ndA, ndB = len(divs(z['A_time'], z['A_MPF'])), len(divs(z['B_time'], z['B_MPF']))

C_DEF, C_EZI, C_SHH, C_MK = '#117a65', '#c0392b', '#7f8c8d', '#8e44ad'
fig, ax = plt.subplots(1, 2, figsize=(15, 5.6))

# ---- A: CyclinD1 + exit ----
a = ax[0]; a2 = a.twinx()
a2.plot(z['A_time'] / 60, z['A_SHH'], ':', color=C_SHH, lw=1.6, label='SHH (Hh, withdrawing)')
a.plot(z['B_time'] / 60, mov(z['B_time'], z['B_Cd']), color=C_EZI, lw=2.4, label='CyclinD1 — EZH2 inhibited')
a.plot(z['A_time'] / 60, mov(z['A_time'], z['A_Cd']), color=C_DEF, lw=2.4, label='CyclinD1 — EZH2 active (default)')
a.axhline(THR, color='k', ls=':', alpha=0.5); a.text(2, THR, ' commitment threshold', fontsize=7.5, va='bottom', color='#555')
for ex, shh, col in [(exitB, shhB, C_EZI), (exitA, shhA, C_DEF)]:
    if ex == ex:
        a.axvline(ex / 60, color=col, ls='--', lw=1.6, alpha=0.8); a.plot(ex / 60, THR, 'v', color=col, ms=10)
        a.text(ex / 60, THR + 0.4, f'exit\nSHH={shh:.2f}', color=col, fontsize=7.5, ha='center', fontweight='bold')
a.set_xlabel('time (h)'); a.set_ylabel('CyclinD1 (cycle-avg)'); a2.set_ylabel('SHH (Hh dose)', color=C_SHH)
a.set_title('(A) Cell-cycle EXIT as Hh withdraws\nEZH2 keeps CyclinD1 lower → exits earlier, at higher Hh', fontweight='bold', fontsize=11)
h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
a.legend(h1 + h2, l1 + l2, fontsize=8, loc='upper right'); a.grid(alpha=0.15)

# ---- B: epigenetic state (reversible in GNP) ----
b = ax[1]; b2 = b.twinx()
b2.plot(z['A_time'] / 60, z['A_SHH'], ':', color=C_SHH, lw=1.4, label='SHH')
b.plot(z['A_time'] / 60, mov(z['A_time'], z['A_EZH2']), color='#5b2c6f', lw=2.2, label='EZH2 (default) — falls with cycle exit')
b.plot(z['A_time'] / 60, mov(z['A_time'], z['A_Mk']), color=C_MK, lw=2.4, label='H3K27me3 — EZH2 active')
b.plot(z['B_time'] / 60, mov(z['B_time'], z['B_Mk']), color=C_MK, lw=2.0, ls='--', alpha=0.8, label='H3K27me3 — EZH2 inhibited (cleared)')
if exitA == exitA:
    b.axvline(exitA / 60, color=C_DEF, ls='--', lw=1.3, alpha=0.6); b.text(exitA / 60, b.get_ylim()[1] * 0.05, ' exit', color=C_DEF, fontsize=8)
b.set_xlabel('time (h)'); b.set_ylabel('EZH2 / H3K27me3 mark'); b2.set_ylabel('SHH', color=C_SHH)
b.set_title('(B) The GNP arrest is REVERSIBLE, not locked\nHh loss → EZH2 falls → the mark relaxes (contrast MB: MYCN sustains EZH2 → lock)', fontweight='bold', fontsize=11)
h1, l1 = b.get_legend_handles_labels(); h2, l2 = b2.get_legend_handles_labels()
b.legend(h1 + h2, l1 + l2, fontsize=7.5, loc='upper right'); b.grid(alpha=0.15)

fig.suptitle('Linear Hedgehog withdrawal: EZH2 keeps CyclinD1 lower so the cell exits the cycle earlier (at higher Hh); the GNP arrest is reversible — EZH2 and its mark fall as the cycle stops',
             fontsize=11, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/fig_v44_cellcycle_exit.png', dpi=160, bbox_inches='tight')
plt.savefig('simulations/fig_v44_cellcycle_exit.pdf', bbox_inches='tight')
plt.close()
print(f'exit (CyclinD1<{THR}) — EZH2 active: t={exitA/60:.1f} h (SHH={shhA:.2f}), {ndA} div;  EZH2 inhibited: t={exitB/60:.1f} h (SHH={shhB:.2f}), {ndB} div')
print('Saved: fig_v44_cellcycle_exit.png / .pdf')
