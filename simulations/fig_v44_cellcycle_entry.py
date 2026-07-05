"""When does a GNP enter the cell cycle as Hedgehog gradually rises — WITH vs WITHOUT EZH2 repression
of CyclinD1? EZH2 raises the CyclinD1 (Hh) threshold to commit, so with the brake intact a cell waits
until Hh climbs to a genuine threshold before it divides; with EZH2 inhibited (de-repressed), CyclinD1
clears the commitment threshold on basal drive and the cell enters early / promiscuously.

SHH is ramped 0 -> 1.0 over ~150 h (from a Hh-free, arrested start). We track CyclinD1 and division
events (MPF peaks) for the default model (EZH2 active) vs EZH2 inhibited (EZH2i=1), and mark the first
division = cell-cycle entry, with the Hh level at entry.

Panels:
  A  CyclinD1 (cycle-avg) vs time for both, with the rising SHH; division ticks; the cell-cycle-entry
     time (first division) marked for each.
  B  Cumulative divisions vs time — the curve lifts off at cell-cycle entry; the horizontal gap between
     the two lift-offs is the entry delay EZH2 imposes.

Run:  ./venv/bin/python simulations/fig_v44_cellcycle_entry.py [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

FRESH = '--fresh' in sys.argv
CACHE = 'simulations/fig_v44_cellcycle_entry_cache.npz'
M = build_model_v44(with_ezh2=True, with_hh=True)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0)
T0, DUR, T_END = 1000.0, 9000.0, 12500.0
shh_fn = lambda t: 0.0 if t < T0 else (1.0 if t > T0 + DUR else 1.0 * (t - T0) / DUR)


def ramp(ezh2i, chunk=60, fine=3.0):
    rr = te.loada(M); rr.integrator.setValue('relative_tolerance', 1e-6); rr.reset()
    for k, v in GNP.items():
        rr[k] = v
    rr['EZH2i'] = ezh2i; rr['SHH'] = 0.0
    T, CD, MP, SH = [], [], [], []
    t = 0.0
    while t < T_END:
        rr['SHH'] = float(shh_fn(t)); rr.saveState('s'); ok = False
        for atol in (1e-8, 1e-7, 1e-6):
            try:
                rr.integrator.setValue('absolute_tolerance', atol)
                r = rr.simulate(t, t + chunk, max(2, int(chunk / fine)), selections=['time', 'Cd', 'MPF', 'SHH']); ok = True; break
            except Exception:
                rr.loadState('s')
        if not ok:
            break
        T.append(r['time'][1:]); CD.append(r['Cd'][1:]); MP.append(r['MPF'][1:]); SH.append(r['SHH'][1:])
        t += chunk
    return np.concatenate(T), np.concatenate(CD), np.concatenate(MP), np.concatenate(SH)


if FRESH or not os.path.exists(CACHE):
    print('ramping Hh 0 -> 1.0; default vs EZH2i ...')
    tA, cA, mA, sA = ramp(0.0)   # default: EZH2 active
    tB, cB, mB, sB = ramp(1.0)   # EZH2 inhibited
    np.savez(CACHE, tA=tA, cA=cA, mA=mA, sA=sA, tB=tB, cB=cB, mB=mB, sB=sB)
    print('  cached ->', CACHE)
z = np.load(CACHE)


def divs(t, mpf):
    pk, _ = find_peaks(mpf, prominence=0.15, distance=int(200 / (t[1] - t[0])))
    return t[pk]


def mov(t, x, win=1320):
    return np.array([x[(t >= T - win / 2) & (t <= T + win / 2)].mean() for T in t])


pkA, pkB = divs(z['tA'], z['mA']), divs(z['tB'], z['mB'])
entryA = pkA[pkA >= T0][0] if np.any(pkA >= T0) else np.nan   # first division after ramp start
entryB = pkB[pkB >= T0][0] if np.any(pkB >= T0) else np.nan
shhA = float(np.interp(entryA, z['tA'], z['sA'])) if entryA == entryA else np.nan
shhB = float(np.interp(entryB, z['tB'], z['sB'])) if entryB == entryB else np.nan

C_DEF, C_EZI, C_SHH = '#117a65', '#c0392b', '#7f8c8d'
fig, ax = plt.subplots(1, 2, figsize=(15, 5.6))

# ---- A: CyclinD1 + entry ----
a = ax[0]; a2 = a.twinx()
a2.plot(z['tA'] / 60, z['sA'], ':', color=C_SHH, lw=1.6, label='SHH (Hh, ramping)')
a.plot(z['tB'] / 60, mov(z['tB'], z['cB']), color=C_EZI, lw=2.4, label='CyclinD1 — EZH2 inhibited')
a.plot(z['tA'] / 60, mov(z['tA'], z['cA']), color=C_DEF, lw=2.4, label='CyclinD1 — EZH2 active (default)')
for entry, shh, col in [(entryB, shhB, C_EZI), (entryA, shhA, C_DEF)]:
    if entry == entry:
        a.axvline(entry / 60, color=col, ls='--', lw=1.6, alpha=0.8)
        a.plot(entry / 60, 0, 'v', color=col, ms=10, clip_on=False)
a.set_xlabel('time (h)'); a.set_ylabel('CyclinD1 (cycle-avg)'); a2.set_ylabel('SHH (Hh dose)', color=C_SHH)
a.set_title('(A) Cell-cycle entry as Hh rises\nEZH2 delays commitment until Hh clears a threshold', fontweight='bold', fontsize=11)
h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
a.legend(h1 + h2, l1 + l2, fontsize=8, loc='upper left'); a.grid(alpha=0.15)

# ---- B: cumulative divisions ----
b = ax[1]; b2 = b.twinx()
b2.plot(z['tA'] / 60, z['sA'], ':', color=C_SHH, lw=1.4)
for t, pk, col, lab, entry, shh in [(z['tB'], pkB, C_EZI, 'EZH2 inhibited', entryB, shhB),
                                     (z['tA'], pkA, C_DEF, 'EZH2 active (default)', entryA, shhA)]:
    pk = pk[pk >= 0]
    steps_t = np.concatenate([[0], np.repeat(pk, 2), [T_END]]) / 60
    steps_n = np.concatenate([[0], np.repeat(np.arange(1, len(pk) + 1), 2)[:-1], [len(pk), len(pk)]])[:len(steps_t)]
    b.step(np.sort(pk) / 60, np.arange(1, len(pk) + 1), where='post', color=col, lw=2.4, label=lab)
    if entry == entry:
        b.annotate(f'entry: t={entry/60:.0f} h\nSHH={shh:.2f}', xy=(entry / 60, 0.5),
                   xytext=(entry / 60 + 8, 1.5 if col == C_DEF else 4), fontsize=8, color=col, fontweight='bold',
                   arrowprops=dict(arrowstyle='->', color=col))
b.set_xlabel('time (h)'); b.set_ylabel('cumulative divisions'); b2.set_ylabel('SHH (Hh dose)', color=C_SHH)
b.set_title('(B) When they enter the cycle\nEZH2i lets cells commit on basal Hh; EZH2 gates entry to higher Hh', fontweight='bold', fontsize=11)
b.legend(fontsize=8, loc='upper left'); b.grid(alpha=0.15)

fig.suptitle('Cell-cycle entry under gradually rising Hedgehog: EZH2 repression of CyclinD1 delays commitment until a Hh threshold; EZH2 inhibition lets cells enter early',
             fontsize=11.5, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/fig_v44_cellcycle_entry.png', dpi=160, bbox_inches='tight')
plt.savefig('simulations/fig_v44_cellcycle_entry.pdf', bbox_inches='tight')
plt.close()
print(f'entry — EZH2 active: t={entryA/60:.1f} h (SHH={shhA:.2f});  EZH2 inhibited: t={entryB/60:.1f} h (SHH={shhB:.2f})')
print('Saved: fig_v44_cellcycle_entry.png / .pdf')
