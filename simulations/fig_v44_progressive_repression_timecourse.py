"""Time-course version of the progressive-repression figure: a GNP experiences GRADUALLY INCREASING
Hedgehog. As Hh climbs, Gli1 and the CyclinD1 drive rise quickly, but the (slow, mitogen-responsive)
EZH2 -> H3K27me3 brake rises behind them and progressively reins CyclinD1 in — the repression develops
in real time.

SHH is ramped linearly from 0.2 -> 2.0 over ~180 h (gradual, slower than the ~day EZH2 timescale).
Everything is cycle-averaged (moving average over ~1 cell cycle) to remove the cell-cycle oscillation.
The un-repressed CyclinD1 drive is recovered self-consistently as Cd_mRNA / R(Mk).

Panels:
  A  Inputs/mediators vs time: SHH (input) and Gli1 (gas) rise fast; EZH2 + H3K27me3 (brake) lag.
  B  CyclinD1 vs time: the raw drive (un-repressed) vs the actual transcript; shaded gap = H3K27me3
     repression, growing over time. Twin axis: surviving fraction R = transcript/drive, falling as the
     brake engages.

Run:  ./venv/bin/python simulations/fig_v44_progressive_repression_timecourse.py [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt, tellurium as te
from src.build_model_v44_heldt import build_model_v44

FRESH = '--fresh' in sys.argv
CACHE = 'simulations/fig_v44_progressive_repression_timecourse_cache.npz'
M = build_model_v44(with_ezh2=True, with_hh=True)
_rr = te.loada(M)
F0, KMK, NMK = _rr['f0_prc2'], _rr['K_prc2'], _rr['n_prc2']
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)
SEL = ['time', 'Cd_mRNA', 'Mk', 'EZH2', 'Gli_act', 'Gli1', 'SHH']

T0, DUR, T_END = 2000.0, 10800.0, 16000.0     # hold 0.2, ramp 0.2->2.0 over 180 h, then hold
shh_fn = lambda t: 0.2 if t < T0 else (2.0 if t > T0 + DUR else 0.2 + 1.8 * (t - T0) / DUR)


def drive_ramp(chunk=60, fine=3.0):
    rr = te.loada(M); rr.integrator.setValue('relative_tolerance', 1e-6); rr.reset()
    for k, v in GNP.items():
        rr[k] = v
    out = {k: [] for k in SEL}
    t = 0.0
    while t < T_END:
        rr['SHH'] = float(shh_fn(t)); rr.saveState('s'); ok = False
        for atol in (1e-8, 1e-7, 1e-6):
            try:
                rr.integrator.setValue('absolute_tolerance', atol)
                r = rr.simulate(t, t + chunk, max(2, int(chunk / fine)), selections=SEL); ok = True; break
            except Exception:
                rr.loadState('s')
        if not ok:
            break
        for k in SEL:
            out[k].append(r[k][1:])
        t += chunk
    return {k: np.concatenate(v) for k, v in out.items()}


if FRESH or not os.path.exists(CACHE):
    print('ramping Hh 0.2 -> 2.0 over 180 h ...')
    d = drive_ramp()
    np.savez(CACHE, **d)
    print('  cached ->', CACHE)
z = np.load(CACHE)
t = z['time']; th = t / 60.0


def mov(x, win=2640):
    return np.array([x[(t >= T - win / 2) & (t <= T + win / 2)].mean() for T in t])


shh = mov(z['SHH']); gli1 = mov(z['Gli_act'] + z['Gli1']); ez = mov(z['EZH2']); mk = mov(z['Mk'])
cdm = mov(z['Cd_mRNA'])
R = F0 + (1 - F0) / (1 + (mk / KMK) ** NMK)
drive = cdm / R

C_GLI, C_EZ, C_CD, C_DRIVE, C_SHH = '#e67e22', '#8e44ad', '#117a65', '#c0392b', '#7f8c8d'
fig, ax = plt.subplots(1, 2, figsize=(15, 5.6))

# ---- A: inputs / mediators vs time ----
a = ax[0]; a2 = a.twinx()
a.plot(th, shh, color=C_SHH, lw=1.6, ls=':', label='SHH (Hh input, ramping)')
a.plot(th, gli1 / gli1.max(), color=C_GLI, lw=2.4, label='Gli1 (gas) — rises fast')
a2.plot(th, ez, color=C_EZ, lw=2.4, label='EZH2 (brake) — lags')
a2.plot(th, mk, color='#5b2c6f', lw=2.0, alpha=0.8, label='H3K27me3 mark')
a.set_xlabel('time (h)'); a.set_ylabel('SHH  /  Gli1 (norm.)'); a2.set_ylabel('EZH2 / H3K27me3', color=C_EZ)
a.set_title('(A) Gradually rising Hh: the gas leads, the brake lags\nGli1 tracks Hh quickly; slow EZH2→H3K27me3 catches up behind', fontweight='bold', fontsize=11)
h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
a.legend(h1 + h2, l1 + l2, fontsize=8, loc='upper left'); a.grid(alpha=0.15)

# ---- B: CyclinD1 drive vs actual + surviving fraction ----
b = ax[1]; b2 = b.twinx()
b.plot(th, drive, '--', color=C_DRIVE, lw=2.2, label='CyclinD1 drive (un-repressed)')
b.plot(th, cdm, color=C_CD, lw=2.6, label='CyclinD1 transcript (actual)')
b.fill_between(th, cdm, drive, color=C_EZ, alpha=0.16, label='repressed by H3K27me3')
b2.plot(th, R, color=C_EZ, lw=1.8, alpha=0.9, label='surviving fraction R')
b2.set_ylim(0, 1.05)
b.set_xlabel('time (h)'); b.set_ylabel('CyclinD1 transcript (a.u.)'); b2.set_ylabel('surviving fraction  R = transcript/drive', color=C_EZ)
b.set_title('(B) CyclinD1 reined in as the brake engages\nthe drive climbs with Hh; H3K27me3 progressively represses it', fontweight='bold', fontsize=11)
h1, l1 = b.get_legend_handles_labels(); h2, l2 = b2.get_legend_handles_labels()
b.legend(h1 + h2, l1 + l2, fontsize=8, loc='upper left'); b.grid(alpha=0.15)

fig.suptitle('Gradually increasing Hedgehog (time course): the mitogen-responsive EZH2→H3K27me3 brake lags the Gli1 gas and progressively represses CyclinD1',
             fontsize=12, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/fig_v44_progressive_repression_timecourse.png', dpi=160, bbox_inches='tight')
plt.savefig('simulations/fig_v44_progressive_repression_timecourse.pdf', bbox_inches='tight')
plt.close()
i0 = np.argmin(np.abs(th - 40)); i1 = np.argmin(np.abs(th - 250))
print(f'R over the ramp: {R[i0]:.2f} -> {R[i1]:.2f};  transcript {cdm[i0]:.2f} -> {cdm[i1]:.2f};  drive {drive[i0]:.2f} -> {drive[i1]:.2f}')
print('Saved: fig_v44_progressive_repression_timecourse.png / .pdf')
