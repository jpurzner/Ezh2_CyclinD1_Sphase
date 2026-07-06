"""PROTOTYPE: cell-cycle duration / growth phase sets the H3K27me3 brake on CyclinD1.

The erasing arm of the mark is replicative dilution (Mk halved every S phase). Its impact is
set by how OFTEN it fires = 1/period. So as the cell cycle lengthens (slower growth / longer G1),
dilution matters less, the mark climbs toward its no-dilution ceiling, and CyclinD1 is more repressed.

Panels:
  (A) Mechanism: Mk(t) sawtooth at a SHORT vs LONG period (mean lines) — longer cycle rebuilds
      more mark between halvings.
  (B) The user's statement, quantified: mean Mk vs period with dilution ON vs OFF; the ON->OFF gap
      = the "dilution relief", which SHRINKS as the cycle lengthens.
  (C) Developmental slowing: growth rate ramped DOWN over developmental time -> period lengthens ->
      H3K27me3 accumulates -> CyclinD1 progressively repressed (feed-forward toward exit).
  (D) Consequence: surviving CyclinD1 fraction & division rate vs period.

Run:  ./venv/bin/python simulations/fig_v44_cellcycle_duration.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

C_MARK = '#7a5aa6'; C_TX = '#b0402f'; C_FAST = '#117a65'; C_SLOW = '#c0392b'; C_OFF = '#888'
plt.rcParams.update({'font.size': 10, 'axes.titlesize': 11, 'axes.labelsize': 10,
                     'axes.titleweight': 'bold', 'figure.dpi': 110})

M_ON = build_model_v44(with_ezh2=True, with_hh=True)
assert 'Mk = 0.5*Mk' in M_ON, 'dilution event string changed'
M_OFF = M_ON.replace('Mk = 0.5*Mk', 'Mk = 1.0*Mk')     # dilution disabled (halving -> no-op)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)


def rr(M):
    r = te.loada(M); r.integrator.setValue('relative_tolerance', 1e-6)
    r.integrator.setValue('absolute_tolerance', 1e-8)
    try: r.integrator.setValue('maximum_num_steps', 60000)
    except Exception: pass
    return r


def setgnp(r, shh=1.0, mu=0.0005):
    r.reset()
    for k, v in GNP.items(): r[k] = v
    r['SHH'] = shh; r['mu'] = float(mu)


def movavg(t, x, win=2400):
    return np.array([x[(t >= T - win / 2) & (t <= T + win / 2)].mean() for T in t])


def measure(M, mu, shh=1.0, T=15000, settle=6000):
    """return (period_h, meanMk, mean Cd_mRNA, divisions/168h) for a cycling cell at growth rate mu."""
    try:
        r = rr(M); setgnp(r, shh, mu)
        d = r.simulate(0, T, int(T / 0.5), selections=['time', 'Mk', 'Cd_mRNA', 'MPF'])
    except Exception:
        return np.nan, np.nan, np.nan, np.nan
    t = d['time']; m = t >= settle
    pk, _ = find_peaks(d['MPF'][m], prominence=0.15, distance=200)
    tp = t[m][pk]
    per = np.median(np.diff(tp)) / 60.0 if len(tp) >= 2 else np.nan
    rate = len(tp) / ((t[m][-1] - t[m][0]) / 60.0) * 168.0 if len(tp) else 0.0
    return per, float(np.nanmean(d['Mk'][m])), float(np.nanmean(d['Cd_mRNA'][m])), rate


fig = plt.figure(figsize=(14, 9.5))
gs = GridSpec(2, 2, figure=fig, hspace=0.34, wspace=0.28)

# ---------- (A) sawtooth at short vs long period ----------
axA = fig.add_subplot(gs[0, 0])
for mu, col, lab in [(0.00090, C_FAST, 'fast cycle'), (0.00033, C_SLOW, 'slow cycle')]:
    r = rr(M_ON); setgnp(r, 1.0, mu)
    d = r.simulate(0, 13000, 26000, selections=['time', 'Mk', 'MPF'])
    t = d['time']; m = t >= 3000; tt = t[m] / 60.0
    axA.plot(tt, d['Mk'][m], color=col, lw=1.3, label=lab)
    mk_mean = d['Mk'][m].mean()
    axA.axhline(mk_mean, color=col, ls='--', lw=1.2, alpha=0.8)
    pk, _ = find_peaks(d['MPF'][m], prominence=0.15, distance=200)
    per = np.median(np.diff(t[m][pk])) / 60.0 if len(pk) >= 2 else np.nan
    axA.text(tt[-1] - 2, mk_mean + 0.008, f'{lab}: {per:.0f} h  ⟨Mk⟩={mk_mean:.2f}',
             color=col, fontsize=8.5, ha='right')
axA.set_xlabel('time (h)'); axA.set_ylabel('H3K27me3  Mk')
axA.set_title('(A) Longer cycle → mark rebuilds more\nbetween S-phase halvings')
axA.legend(fontsize=8.5, loc='lower right'); axA.grid(alpha=0.15)

# ---------- (B)+(D) period sweep, dilution ON vs OFF ----------
MUS = np.geomspace(0.00022, 0.00110, 14)
perON, mkON, cdON, rateON = [], [], [], []
perOFF, mkOFF = [], []
for mu in MUS:
    p, mk, cd, rt = measure(M_ON, mu); perON.append(p); mkON.append(mk); cdON.append(cd); rateON.append(rt)
    p2, mk2, _, _ = measure(M_OFF, mu); perOFF.append(p2); mkOFF.append(mk2)
perON = np.array(perON); mkON = np.array(mkON); cdON = np.array(cdON); rateON = np.array(rateON)
perOFF = np.array(perOFF); mkOFF = np.array(mkOFF)
o = np.argsort(perON); oo = np.argsort(perOFF)

axB = fig.add_subplot(gs[0, 1])
axB.plot(perON[o], mkON[o], '-o', color=C_MARK, lw=2.2, ms=5, label='dilution ON (cycling)')
axB.plot(perOFF[oo], mkOFF[oo], '--', color=C_OFF, lw=2.0, label='dilution OFF (no-dilution ceiling)')
# shade the dilution relief (gap), interpolating OFF onto ON periods
mkOFF_i = np.interp(perON[o], perOFF[oo], mkOFF[oo])
axB.fill_between(perON[o], mkON[o], mkOFF_i, color=C_MARK, alpha=0.12)
axB.annotate('“dilution relief”\nshrinks as the\ncycle lengthens', xy=(perON[o][2], 0.5 * (mkON[o][2] + mkOFF_i[2])),
             xytext=(perON[o][2] + 4, mkON[o][2] - 0.14), fontsize=8.5, color=C_MARK,
             arrowprops=dict(arrowstyle='->', color=C_MARK))
axB.set_xlabel('cell-cycle period (h)'); axB.set_ylabel('mean H3K27me3  Mk')
axB.set_title('(B) Impact of replicative dilution\ndiminishes as duration increases')
axB.legend(fontsize=8, loc='lower right'); axB.grid(alpha=0.15)

# ---------- (C) developmental slowing time-course ----------
axC = fig.add_subplot(gs[1, 0])
mu_hi, mu_lo = 0.00095, 0.00020
t0, D = 2000.0, 16000.0
mu_of_t = lambda tt: mu_hi if tt < t0 else (mu_lo if tt > t0 + D else mu_hi + (mu_lo - mu_hi) * (tt - t0) / D)
r = rr(M_ON); setgnp(r, 1.0, mu_hi)
T, MK, CD, MP = [], [], [], []; t = 0.0; chunk = 200.0; T_end = t0 + D + 3000
while t < T_end:
    r['mu'] = float(mu_of_t(t))
    try:
        s = r.simulate(t, t + chunk, 100, selections=['time', 'Mk', 'Cd', 'MPF'])
        T.append(s['time'][1:]); MK.append(s['Mk'][1:]); CD.append(s['Cd'][1:]); MP.append(s['MPF'][1:])
    except Exception:
        break
    t += chunk
T = np.concatenate(T); MK = np.concatenate(MK); CD = np.concatenate(CD); MP = np.concatenate(MP)
th = T / 60.0
axC.plot(th, MK, color=C_MARK, lw=0.6, alpha=0.35)
axC.plot(th, movavg(T, MK), color=C_MARK, lw=2.4, label='H3K27me3 (cycle-avg)')
axC.set_xlabel('developmental time (h)  —  growth rate falling →'); axC.set_ylabel('H3K27me3  Mk', color=C_MARK)
axC.tick_params(axis='y', colors=C_MARK)
axC.set_title('(C) As the cell cycle slows, the mark accumulates\nand CyclinD1 is progressively repressed')
axCt = axC.twinx()
axCt.plot(th, movavg(T, CD), color=C_TX, lw=2.4, label='CyclinD1 (cycle-avg)')
axCt.set_ylabel('CyclinD1 (Cd)', color=C_TX); axCt.tick_params(axis='y', colors=C_TX)
# division ticks (thinning out as cycle slows)
pk, _ = find_peaks(MP, prominence=0.15, distance=150)
for p in pk:
    axC.plot(th[p], axC.get_ylim()[1], marker='v', color='#555', ms=5, clip_on=False)
# annotate lengthening period
if len(pk) >= 3:
    per_early = (T[pk[1]] - T[pk[0]]) / 60.0; per_late = (T[pk[-1]] - T[pk[-2]]) / 60.0
    axC.text(0.02, 0.06, f'cycle {per_early:.0f} h → {per_late:.0f} h  (▼ = divisions)',
             transform=axC.transAxes, fontsize=8.5, color='#555')

# ---------- (D) CyclinD1 repression & proliferation vs period ----------
axD = fig.add_subplot(gs[1, 1])
axD.plot(perON[o], cdON[o], '-s', color=C_TX, lw=2.2, ms=5, label='mean CyclinD1 transcript')
axD.set_xlabel('cell-cycle period (h)'); axD.set_ylabel('mean CyclinD1 transcript', color=C_TX)
axD.tick_params(axis='y', colors=C_TX)
axD.set_title('(D) The consequence: slower cycle →\nmore mark → less CyclinD1')
axDt = axD.twinx()
axDt.plot(perON[o], rateON[o], '-^', color='#2c3e50', lw=1.8, ms=5, alpha=0.8, label='division rate')
axDt.set_ylabel('divisions / 168 h', color='#2c3e50'); axDt.tick_params(axis='y', colors='#2c3e50')
axD.grid(alpha=0.15)

fig.suptitle('Cell-cycle duration / growth phase sets the H3K27me3 brake on CyclinD1  (replicative-dilution frequency = 1/period)',
             fontsize=13, fontweight='bold', y=0.985)
plt.savefig('simulations/fig_v44_cellcycle_duration.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/fig_v44_cellcycle_duration.pdf', bbox_inches='tight')
plt.close()
print('Saved fig_v44_cellcycle_duration.png')
print('period(h):', np.round(perON[o], 1))
print('meanMk ON:', np.round(mkON[o], 3))
print('meanMk OFF:', np.round(np.interp(perON[o], perOFF[oo], mkOFF[oo]), 3))
print('meanCd_mRNA:', np.round(cdON[o], 2))
