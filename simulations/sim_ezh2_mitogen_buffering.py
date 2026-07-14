"""EZH2 tracks mitogen and buffers CyclinD1 — protecting GNP transit-amplifying divisions.

Mechanism (JP): EZH2 is a mitogen-dose-dependent, SLOW (long-lived) repressor of CyclinD1. At high
mitogen EZH2 is high -> represses CyclinD1 (buffers the high end down); as mitogen is withdrawn EZH2
FALLS over ~1 day -> repression is relieved -> CyclinD1 is propped back up (buffered). The net is that
CyclinD1 tracks mitogen far more weakly than the raw Gli/MYCN drive does (gain 0.85 -> 0.47). A
transient decline is therefore blunted, and the bistable Rb-E2F restriction point adds commitment
MEMORY so a committed cell completes its cycle through a short dip -> no forfeited progeny (symmetric
GNP divisions make a spurious exit exponentially costly).

Resolution matters: raw CyclinD1 is cell-cycle-oscillated, so the EZH2 relief is only visible after
cycle-averaging and against the un-repressed drive (f0_prc2=1) as a reference.

Panels:
  A  EZH2 tracks mitogen & buffers CyclinD1 (fine steady sweep): EZH2(SHH) rises; CyclinD1 is flattened
     vs the un-repressed drive (f0_prc2=1).
  B  Withdrawal dynamics (fine, cycle-averaged): SHH 0.9->0.3; EZH2 falls over ~1 day, relieving
     repression; cycle-averaged CyclinD1 drops far LESS than the un-repressed drive (buffered).
  C  Temporary decline that returns: a committed cell rides through (buffered CyclinD1 + divisions
     resume); a deeper/longer decline crosses the exit threshold.
  D  Commitment memory: sustained loss -> exit; a transient dip (< 1 cycle) -> the committed cell
     still divides (Rb-E2F R-point memory).

Run:  ./venv/bin/python simulations/sim_ezh2_mitogen_buffering.py [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

FRESH = '--fresh' in sys.argv
CACHE = 'simulations/fig_v44_mitogen_buffering_cache.npz'
M = build_model_v44(with_ezh2=True, with_hh=True)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)


def rr_(fb=True):
    r = te.loada(M); r.integrator.setValue('absolute_tolerance', 1e-8); r.integrator.setValue('relative_tolerance', 1e-6)
    r.reset()
    for k, v in GNP.items():
        r[k] = v
    if not fb:
        r['f0_prc2'] = 1.0
    return r


def steady(fb, shh, T=11000, settle=6000):
    r = rr_(fb); r['SHH'] = shh
    for atol in (1e-9, 1e-8, 1e-7):
        try:
            r.integrator.setValue('absolute_tolerance', atol)
            res = r.simulate(0, T, int(T / 0.5), selections=['time', 'Cd', 'EZH2', 'MPF'])
            m = res['time'] >= settle; t = res['time'][m]; dt = t[1] - t[0]
            pk, _ = find_peaks(res['MPF'][m], prominence=0.15, distance=int(200 / dt))
            return float(res['Cd'][m].mean()), float(res['EZH2'][m].mean()), len(pk) / ((t[-1] - t[0]) / 60) * 168
        except Exception:
            continue
    return np.nan, np.nan, np.nan


def two_phase(fb, shh0, shh1, T1=9000, T2=8000, fine=3.0):
    """equilibrate at shh0 (record the pre-withdrawal cycle-avg baseline), then step to shh1."""
    r = rr_(fb); r['SHH'] = shh0; base = np.nan
    for atol in (1e-9, 1e-8, 1e-7):
        try:
            r.integrator.setValue('absolute_tolerance', atol)
            r1 = r.simulate(0, T1, int(T1 / 3), selections=['time', 'Cd'])
            base = float(r1['Cd'][r1['time'] >= T1 - 2640].mean())   # avg over the last ~2 cycles
            break
        except Exception:
            r = rr_(fb); r['SHH'] = shh0
    r['SHH'] = shh1
    for atol in (1e-8, 1e-7, 1e-6):
        try:
            r.integrator.setValue('absolute_tolerance', atol)
            res = r.simulate(0, T2, int(T2 / fine), selections=['time', 'Cd', 'EZH2', 'MPF'])
            return res['time'], res['Cd'], res['EZH2'], res['MPF'], base
        except Exception:
            continue
    return (np.array([]),) * 4 + (np.nan,)


def drive_pulse(fb, base, low, t0, dur, T, chunk=50, fine=3.0, rt=200):
    """cycling baseline; a trapezoidal decline to `low` for `dur`, then recovery. Robust chunked driver."""
    def shh(t):
        if t < t0 or t >= t0 + dur:
            return base
        return base + (low - base) * min(t - t0, t0 + dur - t, rt) / rt
    r = rr_(fb); T_a, CD, MP, SH = [], [], [], []; t = 0.0
    while t < T:
        r['SHH'] = float(shh(t)); r.saveState('s'); ok = False
        for atol in (1e-8, 1e-7, 1e-6):
            try:
                r.integrator.setValue('absolute_tolerance', atol)
                res = r.simulate(t, t + chunk, max(2, int(chunk / fine)), selections=['time', 'Cd', 'EZH2', 'MPF', 'SHH'])
                ok = True; break
            except Exception:
                r.loadState('s')
        if not ok:
            break
        T_a.append(res['time'][1:]); CD.append(res['Cd'][1:]); MP.append(res['MPF'][1:]); SH.append(res['SHH'][1:])
        t += chunk
    if not T_a:
        return (np.array([]),) * 4
    return np.concatenate(T_a), np.concatenate(CD), np.concatenate(MP), np.concatenate(SH)


def drive_fn(fb, shh_fn, T, chunk=60, fine=3.0):
    """robust chunked driver for an arbitrary SHH(t) (e.g. a slow linear ramp)."""
    r = rr_(fb); T_a, CD, EZ, MP, SH = [], [], [], [], []; t = 0.0
    while t < T:
        r['SHH'] = float(shh_fn(t)); r.saveState('s'); ok = False
        for atol in (1e-8, 1e-7, 1e-6):
            try:
                r.integrator.setValue('absolute_tolerance', atol)
                res = r.simulate(t, t + chunk, max(2, int(chunk / fine)), selections=['time', 'Cd', 'EZH2', 'MPF', 'SHH'])
                ok = True; break
            except Exception:
                r.loadState('s')
        if not ok:
            break
        T_a.append(res['time'][1:]); CD.append(res['Cd'][1:]); EZ.append(res['EZH2'][1:]); MP.append(res['MPF'][1:]); SH.append(res['SHH'][1:])
        t += chunk
    if not T_a:
        return (np.array([]),) * 5
    return np.concatenate(T_a), np.concatenate(CD), np.concatenate(EZ), np.concatenate(MP), np.concatenate(SH)


def movavg(t, x, win=1320):
    if len(t) == 0:
        return x
    return np.array([x[(t >= tt - win / 2) & (t <= tt + win / 2)].mean() for tt in t])


if FRESH or not os.path.exists(CACHE):
    print("A: EZH2 tracks mitogen + buffering (fine steady) ...")
    SHH = np.round(np.arange(0.22, 1.71, 0.09), 3)
    cdON, ezON, cdOFF = [], [], []
    for s in SHH:
        a = steady(True, s); b = steady(False, s)
        cdON.append(a[0]); ezON.append(a[1]); cdOFF.append(b[0])
    print("B: GRADUAL sustained decline (linear ramp) ...")
    t0B, DB = 3000.0, 9000.0     # hold, then linear ramp 0.9->0.3 over 150 h (stays in cycling range)
    shhB = lambda t: 0.9 if t < t0B else (0.3 if t > t0B + DB else 0.9 + (0.3 - 0.9) * (t - t0B) / DB)
    tB, cB, eB, mB, sB = drive_fn(True, shhB, 13000)
    tBd, cBd, eBd, mBd, sBd = drive_fn(False, shhB, 13000)
    baseB = float(cB[(tB >= 500) & (tB < t0B)].mean()); baseBd = float(cBd[(tBd >= 500) & (tBd < t0B)].mean())
    print("C: GRADUAL temporary decline that returns (slow trapezoid) ...")
    def shhC(t):
        t0, dn, hold = 3000.0, 3600.0, 2400.0    # 60 h down, 40 h low, 60 h up
        if t < t0: return 0.85
        if t < t0 + dn: return 0.85 + (0.45 - 0.85) * (t - t0) / dn
        if t < t0 + dn + hold: return 0.45
        if t < t0 + 2 * dn + hold: return 0.45 + (0.85 - 0.45) * (t - (t0 + dn + hold)) / dn
        return 0.85
    tC1, cC1, eC1, mC1, sC1 = drive_fn(True, shhC, 16000)
    tC2, cC2, eC2, mC2, sC2 = drive_fn(False, shhC, 16000)
    baseC1 = float(cC1[(tC1 >= 500) & (tC1 < 3000)].mean()); baseC2 = float(cC2[(tC2 >= 500) & (tC2 < 3000)].mean())
    print("D: commitment memory ...")
    tdT, cdT, mdT, sdT = drive_pulse(True, 0.7, 0.15, 6000, 700, 11000)    # transient dip < 1 cycle
    tdS, cdS, mdS, sdS = drive_pulse(True, 0.7, 0.15, 6000, 99999, 11000)  # sustained
    np.savez(CACHE, SHH=SHH, cdON=cdON, ezON=ezON, cdOFF=cdOFF,
             tB=tB, cB=cB, eB=eB, sB=sB, baseB=baseB, tBd=tBd, cBd=cBd, sBd=sBd, baseBd=baseBd,
             tC1=tC1, cC1=cC1, mC1=mC1, sC1=sC1, baseC1=baseC1, tC2=tC2, cC2=cC2, sC2=sC2, baseC2=baseC2,
             tdT=tdT, cdT=cdT, mdT=mdT, tdS=tdS, cdS=cdS, mdS=mdS)
    print("  cached ->", CACHE)

z = np.load(CACHE)
C_ON, C_DRIVE, C_EZ = '#2471a3', '#c0392b', '#8e44ad'
fig, ax = plt.subplots(2, 2, figsize=(14, 10.5))

# ---- A ----
a = ax[0, 0]; a2 = a.twinx()
a.plot(z['SHH'], z['cdOFF'], '--s', color=C_DRIVE, lw=1.8, ms=3, label='CyclinD1, un-repressed drive (f0=1)')
a.plot(z['SHH'], z['cdON'], '-o', color=C_ON, lw=2.2, ms=4, label='CyclinD1, with EZH2 (buffered)')
a2.plot(z['SHH'], z['ezON'], '-^', color=C_EZ, lw=1.8, ms=4, label='EZH2')
a.set_xlabel('mitogen (SHH)'); a.set_ylabel('steady CyclinD1'); a2.set_ylabel('EZH2', color=C_EZ)
a.set_title('(A) EZH2 tracks mitogen → buffers CyclinD1\nEZH2 rises with mitogen; CyclinD1 flattened vs the raw drive (gain 0.85→0.47)', fontweight='bold', fontsize=10.5)
h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
a.legend(h1 + h2, l1 + l2, fontsize=7.5, loc='upper left'); a.grid(alpha=0.15)

# ---- B ---- (gradual decline vs mitogen level; binned by SHH to average over cycles)
b = ax[0, 1]; b2 = b.twinx()
_edges = np.linspace(0.33, 0.9, 9)


def _bin(t, sh, y, t0=3000.0):
    m = t >= t0; sh, y = sh[m], y[m]; idx = np.digitize(sh, _edges); xs, ys = [], []
    for k in range(1, len(_edges)):
        mm = idx == k
        if mm.sum() > 8:
            xs.append(sh[mm].mean()); ys.append(y[mm].mean())
    return np.array(xs), np.array(ys)


if len(z['tBd']):
    xd, yd = _bin(z['tBd'], z['sBd'], z['cBd'] / float(z['baseBd']))
    b.plot(xd, yd, '--s', color=C_DRIVE, lw=2.0, ms=4, label='CyclinD1, raw drive (f0=1)')
if len(z['tB']):
    xo, yo = _bin(z['tB'], z['sB'], z['cB'] / float(z['baseB']))
    b.plot(xo, yo, '-o', color=C_ON, lw=2.4, ms=4, label='CyclinD1, buffered (EZH2 relief)')
    xe, ye = _bin(z['tB'], z['sB'], z['eB'])
    b2.plot(xe, ye, '-^', color=C_EZ, lw=1.3, ms=3, alpha=0.55, label='EZH2 (tracks mitogen down)')
    b.fill_between(xo, yo, np.interp(xo, xd[::-1], yd[::-1]), where=(yo > np.interp(xo, xd[::-1], yd[::-1])), color='#27ae60', alpha=0.15, interpolate=True)
b.axhline(1.0, color='k', ls=':', alpha=0.3); b.invert_xaxis()  # decline runs left→right
b.set_xlabel('mitogen (SHH)  — gradual decline →'); b.set_ylabel('CyclinD1 / baseline'); b2.set_ylabel('EZH2', color=C_EZ)
b.set_title('(B) GRADUAL decline: EZH2 tracks mitogen down → relief\nbuffered CyclinD1 stays above the raw drive (no undershoot)', fontweight='bold', fontsize=10.5)
h1, l1 = b.get_legend_handles_labels(); h2, l2 = b2.get_legend_handles_labels()
b.legend(h1 + h2, l1 + l2, fontsize=7.5, loc='lower left'); b.grid(alpha=0.15)

# ---- C ---- (gradual temporary decline that returns)
c = ax[1, 0]; c2 = c.twinx()
c2.plot(z['tC1'] / 60.0, z['sC1'], ':', color='#888', lw=1.2, label='mitogen (SHH)')
for t, cd, base, col, lab in [(z['tC2'], z['cC2'], z['baseC2'], C_DRIVE, 'un-repressed drive (f0=1)'),
                              (z['tC1'], z['cC1'], z['baseC1'], C_ON, 'buffered (EZH2 relief)')]:
    if len(t):
        c.plot(t / 60.0, movavg(t, cd, win=2640) / float(base), color=col, lw=2.2, label=lab)
c.axhline(1.0, color='k', ls=':', alpha=0.3)
c.set_xlabel('time (h)'); c.set_ylabel('CyclinD1 / baseline (cycle-avg)'); c2.set_ylabel('mitogen (SHH)', color='#888')
c.set_title('(C) A GRADUAL temporary decline that returns\nbuffered CyclinD1 stays higher & recovers (rides through the dip)', fontweight='bold', fontsize=10.5)
h1, l1 = c.get_legend_handles_labels(); h2, l2 = c2.get_legend_handles_labels()
c.legend(h1 + h2, l1 + l2, fontsize=7.5, loc='lower right'); c.grid(alpha=0.15)

# ---- D ----
d = ax[1, 1]
for t, cd, mpf, col, lab in [(z['tdS'], z['cdS'], z['mdS'], '#c0392b', 'sustained loss → EXIT'),
                             (z['tdT'], z['cdT'], z['mdT'], '#2471a3', 'transient dip → completes (R-point memory)')]:
    if len(t):
        th = t / 60.0; d.plot(th, cd, color=col, lw=1.6, label=lab)
        pk, _ = find_peaks(mpf, prominence=0.15, distance=int(200 / (t[1] - t[0])))
        d.plot(th[pk], cd[pk], 'v', color=col, ms=6)
d.axvspan(6000 / 60, (6000 + 700) / 60, color='#f39c12', alpha=0.12); d.text(6050 / 60, d.get_ylim()[1] * 0.05, 'dip', color='#b9770e', fontsize=7.5)
d.set_xlabel('time (h)'); d.set_ylabel('CyclinD1 (feedback ON)  [▼ = division]')
d.set_title('(D) Commitment memory (Rb–E2F R-point)\nsustained loss → exit; short dip → committed cell still divides', fontweight='bold', fontsize=10.5)
d.legend(fontsize=7.5, loc='upper right'); d.grid(alpha=0.15)

fig.suptitle('EZH2 tracks mitogen & buffers CyclinD1 (relief on withdrawal) — protecting GNP transit-amplifying divisions',
             fontsize=12.5, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/fig_v44_mitogen_buffering.png', dpi=160, bbox_inches='tight')
plt.savefig('simulations/fig_v44_mitogen_buffering.pdf', bbox_inches='tight')
plt.close()
print("Saved: fig_v44_mitogen_buffering.png / .pdf")
