"""Triangle mitogen ramp trace, WITH vs WITHOUT H3K27me3 repression on CyclinD1.

Same continuous triangle ramp as sim_ramp_hysteresis_trace, run twice:
  WITH   : f0_mk = 0.233  (mark represses CyclinD1)
  WITHOUT: f0_mk = 1.00   (mark still accumulates/dilutes & EZH2 runs, but R=1 for all Mk -> zero effect on Cd)

Shows directly that the mark is what (a) lowers/gates CyclinD1 and (b) makes entry/exit HISTORY-dependent:
without it the up- and down-legs retrace the same path (entry SHH ~= exit SHH); with it the cell enters
late and exits late (a loop). Caches the two runs for the phase-plot figure.

Run:  ./venv/bin/python simulations/sim_ramp_hysteresis_trace_control.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

CACHE = 'simulations/sim_ramp_hysteresis_trace_control_cache.npz'
LN2 = np.log(2.0)
M = build_model_v44(with_ezh2=True, with_hh=True)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)
DEFAULTS = dict(f0_mk=0.233, mu=0.0005, kTlEZ=0.004, del_mk=0.0007, kDeEZ=0.00015, HU=0.0, k_Cd_translation=0.801)
D_LEG, SHH_LO, SHH_HI, CHUNK = 9000.0, 0.05, 1.10, 60.0
PREEQ = float(min(20000.0, max(3000.0, 4.0 * LN2 / DEFAULTS['kDeEZ'])))
SEL = ['time', 'SHH', 'Cd', 'Mk', 'EZH2', 'mass']

_R = te.loada(M); _R.integrator.setValue('relative_tolerance', 1e-6)
try: _R.integrator.setValue('maximum_num_steps', 100000)
except Exception: pass


def shh_of(t):
    x = t - PREEQ
    if x <= 0: return SHH_LO
    if x <= D_LEG: return SHH_LO + (SHH_HI - SHH_LO) * (x / D_LEG)
    return SHH_HI + (SHH_LO - 0.02 - SHH_HI) * min((x - D_LEG) / D_LEG, 1.0)


def run(f0):
    _R.reset()
    for k, v in GNP.items(): _R[k] = v
    for k, v in DEFAULTS.items(): _R[k] = v
    _R['f0_mk'] = f0
    _R['SHH'] = SHH_LO
    for atol in (1e-8, 1e-7, 1e-6):
        try:
            _R.integrator.setValue('absolute_tolerance', atol); _R.simulate(0, PREEQ, 50, selections=['time']); break
        except Exception: pass
    rows = {s: [] for s in SEL}; t = PREEQ; TEND = PREEQ + 2 * D_LEG + 200
    while t < TEND:
        _R['SHH'] = float(shh_of(t))
        r = None
        for atol in (1e-8, 1e-7, 1e-6):
            try:
                _R.integrator.setValue('absolute_tolerance', atol); r = _R.simulate(t, t + CHUNK, 12, selections=SEL); break
            except Exception: r = None
        if r is None: break
        for s in SEL: rows[s].append(r[s][1:])
        t += CHUNK
    for s in SEL: rows[s] = np.concatenate(rows[s])
    return rows


if not os.path.exists(CACHE):
    W = run(0.233); N = run(1.0)
    np.savez(CACHE, PREEQ=PREEQ, D_LEG=D_LEG,
             **{f'W_{s}': W[s] for s in SEL}, **{f'N_{s}': N[s] for s in SEL})
    print('cached ->', CACHE, flush=True)

z = np.load(CACHE)
D = {'W': {s: z[f'W_{s}'] for s in SEL}, 'N': {s: z[f'N_{s}'] for s in SEL}}
tt = float(z['D_LEG']) / 60.0


def prep(d):
    T = (d['time'] - float(z['PREEQ'])) / 60.0
    dt = np.median(np.diff(d['time'])); w = max(3, int(1400.0 / dt))
    cav = np.convolve(d['Cd'], np.ones(w) / w, mode='same')
    divs = np.where(np.diff(d['mass']) < -0.15 * d['mass'][:-1])[0]
    return T, cav, divs


TW, cavW, dW = prep(D['W']); TN, cavN, dN = prep(D['N'])
up_W = TW <= tt; up_N = TN <= tt
enW = D['W']['SHH'][dW[0]] if len(dW) else np.nan
exW = D['W']['SHH'][dW[(TW[dW] > tt)][-1]] if (TW[dW] > tt).any() else np.nan
enN = D['N']['SHH'][dN[0]] if len(dN) else np.nan
exN = D['N']['SHH'][dN[(TN[dN] > tt)][-1]] if (TN[dN] > tt).any() else np.nan

fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(3, 2, width_ratios=[1.5, 1])
axS = fig.add_subplot(gs[0, 0]); axC = fig.add_subplot(gs[1, 0], sharex=axS); axM = fig.add_subplot(gs[2, 0], sharex=axS)
axL = fig.add_subplot(gs[:, 1])
for a in (axS, axC, axM):
    a.axvspan(0, tt, color='#1f8a4c', alpha=0.05); a.axvspan(tt, TW.max(), color='#c0392b', alpha=0.05); a.axvline(tt, color='k', lw=0.8, ls=':')

# (A) SHH
axS.plot(TW, D['W']['SHH'], color='#333', lw=2)
axS.set_ylabel('SHH'); axS.set_title('(A) Triangle mitogen ramp (identical for both)', fontweight='bold', fontsize=11)
axS.text(tt * 0.5, 1.0, 'UP', color='#1f8a4c', ha='center', fontweight='bold'); axS.text(tt * 1.5, 1.0, 'DOWN', color='#c0392b', ha='center', fontweight='bold')

# (B) CyclinD1 cycle-avg WITH vs WITHOUT
axC.plot(TN, cavN, color='#888', lw=2.4, label='WITHOUT repression (f0_mk=1)')
axC.plot(TW, cavW, color='#0b3d61', lw=2.4, label='WITH mark (f0_mk=0.233)')
axC.plot(TW, D['W']['Cd'], color='#1f6fb0', lw=0.5, alpha=0.3); axC.plot(TN, D['N']['Cd'], color='#aaa', lw=0.5, alpha=0.3)
axC.set_ylabel('CyclinD1 (cycle-avg)'); axC.legend(fontsize=9, loc='upper left')
axC.set_title('(B) CyclinD1: without the mark it is HIGHER and symmetric up vs down; with the mark it is gated & lags', fontweight='bold', fontsize=10)

# (C) mark WITH vs WITHOUT
axM.plot(TN, D['N']['Mk'], color='#c9a3e0', lw=1.6, label='mark (WITHOUT repression) — still runs, no effect on Cd')
axM.plot(TW, D['W']['Mk'], color='#8e44ad', lw=1.8, label='mark (WITH repression)')
axM.set_ylabel('H3K27me3 mark'); axM.set_xlabel('time since ramp start (h)'); axM.legend(fontsize=9, loc='upper left')
axM.set_title('(C) The mark accumulates in arrest & dilutes with divisions in BOTH; only WITH does it feed back on CyclinD1', fontweight='bold', fontsize=10)

# (D) loops
axL.plot(D['W']['SHH'][up_W], cavW[up_W], color='#1f8a4c', lw=2.6, label='WITH — up-leg')
axL.plot(D['W']['SHH'][~up_W], cavW[~up_W], color='#c0392b', lw=2.6, label='WITH — down-leg')
axL.plot(D['N']['SHH'][up_N], cavN[up_N], color='#555', lw=2.0, ls='--', label='WITHOUT — up-leg')
axL.plot(D['N']['SHH'][~up_N], cavN[~up_N], color='#999', lw=2.0, ls=':', label='WITHOUT — down-leg')
for x, c in [(enW, '#1f8a4c'), (exW, '#c0392b')]:
    if np.isfinite(x): axL.axvline(x, color=c, ls='--', lw=1, alpha=0.6)
axL.set_xlabel('mitogen SHH'); axL.set_ylabel('CyclinD1 (cycle-avg)')
axL.set_title(f'(D) Input–output loop\nWITH: entry {enW:.2f} / exit {exW:.2f} (loop)\nWITHOUT: entry {enN:.2f} / exit {exN:.2f} (retraces — nearly no loop)', fontweight='bold', fontsize=10.5)
axL.legend(fontsize=8.5, loc='upper left'); axL.grid(alpha=0.15)

fig.suptitle('Triangle mitogen ramp — WITH vs WITHOUT H3K27me3 repression on CyclinD1', fontsize=13, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/sim_ramp_hysteresis_trace_control.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/sim_ramp_hysteresis_trace_control.pdf', bbox_inches='tight')
plt.close()
print(f'WITH   : entry={enW:.3f} exit={exW:.3f} divisions up={int((TW[dW]<=tt).sum())} down={int((TW[dW]>tt).sum())}')
print(f'WITHOUT: entry={enN:.3f} exit={exN:.3f} divisions up={int((TN[dN]<=tt).sum())} down={int((TN[dN]>tt).sum())}')
print('Saved sim_ramp_hysteresis_trace_control.png')
