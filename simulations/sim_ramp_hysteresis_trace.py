"""Worked example: a single continuous TRIANGLE mitogen ramp (SHH up, then down) showing the cell-cycle
trace and WHY the CyclinD1-vs-SHH loop opens. One run, default timing, fast-ish ramp.

The story in one trace:
 - low SHH (arrest): cell not cycling; no replicative dilution → the H3K27me3 mark ACCUMULATES.
 - ramp UP: SHH rises, but the accumulated mark still represses CyclinD1 → the cell re-enters LATE
   (entry SHH high). Once cycling, each division HALVES the mark (sawtooth down).
 - at the top the mark is low (diluted by cycling) → CyclinD1 high.
 - ramp DOWN: SHH falls, but the mark is still low → CyclinD1 stays high and the cell keeps cycling to
   LOWER SHH than it needed to enter (exit SHH low). => hysteresis.

Panels:
 (A) SHH(t) triangle              (B) CyclinD1(t) sawtooth + H3K27me3 mark(t)
 (C) EZH2(t) + mass(t) w/ division ticks   (D) the CyclinD1-vs-SHH loop traced out (up=green, down=red)

Run:  ./venv/bin/python simulations/sim_ramp_hysteresis_trace.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

LN2 = np.log(2.0)
M = build_model_v44(with_ezh2=True, with_hh=True)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)
DEFAULTS = dict(f0_mk=0.233, mu=0.0005, kTlEZ=0.004, del_mk=0.0007, kDeEZ=0.00015, HU=0.0, k_Cd_translation=0.801)
_R = te.loada(M); _R.integrator.setValue('relative_tolerance', 1e-6)
try: _R.integrator.setValue('maximum_num_steps', 100000)
except Exception: pass

for k, v in GNP.items(): _R[k] = v
for k, v in DEFAULTS.items(): _R[k] = v

D_LEG = 9000.0            # 150 h per leg (fast enough for a clear loop, ~6-8 divisions/leg)
SHH_LO, SHH_HI = 0.05, 1.10
PREEQ = float(min(20000.0, max(3000.0, 4.0 * LN2 / DEFAULTS['kDeEZ'])))   # load the mark before we start
CHUNK = 60.0
SEL = ['time', 'SHH', 'Cd', 'Mk', 'EZH2', 'mass', 'Cdc20']


def shh_of(t):
    x = t - PREEQ
    if x <= 0:
        return SHH_LO
    if x <= D_LEG:
        return SHH_LO + (SHH_HI - SHH_LO) * (x / D_LEG)          # up
    y = x - D_LEG
    return SHH_HI + (SHH_LO - 0.02 - SHH_HI) * min(y / D_LEG, 1.0)  # down (to just below start)


# preequilibrate arrested at low SHH
_R['SHH'] = SHH_LO
for atol in (1e-8, 1e-7, 1e-6):
    try:
        _R.integrator.setValue('absolute_tolerance', atol); _R.simulate(0, PREEQ, 50, selections=['time']); break
    except Exception:
        pass

rows = {s: [] for s in SEL}
t = PREEQ; TEND = PREEQ + 2 * D_LEG + 200
while t < TEND:
    _R['SHH'] = float(shh_of(t))
    for atol in (1e-8, 1e-7, 1e-6):
        try:
            _R.integrator.setValue('absolute_tolerance', atol)
            r = _R.simulate(t, t + CHUNK, 12, selections=SEL); break
        except Exception:
            r = None
    if r is None:
        break
    for s in SEL:
        rows[s].append(r[s][1:])
    t += CHUNK

for s in SEL:
    rows[s] = np.concatenate(rows[s])
T = (rows['time'] - PREEQ) / 60.0                                # h since ramp start
SH, CD, MK, EZ, MASS, C20 = rows['SHH'], rows['Cd'], rows['Mk'], rows['EZH2'], rows['mass'], rows['Cdc20']
t_top = D_LEG / 60.0                                             # leg boundary (h)
up = T <= t_top


def movavg(x, win=1400.0):
    dt = np.median(np.diff(rows['time']))
    w = max(3, min(int(win / dt), max(3, len(x) // 2)))
    return np.convolve(x, np.ones(w) / w, mode='same')


# divisions = mass drops (halving)
dmass = np.diff(MASS)
divs = np.where(dmass < -0.15 * MASS[:-1])[0]
div_t = T[divs]
n_up = int((div_t <= t_top).sum()); n_dn = int((div_t > t_top).sum())
shh_entry = SH[divs[0]] if len(divs) else np.nan
dn_divs = div_t[div_t > t_top]
shh_exit = SH[divs[div_t > t_top][-1]] if len(dn_divs) else np.nan

fig, ax = plt.subplots(3, 2, figsize=(16, 11), gridspec_kw={'width_ratios': [1.55, 1]})
gs = ax[0, 1].get_gridspec()
for a in ax[:, 1]:
    a.remove()
axD = fig.add_subplot(gs[:, 1])


def shade(a):
    a.axvspan(0, t_top, color='#1f8a4c', alpha=0.05); a.axvspan(t_top, T.max(), color='#c0392b', alpha=0.05)
    a.axvline(t_top, color='k', lw=0.8, ls=':')


# (A) SHH
a = ax[0, 0]; shade(a)
a.plot(T, SH, color='#333', lw=2)
a.set_ylabel('SHH (mitogen)'); a.set_title('(A) Triangle mitogen ramp — UP (green) then DOWN (red)', fontweight='bold', fontsize=11)
a.text(t_top * 0.5, SHH_HI * 0.95, 'ramp UP', color='#1f8a4c', ha='center', fontweight='bold')
a.text(t_top * 1.5, SHH_HI * 0.95, 'ramp DOWN', color='#c0392b', ha='center', fontweight='bold')
a.grid(alpha=0.15)

# (B) CyclinD1 sawtooth + mark
a = ax[1, 0]; shade(a)
a.plot(T, CD, color='#1f6fb0', lw=1.3, label='CyclinD1 (Cd) — cycling sawtooth')
a.plot(T, movavg(CD), color='#0b3d61', lw=2.2, label='CyclinD1 cycle-averaged')
a.set_ylabel('CyclinD1', color='#0b3d61'); a.tick_params(axis='y', labelcolor='#0b3d61')
for d in div_t:
    a.axvline(d, color='#1f6fb0', lw=0.5, alpha=0.25)
a2 = a.twinx()
a2.plot(T, MK, color='#8e44ad', lw=2.2, label='H3K27me3 mark (Mk)')
a2.set_ylabel('H3K27me3 mark', color='#8e44ad'); a2.tick_params(axis='y', labelcolor='#8e44ad')
a.set_title('(B) CyclinD1 sawtooth + H3K27me3 mark: mark accumulates in arrest, diluted by divisions, lags SHH', fontweight='bold', fontsize=10.5)
h1, l1 = a.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
a.legend(h1 + h2, l1 + l2, fontsize=8.5, loc='upper left')
a.grid(alpha=0.12)

# (C) EZH2 + mass + division ticks
a = ax[2, 0]; shade(a)
a.plot(T, EZ, color='#c47f17', lw=2.2, label='EZH2 (writer)')
a.set_ylabel('EZH2', color='#c47f17'); a.tick_params(axis='y', labelcolor='#c47f17')
a.set_xlabel('time since ramp start (h)')
a3 = a.twinx()
a3.plot(T, MASS, color='#666', lw=1.1, label='cell mass (halves at division)')
a3.plot(div_t, np.full_like(div_t, MASS.max() * 1.05), 'v', color='k', ms=6, label='division')
a3.set_ylabel('mass', color='#666'); a3.tick_params(axis='y', labelcolor='#666')
a.set_title(f'(C) EZH2 writer + divisions — {n_up} divisions on UP-leg, {n_dn} on DOWN-leg', fontweight='bold', fontsize=10.5)
h1, l1 = a.get_legend_handles_labels(); h2, l2 = a3.get_legend_handles_labels()
a.legend(h1 + h2, l1 + l2, fontsize=8.5, loc='upper right')
a.grid(alpha=0.12)

# (D) the loop traced by this single run
cav = movavg(CD)
axD.plot(SH[up], cav[up], color='#1f8a4c', lw=2.6, label='ramp UP (mark low → CyclinD1 high)')
axD.plot(SH[~up], cav[~up], color='#c0392b', lw=2.6, label='ramp DOWN (mark high → CyclinD1 low)')
axD.plot(SH[up], CD[up], color='#1f8a4c', lw=0.6, alpha=0.3)
axD.plot(SH[~up], CD[~up], color='#c0392b', lw=0.6, alpha=0.3)
if np.isfinite(shh_entry):
    axD.axvline(shh_entry, color='#1f8a4c', ls='--', lw=1.2, alpha=0.7)
    axD.text(shh_entry, axD.get_ylim()[1] * 0.05, f' entry\n SHH={shh_entry:.2f}', color='#1f8a4c', fontsize=8.5)
if np.isfinite(shh_exit):
    axD.axvline(shh_exit, color='#c0392b', ls='--', lw=1.2, alpha=0.7)
    axD.text(shh_exit, axD.get_ylim()[1] * 0.30, f'exit\nSHH={shh_exit:.2f} ', color='#c0392b', fontsize=8.5, ha='right')
axD.set_xlabel('mitogen SHH'); axD.set_ylabel('CyclinD1')
axD.set_title('(D) The hysteresis loop traced by this ONE run\n(faint = raw sawtooth, bold = cycle-averaged)', fontweight='bold', fontsize=11)
axD.legend(fontsize=9, loc='upper left'); axD.grid(alpha=0.15)

fig.suptitle(f'Worked cell-cycle trace through a triangle mitogen ramp ({D_LEG/60:.0f} h/leg, default timing) — how the CyclinD1–SHH hysteresis loop is drawn',
             fontsize=13, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/sim_ramp_hysteresis_trace.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/sim_ramp_hysteresis_trace.pdf', bbox_inches='tight')
plt.close()
print(f'entry SHH={shh_entry:.3f}  exit SHH={shh_exit:.3f}  divisions up={n_up} down={n_dn}')
print('Saved sim_ramp_hysteresis_trace.png')
