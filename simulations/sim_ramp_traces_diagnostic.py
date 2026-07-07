"""DIAGNOSTIC: complete traces of the mitogen ramp (entry & withdrawal) across the parameter corners,
to see the mechanism directly (where CyclinD1 crosses, how the H3K27me3 mark accumulates/dilutes, when
the cell divides) rather than relying on a single threshold number.

4 conditions = corners of (H3K27me3 strength × cell-cycle duration); columns = entry ramp (SHH up) and
withdrawal ramp (SHH down). Each panel overlays SHH (input), CyclinD1 (Cd, the readout), the H3K27me3
mark (Mk), and division events (MPF peaks).

Run:  ./venv/bin/python simulations/sim_ramp_traces_diagnostic.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

D, PREEQ, CHUNK = 9000.0, 3000.0, 150.0
M = build_model_v44(with_ezh2=True, with_hh=True)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)
DEFAULTS = dict(f0_mk=0.233, mu=0.0005, kTlEZ=0.004, del_mk=0.0007, HU=0.0, k_Cd_translation=0.801)
_R = te.loada(M); _R.integrator.setValue('relative_tolerance', 1e-6)
try: _R.integrator.setValue('maximum_num_steps', 100000)
except Exception: pass


def ramp_trace(f0, mu, entry):
    _R.reset()
    for k, v in GNP.items(): _R[k] = v
    for k, v in DEFAULTS.items(): _R[k] = v
    _R['f0_mk'] = f0; _R['mu'] = mu
    shh0, shh1 = (0.05, 1.10) if entry else (1.10, 0.03)
    _R['SHH'] = shh0
    sel = ['time', 'SHH', 'Cd', 'Mk', 'MPF']
    t = 0.0
    for atol in (1e-8, 1e-7, 1e-6):
        try:
            _R.integrator.setValue('absolute_tolerance', atol); _R.simulate(0, PREEQ, 200, selections=sel); break
        except Exception:
            pass
    t = PREEQ
    shh_fn = lambda tt: shh0 + (shh1 - shh0) * min(max((tt - PREEQ) / D, 0.0), 1.0)
    T, SH, CD, MK, MP = [], [], [], [], []
    while t < PREEQ + D + 1500:
        _R['SHH'] = float(shh_fn(t)); ok = False
        for atol in (1e-8, 1e-7, 1e-6):
            try:
                _R.integrator.setValue('absolute_tolerance', atol)
                r = _R.simulate(t, t + CHUNK, max(30, int(CHUNK / 5)), selections=sel); ok = True; break
            except Exception:
                pass
        if not ok:
            break
        T.append(r['time'][1:]); SH.append(r['SHH'][1:]); CD.append(r['Cd'][1:]); MK.append(r['Mk'][1:]); MP.append(r['MPF'][1:])
        t += CHUNK
    if not T:
        return None
    return (np.concatenate(T) - PREEQ) / 60.0, np.concatenate(SH), np.concatenate(CD), np.concatenate(MK), np.concatenate(MP)


CONDS = [
    ('weak mark · short cycle',   0.60, 0.00090),
    ('strong mark · short cycle', 0.10, 0.00090),
    ('weak mark · long cycle',    0.60, 0.00030),
    ('strong mark · long cycle',  0.10, 0.00030),
]
THR = 2.0
fig, axes = plt.subplots(len(CONDS), 2, figsize=(15, 12), sharex='col')
for row, (name, f0, mu) in enumerate(CONDS):
    for col, entry in enumerate([True, False]):
        a = axes[row, col]
        tr = ramp_trace(f0, mu, entry)
        if tr is None:
            a.text(0.5, 0.5, 'integration failed', transform=a.transAxes, ha='center'); continue
        th, sh, cd, mk, mp = tr
        a.plot(th, cd, color='#b0402f', lw=1.4, label='CyclinD1 (Cd)')
        a.axhline(THR, color='#b0402f', ls=':', lw=0.9, alpha=0.6)
        a.set_ylim(0, max(6, np.nanmax(cd) * 1.05)); a.set_ylabel('CyclinD1', color='#b0402f', fontsize=9)
        a.tick_params(axis='y', colors='#b0402f')
        at = a.twinx()
        at.plot(th, sh, color='#888', lw=1.3, label='SHH (input)')
        at.plot(th, mk, color='#7a5aa6', lw=1.3, label='H3K27me3 (Mk)')
        at.set_ylim(0, 1.2); at.set_ylabel('SHH / Mk', fontsize=9)
        pk, _ = find_peaks(mp, prominence=0.15, distance=200)
        for p in pk:
            a.plot(th[p], a.get_ylim()[1] * 0.97, marker='v', color='#2ca25f', ms=6, clip_on=False)
        # entry/exit SHH annotation
        if len(pk):
            key = th[pk[0]] if entry else th[pk[-1]]
            skey = sh[pk[0]] if entry else sh[pk[-1]]
            a.axvline(key, color='#2ca25f', lw=1.4, alpha=0.6)
            a.text(key, a.get_ylim()[1] * 0.55, f'{"entry" if entry else "exit"}\nSHH {skey:.2f}', color='#1a7a43', fontsize=8, ha='center')
        if row == 0:
            a.set_title(('ENTRY ramp (SHH ↑)' if entry else 'WITHDRAWAL ramp (SHH ↓)'), fontweight='bold')
        if col == 0:
            a.text(-0.14, 0.5, name, transform=a.transAxes, rotation=90, va='center', fontweight='bold', fontsize=10)
        if row == len(CONDS) - 1:
            a.set_xlabel('time along ramp (h)')
        if row == 0 and col == 0:
            l1, lab1 = a.get_legend_handles_labels(); l2, lab2 = at.get_legend_handles_labels()
            a.legend(l1 + l2 + [plt.Line2D([], [], marker='v', color='#2ca25f', ls='none')], lab1 + lab2 + ['division'], fontsize=7.5, loc='upper left')

fig.suptitle('Mitogen-ramp TRACE DIAGNOSTIC — CyclinD1 (red), SHH input (grey), H3K27me3 mark (purple), divisions (▼)',
             fontsize=13, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/sim_ramp_traces_diagnostic.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/sim_ramp_traces_diagnostic.pdf', bbox_inches='tight')
plt.close()
print('Saved sim_ramp_traces_diagnostic.png')
