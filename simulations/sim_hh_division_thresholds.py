"""Hh threshold for cell division — REPRESSED (f0_mk=0.233) vs NOT repressed (f0_mk=1.0), two scenarios:

  Figure 1  RAMP UP     : start ARRESTED at low Hh (mark loaded by arrest), ramp Hh up. Mark the Hh at
                          which cells START dividing (entry threshold).
  Figure 2  WITHDRAWAL  : start ESTABLISHED-CYCLING at high Hh, withdraw Hh. Mark the Hh at which cells
                          STOP dividing (exit threshold).

Each figure shows only: the CyclinD1 context curve, the individual division events (rug), and the Hh
threshold for each condition — so the repressed-vs-not shift is the whole story.

Run:  ./venv/bin/python simulations/sim_hh_division_thresholds.py [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

FRESH = '--fresh' in sys.argv
CACHE = 'simulations/sim_hh_division_thresholds_cache.npz'
LN2 = np.log(2.0)
M = build_model_v44(with_ezh2=True, with_hh=True)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)
DEFAULTS = dict(f0_mk=0.233, mu=0.0005, kTlEZ=0.004, del_mk=0.0007, kDeEZ=0.00015, HU=0.0, k_Cd_translation=0.801)
D_RAMP, SHH_LO, SHH_HI, CHUNK = 9000.0, 0.05, 1.10, 60.0
PREEQ = float(min(20000.0, max(3000.0, 4.0 * LN2 / DEFAULTS['kDeEZ'])))
SEL = ['time', 'SHH', 'Cd', 'mass']
_R = te.loada(M); _R.integrator.setValue('relative_tolerance', 1e-6)
try: _R.integrator.setValue('maximum_num_steps', 100000)
except Exception: pass


def ramp(f0, entry):
    _R.reset()
    for k, v in GNP.items(): _R[k] = v
    for k, v in DEFAULTS.items(): _R[k] = v
    _R['f0_mk'] = f0
    shh0, shh1 = (SHH_LO, SHH_HI) if entry else (SHH_HI, 0.02)
    _R['SHH'] = shh0
    for atol in (1e-8, 1e-7, 1e-6):
        try:
            _R.integrator.setValue('absolute_tolerance', atol); _R.simulate(0, PREEQ, 50, selections=['time']); break
        except Exception: pass
    t = PREEQ; rows = {s: [] for s in SEL}
    shh_fn = lambda tt: shh0 + (shh1 - shh0) * min(max((tt - PREEQ) / D_RAMP, 0.0), 1.0)
    while t < PREEQ + D_RAMP + 1500:
        _R['SHH'] = float(shh_fn(t)); r = None
        for atol in (1e-8, 1e-7, 1e-6):
            try:
                _R.integrator.setValue('absolute_tolerance', atol); r = _R.simulate(t, t + CHUNK, 12, selections=SEL); break
            except Exception: r = None
        if r is None: break
        for s in SEL: rows[s].append(r[s][1:])
        t += CHUNK
    for s in SEL: rows[s] = np.concatenate(rows[s])
    m = rows['time'] >= PREEQ
    return {s: rows[s][m] for s in SEL}


if FRESH or not os.path.exists(CACHE):
    data = {}
    for tag, entry in [('up', True), ('dn', False)]:
        for cond, f0 in [('W', 0.233), ('N', 1.0)]:
            r = ramp(f0, entry)
            for s in SEL: data[f'{tag}_{cond}_{s}'] = r[s]
            print(f'  {tag} {cond} done', flush=True)
    np.savez(CACHE, **data); print('cached ->', CACHE, flush=True)

z = np.load(CACHE)


def movavg(x, win=1400.0):
    w = max(3, int(win / CHUNK / (1 / 12.0)))            # ~cycle window; SEL sampled 12/chunk
    w = max(3, min(w, max(3, len(x) // 2)))
    return np.convolve(x, np.ones(w) / w, mode='same')


def divisions_shh(mass, shh):
    idx = np.where(np.diff(mass) < -0.15 * mass[:-1])[0]
    return shh[idx]


def get(tag, cond):
    SH = z[f'{tag}_{cond}_SHH']; CD = z[f'{tag}_{cond}_Cd']; MS = z[f'{tag}_{cond}_mass']
    dv = divisions_shh(MS, SH)
    return SH, movavg(CD), dv


def make(tag, entry, fname, ttl):
    SHW, cdW, dvW = get(tag, 'W'); SHN, cdN, dvN = get(tag, 'N')
    thrW = np.nanmin(dvW) if len(dvW) else np.nan                # lowest-Hh division = the threshold (start=up / stop=down)
    thrN = np.nanmin(dvN) if len(dvN) else np.nan
    fig, ax = plt.subplots(figsize=(11, 6))
    # CyclinD1 context
    ax.plot(SHW, cdW, color='#0b3d61', lw=2.2, label='CyclinD1 — REPRESSED (mark on)')
    ax.plot(SHN, cdN, color='#c0392b', lw=2.2, label='CyclinD1 — NOT repressed (mark off)')
    ymax = max(np.nanmax(cdW), np.nanmax(cdN)) * 1.12
    ax.set_ylim(0, ymax)
    # division rug (each event at its Hh)
    yW, yN = ymax * 0.97, ymax * 0.90
    ax.plot(dvW, np.full_like(dvW, yW), '|', color='#0b3d61', ms=16, mew=2.2)
    ax.plot(dvN, np.full_like(dvN, yN), '|', color='#c0392b', ms=16, mew=2.2)
    ax.text(1.11, yW, ' divisions (repressed)', color='#0b3d61', va='center', fontsize=8.5)
    ax.text(1.11, yN, ' divisions (not repressed)', color='#c0392b', va='center', fontsize=8.5)
    # threshold lines + shaded shift
    verb = 'START' if entry else 'STOP'
    for thr, c in [(thrW, '#0b3d61'), (thrN, '#c0392b')]:
        if np.isfinite(thr): ax.axvline(thr, color=c, ls='--', lw=1.8, alpha=0.8)
    if np.isfinite(thrW) and np.isfinite(thrN):
        ax.axvspan(min(thrW, thrN), max(thrW, thrN), color='gold', alpha=0.25)
        ax.text((thrW + thrN) / 2, ymax * 0.66, f'ΔHh\n{abs(thrW - thrN):.2f}', ha='center', fontsize=10, fontweight='bold')
    ax.text(thrW, -ymax * 0.07, f'repressed\nHh={thrW:.2f}', color='#0b3d61', ha='center', fontsize=8.5, fontweight='bold')
    ax.text(thrN, -ymax * 0.14, f'not repressed\nHh={thrN:.2f}', color='#c0392b', ha='center', fontsize=8.5, fontweight='bold')
    ax.set_xlabel('Hedgehog (SHH) — rising →' if entry else '← Hedgehog (SHH) withdrawn')
    ax.set_ylabel('CyclinD1 (cycle-avg)')
    ax.set_title(f'{ttl}\ncells {verb} dividing at Hh = {thrW:.2f} (repressed) vs {thrN:.2f} (not repressed)', fontweight='bold', fontsize=12)
    ax.set_xlim(0, 1.20); ax.legend(fontsize=9, loc=('lower right' if entry else 'lower left')); ax.grid(alpha=0.15)
    if not entry: ax.invert_xaxis()                             # withdrawal reads high→low left-to-right
    plt.tight_layout()
    plt.savefig(f'simulations/{fname}.png', dpi=150, bbox_inches='tight')
    plt.savefig(f'simulations/{fname}.pdf', bbox_inches='tight')
    plt.close()
    print(f'{fname}: repressed Hh_thr={thrW:.3f}  not-repressed Hh_thr={thrN:.3f}  (n_div W={len(dvW)} N={len(dvN)})')


make('up', True, 'sim_hh_entry_threshold', 'RAMP UP — Hh threshold to START dividing (entry)')
make('dn', False, 'sim_hh_withdrawal_threshold', 'Hh WITHDRAWAL — Hh threshold to STOP dividing (exit)')
print('done')
