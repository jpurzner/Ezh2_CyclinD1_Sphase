"""Parameter exploration of the EZH2->CyclinD1 circuit against two readouts:
  (1) the mitogen (Hh/SHH) THRESHOLD for the first division (cell-cycle entry), and
  (2) the number of divisions during a linear Hh WITHDRAWAL (cell-cycle exit).

Three knobs:
  - EZH2 mitogen-responsiveness  = kEZE2f   (default 0.022; how strongly EZH2 production tracks cycling/mitogen)
  - CyclinD1 inhibition from EZH2 = repression depth (1 - f0_mk)  (default f0_mk 0.233 -> depth 0.77; higher = more inhibition)
  - rate of Hh decline           = withdrawal ramp duration D (fast->slow); only affects readout (2)

Panels:
  A  Entry mitogen threshold over (EZH2 responsiveness x inhibition depth) — heatmap.
  B  Withdrawal divisions over (EZH2 responsiveness x inhibition depth), fixed decline rate — heatmap.
  C  Withdrawal divisions vs Hh decline rate, for a few EZH2 settings.

Run:  ./venv/bin/python simulations/fig_v44_param_exploration.py [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

FRESH = '--fresh' in sys.argv
CACHE = 'simulations/fig_v44_param_exploration_cache.npz'
M = build_model_v44(with_ezh2=True, with_hh=True)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)
KEZ0, F00 = 0.022, 0.233     # defaults


def _rr():
    r = te.loada(M); r.integrator.setValue('relative_tolerance', 1e-6); return r


def _set(r, kez, f0):
    r.reset()
    for k, v in GNP.items():
        r[k] = v
    r['kEZE2f'] = kez; r['f0_mk'] = f0


def divrate(r, shh, T=8000, settle=3500):
    r['SHH'] = shh
    for atol in (1e-9, 1e-8, 1e-7):
        try:
            r.integrator.setValue('absolute_tolerance', atol)
            res = r.simulate(0, T, int(T / 0.5), selections=['time', 'MPF']); m = res['time'] >= settle; t = res['time'][m]
            pk, _ = find_peaks(res['MPF'][m], prominence=0.15, distance=int(200 / (t[1] - t[0])))
            return len(pk) / ((t[-1] - t[0]) / 60) * 168.0
        except Exception:
            continue
    return np.nan


def entry_threshold(r, kez, f0, lo=0.02, hi=2.0, it=7):
    """bisection: smallest SHH that gives >0.5 div/168h (quasi-static)."""
    _set(r, kez, f0)
    if divrate(r, hi) < 0.5:
        return np.nan     # never cycles even at hi
    _set(r, kez, f0)
    if divrate(r, lo) > 0.5:
        return lo         # cycles even at lo
    for _ in range(it):
        mid = 0.5 * (lo + hi); _set(r, kez, f0)
        if divrate(r, mid) > 0.5:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def withdrawal_divs(r, kez, f0, D, s0=1.0, s1=0.03, t0=1500.0, chunk=100.0):
    """ramp SHH s0->s1 over D (min); count divisions after t0."""
    _set(r, kez, f0); r['SHH'] = s0
    shh = lambda t: s0 if t < t0 else (s1 if t > t0 + D else s0 + (s1 - s0) * (t - t0) / D)
    T_end = t0 + D + 2000; tt, mp = [], []; t = 0.0
    while t < T_end:
        r['SHH'] = float(shh(t))
        try:
            r.integrator.setValue('absolute_tolerance', 1e-8)
            res = r.simulate(t, t + chunk, 20, selections=['time', 'MPF']); tt.append(res['time'][1:]); mp.append(res['MPF'][1:])
        except Exception:
            try:
                r.integrator.setValue('absolute_tolerance', 1e-6)
                res = r.simulate(t, t + chunk, 20, selections=['time', 'MPF']); tt.append(res['time'][1:]); mp.append(res['MPF'][1:])
            except Exception:
                break
        t += chunk
    if not tt:
        return np.nan
    tt = np.concatenate(tt); mp = np.concatenate(mp)
    pk, _ = find_peaks(mp, prominence=0.15, distance=int(200 / (tt[1] - tt[0])))
    return int(np.sum(tt[pk] >= t0))


if FRESH or not os.path.exists(CACHE):
    KEZ = np.array([0.011, 0.0165, 0.022, 0.033, 0.044])
    F0 = np.array([0.05, 0.15, 0.233, 0.35, 0.5])          # inhibition: depth = 1 - f0
    DEPTH = 1 - F0
    r = _rr()
    print('A: entry threshold heatmap ...')
    ZA = np.full((len(F0), len(KEZ)), np.nan)
    for i, f0 in enumerate(F0):
        for j, kz in enumerate(KEZ):
            ZA[i, j] = entry_threshold(r, kz, f0)
        print(f'  f0={f0:.2f} row done')
    print('B: withdrawal divisions heatmap (D=150h) ...')
    Dfix = 9000.0
    ZB = np.full((len(F0), len(KEZ)), np.nan)
    for i, f0 in enumerate(F0):
        for j, kz in enumerate(KEZ):
            ZB[i, j] = withdrawal_divs(r, kz, f0, Dfix)
        print(f'  f0={f0:.2f} row done')
    print('C: withdrawal divisions vs decline rate ...')
    DRATES = np.array([3000.0, 6000.0, 9000.0, 15000.0, 24000.0])   # 50..400 h
    SETTINGS = [('no EZH2 (f0=1)', KEZ0, 1.0), ('default', KEZ0, F00), ('strong inhib (f0=0.1)', KEZ0, 0.1),
                ('high EZH2 resp (2x)', 0.044, F00)]
    ZC = np.full((len(SETTINGS), len(DRATES)), np.nan)
    for si, (lab, kz, f0) in enumerate(SETTINGS):
        for di, D in enumerate(DRATES):
            ZC[si, di] = withdrawal_divs(r, kz, f0, D)
        print(f'  setting "{lab}" done')
    np.savez(CACHE, KEZ=KEZ, F0=F0, DEPTH=DEPTH, ZA=ZA, ZB=ZB, DRATES=DRATES, ZC=ZC,
             SET_LAB=[s[0] for s in SETTINGS])
    print('  cached ->', CACHE)

z = np.load(CACHE, allow_pickle=True)
KEZ, F0, DEPTH, ZA, ZB, DRATES, ZC = z['KEZ'], z['F0'], z['DEPTH'], z['ZA'], z['ZB'], z['DRATES'], z['ZC']
SET_LAB = list(z['SET_LAB'])
fig, ax = plt.subplots(1, 3, figsize=(18, 5.4))


def heat(a, Z, title, clabel, cmap):
    pm = a.pcolormesh(KEZ, DEPTH, Z, shading='nearest', cmap=cmap)
    a.set_xlabel('EZH2 mitogen-responsiveness  (kEZE2f)'); a.set_ylabel('CyclinD1 inhibition depth  (1 − f0)')
    a.axvline(0.022, color='w', ls=':', lw=1); a.axhline(1 - 0.233, color='w', ls=':', lw=1)
    a.set_title(title, fontweight='bold', fontsize=11)
    fig.colorbar(pm, ax=a, label=clabel)
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            if Z[i, j] == Z[i, j]:
                a.text(KEZ[j], DEPTH[i], f'{Z[i, j]:.2f}' if Z[i, j] < 5 else f'{Z[i, j]:.0f}', ha='center', va='center', fontsize=6.5, color='w')


heat(ax[0], ZA, '(A) Mitogen threshold for FIRST DIVISION\nstronger/ more-responsive EZH2 → higher Hh needed', 'entry SHH threshold', 'viridis')
ax[0].text(0.033, 1 - 0.233, ' default', color='w', fontsize=7)
heat(ax[1], ZB, '(B) Divisions during Hh WITHDRAWAL (150 h)\nstronger/ more-responsive EZH2 → fewer divisions', 'withdrawal divisions', 'magma')

c = ax[2]
for si, lab in enumerate(SET_LAB):
    c.plot(DRATES / 60.0, ZC[si], '-o', lw=2.2, ms=4, label=lab)
c.set_xlabel('Hh decline time (h)  [fast → slow]'); c.set_ylabel('divisions during withdrawal')
c.set_title('(C) Withdrawal divisions vs Hh decline rate\nslower decline → more divisions; EZH2 caps them', fontweight='bold', fontsize=11)
c.legend(fontsize=8, loc='upper left'); c.grid(alpha=0.15)

fig.suptitle('Parameter exploration: EZH2 mitogen-responsiveness (kEZE2f) × CyclinD1 inhibition (1−f0) × Hh decline rate → entry threshold & withdrawal divisions',
             fontsize=11.5, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/fig_v44_param_exploration.png', dpi=160, bbox_inches='tight')
plt.savefig('simulations/fig_v44_param_exploration.pdf', bbox_inches='tight')
plt.close()
print('Saved: fig_v44_param_exploration.png / .pdf')
