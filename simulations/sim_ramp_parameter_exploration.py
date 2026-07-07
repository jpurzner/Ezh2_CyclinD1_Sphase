"""How H3K27me3 repression STRENGTH, cell-cycle DURATION, and other features shift the mitogen (Hh)
response for cell-cycle ENTRY (rising Hh) and WITHDRAWAL (falling Hh).

Robust metrics (single representative GNP cell), reusing the proven approach from sweep_param_exploration:
  ENTRY threshold  = smallest SHH giving sustained cycling (bisection on settled division rate).
  WITHDRAWAL divisions = # divisions during a linear SHH ramp-down from cycling (how long the cell coasts).

Panels:
  (A) entry threshold over (H3K27me3 strength [1-f0_mk] × cell-cycle duration [growth rate mu])
  (B) withdrawal divisions over the same plane
  (C) ENTRY-threshold sensitivity tornado (each feature low→high)
  (D) WITHDRAWAL-divisions sensitivity tornado

Run:  ./venv/bin/python simulations/sim_ramp_parameter_exploration.py [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

FRESH = '--fresh' in sys.argv
CACHE = 'simulations/sim_ramp_parameter_exploration_cache.npz'

M = build_model_v44(with_ezh2=True, with_hh=True)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)
_R = te.loada(M); _R.integrator.setValue('relative_tolerance', 1e-6)
try: _R.integrator.setValue('maximum_num_steps', 80000)
except Exception: pass


# roadrunner reset() resets SPECIES only, NOT parameters — so every call must re-establish a full
# baseline for every parameter any sweep touches, otherwise values leak between calls.
DEFAULTS = dict(f0_mk=0.233, mu=0.0005, kTlEZ=0.004, del_mk=0.0007, HU=0.0, p18=0.464, k_Cd_translation=0.801)


def _set(ov):
    _R.reset()
    for k, v in GNP.items(): _R[k] = v
    for k, v in DEFAULTS.items(): _R[k] = v          # clean baseline (defeats parameter leakage)
    for k, v in ov.items(): _R[k] = v


def divrate(ov, shh, T=8000, settle=3500):
    _set(ov); _R['SHH'] = float(shh)
    for atol in (1e-8, 1e-7, 1e-6):
        try:
            _R.integrator.setValue('absolute_tolerance', atol)
            r = _R.simulate(0, T, int(T / 0.5), selections=['time', 'MPF']); m = r['time'] >= settle
            t = r['time'][m]; pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=int(200 / (t[1] - t[0])))
            return len(pk) / ((t[-1] - t[0]) / 60.0) * 168.0
        except Exception:
            continue
    return np.nan


def entry_threshold(ov, lo=0.02, hi=1.5, it=7):
    if not (divrate(ov, hi) > 0.5):
        return np.nan
    if divrate(ov, lo) > 0.5:
        return lo
    for _ in range(it):
        mid = 0.5 * (lo + hi)
        if divrate(ov, mid) > 0.5: hi = mid
        else: lo = mid
    return 0.5 * (lo + hi)


def movavg(t, x, win=1320):
    return np.array([x[(t >= T - win / 2) & (t <= T + win / 2)].mean() for T in t])


def exit_threshold(ov, D=9000.0, s0=1.0, s1=0.03, t0=1500.0):
    """EXIT Hh = the SHH at which cycle-averaged CyclinD1 drops below HALF its cycling value, during a
    linear SHH ramp-down from cycling. Scale-invariant (relative to each cell's own cycling CyclinD1),
    so it isolates the mitogen threshold rather than confounding with cycle speed."""
    _set(ov); _R['SHH'] = s0
    shh = lambda tt: s0 if tt < t0 else (s1 if tt > t0 + D else s0 + (s1 - s0) * (tt - t0) / D)
    chunk = max(150.0, D / 90.0); T_end = t0 + D + 800
    T, CD, SH = [], [], []; t = 0.0; fails = 0
    while t < T_end:
        _R['SHH'] = float(shh(t)); ok = False
        for atol in (1e-8, 1e-7, 1e-6):
            try:
                _R.integrator.setValue('absolute_tolerance', atol)
                res = _R.simulate(t, t + chunk, max(40, int(chunk / 4)), selections=['time', 'Cd', 'SHH']); ok = True; break
            except Exception:
                pass
        if not ok:
            fails += 1
            if fails > 3: break
            t += chunk; continue
        T.append(res['time'][1:]); CD.append(res['Cd'][1:]); SH.append(res['SHH'][1:]); t += chunk
    if not T:
        return np.nan
    T = np.concatenate(T); CD = np.concatenate(CD); SH = np.concatenate(SH)
    cav = movavg(T, CD)
    cyc = np.median(cav[(T >= t0 - 800) & (T < t0)])          # cycling CyclinD1 before withdrawal
    above = np.where((cav >= 0.5 * cyc) & (T >= t0))[0]        # last SHH with Cd above half-cycling
    return float(SH[above[-1]]) if len(above) else s1


FEATURES = {
    'H3K27me3 strength':   ('f0_mk', 0.70, 0.05),
    'cell-cycle duration': ('mu', 0.00095, 0.00030),
    'EZH2 level':          ('kTlEZ', 0.002, 0.008),
    'mark turnover':       ('del_mk', 0.0021, 0.00035),
    'S-phase length (HU)': ('HU', 0.0, 1.0),
    'CKI brake (p18)':     ('p18', 0.23, 1.40),
    'CyclinD1 drive':      ('k_Cd_translation', 0.40, 1.60),
}

if FRESH or not os.path.exists(CACHE):
    F0 = np.round(np.linspace(0.70, 0.05, 6), 3)
    MU = np.round(np.geomspace(0.00030, 0.00095, 6), 6)
    ENT = np.full((len(MU), len(F0)), np.nan); EXI = np.full((len(MU), len(F0)), np.nan)
    for i, mu in enumerate(MU):
        for j, f0 in enumerate(F0):
            ENT[i, j] = entry_threshold({'f0_mk': f0, 'mu': mu})
            EXI[i, j] = exit_threshold({'f0_mk': f0, 'mu': mu})
        print(f"  [plane] mu={mu:.5f}  entry={np.round(ENT[i],2)}  exit={np.round(EXI[i],2)}", flush=True)
    names = list(FEATURES); ent_lo = []; ent_hi = []; exi_lo = []; exi_hi = []
    for nm in names:
        p, lo, hi = FEATURES[nm]
        ent_lo.append(entry_threshold({p: lo})); ent_hi.append(entry_threshold({p: hi}))
        exi_lo.append(exit_threshold({p: lo})); exi_hi.append(exit_threshold({p: hi}))
        print(f"  [tornado] {nm:22s} entry {np.round(ent_lo[-1],2)}->{np.round(ent_hi[-1],2)}  exit {np.round(exi_lo[-1],2)}->{np.round(exi_hi[-1],2)}", flush=True)
    np.savez(CACHE, F0=F0, MU=MU, ENT=ENT, EXI=EXI, names=names,
             ent_lo=ent_lo, ent_hi=ent_hi, exi_lo=exi_lo, exi_hi=exi_hi)
    print('cached ->', CACHE, flush=True)

z = np.load(CACHE, allow_pickle=True)
F0, MU, ENT, EXI = z['F0'], z['MU'], z['ENT'], z['EXI']
names = list(z['names']); ent_lo, ent_hi = z['ent_lo'].astype(float), z['ent_hi'].astype(float)
exi_lo, exi_hi = z['exi_lo'].astype(float), z['exi_hi'].astype(float)
DEPTH = 1 - F0; PERIOD_LABEL = [f'{np.log(2)/mu/60:.0f}' for mu in MU]   # ~ln2/mu as a duration proxy (h)

fig, ax = plt.subplots(2, 2, figsize=(15, 11))

# (A) entry threshold heatmap
a = ax[0, 0]
pm = a.pcolormesh(DEPTH, np.arange(len(MU)), np.ma.masked_invalid(ENT), cmap='viridis', shading='nearest')
a.set_yticks(np.arange(len(MU))); a.set_yticklabels([f'{m:.5f}\n(~{p}h)' for m, p in zip(MU, PERIOD_LABEL)], fontsize=7)
a.set_xlabel('H3K27me3 strength (1 − f0_mk)'); a.set_ylabel('growth rate mu  (↑ = shorter cycle)')
a.set_title('(A) ENTRY Hh threshold\nstronger mark + longer cycle → higher threshold', fontweight='bold', fontsize=11)
fig.colorbar(pm, ax=a, label='entry SHH threshold')
for i in range(ENT.shape[0]):
    for j in range(ENT.shape[1]):
        if np.isfinite(ENT[i, j]): a.text(DEPTH[j], i, f'{ENT[i,j]:.2f}', ha='center', va='center', fontsize=7.5, color='w')

# (B) exit threshold heatmap (same units/colormap as entry)
a = ax[0, 1]
pm = a.pcolormesh(DEPTH, np.arange(len(MU)), np.ma.masked_invalid(EXI), cmap='viridis', shading='nearest')
a.set_yticks(np.arange(len(MU))); a.set_yticklabels([f'{m:.5f}' for m in MU], fontsize=7)
a.set_xlabel('H3K27me3 strength (1 − f0_mk)'); a.set_ylabel('growth rate mu')
a.set_title('(B) EXIT Hh threshold (on withdrawal)\nstronger mark + longer cycle → exit at higher Hh', fontweight='bold', fontsize=11)
fig.colorbar(pm, ax=a, label='exit SHH threshold')
for i in range(EXI.shape[0]):
    for j in range(EXI.shape[1]):
        if np.isfinite(EXI[i, j]): a.text(DEPTH[j], i, f'{EXI[i,j]:.2f}', ha='center', va='center', fontsize=7.5, color='w')

# tornados
for a, lo, hi, ttl, unit in [(ax[1, 0], ent_lo, ent_hi, '(C) ENTRY-threshold sensitivity', 'entry SHH'),
                             (ax[1, 1], exi_lo, exi_hi, '(D) EXIT-threshold sensitivity', 'exit SHH')]:
    lo = np.nan_to_num(lo, nan=np.nanmax(np.concatenate([lo, hi])))
    rng = np.abs(hi - lo); order = np.argsort(rng)
    for k, idx in enumerate(order):
        a.plot([lo[idx], hi[idx]], [k, k], color='#999', lw=1.6, zorder=1)
        a.plot(lo[idx], k, 'o', color='#4a90c2', ms=9, zorder=2)
        a.plot(hi[idx], k, 'o', color='#c0392b', ms=9, zorder=2)
    a.set_yticks(range(len(names))); a.set_yticklabels([names[i] for i in order], fontsize=9)
    a.set_xlabel(unit); a.set_title(ttl + '   (● low feature   ● high feature)', fontweight='bold', fontsize=11)
    a.grid(alpha=0.15, axis='x')

fig.suptitle('What shifts the Hh response for ENTRY & WITHDRAWAL: H3K27me3 strength, cell-cycle duration, and other features',
             fontsize=13, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/sim_ramp_parameter_exploration.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/sim_ramp_parameter_exploration.pdf', bbox_inches='tight')
plt.close()
print('Saved sim_ramp_parameter_exploration.png')
