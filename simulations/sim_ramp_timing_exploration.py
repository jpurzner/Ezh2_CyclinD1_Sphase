"""TIMING sweep (FINE GRID): how the H3K27me3→CyclinD1 kinetics (mark turnover del_mk, EZH2 half-life
kDeEZ) and REPLICATIVE-DILUTION frequency (cell-cycle period via mu) set the entry/exit HYSTERESIS under
mitogen ramps, as a function of ramp SPEED.

The mark is a slow integrator: on a fast ramp it lags the mitogen and a hysteresis opens between the
ramp-up (entry) and ramp-down (exit) Hh thresholds; on a slow (quasi-static) ramp the lag vanishes. We
grid ramp duration D (fast→slow) against each timing variable and read out the hysteresis gap
|entry−exit| (thresholds = SHH where cycle-averaged CyclinD1 crosses half its cycling value).

Panels (finely gridded heatmaps + one line panel):
  (A) gap over (ramp duration × mark turnover t½)     — del_mk
  (B) gap over (ramp duration × EZH2 t½)              — kDeEZ
  (C) gap over (ramp duration × cell-cycle period)    — mu (dilution frequency)
  (D) entry & exit threshold vs ramp duration (default timing) — the hysteresis loop, fine

Run:  ./venv/bin/python simulations/sim_ramp_timing_exploration.py [--fresh]  (~25 min cold; re-plots from cache)
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
CACHE = 'simulations/sim_ramp_timing_exploration_cache.npz'
PREEQ, CHUNK = 3000.0, 120.0
LN2 = np.log(2.0); hl = lambda r: LN2 / r / 60.0

M = build_model_v44(with_ezh2=True, with_hh=True)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)
DEFAULTS = dict(f0_mk=0.233, mu=0.0005, kTlEZ=0.004, del_mk=0.0007, kDeEZ=0.00015, HU=0.0, k_Cd_translation=0.801)
_R = te.loada(M); _R.integrator.setValue('relative_tolerance', 1e-6)
try: _R.integrator.setValue('maximum_num_steps', 100000)
except Exception: pass


def _set(ov):
    _R.reset()
    for k, v in GNP.items(): _R[k] = v
    for k, v in DEFAULTS.items(): _R[k] = v
    for k, v in ov.items(): _R[k] = v


def movavg(t, x, win=2000.0):
    dt = np.median(np.diff(t)) if len(t) > 1 else 1.0        # efficient box filter (~1.5 cell cycles)
    w = max(3, min(int(win / dt), max(3, len(x) // 2)))      # cap window at half the trace (mode='same' needs w<=len)
    return np.convolve(x, np.ones(w) / w, mode='same')


def ramp_thresh(ov, D, entry):
    _set(ov)
    shh0, shh1 = (0.05, 1.10) if entry else (1.10, 0.03)
    _R['SHH'] = shh0
    for atol in (1e-8, 1e-7, 1e-6):
        try:
            _R.integrator.setValue('absolute_tolerance', atol); _R.simulate(0, PREEQ, 100, selections=['time']); break
        except Exception:
            pass
    t = PREEQ
    shh_fn = lambda tt: shh0 + (shh1 - shh0) * min(max((tt - PREEQ) / D, 0.0), 1.0)
    T, SH, CD = [], [], []
    while t < PREEQ + D + 1500:
        _R['SHH'] = float(shh_fn(t)); ok = False
        for atol in (1e-8, 1e-7, 1e-6):
            try:
                _R.integrator.setValue('absolute_tolerance', atol)
                r = _R.simulate(t, t + CHUNK, max(20, int(CHUNK / 6)), selections=['time', 'SHH', 'Cd']); ok = True; break
            except Exception:
                pass
        if not ok:
            break
        T.append(r['time'][1:]); SH.append(r['SHH'][1:]); CD.append(r['Cd'][1:]); t += CHUNK
    if not T:
        return np.nan
    T = np.concatenate(T); SH = np.concatenate(SH); CD = np.concatenate(CD)
    m = T >= PREEQ; T, SH, CD = T[m], SH[m], CD[m]
    cav = movavg(T, CD)
    hi = SH >= 0.85 * SH.max()                                  # cycling plateau (high-SHH end of either ramp)
    ref = np.nanmedian(cav[hi]) if hi.any() else np.nanmax(cav)
    thr = 0.5 * ref

    def interp(k):                                              # SHH where cav crosses thr between samples k,k+1
        d = cav[k + 1] - cav[k]
        f = (thr - cav[k]) / d if abs(d) > 1e-9 else 0.0
        return float(SH[k] + f * (SH[k + 1] - SH[k]))

    if entry:                                                   # rising SHH → first up-crossing of thr
        cr = np.where((cav[:-1] < thr) & (cav[1:] >= thr))[0]
        return interp(cr[0]) if len(cr) else (float(SH[np.argmax(cav >= thr)]) if (cav >= thr).any() else np.nan)
    cr = np.where((cav[:-1] >= thr) & (cav[1:] < thr))[0]       # falling SHH → last down-crossing of thr
    return interp(cr[-1]) if len(cr) else (float(SH[len(SH) - 1 - np.argmax((cav >= thr)[::-1])]) if (cav >= thr).any() else np.nan)


NPD, NPV = 20, 16                                                # ramp-duration × timing-value grid (dense)
D_RANGE = np.round(np.geomspace(2500.0, 24000.0, NPD), 0)      # ~42 h .. 400 h (≥~2 cycles: cycle-avg valid)
VARS = {
    'del_mk': np.geomspace(0.0021, 0.00023, NPV),               # fast→slow mark turnover
    'kDeEZ':  np.geomspace(0.00045, 0.00005, NPV),              # fast→slow EZH2 half-life
    'mu':     np.geomspace(0.00095, 0.00028, NPV),              # short→long cycle (dilution freq)
}

if FRESH or not os.path.exists(CACHE):
    GAP = {}; ENTd = {}; EXId = {}
    for param, vals in VARS.items():
        ent = np.full((NPV, NPD), np.nan); exi = np.full((NPV, NPD), np.nan)
        for i, v in enumerate(vals):
            for j, D in enumerate(D_RANGE):
                ent[i, j] = ramp_thresh({param: v}, D, True)
                exi[i, j] = ramp_thresh({param: v}, D, False)
            print(f"  {param:8s} {i+1}/{NPV} v={v:.5f}  gap={np.round(np.abs(ent[i]-exi[i]),2)}", flush=True)
        GAP[param] = np.abs(ent - exi); ENTd[param] = ent; EXId[param] = exi
    np.savez(CACHE, D_RANGE=D_RANGE, delmk_vals=VARS['del_mk'], kdeez_vals=VARS['kDeEZ'], mu_vals=VARS['mu'],
             GAP_delmk=GAP['del_mk'], GAP_kdeez=GAP['kDeEZ'], GAP_mu=GAP['mu'],
             ENT_delmk=ENTd['del_mk'], EXI_delmk=EXId['del_mk'])
    print('cached ->', CACHE, flush=True)

z = np.load(CACHE)
D_RANGE = z['D_RANGE']; Dh = D_RANGE / 60.0
delmk_vals, kdeez_vals, mu_vals = z['delmk_vals'], z['kdeez_vals'], z['mu_vals']
GAP = {'del_mk': z['GAP_delmk'], 'kDeEZ': z['GAP_kdeez'], 'mu': z['GAP_mu']}


def _edges(c):
    lc = np.log10(np.asarray(c, float)); e = np.empty(len(c) + 1)
    e[1:-1] = 10 ** (0.5 * (lc[:-1] + lc[1:])); e[0] = 10 ** (lc[0] - 0.5 * (lc[1] - lc[0])); e[-1] = 10 ** (lc[-1] + 0.5 * (lc[-1] - lc[-2])); return e


XE = _edges(Dh)
vmax = max(np.nanmax(GAP['del_mk']), np.nanmax(GAP['kDeEZ']), np.nanmax(GAP['mu']))
fig, ax = plt.subplots(2, 2, figsize=(15, 11))

for a, key, yvals, ylab, ttl in [
    (ax[0, 0], 'del_mk', hl(delmk_vals), 'mark turnover t½ (h)', '(A) del_mk (mark turnover)'),
    (ax[0, 1], 'kDeEZ', hl(kdeez_vals), 'EZH2 t½ (h)', '(B) kDeEZ (EZH2 half-life)'),
    (ax[1, 0], 'mu', LN2 / mu_vals / 60.0, 'cell-cycle timescale ln2/mu (h)', '(C) mu (dilution frequency / cycle period)')]:
    YE = _edges(yvals)
    pm = a.pcolormesh(XE, YE, GAP[key], cmap='magma', shading='flat', vmin=0, vmax=vmax)
    a.set_xscale('log'); a.set_yscale('log')
    a.set_xlabel('ramp duration (h)  — faster ← → slower'); a.set_ylabel(ylab)
    a.set_title(ttl + '\nhysteresis |entry − exit| (SHH): biggest at FAST ramp + SLOW kinetics', fontweight='bold', fontsize=11)
    fig.colorbar(pm, ax=a, label='hysteresis gap (SHH)')

# (D) entry & exit vs ramp duration at default timing (nearest del_mk row to 0.0007)
di = int(np.argmin(np.abs(delmk_vals - 0.0007)))
ent = z['ENT_delmk'][di]; exi = z['EXI_delmk'][di]
a = ax[1, 1]
a.plot(Dh, ent, '-o', color='#1f8a4c', lw=2.4, ms=5, label='ENTRY (ramp-up)')
a.plot(Dh, exi, '-s', color='#c0392b', lw=2.4, ms=5, label='EXIT (ramp-down)')
a.fill_between(Dh, ent, exi, color='gold', alpha=0.25)
a.set_xscale('log'); a.set_xlabel('ramp duration (h)'); a.set_ylabel('Hh (SHH) threshold')
a.set_title('(D) Entry & exit threshold vs ramp speed (default timing)\nhysteresis widest at fast ramps, closing toward quasi-static', fontweight='bold', fontsize=11)
a.legend(fontsize=9); a.grid(alpha=0.15, which='both')

fig.suptitle('TIMING sweep (fine grid) — H3K27me3 kinetics & replicative-dilution frequency vs mitogen-ramp speed → entry/exit hysteresis',
             fontsize=12.5, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/sim_ramp_timing_exploration.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/sim_ramp_timing_exploration.pdf', bbox_inches='tight')
plt.close()
print('Saved sim_ramp_timing_exploration.png')
