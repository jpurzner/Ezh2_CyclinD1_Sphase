"""CONTROL for sim_ramp_hysteresis_area: is the CyclinD1-vs-SHH hysteresis loop actually CAUSED by the
H3K27me3 mark, or would the cell-cycle machinery produce it anyway?

Identical model, two conditions at DEFAULT timing:
  WITH mark   : f0_mk = 0.233 (default; mark represses CyclinD1)
  NO repression: f0_mk = 1.00  (mark still accumulates/dilutes & EZH2 runs, but exerts ZERO effect on
                 CyclinD1 transcription — R = f0_mk + (1-f0_mk)/(1+(Mk/K)^n) = 1 for all Mk).
Only the mark→CyclinD1 coupling differs. Residual loop in the NO-repression arm = whatever the Rb-E2f /
Skp2-p27 commitment switch contributes on its own.

Panels:
  (A) the loops themselves at a FAST ramp — WITH (green/red, gold fill) vs NO repression (grey, hatched)
  (B) normalized loop AREA vs ramp duration for both conditions — the mark's contribution across speeds

Run:  ./venv/bin/python simulations/sim_ramp_hysteresis_mark_control.py [--fresh]
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
CACHE = 'simulations/sim_ramp_hysteresis_mark_control_cache.npz'
PREEQ, CHUNK = 3000.0, 120.0
LN2 = np.log(2.0)

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
    dt = np.median(np.diff(t)) if len(t) > 1 else 1.0
    w = max(3, min(int(win / dt), max(3, len(x) // 2)))
    return np.convolve(x, np.ones(w) / w, mode='same')


def _preeq_for(ov):
    dm = ov.get('del_mk', DEFAULTS['del_mk']); ke = ov.get('kDeEZ', DEFAULTS['kDeEZ'])
    tau = max(LN2 / dm, LN2 / ke)
    return float(min(20000.0, max(PREEQ, 4.0 * tau)))


def ramp_curve(ov, D, entry):
    _set(ov)
    preeq = _preeq_for(ov)
    shh0, shh1 = (0.05, 1.10) if entry else (1.10, 0.03)
    _R['SHH'] = shh0
    for atol in (1e-8, 1e-7, 1e-6):
        try:
            _R.integrator.setValue('absolute_tolerance', atol); _R.simulate(0, preeq, 100, selections=['time']); break
        except Exception:
            pass
    t = preeq
    shh_fn = lambda tt: shh0 + (shh1 - shh0) * min(max((tt - preeq) / D, 0.0), 1.0)
    T, SH, CD = [], [], []
    while t < preeq + D + 1500:
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
        return None
    T = np.concatenate(T); SH = np.concatenate(SH); CD = np.concatenate(CD)
    m = T >= preeq; T, SH, CD = T[m], SH[m], CD[m]
    cav = movavg(T, CD)
    o = np.argsort(SH)
    return SH[o], cav[o]


GRID = np.linspace(0.10, 1.00, 60)


def loop(ov, D):
    """returns (up_sh, up_cd, dn_sh, dn_cd, normalized_area)."""
    up = ramp_curve(ov, D, True); dn = ramp_curve(ov, D, False)
    if up is None or dn is None:
        return None
    cu = np.interp(GRID, up[0], up[1]); cd = np.interp(GRID, dn[0], dn[1])
    ref = 0.5 * (np.nanmax(up[1]) + np.nanmax(dn[1]))
    area = float(np.trapz(cu - cd, GRID) / (ref + 1e-9))
    return up, dn, area


COND = {'with': {}, 'nomark': {'f0_mk': 1.0}}
D_LIST = np.round(np.geomspace(2500.0, 24000.0, 8), 0)     # ~42 h .. 400 h
DFAST = D_LIST[1]

if FRESH or not os.path.exists(CACHE):
    areas = {c: np.full(len(D_LIST), np.nan) for c in COND}
    loops = {}
    for c, ov in COND.items():
        for j, D in enumerate(D_LIST):
            r = loop(ov, D)
            if r is None:
                continue
            areas[c][j] = r[2]
            if D == DFAST:
                loops[c] = (r[0][0], r[0][1], r[1][0], r[1][1])
            print(f"  {c:7s} D={D/60:6.1f}h  area={areas[c][j]:.3f}", flush=True)
    np.savez(CACHE, D_LIST=D_LIST, DFAST=DFAST,
             area_with=areas['with'], area_nomark=areas['nomark'],
             w_ush=loops['with'][0], w_ucd=loops['with'][1], w_dsh=loops['with'][2], w_dcd=loops['with'][3],
             n_ush=loops['nomark'][0], n_ucd=loops['nomark'][1], n_dsh=loops['nomark'][2], n_dcd=loops['nomark'][3])
    print('cached ->', CACHE, flush=True)

z = np.load(CACHE)
D_LIST = z['D_LIST']; Dh = D_LIST / 60.0; DFAST = float(z['DFAST'])
area_with, area_nomark = z['area_with'], z['area_nomark']

fig, ax = plt.subplots(1, 2, figsize=(15, 6))

# (A) the loops at the fast ramp
a = ax[0]
a.plot(z['w_ush'], z['w_ucd'], color='#1f8a4c', lw=2.6, label='WITH mark — ramp-UP')
a.plot(z['w_dsh'], z['w_dcd'], color='#c0392b', lw=2.6, label='WITH mark — ramp-DOWN')
gi = np.interp(GRID, z['w_ush'], z['w_ucd']); gd = np.interp(GRID, z['w_dsh'], z['w_dcd'])
a.fill_between(GRID, gd, gi, where=gi >= gd, color='gold', alpha=0.35, label=f'WITH-mark loop (area {area_with[1]:.2f})')
a.plot(z['n_ush'], z['n_ucd'], color='#555555', lw=2.0, ls='--', label='NO repression — ramp-UP')
a.plot(z['n_dsh'], z['n_dcd'], color='#999999', lw=2.0, ls=':', label='NO repression — ramp-DOWN')
ni = np.interp(GRID, z['n_ush'], z['n_ucd']); nd = np.interp(GRID, z['n_dsh'], z['n_dcd'])
a.fill_between(GRID, nd, ni, where=ni >= nd, color='#888888', alpha=0.25, hatch='//', label=f'NO-repression loop (area {area_nomark[1]:.2f})')
a.set_xlabel('mitogen SHH'); a.set_ylabel('cycle-averaged CyclinD1')
a.set_title(f'(A) CyclinD1–SHH loop at a fast ramp ({DFAST/60:.0f} h)\nthe mark OPENS the loop; without it up & down nearly coincide', fontweight='bold', fontsize=11)
a.legend(fontsize=8.5, loc='upper left'); a.grid(alpha=0.15)

# (B) normalized loop area vs ramp duration, both conditions
a = ax[1]
a.plot(Dh, area_with, '-o', color='#1f8a4c', lw=2.6, ms=6, label='WITH mark (H3K27me3 → CyclinD1)')
a.plot(Dh, area_nomark, '-s', color='#888888', lw=2.6, ms=6, label='NO repression (f0_mk = 1)')
a.fill_between(Dh, area_nomark, area_with, where=area_with >= area_nomark, color='gold', alpha=0.25, label='mark-attributable hysteresis')
a.set_xscale('log'); a.set_xlabel('ramp duration (h)  — faster ← → slower'); a.set_ylabel('normalized loop area')
a.set_title('(B) Loop area vs ramp speed: with vs without the mark\nthe gold gap is the hysteresis the H3K27me3 mark contributes', fontweight='bold', fontsize=11)
a.legend(fontsize=9.5, loc='upper right'); a.grid(alpha=0.15, which='both'); a.set_ylim(bottom=min(0, np.nanmin(area_nomark) * 1.1))

fig.suptitle('Is the CyclinD1–SHH hysteresis caused by the H3K27me3 mark? — WITH vs WITHOUT repression control',
             fontsize=13, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('simulations/sim_ramp_hysteresis_mark_control.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/sim_ramp_hysteresis_mark_control.pdf', bbox_inches='tight')
plt.close()
print('Saved sim_ramp_hysteresis_mark_control.png')
print(f'  area WITH mark  : {np.round(area_with,3)}')
print(f'  area NO repress : {np.round(area_nomark,3)}')
