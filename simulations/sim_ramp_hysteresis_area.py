"""TIMING sweep, AREA metric: measure the CyclinD1-vs-SHH hysteresis LOOP AREA (up-ramp minus down-ramp)
instead of a threshold crossing. Continuous & smooth — no dependence on which discrete cycle divides last,
so it eliminates the exit-threshold jaggedness of sim_ramp_timing_exploration.

The mark lags the mitogen: on the UP ramp the mark is lower-than-steady (less repression → CyclinD1 HIGH);
on the DOWN ramp it's higher-than-steady (CyclinD1 LOW). The area between the two cycle-averaged
CyclinD1-vs-SHH curves is the hysteresis, normalized by the cycling CyclinD1 (dimensionless, scale-free).

Result: the loop area peaks at an INTERMEDIATE kinetic timescale, NOT the slowest. Too-fast marks track
the ramp (no lag → no loop); too-slow marks over-accumulate in the arrested low-SHH state (no replicative
dilution there to clear them), so the up-ramp starts already repressed and can't re-enter → loop collapses.
A mild fast>slow-ramp tilt rides on top. (Pre-equilibration is scaled to the kinetic timescale — see
_preeq_for — so a fixed preeq can't under-load slow marks and fake this structure.)

Panels:
  (A) area over (ramp duration × mark turnover t½)   — del_mk
  (B) area over (ramp duration × EZH2 t½)            — kDeEZ
  (C) area over (ramp duration × cell-cycle period)  — mu (dilution frequency)
  (D) the CyclinD1-vs-SHH loop itself (up vs down) at a fast ramp, default timing — what the area measures

Run:  ./venv/bin/python simulations/sim_ramp_hysteresis_area.py [--fresh]
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
CACHE = 'simulations/sim_ramp_hysteresis_area_cache.npz'
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
    dt = np.median(np.diff(t)) if len(t) > 1 else 1.0
    w = max(3, min(int(win / dt), max(3, len(x) // 2)))
    return np.convolve(x, np.ones(w) / w, mode='same')


def _preeq_for(ov):
    """pre-equilibrate for 4 turnover-times of the SLOWEST kinetic species (mark or EZH2), so the
    down-ramp starts from a FULLY-loaded (repressed) mark. A fixed PREEQ under-loads slow marks and
    artificially collapses the loop at the slow-kinetics corner (finite-preequilibration confound)."""
    dm = ov.get('del_mk', DEFAULTS['del_mk']); ke = ov.get('kDeEZ', DEFAULTS['kDeEZ'])
    tau = max(LN2 / dm, LN2 / ke)                               # slowest half-life (min)
    return float(min(20000.0, max(PREEQ, 4.0 * tau)))          # 4 t½, capped at 20000 min (~333 h)


def ramp_curve(ov, D, entry):
    """cycle-averaged CyclinD1 vs SHH along the ramp, sorted by SHH increasing."""
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


def hysteresis_area(ov, D):
    up = ramp_curve(ov, D, True); dn = ramp_curve(ov, D, False)
    if up is None or dn is None:
        return np.nan
    cu = np.interp(GRID, up[0], up[1]); cd = np.interp(GRID, dn[0], dn[1])
    ref = 0.5 * (np.nanmax(up[1]) + np.nanmax(dn[1]))
    return float(np.trapz(cu - cd, GRID) / (ref + 1e-9))          # dimensionless loop area


NPD, NPV = 14, 12
D_RANGE = np.round(np.geomspace(2500.0, 24000.0, NPD), 0)
VARS = {
    'del_mk': np.geomspace(0.0021, 0.00023, NPV),
    'kDeEZ':  np.geomspace(0.00045, 0.00005, NPV),
    'mu':     np.geomspace(0.00095, 0.00028, NPV),
}

if FRESH or not os.path.exists(CACHE):
    A = {}
    for param, vals in VARS.items():
        grid = np.full((NPV, NPD), np.nan)
        for i, v in enumerate(vals):
            for j, D in enumerate(D_RANGE):
                grid[i, j] = hysteresis_area({param: v}, D)
            print(f"  {param:8s} {i+1}/{NPV} v={v:.5f}  area={np.round(grid[i],3)}", flush=True)
        A[param] = grid
    # representative loop for panel D: fast ramp, default timing
    Dfast = D_RANGE[2]
    up = ramp_curve({}, Dfast, True); dn = ramp_curve({}, Dfast, False)
    np.savez(CACHE, D_RANGE=D_RANGE, delmk_vals=VARS['del_mk'], kdeez_vals=VARS['kDeEZ'], mu_vals=VARS['mu'],
             A_delmk=A['del_mk'], A_kdeez=A['kDeEZ'], A_mu=A['mu'],
             up_sh=up[0], up_cd=up[1], dn_sh=dn[0], dn_cd=dn[1], Dfast=Dfast)
    print('cached ->', CACHE, flush=True)

z = np.load(CACHE)
D_RANGE = z['D_RANGE']; Dh = D_RANGE / 60.0
A = {'del_mk': z['A_delmk'], 'kDeEZ': z['A_kdeez'], 'mu': z['A_mu']}
delmk_vals, kdeez_vals, mu_vals = z['delmk_vals'], z['kdeez_vals'], z['mu_vals']


def _edges(c):
    lc = np.log10(np.asarray(c, float)); e = np.empty(len(c) + 1)
    e[1:-1] = 10 ** (0.5 * (lc[:-1] + lc[1:])); e[0] = 10 ** (lc[0] - 0.5 * (lc[1] - lc[0])); e[-1] = 10 ** (lc[-1] + 0.5 * (lc[-1] - lc[-2])); return e


XE = _edges(Dh)
vmax = max(np.nanmax(A['del_mk']), np.nanmax(A['kDeEZ']), np.nanmax(A['mu']))
fig, ax = plt.subplots(2, 2, figsize=(15, 11))
for a, key, yvals, ylab, ttl in [
    (ax[0, 0], 'del_mk', hl(delmk_vals), 'mark turnover t½ (h)', '(A) del_mk (mark turnover)'),
    (ax[0, 1], 'kDeEZ', hl(kdeez_vals), 'EZH2 t½ (h)', '(B) kDeEZ (EZH2 half-life)'),
    (ax[1, 0], 'mu', LN2 / mu_vals / 60.0, 'cell-cycle timescale ln2/mu (h)', '(C) mu (dilution frequency / cycle period)')]:
    YE = _edges(yvals)
    pm = a.pcolormesh(XE, YE, A[key], cmap='magma', shading='flat', vmin=0, vmax=vmax)
    a.set_xscale('log'); a.set_yscale('log')
    a.set_xlabel('ramp duration (h)  — faster ← → slower'); a.set_ylabel(ylab)
    a.set_title(ttl + '\nloop AREA peaks at INTERMEDIATE timescale; mild fast>slow-ramp tilt', fontweight='bold', fontsize=11)
    fig.colorbar(pm, ax=a, label='normalized loop area')

# (D) the CyclinD1-vs-SHH loop (up vs down) at a fast ramp, default timing
a = ax[1, 1]
a.plot(z['up_sh'], z['up_cd'], color='#1f8a4c', lw=2.4, label='ramp-UP (mark lags low → CyclinD1 high)')
a.plot(z['dn_sh'], z['dn_cd'], color='#c0392b', lw=2.4, label='ramp-DOWN (mark lags high → CyclinD1 low)')
gi = np.interp(GRID, z['up_sh'], z['up_cd']); gd = np.interp(GRID, z['dn_sh'], z['dn_cd'])
a.fill_between(GRID, gd, gi, where=gi >= gd, color='gold', alpha=0.3, label='hysteresis area')
a.set_xlabel('mitogen SHH'); a.set_ylabel('cycle-averaged CyclinD1')
a.set_title(f'(D) The CyclinD1–SHH hysteresis loop (fast ramp {z["Dfast"]/60:.0f} h, default timing)\nthe area between up- & down-ramp is the smooth metric', fontweight='bold', fontsize=11)
a.legend(fontsize=8.5, loc='upper left'); a.grid(alpha=0.15)

fig.suptitle('CyclinD1–SHH hysteresis AREA (smooth metric) — H3K27me3 kinetics & dilution frequency vs mitogen-ramp speed',
             fontsize=12.5, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/sim_ramp_hysteresis_area.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/sim_ramp_hysteresis_area.pdf', bbox_inches='tight')
plt.close()
print('Saved sim_ramp_hysteresis_area.png')
