"""Phase diagram: does the H3K27me3 mark ACCUMULATE or RELAX in a persistently arrested G0 cell?

The outcome is a race between two timescales (Purzner, this thread):
  - EZH2-writer PERSISTENCE  (protein half-life, set by kDeEZ) — a slower-decaying writer keeps writing
    in G0 AND settles to a HIGHER basal level (EZH2_basal = kTlEZ*kEZbas/(kDeEZm*kDeEZ)), pushing the
    arrested mark UP.
  - MARK TURNOVER            (demethylation half-life, set by del_mk) — faster turnover pulls the mark
    DOWN toward whatever the (falling) writer can sustain.
Because arrest removes replicative dilution, the arrested steady mark is set purely by writer-vs-turnover;
the cycling mark is additionally diluted. So the sign of (arrested - cycling) can go EITHER way.

For each (EZH2 half-life x mark half-life) we simulate a CYCLING GNP (SHH=1.0; mean Mk over the sawtooth)
and an ARRESTED cell (SHH=0.1; steady Mk after ~2 weeks, no S phase -> no dilution), and color the
difference. The measured EZH2 half-life (G2/G0=1.48x IF -> kDeEZ=1.5e-4 ~77 h; range ~46-580 h) and the
mark half-life (de-repression t1/2~15 h -> del_mk) are drawn as the operating point + uncertainty box.

Run:  ./venv/bin/python simulations/fig_v44_g0_mark_regime.py [--mini] [--fresh]
Cached in simulations/fig_v44_g0_mark_regime_cache.npz (~20 min cold; re-plots instantly).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

MINI = '--mini' in sys.argv
FRESH = '--fresh' in sys.argv or MINI
CACHE = 'simulations/fig_v44_g0_mark_regime_cache.npz' if not MINI else 'simulations/fig_v44_g0_mark_regime_mini.npz'
LN2 = np.log(2.0)
KDEEZ0, DELMK0 = 0.00015, 0.00070                    # operating point (baked)
hl = lambda rate: LN2 / rate / 60.0                  # 1/min -> half-life in hours

M = build_model_v44(with_ezh2=True, with_hh=True)
GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, HHi=0, EZH2i=0)


def rr():
    r = te.loada(M); r.integrator.setValue('relative_tolerance', 1e-6)
    try: r.integrator.setValue('maximum_num_steps', 200000)
    except Exception: pass
    return r


def _run(r, shh, kdeez, delmk, T, settle):
    r.reset()
    for k, v in GNP.items(): r[k] = v
    r['SHH'] = shh; r['kDeEZ'] = float(kdeez); r['del_mk'] = float(delmk)
    for atol in (1e-8, 1e-7, 1e-6):
        try:
            r.integrator.setValue('absolute_tolerance', atol)
            d = r.simulate(0, T, int(T / 2), selections=['time', 'Mk', 'MPF', 'EZH2'])
            m = d['time'] >= settle
            pk, _ = find_peaks(d['MPF'][m], prominence=0.15, distance=200)
            return float(np.nanmean(d['Mk'][m])), len(pk), float(np.nanmean(d['EZH2'][m]))
        except Exception:
            continue
    return np.nan, -1, np.nan


if FRESH or not os.path.exists(CACHE):
    if MINI:
        KDEEZ = np.geomspace(6e-5, 5e-4, 6); DELMK = np.geomspace(3e-4, 2.5e-3, 6)
    else:
        KDEEZ = np.geomspace(5e-5, 6e-4, 16)         # EZH2 half-life ~19-231 h
        DELMK = np.geomspace(2e-4, 3e-3, 16)         # mark half-life ~3.9-58 h
    r = rr()
    CYC = np.full((len(DELMK), len(KDEEZ)), np.nan)   # cycling mean Mk
    ARR = np.full((len(DELMK), len(KDEEZ)), np.nan)   # arrested steady Mk
    NDIV_C = np.full((len(DELMK), len(KDEEZ)), np.nan)  # cycling divisions (flag mark-arrest)
    import time as _t
    for i, dm in enumerate(DELMK):
        for j, ke in enumerate(KDEEZ):
            CYC[i, j], NDIV_C[i, j], _ = _run(r, 1.0, ke, dm, 15000, 6000)
            ARR[i, j], _, _ = _run(r, 0.10, ke, dm, 26000, 20000)   # long: reach G0 steady state
        np.savez(CACHE, KDEEZ=KDEEZ, DELMK=DELMK, CYC=CYC, ARR=ARR, NDIV_C=NDIV_C)
        print(f'  row {i+1}/{len(DELMK)}  del_mk={dm:.1e} (t1/2 {hl(dm):.0f}h)', flush=True)
    print('cached ->', CACHE, flush=True)

z = np.load(CACHE)
KDEEZ, DELMK, CYC, ARR, NDIV_C = z['KDEEZ'], z['DELMK'], z['CYC'], z['ARR'], z['NDIV_C']
DIFF = ARR - CYC                                      # >0 accumulate, <0 relax
HLE = hl(KDEEZ); HLM = hl(DELMK)                      # half-life axes (hours)


def _edges(c):
    lc = np.log10(c); e = np.empty(len(c) + 1)
    e[1:-1] = 10 ** (0.5 * (lc[:-1] + lc[1:])); e[0] = 10 ** (lc[0] - 0.5 * (lc[1] - lc[0]))
    e[-1] = 10 ** (lc[-1] + 0.5 * (lc[-1] - lc[-2])); return e


XE, YE = _edges(HLE), _edges(HLM)
fig, ax = plt.subplots(figsize=(9.2, 7.2))
vlim = np.nanmax(np.abs(DIFF))
pm = ax.pcolormesh(XE, YE, DIFF, cmap='RdBu_r', norm=TwoSlopeNorm(vcenter=0, vmin=-vlim, vmax=vlim), shading='flat')
# zero boundary (accumulate vs relax)
try:
    cs = ax.contour(HLE, HLM, DIFF, levels=[0], colors='k', linewidths=2.2)
    ax.clabel(cs, fmt={0: 'arrested = cycling'}, fontsize=9)
except Exception:
    pass
# mark-arrested region (SHH=1.0 fails to cycle -> strong-brake regime): hatch
mask = NDIV_C < 2
if mask.any():
    ax.pcolor(XE, YE, np.ma.masked_where(~mask, mask), hatch='xx', alpha=0, edgecolor='k', linewidth=0)

ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('EZH2 protein half-life (h)  →  more persistent writer', fontsize=11)
ax.set_ylabel('H3K27me3 turnover half-life (h)  →  slower erasing', fontsize=11)

# operating point + measured uncertainty box
ax.plot(hl(KDEEZ0), hl(DELMK0), 'o', ms=13, mfc='yellow', mec='k', mew=1.8, zorder=6, label='operating point')
ax.add_patch(plt.Rectangle((46, hl(DELMK0 * 2.0)), 580 - 46, hl(DELMK0 * 0.5) - hl(DELMK0 * 2.0),
                           fill=False, ec='k', ls='--', lw=1.5, zorder=5))
ax.text(hl(KDEEZ0) * 1.06, hl(DELMK0) * 1.06, ' measured\n EZH2 t½ 77 h\n (46–580 h)', fontsize=8.5, zorder=6)

# regime labels
ax.text(0.06, 0.93, 'RELAX\n(arrested mark < cycling —\nwriter collapse wins,\nG0 reversible)', transform=ax.transAxes,
        fontsize=10, va='top', color='#1a4a7a', fontweight='bold')
ax.text(0.60, 0.14, 'ACCUMULATE\n(arrested mark > cycling —\ndilution loss wins,\nG0 ratchets)', transform=ax.transAxes,
        fontsize=10, va='top', color='#8a2020', fontweight='bold')

cb = fig.colorbar(pm, ax=ax, label='arrested-G0 mark  −  cycling mean mark  (Mk units)')
ax.legend(loc='lower left', fontsize=9)
ax.set_title('Does H3K27me3 accumulate or relax in arrested G0?\nA race between EZH2-writer persistence and mark turnover',
             fontsize=12.5, fontweight='bold')
plt.tight_layout()
plt.savefig('simulations/fig_v44_g0_mark_regime.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/fig_v44_g0_mark_regime.pdf', bbox_inches='tight')
plt.close()
print('Saved fig_v44_g0_mark_regime.png')
op_i = int(np.argmin(np.abs(DELMK - DELMK0))); op_j = int(np.argmin(np.abs(KDEEZ - KDEEZ0)))
print(f'at operating point: arrested {ARR[op_i,op_j]:.3f}  cycling {CYC[op_i,op_j]:.3f}  diff {DIFF[op_i,op_j]:+.3f}')
