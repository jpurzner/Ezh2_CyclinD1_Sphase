"""Population consequence of the two-compartment transient-G0 amplifier (JP 2026-08-09).

The amplifier does NOT lower the CDKI arrest threshold (it engages only AFTER a cell stops cycling); its
population signature is the DWELL: a cell that transiently gets high CDKI (a stochastic excursion) arrests,
the reservoir re-loads the proximal promoter, and it DWELLS in G0 for many cycle-times before re-entering
when CDKI returns -- instead of the quick (~1 cycle) re-entry without the amplifier. So the transient-G0
subpopulation (JP's Atoh1+ G0 cells) has LONGER dwells and a higher instantaneous G0 fraction.

Protocol per cell: cycle at 1x CDKI, apply a high-CDKI EXCURSION (mult for t_exc), return to 1x, measure
the DWELL = time to first re-entry division. Compares amplifier OFF vs TRANSIENT (a_P0.10/f_P0.80) vs
PERMANENT (a_P0.50/f_P0.30). Population = lognormal excursion magnitudes.

Run:  PYTHONPATH=. ./venv/bin/python simulations/population_amplifier.py   (~10-15 min)
Writes fig_v44_population_amplifier.png/.pdf + prints a summary.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.build_model_v44_heldt import build_model_v44

BP = dict(kPhRbCd=0.1748299602292031, K_CdRb=0.4025007594938155, w_ink4=4.646376316224687,
          kSyP21=0.0004043793202695936, K_cdk6_sink=3.879643001142466, w_p18=1.003950915185932,
          w_p19=0.6655007434943913, k_Cd_tx_Gli_max=0.24197652013556434, k_Cd_tx_MYCN=0.11031706915457935)
CONDS = [('amplifier OFF', None, '#888888'),
         ('transient (C)', dict(a_P=0.10, f_P=0.80), '#b0402f'),
         ('permanent (A)', dict(a_P=0.50, f_P=0.30), '#a02c24')]
T_EXC = 22000       # excursion duration (min) ~ long enough to load P_prox
T_POST = 90000      # post-excursion window to detect re-entry


def make(amp):
    P = dict(BP)
    if amp: P.update(amp)
    rr = te.loada(build_model_v44(with_p27_optionB=True, with_cdk6_gli=True, with_proximal_distal=True, params=P))
    rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
    rr.integrator.setValue('maximum_num_steps', 15000000)
    return rr, P


def setp(rr, P, mult):
    for k, v in P.items():
        try: rr[k] = v
        except Exception: pass
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
    try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
    except Exception: pass
    rr['p18'] = 1.73 * mult; rr['p19'] = 0.58 * mult; rr['kSyP21'] = BP['kSyP21'] * mult


def dwell(rr, P, mult):
    """excursion to `mult` CDKI, then back to 1x; return dwell-to-re-entry (h) or None (permanent)."""
    setp(rr, P, mult); rr.reset()
    rr.simulate(0, T_EXC, T_EXC // 10, selections=['time', 'Dna'])
    setp(rr, P, 1.0)  # return to baseline CDKI (species state carries the loaded reservoir/P_prox)
    r2 = rr.simulate(T_EXC, T_EXC + T_POST, T_POST // 10, selections=['time', 'Dna'])
    dna = r2['Dna']; dd = np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]
    return float((r2['time'][dd[0]] - T_EXC) / 60) if len(dd) else None


# ---- (1) dwell vs excursion magnitude ----
mags = [2.0, 3.0, 4.0, 5.0, 6.0]
dwell_curve = {}
for name, amp, _ in CONDS:
    rr, P = make(amp)
    dwell_curve[name] = [dwell(rr, P, m) for m in mags]
    print(f'{name:14s} dwell(h) vs excursion {mags}: {dwell_curve[name]}')

# ---- (2) population of excursion magnitudes (lognormal, those that reach arrest) ----
rng = np.random.default_rng(0); N = 30
cv = 0.5; sd = np.sqrt(np.log(1 + cv**2)); muL = np.log(4.0) - 0.5 * sd**2   # centered ~4x (arrest-competent tail)
pop_mags = np.clip(rng.lognormal(muL, sd, N), 2.0, 9.0)
pop = {}
for name, amp, _ in CONDS:
    rr, P = make(amp)
    ds = [dwell(rr, P, float(m)) for m in pop_mags]
    fin = [d for d in ds if d is not None]
    perm = 100.0 * (len(ds) - len(fin)) / len(ds)
    pop[name] = dict(dwells=ds, median=(float(np.median(fin)) if fin else None), pct_permanent=perm)
    print(f'{name:14s} population: median dwell={pop[name]["median"]}h  permanent={perm:.0f}%')

# ---- figure ----
fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 5.2))
for name, amp, col in CONDS:
    y = [d if d is not None else np.nan for d in dwell_curve[name]]
    axA.plot(mags, y, 'o-', color=col, lw=1.7, ms=7, label=name)
    for xm, yv, dv in zip(mags, y, dwell_curve[name]):
        if dv is None: axA.annotate('lock', (xm, 0), color=col, fontsize=7, ha='center', va='bottom')
axA.set_xlabel('excursion CDKI (× baseline)'); axA.set_ylabel('transient-G0 dwell after re-entry (h)')
axA.set_title('(A) Dwell vs CDKI excursion — amplifier prolongs G0', fontweight='bold')
axA.axhline(24, color='k', ls=':', lw=0.8); axA.text(2.05, 26, 'one cycle (~24 h)', fontsize=7.5)
axA.legend(fontsize=8.5); axA.set_ylim(bottom=0)

off = [d for d in pop['amplifier OFF']['dwells'] if d is not None]
tr = [d for d in pop['transient (C)']['dwells'] if d is not None]
bins = np.linspace(0, 260, 14)
axB.hist(off, bins=bins, alpha=0.6, color='#888888', label=f"OFF (median {pop['amplifier OFF']['median']:.0f}h)")
axB.hist(tr, bins=bins, alpha=0.6, color='#b0402f', label=f"transient (median {pop['transient (C)']['median']:.0f}h)")
axB.set_xlabel('transient-G0 dwell (h)'); axB.set_ylabel('cells')
axB.set_title('(B) Population dwell distribution (excursion CV 0.5, N=30)', fontweight='bold')
axB.legend(fontsize=8.5)
fig.suptitle('Population effect of the reservoir-sustained transient-G0 amplifier (with_proximal_distal)',
             fontsize=12.5, fontweight='bold')
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_v44_population_amplifier')
fig.savefig(out + '.png', dpi=140, bbox_inches='tight'); fig.savefig(out + '.pdf', bbox_inches='tight')
print('wrote', out + '.png/.pdf')
