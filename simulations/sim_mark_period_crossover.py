"""Clean μ-only sweep: the H3K27me3 dilution-vs-restoration RACE as a smooth curve, and the CROSSOVER period.

Controlled version of sim_withdrawal_mark_dilution panel D: hold ALL abundance axes at their model defaults
(k_Cd_translation, kTlEZ, P21_div) and vary ONLY the growth rate mu → cell-cycle period becomes a clean
independent variable, so mark-vs-period shows as a curve instead of a scatter cloud. Constant Hh (SHH=1.0),
single deterministic cell per point (no ensemble/noise).

Two conditions:
  PASSIVE  (f0_mk=1.0): the mark accumulates/dilutes but does NOT repress CyclinD1 → period is set purely
                        by mu, so the mark is a clean readout of the dilution-vs-restoration balance.
  ACTIVE   (f0_mk=0.233): the mark represses CyclinD1 (feedback on) → shows how the feedback bends the race.

Expectation (the regime model): short period → dilution outruns restoration → LOW mark (regime 1);
long period → restoration reaches the (1−Mk) plateau → HIGH mark (regime 2/3). The CROSSOVER = the period
where the mark is halfway up its transition (dilution ≈ restoration). NB restoration (EZH2 methylation) is
itself E2f/cycle-gated (panel C) — that's why the curve isn't a pure 1/period dilution law.

Run:  ./venv/bin/python simulations/sim_mark_period_crossover.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

GNP = dict(MYCN_amplification=1.0, Ptch1_copy_number=1.0, p16=0.0, p18=0.464, kSyP21=0.002, EZH2i=0)
SEL = ['time', 'Dna', 'Mk', 'Cd', 'EZH2']
T_END, N_PTS, SETTLE = 16000, 32000, 8000
MU_GRID = np.geomspace(0.00090, 0.00024, 34)          # fast (short period) -> slow (long period)

_R = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
_R.integrator.setValue('absolute_tolerance', 1e-9); _R.integrator.setValue('relative_tolerance', 1e-6)
try: _R.integrator.setValue('maximum_num_steps', 400000)
except Exception: pass


def run(mu, f0):
    _R.reset(); _R['SHH'] = 1.0
    for k, v in GNP.items(): _R[k] = v
    _R['mu'] = float(mu); _R['f0_mk'] = f0
    try:
        r = _R.simulate(0, T_END, N_PTS, selections=SEL)
    except Exception:
        return np.nan, np.nan, np.nan, np.nan
    t = r['time']; m = t >= SETTLE; tt = t[m]
    Dna = r['Dna'][m]
    dv = np.where((Dna[:-1] > 0.9) & (Dna[1:] < 0.1))[0]
    period = float(np.mean(np.diff(tt[dv])) / 60.0) if len(dv) >= 2 else np.nan
    return period, float(np.mean(r['Mk'][m])), float(np.mean(r['Cd'][m])), float(np.mean(r['EZH2'][m]))


res = {'P': {}, 'A': {}}
for tag, f0 in [('P', 1.0), ('A', 0.233)]:
    per, mk, cd, ez = [], [], [], []
    for mu in MU_GRID:
        p, m, c, e = run(mu, f0)
        per.append(p); mk.append(m); cd.append(c); ez.append(e)
        print(f'  {"PASSIVE" if tag=="P" else "ACTIVE "}  mu={mu:.5f}  period={p if np.isfinite(p) else float("nan"):.1f}h  Mk={m:.3f}  Cd={c:.2f}  EZH2={e:.3f}', flush=True)
    res[tag] = dict(per=np.array(per), mk=np.array(mk), cd=np.array(cd), ez=np.array(ez))


def crossover(per, mk):
    ok = np.isfinite(per) & np.isfinite(mk)
    p, m = per[ok], mk[ok]; o = np.argsort(p); p, m = p[o], m[o]
    if len(p) < 4: return np.nan
    half = 0.5 * (np.nanmin(m) + np.nanmax(m))
    cr = np.where((m[:-1] - half) * (m[1:] - half) < 0)[0]
    if not len(cr): return np.nan
    k = cr[0]; f = (half - m[k]) / (m[k + 1] - m[k]) if m[k + 1] != m[k] else 0.0
    return float(p[k] + f * (p[k + 1] - p[k]))


xoP = crossover(res['P']['per'], res['P']['mk'])
xoA = crossover(res['A']['per'], res['A']['mk'])

fig, ax = plt.subplots(2, 2, figsize=(14, 10))
COND = [('P', '#c47f17', 'PASSIVE (no feedback, f0=1)'), ('A', '#8b1a1a', 'ACTIVE (mark represses CyclinD1)')]

# (A) the race: mark vs period
a = ax[0, 0]
for tag, col, nm in COND:
    o = np.argsort(res[tag]['per']); a.plot(res[tag]['per'][o], res[tag]['mk'][o], '-o', color=col, lw=2.6, ms=5, label=nm)
for xo, col in [(xoP, '#c47f17'), (xoA, '#8b1a1a')]:
    if np.isfinite(xo): a.axvline(xo, color=col, ls='--', lw=1.4, alpha=0.7)
a.text(0.02, 0.96, f'crossover period ≈ {xoP:.0f} h (passive)', transform=a.transAxes, va='top', fontsize=9, color='#7a5210')
a.set_xlabel('cell-cycle period (h)'); a.set_ylabel('steady-state H3K27me3 mark (Mk)')
a.set_title('(A) The dilution-vs-restoration RACE\nshort period → diluted-LOW; long → plateau; dashed = crossover', fontweight='bold', fontsize=11.5)
a.legend(fontsize=9, loc='center right'); a.grid(alpha=0.15)

# (B) CyclinD1 vs period
a = ax[0, 1]
for tag, col, nm in COND:
    o = np.argsort(res[tag]['per']); a.plot(res[tag]['per'][o], res[tag]['cd'][o], '-o', color=col, lw=2.6, ms=5, label=nm)
a.set_xlabel('cell-cycle period (h)'); a.set_ylabel('mean CyclinD1')
a.set_title('(B) CyclinD1 vs period\nACTIVE: the mark represses CyclinD1 more as period lengthens (mark builds)', fontweight='bold', fontsize=11)
a.legend(fontsize=9); a.grid(alpha=0.15)

# (C) EZH2 (writer) vs period — restoration is cycle-gated
a = ax[1, 0]
for tag, col, nm in COND:
    o = np.argsort(res[tag]['per']); a.plot(res[tag]['per'][o], res[tag]['ez'][o], '-o', color=col, lw=2.6, ms=5, label=nm)
a.set_xlabel('cell-cycle period (h)'); a.set_ylabel('mean EZH2 (the writer)')
a.set_title('(C) EZH2 vs period — the RESTORATION capacity\n(E2f/cycle-gated, so it is not constant across period)', fontweight='bold', fontsize=11)
a.legend(fontsize=9); a.grid(alpha=0.15)

# (D) the balance: dilution removal-rate vs restoration, both vs period (PASSIVE)
a = ax[1, 1]
p = res['P']['per']; mk = res['P']['mk']; ez = res['P']['ez']
ok = np.isfinite(p); o = np.argsort(p[ok])
pp = p[ok][o]; mkp = mk[ok][o]; ezp = ez[ok][o]
dil = mkp * 0.5 / pp                                   # avg mark HALVED per division, per hour  (removal)
a.plot(pp, dil, '-o', color='#2b6cb0', lw=2.6, ms=5, label='dilution removal  (Mk/2 ÷ period)')
a.plot(pp, res['P']['mk'][ok][o] * 0, alpha=0)         # keep scale
delmk_h = 0.00070 * 60.0                               # del_mk turnover per hour
a.plot(pp, delmk_h * mkp, '-s', color='#c0392b', lw=2.2, ms=4, label='turnover  (del_mk·Mk)')
if np.isfinite(xoP): a.axvline(xoP, color='#c47f17', ls='--', lw=1.4, alpha=0.7)
a.set_xlabel('cell-cycle period (h)'); a.set_ylabel('mark removal rate (per h)')
a.set_title('(D) Removal terms vs period (PASSIVE): dilution dominates at\nSHORT period, turnover-limited plateau at LONG period', fontweight='bold', fontsize=11)
a.legend(fontsize=9); a.grid(alpha=0.15)

fig.suptitle('Clean μ-only sweep — the H3K27me3 dilution-vs-restoration race and the crossover period (constant Hh=1.0)',
             fontsize=13, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/sim_mark_period_crossover.png', dpi=150, bbox_inches='tight')
plt.savefig('simulations/sim_mark_period_crossover.pdf', bbox_inches='tight')
plt.close()
print(f'\ncrossover period (half-max mark):  passive {xoP:.1f} h   active {xoA:.1f} h')
print('Saved sim_mark_period_crossover.png')
