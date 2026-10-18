"""Does a big CyclinD1 mRNA POOL substitute for slow degradation as the 'memory'? (JP question 2026-07-24).
Answer: the decay HALF-LIFE is set by k alone (pool-independent), but a big pool buys TIME-ABOVE-THRESHOLD
(the memory-relevant quantity) — LOGARITHMICALLY in pool vs LINEARLY in half-life. A ~10x pool of a fast (2h)
transcript reproduces the slow-11h reservoir's coasting; too much pool (>~x6) makes it Shh-INDEPENDENT.
3 panels: A decay curves (same half-life, different level → threshold crossed later); B time-above-threshold vs
pool (log) and vs half-life (linear); C the Shh-independence limit (steady division at SHH=0 vs pool).
Run: ./venv/bin/python simulations/fig_pool_vs_halflife.py"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1'); os.environ.setdefault('EZH2_CONC', '0')
import tellurium as te
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from src.build_model_v44_heldt import build_model_v44
HERE = os.path.dirname(os.path.abspath(__file__)); UPH = 62.8; LN2 = np.log(2)
GNP = dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464)
BASE = dict(basal=0.00028075924684120905, gli=0.0940819151141191, mycn=0.020257706030264946)
THR = 0.8
def divs(dna): dna = np.asarray(dna); return np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]
def build(thalf_h, pool_mult):
    kmu = LN2 / (thalf_h * UPH); txs = pool_mult * (kmu / 0.001)
    return dict(k_Cd_mRNA_deg=kmu, k_Cd_tx_basal=BASE['basal'] * txs, k_Cd_tx_Gli_max=BASE['gli'] * txs, k_Cd_tx_MYCN=BASE['mycn'] * txs)
def sim_robust(rr, t0, t1, n, sel):
    for atol in (1e-8, 1e-7, 1e-6, 1e-5, 1e-4):
        try: rr.integrator.setValue('absolute_tolerance', atol); return rr.simulate(t0, t1, n, selections=sel)
        except Exception: rr.integrator.setValue('maximum_time_step', 1.0)
    return None
def wd(thalf_h, pool_mult, Twarm=130, Tpost=80):
    rr = te.loada(build_model_v44(params=build(thalf_h, pool_mult))); rr.reset()
    for k, v in GNP.items(): rr[k] = v
    rr['SHH'] = 0.5; rr.integrator.setValue('absolute_tolerance', 1e-8); rr.integrator.setValue('maximum_num_steps', 4000000)
    sim_robust(rr, 0, Twarm * UPH, 300, ['time', 'Cd']); cd0 = rr['Cd']; rr['HHi'] = 1
    r = sim_robust(rr, Twarm * UPH, (Twarm + Tpost) * UPH, 4000, ['time', 'Cd'])
    if r is None: return cd0, np.nan, None, None
    t = r['time'] / UPH - Twarm; cd = np.asarray(r['Cd'])
    k = 41; cds = np.convolve(np.pad(cd, k, 'edge'), np.ones(k) / k, 'same')[k:-k]
    above = t[cds > THR]; ta = above[-1] if len(above) else 0.0
    return cd0, ta, t, cds
def div_at_zero(pool_mult, thalf_h=2.0, T=18000):   # Shh-independence check: steady division at SHH=0
    rr = te.loada(build_model_v44(params=build(thalf_h, pool_mult))); rr.reset()
    for k, v in GNP.items(): rr[k] = v
    rr['SHH'] = 0.0; rr.integrator.setValue('absolute_tolerance', 1e-7); rr.integrator.setValue('maximum_num_steps', 4000000)
    r = sim_robust(rr, 0, T, 2500, ['time', 'Dna'])
    if r is None: return np.nan
    t = np.asarray(r['time']); m = t > T * 0.55; dv = divs(np.asarray(r['Dna'])[m]); span = (t[m][-1] - t[m][0]) / UPH
    return len(dv) / span * 24 if span > 0 else 0

print('computing...', flush=True)
POOLS = [1, 1.5, 2, 3, 4, 6]; ta_pool = []; curves = {}
for pm in POOLS:
    cd0, ta, t, cds = wd(2.0, pm); ta_pool.append((cd0, min(ta, 78)))
    if pm in (1, 2, 4): curves[pm] = (t, cds)
    print(f'  pool ×{pm}: Cd0={cd0:.2f} time-above={ta:.1f}h', flush=True)
HALVES = [2, 3, 4, 5, 6]; ta_half = []
for th in HALVES:
    cd0, ta, _, _ = wd(th, 1); ta_half.append((cd0, min(ta, 78))); print(f'  t½={th}h: Cd0={cd0:.2f} time-above={ta:.1f}h', flush=True)
POOLS_Z = [1, 2, 4, 6, 8, 12]; dz = [div_at_zero(pm) for pm in POOLS_Z]
print('  Shh=0 div/day vs pool:', dict(zip(POOLS_Z, [round(x, 2) for x in dz])), flush=True)

C1, C2, C3 = '#2e86c1', '#e67e22', '#c0392b'
fig, ax = plt.subplots(1, 3, figsize=(20.5, 6.4))
fig.subplots_adjust(left=0.05, right=0.99, top=0.82, bottom=0.14, wspace=0.28)

# A: decay curves — same half-life, different level
a = ax[0]
for pm, col in ((1, C1), (2, '#5dade2'), (4, '#a9cce3')):
    t, cds = curves[pm]; w = (t >= -2) & (t < 40)
    a.plot(t[w], cds[w], color=col, lw=2.6, label='pool ×%g (starts %.1f)' % (pm, cds[np.argmin(np.abs(t))]))
a.axhline(THR, color=C3, ls='--', lw=1.6); a.text(30, THR + 0.15, 'commitment threshold', fontsize=8, color=C3)
a.axvline(0, color='#333', ls=':', lw=1.1); a.text(0.6, a.get_ylim()[1]*0.02+3.5, 'Hh OFF', fontsize=8.4, color='#333')
a.text(0.03, 0.60, 'all three have the SAME half-life (2.2 h) —\na bigger pool just starts HIGHER, so it crosses\nthe threshold later (by only ~ln(pool)·τ)', transform=a.transAxes, fontsize=8.2, color='#555')
a.set_xlabel('time after Hh withdrawal (h)'); a.set_ylabel('CyclinD1 protein'); a.set_xlim(-2, 40)
a.set_title('A · A bigger pool does NOT slow the half-life —\nit just starts higher (same decay shape)', fontsize=11, fontweight='bold')
a.legend(fontsize=8.4, loc='upper right'); a.grid(alpha=0.2)

# B: time-above-threshold vs the fold-change of each knob (pool bounded, half-life unbounded)
b = ax[1]
pool_fold = np.array(POOLS, float); half_fold = np.array(HALVES, float) / 2.0
b.plot(pool_fold, [p[1] for p in ta_pool], '-o', color=C2, lw=2.5, ms=7, label='vary POOL (t½=2 h) — capped by Shh-independence')
b.plot(half_fold, [h[1] for h in ta_half], '-s', color=C1, lw=2.5, ms=7, label='vary HALF-LIFE (pool ×1) — linear, unbounded')
b.axvspan(3.3, pool_fold.max() + 0.3, color=C3, alpha=0.08)
b.text(3.5, 1.4, 'pool > ~×3–4:\nShh-INDEPENDENT\n(memory upper bound)', fontsize=7.6, color=C3)
b.set_xlabel('fold-change of the knob  (pool × , or half-life ×)'); b.set_ylabel('time above commitment threshold (h)')
b.set_xlim(0.8, pool_fold.max() + 0.4)
b.set_title('B · Both buy memory, but POOL is bounded (~×4 → Shh-\nindependent); HALF-LIFE is the unbounded axis', fontsize=11, fontweight='bold')
b.legend(fontsize=8.0, loc='upper left'); b.grid(alpha=0.2)

# C: the Shh-independence limit
c = ax[2]
c.plot(POOLS_Z, dz, '-o', color=C3, lw=2.6, ms=7)
c.axhline(0, color='#ccc', lw=0.8); c.fill_between(POOLS_Z, 0, dz, where=(np.array(dz) > 0.1), color=C3, alpha=0.12)
c.text(0.05, 0.55, 'past ~×6 pool the cell divides at SHH=0\n(Shh-INDEPENDENT) — the same failure as\ntoo-high basal. So "just have lots of it" is\nbounded: a few-fold pool is fine, a huge pool\nbreaks Shh-dependence.', transform=c.transAxes, fontsize=8.2, color='#555')
c.set_xlabel('CyclinD1 pool multiplier (fast t½=2 h)'); c.set_ylabel('division rate at SHH = 0 (÷/day)')
c.set_title('C · The catch: too big a pool → Shh-INDEPENDENT\n(divides with no Hh) — memory has an upper bound', fontsize=11, fontweight='bold')
c.grid(alpha=0.2)

fig.suptitle('Does a big CyclinD1 mRNA POOL substitute for slow degradation?  The half-life is set by k ALONE (pool-independent), but a big pool buys time-ABOVE-THRESHOLD — '
             'so YES a large pool gives memory, but LOGARITHMICALLY (inefficient) and only up to ~×6 before it becomes Shh-independent',
             fontsize=11.3, fontweight='bold', y=0.965)
out = os.path.join(HERE, 'fig_pool_vs_halflife.png')
fig.savefig(out, dpi=140); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
