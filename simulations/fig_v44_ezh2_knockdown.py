"""How does the cell cycle respond to DECREASING EZH2 EXPRESSION (protein level), as opposed to
catalytic inhibition? Models the developmental EZH2 decline / shRNA knockdown by scaling EZH2
translation (kTlEZ). Question: does lowering EZH2 give a 'more rapid' cell cycle?

We separate the senses of 'more rapid':
  - cell-cycle PERIOD (length of one cycle)        -> growth-limited in this model
  - DIVISION RATE / fraction cycling (net output)  -> set by CyclinD1 vs threshold
  - transient G0 (time spent quiescent)            -> the G0/G1 commitment delay

Part A: steady-state knockdown sweep (EZH2 expression 100% -> 0%) in GNP and MB.
Part B: a dynamic step-down (WT, then drop EZH2 expression mid-run) to show the live response.

Run:  ./venv/bin/python simulations/fig_v44_ezh2_knockdown.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44
from simulations.validate_v44 import count_divisions, mean_settled, classify, P27_THR, SEL

KTLEZ0 = 0.004                       # WT EZH2 translation rate
KD = [1.0, 0.75, 0.50, 0.30, 0.15, 0.05, 0.0]   # EZH2-expression knockdown factor (1=WT, 0=null)
T_END, N_PTS = 10080, 20160
_M = build_model_v44()


def run(kd, shh=0.5, mycn=1.0, ptch1=1.0):
    for atol in (1e-9, 1e-8, 1e-7):
        rr = te.loada(_M); rr.integrator.setValue("absolute_tolerance", atol)
        rr.integrator.setValue("relative_tolerance", 1e-6)
        rr['SHH'] = shh; rr['MYCN_amplification'] = mycn; rr['Ptch1_copy_number'] = ptch1
        rr['kTlEZ'] = KTLEZ0 * kd
        try:
            return rr.simulate(0, T_END, N_PTS, selections=SEL)
        except Exception:
            continue
    return None


ctxs = {'GNP': dict(shh=0.5, mycn=1.0, ptch1=1.0), 'MB': dict(shh=0.5, mycn=2.8, ptch1=0.1)}
D = {c: dict(ez=[], cd=[], ndiv=[], per=[], g0=[]) for c in ctxs}
print("=" * 74)
print("DECREASING EZH2 EXPRESSION (translation knockdown) -> cell cycle")
for c, cond in ctxs.items():
    print(f"\n{c}:   EZH2%   EZH2prot  CyclinD1   div/168h   period   G0%")
    for kd in KD:
        r = run(kd, **cond)
        n, _, per = count_divisions(r); f, _, _ = classify(r, P27_THR)   # f = COUNT-fraction (flow quiescent G0)
        ez = mean_settled(r, 'EZH2'); cd = mean_settled(r, 'Cd_mRNA')
        pp = np.mean(per) if len(per) else np.nan
        D[c]['ez'].append(ez); D[c]['cd'].append(cd); D[c]['ndiv'].append(n)
        D[c]['per'].append(pp); D[c]['g0'].append(f['G0'])
        print(f"        {kd*100:5.0f}%  {ez:8.3f}  {cd:8.3f}   {n:8d}   {pp:6.1f}h  {f['G0']:5.1f}")

# ---- Part B: dynamic step-down (GNP): WT then knock EZH2 expression down at t0 ----
t0 = 6000
rr = te.loada(_M); rr.integrator.setValue("absolute_tolerance", 1e-9); rr.integrator.setValue("relative_tolerance", 1e-6)
rr['SHH'] = 0.5
r1 = rr.simulate(0, t0, 12000, selections=SEL)
rr['kTlEZ'] = KTLEZ0 * 0.1            # drop EZH2 expression to 10%
r2 = rr.simulate(t0, 18000, 24000, selections=SEL)
tt = np.concatenate([r1['time'], r2['time']]) / 60.0
EZ = np.concatenate([r1['EZH2'], r2['EZH2']])
CD = np.concatenate([r1['Cd_mRNA'], r2['Cd_mRNA']])
MPF = np.concatenate([r1['MPF'], r2['MPF']])
P21 = np.concatenate([r1['P21'], r2['P21']])

# ================= figure =================
fig = plt.figure(figsize=(14, 9))
gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.38)
xpct = [k * 100 for k in KD]

ax = fig.add_subplot(gs[0, 0])
for c, col in [('GNP', '#1b9e77'), ('MB', '#762A83')]:
    ax.plot(xpct, D[c]['cd'], col, lw=2.2, marker='o', label=c)
ax.set_xlabel('EZH2 expression (% of WT)'); ax.set_ylabel('CyclinD1 transcript'); ax.invert_xaxis()
ax.set_title('(a) Less EZH2 -> more CyclinD1', fontweight='bold', fontsize=11); ax.legend(fontsize=8); ax.grid(alpha=0.2)

ax = fig.add_subplot(gs[0, 1])
for c, col in [('GNP', '#1b9e77'), ('MB', '#762A83')]:
    ax.plot(xpct, D[c]['ndiv'], col, lw=2.2, marker='o', label=c)
ax.set_xlabel('EZH2 expression (% of WT)'); ax.set_ylabel('divisions / 168 h'); ax.invert_xaxis()
ax.set_title('(b) Division RATE (net output) rises', fontweight='bold', fontsize=11); ax.legend(fontsize=8); ax.grid(alpha=0.2)

ax = fig.add_subplot(gs[0, 2])
for c, col in [('GNP', '#1b9e77'), ('MB', '#762A83')]:
    ax.plot(xpct, D[c]['per'], col, lw=2.2, marker='o', label=c)
ax.set_xlabel('EZH2 expression (% of WT)'); ax.set_ylabel('cell-cycle period (h)'); ax.invert_xaxis()
ax.set_ylim(0, 30)
ax.set_title('(c) PERIOD is ~flat (growth-limited)\n-> not a faster individual cycle', fontweight='bold', fontsize=10.5)
ax.legend(fontsize=8); ax.grid(alpha=0.2)

ax = fig.add_subplot(gs[1, 0])
for c, col in [('GNP', '#1b9e77'), ('MB', '#762A83')]:
    ax.plot(xpct, D[c]['g0'], col, lw=2.2, marker='s', label=c)
ax.set_xlabel('EZH2 expression (% of WT)'); ax.set_ylabel('transient G0 (p27) %'); ax.invert_xaxis()
ax.set_title('(d) Transient G0 shrinks (faster G0->G1)', fontweight='bold', fontsize=11); ax.legend(fontsize=8); ax.grid(alpha=0.2)

# dynamic step-down (spans bottom-middle + bottom-right)
ax = fig.add_subplot(gs[1, 1:])
ax.axvline(t0 / 60.0, color='k', ls='--', alpha=0.6)
ax.text(t0/60.0 + 2, ax.get_ylim()[1] if False else 4.6, 'EZH2 expression\n dropped to 10%', fontsize=8, color='#b03a2e')
ax.plot(tt, EZ, color='#8e44ad', lw=1.6, label='EZH2 protein')
ax.plot(tt, CD, color='#117a65', lw=1.6, label='CyclinD1 transcript')
ax.plot(tt, MPF * 8, color='#e67e22', lw=0.9, alpha=0.8, label='divisions (MPF x8)')
ax.set_xlabel('time (h)'); ax.set_ylabel('level (a.u.)')
ax.set_title('(e) Dynamic step-down (GNP): drop EZH2 mid-run -> CyclinD1 rises, divisions speed up',
             fontweight='bold', fontsize=10.5)
ax.legend(fontsize=8, loc='upper left'); ax.grid(alpha=0.2)

# report period before/after in the step-down
def per_window(t, mpf, a, b):
    m = (t >= a) & (t < b); tt2 = t[m]; dt = tt2[1]-tt2[0]
    pk, _ = find_peaks(mpf[m], prominence=0.1, distance=int(150/dt))
    iv = np.diff(tt2[pk]) if len(pk) > 1 else []
    return (np.mean(iv), len(pk))
pa = per_window(tt, MPF, 60, t0/60.0); pb = per_window(tt, MPF, t0/60.0 + 20, 300)
print(f"\nStep-down (GNP): before drop period~{pa[0]:.1f}h ({pa[1]} div); after drop period~{pb[0]:.1f}h ({pb[1]} div)")

fig.suptitle('Decreasing EZH2 EXPRESSION: more divisions & shorter G0, but ~same cycle period',
             fontsize=13, fontweight='bold', y=0.99)
plt.savefig('simulations/fig_v44_ezh2_knockdown.png', dpi=170, bbox_inches='tight')
plt.savefig('simulations/fig_v44_ezh2_knockdown.pdf', bbox_inches='tight')
plt.close()
print("\nSaved: fig_v44_ezh2_knockdown.png")
