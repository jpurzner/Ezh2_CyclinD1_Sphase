"""Within-cycle anatomy: where is the transient G0, how long is it, and how does CyclinD1 set it?

Plots ~2.5 cycles of the key species (CyclinD1 protein, phospho-Rb, p27[=P21], cell mass, Skp2,
replication aRc/Dna, MPF) and measures, per cycle:
  - transient G0 by the p27 marker:        time from division until P21 < P21_lo  (p27 cleared)
  - transient G0 by the phospho-Rb marker: time from division until pRb > 0.5*max  (Rb committed)
  - G1 (committed, pre-S) and the growth time to the commitment size gate
across three CyclinD1 levels (low / default / high, set via k_Cd_translation).

Run:  ./venv/bin/python simulations/fig_v44_cycle_anatomy.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
import tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

SEL = ["time", "Cd", "pRb", "P21", "Skp2", "aRc", "Dna", "MPF", "mass", "E2f"]
T_END, N_PTS, SETTLE = 14000, 70000, 6000
P21_LO = 0.1


def simulate(kTl=0.26, shh=0.5):
    rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True,
                                  params={"k_Cd_translation": kTl}))
    rr.integrator.setValue("absolute_tolerance", 1e-9)
    rr.integrator.setValue("relative_tolerance", 1e-6)
    rr['SHH'] = shh
    return rr.simulate(0, T_END, N_PTS, selections=SEL)


def division_times(res):
    """Divisions = sharp drops in mass (halving event)."""
    t = res['time']; mass = res['mass']
    drops = np.where((mass[:-1] - mass[1:]) > 0.2 * mass[:-1])[0]
    # collapse consecutive indices
    if len(drops) == 0:
        return np.array([])
    keep = [drops[0]]
    for d in drops[1:]:
        if d - keep[-1] > 50:
            keep.append(d)
    return t[keep]


def g0_metrics(res):
    """Alignment-free transient-G0 by both markers = fraction(pre-S & marker-G0) x period.
    With the two-step Rb, pRb(hyper) is low in G0, so G0(p27) and G0(pRb-low) should agree."""
    t = res['time']; m = t >= SETTLE; tt = t[m]; dt = tt[1] - tt[0]
    P21 = res['P21'][m]; pRb = res['pRb'][m]; aRc = res['aRc'][m]; Dna = res['Dna'][m]
    pk, _ = find_peaks(res['MPF'][m], prominence=0.1, distance=int(150 / dt))
    per = float(np.mean(np.diff(tt[pk]))) / 60 if len(pk) > 1 else np.nan
    in_S = (aRc > 0.05) & (Dna < 0.98); in_G2 = Dna >= 0.98; preS = ~in_S & ~in_G2
    tot = (preS | in_S | in_G2).mean() or 1
    prb_thr = 0.5 * pRb.max()
    g0_p27 = per * (preS & (P21 > P21_LO)).mean() / tot
    g0_prb = per * (preS & (pRb < prb_thr)).mean() / tot
    return per, g0_p27, g0_prb


LEVELS = [("low (kTl=0.30)", 0.30), ("default (kTl=0.26)", 0.26), ("high (kTl=0.40)", 0.40)]
sims = {lab: simulate(kTl=k) for lab, k in LEVELS}

print("=" * 78)
print("Transient-G0 duration vs CyclinD1 protein level")
print(f"{'level':>20} {'meanCd':>7} {'period':>7} {'G0(p27)':>8} {'G0(pRb)':>8}")
for lab, k in LEVELS:
    r = sims[lab]; m = r['time'] >= SETTLE
    per, g0p27, g0prb = g0_metrics(r)
    print(f"{lab:>20} {np.mean(r['Cd'][m]):>7.3f} {per:>6.1f}h {g0p27:>7.2f}h {g0prb:>7.2f}h")

# ---- figure: 2.5-cycle anatomy of the DEFAULT level ----
r = sims["default (kTl=0.26)"]
t = r['time'] / 60.0
divs = division_times(r) / 60.0
divs = divs[(divs >= SETTLE/60.0)]
t0 = divs[0]; t1 = divs[min(3, len(divs)-1)]
m = (t >= t0 - 1) & (t <= t1 + 1)
tt = t[m] - t0

fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
prb_thr = 0.5 * r['pRb'][r['time'] >= SETTLE].max()

ax = axes[0]
ax.plot(tt, r['Cd'][m], color='#117a65', lw=1.8, label='CyclinD1 protein (Cd)')
ax.plot(tt, r['mass'][m], color='#607d8b', lw=1.4, ls='--', label='cell mass')
ax.axhline(1.3, color='#607d8b', ls=':', alpha=0.7, label='M_commit (size gate)')
ax.set_ylabel('a.u.'); ax.legend(fontsize=8, loc='upper right'); ax.grid(alpha=0.2)
ax.set_title('Within-cycle anatomy (GNP, default CyclinD1): the transient G0 is the p27-high / phospho-Rb-low window after division',
             fontsize=10, fontweight='bold')

ax = axes[1]
ax.plot(tt, r['pRb'][m], color='#c0392b', lw=1.8, label='phospho-Rb (pRb)')
ax.plot(tt, r['P21'][m], color='#b9770e', lw=1.8, label='p27 (=P21)')
ax.plot(tt, r['Skp2'][m], color='#2980b9', lw=1.4, label='Skp2')
ax.axhline(prb_thr, color='#c0392b', ls=':', alpha=0.6, label='pRb commit thr (0.5 max)')
ax.axhline(P21_LO, color='#b9770e', ls=':', alpha=0.6, label='p27-low thr')
ax.set_ylabel('a.u.'); ax.legend(fontsize=8, loc='upper right'); ax.grid(alpha=0.2)

ax = axes[2]
ax.plot(tt, r['aRc'][m], color='#17a2b8', lw=1.6, label='aRc (active forks, S)')
ax.plot(tt, r['Dna'][m], color='#8e44ad', lw=1.6, label='Dna (0->1)')
ax.plot(tt, r['MPF'][m], color='#e67e22', lw=1.6, label='MPF (mitosis)')
ax.set_ylabel('a.u.'); ax.set_xlabel('time since a division (h)')
ax.legend(fontsize=8, loc='upper right'); ax.grid(alpha=0.2)

# shade the transient G0 (p27-high) of the FIRST displayed cycle on all axes
p27seg = r['P21'][m]
cyc1_len = (divs[1] - t0) if len(divs) > 1 else (t1 - t0)
g0_region = (tt >= 0) & (tt < cyc1_len) & (p27seg > P21_LO)
g0_end = float(tt[g0_region].max()) if g0_region.any() else 0.0
for ax in axes:
    ax.axvspan(0, g0_end, color='#fbeee6', alpha=0.5, zorder=0)
    ax.axvline(0, color='k', ls='-', alpha=0.3, lw=0.8)
axes[0].text(g0_end/2, axes[0].get_ylim()[1]*0.9, 'transient\nG0', ha='center', fontsize=9,
             fontweight='bold', color='#b9770e')

plt.tight_layout()
plt.savefig('simulations/fig_v44_cycle_anatomy.png', dpi=170, bbox_inches='tight')
plt.savefig('simulations/fig_v44_cycle_anatomy.pdf', bbox_inches='tight')
plt.close()
print("\nSaved: fig_v44_cycle_anatomy.png")
