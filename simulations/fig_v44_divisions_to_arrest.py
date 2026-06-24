"""Figure: MB's reduced vismodegib sensitivity — a Gli1 residual + a resistant proliferating fraction.

After vismodegib, a GNP's Gli1 collapses (Smo blocked, nothing sustains it), but MB retains a substantial
Gli1 RESIDUAL: the Gli autoregulatory epigenetic memory (Gli1_epi), charged while the pathway was active,
discharges only over ~a day. Data (Purzner-lab bulk RNA-seq, 24 h vismo): MB Gli1 413/17881 = 2.3% of
baseline vs GNP 10/2591 = 0.4% -- a 41x larger residual. MYCN floors CyclinD1 but has NO path to Gli1, so
this residual is Gli-intrinsic and only the memory reproduces it (the no-memory model gives Gli1~0 by 24 h).

Because the residual Gli1 is small, CyclinD1 after vismo stays MYCN-floored, so proliferation is governed
by the threshold/heterogeneity (a RESISTANT FRACTION, ~matching pRb 1/4), not by the Gli1 residual -- the
memory explains the transcript residual without inflating the cycling fraction.

Panel A: Gli1 after vismo (GNP vs MB) vs the 24 h data residuals. Panel B: % of cells still proliferating.
Run (background):  ./venv/bin/python simulations/fig_v44_divisions_to_arrest.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

N = 120
T1, TEND = 8000, 30000
KTL0 = 0.75
rng = np.random.default_rng(7)
cd_scale = np.clip(np.exp(rng.normal(0.0, 0.68, N)), 0.25, 4.0)
p21_div = np.clip(np.exp(rng.normal(np.log(0.85), 0.15, N)), 0.5, 2.0)

CTX = {
    'GNP': dict(f=1.0, mycn=1.0, p16=0.0,  p18=0.4, ksy=0.0020, col='#1b9e77'),
    'MB':  dict(f=0.1, mycn=2.8, p16=0.15, p18=1.5, ksy=0.0040, col='#762A83'),
}
DATA_RESID = {'GNP': 10 / 2591, 'MB': 413 / 17881}   # 24 h vismo Gli1 / baseline (RNA-seq)
_RR = te.loada(build_model_v44())
_RR.integrator.setValue("relative_tolerance", 1e-6)


def withdraw(c, ktl, p21d, sel):
    for atol in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
        _RR.reset()
        _RR['SHH'] = 0.5; _RR['Ptch1_copy_number'] = c['f']; _RR['MYCN_amplification'] = c['mycn']
        _RR['p16'] = c['p16']; _RR['p18'] = c['p18']; _RR['kSyP21'] = c['ksy']
        _RR['EZH2i'] = 0; _RR['HHi'] = 0.0; _RR['P21_div'] = p21d; _RR['k_Cd_translation'] = ktl
        _RR.integrator.setValue("absolute_tolerance", atol)
        try: _RR.integrator.setValue("maximum_num_steps", 300000)
        except Exception: pass
        try:
            s1 = _RR.simulate(0, T1, T1 * 3, selections=sel)
            _RR['HHi'] = 1.0
            s2 = _RR.simulate(T1, TEND, (TEND - T1) * 3, selections=sel)
        except Exception:
            continue
        return {k: np.concatenate([s1[k], s2[k]]) for k in sel}
    return None


def resistant_fraction(c):
    res = []
    for i in range(N):
        r = withdraw(c, KTL0 * cd_scale[i], p21_div[i], ["time", "MPF"])
        if r is None:
            continue
        dt = np.median(np.diff(r['time'])); pk, _ = find_peaks(r['MPF'], prominence=0.15, distance=int(200 / dt))
        res.append(bool(np.any(r['time'][pk] >= TEND - 2200)))
    return 100 * np.mean(res) if res else 0.0


fracs = {}
for name, c in CTX.items():
    fracs[name] = resistant_fraction(c)
    print(f"{name}: resistant (still proliferating) = {fracs[name]:.0f}%")

# representative Gli1 trajectory + its baseline, per context
traj, g0 = {}, {}
for name, c in CTX.items():
    r = withdraw(c, KTL0 * 1.1, 0.85, ["time", "Gli1"])
    traj[name] = r
    g0[name] = float(np.mean(r['Gli1'][(r['time'] >= 5000) & (r['time'] <= T1)])) if r else np.nan
    if r is not None:
        i24 = np.argmin(np.abs(r['time'] - (T1 + 1440)))
        print(f"  {name} Gli1(24h)/baseline = {r['Gli1'][i24] / g0[name]:.4f}  (data {DATA_RESID[name]:.4f})")

# ======================================================================
fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.5, 5.0), gridspec_kw=dict(width_ratios=[1.7, 1.0]))

# Panel A: Gli1 (fraction of baseline) after vismodegib, log scale
for name, c in CTX.items():
    r = traj[name]
    if r is None:
        continue
    th = (r['time'] - T1) / 60.0; m = th >= -30
    axA.plot(th[m], np.clip(r['Gli1'][m] / g0[name], 1e-4, None), color=c['col'], lw=2.0, label=f"{name} (model)", zorder=3)
    axA.scatter([24], [DATA_RESID[name]], color=c['col'], s=90, edgecolor='k', zorder=5, marker='D')
axA.axvline(0, color="#444", ls=":", lw=1.2); axA.axvline(24, color="#999", ls=":", lw=1.0)
axA.text(2, 0.0008, "vismodegib", fontsize=8.5, color="#444")
axA.text(24.5, 0.0008, "24 h\n(data)", fontsize=8, color="#666")
axA.set_yscale("log"); axA.set_xlim(-30, 120); axA.set_ylim(3e-4, 1.5)
axA.set_xlabel("time since vismodegib (h)", fontsize=10.5)
axA.set_ylabel("Gli1 / baseline (log)", fontsize=10.5)
axA.set_title("A   MB retains a Gli1 residual after vismo (Gli memory);\nGNP Gli1 collapses  (◆ = 24 h RNA-seq)",
              fontsize=10.5, fontweight="bold", loc="left")
axA.legend(fontsize=9.5, loc="upper right"); axA.grid(alpha=0.15, which="both")

# Panel B: resistant proliferating fraction
xb = np.arange(2); vals = [fracs['GNP'], fracs['MB']]
axB.bar(xb, vals, color=[CTX['GNP']['col'], CTX['MB']['col']], edgecolor='k', width=0.6)
for x, v in zip(xb, vals):
    axB.text(x, v + 1.5, f"{v:.0f}%", ha='center', fontweight='bold', fontsize=12)
axB.axhline(25, color='#c0392b', ls='--', lw=1.2)
axB.text(1.4, 27, "pRb ~1/4", color='#c0392b', fontsize=8.5, ha='right')
axB.set_xticks(xb); axB.set_xticklabels(['GNP', 'MB'], fontsize=11, fontweight='bold')
axB.set_ylabel("% still proliferating after vismodegib", fontsize=10.5); axB.set_ylim(0, 100)
axB.set_title("B   Proliferation is MYCN-floored\n(resistant fraction, not the Gli1 residual)",
              fontsize=10.5, fontweight="bold", loc="left")
axB.grid(axis='y', alpha=0.2)

fig.suptitle("MB's reduced vismodegib sensitivity: a Gli-intrinsic Gli1 residual + a resistant fraction (v44)",
             fontsize=12, fontweight='bold', y=1.0)
fig.text(0.5, -0.02,
         "MYCN floors CyclinD1 but cannot touch Gli1, so the 41x-larger MB Gli1 residual is Gli-intrinsic "
         "(the autoregulatory memory) and matches the 24 h data. Because that residual is small, CyclinD1 "
         "stays MYCN-floored -> proliferation is a resistant fraction (~pRb 1/4), set by threshold/heterogeneity.",
         ha='center', fontsize=8.0, color='#555', style='italic')
plt.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(os.path.dirname(__file__), "fig_v44_divisions_to_arrest")
fig.savefig(out + ".png", dpi=160, bbox_inches='tight')
fig.savefig(out + ".pdf", bbox_inches='tight')
print("\nwrote", out + ".png /.pdf")
