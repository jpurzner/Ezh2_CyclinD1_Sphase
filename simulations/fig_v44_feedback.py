"""Supp Fig 8 A-H: impact of the EZH2 -> CyclinD1 feedback (v44 model).

Compares the model WITH feedback (EZH2 represses CyclinD1) vs WITHOUT (repression removed by
K_EZH2_repression -> 1e6) in GNP+SHH and Ptch+/- MB contexts. Panels:
  A,B  Cyclin B traces +/- feedback (GNP, MB)
  C,D  CyclinD1 mRNA +/- feedback   (GNP, MB)
  E,F  EZH2 protein +/- feedback    (GNP, MB)
  G,H  normalized CyclinD1 / EZH2 / CycB phase relationship (with feedback, last 70 h)
Also a bar summary (period / CyclinD1 / divisions).

NOTE (v44 vs v42): in v44 the feedback strongly sets CyclinD1 LEVEL but, because cell-cycle
phases are largely uncoupled (explicit-replication core), it does NOT lengthen the period the
way the eps-clock v42 did. This script reports whatever the calibrated model produces.

Run:  ./venv/bin/python simulations/fig_v44_feedback.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.signal import find_peaks
import tellurium as te
from src.build_model_v44_heldt import build_model_v44
from simulations.validate_v44 import PARAMS

T_END, N_PTS = 10080, 20160      # 168 h

CTX = {
    'GNP + SHH':    dict(SHH=0.5, Ptch1_copy_number=1.0, MYCN_amplification=1.0),
    'Ptch+/- MB':   dict(SHH=0.5, Ptch1_copy_number=0.1, MYCN_amplification=2.8),
}


def run(ctx, feedback=True):
    p = dict(PARAMS)
    if not feedback:
        p["K_EZH2_repression"] = 1e6        # ablate EZH2 -> CyclinD1 repression
    rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True, params=p or None))
    rr.integrator.setValue("absolute_tolerance", 1e-9)
    rr.integrator.setValue("relative_tolerance", 1e-6)
    for k, v in CTX[ctx].items():
        rr[k] = v
    return rr.simulate(0, T_END, N_PTS,
                       selections=["time", "Cb", "Cd_mRNA", "EZH2", "MPF"])


def metrics(r):
    t = r['time']; m = t >= 3000; tt = t[m]; dt = tt[1] - tt[0]
    pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=int(200 / dt))
    per = np.diff(tt[pk] / 60.0) if len(pk) > 1 else np.array([])
    idx = -3000
    return dict(ndiv=len(pk), period=(np.mean(per) if len(per) else None),
                cd_mean=float(np.mean(r['Cd_mRNA'][idx:])))


sims = {(c, fb): run(c, fb) for c in CTX for fb in (True, False)}
mets = {k: metrics(v) for k, v in sims.items()}

ctxs = list(CTX)
col_with, col_without = '#2c3e50', '#e74c3c'

fig = plt.figure(figsize=(15, 14))
gs = GridSpec(4, 2, figure=fig, hspace=0.45, wspace=0.28)
panel = [['A', 'B'], ['C', 'D'], ['E', 'F'], ['G', 'H']]

cd_ymax = max(sims[(c, fb)]['Cd_mRNA'].max() for c in ctxs for fb in (True, False)) * 1.08

for col, ctx in enumerate(ctxs):
    rw, ro = sims[(ctx, True)], sims[(ctx, False)]
    mw, mo = mets[(ctx, True)], mets[(ctx, False)]
    th = rw['time'] / 60.0

    ax = fig.add_subplot(gs[0, col])
    lw = f"with fb ({mw['ndiv']} div" + (f", {mw['period']:.1f} h)" if mw['period'] else ")")
    lo = f"no fb ({mo['ndiv']} div" + (f", {mo['period']:.1f} h)" if mo['period'] else ")")
    ax.plot(th, rw['Cb'], col_with, lw=1.0, label=lw)
    ax.plot(ro['time'] / 60.0, ro['Cb'], col_without, lw=1.0, ls='--', label=lo)
    ax.set_title(f"({panel[0][col]}) {ctx}: Cyclin B", fontsize=11, fontweight='bold')
    ax.set_ylabel('Cyclin B'); ax.set_xlabel('Time (h)'); ax.legend(fontsize=8); ax.grid(alpha=0.2)

    ax = fig.add_subplot(gs[1, col])
    ax.plot(th, rw['Cd_mRNA'], col_with, lw=1.0, label='with feedback')
    ax.plot(ro['time'] / 60.0, ro['Cd_mRNA'], col_without, lw=1.0, ls='--', label='no feedback')
    ax.set_title(f"({panel[1][col]}) {ctx}: CyclinD1 mRNA  "
                 f"({mo['cd_mean']/mw['cd_mean']:.1f}x without fb)", fontsize=11, fontweight='bold')
    ax.set_ylabel('CyclinD1 mRNA'); ax.set_xlabel('Time (h)'); ax.set_ylim(0, cd_ymax)
    ax.legend(fontsize=8); ax.grid(alpha=0.2)

    ax = fig.add_subplot(gs[2, col])
    ax.plot(th, rw['EZH2'], col_with, lw=1.0, label='with feedback')
    ax.plot(ro['time'] / 60.0, ro['EZH2'], col_without, lw=1.0, ls='--', label='no feedback')
    ax.set_title(f"({panel[2][col]}) {ctx}: EZH2 protein", fontsize=11, fontweight='bold')
    ax.set_ylabel('EZH2 protein'); ax.set_xlabel('Time (h)'); ax.set_ylim(bottom=0)
    ax.legend(fontsize=8); ax.grid(alpha=0.2)

    ax = fig.add_subplot(gs[3, col])
    zoom = 70
    mk = rw['time'] / 60.0 >= (th[-1] - zoom)
    tp = th[mk]
    def nrm(a): a = a[mk]; return (a - a.min()) / (a.max() - a.min() + 1e-10)
    ax.plot(tp, nrm(rw['Cd_mRNA']), '#e67e22', lw=1.4, label='CyclinD1 mRNA')
    ax.plot(tp, nrm(rw['EZH2']), '#8e44ad', lw=1.4, label='EZH2')
    ax.plot(tp, nrm(rw['Cb']), '#2980b9', lw=1.0, alpha=0.6, label='CycB')
    ax.set_title(f"({panel[3][col]}) {ctx}: phase relationship (last {zoom} h)",
                 fontsize=11, fontweight='bold')
    ax.set_ylabel('Normalized (0-1)'); ax.set_xlabel('Time (h)'); ax.legend(fontsize=8); ax.grid(alpha=0.2)

fig.suptitle('Supp Fig 8 A-H: EZH2 -> CyclinD1 feedback in the v44 model',
             fontsize=14, fontweight='bold', y=1.005)
plt.savefig('simulations/fig_v44_feedback.png', dpi=200, bbox_inches='tight')
plt.savefig('simulations/fig_v44_feedback.pdf', bbox_inches='tight')
plt.close()
print("Saved: fig_v44_feedback.png/pdf")

print("\nEZH2 feedback impact (v44):")
for c in ctxs:
    mw, mo = mets[(c, True)], mets[(c, False)]
    pw = f"{mw['period']:.1f}h" if mw['period'] else "arrest"
    po = f"{mo['period']:.1f}h" if mo['period'] else "arrest"
    print(f"  {c:14s}: period {pw} -> {po} (no fb);  "
          f"CyclinD1 {mw['cd_mean']:.2f} -> {mo['cd_mean']:.2f} a.u. "
          f"({mo['cd_mean']/mw['cd_mean']:.2f}x);  div {mw['ndiv']} -> {mo['ndiv']}")
