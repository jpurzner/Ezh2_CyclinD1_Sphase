"""EZH2-inhibition dose-response: how de-repression of CyclinD1 drives the cell cycle.

EZH2i is a graded input (0 = no inhibition, 1 = full): the CyclinD1 repression factor is
K/(K + EZH2*(1-EZH2i)), so intermediate EZH2i = partial catalytic inhibition (a Tazemetostat dose).
We sweep EZH2i in four contexts and trace CyclinD1 (transcript + protein) -> proliferation:

  GNP            : SHH-driven, already cycling      -> effect on a cycling cell
  MB             : Ptch+/-, MYCN-amplified, cycling -> effect on a tumour that already cycles
  MB + HHi       : vismodegib-arrested              -> the RESCUE: threshold to re-enter the cycle
  MB + CDK4/6i   : palbociclib-arrested (downstream)-> control: cannot be rescued at any dose

Outputs (fig_v44_ezh2i_doseresponse.png) + a printed table:
  (a) CyclinD1 transcript & protein vs EZH2i dose
  (b) proliferation (divisions / 168 h) vs EZH2i dose  -> the de-repression THRESHOLD that the paper
      predicts (low dose may be insufficient or harmful; a critical dose is needed to rescue)
  (c) transient-G0 fraction (p27) vs EZH2i dose
  (d) where each context sits relative to the cycling threshold (CyclinD1 protein vs divisions)

Run:  ./venv/bin/python simulations/fig_v44_ezh2i_doseresponse.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, matplotlib.pyplot as plt
import tellurium as te
from simulations.validate_v44 import count_divisions, mean_settled, classify, P27_THR, SEL, _MODEL

T_END, N_PTS = 10080, 20160          # 168 h
DOSES = np.round(np.linspace(0, 1, 21), 3)


def run(shh=0.5, ptch1_cn=1.0, gdc=0.0, ezh2i=0.0, mycn_amp=1.0, cdk46i=False,
        p16=0.0, p18=None, ksyp21=None, t_end=T_END, n_pts=N_PTS):
    """Robust single-condition run (retry on the flaky CVODE init error with looser tol)."""
    for te_end, te_pts in ((t_end, n_pts), (8000, 32000), (6000, 24000)):
        for atol in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
            rr = te.loada(_MODEL)
            rr.integrator.setValue("absolute_tolerance", atol)
            rr.integrator.setValue("relative_tolerance", 1e-6)
            try: rr.integrator.setValue("maximum_num_steps", 300000)
            except Exception: pass
            rr['SHH'] = shh; rr['Ptch1_copy_number'] = ptch1_cn; rr['GDC0449'] = gdc
            rr['EZH2i'] = ezh2i; rr['MYCN_amplification'] = mycn_amp
            rr['p16'] = p16                                  # INK4 (Cdkn2a) competitive CDK4/6 brake
            if p18 is not None: rr['p18'] = p18              # INK4 (Cdkn2c), GNP default 0.4
            if ksyp21 is not None: rr['kSyP21'] = ksyp21     # CIP/KIP (p21/p27) CDK2 brake
            if cdk46i:
                rr['kPhRbCd'] = 0.0
            try:
                return rr.simulate(0, te_end, te_pts, selections=SEL)
            except Exception:
                continue
    return None     # very high-CyclinD1 (high EZH2i in MB) can be too stiff -> caller carries NaN

_MB = dict(p16=0.15, p18=1.5, ksyp21=0.004)              # MB CDK-inhibitor brake (INK4 + CIP/KIP)
CONTEXTS = {
    'GNP':          dict(shh=0.5, ptch1_cn=1.0, gdc=0.0, mycn_amp=1.0),
    'MB':           dict(shh=0.5, ptch1_cn=0.1, gdc=0.0, mycn_amp=2.8, **_MB),
    'MB + HHi':     dict(shh=0.5, ptch1_cn=0.1, gdc=1.0, mycn_amp=2.8, **_MB),
    'MB + CDK4/6i': dict(shh=0.5, ptch1_cn=0.1, gdc=0.0, mycn_amp=2.8, cdk46i=True, **_MB),
}
COL = {'GNP': '#1b9e77', 'MB': '#762A83', 'MB + HHi': '#E08214', 'MB + CDK4/6i': '#C2185B'}

data = {c: dict(cdm=[], cd=[], ndiv=[], period=[], g0=[]) for c in CONTEXTS}
for cname, cond in CONTEXTS.items():
    for d in DOSES:
        r = run(**cond, ezh2i=float(d), t_end=T_END, n_pts=N_PTS)
        if r is None:                                  # stiff (very high CyclinD1) -> record NaN
            for k in ('cdm', 'cd', 'period', 'g0'):
                data[cname][k].append(np.nan)
            data[cname]['ndiv'].append(np.nan)
            continue
        n, _, per = count_divisions(r)
        f, _, _ = classify(r, P27_THR)   # f = COUNT-fraction (flow quiescent G0)
        data[cname]['cdm'].append(mean_settled(r, 'Cd_mRNA'))
        data[cname]['cd'].append(mean_settled(r, 'Cd'))
        data[cname]['ndiv'].append(n)
        data[cname]['period'].append(np.mean(per) if len(per) else np.nan)
        data[cname]['g0'].append(f['G0'])

# threshold (first dose that cycles) for the arrested contexts
def threshold(c):
    nd = data[c]['ndiv']
    for d, n in zip(DOSES, nd):
        if n > 0:
            return d
    return None

print("=" * 78)
print("EZH2i dose-response (divisions over 168 h; CyclinD1 transcript fold vs EZH2i=0)")
for c in CONTEXTS:
    base = data[c]['cdm'][0]
    print(f"\n{c}:")
    print(f"  {'EZH2i':>6} {'CycD1mRNA':>10} {'fold':>6} {'CycD1prot':>10} {'div/168h':>9} {'G0%':>6}")
    for i, d in enumerate(DOSES):
        if i % 2:  # print every other for brevity
            continue
        print(f"  {d:>6.2f} {data[c]['cdm'][i]:>10.3f} {data[c]['cdm'][i]/base:>6.2f} "
              f"{data[c]['cd'][i]:>10.3f} {data[c]['ndiv'][i]:>9d} {data[c]['g0'][i]:>6.1f}")
    th = threshold(c)
    print(f"  -> cycling threshold: EZH2i = {th}" if th is not None and data[c]['ndiv'][0] == 0
          else ("  -> cycles at all doses" if data[c]['ndiv'][0] > 0 else "  -> never cycles (no rescue)"))

# ---------------- figure ----------------
fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))

# (a) CyclinD1 transcript (solid) + protein (dashed), fold vs EZH2i=0
for c in CONTEXTS:
    base_m = data[c]['cdm'][0]; base_p = data[c]['cd'][0]
    ax[0, 0].plot(DOSES, np.array(data[c]['cdm']) / base_m, color=COL[c], lw=2, label=c)
    ax[0, 0].plot(DOSES, np.array(data[c]['cd']) / base_p, color=COL[c], lw=1.3, ls='--', alpha=0.7)
ax[0, 0].set_xlabel('EZH2i dose (0 = none, 1 = full)'); ax[0, 0].set_ylabel('CyclinD1 fold (vs EZH2i=0)')
ax[0, 0].set_title('(a) EZH2i de-represses CyclinD1\n(solid = transcript, dashed = protein)', fontweight='bold', fontsize=11)
ax[0, 0].legend(fontsize=8); ax[0, 0].grid(alpha=0.2)

# (b) proliferation vs EZH2i dose -> the rescue threshold
for c in CONTEXTS:
    ax[0, 1].plot(DOSES, data[c]['ndiv'], color=COL[c], lw=2, marker='o', ms=3, label=c)
ax[0, 1].set_xlabel('EZH2i dose'); ax[0, 1].set_ylabel('divisions / 168 h')
ax[0, 1].set_title('(b) Proliferation vs EZH2i dose\n(MB+HHi: threshold rescue; CDK4/6i: no rescue)', fontweight='bold', fontsize=11)
ax[0, 1].legend(fontsize=8); ax[0, 1].grid(alpha=0.2)
th = threshold('MB + HHi')
if th is not None and data['MB + HHi']['ndiv'][0] == 0:
    ax[0, 1].axvline(th, color=COL['MB + HHi'], ls=':', alpha=0.6)
    ax[0, 1].text(th, ax[0, 1].get_ylim()[1]*0.5, f' rescue\n threshold\n EZH2i={th:g}', fontsize=8, color=COL['MB + HHi'])

# (c) transient-G0 (p27) vs EZH2i dose
for c in CONTEXTS:
    ax[1, 0].plot(DOSES, data[c]['g0'], color=COL[c], lw=2, marker='s', ms=3, label=c)
ax[1, 0].set_xlabel('EZH2i dose'); ax[1, 0].set_ylabel('transient G0 (p27) %')
ax[1, 0].set_title('(c) EZH2i shrinks the transient G0\n(cells leave quiescence)', fontweight='bold', fontsize=11)
ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=0.2)

# (d) CyclinD1 protein vs divisions -> the cycling threshold the cell must cross
for c in CONTEXTS:
    ax[1, 1].scatter(data[c]['cd'], data[c]['ndiv'], color=COL[c], s=28, label=c, edgecolor='none')
ax[1, 1].set_xlabel('CyclinD1 protein (mean, a.u.)'); ax[1, 1].set_ylabel('divisions / 168 h')
ax[1, 1].set_title('(d) Proliferation is set by CyclinD1 protein\n(EZH2i moves cells along this axis)', fontweight='bold', fontsize=11)
ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=0.2)

fig.suptitle('Impact of EZH2 inhibition (CyclinD1 de-repression) on the cell cycle — v44 model',
             fontsize=13, fontweight='bold', y=1.0)
plt.tight_layout()
plt.savefig('simulations/fig_v44_ezh2i_doseresponse.png', dpi=170, bbox_inches='tight')
plt.savefig('simulations/fig_v44_ezh2i_doseresponse.pdf', bbox_inches='tight')
plt.close()
print("\nSaved: fig_v44_ezh2i_doseresponse.png")
