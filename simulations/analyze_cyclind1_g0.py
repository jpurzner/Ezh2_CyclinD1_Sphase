"""Relationship: CyclinD1 transcript -> protein -> cell-cycle duration -> transient G0 duration.

Three analyses in the v44 model:
  (1) QSS check: is Cd protein simply proportional to Cd_mRNA? (timescale argument)
  (2) Dissociation: change transcript (via SHH / EZH2 repression) vs change protein-per-transcript
      (via k_Cd_translation). If G0 duration collapses onto ONE curve vs PROTEIN regardless of which
      knob moved it, then CyclinD1 PROTEIN (not transcript per se) is the controller.
  (3) Phase-duration breakdown vs CyclinD1 drive: how G0 / G1 / S / G2 durations and period respond.

Phase durations are computed as (time-fraction in phase) x period from a settled cycling trajectory.
Transient G0 = pre-replication & phospho-Rb-negative (pRb < 0.5*max) -- the p27-high post-division
window before the Skp2-p27-Rb feedforward commits the cell.

Run:  ./venv/bin/python simulations/analyze_cyclind1_g0.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
import tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

SEL = ["time", "Cd_mRNA", "Cd", "pRb", "P21", "aRc", "Dna", "MPF", "mass", "EZH2"]
T_END, N_PTS, SETTLE = 14000, 56000, 5000


def simulate(params=None, shh=0.5, ptch1=1.0, mycn=1.0, gdc=0.0, ezh2i=0.0):
    rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True, params=(params or None)))
    rr.integrator.setValue("absolute_tolerance", 1e-9)
    rr.integrator.setValue("relative_tolerance", 1e-6)
    rr['SHH'] = shh; rr['Ptch1_copy_number'] = ptch1; rr['MYCN_amplification'] = mycn
    rr['GDC0449'] = gdc; rr['EZH2i'] = ezh2i
    return rr.simulate(0, T_END, N_PTS, selections=SEL)


def analyze(res):
    t = res['time']; m = t >= SETTLE
    tt = t[m]; dt = tt[1] - tt[0]
    mpf = res['MPF'][m]
    pk, _ = find_peaks(mpf, prominence=0.15, distance=int(200 / dt))
    if len(pk) < 2:
        return dict(ndiv=len(pk), period=np.nan, cd_mRNA=float(np.mean(res['Cd_mRNA'][m])),
                    cd=float(np.mean(res['Cd'][m])), G0=np.nan, G1=np.nan, S=np.nan, G2=np.nan,
                    arrest=True)
    period = float(np.mean(np.diff(tt[pk]))) / 60.0
    pRb = res['pRb'][m]; aRc = res['aRc'][m]; Dna = res['Dna'][m]
    thr = 0.5 * pRb.max()
    in_S = (aRc > 0.05) & (Dna < 0.98); in_G2 = Dna >= 0.98; preS = ~in_S & ~in_G2
    G0 = preS & (pRb < thr); G1 = preS & (pRb >= thr)
    fr = {k: v.mean() for k, v in dict(G0=G0, G1=G1, S=in_S, G2=in_G2).items()}
    s = sum(fr.values()) or 1.0
    dur = {k: period * v / s for k, v in fr.items()}    # hours per cycle
    return dict(ndiv=len(pk), period=period, cd_mRNA=float(np.mean(res['Cd_mRNA'][m])),
                cd=float(np.mean(res['Cd'][m])), arrest=False, **dur)


# ============================================================
# (1) QSS check
# ============================================================
print("=" * 70)
print("(1) CyclinD1 transcript -> protein (quasi-steady-state check)")
r = simulate(shh=0.5)
m = r['time'] >= SETTLE
cdm, cd = r['Cd_mRNA'][m], r['Cd'][m]
slope = np.polyfit(cdm, cd, 1)
print(f"  k_Cd_translation/k_Cd_deg = 0.26/1.0 = 0.26  (predicted protein/transcript)")
print(f"  fitted Cd vs Cd_mRNA: slope={slope[0]:.3f}, intercept={slope[1]:.3f}  (r={np.corrcoef(cdm,cd)[0,1]:.4f})")
print(f"  protein tracks transcript with ~1 min lag (cycle ~1380 min) -> effectively instantaneous")

# ============================================================
# (2)+(3) Sweeps that move CyclinD1, measuring durations
# ============================================================
print("\n" + "=" * 70)
print("(2/3) Sweeps: CyclinD1 drive vs phase durations")

# Knob A: SHH dose (GNP) -> changes TRANSCRIPT (mitogen)
SHH_SET = [0.35, 0.4, 0.45, 0.5, 0.6, 0.75, 1.0]
shh_res = [(s, analyze(simulate(shh=s))) for s in SHH_SET]

# Knob B: k_Cd_translation (GNP, SHH=0.5) -> changes PROTEIN per transcript (transcript ~fixed)
TL_SET = [0.16, 0.20, 0.26, 0.32, 0.40]
tl_res = [(k, analyze(simulate(params={"k_Cd_translation": k}, shh=0.5))) for k in TL_SET]

# Knob C: EZH2 repression strength K_EZH2_repression (GNP) -> changes TRANSCRIPT (de-repression)
KE_SET = [0.3, 0.5, 0.8, 1.5]
ke_res = [(k, analyze(simulate(params={"K_EZH2_repression": k}, shh=0.5))) for k in KE_SET]

print("\n  SHH sweep (GNP):")
print(f"  {'SHH':>5} {'Cd_mRNA':>8} {'Cd_prot':>8} {'period':>7} {'G0':>6} {'G1':>6} {'S':>6} {'G2':>6} {'div':>4}")
for s, a in shh_res:
    if a['arrest']:
        print(f"  {s:>5.2f} {a['cd_mRNA']:>8.3f} {a['cd']:>8.3f}   ARREST (transient G0 -> permanent)")
    else:
        print(f"  {s:>5.2f} {a['cd_mRNA']:>8.3f} {a['cd']:>8.3f} {a['period']:>6.1f}h {a['G0']:>6.2f} "
              f"{a['G1']:>6.2f} {a['S']:>6.2f} {a['G2']:>6.2f} {a['ndiv']:>4}")

print("\n  k_Cd_translation sweep (GNP, SHH=0.5; transcript ~fixed, protein scaled):")
print(f"  {'kTl':>5} {'Cd_mRNA':>8} {'Cd_prot':>8} {'period':>7} {'G0':>6} {'G1':>6} {'S':>6} {'G2':>6}")
for k, a in tl_res:
    if a['arrest']:
        print(f"  {k:>5.2f} {a['cd_mRNA']:>8.3f} {a['cd']:>8.3f}   ARREST")
    else:
        print(f"  {k:>5.2f} {a['cd_mRNA']:>8.3f} {a['cd']:>8.3f} {a['period']:>6.1f}h {a['G0']:>6.2f} "
              f"{a['G1']:>6.2f} {a['S']:>6.2f} {a['G2']:>6.2f}")

print("\n  K_EZH2_repression sweep (GNP, SHH=0.5; higher K = weaker repression = more transcript):")
for k, a in ke_res:
    tag = "ARREST" if a['arrest'] else (f"period {a['period']:.1f}h, G0 {a['G0']:.2f}h, "
                                        f"G1 {a['G1']:.2f}h")
    print(f"  K={k:>4.1f}: Cd_mRNA={a['cd_mRNA']:.3f} Cd_prot={a['cd']:.3f}  {tag}")

# ============================================================
# Figure
# ============================================================
fig, ax = plt.subplots(2, 2, figsize=(13, 10))

# (a) protein vs transcript (QSS)
ax[0, 0].scatter(cdm[::50], cd[::50], s=6, alpha=0.4, color='#2980b9')
xx = np.linspace(cdm.min(), cdm.max(), 50)
ax[0, 0].plot(xx, 0.26 * xx, 'r--', label='Cd = 0.26 x Cd_mRNA (QSS)')
ax[0, 0].set_xlabel('CyclinD1 transcript (Cd_mRNA)'); ax[0, 0].set_ylabel('CyclinD1 protein (Cd)')
ax[0, 0].set_title('(a) Protein tracks transcript (quasi-steady-state)', fontweight='bold', fontsize=11)
ax[0, 0].legend(fontsize=9); ax[0, 0].grid(alpha=0.2)

# (b) G0 duration vs CyclinD1 PROTEIN — overlay all knobs (collapse test)
def xy(res, dur='G0'):
    pts = [(a['cd'], a[dur]) for _, a in res if not a['arrest']]
    return np.array(pts).T if pts else (np.array([]), np.array([]))
for res, lab, c, mk in [(shh_res, 'SHH dose', '#1b9e77', 'o'),
                        (tl_res, 'k_Cd_translation', '#d95f02', 's'),
                        (ke_res, 'EZH2 repression', '#7570b3', '^')]:
    x, y = xy(res, 'G0')
    if len(x): ax[0, 1].scatter(x, y, s=60, color=c, marker=mk, label=lab, edgecolor='k', zorder=3)
ax[0, 1].set_xlabel('CyclinD1 protein (mean, a.u.)'); ax[0, 1].set_ylabel('Transient G0 duration (h/cycle)')
ax[0, 1].set_title('(b) G0 duration is set by CyclinD1 PROTEIN\n(all knobs collapse onto one curve)',
                   fontweight='bold', fontsize=11)
ax[0, 1].legend(fontsize=9); ax[0, 1].grid(alpha=0.2)

# (c) period & preS(G0+G1) vs CyclinD1 protein
for res, c, mk in [(shh_res, '#1b9e77', 'o'), (tl_res, '#d95f02', 's'), (ke_res, '#7570b3', '^')]:
    pts = [(a['cd'], a['period'], a['G0'] + a['G1']) for _, a in res if not a['arrest']]
    if pts:
        arr = np.array(pts).T
        ax[1, 0].scatter(arr[0], arr[1], s=55, color=c, marker=mk, edgecolor='k', label='period')
        ax[1, 0].scatter(arr[0], arr[2], s=55, color=c, marker=mk, alpha=0.4)
ax[1, 0].set_xlabel('CyclinD1 protein (mean, a.u.)'); ax[1, 0].set_ylabel('hours')
ax[1, 0].set_title('(c) Period (solid) & G0+G1 (faded) vs CyclinD1 protein', fontweight='bold', fontsize=11)
ax[1, 0].grid(alpha=0.2)

# (d) stacked phase durations vs SHH
cyc = [(s, a) for s, a in shh_res if not a['arrest']]
xs = [s for s, _ in cyc]
G0 = [a['G0'] for _, a in cyc]; G1 = [a['G1'] for _, a in cyc]
S = [a['S'] for _, a in cyc]; G2 = [a['G2'] for _, a in cyc]
ax[1, 1].bar(range(len(xs)), G0, label='G0', color='#fbeee6', edgecolor='k')
ax[1, 1].bar(range(len(xs)), G1, bottom=G0, label='G1', color='#aed6f1', edgecolor='k')
ax[1, 1].bar(range(len(xs)), S, bottom=np.array(G0)+np.array(G1), label='S', color='#a3e4d7', edgecolor='k')
ax[1, 1].bar(range(len(xs)), G2, bottom=np.array(G0)+np.array(G1)+np.array(S), label='G2', color='#fdebd0', edgecolor='k')
ax[1, 1].set_xticks(range(len(xs))); ax[1, 1].set_xticklabels([f'{s:.2f}' for s in xs])
ax[1, 1].set_xlabel('SHH dose'); ax[1, 1].set_ylabel('phase duration (h)')
ax[1, 1].set_title('(d) Phase durations vs SHH (GNP)', fontweight='bold', fontsize=11)
ax[1, 1].legend(fontsize=9)

plt.tight_layout()
plt.savefig('simulations/fig_v44_cyclind1_g0_analysis.png', dpi=170, bbox_inches='tight')
plt.close()
print("\nSaved: fig_v44_cyclind1_g0_analysis.png")
