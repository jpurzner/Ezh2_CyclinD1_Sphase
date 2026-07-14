"""Two-step Rb switch restores hyper-P Rb as a faithful G0/R-point marker (v44).

The calibrated v44 DEFAULT is a TWO-STEP Rb switch:
  * CyclinD-CDK4/6 PRIMES Rb -> Rbm (mono-phospho), a size-gated, E2f-still-bound state that
    accumulates through the pre-commitment (p27-high) window WITHOUT firing the switch.
  * CyclinE/A-CDK2 COMMITS Rbm -> pRb (hyper-phospho), releasing E2f at the R-point.
The old single-step framing (with_two_step_rb=False) collapsed both into one pRb pool, so any kinase
drove pRb to max within ~20 min of anaphase -- long before functional commitment -- making hyper-P Rb
a USELESS G0 marker. The two-step default keeps pRb LOW through the p27-high G0 window and raises it
only at commitment, so hyper-P-Rb-LOW and p27-HIGH coincide in G0 (a faithful marker again).

(A) one-cycle GNP anatomy on the two-step default: Rb(hypo) -> Rbm(mono-P) -> pRb(hyper-P), with
    p27 and the pre-commitment p27-high window shaded.
(B) single-step vs two-step overlay of hyper-P Rb (pRb) across a cycle: single-step SATURATES within
    ~20 min of anaphase; two-step stays LOW through the p27-high window and rises only at commitment.
(C) marker agreement: fraction of the pre-commit (p27-high) window during which hyper-P Rb is LOW --
    high under two-step (markers AGREE), near zero under single-step (pRb already high; markers DISAGREE).

Run:  ./venv/bin/python simulations/fig_v44_two_step_rpoint.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

T_END, N_PTS, SETTLE = 15000, 60000, 6000
P27_THR = 0.1                    # p27 (P21) commitment threshold (validate classify)
# GNP background (SHH=0.5). The deterministic GNP cell AT the setpoint (P21_div=0.6) is born committed
# (Spencer carryover keeps p27<thr, pRb high) -> it has NO p27-high G0 window to mark. To EXHIBIT the
# marker (i.e. a cell that actually dwells in G0 before committing) we use the transient-G0 birth p27
# used in fig_v44_transient_g0_definition (P21_div=1.6): a GNP cell born p27-high that dwells, then
# commits. Same calibrated defaults; only the Rb-switch toggle differs between the two models compared.
P21_DIV_G0 = 1.6
SEL_TWO = ["time", "Cd", "P21", "E2f", "Ce", "Ca", "aRc", "Dna", "MPF", "Rb", "Rbm", "pRb", "mass"]
SEL_ONE = ["time", "Cd", "P21", "E2f", "Ce", "Ca", "aRc", "Dna", "MPF", "Rb", "pRb", "mass"]


def sim(two_step=True):
    """Fresh model per condition (reset() does NOT restore parameters -> re-build per condition).
    GNP background (SHH=0.5) with the transient-G0 birth p27 (P21_DIV_G0) so the cell has a real
    p27-high G0 window to mark."""
    rr = te.loada(build_model_v44(with_two_step_rb=two_step))
    rr.integrator.setValue("absolute_tolerance", 1e-8)
    rr.integrator.setValue("relative_tolerance", 1e-6)
    try:
        rr.integrator.setValue("maximum_time_step", 5.0)
        rr.integrator.setValue("maximum_num_steps", 2000000)
    except Exception:
        pass
    rr['SHH'] = 0.5
    rr['P21_div'] = P21_DIV_G0
    sel = SEL_TWO if two_step else SEL_ONE
    return rr.simulate(0, T_END, N_PTS, selections=sel)


def one_cycle(res):
    """Extract one representative cycle aligned to birth (anaphase). Division = Dna 1->0 reset
    (NOT MPF peaks -- too sharp for coarse sampling)."""
    t = res['time']; Dna = res['Dna']
    div = np.where((Dna[:-1] > 0.9) & (Dna[1:] < 0.5))[0] + 1
    div = div[t[div] >= SETTLE]
    if len(div) < 3:
        raise RuntimeError(f"need >=3 divisions after settle, got {len(div)}")
    a0, a1 = div[1], div[2]         # a well-settled interior cycle
    seg = slice(a0, a1)
    th = (t[seg] - t[a0]) / 60.0    # hours since anaphase/birth
    return th, seg


def commit_time(th, P21seg):
    """Pre-commitment p27-high window: birth -> first time p27 drops below threshold (Skp2-p27 toggle
    flips, CDK2 freed). Returns commit time in h (or last time if it never drops)."""
    below = np.where(P21seg < P27_THR)[0]
    return float(th[below[0]]) if len(below) else float(th[-1])


def t_half(th, x):
    """Time (h) for x to first reach 50% of its cycle max after birth."""
    xm = np.nanmax(x)
    if xm <= 0:
        return np.nan
    hit = np.where(x >= 0.5 * xm)[0]
    return float(th[hit[0]]) if len(hit) else np.nan


# ---- simulate both models ----
r2 = sim(two_step=True)
r1 = sim(two_step=False)

th2, s2 = one_cycle(r2)
th1, s1 = one_cycle(r1)

# two-step cycle traces
Rb2, Rbm2, pRb2 = r2['Rb'][s2], r2['Rbm'][s2], r2['pRb'][s2]
P21_2, E2f2, aRc2, Cd2 = r2['P21'][s2], r2['E2f'][s2], r2['aRc'][s2], r2['Cd'][s2]

# single-step cycle traces
Rb1, pRb1 = r1['Rb'][s1], r1['pRb'][s1]
P21_1, E2f1, aRc1 = r1['P21'][s1], r1['E2f'][s1], r1['aRc'][s1]

# commitment / p27-high windows
tc2 = commit_time(th2, P21_2)
tc1 = commit_time(th1, P21_1)

# pRb half-rise times (post-anaphase) -- the crux number
th_pRb2 = t_half(th2, pRb2)
th_pRb1 = t_half(th1, pRb1)

def t_frac(th, x, f):
    xm = np.nanmax(x)
    hit = np.where(x >= f * xm)[0]
    return float(th[hit[0]]) if len(hit) and xm > 0 else np.nan
t90_1 = t_frac(th1, pRb1, 0.9)   # single-step: when hyper-P Rb is ~saturated post-anaphase

# marker-agreement metric: fraction of the p27-high pre-commit window with hyper-P Rb LOW (< half-max)
def agreement(th, pRb, tc):
    win = th <= tc
    if win.sum() == 0:
        return np.nan
    low = pRb < 0.5 * np.nanmax(pRb)
    return 100.0 * float(np.mean(low[win]))

agr2 = agreement(th2, pRb2, tc2)
agr1 = agreement(th1, pRb1, tc1)

# ============================ FIGURE ============================
fig, ax = plt.subplots(1, 3, figsize=(17, 5.4))
fig.suptitle("Two-step Rb switch restores hyper-P Rb as a faithful G0 / R-point marker (v44 default)",
             fontsize=14, fontweight="bold", y=1.00)
fig.text(0.5, 0.925, "CyclinD-CDK4/6 primes Rb→Rbm (mono-P); CyclinE/A-CDK2 commits Rbm→pRb (hyper-P) at "
         "the R-point.  GNP transient-G0 cell.  Single-step pRb is high from anaphase (bad marker); "
         "two-step pRb stays low through G0 (good marker).",
         ha="center", fontsize=9.0, color="#555", style="italic")

# ---- (a) one-cycle two-step anatomy: Rb -> Rbm -> pRb ----
axA = ax[0]
axA.plot(th2, Rb2,  color="#2471a3", lw=2.0, label="Rb (hypo-P, active)")
axA.plot(th2, Rbm2, color="#e67e22", lw=2.2, label="Rbm (mono-P, primed)")
axA.plot(th2, pRb2, color="#c0392b", lw=2.2, label="pRb (hyper-P, committed)")
axA.plot(th2, P21_2, color="#8e44ad", lw=1.8, ls="--", label="p27 (=P21)")
axA.axhline(P27_THR, color="#8e44ad", lw=0.8, ls=":", alpha=0.6)
axA.axvspan(0, tc2, color="#f9e79f", alpha=0.55, zorder=0)
ymax = max(np.nanmax(Rb2), np.nanmax(Rbm2), np.nanmax(pRb2)) * 1.05
axA.text(tc2/2, ymax*0.93, f"pre-commit G0\n(p27-high, ~{tc2:.1f} h)", ha="center",
         fontsize=8.6, color="#8a6d00", fontweight="bold")
axA.axvline(tc2, color="#7f8c8d", lw=1.0, ls="-.", alpha=0.7)
axA.text(tc2, ymax*0.55, " commitment\n (R-point)", fontsize=8, color="#7f8c8d", fontweight="bold")
axA.set_xlabel("time since anaphase (h)"); axA.set_ylabel("level (a.u.)")
axA.set_ylim(-0.03, ymax)
axA.set_title("(a) One GNP cycle (two-step default):\nRb(hypo) → Rbm(mono-P) → pRb(hyper-P)",
              fontsize=10.5, fontweight="bold")
axA.legend(fontsize=8, loc="center right"); axA.grid(alpha=0.2)

# ---- (b) single-step vs two-step hyper-P Rb overlay ----
axB = ax[1]
axB.plot(th1, pRb1/np.nanmax(pRb1), color="#c0392b", lw=2.2, ls="--",
         label="single-step pRb (saturates!)")
axB.plot(th2, pRb2/np.nanmax(pRb2), color="#c0392b", lw=2.4,
         label="two-step pRb (stays low in G0)")
# shade the two-step p27-high pre-commit window as the "should be pRb-low" zone
axB.axvspan(0, tc2, color="#f9e79f", alpha=0.5, zorder=0)
axB.text(tc2/2, 0.06, "p27-high G0\nwindow", ha="center", fontsize=8.4,
         color="#8a6d00", fontweight="bold")
# annotate half-rise times
axB.axhline(0.5, color="#95a5a6", lw=0.8, ls=":", alpha=0.7)
if not np.isnan(t90_1):
    axB.annotate(f"single-step ~sat\n@ {t90_1*60:.0f} min", xy=(t90_1, 0.9),
                 xytext=(t90_1+2.2, 0.55), fontsize=8, color="#c0392b",
                 arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.0))
if not np.isnan(th_pRb2):
    axB.annotate(f"two-step 50%\n@ {th_pRb2:.1f} h", xy=(th_pRb2, 0.5),
                 xytext=(th_pRb2-4.5, 0.72), fontsize=8, color="#7b241c",
                 arrowprops=dict(arrowstyle="->", color="#7b241c", lw=1.0))
axB.set_xlabel("time since anaphase (h)"); axB.set_ylabel("hyper-P Rb (pRb, normalized)")
axB.set_ylim(-0.03, 1.08)
axB.set_title("(b) Hyper-P Rb across a cycle:\nsingle-step saturates ~20 min post-anaphase vs two-step",
              fontsize=10.5, fontweight="bold")
axB.legend(fontsize=8.4, loc="center right"); axB.grid(alpha=0.2)

# ---- (c) marker agreement ----
axC = ax[2]
labels = ["single-step\n(with_two_step_rb=False)", "two-step\n(default)"]
vals = [agr1, agr2]
cols = ["#c0392b", "#1e8449"]
bars = axC.bar([0, 1], vals, color=cols, edgecolor="k", width=0.6)
for i, v in enumerate(vals):
    axC.text(i, v + 2, f"{v:.0f}%", ha="center", fontsize=12, fontweight="bold", color=cols[i])
axC.axhline(50, color="#95a5a6", lw=0.8, ls=":")
axC.set_xticks([0, 1]); axC.set_xticklabels(labels, fontsize=9)
axC.set_ylim(0, 108)
axC.set_ylabel("% of p27-high G0 window with hyper-P Rb LOW")
axC.set_title("(c) Marker agreement in the pre-commit G0 window:\np27-HIGH ↔ hyper-P-Rb-LOW",
              fontsize=10.5, fontweight="bold")
axC.text(0, 12, "pRb already HIGH\n→ markers DISAGREE", ha="center", fontsize=8.2,
         color="#7b241c", fontweight="bold")
axC.text(1, agr2*0.5, "pRb LOW while\np27 HIGH\n→ markers AGREE", ha="center", fontsize=8.2,
         color="#0e5a2b", fontweight="bold")
axC.grid(alpha=0.2, axis="y")

plt.tight_layout(rect=[0, 0, 1, 0.92])
# NB: fixed figsize + tight_layout, NO bbox_inches='tight' (a near-arrest lineage could place an
# off-axis artist -> tight-bbox would blow up the canvas).
out = os.path.join(os.path.dirname(__file__), "fig_v44_two_step_rpoint")
plt.savefig(out + ".png", dpi=170)
plt.savefig(out + ".pdf")
plt.close()

print("Saved: fig_v44_two_step_rpoint.png /.pdf")
print(f"  cycle length                : two-step {th2[-1]:.1f} h | single-step {th1[-1]:.1f} h")
print(f"  pre-commit p27-high window  : two-step {tc2:.2f} h | single-step {tc1:.2f} h")
print(f"  birth hyper-P Rb (frac max)  : two-step {pRb2[0]/np.nanmax(pRb2):.2f} | "
      f"single-step {pRb1[0]/np.nanmax(pRb1):.2f}")
print(f"  hyper-P Rb 50%-rise time    : two-step {th_pRb2:.2f} h ({th_pRb2*60:.0f} min) | "
      f"single-step {th_pRb1*60:.0f} min")
print(f"  single-step ~saturation (90%): {t90_1*60:.0f} min post-anaphase")
print(f"  marker agreement (% of p27-high window with hyper-P Rb LOW): "
      f"two-step {agr2:.0f}%  vs  single-step {agr1:.0f}%")
