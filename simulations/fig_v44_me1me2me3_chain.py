"""Exploratory extension (NOT baked): explicit serial H3K27 methylation me0->me1->me2->me3 vs the
lumped single-Mk mark. Shows how the slow, rate-limiting me2->me3 step changes the KINETICS of me3
accumulation: a lag, and a strong dependence on cell-cycle period that the lumped mark does not have.

Fractions m0,m1,m2,m3 (sum=1) of K27 residues in a Ccnd1 domain. PRC2 catalyses the chain forward
(activity P); KDM/turnover reverses it. Fast me0->me1, me1->me2; SLOW me2->me3 (rate-limiting).
Replication: at each division every T hours, half the histones are new (me0) -> m_{1,2,3} *= 1/2.

Pure ODE (RK4), controllable period T -> isolates the chain kinetics. Not wired into v44's cell cycle;
this is the methylation MODULE that would replace the lumped Mk. Run: ./venv/bin/python simulations/fig_v44_me1me2me3_chain.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ---- chain rates (per hour). me0->me1, me1->me2 FAST; me2->me3 SLOW (rate-limiting) ----
K1, K2, K3 = 0.5, 0.5, 0.02          # forward methylation (x PRC2 activity P)
D1, D2, D3 = 0.02, 0.02, 0.01        # reverse demethylation / turnover
P0 = 1.0
# lumped one-step comparator: dMk/dt = P0*KW*(1-Mk) - DW*Mk  (matched plateau ~0.67, fast)
KW, DW = 0.30, 0.15
DT = 0.1

def _chain_deriv(y, P):
    m1, m2, m3 = y
    m0 = 1.0 - m1 - m2 - m3
    dm1 = P * K1 * m0 - P * K2 * m1 + D2 * m2 - D1 * m1
    dm2 = P * K2 * m1 - P * K3 * m2 + D3 * m3 - D2 * m2
    dm3 = P * K3 * m2 - D3 * m3
    return np.array([dm1, dm2, dm3])

def simulate(T=24.0, tmax=400.0, divide=True, chain=True):
    """Return t, and either (m1,m2,m3) if chain else Mk. Dilution every T h if divide."""
    n = int(tmax / DT)
    t = np.arange(n) * DT
    if chain:
        Y = np.zeros((n, 3)); y = np.array([0.0, 0.0, 0.0])
    else:
        Y = np.zeros(n); y = 0.0
    next_div = T
    for i in range(n):
        if chain:
            Y[i] = y
        else:
            Y[i] = y
        # RK4 step
        if chain:
            k1 = _chain_deriv(y, P0); k2 = _chain_deriv(y + 0.5 * DT * k1, P0)
            k3 = _chain_deriv(y + 0.5 * DT * k2, P0); k4 = _chain_deriv(y + DT * k3, P0)
            y = y + DT / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
            y = np.clip(y, 0, 1)
        else:
            def dl(mk): return P0 * KW * (1 - mk) - DW * mk
            k1 = dl(y); k2 = dl(y + 0.5 * DT * k1); k3 = dl(y + 0.5 * DT * k2); k4 = dl(y + DT * k3)
            y = y + DT / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        if divide and (i + 1) * DT >= next_div:
            if chain:
                y = y * 0.5
            else:
                y = y * 0.5
            next_div += T
    return (t, Y) if chain else (t, Y)

def settled_mean(T, chain=True, ncyc=40):
    tmax = max(30 * T, 1500)
    t, Y = simulate(T, tmax, divide=True, chain=chain)
    s = t >= tmax - 6 * T
    if chain:
        return Y[s].mean(axis=0)      # (m1,m2,m3)
    return Y[s].mean()

# =================== figure ===================
C1, C2, C3, CMK, GREY = "#f39c12", "#e67e22", "#7b241c", "#2471a3", "#7f8c8d"
fig = plt.figure(figsize=(16.5, 11))
gs = GridSpec(2, 2, figure=fig, hspace=0.32, wspace=0.24, left=0.06, right=0.975, top=0.9, bottom=0.07)
fig.suptitle("Serial H3K27 methylation (me1→me2→me3) vs a lumped mark — how the slow last step changes me3 kinetics",
             fontsize=14.5, fontweight="bold", y=0.965)
fig.text(0.5, 0.925, "me0→me1→me2→me3 with FAST me1/me2 and a SLOW, rate-limiting me2→me3. Every division halves the modified fractions (new histones = me0). "
         "Exploratory module (not baked).", ha="center", fontsize=9.5, color="#555", style="italic")

# ---- A: scheme ----
axA = fig.add_subplot(gs[0, 0]); axA.axis("off"); axA.set_xlim(0, 10); axA.set_ylim(0, 10)
axA.set_title("A.  The serial scheme", fontsize=11, fontweight="bold", loc="left")
xs = [0.6, 3.0, 5.4, 7.8]; labs = ["me0", "me1", "me2", "me3"]; cols = ["#d5dbdb", C1, C2, C3]
for x, l, c in zip(xs, labs, cols):
    axA.add_patch(FancyBboxPatch((x, 6.2), 1.6, 1.1, boxstyle="round,pad=0.05,rounding_size=0.1", fc=c, ec="#2c3e50", lw=1.4))
    axA.text(x + 0.8, 6.75, l, ha="center", va="center", fontsize=11, fontweight="bold",
             color="white" if l in ("me2", "me3") else "#2c3e50")
frates = ["k₁·PRC2\n(fast)", "k₂·PRC2\n(fast)", "k₃·PRC2\n(SLOW,\nrate-limiting)"]
for i in range(3):
    axA.add_patch(FancyArrowPatch((xs[i] + 1.65, 7.0), (xs[i + 1] - 0.05, 7.0), arrowstyle="->",
                  color=("#c0392b" if i == 2 else "#27ae60"), lw=(2.4 if i == 2 else 1.8), mutation_scale=16))
    axA.text((xs[i] + xs[i + 1]) / 2 + 0.8, 7.75, frates[i], ha="center", fontsize=7.6,
             color=("#c0392b" if i == 2 else "#27ae60"), fontweight="bold")
    axA.add_patch(FancyArrowPatch((xs[i + 1] - 0.05, 6.4), (xs[i] + 1.65, 6.4), arrowstyle="->",
                  color=GREY, lw=1.2, mutation_scale=11))
axA.text(5.4, 5.55, "reverse: KDM6B / turnover", ha="center", fontsize=7.6, color=GREY, style="italic")
axA.text(0.4, 4.4, "At each division (S phase): new histones are unmodified →\n"
         "m₁, m₂, m₃  →  ½·m₁, ½·m₂, ½·m₃   (the lumped model's Mk→½Mk)",
         fontsize=9, color="#2c3e50")
axA.text(0.4, 3.0, "Why me3 is special:\n"
         "• me1, me2 equilibrate fast (~2 h) → track current PRC2 activity, high turnover.\n"
         "• me3 sits behind the slow k₃ bottleneck → builds over ~tens of hours,\n"
         "   integrates activity over time → the heritable 'memory' state.",
         fontsize=8.6, color="#333")
axA.add_patch(plt.Rectangle((0.3, 2.6), 9.4, 5.1, fill=False, ec="#ccc", lw=1))

# ---- B: accumulation lag (no division) ----
axB = fig.add_subplot(gs[0, 1])
t, Y = simulate(T=24, tmax=180, divide=False, chain=True)
axB.plot(t, Y[:, 0], color=C1, lw=2, label="me1 (fast)")
axB.plot(t, Y[:, 1], color=C2, lw=2, label="me2")
axB.plot(t, Y[:, 2], color=C3, lw=2.4, label="me3 (slow, lags)")
tm, Mk = simulate(T=24, tmax=180, divide=False, chain=False)
axB.plot(tm, Mk, color=CMK, lw=2, ls="--", label="lumped Mk (one-step)")
axB.set_xlabel("time (h)  —  no division, constant PRC2"); axB.set_ylabel("fraction methylated")
axB.set_ylim(-0.02, 1.0); axB.grid(alpha=0.25); axB.legend(fontsize=8.5, loc="center right")
axB.set_title("B.  Accumulation LAG: me3 builds slowly, sigmoidally", fontsize=10.6, fontweight="bold")
axB.text(0.03, 0.93, "lumped Mk jumps to plateau (exponential);\nme3 lags ~tens of h behind me1/me2",
         transform=axB.transAxes, fontsize=8, va="top", color="#555", style="italic")

# ---- C: across divisions (T=24h) ----
axC = fig.add_subplot(gs[1, 0])
T0 = 22.0
t, Y = simulate(T=T0, tmax=300, divide=True, chain=True)
for x in np.arange(T0, 300, T0):
    axC.axvline(x, color=GREY, lw=0.5, ls=":", alpha=0.5)
axC.plot(t, Y[:, 0], color=C1, lw=1.6, label="me1")
axC.plot(t, Y[:, 1], color=C2, lw=1.6, label="me2")
axC.plot(t, Y[:, 2], color=C3, lw=2.2, label="me3")
tm, Mk = simulate(T=T0, tmax=300, divide=True, chain=False)
axC.plot(tm, Mk, color=CMK, lw=1.8, ls="--", label="lumped Mk")
axC.set_xlabel(f"time (h)  —  dividing every {T0:.0f} h (dotted)"); axC.set_ylabel("fraction methylated")
axC.set_ylim(-0.02, 1.0); axC.grid(alpha=0.25); axC.legend(fontsize=8.5, loc="upper right")
axC.set_title("C.  Under division: me1/me2 recover; me3 held down by dilution", fontsize=10.6, fontweight="bold")
axC.text(0.03, 0.5, "me1/me2 refill each cycle (fast);\nme3 can't finish the slow climb\nbefore the next ÷2 → lower plateau",
         transform=axC.transAxes, fontsize=8, color="#555", style="italic")

# ---- D: cycle-period dependence (payoff) ----
axD = fig.add_subplot(gs[1, 1])
Ts = np.linspace(6, 96, 16)
M = np.array([settled_mean(T, chain=True) for T in Ts])
Mk = np.array([settled_mean(T, chain=False) for T in Ts])
axD.plot(Ts, M[:, 0], "-o", color=C1, ms=4, lw=1.8, label="me1")
axD.plot(Ts, M[:, 1], "-o", color=C2, ms=4, lw=1.8, label="me2")
axD.plot(Ts, M[:, 2], "-o", color=C3, ms=5, lw=2.4, label="me3 (strongly cycle-dependent)")
axD.plot(Ts, Mk, "--s", color=CMK, ms=4, lw=1.8, label="lumped Mk (nearly flat)")
axD.set_xlabel("cell-cycle period (h)  →  slower cycling"); axD.set_ylabel("settled mean fraction")
axD.set_ylim(-0.02, 1.0); axD.grid(alpha=0.25); axD.legend(fontsize=8.3, loc="center right")
axD.set_title("D.  me3 depends strongly on cycle speed — the lumped mark does not", fontsize=10.6, fontweight="bold")
axD.text(0.03, 0.93, "the slow me2→me3 step makes me3 a DIVISION-TIMING readout:\nfast cycling → me3 diluted out; slow cycling → me3 accumulates",
         transform=axD.transAxes, fontsize=8, va="top", color="#7b241c", fontweight="bold")

out = os.path.join(os.path.dirname(__file__), "fig_v44_me1me2me3_chain")
fig.savefig(out + ".png", dpi=150, bbox_inches="tight")
fig.savefig(out + ".pdf", bbox_inches="tight")
print("wrote", out + ".png /.pdf")
print("steady (no div) me1/me2/me3 =", np.round(simulate(T=24, tmax=400, divide=False)[1][-1], 3))
