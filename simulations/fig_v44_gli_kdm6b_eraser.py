"""Figure: DETAIL wiring + equations of the Gli1 -> Kdm6b/Jmjd3 -> H3K27me3 eraser arm (v44.1).

Focused schematic of the mitogen-gated eraser that opposes the EZH2/PRC2 writer at the Ccnd1
H3K27me3 domain. Shows (top) the wiring from Hedgehog/Smo through Gli1 (incl. the Gli1_epi
slow-memory capacitor) to the demethylase and the mark, and (bottom) the exact model equations.

IMPORTANT modeling detail made explicit here: Kdm6b/Jmjd3 is NOT an explicit species. Its level
is taken quasi-static and proportional to Gli1, so the demethylation flux is the lumped term
k_jmjd3_gli * Gli1 * Mk. (Shi 2014, Nat Commun 5:5425 — Gli activators recruit KDM6B/JMJD3.)

Generates simulations/fig_v44_gli_kdm6b_eraser.{png,pdf}. Pure schematic (no simulation).
Run:  ./venv/bin/python simulations/fig_v44_gli_kdm6b_eraser.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

C = dict(
    hh="#d4edda", hh_b="#28a745",
    gli="#eafaf1", gli_b="#1e8449",
    epi="#fef9e7", epi_b="#b7950b",
    erase="#d5f5e3", erase_b="#148f77",
    ezh2="#e8daef", ezh2_b="#8e44ad",
    mark="#e8daef", mark_b="#8e44ad",
    cd="#fdfefe", cd_b="#117a65",
    node="#ffffff",
    act="#27ae60", inh="#c0392b", drug="#e74c3c", rep_b="#17a2b8",
    eqbg="#f7f9fb", eqbg2="#f4ecf7", eqbg3="#eafaf1", eqbg4="#fef9e7",
)

fig, ax = plt.subplots(figsize=(15.5, 11.6))
ax.set_xlim(0, 15.5); ax.set_ylim(0, 11.6); ax.axis("off")


def node(x, y, text, w=1.7, h=0.66, fc=C["node"], ec="#2c3e50", fs=9.2, bold=False, ls="-", italic=False):
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h, boxstyle="round,pad=0.03,rounding_size=0.09",
                                fc=fc, ec=ec, lw=1.4, ls=ls, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, zorder=4,
            fontweight="bold" if bold else "normal", fontstyle="italic" if italic else "normal")
    return (x, y, w, h)


def arrow(p1, p2, color=C["act"], style="->", lw=1.9, rad=0.0, ls="-", z=2, shrink=8):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle=style, color=color, lw=lw, ls=ls,
                                 mutation_scale=15, shrinkA=shrink, shrinkB=shrink,
                                 connectionstyle=f"arc3,rad={rad}", zorder=z))


def inhibit(p1, p2, color=C["inh"], lw=1.9, rad=0.0, z=2, shrink=8):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-[", color=color, lw=lw,
                                 mutation_scale=10, shrinkA=shrink, shrinkB=shrink,
                                 connectionstyle=f"arc3,rad={rad}", zorder=z))


ax.text(7.75, 11.28, "Gli1 → Kdm6b/Jmjd3 eraser of H3K27me3 at $\\it{Ccnd1}$ — detailed wiring + equations (v44.1)",
        fontsize=14.5, fontweight="bold", ha="center")
ax.text(7.75, 10.92, "mitogen-gated ERASER (Gli1 → KDM6B/JMJD3 demethylase)  vs  cycle-gated WRITER (EZH2/PRC2)  —  the per-division race that sets the mark",
        fontsize=9.2, ha="center", color="#555", style="italic")

# ================= WIRING (top) =================
# --- Hedgehog / Smo input ---
smo = node(1.35, 9.75, "Smo·active\n(Hh input;\nvismo target)", w=1.9, h=1.0, fc=C["hh"], ec=C["hh_b"], fs=8.2, bold=True)

# --- Gli activator/repressor switch ---
grep_ = node(3.35, 10.15, "Gli_rep", w=1.5, h=0.56, fc="#fdecea", ec=C["inh"], fs=8.5)
gact = node(3.35, 9.2, "Gli_act", w=1.5, h=0.56, fc=C["gli"], ec=C["gli_b"], fs=8.5, bold=True)
ax.annotate("", xy=(gact[0]-0.0, gact[1]+0.30), xytext=(grep_[0]-0.0, grep_[1]-0.30),
            arrowprops=dict(arrowstyle="<|-|>", color="#7f8c8d", lw=1.5), zorder=2)
ax.text(4.15, 9.68, "Smo-gated\nswitch", fontsize=6.6, color="#7f8c8d", ha="left", style="italic")
arrow((smo[0]+0.55, smo[1]+0.25), (grep_[0]-0.78, grep_[1]-0.05), C["hh_b"], rad=-0.15, lw=1.5)

# --- Gli1_epi slow-memory capacitor ---
epi = node(1.35, 7.7, "Gli1_epi\n(slow chromatin\nmemory capacitor)", w=2.1, h=0.98, fc=C["epi"], ec=C["epi_b"], fs=7.6, bold=True, ls=(0,(4,2)))
arrow((smo[0], smo[1]-0.52), (epi[0], epi[1]+0.5), C["epi_b"], lw=1.6)                       # Smo charges epi
ax.text(1.62, 9.0, "charge ∝ Smo\n(k_epi_on)", fontsize=6.6, color=C["epi_b"], ha="left", style="italic")
# self-discharge loop
ax.add_patch(FancyArrowPatch((epi[0]-1.02, epi[1]-0.12), (epi[0]-1.02, epi[1]+0.12), arrowstyle="->",
             color=C["epi_b"], lw=1.3, mutation_scale=11, connectionstyle="arc3,rad=-1.8", shrinkA=2, shrinkB=2, zorder=2))
ax.text(0.12, 7.7, "slow\ndischarge\n(k_epi_off)", fontsize=6.3, color=C["epi_b"], ha="left", va="center", style="italic")

# --- Gli1 mRNA -> Gli1 protein ---
g1m = node(5.55, 8.95, "Gli1 mRNA", w=1.7, h=0.6, fc=C["gli"], ec=C["gli_b"], fs=8.4)
g1 = node(7.6, 8.95, "Gli1\n(protein)", w=1.6, h=0.8, fc=C["gli"], ec=C["gli_b"], fs=9.6, bold=True)
arrow((gact[0]+0.78, gact[1]+0.05), (g1m[0]-0.78, g1m[1]-0.18), C["act"], rad=0.1)            # Gli_act -> tx
inhibit((grep_[0]+0.78, grep_[1]-0.1), (g1m[0]-0.7, g1m[1]+0.2), C["inh"], rad=-0.15, shrink=4)  # Gli_rep -| tx
arrow((epi[0]+1.08, epi[1]+0.2), (g1m[0]-0.72, g1m[1]-0.05), C["epi_b"], rad=-0.28, lw=1.6)   # epi -> tx (autoreg)
ax.text(4.35, 8.05, "autoregulation\n(k_Gli1_auto·Gli1_epi)", fontsize=6.6, color=C["epi_b"], ha="center", style="italic")
arrow((g1m[0]+0.78, g1m[1]), (g1[0]-0.72, g1[1]), C["act"])                                    # translation

# --- Kdm6b/Jmjd3 (implicit) ---
kdm = node(9.9, 8.95, "KDM6B / JMJD3\n(H3K27 demethylase)\nimplicit: [KDM6B] ∝ Gli1", w=2.7, h=1.0,
           fc=C["erase"], ec=C["erase_b"], fs=7.8, bold=True, ls=(0,(4,2)))
arrow((g1[0]+0.72, g1[1]), (kdm[0]-1.38, kdm[1]), C["erase_b"], lw=2.0)                        # Gli1 induces demethylase
ax.text(8.75, 9.42, "induces\n(recruit)", fontsize=6.8, color=C["erase_b"], ha="center", style="italic")

# --- H3K27me3 (Mk) mark ---
mk = node(12.75, 8.35, "H3K27me3\nat $\\it{Ccnd1}$  (Mk)\n[0,1] occupancy", w=2.4, h=1.05, fc=C["mark"], ec=C["mark_b"], fs=8.6, bold=True)
# ERASER: Kdm6b -| Mk
inhibit((kdm[0]+0.9, kdm[1]-0.28), (mk[0]-1.05, mk[1]+0.32), C["erase_b"], rad=-0.12, lw=2.2, shrink=5)
ax.text(11.45, 9.12, "ERASE\nk_jmjd3_gli·Gli1·Mk", fontsize=7.4, color=C["erase_b"], ha="center", fontweight="bold")

# --- WRITER arm: EZH2 -> PRC2 -> Mk ---
ezh2 = node(10.4, 7.55, "EZH2\n(E2F target;\ncycle-gated)", w=1.7, h=0.9, fc=C["ezh2"], ec=C["ezh2_b"], fs=7.8, bold=True)
prc2 = node(11.35, 6.55, "PRC2 complex\n(occupancy)", w=2.0, h=0.8, fc="#d7bde2", ec=C["ezh2_b"], fs=8.2, bold=True)
arrow((ezh2[0]+0.35, ezh2[1]-0.46), (prc2[0]-0.25, prc2[1]+0.42), C["ezh2_b"], lw=1.8)
arrow((prc2[0]+0.35, prc2[1]+0.44), (mk[0]-0.55, mk[1]-0.55), C["ezh2_b"], rad=-0.18, lw=2.0)   # PRC2 writes Mk
ax.text(12.6, 7.5, "WRITE\n(methylate)\nPRC2·(1−Mk)", fontsize=7.1, color=C["ezh2_b"], ha="center", fontweight="bold")
# read-write: Mk -> PRC2
arrow((mk[0]-0.78, mk[1]-0.46), (prc2[0]+0.6, prc2[1]+0.04), C["ezh2_b"], rad=-0.3, lw=1.6, ls=(0,(4,2)))
ax.text(11.2, 7.42, "read-write\n(a_rw·Mk)", fontsize=6.8, color=C["ezh2_b"], ha="center", style="italic")
ezi = node(9.15, 6.72, "EZH2i (taz)", w=1.5, h=0.5, fc="#fdedeb", ec=C["drug"], fs=7.5, bold=True)
inhibit((ezi[0]+0.35, ezi[1]-0.02), (prc2[0]-0.85, prc2[1]-0.02), C["drug"], rad=0.12, shrink=4)

# --- Mk turnover + replicative dilution + downstream ---
ax.add_patch(FancyArrowPatch((mk[0]+1.02, mk[1]+0.28), (mk[0]+1.02, mk[1]-0.28), arrowstyle="->",
             color=C["rep_b"], lw=1.4, mutation_scale=12, connectionstyle="arc3,rad=-1.7", shrinkA=2, shrinkB=2, zorder=2))
ax.text(14.05, 8.35, "÷2 at\nearly S\n+ turnover\n(del_mk)", fontsize=6.4, color=C["rep_b"], ha="left", va="center", style="italic")
cd = node(14.15, 6.5, "CyclinD1\ntranscription", w=1.9, h=0.72, fc=C["cd"], ec=C["cd_b"], fs=8.1, bold=True)
inhibit((mk[0]+0.55, mk[1]-0.52), (cd[0]-0.15, cd[1]+0.4), C["ezh2_b"], rad=-0.15, lw=1.8, shrink=4)
ax.text(13.05, 7.2, "PRC2 occ.\n⊣ Ccnd1", fontsize=6.5, color=C["ezh2_b"], ha="center", style="italic")

# race banner
ax.add_patch(FancyBboxPatch((4.6, 6.5), 3.95, 1.2, boxstyle="round,pad=0.05,rounding_size=0.12",
                            fc="#fbfcfc", ec="#7f8c8d", lw=1.3, ls=(0,(3,2)), zorder=1))
ax.text(6.57, 7.38, "the per-division RACE", fontsize=9.0, fontweight="bold", ha="center", color="#34495e")
ax.text(6.57, 7.02, "WRITER  EZH2/PRC2  (↑ in S, cycle-gated)", fontsize=7.4, ha="center", color=C["ezh2_b"])
ax.text(6.57, 6.73, "ERASER  Gli1/KDM6B  (↑ with mitogen, MB≈6.9×)", fontsize=7.4, ha="center", color=C["erase_b"])

# divider
ax.plot([0.3, 15.2], [5.95, 5.95], color="#d5dbdb", lw=1.0, zorder=0)

# ================= EQUATIONS (bottom) =================
def eqbox(x, y, w, h, fc, ec, title):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.1",
                                fc=fc, ec=ec, lw=1.4, alpha=0.9, zorder=1))
    ax.text(x + 0.18, y + h - 0.28, title, fontsize=9.4, fontweight="bold", color=ec, zorder=2)

def eqline(x, y, s, fs=9.0, color="#1b2631", bold=False):
    ax.text(x, y, s, fontsize=fs, family="monospace", color=color, zorder=2, va="center",
            fontweight="bold" if bold else "normal")

# --- Panel A: Gli1 axis ---
eqbox(0.3, 3.35, 9.35, 2.15, C["eqbg3"], C["gli_b"], "① Gli1 axis  (mitogen input, with slow autoregulatory memory)")
eqline(0.55, 4.86, "Gli_rep  ⇌  Gli_act        (Smo-gated activator/repressor switch)", 8.6, "#5d6d7e")
eqline(0.55, 4.50, "d[Gli1_mRNA]/dt = Vmax·Gli_act²/(K_GA²+Gli_act²)·(1 − Gli_rep²/(K_GR²+Gli_rep²))·γ_Ptch", 8.4)
eqline(0.55, 4.20, "                  + k_Gli1_auto·Gli1_epi  −  k_mdeg·Gli1_mRNA", 8.4, C["epi_b"])
eqline(0.55, 3.84, "d[Gli1_epi]/dt  = k_epi_on·(1−Ptch)·Smo⁴/(K_auto⁴+Smo⁴)·(1−Gli1_epi) − k_epi_off·Gli1_epi", 8.4, C["epi_b"])
eqline(0.55, 3.56, "d[Gli1]/dt      = k_Gli1_tl·Gli1_mRNA − k_Gli1_deg·Gli1", 8.4)

# --- Panel B: eraser ---
eqbox(9.85, 3.35, 5.35, 2.15, C["eqbg"], C["erase_b"], "② Eraser  (KDM6B/JMJD3, implicit)")
eqline(10.08, 4.82, "KDM6B is NOT a state variable.", 8.2, "#5d6d7e")
eqline(10.08, 4.54, "Quasi-static:  [KDM6B] ∝ Gli1,", 8.2, "#5d6d7e")
eqline(10.08, 4.28, "folded into k_jmjd3_gli.", 8.2, "#5d6d7e")
eqline(10.08, 3.88, "erase flux =", 8.6, C["erase_b"], bold=True)
eqline(10.08, 3.60, "  k_jmjd3_gli · Gli1 · Mk", 9.0, C["erase_b"], bold=True)

# --- Panel C: Mk balance ---
eqbox(0.3, 0.95, 9.35, 2.2, C["eqbg2"], C["mark_b"], "③ H3K27me3 (Mk) balance at the $Ccnd1$ domain")
eqline(0.55, 2.55, "d[Mk]/dt =  PRC2·(1 − Mk)   −   del_mk·Mk   −   k_jmjd3_gli·Gli1·Mk", 8.9, "#1b2631", bold=True)
eqline(0.55, 2.28, "            └ WRITE ────────┘  └ turnover ┘   └ ERASE ─────────┘", 8.2, "#7f8c8d")
eqline(0.55, 1.98, "            (EZH2/PRC2, cycle)               (Gli1/KDM6B, mitogen)", 7.8, "#7f8c8d")
eqline(0.55, 1.64, "discrete:  at Dna > 0.05  →  Mk ← ½·Mk        (replicative dilution, once per S)", 8.4, C["rep_b"])
eqline(0.55, 1.30, "PRC2 = EZH2·(1−EZH2i)·(a0 + a_rw·Mk)·(1 − g·Cd_mRNA^p/(K_tx^p+Cd_mRNA^p))", 8.2, C["ezh2_b"])

# --- Panel D: params + net ---
eqbox(9.85, 0.95, 5.35, 2.2, C["eqbg4"], C["epi_b"], "④ Baked parameters (v44.1) + net effect")
eqline(10.08, 2.58, "k_jmjd3_gli 0.0478   del_mk 0.0015", 7.9)
eqline(10.08, 2.32, "a0 3.78e-4  a_rw 4.30e-3 (a_rw/a0 11.4)", 7.9, C["ezh2_b"])
eqline(10.08, 2.06, "g_prc2 0.0424  K_prc2 3.5e-3  n 3.39", 7.9, C["ezh2_b"])
eqline(10.08, 1.80, "f0_prc2 0.0505  K_tx 5.12  p_tx 3.27", 7.9, C["ezh2_b"])
eqline(10.08, 1.52, "Gli: Vmax 1.27  auto 0.012  on/off 0.02/8e-4", 7.9, C["gli_b"])
ax.text(10.05, 1.16, "MB Gli1 ≈ 6.9× GNP → erase ≈ 6.9× stronger →\nChIP $Ccnd1$ mark MB/GNP ≈ 0.55 (MB < GNP)",
        fontsize=7.4, color=C["erase_b"], ha="left", va="center", style="italic", fontweight="bold")

# bottom takeaway
ax.text(7.75, 0.5, "Net over divisions:  Mk accumulates when  WRITE  >  (turnover + ERASE + ½-dilution)  each cycle; "
        "MB's high Gli1 tilts the race toward erasure despite higher EZH2.",
        fontsize=8.4, ha="center", color="#34495e", style="italic")

# legend
leg = [Line2D([0],[0], color=C["act"], lw=2, label="activation / induction"),
       Line2D([0],[0], color=C["inh"], lw=2, label="repression (⊣)"),
       Line2D([0],[0], color=C["erase_b"], lw=2.2, label="Gli1→KDM6B eraser (demethylate)"),
       Line2D([0],[0], color=C["ezh2_b"], lw=2, label="EZH2/PRC2 writer"),
       Line2D([0],[0], color=C["ezh2_b"], lw=1.6, ls=(0,(4,2)), label="H3K27me3 read-write recruit"),
       Line2D([0],[0], color=C["rep_b"], lw=1.5, label="replicative dilution / turnover"),
       Line2D([0],[0], color="#7f8c8d", lw=1.3, ls=(0,(4,2)), label="implicit / quasi-static node")]
ax.legend(handles=leg, loc="upper left", bbox_to_anchor=(0.002, 0.585), fontsize=7.8, framealpha=0.92, ncol=1).set_zorder(10)

plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), "fig_v44_gli_kdm6b_eraser")
fig.savefig(out + ".png", dpi=160, bbox_inches="tight")
fig.savefig(out + ".pdf", bbox_inches="tight")
print("wrote", out + ".png /.pdf")
