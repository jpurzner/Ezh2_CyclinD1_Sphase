"""Figure: v44.1 model wiring diagram — modules, cell-cycle flow, and key feedback loops.

CyclinD1 is the convergence node: the Hedgehog/MYCN mitogenic drive and the EZH2/PRC2 epigenetic
brake both act on it, and it hands off to the cell-cycle engine's restriction point.

v44.1 mechanism (this diagram reflects the baked, 28/28 read-write model):
  * Repression flows through a FORMAL PRC2 COMPLEX occupancy at Ccnd1, NOT the mark level directly.
    PRC2 = EZH2(1-EZH2i) x (accessory a0 + read-write a_rw*Mk) x (1 - nascent-Ccnd1 eviction).
  * H3K27me3 (Mk) is a READ-WRITE AMPLIFIER: the same PRC2 writes Mk and Mk recruits more PRC2.
  * Writer (EZH2/PRC2) is cycle-gated (E2F); eraser (Gli->Jmjd3/Kdm6b) is mitogen-gated ->
    high-Gli MB strips the mark (ChIP MB<GNP) while high EZH2 keeps occupancy (5.07 dose fold).
  * Feature C: commitment carryover across division (f_commit_carry) — committed daughters inherit
    phospho-Rb / CyclinE / low-p27 and re-enter G1 instead of resetting to deep G0.

Generates simulations/fig_v44_wiring.{png,pdf}. Pure schematic (no simulation).
Run:  ./venv/bin/python simulations/fig_v44_wiring.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

# ---------------- palette ----------------
C = dict(
    hh="#d4edda", hh_b="#28a745",          # Hedgehog / Gli (green)
    mycn="#fff3cd", mycn_b="#e67e22",       # MYCN (orange)
    cd="#fdfefe", cd_b="#117a65",           # CyclinD1 module (green-teal)
    cc="#d6eaf8", cc_b="#2980b9",           # cell-cycle engine (blue)
    rep="#d1f2eb", rep_b="#17a2b8",         # DNA replication (teal)
    ezh2="#e8daef", ezh2_b="#8e44ad",       # EZH2 / PRC2 / H3K27me3 (purple)
    growth="#eceff1", growth_b="#607d8b",   # growth (slate)
    chk="#fdebd0", chk_b="#e67e22",         # checkpoint (amber)
    node="#ffffff",
    act="#27ae60", inh="#c0392b", fb="#8e44ad", drug="#e74c3c", erase="#1f8a70",
)

fig, ax = plt.subplots(figsize=(18, 10.6))
ax.set_xlim(0, 17.8); ax.set_ylim(0, 10.6); ax.set_aspect("equal"); ax.axis("off")


def region(x, y, w, h, fc, ec, label, lx=None, ly=None):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.15",
                                fc=fc, ec=ec, lw=1.6, alpha=0.55, zorder=0))
    if label:
        ax.text(lx if lx is not None else x + 0.15, ly if ly is not None else y + h - 0.32,
                label, fontsize=10.5, fontweight="bold", color=ec, zorder=1, ha="left")


def node(x, y, text, w=1.55, h=0.6, fc=C["node"], ec="#2c3e50", fs=9, bold=False, italic=False):
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h, boxstyle="round,pad=0.03,rounding_size=0.08",
                                fc=fc, ec=ec, lw=1.3, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, zorder=4,
            fontweight="bold" if bold else "normal", fontstyle="italic" if italic else "normal")
    return (x, y, w, h)


def arrow(p1, p2, color=C["act"], style="->", lw=1.8, rad=0.0, ls="-", z=2, shrink=8):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle=style, color=color, lw=lw, ls=ls,
                                 mutation_scale=15, shrinkA=shrink, shrinkB=shrink,
                                 connectionstyle=f"arc3,rad={rad}", zorder=z))


def inhibit(p1, p2, color=C["inh"], lw=1.8, rad=0.0, z=2, shrink=8):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-[", color=color, lw=lw,
                                 mutation_scale=9, shrinkA=shrink, shrinkB=shrink,
                                 connectionstyle=f"arc3,rad={rad}", zorder=z))


ax.text(8.6, 10.32, "v44.1 cell-cycle model — module wiring", fontsize=15, fontweight="bold", ha="center")
ax.text(8.6, 9.82, "Heldt 2018 core + explicit replication + mitotic switch + growth-gated R-point + Skp2–p27 feedforward"
        "  ·  CyclinD1 = convergence node  ·  PRC2/H3K27me3 read-write brake  ·  Gli→Jmjd3 eraser",
        fontsize=8.8, ha="center", color="#555")

# ===================== module regions =====================
region(0.2, 6.7, 4.5, 2.7, C["hh"], C["hh_b"], "Hedgehog / MYCN  (mitogenic drive)")
region(0.2, 4.5, 4.5, 1.95, C["cd"], C["cd_b"], "CyclinD1 module  (convergence node)")
region(0.2, 0.4, 4.5, 3.9, C["growth"], C["growth_b"], "")
ax.text(0.35, 3.9, "Cell growth & size control", fontsize=10.5, fontweight="bold", color=C["growth_b"])
region(4.9, 3.0, 8.2, 6.2, C["cc"], C["cc_b"], "Cell-cycle engine", lx=8.25)
region(4.9, 0.4, 8.2, 2.3, C["rep"], C["rep_b"], "Explicit DNA replication  (fork speed kSyDna)")
region(13.05, 3.7, 4.55, 5.55, C["ezh2"], C["ezh2_b"], "EZH2 / PRC2 layer")
region(13.3, 1.5, 3.5, 2.3, C["chk"], C["chk_b"], "Mitotic switch + intra-S checkpoint")

# ===================== Hedgehog / MYCN  (with Gli→Ptch1 negative feedback) =====================
shh = node(0.85, 8.45, "SHH", w=0.95, fc="#fef9e7", fs=8.5)
smo = node(2.55, 8.45, "Smo → GliA", w=1.7, fs=8)
ptch = node(0.95, 7.25, "Ptch1\n(Gli target)", w=1.35, h=0.62, fs=7.2)
mycn = node(3.55, 7.15, "MYCN", w=1.2, h=0.6, fc=C["mycn"], ec=C["mycn_b"], fs=8.5, bold=True)
# SHH relieves Ptch1's inhibition of Smo
arrow(shh[:2], (smo[0]-0.75, smo[1]), C["act"])
# --- canonical negative feedback: GliA → Ptch1 ⊣ Smo  (intact in GNP; BROKEN in MB) ---
arrow((smo[0]-0.6, smo[1]-0.3), (ptch[0]+0.4, ptch[1]+0.28), C["hh_b"], rad=-0.30, lw=1.5)      # GliA induces Ptch1
inhibit((ptch[0]+0.0, ptch[1]+0.33), (smo[0]-0.67, smo[1]-0.25), C["inh"], rad=0.30, shrink=3)  # functional Ptch1 ⊣ Smo
ax.text(0.30, 6.86, "neg. feedback GliA→Ptch1⊣Smo  ·  Ptch1 f: GNP 1.0 → MB ~0.1 (broken → Gli↑)",
        fontsize=6.3, color=C["hh_b"], ha="left", fontweight="bold")
# drug HHi ⊣ Smo
hhi = node(2.55, 9.15, "HHi (vismo)", w=1.6, h=0.5, fc="#fdedeb", ec=C["drug"], fs=7.5, bold=True)
inhibit((hhi[0], hhi[1]-0.25), (smo[0], smo[1]+0.3), C["drug"], shrink=4)

# ===================== CyclinD1 module =====================
cdm = node(1.55, 5.4, "Ccnd1 mRNA", w=1.75, h=0.62, fc="#eafaf1", ec=C["cd_b"], fs=8)
cd = node(3.6, 5.4, "CyclinD1\n(Cd)", w=1.55, h=0.7, fc=C["cd"], ec=C["cd_b"], fs=9, bold=True)
arrow((cdm[0]+0.55, cdm[1]), (cd[0]-0.55, cd[1]), C["act"])
# drive in: Smo→GliA and MYCN drive Ccnd1 transcription
arrow((smo[0]-0.2, smo[1]-0.32), (cdm[0]+0.15, cdm[1]+0.35), C["act"], rad=0.18)   # GliA → Ccnd1 tx
arrow((mycn[0], mycn[1]-0.32), (cd[0]+0.1, cd[1]+0.38), C["mycn_b"], rad=-0.12)    # MYCN → Ccnd1 tx
ax.text(0.55, 4.62, "Cd ≈ (basal + Gli + MYCN drive) × R(PRC2 occupancy)", fontsize=6.5, color=C["cd_b"], ha="left", style="italic")

# ===================== growth =====================
mass = node(2.4, 3.2, "cell mass\n(grows, ÷2 at division)", w=2.6, h=0.7, fc="#ffffff", ec=C["growth_b"], fs=8)
ax.text(2.45, 2.35, "size gates:", fontsize=8.5, ha="center", color=C["growth_b"], fontstyle="italic")
node(1.5, 1.65, "M_commit\n(G0→G1)", w=1.5, h=0.6, fc="#ffffff", ec=C["growth_b"], fs=7.5)
msize = node(3.4, 1.65, "M_size\n(→ S-entry)", w=1.5, h=0.6, fc="#ffffff", ec=C["growth_b"], fs=7.5)

# ===================== cell-cycle engine (flow) =====================
ax.add_patch(FancyBboxPatch((5.6, 6.55), 3.0, 2.1, boxstyle="round,pad=0.03,rounding_size=0.1",
                            fc="#eaf2f8", ec=C["cc_b"], lw=1.4, ls="--", zorder=1))
ax.text(7.1, 8.45, "Restriction point", fontsize=9, fontweight="bold", color=C["cc_b"], ha="center", style="italic")
p27 = node(6.1, 7.85, "p27↑\n(=P21)", w=1.1, h=0.62, fc="#fbeee6", ec="#b9770e", fs=8, bold=True)
skp2 = node(6.1, 6.95, "Skp2↓", w=1.0, h=0.55, fs=8)
rb = node(7.5, 7.85, "Rb / pRb", w=1.3, fs=8.5)
e2f = node(7.9, 6.9, "E2F", w=1.0, fc="#fdfefe", ec="#117a65", fs=9, bold=True)
inhibit((p27[0]+0.45, p27[1]-0.3), (skp2[0], skp2[1]+0.28), C["inh"], shrink=3)
inhibit((skp2[0]+0.4, skp2[1]+0.1), (p27[0]+0.2, p27[1]-0.3), C["inh"], rad=0.3, shrink=3)
arrow((rb[0], rb[1]-0.3), (e2f[0]-0.2, e2f[1]+0.25), C["act"], rad=-0.1)
arrow((e2f[0]+0.2, e2f[1]+0.3), (skp2[0]+0.4, skp2[1]+0.3), C["act"], rad=-0.4, lw=1.4)  # E2F→Skp2
inhibit((p27[0]+0.5, p27[1]), (rb[0]-0.55, rb[1]), C["inh"], shrink=3)  # p27 ⊣ CDK2/Rb
# CyclinD1 module → Rb (via CDK4/6) — the hand-off from the convergence node to the engine
arrow((cd[0]+0.55, cd[1]+0.2), (rb[0]-0.6, rb[1]-0.2), C["cd_b"], rad=-0.18, lw=2.0)
ax.text(4.55, 6.95, "→ CDK4/6 → Rb", fontsize=6.8, color=C["cd_b"], ha="center", style="italic", zorder=5)
# --- INK4 (p16 + p18) COMPETITIVE brake + palbociclib (Vmax block) on the CyclinD1→Rb arm ---
p16 = node(5.05, 8.72, "p16 + p18\n(INK4)", w=1.5, h=0.55, fc="#f4ecf7", ec=C["ezh2_b"], fs=7.2, bold=True)
inhibit((p16[0]-0.1, p16[1]-0.3), (5.0, 7.84), C["inh"], shrink=3)
ax.text(4.2, 9.18, "INK4 CDK4/6 brake: p16 (GNP 0→MB) + p18 (GNP-expressed)", fontsize=6.2,
        color=C["ezh2_b"], ha="left", style="italic")
palbo = node(5.15, 6.78, "CDK4/6i\n(palbo)", w=1.35, h=0.6, fc="#fdedeb", ec=C["drug"], fs=7.2, bold=True)
inhibit((palbo[0]+0.05, palbo[1]+0.32), (5.1, 7.66), C["drug"], shrink=3)
ax.text(5.15, 6.36, "Vmax block (not surmountable by Cd)", fontsize=6.0, color=C["drug"], ha="center", style="italic")
# size gate on commitment
arrow((mass[0]+0.9, mass[1]+0.35), (5.7, 7.0), C["growth_b"], rad=-0.25, ls=":", lw=1.4)

# CDK2 cyclins
ce = node(9.8, 7.5, "CyclinE/A\n(CDK2)", w=1.5, h=0.7, fc="#d6eaf8", ec=C["cc_b"], fs=8.5, bold=True)
arrow((e2f[0]+0.4, e2f[1]), (ce[0]-0.7, ce[1]), C["act"])

# origins / replication (size-gated S entry)
ori = node(9.8, 5.0, "origins fire\n(Rc→aRc)", w=1.6, h=0.7, fc="#ffffff", ec=C["rep_b"], fs=8)
arrow((ce[0], ce[1]-0.35), (ori[0], ori[1]+0.35), C["act"])
arrow((msize[0]+0.5, msize[1]+0.2), (ori[0]-0.55, ori[1]-0.35), C["growth_b"], rad=0.2, ls=":", lw=1.4)  # M_size gate
dna = node(8.0, 1.5, "DNA replication\nDna: 0 → 1", w=2.0, h=0.7, fc="#ffffff", ec=C["rep_b"], fs=8.5, bold=True)
arrow((ori[0], ori[1]-0.35), (dna[0]+0.7, dna[1]+0.4), C["act"], rad=0.2)
hu = node(5.7, 1.5, "HU", w=0.9, fc="#fdedeb", ec=C["drug"], fs=8.5, bold=True)
inhibit((hu[0]+0.45, hu[1]), (dna[0]-1.0, dna[1]), C["drug"], shrink=3)
ax.text(6.5, 1.02, "↓ fork speed → S lengthens", fontsize=7.5, color=C["drug"], ha="center", style="italic")

# mitotic switch + checkpoint
mpf = node(14.9, 2.95, "CyclinB/CDK1\n(MPF)\nCdc25 / Wee1", w=2.4, h=1.0, fc="#ffffff", ec=C["chk_b"], fs=8, bold=True)
chk = node(14.4, 1.9, "CHK1\n(unfinished S)", w=1.7, h=0.55, fc="#fdebd0", ec=C["chk_b"], fs=7.5)
arrow((dna[0]+1.0, dna[1]+0.2), (mpf[0]-1.2, mpf[1]-0.2), C["act"], rad=-0.15)  # Dna done → MPF (G2)
inhibit((chk[0]+0.55, chk[1]+0.25), (mpf[0]-0.7, mpf[1]-0.45), C["inh"], shrink=3)  # CHK1 ⊣ MPF
arrow((dna[0]+1.0, dna[1]-0.1), (chk[0]-0.85, chk[1]), C["act"], rad=0.15, lw=1.3)  # forks → CHK1
# mitosis → division → back to G0/G1 (Feature C carryover)
arrow((mpf[0], mpf[1]+0.55), (8.0, 8.9), C["cc_b"], rad=0.35, lw=2.2)
ax.text(11.4, 9.15, "mitosis → division (÷2)", fontsize=8.5, color=C["cc_b"], ha="center", fontweight="bold")
arrow((8.0, 8.9), (p27[0], p27[1]+0.35), C["cc_b"], rad=0.15, lw=2.2)  # back toward G0/G1
ax.text(9.55, 8.62, "Feature C: commitment carryover (f_commit_carry≈0.57) —\ncommitted daughters keep pRb / CyclinE / low-p27 → re-enter G1, not deep G0",
        fontsize=6.3, color=C["cc_b"], ha="center", style="italic", zorder=6)

# ===================== EZH2 / PRC2 / H3K27me3 layer =====================
ezi = node(16.5, 9.0, "EZH2i (taz)", w=1.5, h=0.48, fc="#fdedeb", ec=C["drug"], fs=7.5, bold=True)
ezh2m = node(14.35, 8.38, "EZH2 mRNA\n(E2F target)", w=2.15, h=0.6, fc="#ffffff", ec=C["ezh2_b"], fs=7.8)
ezh2 = node(14.05, 7.28, "EZH2 protein", w=1.95, h=0.62, fc="#f4ecf7", ec=C["ezh2_b"], fs=8.2, bold=True)
prc2 = node(14.05, 5.95, "PRC2 complex\n(occupancy @ Ccnd1)", w=2.15, h=0.92, fc="#d7bde2", ec=C["ezh2_b"], fs=8, bold=True)
mark = node(16.45, 5.95, "H3K27me3\n(Mk)", w=1.6, h=0.92, fc="#e8daef", ec=C["ezh2_b"], fs=8.4, bold=True)
jmjd3 = node(16.45, 7.35, "Jmjd3 / Kdm6b\n(Gli1-induced)", w=1.6, h=0.72, fc="#d5f5e3", ec=C["erase"], fs=7.2, bold=True)

# writer chain: mRNA → protein → PRC2 (catalytic core)
arrow((ezh2m[0], ezh2m[1]-0.32), (ezh2[0]+0.1, ezh2[1]+0.34), C["act"])
arrow((ezh2[0], ezh2[1]-0.34), (prc2[0], prc2[1]+0.48), C["ezh2_b"], lw=1.8)
ax.text(13.0, 6.62, "catalytic\ncore", fontsize=6.2, color=C["ezh2_b"], ha="center", style="italic")
# EZH2 mRNA drivers: E2F (cycle) + mitogen-dose basal floor
ax.text(14.35, 8.82, "+ basal floor ∝ mitogen (kEZbas_Cd·Cd)", fontsize=5.9, color=C["ezh2_b"], ha="center", style="italic")
# EZH2i ⊣ PRC2 activity
inhibit((ezi[0]-0.3, ezi[1]-0.24), (prc2[0]+0.75, prc2[1]+0.42), C["drug"], rad=0.28, shrink=3)
# --- READ-WRITE LOOP: PRC2 writes Mk, Mk recruits PRC2 (amplifier) ---
arrow((prc2[0]+1.02, prc2[1]+0.14), (mark[0]-0.78, mark[1]+0.14), C["ezh2_b"], lw=1.9)   # write
ax.text(15.55, 6.34, "methylate (write)", fontsize=6.2, color=C["ezh2_b"], ha="center", style="italic")
arrow((mark[0]-0.78, mark[1]-0.16), (prc2[0]+1.02, prc2[1]-0.16), C["ezh2_b"], lw=1.9, ls=(0,(4,2)))  # read-write recruit
ax.text(15.55, 5.5, "read-write recruit\n(a_rw·Mk)", fontsize=6.2, color=C["ezh2_b"], ha="center", style="italic")
# eraser: Jmjd3 ⊣ Mk
inhibit((jmjd3[0], jmjd3[1]-0.36), (mark[0], mark[1]+0.48), C["erase"], shrink=3)
ax.text(17.25, 6.66, "erase", fontsize=6.4, color=C["erase"], ha="left", style="italic")
# Gli1 → Jmjd3 (long, mitogen-gated eraser arm)
arrow((smo[0]+0.35, smo[1]+0.28), (jmjd3[0]-0.6, jmjd3[1]+0.36), C["erase"], rad=-0.20, lw=1.5, ls=(0,(5,2)))
ax.text(8.7, 9.36, "Gli1 → Jmjd3/Kdm6b : mitogen-gated eraser (strips the Ccnd1 mark; high in MB → ChIP MB<GNP)",
        fontsize=7.2, color=C["erase"], ha="center", fontweight="bold")
# Mk turnover + replicative dilution
arrow((mark[0]+0.78, mark[1]+0.22), (mark[0]+0.78, mark[1]-0.22), C["rep_b"], rad=-1.6, lw=1.3, shrink=2)
ax.text(17.32, 5.95, "÷2 at\nearly S\n+ turnover", fontsize=5.8, color=C["rep_b"], ha="left", va="center", style="italic")
# nascent Ccnd1 RNA ⊣ PRC2 (g_prc2 eviction; double-negative)
inhibit((13.1, 4.95), (prc2[0]-0.5, prc2[1]-0.42), C["cd_b"], rad=0.22, shrink=3)
ax.text(13.0, 4.62, "nascent Ccnd1 RNA ⊣ PRC2 (g_prc2 eviction)", fontsize=6.0, color=C["cd_b"], ha="left", style="italic")

# ===================== central feedback: PRC2 occupancy ⊣ CyclinD1 =====================
inhibit((prc2[0]-1.02, prc2[1]+0.15), (cd[0]+0.2, cd[1]+0.42), C["fb"], rad=0.52, lw=2.4, shrink=4)
ax.text(8.6, 4.28, "PRC2 occupancy ⊣ Ccnd1 transcription  (Hill on PRC2, leaky floor f0_prc2;  H3K27me3 amplifies occupancy, not a standalone repressor)",
        fontsize=8.0, color=C["fb"], ha="center", fontweight="bold")

# ===================== phase bar (top-right) =====================
pm = [("G0", "#fbeee6", "#b9770e"), ("G1", "#eaf2f8", C["cc_b"]),
      ("S", "#d1f2eb", C["rep_b"]), ("G2/M", "#fdebd0", C["chk_b"])]
bx = 12.0
for i, (name, col, ec) in enumerate(pm):
    ax.add_patch(FancyBboxPatch((bx + i*1.22, 10.02), 1.18, 0.32, boxstyle="round,pad=0.01",
                                fc=col, ec=ec, lw=1.2, zorder=5))
    ax.text(bx + i*1.22 + 0.59, 10.18, name, ha="center", va="center", fontsize=8.5, fontweight="bold", zorder=6)
ax.text(bx - 0.15, 10.18, "phases:", ha="right", va="center", fontsize=8.5, fontstyle="italic", color="#555")
# legend
leg = [Line2D([0],[0], color=C["act"], lw=2, label="activation"),
       Line2D([0],[0], color=C["inh"], lw=2, label="inhibition (⊣)"),
       Line2D([0],[0], color=C["fb"], lw=2.4, label="PRC2 occupancy ⊣ CyclinD1 (epigenetic brake)"),
       Line2D([0],[0], color=C["ezh2_b"], lw=1.9, ls=(0,(4,2)), label="H3K27me3 read-write amplifier loop"),
       Line2D([0],[0], color=C["erase"], lw=1.6, ls=(0,(5,2)), label="Gli→Jmjd3 mark eraser (mitogen-gated)"),
       Line2D([0],[0], color=C["cd_b"], lw=2.0, label="CyclinD1 hand-off / PRC2 eviction"),
       Line2D([0],[0], color=C["growth_b"], lw=1.6, ls=":", label="size gate"),
       Line2D([0],[0], color=C["drug"], lw=2, label="drug / input")]
ax.legend(handles=leg, loc="lower right", fontsize=7.8, framealpha=0.92, ncol=2,
          columnspacing=1.1, handlelength=1.8).set_zorder(10)

plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), "fig_v44_wiring")
fig.savefig(out + ".png", dpi=160, bbox_inches="tight")
fig.savefig(out + ".pdf", bbox_inches="tight")
print("wrote", out + ".png /.pdf")
