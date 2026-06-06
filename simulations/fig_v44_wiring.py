"""Figure: v44 model wiring diagram — modules, cell-cycle flow, and key feedback loops.

Generates simulations/fig_v44_wiring.{png,pdf}. Pure schematic (no simulation).
Run:  ./venv/bin/python simulations/fig_v44_wiring.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D

# ---------------- palette ----------------
C = dict(
    hh="#d4edda", hh_b="#28a745",          # Hedgehog (green)
    mycn="#fff3cd", mycn_b="#e67e22",       # MYCN (orange)
    cc="#d6eaf8", cc_b="#2980b9",           # cell-cycle engine (blue)
    rep="#d1f2eb", rep_b="#17a2b8",         # DNA replication (teal)
    ezh2="#e8daef", ezh2_b="#8e44ad",       # EZH2 (purple)
    growth="#eceff1", growth_b="#607d8b",   # growth (slate)
    chk="#fdebd0", chk_b="#e67e22",         # checkpoint (amber)
    node="#ffffff",
    act="#27ae60", inh="#c0392b", fb="#8e44ad", drug="#e74c3c",
)

fig, ax = plt.subplots(figsize=(17, 10.5))
ax.set_xlim(0, 17); ax.set_ylim(0, 10.5); ax.set_aspect("equal"); ax.axis("off")


def region(x, y, w, h, fc, ec, label, lx=None, ly=None):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.15",
                                fc=fc, ec=ec, lw=1.6, alpha=0.55, zorder=0))
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


ax.text(8.5, 10.2, "v44 cell-cycle model — module wiring", fontsize=15, fontweight="bold", ha="center")
ax.text(8.5, 9.82, "Heldt 2018 core + explicit replication + mitotic switch + growth-gated R-point + Skp2–p27 feedforward + EZH2"
        "  ·  Gli→Ptch1 feedback  ·  p16 CDK4/6 brake",
        fontsize=9, ha="center", color="#555")

# ===================== module regions =====================
region(0.2, 5.85, 4.5, 3.35, C["hh"], C["hh_b"], "Hedgehog / MYCN  →  CyclinD1")
region(0.2, 0.4, 4.5, 5.3, C["growth"], C["growth_b"], "")
ax.text(0.35, 5.4, "Cell growth & size control", fontsize=10.5, fontweight="bold", color=C["growth_b"])
region(4.9, 3.0, 8.2, 6.2, C["cc"], C["cc_b"], "Cell-cycle engine", lx=8.25)
region(4.9, 0.4, 8.2, 2.3, C["rep"], C["rep_b"], "Explicit DNA replication  (fork speed kSyDna)")
region(13.3, 4.4, 3.5, 4.8, C["ezh2"], C["ezh2_b"], "EZH2 layer")
region(13.3, 0.4, 3.5, 3.6, C["chk"], C["chk_b"], "Mitotic switch + intra-S checkpoint")

# ===================== Hedgehog / MYCN  (with Gli→Ptch1 negative feedback) =====================
shh = node(0.85, 8.55, "SHH", w=0.95, fc="#fef9e7", fs=8.5)
smo = node(2.45, 8.55, "Smo → GliA", w=1.6, fs=8)
ptch = node(0.95, 7.05, "Ptch1\n(Gli target)", w=1.35, h=0.62, fs=7.2)
mycn = node(2.5, 6.6, "MYCN", w=1.2, fc=C["mycn"], ec=C["mycn_b"], fs=8.5, bold=True)
cd = node(4.0, 7.65, "CyclinD1\n(Cd)", w=1.4, h=0.7, fc="#fdfefe", ec="#117a65", fs=9, bold=True)
# SHH relieves Ptch1's inhibition of Smo; Smo→GliA drives CyclinD1 (& Gli1); MYCN also drives Cd
arrow(shh[:2], (smo[0]-0.7, smo[1]), C["act"])
arrow((smo[0]+0.3, smo[1]-0.3), (cd[0]-0.25, cd[1]+0.28), C["act"], rad=-0.12)
arrow(mycn[:2], (cd[0]-0.25, cd[1]-0.28), C["act"], rad=0.12)
# --- canonical negative feedback: GliA → Ptch1 ⊣ Smo  (intact in GNP; BROKEN in MB) ---
arrow((smo[0]-0.55, smo[1]-0.3), (ptch[0]+0.4, ptch[1]+0.28), C["hh_b"], rad=-0.30, lw=1.5)      # GliA induces Ptch1
inhibit((ptch[0]+0.0, ptch[1]+0.33), (smo[0]-0.62, smo[1]-0.25), C["inh"], rad=0.30, shrink=3)   # functional Ptch1 ⊣ Smo
ax.text(0.30, 6.32, "neg. feedback  GliA→Ptch1⊣Smo", fontsize=7.0, color=C["hh_b"], ha="left", fontweight="bold")
ax.text(0.30, 6.06, "functional Ptch1 f: GNP 1.0 → MB ~0.1 (broken → Gli↑)", fontsize=6.3, color="#555", ha="left", style="italic")
# drug GDC0449 ⊣ Smo
gdc = node(2.45, 9.6, "GDC0449", w=1.3, fc="#fdedeb", ec=C["drug"], fs=8, bold=True)
inhibit((gdc[0], gdc[1]-0.3), (smo[0], smo[1]+0.3), C["drug"], shrink=4)

# ===================== growth =====================
mass = node(2.4, 4.2, "cell mass\n(grows, ÷2 at division)", w=2.6, h=0.7, fc="#ffffff", ec=C["growth_b"], fs=8)
ax.text(2.45, 3.2, "size gates:", fontsize=8.5, ha="center", color=C["growth_b"], fontstyle="italic")
node(1.5, 2.5, "M_commit\n(G0→G1)", w=1.5, h=0.62, fc="#ffffff", ec=C["growth_b"], fs=7.5)
node(3.4, 2.5, "M_size\n(→ S-entry)", w=1.5, h=0.62, fc="#ffffff", ec=C["growth_b"], fs=7.5)

# ===================== cell-cycle engine (flow) =====================
# R-point box: Skp2-p27 feedforward
rp_x, rp_y = 6.6, 7.6
ax.add_patch(FancyBboxPatch((5.6, 6.55), 3.0, 2.1, boxstyle="round,pad=0.03,rounding_size=0.1",
                            fc="#eaf2f8", ec=C["cc_b"], lw=1.4, ls="--", zorder=1))
ax.text(7.1, 8.45, "Restriction point", fontsize=9, fontweight="bold", color=C["cc_b"], ha="center", style="italic")
p27 = node(6.1, 7.85, "p27↑\n(=P21)", w=1.1, h=0.62, fc="#fbeee6", ec="#b9770e", fs=8, bold=True)
skp2 = node(6.1, 6.95, "Skp2↓", w=1.0, h=0.55, fs=8)
rb = node(7.5, 7.85, "Rb / pRb", w=1.3, fs=8.5)
e2f = node(7.9, 6.9, "E2F", w=1.0, fc="#fdfefe", ec="#117a65", fs=9, bold=True)
# feedforward arrows
inhibit((p27[0]+0.45, p27[1]-0.3), (skp2[0], skp2[1]+0.28), C["inh"], shrink=3)
inhibit((skp2[0]+0.4, skp2[1]+0.1), (p27[0]+0.2, p27[1]-0.3), C["inh"], rad=0.3, shrink=3)
arrow((rb[0], rb[1]-0.3), (e2f[0]-0.2, e2f[1]+0.25), C["act"], rad=-0.1)
arrow((e2f[0]+0.2, e2f[1]+0.3), (skp2[0]+0.4, skp2[1]+0.3), C["act"], rad=-0.4, lw=1.4)  # E2F→Skp2
inhibit((p27[0]+0.5, p27[1]), (rb[0]-0.55, rb[1]), C["inh"], shrink=3)  # p27 ⊣ CDK2/Rb
arrow((cd[0]+0.5, cd[1]), (rb[0]-0.6, rb[1]+0.15), C["act"], rad=-0.15)  # CyclinD → Rb (via CDK4/6)
ax.text(5.2, 7.42, "CDK4/6", fontsize=6.8, color="#117a65", ha="center", style="italic", zorder=5)
# --- p16 (Cdkn2a) COMPETITIVE brake + palbociclib (Vmax block) on the CyclinD1→Rb arm ---
p16 = node(5.05, 8.72, "p16 (Cdkn2a)", w=1.5, h=0.5, fc="#f4ecf7", ec=C["ezh2_b"], fs=7.3, bold=True)
inhibit((p16[0]-0.1, p16[1]-0.27), (5.0, 7.84), C["inh"], shrink=3)                 # competitive: ↑Cd→Rb half-max
ax.text(4.25, 9.12, "EZH2-silenced (GNP 0) · high in MB", fontsize=6.3, color=C["ezh2_b"],
        ha="left", style="italic")
palbo = node(5.15, 6.78, "CDK4/6i\n(palbo)", w=1.35, h=0.6, fc="#fdedeb", ec=C["drug"], fs=7.2, bold=True)
inhibit((palbo[0]+0.05, palbo[1]+0.32), (5.1, 7.66), C["drug"], shrink=3)            # Vmax block (kPhRbCd=0)
ax.text(5.15, 6.36, "Vmax block (not surmountable by Cd)", fontsize=6.0, color=C["drug"], ha="center", style="italic")
# size gate on commitment
arrow((mass[0]+0.9, mass[1]+0.45), (5.7, 7.2), C["growth_b"], rad=-0.25, ls=":", lw=1.4)

# CDK2 cyclins
ce = node(9.8, 7.5, "CyclinE/A\n(CDK2)", w=1.5, h=0.7, fc="#d6eaf8", ec=C["cc_b"], fs=8.5, bold=True)
arrow((e2f[0]+0.4, e2f[1]), (ce[0]-0.7, ce[1]), C["act"])

# origins / replication (size-gated S entry)
ori = node(9.8, 5.0, "origins fire\n(Rc→aRc)", w=1.6, h=0.7, fc="#ffffff", ec=C["rep_b"], fs=8)
arrow((ce[0], ce[1]-0.35), (ori[0], ori[1]+0.35), C["act"])
arrow((3.4, 2.85), (ori[0]-0.5, ori[1]-0.35), C["growth_b"], rad=0.2, ls=":", lw=1.4)  # M_size gate
dna = node(8.0, 1.5, "DNA replication\nDna: 0 → 1", w=2.0, h=0.7, fc="#ffffff", ec=C["rep_b"], fs=8.5, bold=True)
arrow((ori[0], ori[1]-0.35), (dna[0]+0.7, dna[1]+0.4), C["act"], rad=0.2)
hu = node(5.7, 1.5, "HU", w=0.9, fc="#fdedeb", ec=C["drug"], fs=8.5, bold=True)
inhibit((hu[0]+0.45, hu[1]), (dna[0]-1.0, dna[1]), C["drug"], shrink=3)
ax.text(6.5, 1.05, "↓ fork speed → S lengthens", fontsize=7.5, color=C["drug"], ha="center", style="italic")

# mitotic switch + checkpoint
mpf = node(14.9, 2.6, "CyclinB/CDK1\n(MPF)\nCdc25 / Wee1", w=2.4, h=1.0, fc="#ffffff", ec=C["chk_b"], fs=8, bold=True)
chk = node(14.4, 1.1, "CHK1\n(unfinished S)", w=1.7, h=0.62, fc="#fdebd0", ec=C["chk_b"], fs=7.5)
arrow((dna[0]+1.0, dna[1]+0.2), (mpf[0]-1.2, mpf[1]-0.2), C["act"], rad=-0.15)  # Dna done → MPF (G2)
inhibit((chk[0]+0.55, chk[1]+0.25), (mpf[0]-0.7, mpf[1]-0.45), C["inh"], shrink=3)  # CHK1 ⊣ MPF
arrow((dna[0]+1.0, dna[1]-0.1), (chk[0]-0.85, chk[1]), C["act"], rad=0.15, lw=1.3)  # forks → CHK1
# mitosis → division → back to G0
arrow((mpf[0], mpf[1]+0.55), (8.0, 8.9), C["cc_b"], rad=0.35, lw=2.2)
ax.text(11.4, 9.1, "mitosis → division (÷2)", fontsize=8.5, color=C["cc_b"], ha="center", fontweight="bold")
arrow((8.0, 8.9), (p27[0], p27[1]+0.35), C["cc_b"], rad=0.15, lw=2.2)  # back to G0 (p27 high)

# ===================== EZH2 layer =====================
ezh2 = node(15.0, 7.0, "EZH2 protein\n(integrates S)", w=2.4, h=0.8, fc="#f4ecf7", ec=C["ezh2_b"], fs=8.5, bold=True)
ezh2m = node(15.0, 8.2, "EZH2 mRNA\n(E2F target)", w=2.4, h=0.7, fc="#ffffff", ec=C["ezh2_b"], fs=8)
arrow((ezh2m[0], ezh2m[1]-0.35), (ezh2[0], ezh2[1]+0.4), C["act"])
# EZH2 is an E2F target: E2F (released by pRb phosphorylation) drives EZH2 transcription;
# CyclinE/A (CDK2) only gate the S-window in the rate law E2f*(CycE+CycA).
arrow((e2f[0]+0.4, e2f[1]+0.2), (ezh2m[0]-1.25, ezh2m[1]-0.05), C["ezh2_b"], rad=-0.32, lw=1.6, ls="-")
ax.text(11.6, 8.7, "E2F target", fontsize=7.5, color=C["ezh2_b"], ha="center", style="italic")
ax.text(13.05, 7.95, "(S-window\ngated by CDK2)", fontsize=6.5, color="#888", ha="center", style="italic")
ezi = node(15.0, 9.0, "EZH2i", w=1.1, fc="#fdedeb", ec=C["drug"], fs=8, bold=True)
# central feedback: EZH2 ⊣ CyclinD1
inhibit((ezh2[0]-1.2, ezh2[1]-0.1), (cd[0]+0.1, cd[1]+0.4), C["fb"], rad=0.32, lw=2.4, shrink=4)
ax.text(9.5, 9.48, "EZH2 ⊣ CyclinD1  (central feedback: longer S → more EZH2 → ↓CyclinD1)",
        fontsize=8.5, color=C["fb"], ha="center", fontweight="bold")

# ===================== phase bar (top-right) =====================
pm = [("G0", "#fbeee6", "#b9770e", "p27⁺ / phospho-Rb⁻"),
      ("G1", "#eaf2f8", C["cc_b"], "phospho-Rb⁺, p27⁻"),
      ("S", "#d1f2eb", C["rep_b"], "0 < Dna < 1"),
      ("G2/M", "#fdebd0", C["chk_b"], "Dna≈1, CyclinB-CDK1")]
bx = 11.7
for i, (name, col, ec, sub) in enumerate(pm):
    ax.add_patch(FancyBboxPatch((bx + i*1.25, 9.95), 1.2, 0.32, boxstyle="round,pad=0.01",
                                fc=col, ec=ec, lw=1.2, zorder=5))
    ax.text(bx + i*1.25 + 0.6, 10.11, name, ha="center", va="center", fontsize=8.5, fontweight="bold", zorder=6)
ax.text(bx - 0.15, 10.11, "phases:", ha="right", va="center", fontsize=8.5, fontstyle="italic", color="#555")
# legend
leg = [Line2D([0],[0], color=C["act"], lw=2, label="activation"),
       Line2D([0],[0], color=C["inh"], lw=2, label="inhibition (⊣)"),
       Line2D([0],[0], color=C["fb"], lw=2.4, label="EZH2→CyclinD1 feedback"),
       Line2D([0],[0], color=C["growth_b"], lw=1.6, ls=":", label="size gate"),
       Line2D([0],[0], color=C["drug"], lw=2, label="drug / input")]
ax.legend(handles=leg, loc="lower right", fontsize=8.5, framealpha=0.9, ncol=1).set_zorder(10)

plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), "fig_v44_wiring")
fig.savefig(out + ".png", dpi=160, bbox_inches="tight")
fig.savefig(out + ".pdf", bbox_inches="tight")
print("wrote", out + ".png /.pdf")
