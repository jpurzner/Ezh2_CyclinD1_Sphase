"""Figure: v44 model wiring diagram — CyclinD1-centric, with the Ccnd1 gene split into a separate
PROMOTER (epigenetic locus) and TRANSCRIPT.

CyclinD1 PROTEIN is the central convergence node. The gene-expression axis runs into it:
    Ccnd1 PROMOTER (H3K27me3 / PRC2 locus, gated by writer vs eraser)  ->  Ccnd1 mRNA (transcript)  ->  CyclinD1
Mitogen (Hedgehog/Gli + MYCN) activates the promoter; the EZH2/PRC2 H3K27me3 mark represses it. CyclinD1 then
hands off to the two-step Rb restriction point and the cell-cycle engine.

Reflects the current default model:
  * TWO-STEP Rb: CyclinD-CDK4/6 mono-phosphorylates (Rb -> Rbm), CyclinE/A-CDK2 hyper-phosphorylates (Rbm -> pRb)
    and releases E2F -- the R-point.
  * SERIAL H3K27 chain me1 -> me2 -> me3 at the Ccnd1 promoter (me2->me3 rate-limiting); Mk = me3 represses.
  * PRC2 = EZH2 x (accessory + read-write a_rw*Mk); writer cycle-gated (E2F target), eraser Gli->Jmjd3/Kdm6b
    mitogen-gated; nascent-Ccnd1 RNA evicts PRC2 (double-negative).
  * EZH2 is a CONCENTRATION (preserved at division, reset by ~10h turnover -- not amount-halved).
  * MB-specific BIRTH p27 (P21_div: GNP 0.6 -> MB 1.8) makes MB p27-high/G0-rich; Skp2-p27 toggle is the switch.
  * Feature C: commitment carryover across division (committed daughters re-enter G1, not deep G0).

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
    hh="#d4edda", hh_b="#28a745",           # Hedgehog / Gli (green)
    mycn="#fff3cd", mycn_b="#e67e22",        # MYCN (orange)
    cd="#d0ece7", cd_b="#0e6251",            # CyclinD1 (deep teal, central)
    prom="#fef5e7", prom_b="#b9770e",        # Ccnd1 promoter / locus (amber-gold)
    tx="#eafaf1", tx_b="#117a65",            # transcript (green-teal)
    cc="#d6eaf8", cc_b="#2471a3",            # cell-cycle engine (blue)
    rep="#d1f2eb", rep_b="#17a2b8",          # DNA replication (teal)
    ezh2="#e8daef", ezh2_b="#7d3c98",        # EZH2 / PRC2 / H3K27me3 (purple)
    growth="#eceff1", growth_b="#607d8b",    # growth (slate)
    chk="#fdebd0", chk_b="#ca6f1e",          # checkpoint (amber)
    rb="#fadbd8", rb_b="#a93226",            # Rb states (rose/crimson)
    node="#ffffff",
    act="#229954", inh="#c0392b", fb="#7d3c98", drug="#e74c3c", erase="#1f8a70",
)

fig, ax = plt.subplots(figsize=(18.4, 11.0))
ax.set_xlim(0, 18.4); ax.set_ylim(0, 11.0); ax.set_aspect("equal"); ax.axis("off")


def region(x, y, w, h, fc, ec, label, lx=None, ly=None, fs=10.5):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.15",
                                fc=fc, ec=ec, lw=1.5, alpha=0.35, zorder=0))
    if label:
        ax.text(lx if lx is not None else x + 0.15, ly if ly is not None else y + h - 0.30,
                label, fontsize=fs, fontweight="bold", color=ec, zorder=1, ha="left")


def node(x, y, text, w=1.55, h=0.6, fc=C["node"], ec="#2c3e50", fs=9, bold=False, italic=False, lw=1.3):
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h, boxstyle="round,pad=0.03,rounding_size=0.08",
                                fc=fc, ec=ec, lw=lw, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, zorder=4,
            fontweight="bold" if bold else "normal", fontstyle="italic" if italic else "normal")
    return (x, y, w, h)


def edge_pt(nd, side):
    x, y, w, h = nd
    return {"l": (x - w/2, y), "r": (x + w/2, y), "t": (x, y + h/2), "b": (x, y - h/2)}[side]


def arrow(p1, p2, color=C["act"], style="->", lw=1.8, rad=0.0, ls="-", z=2, shrink=7):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle=style, color=color, lw=lw, ls=ls,
                                 mutation_scale=15, shrinkA=shrink, shrinkB=shrink,
                                 connectionstyle=f"arc3,rad={rad}", zorder=z))


def inhibit(p1, p2, color=C["inh"], lw=1.8, rad=0.0, z=2, shrink=7):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-[", color=color, lw=lw,
                                 mutation_scale=9, shrinkA=shrink, shrinkB=shrink,
                                 connectionstyle=f"arc3,rad={rad}", zorder=z))


# ===================== title =====================
ax.text(9.2, 10.68, "v44 cell-cycle model — CyclinD1-centric wiring", fontsize=16, fontweight="bold", ha="center")
ax.text(9.2, 10.22, "Ccnd1 PROMOTER (epigenetic locus) -> Ccnd1 mRNA (transcript) -> CyclinD1 protein (central node) -> two-step Rb restriction point -> cycle",
        fontsize=9.0, ha="center", color="#555", style="italic")

# ===================== module regions =====================
region(0.25, 6.35, 5.9, 3.9, C["hh"], C["hh_b"], "Mitogenic drive  —  Hedgehog / Gli + MYCN")
region(0.25, 0.35, 6.35, 4.7, C["ezh2"], C["ezh2_b"], "EZH2 / PRC2 / H3K27me3")
region(11.15, 5.55, 7.0, 4.7, C["cc"], C["cc_b"], "Cell-cycle engine  —  two-step R-point → S → mitosis", lx=11.35)
region(11.15, 0.35, 7.0, 4.9, C["rep"], C["rep_b"], "DNA replication + mitotic switch / checkpoint")

# ===================== central node: CyclinD1 protein =====================
cd = node(8.55, 5.55, "CyclinD1\n(Cd protein)", w=2.35, h=1.15, fc=C["cd"], ec=C["cd_b"], fs=12, bold=True, lw=2.4)

# ===================== gene axis: PROMOTER -> mRNA -> CyclinD1 =====================
prom = node(3.05, 5.55, "Ccnd1 PROMOTER\n(H3K27me3 locus)", w=2.55, h=1.0, fc=C["prom"], ec=C["prom_b"], fs=8.6, bold=True, lw=1.8)
mrna = node(6.0, 5.55, "Ccnd1 mRNA\n(transcript)", w=1.85, h=0.72, fc=C["tx"], ec=C["tx_b"], fs=8.4, bold=True)
arrow(edge_pt(prom, "r"), edge_pt(mrna, "l"), C["tx_b"], lw=2.0)
ax.text(4.52, 5.94, "transcribe", fontsize=7.0, color=C["tx_b"], ha="center", style="italic")
ax.text(4.52, 5.16, "(rate = drive × open)", fontsize=6.0, color=C["tx_b"], ha="center", style="italic")
arrow(edge_pt(mrna, "r"), edge_pt(cd, "l"), C["cd_b"], lw=2.2)
ax.text(7.35, 5.9, "translate\n(k_Cd_translation)", fontsize=6.6, color=C["cd_b"], ha="center", style="italic")

# ===================== Hedgehog / MYCN drive (into the promoter) =====================
shh = node(0.9, 9.55, "SHH", w=0.9, h=0.5, fc="#fef9e7", fs=8.2)
smo = node(2.75, 9.55, "Smo → Gli(A)", w=1.95, h=0.58, fc=C["hh"], ec=C["hh_b"], fs=8.2, bold=True)
ptch = node(1.05, 8.5, "Ptch1\n(Gli target)", w=1.35, h=0.62, fs=7.0)
mycn = node(4.9, 9.35, "MYCN", w=1.2, h=0.58, fc=C["mycn"], ec=C["mycn_b"], fs=8.5, bold=True)
hhi = node(2.75, 10.28, "HHi (vismo)", w=1.55, h=0.46, fc="#fdedeb", ec=C["drug"], fs=7.4, bold=True)
arrow(edge_pt(shh, "r"), edge_pt(smo, "l"), C["act"])
arrow((smo[0]-0.55, smo[1]-0.29), (ptch[0]+0.45, ptch[1]+0.30), C["hh_b"], rad=-0.30, lw=1.4)  # Gli→Ptch1
inhibit((ptch[0]+0.15, ptch[1]+0.33), (smo[0]-0.62, smo[1]-0.26), C["inh"], rad=0.30, shrink=3)  # Ptch1⊣Smo
inhibit((hhi[0], hhi[1]-0.23), (smo[0], smo[1]+0.29), C["drug"], shrink=4)
ax.text(0.4, 7.95, "neg-fb Gli→Ptch1⊣Smo (broken in MB: Ptch1 f 1.0→~0.1 → Gli↑)", fontsize=6.1, color=C["hh_b"], ha="left", fontweight="bold")
# Gli + MYCN activate the PROMOTER
arrow((smo[0]-0.3, smo[1]-0.32), (prom[0]-0.2, prom[1]+0.52), C["act"], rad=0.16, lw=1.7)   # Gli → promoter
arrow((mycn[0], mycn[1]-0.31), (prom[0]+0.7, prom[1]+0.52), C["mycn_b"], rad=-0.12, lw=1.7)  # MYCN → promoter
ax.text(3.05, 6.72, "activate transcription", fontsize=6.6, color=C["hh_b"], ha="center", style="italic")

# ===================== EZH2 / PRC2 / H3K27me3 (writes the promoter mark) =====================
ezh2m = node(1.35, 4.15, "EZH2 mRNA\n(E2F target)", w=1.9, h=0.6, fc="#ffffff", ec=C["ezh2_b"], fs=7.6)
ezh2 = node(1.35, 3.4, "EZH2 protein\n(concentration:\npreserved at ÷)", w=1.9, h=0.86, fc="#f4ecf7", ec=C["ezh2_b"], fs=7.2, bold=True)
prc2 = node(1.55, 1.55, "PRC2 complex\n(EZH2 × [a0 + a_rw·me3])", w=2.5, h=0.82, fc="#d7bde2", ec=C["ezh2_b"], fs=7.6, bold=True)
mark = node(4.55, 2.55, "H3K27me3 @ Ccnd1\nme1 → me2 → me3\n(me2→me3 rate-limiting)", w=2.9, h=1.05, fc="#e8daef", ec=C["ezh2_b"], fs=7.8, bold=True)
jmjd3 = node(4.7, 4.35, "Jmjd3 / Kdm6b\n(Gli1-induced eraser)", w=2.2, h=0.66, fc="#d5f5e3", ec=C["erase"], fs=7.2, bold=True)
ezi = node(0.95, 0.72, "EZH2i (taz)", w=1.5, h=0.46, fc="#fdedeb", ec=C["drug"], fs=7.4, bold=True)
arrow(edge_pt(ezh2m, "b"), edge_pt(ezh2, "t"), C["act"])
arrow(edge_pt(ezh2, "b"), (prc2[0]-0.2, prc2[1]+0.42), C["ezh2_b"], rad=0.12, lw=1.7)
inhibit((ezi[0]+0.35, ezi[1]+0.22), (prc2[0]-0.55, prc2[1]-0.30), C["drug"], rad=-0.2, shrink=3)
# write / read-write between PRC2 and the mark chain
arrow((prc2[0]+1.15, prc2[1]+0.18), (mark[0]-1.15, mark[1]-0.30), C["ezh2_b"], rad=-0.12, lw=1.9)  # write
ax.text(3.1, 1.72, "methylate (write)", fontsize=6.4, color=C["ezh2_b"], ha="center", style="italic")
arrow((mark[0]-1.25, mark[1]-0.05), (prc2[0]+1.1, prc2[1]+0.35), C["ezh2_b"], rad=-0.35, lw=1.6, ls=(0,(4,2)))  # read-write
ax.text(2.35, 2.62, "read-write\nrecruit (a_rw·me3)", fontsize=6.1, color=C["ezh2_b"], ha="center", style="italic")
inhibit((jmjd3[0]+0.2, jmjd3[1]-0.35), (mark[0]+0.2, mark[1]+0.55), C["erase"], shrink=3)   # Jmjd3 ⊣ mark
ax.text(5.85, 3.5, "erase", fontsize=6.5, color=C["erase"], ha="left", style="italic")
arrow((smo[0]+0.9, smo[1]-0.2), (jmjd3[0]+0.2, jmjd3[1]+0.36), C["erase"], rad=-0.28, lw=1.5, ls=(0,(5,2)))  # Gli→Jmjd3
ax.text(6.0, 5.15, "Gli→Jmjd3: mitogen-gated eraser\n(high-Gli MB strips mark → ChIP MB<GNP)", fontsize=6.3, color=C["erase"], ha="left", fontweight="bold")
# H3K27me3 basal turnover loop (the S-phase ÷2 arrow from DNA replication is drawn later, once `dna` exists)
arrow((mark[0]-1.5, mark[1]+0.15), (mark[0]-1.5, mark[1]-0.15), C["rep_b"], rad=-1.7, lw=1.3, shrink=2)   # + turnover loop
ax.text(2.75, 2.55, "+ turnover", fontsize=5.8, color=C["rep_b"], ha="center", va="center", style="italic")

# ===================== the epigenetic gate on the PROMOTER =====================
inhibit((mark[0]-0.4, mark[1]+0.55), (prom[0]+0.0, prom[1]-0.52), C["fb"], rad=-0.15, lw=2.4, shrink=4)  # me3 ⊣ promoter
ax.text(6.75, 4.62, "H3K27me3 closes the Ccnd1 promoter\n(↓ transcription; Hill on PRC2, floor f0)",
        fontsize=6.6, color=C["fb"], ha="center", fontweight="bold")
# nascent Ccnd1 RNA ⊣ PRC2 (eviction; double-negative)
inhibit((mrna[0]-0.35, mrna[1]-0.36), (prc2[0]+1.0, prc2[1]+0.42), C["tx_b"], rad=0.45, shrink=3)
ax.text(6.6, 3.95, "nascent RNA ⊣ PRC2 (eviction)", fontsize=6.0, color=C["tx_b"], ha="left", style="italic")

# ===================== cell-cycle engine: two-step Rb R-point =====================
# CyclinD1 -> CDK4/6 -> Rb->Rbm (mono-P);  CDK2 -> Rbm->pRb (hyper-P) -> E2F
rb = node(11.9, 8.9, "Rb", w=0.95, h=0.55, fc=C["rb"], ec=C["rb_b"], fs=8.5, bold=True)
rbm = node(13.35, 8.9, "Rbm\n(mono-P)", w=1.25, h=0.66, fc="#f1948a", ec=C["rb_b"], fs=7.6, bold=True)
prb = node(14.9, 8.9, "pRb\n(hyper-P)", w=1.3, h=0.66, fc="#cb4335", ec=C["rb_b"], fs=7.6, bold=True)
e2f = node(15.05, 7.55, "E2F", w=1.05, h=0.58, fc="#fdfefe", ec=C["cc_b"], fs=9.5, bold=True)
ce = node(16.9, 8.5, "CyclinE/A\n(CDK2)", w=1.6, h=0.72, fc=C["cc"], ec=C["cc_b"], fs=8.2, bold=True)
arrow(edge_pt(rb, "r"), edge_pt(rbm, "l"), C["rb_b"], lw=1.7)
arrow(edge_pt(rbm, "r"), edge_pt(prb, "l"), C["rb_b"], lw=1.7)
arrow((prb[0]-0.1, prb[1]-0.33), (e2f[0]+0.1, e2f[1]+0.30), C["act"], rad=-0.1)   # pRb releases E2F
arrow(edge_pt(e2f, "r"), (ce[0]-0.55, ce[1]-0.3), C["act"], rad=-0.15)             # E2F → CyclinE/A
# CyclinD1 → CDK4/6 → Rb->Rbm (the mono-P priming step)
arrow((cd[0]+1.0, cd[1]+0.45), (rb[0]-0.5, rb[1]-0.35), C["cd_b"], rad=-0.28, lw=2.2)
ax.text(10.55, 7.55, "CDK4/6\n(mono-P)", fontsize=6.8, color=C["cd_b"], ha="center", style="italic", fontweight="bold")
# CDK2 → Rbm->pRb (the hyper-P commitment step)
arrow((ce[0]-0.2, ce[1]-0.38), (rbm[0]+0.35, rbm[1]+0.35), C["cc_b"], rad=0.3, lw=1.7)
ax.text(16.1, 8.98, "CDK2 (hyper-P)", fontsize=6.6, color=C["cc_b"], ha="center", style="italic")
ax.text(13.4, 9.72, "TWO-STEP Rb R-point:  CyclinD-CDK4/6 primes (mono-P) → CyclinE/A-CDK2 commits (hyper-P) → E2F released",
        fontsize=7.4, color=C["rb_b"], ha="center", fontweight="bold")

# INK4 + palbo brake on CDK4/6
ink4 = node(10.35, 9.35, "p16 + p18\n(INK4)", w=1.5, h=0.6, fc="#f4ecf7", ec=C["ezh2_b"], fs=7.2, bold=True)
inhibit((ink4[0]+0.1, ink4[1]-0.33), (rb[0]-0.35, rb[1]+0.30), C["inh"], rad=-0.1, shrink=3)
palbo = node(9.7, 8.15, "CDK4/6i\n(palbo)", w=1.3, h=0.58, fc="#fdedeb", ec=C["drug"], fs=7.0, bold=True)
inhibit((palbo[0]+0.35, palbo[1]+0.25), (10.55, 7.9), C["drug"], shrink=3)
ax.text(9.7, 7.72, "Vmax block", fontsize=6.0, color=C["drug"], ha="center", style="italic")

# ===================== Skp2-p27 toggle + MB birth-p27 =====================
p27 = node(11.75, 6.6, "p27 (=P21)", w=1.5, h=0.58, fc="#fbeee6", ec="#b9770e", fs=8.2, bold=True)
skp2 = node(13.4, 6.55, "Skp2", w=1.05, h=0.52, fs=8.2)
inhibit((p27[0]+0.55, p27[1]-0.15), (skp2[0]-0.5, skp2[1]+0.1), C["inh"], rad=-0.15, shrink=3)  # p27 side
inhibit((skp2[0]-0.35, skp2[1]-0.2), (p27[0]+0.45, p27[1]-0.28), C["inh"], rad=0.3, shrink=3)   # Skp2 ⊣ p27 (degrade)
arrow((e2f[0]-0.2, e2f[1]+0.28), (skp2[0]+0.2, skp2[1]-0.28), C["act"], rad=-0.35, lw=1.3)       # E2F → Skp2
inhibit((p27[0]+0.4, p27[1]+0.28), (rbm[0]-0.3, rbm[1]-0.34), C["inh"], rad=-0.2, shrink=3)      # p27 ⊣ CDK2
ax.text(11.75, 5.95, "birth p27:  GNP 0.6 → MB 1.8\n→ MB p27-high / G0-rich (fold-safe)", fontsize=6.6, color="#b9770e", ha="center", fontweight="bold")

# ===================== S-phase / replication / mitosis =====================
ori = node(16.7, 6.35, "origins fire\n(Rc→aRc)", w=1.55, h=0.66, fc="#ffffff", ec=C["rep_b"], fs=7.6)
arrow((ce[0]+0.3, ce[1]-0.38), (ori[0]+0.1, ori[1]+0.35), C["act"], rad=-0.2)
dna = node(14.3, 3.3, "DNA replication\nDna: 0 → 1", w=2.05, h=0.72, fc="#ffffff", ec=C["rep_b"], fs=8.2, bold=True)
arrow(edge_pt(ori, "b"), (dna[0]+0.7, dna[1]+0.4), C["act"], rad=0.25)
# DNA replication HALVES the H3K27me3 mark (replicative dilution at S)
arrow((dna[0]-1.05, dna[1]-0.35), (mark[0]+1.5, mark[1]-0.32), C["rep_b"], rad=0.45, lw=1.6, ls=(0,(5,2)))
ax.text(9.35, 1.72, "S-PHASE ÷2:  DNA replication halves H3K27me3\n(DNA doubles, new histones unmethylated → dilution)",
        fontsize=6.3, color=C["rep_b"], ha="center", va="center", style="italic", fontweight="bold")
hu = node(11.95, 3.3, "HU", w=0.85, h=0.5, fc="#fdedeb", ec=C["drug"], fs=8.2, bold=True)
inhibit(edge_pt(hu, "r"), (dna[0]-1.05, dna[1]), C["drug"], shrink=3)
ax.text(12.05, 2.7, "↓ fork speed → S↑", fontsize=6.6, color=C["drug"], ha="center", style="italic")
mpf = node(16.75, 2.35, "CyclinB/CDK1 (MPF)\nCdc25 / Wee1", w=2.5, h=0.9, fc="#ffffff", ec=C["chk_b"], fs=7.8, bold=True)
chk = node(15.9, 1.15, "CHK1\n(unfinished S)", w=1.7, h=0.55, fc="#fdebd0", ec=C["chk_b"], fs=7.2)
arrow((dna[0]+1.0, dna[1]-0.1), (mpf[0]-1.3, mpf[1]+0.05), C["act"], rad=-0.12)   # Dna done → MPF
inhibit((chk[0]+0.5, chk[1]+0.27), (mpf[0]-0.75, mpf[1]-0.42), C["inh"], shrink=3)
arrow((dna[0]+0.6, dna[1]-0.38), (chk[0]-0.55, chk[1]+0.15), C["act"], rad=0.2, lw=1.2)  # forks → CHK1
# ===================== MITOSIS / division event (E_div) =====================
mit = node(12.35, 1.35, "MITOSIS ÷ division  (E_div)", w=2.85, h=0.6, fc="#f5b7b1", ec="#922b21", fs=8.0, bold=True, lw=1.9)
arrow((mpf[0]-1.3, mpf[1]-0.3), (mit[0]+1.35, mit[1]+0.02), "#922b21", rad=0.2, lw=1.9)             # MPF crosses threshold → mitosis
arrow((mit[0]+1.1, mit[1]+0.25), (dna[0]-0.15, dna[1]-0.38), "#922b21", rad=-0.15, lw=1.3, ls=(0,(3,2)))  # Dna → 0 (re-license)
arrow((mit[0]-0.6, mit[1]+0.28), (cd[0]+0.85, cd[1]-0.56), C["cc_b"], rad=-0.28, lw=2.2)            # daughter restarts → CyclinD1 / G1
ax.text(10.15, 4.05, "daughter\nre-enters G1", fontsize=6.3, color=C["cc_b"], ha="center", style="italic", fontweight="bold")
ax.text(8.15, 0.5, "at division (E_div):  Dna→0 · mass÷2 · MPF/Cdc20→0 · cyclins→low · p27→BIRTH (P21_div; MB 1.8)\n"
        "PRESERVED across ÷:  EZH2 = concentration (NOT halved) · Feature-C carryover (pRb / CycE / low-p27)",
        fontsize=6.1, color="#922b21", ha="center", fontweight="bold")

# ===================== growth / size gates =====================
region(6.75, 0.35, 4.15, 3.15, C["growth"], C["growth_b"], "Cell growth & size", fs=9.5)
mass = node(8.8, 2.55, "cell mass\n(grows, ÷2 at ÷)", w=2.2, h=0.68, fc="#ffffff", ec=C["growth_b"], fs=7.8)
mcommit = node(7.75, 1.25, "M_commit\n(G0→G1)", w=1.55, h=0.6, fc="#ffffff", ec=C["growth_b"], fs=7.2)
msize = node(9.85, 1.25, "M_size\n(→ S-entry)", w=1.55, h=0.6, fc="#ffffff", ec=C["growth_b"], fs=7.2)
arrow((mass[0]-0.5, mass[1]-0.3), edge_pt(mcommit, "t"), C["growth_b"], ls=":", lw=1.3, rad=0.15)
arrow((mass[0]+0.5, mass[1]-0.3), edge_pt(msize, "t"), C["growth_b"], ls=":", lw=1.3, rad=-0.15)
arrow((mit[0]-1.35, mit[1]+0.12), (mass[0]+0.55, mass[1]-0.34), "#922b21", rad=0.22, lw=1.3, ls=(0,(3,2)))  # mitosis: mass ÷2
arrow((mcommit[0]+0.2, mcommit[1]+0.32), (10.9, 6.9), C["growth_b"], ls=":", lw=1.4, rad=-0.2)  # commit gate → R-point
arrow(edge_pt(msize, "r"), (ori[0]-0.6, ori[1]-0.3), C["growth_b"], ls=":", lw=1.4, rad=-0.25)   # size gate → S entry

# ===================== E2F → EZH2 writer feedback (cycle-gated) =====================
arrow((e2f[0]-0.35, e2f[1]-0.2), (ezh2m[0]+0.6, ezh2m[1]+0.28), C["ezh2_b"], rad=0.30, lw=1.5, ls=(0,(6,2)))
ax.text(11.7, 4.45, "E2F drives EZH2\n(writer is cycle-gated)", fontsize=6.6, color=C["ezh2_b"], ha="center", style="italic", fontweight="bold")

# ===================== phase bar + legend =====================
pm = [("G0", "#fbeee6", "#b9770e"), ("G1", "#eaf2f8", C["cc_b"]), ("S", "#d1f2eb", C["rep_b"]), ("G2/M", "#fdebd0", C["chk_b"])]
bx = 14.6
for i, (name, col, ec) in enumerate(pm):
    ax.add_patch(FancyBboxPatch((bx + i*0.9, 10.42), 0.86, 0.30, boxstyle="round,pad=0.01", fc=col, ec=ec, lw=1.2, zorder=5))
    ax.text(bx + i*0.9 + 0.43, 10.57, name, ha="center", va="center", fontsize=8.0, fontweight="bold", zorder=6)
ax.text(bx - 0.15, 10.57, "phases:", ha="right", va="center", fontsize=8.0, fontstyle="italic", color="#555")

leg = [Line2D([0],[0], color=C["act"], lw=2, label="activation"),
       Line2D([0],[0], color=C["inh"], lw=2, label="inhibition (⊣)"),
       Line2D([0],[0], color=C["fb"], lw=2.4, label="H3K27me3 closes the Ccnd1 promoter"),
       Line2D([0],[0], color=C["ezh2_b"], lw=1.6, ls=(0,(4,2)), label="PRC2 read-write / E2F→EZH2 writer"),
       Line2D([0],[0], color=C["erase"], lw=1.6, ls=(0,(5,2)), label="Gli→Jmjd3 mark eraser"),
       Line2D([0],[0], color=C["cd_b"], lw=2.2, label="CyclinD1 axis / hand-off"),
       Line2D([0],[0], color=C["rb_b"], lw=1.7, label="two-step Rb (mono-P → hyper-P)"),
       Line2D([0],[0], color=C["growth_b"], lw=1.6, ls=":", label="size gate"),
       Line2D([0],[0], color=C["drug"], lw=2, label="drug / input")]
ax.legend(handles=leg, loc="lower right", fontsize=7.6, framealpha=0.93, ncol=2,
          columnspacing=1.1, handlelength=1.8).set_zorder(10)

plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), "fig_v44_wiring")
fig.savefig(out + ".png", dpi=160, bbox_inches="tight")
fig.savefig(out + ".pdf", bbox_inches="tight")
print("wrote", out + ".png /.pdf")
