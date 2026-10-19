"""Wireframe — Module 2: the Ccnd1 / Cdk6 loci (EZH2 → PRC2 → H3K27me3 reservoir; proximal PRC2; eraser; eviction;
S-phase dilution; leaky repression of BOTH Cyclin D1 and CDK6).

Boxes = species in the 2026-10 default branch of src/build_model_v44_heldt.py (EZH2_CORE_BLOCK with Cdh1-gated
degradation; serial me-chain with cooperative read-write, continuous fork-gated dilution and H3.3 exchange;
PRC2 catalytic/occupancy split; with_proximal_distal; with_cdk6_gli). No equations.
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from wf import SVG, PAL

W, H = 1400, 800


def build():
    s = SVG(W, H)
    s.text(24, 34, "Module 2 — The Ccnd1 / Cdk6 loci: PRC2 / H3K27me3 governor of the mitogen-drive axis", fs=20, bold=True)
    s.text(24, 56, "EZH2 (cell-cycle-coupled writer) → PRC2 → me1 → me2 → me3 reservoir; proximal PRC2 loads only when Pol II stops; "
                   "transcription evicts PRC2; replication dilutes; repression is leaky and governs Cyclin D1 AND CDK6", fs=11, fill=PAL["grey"], italic=True)

    P, PL = PAL["ezh2_fill"], PAL["ezh2_line"]
    ER = "#1f8a70"

    # ---- EZH2 axis (left column) ----
    e2f_in = s.io(150, 100, "E2F × Cyclin E/A", sub="Module 3 · S/G2 window", w=210)
    cd_in = s.io(150, 170, "Cyclin D1 (mitogen dose)", sub="Module 1/3", w=210)
    ezm = s.node(150, 250, 150, 40, "EZH2 mRNA", fs=12, fill=P, stroke=PL)
    ez = s.node(150, 340, 190, 64, "EZH2 protein", sub="stable · banked in S/G2 and in arrest\nhalved at division", fs=14, fill=P, stroke=PL, sw=2)
    s.edge([e2f_in["b"], ezm["t"]], kind="act", width=1.8)
    s.edge([cd_in["r"], (300, 170), (300, 250), ezm["r"]], kind="act", width=1.4)
    s.label(320, 214, "dose scaling (saturating)", fs=9.5, anchor="start")
    s.edge([ezm["b"], ez["t"]], kind="act", width=1.6)
    s.badge(275, 86, 1)
    s.label(292, 82, "cycle-coupled writer: EZH2 rises in S/G2,\nis degraded in G0/G1, and is diluted at division", fs=9.5, anchor="start", color=PAL["fb"], bold=True)
    cdh1 = s.io(150, 440, "APC/C–Cdh1 (active in G0/G1)", sub="Module 4", w=210, h=40, fs=11)
    s.edge([cdh1["t"], ez["b"]], kind="inh", width=1.6)
    s.label(265, 408, "degrades EZH2 in G0/G1", fs=9, anchor="start")

    # ---- PRC2 ----
    prc2 = s.node(430, 340, 160, 72, "PRC2", sub="EZH2 · EED · SUZ12\nrecruited by sequence + by me3", fs=15, fill="#d7bde2", stroke=PL, sw=2.2)
    s.edge([ez["r"], prc2["l"]], kind="act", width=2.2, color=PL)
    s.label(335, 332, "assembles", fs=9.5)

    # EZH2i: catalytic inhibitor — blocks the SET domain (writing) only; PRC2 stays bound
    ezi = s.pill(430, 240, "tazemetostat (EZH2i)", w=150)
    s.edge([ezi["b"], prc2["t"]], kind="drug", width=2.0)
    s.label(515, 232, "blocks the SET domain (WRITING) only — PRC2 stays bound,\nexisting mark decays over ~1 day (turnover + eraser + dilution)",
            fs=9.5, color=PAL["drug"], anchor="start")

    # ---- methylation chain = the DISTAL RESERVOIR ----
    CY = 470
    xs = [440, 580, 720, 860]
    names = ["me0", "H3K27me1", "H3K27me2", "H3K27me3"]
    subs = ["unmethylated", None, None, "repressive · heritable"]
    fills = ["#ffffff", "#f4ecf7", "#e8daef", "#c39bd3"]
    ch = [s.node(x, CY, 100 if i else 76, 40, n, sub=sub, fs=11.5, fill=f, stroke=PL, sub_fs=9)
          for i, (x, n, sub, f) in enumerate(zip(xs, names, subs, fills))]
    for i in range(3):
        s.edge([(ch[i]["r"][0], CY - 6), (ch[i + 1]["l"][0], CY - 6)], kind="act", width=1.8, color=PL)
        s.edge([(ch[i + 1]["l"][0], CY + 8), (ch[i]["r"][0], CY + 8)], kind="plain", width=1.2)
    for i, lab in enumerate(["fast", "fast", "SLOW · rate-limiting"]):
        s.label((xs[i] + xs[i + 1]) / 2, CY - 14, lab, fs=9.5, color=PL, bold=(i == 2))
    s.line([(500, 428), (500, 420), (812, 420), (812, 428)], color=PL, width=1.4)
    s.edge([prc2["b"], (430, 420), (500, 420)], kind="act", width=1.8, color=PL)
    s.label(656, 412, "PRC2 writes each step", fs=9.5, color=PL, bold=True)
    s.label(656, 441, "DISTAL RESERVOIR = the broad H3K27me3 domain = the memory", fs=9, color=PL, bold=True)

    # erase bracket below
    s.line([(520, 494), (520, 502), (900, 502), (900, 494)], color=ER, width=1.4)
    jm = s.node(1000, 605, 170, 50, "Jmjd3 / Kdm6b", sub="H3K27 demethylase (weak)", fs=12.5, fill="#d5f5e3", stroke=ER)
    gli_in = s.io(1240, 605, "Gli1", sub="Module 1 · mitogen-gated", w=150)
    s.edge([gli_in["l"], jm["r"]], kind="act", width=1.6, color=ER)
    s.edge([jm["t"], (1000, 502), (900, 502)], kind="act", width=1.6, color=ER)
    s.label(760, 517, "eraser strips me3 → me2 → me1 (weak; Gli1-induced)", fs=9.5, color=ER, bold=True)
    s.label(1110, 660, "writer follows the CYCLE, eraser follows MITOGEN", fs=9.5, color=ER)

    # read-write feedback: me3 -> PRC2
    s.edge([ch[3]["t"], prc2["tr"]], kind="fb", width=1.6, curve=(700, 300))
    s.label(760, 292, "EED reads me3 → recruits more PRC2 (cooperative read-write, self-reinforcing)", fs=9.5, color=PAL["fb"], anchor="start")

    # replicative dilution (continuous, fork-gated) + H3.3 exchange
    rep = s.io(1240, 470, "DNA replication (S phase)", sub="Module 4 · fork-gated · longer S dilutes more", w=250)
    s.edge([rep["l"], ch[3]["r"]], kind="fb", width=1.8, color=PAL["blue"])
    s.label(925, 452, "dilutes me1/me2/me3 through S", fs=9.5, color=PAL["blue"], bold=True, anchor="start")
    s.label(925, 464, "(new histones unmethylated)", fs=9, color=PAL["blue"], anchor="start")
    s.badge(1358, 440, 2)
    s.label(1350, 420, "faster cycling dilutes the brake", fs=9.5, color=PAL["fb"], bold=True, anchor="end")
    s.label(1012, 508, "+ H3.3/HIRA exchange: replication-INDEPENDENT\n   me3 loss where the gene is transcribed", fs=9, color=PAL["blue"], anchor="start")

    # ---- promoter / transcription (proximal compartment) ----
    loc = s.locus(700, 640, w=260, n_nuc=5, marked=3)
    # proximal PRC2 occupancy (P_prox), dashed because ≈0 in cycling cells
    px, py = loc["pol2"][0] + 24, loc["y"] - 40
    s.body.append(f'<ellipse cx="{px}" cy="{py}" rx="19" ry="11" fill="#d7bde2" stroke="{PL}" stroke-width="1.4" stroke-dasharray="4,3" opacity="0.9"/>')
    s.text(px, py + 3.5, "PRC2", fs=8, anchor="middle", bold=True, fill=PL)
    s.label(px, py - 16, "proximal", fs=8.5, color=PL)
    d_in = s.io(480, 640, "transcriptional drive", sub="Gli + MYCN · Module 1", w=160, h=44, fs=11.5)
    s.edge([d_in["r"], loc["l"]], kind="act", width=3.0, color=PAL["prom_line"])
    # bound PRC2 represses (leaky) — shared by Ccnd1 and Cdk6
    s.edge([(prc2["bl"][0] + 20, prc2["bl"][1]), (370, 560), (640, 560), (640, 596)], kind="inh", width=2.2, color=PAL["fb"])
    s.dot(370, 560, color=PAL["fb"])
    s.label(455, 550, "bound PRC2 represses (occupancy, not the mark)", fs=10, color=PAL["fb"], bold=True, anchor="start")
    s.label(455, 575, "LEAKY: dampens, never silences (Pol II stays)", fs=9.5, color=PAL["fb"], anchor="start")
    # transcription evicts PRC2
    s.edge([loc["rna"], (780, 540), (490, 540), (490, prc2["b"][1])], kind="inh", width=1.4, color=PAL["teal"], dashed=True)
    s.label(560, 533, "nascent transcription evicts PRC2", fs=9.5, color=PAL["teal"], anchor="start")
    s.badge(830, 532, 3)
    s.label(846, 536, "double-negative: active locus stays poised; vismodegib arrest self-reinforces", fs=9.5, color=PAL["fb"], bold=True, anchor="start")
    # proximal compartment explanation
    s.label(560, 702, "PROXIMAL PRC2 at the TSS (the actual repressor in arrest):", fs=9.5, color=PL, anchor="start", bold=True)
    s.label(560, 715, "loaded from the reservoir only while Pol II is NOT elongating (E2F low);\n"
                      "evicted by elongation → ≈0 in cycling cells (marked-but-expressed);\n"
                      "re-loads in sustained arrest → deeper, REVERSIBLE repression of Ccnd1 + Cdk6", fs=9, color=PL, anchor="start")

    # ---- CDK6: the second governed arm ----
    cdk6 = s.node(370, 720, 160, 52, "CDK6", sub="Gli-driven · slow reservoir", fs=14, fill="#d6eaf8", stroke=PAL["blue"], sw=1.8)
    s.edge([(370, 560), cdk6["t"]], kind="inh", width=2.0, color=PAL["fb"])
    gli6 = s.io(150, 720, "Gli1 / Gli-A", sub="Module 1", w=150, h=44, fs=11.5)
    s.edge([gli6["r"], cdk6["l"]], kind="act", width=1.8)
    s.label(290, 790, "same PRC2 repression as Ccnd1 (leaky floor) → CDK6 is the 2nd governed arm → Module 3 (INK4 sink)", fs=9, color=PAL["blue"], anchor="start")

    # output
    mrna = s.node(1000, 725, 150, 42, "Ccnd1 mRNA", fs=12.5, fill="#eafaf1", stroke=PAL["teal"])
    cd = s.node(1240, 725, 190, 60, "Cyclin D1", sub="→ Module 3", fs=17, fill=PAL["cd_fill"], stroke=PAL["cd_line"], sw=2.6)
    s.edge([(loc["r"][0], 640), (900, 640), (900, 725), mrna["l"]], kind="act", width=3.0, color=PAL["teal"])
    s.edge([mrna["r"], cd["l"]], kind="act", width=3.0, color=PAL["teal"])

    s.legend(40, 560, [("act", "activation / conversion"), ("inh", "inhibition"), ("drug", "drug"),
                       ("fb", "feedback / slow / hand-off")])
    return s


if __name__ == "__main__":
    out = os.path.join(os.path.dirname(__file__), "module2_ccnd1_locus.svg")
    build().save(out)
    print("wrote", out)
