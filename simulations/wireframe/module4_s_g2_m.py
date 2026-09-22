"""Wireframe — Module 4: S / G2 / M engine (origin firing, explicit replication, CHK1, Cyclin B–CDK1 switch, APC/C, division).

Boxes = Heldt replication core + MITOSIS_BLOCK / HU_BLOCK / GROWTH_BLOCK of src/build_model_v44_heldt.py. No equations.
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from wf import SVG, PAL

W, H = 1400, 800


def build():
    s = SVG(W, H)
    s.text(24, 34, "Module 4 — S / G2 / M engine: explicit replication, checkpoint, mitotic switch, division", fs=20, bold=True)
    s.text(24, 56, "CDK2 + size gate fire origins → forks synthesise DNA at a fork speed (HU slows it) → CHK1 holds mitosis until "
                   "replication completes → Cyclin B–CDK1 switch → APC/C → division resets the daughter", fs=11, fill=PAL["grey"], italic=True)

    T, TL = "#d1f2eb", "#17a2b8"
    A, AL = "#fdebd0", "#ca6f1e"
    RY = 330

    # ---- inputs ----
    cdk2 = s.io(150, RY, "Cyclin E/A–CDK2", sub="Module 3", w=180, h=50, fs=12.5)
    sg = s.io(150, 430, "size gate (cell mass)", sub="Module 3", w=180, h=44, fs=11.5)
    hu = s.pill(330, 200, "hydroxyurea (HU)", w=140)

    # ---- origin firing ----
    ori = s.node(420, RY, 190, 64, "origin firing", sub="licensed → fired → active forks\n(+ PCNA)", fs=13, fill=T, stroke=TL)
    s.edge([cdk2["r"], ori["l"]], kind="act", width=2.4, color=PAL["blue"])
    s.edge([sg["r"], (300, 430), (300, RY + 16), (ori["l"][0], RY + 16)], kind="plain", width=1.4, dashed=True, color="#607d8b")
    s.gate(300, 380, color="#607d8b", label="size", lx=282)
    s.edge([(hu["l"][0] + 20, hu["b"][1]), (400, ori["t"][1])], kind="drug", width=1.6)
    s.label(345, 245, "blocks new firing\n(replication stress)", fs=9, color=PAL["drug"], anchor="end")

    # ---- DNA synthesis ----
    dna = s.node(680, RY, 200, 64, "DNA replication", sub="forks copy the genome\nat a set fork speed", fs=13, fill=T, stroke=TL, sw=2)
    s.edge([ori["r"], dna["l"]], kind="act", width=2.4, color=TL)
    s.edge([(hu["r"][0] - 20, hu["b"][1]), (660, dna["t"][1])], kind="drug", width=1.6)
    s.label(618, 290, "slows forks → S lengthens (G1 unchanged)", fs=9, color=PAL["drug"], anchor="end")
    p27 = s.io(560, 190, "p27 (Module 3) ⊣ PCNA·forks", w=200, h=34, fs=10.5)
    s.edge([p27["b"], (640, dna["t"][1])], kind="inh", width=1.2, opacity=0.7)
    # outputs of replication to Module 2
    rep_out = s.io(680, 470, "halves H3K27me3 at Ccnd1", sub="Module 2 · replicative dilution", w=230)
    s.edge([dna["b"], rep_out["t"]], kind="fb", width=1.8, color=PAL["blue"])

    # ---- CHK1 ----
    chk = s.node(900, 200, 150, 48, "CHK1", sub="active while forks run", fs=13, fill=A, stroke=AL)
    s.edge([(dna["tr"][0] - 30, dna["tr"][1]), (750, 200), chk["l"]], kind="act", width=1.6, color=AL)

    # ---- Cyclin B / MPF switch ----
    preb = s.node(960, RY, 180, 64, "Cyclin B–CDK1 (inactive)", sub="Tyr15-P · made as replication\ncompletes (real G2)", fs=11.5, fill="#ffffff", stroke=AL)
    s.edge([dna["r"], preb["l"]], kind="act", width=2.2, color=TL)
    s.label(920, 300, "DNA complete", fs=9.5, color=TL)
    sw = s.node(1180, 230, 220, 56, "Cdc25 ⇄ Wee1", sub="bistable Tyr15 switch", fs=13, fill=A, stroke=AL)
    mpf = s.node(1230, RY + 100, 180, 64, "Cyclin B–CDK1 (MPF)", sub="active · mitotic entry", fs=12.5, fill="#f5b7b1", stroke="#922b21", sw=2.2)
    s.edge([preb["r"], (1180, RY), (1180, mpf["t"][1])], kind="act", width=2.4, color="#922b21")
    s.edge([sw["b"], (1180, RY - 8)], kind="act", width=1.6, color=AL)
    s.label(1172, 290, "Cdc25 activates · Wee1 inactivates", fs=9, color=AL, anchor="end")
    s.edge([(1270, mpf["t"][1]), (1270, sw["b"][1])], kind="fb", width=1.4, color="#922b21")
    s.label(1282, 322, "MPF → Cdc25\nMPF ⊣ Wee1\n(hysteresis)", fs=9, color="#922b21", anchor="start")
    s.edge([chk["r"], sw["l"]], kind="inh", width=2.0)
    s.label(1020, 168, "unfinished S holds the switch OFF", fs=9.5, bold=True, anchor="middle")

    # ---- APC/C-Cdc20 ----
    cdc20 = s.node(1000, RY + 100, 170, 52, "APC/C–Cdc20", sub="degrades Cyclin A / B", fs=12.5, fill="#ffffff", stroke="#922b21")
    s.edge([(mpf["l"][0], RY + 112), (cdc20["r"][0], RY + 112)], kind="act", width=1.6, color="#922b21")
    s.edge([(cdc20["r"][0], RY + 88), (mpf["l"][0], RY + 88)], kind="inh", width=1.6)
    s.label(1060, 478, "MPF activates Cdc20; Cdc20 destroys cyclins → mitotic exit", fs=9, anchor="middle")

    # ---- division event ----
    div = s.callout(880, 520, 470, [
        "DNA → 0, origins re-licensed",
        "Rb → hypo-P, E2F re-bound (partial carry-over if committed)",
        "p27 → BIRTH level (GNP low, MB higher) → transient G0",
        "Skp2, Cyclin E/A, MPF, Cdc20 → low",
        "cell mass ÷ 2",
        "INHERITED: H3K27me3 at Ccnd1 · EZH2 concentration",
    ], title="DIVISION  (MPF crosses threshold)", color="#922b21", fill="#fdf2f0", fs=10)
    s.edge([mpf["b"], (1230, 520)], kind="act", width=2.6, color="#922b21")
    out1 = s.io(360, 545, "daughters re-enter G0/G1", sub="Module 3 · p27-high transient G0", w=260)
    s.edge([(div["l"][0], 570), out1["r"]], kind="fb", width=1.8, color="#922b21")
    out2 = s.io(360, 605, "mark + EZH2 carried into next cycle", sub="Module 2", w=260)
    s.edge([(div["l"][0], 600), out2["r"]], kind="fb", width=1.6, color="#922b21")

    # ---- phase bar ----
    s.phase_bar(300, 1340, 730, [("S-entry", 1, "#eaf2f8", PAL["blue"]), ("S phase", 2.2, T, TL),
                                 ("G2", 1.4, A, AL), ("M → division", 1.4, "#f5b7b1", "#922b21")])
    s.label(300, 770, "S-phase length is emergent (fork speed × genome), so it can be stretched by HU without touching G1 — "
                      "and a longer S means more EZH2 (Module 2).", fs=9.5, anchor="start")

    s.legend(40, 660, [("act", "activation / conversion"), ("inh", "inhibition"), ("drug", "drug"),
                       ("fb", "hand-off / feedback")])
    return s


if __name__ == "__main__":
    out = os.path.join(os.path.dirname(__file__), "module4_s_g2_m.svg")
    build().save(out)
    print("wrote", out)
