"""Wireframe — Module 1: Mitogenic drive (SHH ⊣ Ptch1 ⊣ Smo → Gli → {Ptch1, Gli1, MYCN} → Ccnd1 drive).

Boxes = species in HH_MYCN_BLOCK of src/build_model_v44_heldt.py; arrows = fluxes. No equations.
Run:  python3 simulations/wireframe/build_all.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from wf import SVG, PAL

W, H = 1400, 800


def build():
    s = SVG(W, H)
    YM, HM = 112, 26
    s.membrane(YM, HM, x0=0, x1=W)
    s.text(24, 34, "Module 1 — Mitogenic drive", fs=20, bold=True)
    s.text(24, 56, "SHH ⊣ Ptch1 ⊣ Smo → Gli → {Ptch1, Gli1, MYCN} → Ccnd1 transcriptional drive", fs=11,
           fill=PAL["grey"], italic=True)

    # ---- membrane proteins ----
    shh = s.ligand(300, 82, "SHH  (ligand input)", "#f6d55c")
    ptch = s.receptor(300, YM, HM, 12, "Ptch1", PAL["hh_line"])
    s.receptor(130, YM, HM, 12, "SHH·Ptch1 (inactive)", PAL["hh_line"], faded=True)
    s.body.append(f'<circle cx="130" cy="{YM - 2}" r="9" fill="#f6d55c" opacity="0.5" stroke="{PAL["ink"]}" stroke-width="0.8"/>')
    smo = s.receptor(560, YM, HM, 7, "Smo", PAL["hh_line"])
    vismo = s.pill(560, 48, "vismodegib", w=104)

    s.edge([shh["b"], ptch["t"]], kind="act", width=1.6)
    s.edge([(255, 190), (185, 190)], kind="act", width=1.4)
    s.edge([(185, 202), (255, 202)], kind="act", width=1.2, opacity=0.6)
    s.label(220, 183, "SHH binds", fs=9.5)
    s.label(220, 214, "release", fs=9.5)
    s.edge([(130, 178), (130, 222)], kind="act", width=1.2, opacity=0.6)
    s.label(130, 238, "degraded", fs=9.5)
    s.label(130, 262, "bound Ptch1 no longer\nrepresses Smo", fs=9.5)

    s.edge([ptch["r"], smo["l"]], kind="inh", width=2.0)
    s.label(430, 188, "Ptch1 brake ∝ FUNCTIONAL Ptch1", fs=10, color=PAL["hh_line"], bold=True)
    s.label(430, 203, "GNP: intact  ·  MB: lost (Ptch1 mutant)", fs=9.5, color=PAL["hh_line"])
    s.edge([vismo["b"], smo["t"]], kind="drug", width=2.0)
    s.label(650, 52, "Smo antagonist", fs=9.5, color=PAL["drug"], anchor="start")

    # ---- Gli toggle ----
    gliA = s.tf(480, 300, "Gli-A", PAL["hh_fill"], PAL["hh_line"], w=92, h=42, sub="activator")
    gliR = s.tf(640, 300, "Gli-R", "#f4f6f7", PAL["grey"], w=92, h=42, sub="repressor")
    s.edge([(gliR["l"][0], 291), (gliA["r"][0], 291)], kind="act", width=1.8)
    s.edge([(gliA["r"][0], 311), (gliR["l"][0], 311)], kind="plain", width=1.4)
    s.label(560, 283, "Smo-driven", fs=9.5, color=PAL["act"])
    s.label(560, 327, "reverts", fs=9.5)
    s.edge([smo["b"], (560, 270)], kind="act", width=2.0)
    s.label(560, 352, "one Gli pool, two states", fs=9.5)

    # ---- Ptch1 mRNA + negative feedback ----
    pm = s.node(300, 300, 124, 40, "Ptch1 mRNA", fs=12, fill=PAL["hh_fill"], stroke=PAL["hh_line"])
    s.edge([pm["t"], ptch["b"]], kind="act", width=1.6)
    s.label(312, 250, "translated", fs=9.5, anchor="start")
    s.edge([gliA["l"], pm["r"]], kind="act", width=1.8)
    s.label(397, 282, "Gli → Ptch1", fs=10, color=PAL["act"], bold=True)
    s.label(397, 320, "negative feedback", fs=9.5, color=PAL["act"])
    s.label(300, 345, "loop LOST in MB: Ptch1 still transcribed\n(high Ptch1 mRNA = SHH-MB marker) but non-functional",
            fs=9.5, color=PAL["hh_line"], anchor="middle")

    # ---- Gli1 ----
    g1m = s.node(480, 440, 124, 40, "Gli1 mRNA", fs=12, fill=PAL["hh_fill"], stroke=PAL["hh_line"])
    g1 = s.node(480, 530, 124, 44, "Gli1", sub="pathway readout", fs=13, fill=PAL["hh_fill"], stroke=PAL["hh_line"])
    s.edge([gliA["b"], g1m["t"]], kind="act", width=1.8)
    s.edge([g1m["b"], g1["t"]], kind="act", width=1.4)
    s.label(580, 436, "Gli-A activates, Gli-R represses", fs=9.5, anchor="start")
    s.label(580, 532, "Gli1  MB/GNP ≈ 7×", fs=10, color=PAL["hh_line"], anchor="start", bold=True)

    mem = s.node(250, 500, 190, 62, "Gli1 locus memory", sub="MB only · charged by Smo · slow (~14 h)",
                 fs=12, fill=PAL["mem_fill"], stroke=PAL["mem_line"], dashed=True)
    s.edge([mem["r"], g1m["l"]], kind="fb", width=1.4, curve=(400, 470))
    s.label(250, 555, "optional: explains the Gli1 residual\nthat persists ~1 day after vismodegib", fs=9.5)

    # ---- G junction + bus ----
    GY = 600
    G = s.node(480, GY, 150, 34, "Gli-A + Gli1", fs=11.5, fill="#ffffff", stroke=PAL["hh_line"], bold=False)
    s.edge([g1["b"], G["t"]], kind="act", width=1.6)
    s.edge([gliA["br"], (560, 350), (560, GY - 24), (G["tr"][0] - 14, G["tr"][1])], kind="act", width=1.2, opacity=0.7)
    s.label(566, 575, "Gli-A", fs=9, color=PAL["act"], anchor="start")

    prom = s.node(1150, 395, 190, 60, "Ccnd1 promoter", sub="hand-off to Module 2 (chromatin)",
                  fs=14, fill=PAL["prom_fill"], stroke=PAL["prom_line"], sw=2.2)
    s.dna(1065, 340, 170)
    s.edge([G["r"], (1100, GY), (1100, prom["b"][1])], kind="act", width=3.0)
    s.label(620, GY - 8, "Gli drive  (≈ 86 % of MB Cyclin D1)", fs=10, color=PAL["act"], bold=True, anchor="start")
    s.edge([gliR["t"], (gliR["t"][0], 222), (1020, 222), (1020, 378), (prom["l"][0], 378)], kind="inh", width=1.2, opacity=0.7)
    s.label(820, 216, "Gli-R dampens", fs=9.5)
    s.edge([(1300, 395), prom["r"]], kind="act", width=1.4, color=PAL["prom_line"])
    s.label(1306, 399, "basal", fs=10, anchor="start")

    # ---- MYCN ----
    MX = 990
    mycn = s.tf(MX, 470, "MYCN", PAL["mycn_fill"], PAL["mycn_line"], w=100, h=46, sub="TF")
    s.dot(MX, GY)
    s.edge([(MX, GY), mycn["b"]], kind="act", width=1.6)
    s.label(MX + 8, 560, "Gli target", fs=9.5, anchor="start")
    s.edge([(MX, 400), mycn["t"]], kind="act", width=1.6, color=PAL["mycn_line"])
    s.label(MX - 12, 400, "basal × amplification", fs=10, anchor="end", color=PAL["mycn_line"])
    s.edge([mycn["r"], (1065, 470), (1065, prom["b"][1])], kind="act", width=2.0, color=PAL["mycn_line"])
    s.label(1140, 503, "MYCN drive", fs=9.5, color=PAL["mycn_line"], bold=True, anchor="end")
    s.label(940, 528, "mitogen-INDEPENDENT floor · amplified in MB", fs=9.5, color=PAL["mycn_line"], anchor="end", bold=True)
    s.label(940, 544, "→ vismodegib is only a PARTIAL withdrawal in MB", fs=9.5, color=PAL["mycn_line"], anchor="end")

    # ---- output ----
    mrna = s.node(1150, 650, 150, 42, "Ccnd1 mRNA", fs=12.5, fill="#eafaf1", stroke=PAL["teal"])
    s.edge([prom["b"], mrna["t"]], kind="act", width=3.0, color=PAL["teal"])
    s.rect(1170, 515, 170, 26, "#ffffff", PAL["fb"], rx=6, sw=1.2, dashed=True)
    s.text(1255, 532, "× repression (Module 2)", fs=10.5, anchor="middle", fill=PAL["fb"], bold=True)
    cd = s.node(1150, 745, 190, 64, "Cyclin D1", sub="protein — tracks its mRNA within minutes", fs=17,
                fill=PAL["cd_fill"], stroke=PAL["cd_line"], sw=2.6)
    s.edge([mrna["b"], cd["t"]], kind="act", width=3.0, color=PAL["teal"])

    # ---- condition table ----
    tx, ty = 40, 640
    s.text(tx, ty, "Context inputs per condition", fs=12, bold=True)
    rows = [("", "SHH", "vismodegib", "functional Ptch1", "MYCN copy"),
            ("GNP", "on", "–", "1.0", "1×"),
            ("Ptch1+/−", "on", "–", "0.5", "1×"),
            ("MB", "on", "– / +", "0.1", "2.8×")]
    cols = [0, 90, 140, 230, 350]
    for r, row in enumerate(rows):
        yy = ty + 22 + r * 20
        for c, cell in enumerate(row):
            s.text(tx + cols[c], yy, cell, fs=11, bold=(r == 0 or c == 0))
    s.line([(tx, ty + 28), (tx + 420, ty + 28)], color="#d1d5db")
    s.label(tx, ty + 112, "Whole module is fast (~minutes) except the Gli1 locus memory: it acts as a\n"
                          "quasi-static input that sets the Ccnd1 transcriptional drive.", fs=10, anchor="start")

    s.legend(560, 660, [("act", "activation / conversion"), ("inh", "inhibition"),
                        ("drug", "drug (input)"), ("fb", "slow / optional (dashed)")])
    return s


if __name__ == "__main__":
    out = os.path.join(os.path.dirname(__file__), "module1_mitogenic_drive.svg")
    build().save(out)
    print("wrote", out)
