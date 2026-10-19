"""Wireframe — Module 3: the restriction point (two-step Rb, E2F, dynamic CDKIs, p27 buffer, CDK6 INK4 sink,
Skp2–p27 commitment toggle, CDK4/6-alone escape, growth gates).

Boxes = Heldt core species as rewired in the 2026-10 default of src/build_model_v44_heldt.py (_apply_two_step_rb,
SKP2_BLOCK, GROWTH_BLOCK, with_cdki_species + with_cdki_mrna, with_p27_optionB, with_cdk6_gli, with_cd_hyper_escape).
No equations.
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from wf import SVG, PAL

W, H = 1400, 800


def build():
    s = SVG(W, H)
    s.text(24, 34, "Module 3 — Restriction point: two-step Rb, E2F, the CDKI brakes and the Skp2–p27 commitment switch", fs=20, bold=True)
    s.text(24, 56, "Cyclin D1/D2–CDK4/6 primes Rb (mono-P); Cyclin E/A–CDK2 commits (hyper-P) and releases E2F; INK4 (p18/p19) and p27 raise "
                   "the threshold; free CDK6 sinks INK4; Skp2 and KPC clear p27; weak CDK4/6 hyper-P lets high-drive cells escape G0",
           fs=11, fill=PAL["grey"], italic=True)

    B, BL = "#d6eaf8", PAL["blue"]
    R, RL = "#fadbd8", "#a93226"
    GC = "#607d8b"
    RY = 330

    # ---- main axis ----
    cd = s.node(120, RY, 150, 56, "Cyclin D1", sub="from Module 2", fs=15, fill=PAL["cd_fill"], stroke=PAL["cd_line"], sw=2.4)
    cdk46 = s.node(320, RY, 160, 56, "Cyclin D–CDK4/6", sub="kinase, bounded Vmax", fs=12.5, fill=PAL["cd_fill"], stroke=PAL["cd_line"])
    rb = s.node(520, RY, 90, 48, "Rb", sub="hypo-P", fs=13, fill=R, stroke=RL)
    rbm = s.node(690, RY, 110, 48, "Rb mono-P", sub="E2F still bound", fs=12, fill="#f1948a", stroke=RL)
    prb = s.node(870, RY, 110, 48, "Rb hyper-P", sub="E2F released", fs=12, fill="#cb4335", stroke=RL)
    e2f = s.node(1050, RY, 100, 52, "E2F", sub="free", fs=15, fill="#ffffff", stroke=BL, sw=2)
    s.edge([cd["r"], cdk46["l"]], kind="act", width=3.0, color=PAL["cd_line"])
    cd2 = s.io(120, 425, "Cyclin D2 (Module 1)", sub="Hh-buffered · low weight", w=160, h=44, fs=11.5)
    s.edge([cd2["r"], (265, 425), (265, cdk46["b"][1])], kind="act", width=1.4, color=PAL["cd_line"])
    s.edge([cdk46["r"], rb["l"]], kind="act", width=2.2, color=PAL["cd_line"])
    s.label(438, 300, "phosphorylates", fs=9.5, color=PAL["cd_line"])
    s.edge([rb["r"], rbm["l"]], kind="act", width=2.4, color=RL)
    s.edge([(rbm["l"][0], RY + 14), (rb["r"][0], RY + 14)], kind="plain", width=1.0)
    s.label(600, 300, "STEP 1 · mono-P", fs=10, color=RL, bold=True)
    s.edge([rbm["r"], prb["l"]], kind="act", width=2.4, color=RL)
    s.edge([(prb["l"][0], RY + 14), (rbm["r"][0], RY + 14)], kind="plain", width=1.0)
    s.label(780, 300, "STEP 2 · hyper-P", fs=10, color=RL, bold=True)
    s.edge([prb["r"], e2f["l"]], kind="act", width=2.4, color=BL)
    s.label(965, 300, "releases", fs=9.5, color=BL)
    s.label(605, 372, "Rb·E2F and Rb-mono-P·E2F complexes hold E2F: the mono-P state is the\ntransient-G0 buffer; phospho-Rb(Ser807/811)⁺ = the hyper-P species", fs=9.5, anchor="middle")

    # CDK4/6-alone escape: weak direct hyper-P
    s.edge([(cdk46["tr"][0] - 10, cdk46["tr"][1]), (390, 278), (835, 278), (835, RY - 14)], kind="act", width=1.2, color=PAL["cd_line"], dashed=True)
    s.label(640, 272, "weak DIRECT hyper-P by CDK4/6 → high-drive (MB) cells escape the G0 lock (reversible G0)", fs=9, color=PAL["cd_line"])

    # E2F outputs
    ce = s.node(1050, 470, 170, 56, "Cyclin E/A–CDK2", sub="S-phase kinases", fs=12.5, fill=B, stroke=BL)
    s.edge([e2f["b"], ce["t"]], kind="act", width=2.2, color=BL)
    s.label(1062, 410, "transcribes", fs=9.5, color=BL, anchor="start")
    s.edge([ce["l"], (835, 470), (835, RY + 16)], kind="act", width=2.0, color=BL)
    s.label(840, 430, "CDK2 hyper-phosphorylates → COMMITMENT", fs=10, color=BL, bold=True, anchor="start")
    s.edge([(e2f["t"][0] + 20, e2f["t"][1]), (e2f["t"][0] + 20, 260), (e2f["t"][0] - 20, 260), (e2f["t"][0] - 20, e2f["t"][1])], kind="act", width=1.2, color=BL)
    s.label(1050, 250, "auto-activation", fs=9, color=BL)
    ezh = s.io(1280, 250, "EZH2 transcription", sub="Module 2 (with Cyclin E/A + Cyclin D1 dose)", w=220)
    s.edge([e2f["r"], (1180, RY), (1180, 250), ezh["l"]], kind="fb", width=1.6)
    skp2 = s.node(1280, 470, 160, 56, "Skp2", sub="SCF ligase · E2F target\n+ Cyclin D1 dose (MB ×2.3)", fs=13, fill="#ffffff", stroke=PAL["ink"], sub_fs=9)
    s.edge([(1180, RY), (1180, 470), skp2["l"]], kind="act", width=1.8, color=BL)
    s.label(1192, 400, "Skp2 degraded by\nAPC/C–Cdh1 in G0/G1", fs=9, anchor="start")

    # ---- brakes on CDK4/6: INK4 (p18/p19), CDK6 sink, palbociclib ----
    ink4 = s.node(320, 200, 160, 54, "p18 / p19 (INK4)", sub="dynamic proteins (mRNA → protein)\nMB ×3.7 / ×1.6", fs=12, fill="#f4ecf7", stroke=PAL["ezh2_line"], sub_fs=9)
    s.edge([ink4["b"], cdk46["t"]], kind="inh", width=2.0)
    s.label(232, 262, "competitive: raise the\nCyclin D needed", fs=9, anchor="end")
    palbo = s.pill(140, 200, "palbociclib", w=110)
    s.edge([palbo["r"], (245, 200), cdk46["tl"]], kind="drug", width=2.0)
    s.label(140, 228, "Vmax block (not surmountable)", fs=9, color=PAL["drug"])
    cdk6 = s.node(560, 200, 160, 54, "CDK6 (free pool)", sub="Module 2 · Gli-driven, PRC2-governed\nMB ≈ 5× GNP", fs=12, fill=B, stroke=BL, sub_fs=9)
    s.edge([cdk6["l"], ink4["r"]], kind="inh", width=2.0, color=BL)
    s.label(440, 192, "sinks p18/p19", fs=9, color=BL, bold=True)
    s.label(560, 244, "MB's huge CDK6 overwhelms its ×3 INK4", fs=8.5, color=BL)

    # ---- p27, its Cyclin D buffer, and the Skp2 toggle ----
    p27 = s.node(320, 600, 170, 60, "p27 (CIP/KIP)", sub="functional nuclear pool · born HIGH\n(partial carry-over) · mRNA → protein", fs=13, fill="#fbeee6", stroke="#b9770e", sw=2, sub_fs=9)
    buf = s.node(320, 470, 180, 50, "p27 · Cyclin D1–CDK4/6", sub="INACTIVE buffer (consumes Cyclin D1)", fs=11.5, fill="#fdf2e9", stroke="#b9770e", sub_fs=9)
    s.edge([p27["t"], buf["b"]], kind="act", width=1.6, color="#b9770e")
    s.label(332, 552, "binds", fs=9, anchor="start")
    s.edge([buf["t"], cdk46["b"]], kind="inh", width=2.0, color="#b9770e")
    s.label(332, 420, "removes Cyclin D1 from the drive → commitment set by\nthe Cyclin D1 : p27 ratio (Fan–Meyer); INK4 releases p27", fs=9, anchor="start")
    s.edge([p27["r"], (1050, 600), ce["b"]], kind="inh", width=2.0)
    s.label(700, 592, "free p27 binds and inhibits Cyclin E/A–CDK2   (p21: low · p57: off)", fs=10, bold=True)
    s.edge([skp2["b"], (1280, 680), (320, 680), p27["b"]], kind="inh", width=2.0)
    s.label(800, 672, "SCF–Skp2 degrades p27 (plus constant KPC clearance of the free pool) — E2F → Skp2 → ⊣ p27 → CDK2 → more E2F", fs=10, bold=True)
    s.callout(560, 700, 560, ["G0: E2F off → Skp2 low → p27 high → CDK2 off (self-sustaining; pRb-hypo, p27⁺)",
                              "Commit: mono-P leaks E2F → Skp2 up → p27 cleared → hyper-P → E2F on; threshold set by",
                              "Cyclin D1/D2 vs INK4 + p27 (growth-independent); high drive can also escape via CDK4/6 alone"],
              title="Skp2–p27 feedforward = bistable commitment switch", color=BL, fs=10)

    # ---- growth / size gate ----
    mass = s.node(560, 470, 180, 56, "cell growth (mass)", sub="rate slowed by INK4 (MB) and by\nlow Cyclin D1 · halves at division", fs=12.5, fill="#eceff1", stroke=GC, sub_fs=9)
    sg = s.io(560, 560, "S-entry size gate → origin firing", sub="Module 4 (≡ Rb dilution by growth)", w=250, h=40, fs=11)
    s.edge([mass["b"], sg["t"]], kind="plain", width=1.4, dashed=True, color=GC)
    s.label(680, 430, "GNP ≈ 17 h cycle · MB ≈ 24–26 h\n(emerges from INK4 tone)", fs=9, color=GC, anchor="start")

    s.legend(40, 690, [("act", "activation / conversion"), ("inh", "inhibition"), ("drug", "drug"),
                       ("fb", "hand-off to another module")])
    s.label(1200, 795, "Phase read-outs:  G0 = p27-high / hyper-P-Rb-negative  ·  G1 = hyper-P-Rb⁺, p27⁻, pre-S  ·  total p27 (IHC) ≠ functional p27", fs=9.5, anchor="end")
    return s


if __name__ == "__main__":
    out = os.path.join(os.path.dirname(__file__), "module3_restriction_point.svg")
    build().save(out)
    print("wrote", out)
