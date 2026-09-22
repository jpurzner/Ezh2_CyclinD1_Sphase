"""Wireframe — Module 3: the restriction point (two-step Rb, E2F, CKI brakes, Skp2–p27 commitment toggle, growth gates).

Boxes = Heldt core species as rewired in src/build_model_v44_heldt.py (_apply_two_step_rb, SKP2_BLOCK, GROWTH_BLOCK).
No equations.
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from wf import SVG, PAL

W, H = 1400, 800


def build():
    s = SVG(W, H)
    s.text(24, 34, "Module 3 — Restriction point: two-step Rb, E2F and the Skp2–p27 commitment switch", fs=20, bold=True)
    s.text(24, 56, "Cyclin D1–CDK4/6 primes Rb (mono-P, growth-gated); Cyclin E/A–CDK2 commits (hyper-P) and releases E2F; "
                   "INK4 and p27 raise the threshold; Skp2 clears p27", fs=11, fill=PAL["grey"], italic=True)

    B, BL = PAL["cc_fill"] if "cc_fill" in PAL else "#d6eaf8", PAL["blue"]
    R, RL = "#fadbd8", "#a93226"
    RY = 330

    # ---- main axis ----
    cd = s.node(120, RY, 150, 56, "Cyclin D1", sub="from Module 2", fs=15, fill=PAL["cd_fill"], stroke=PAL["cd_line"], sw=2.4)
    cdk46 = s.node(320, RY, 160, 56, "Cyclin D1–CDK4/6", sub="kinase, bounded Vmax", fs=12.5, fill=PAL["cd_fill"], stroke=PAL["cd_line"])
    rb = s.node(520, RY, 90, 48, "Rb", sub="hypo-P", fs=13, fill=R, stroke=RL)
    rbm = s.node(690, RY, 110, 48, "Rb mono-P", sub="E2F still bound", fs=12, fill="#f1948a", stroke=RL)
    prb = s.node(870, RY, 110, 48, "Rb hyper-P", sub="E2F released", fs=12, fill="#cb4335", stroke=RL)
    e2f = s.node(1050, RY, 100, 52, "E2F", sub="free", fs=15, fill="#ffffff", stroke=BL, sw=2)
    s.edge([cd["r"], cdk46["l"]], kind="act", width=3.0, color=PAL["cd_line"])
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
    s.label(605, 372, "Rb·E2F and Rb-mono-P·E2F complexes hold E2F: the mono-P state is the\ntransient-G0 buffer while the cell grows", fs=9.5, anchor="middle")

    # E2F outputs
    ce = s.node(1050, 470, 170, 56, "Cyclin E/A–CDK2", sub="S-phase kinases", fs=12.5, fill=B, stroke=BL)
    s.edge([e2f["b"], ce["t"]], kind="act", width=2.2, color=BL)
    s.label(1062, 410, "transcribes", fs=9.5, color=BL, anchor="start")
    # CDK2 -> hyper step (commits)
    s.edge([ce["l"], (835, 470), (835, RY + 16)], kind="act", width=2.0, color=BL)
    s.label(840, 452, "CDK2 hyper-phosphorylates → COMMITMENT", fs=10, color=BL, bold=True, anchor="start")
    # E2F positive feedback
    s.edge([(e2f["t"][0] + 20, e2f["t"][1]), (e2f["t"][0] + 20, 260), (e2f["t"][0] - 20, 260), (e2f["t"][0] - 20, e2f["t"][1])], kind="act", width=1.2, color=BL)
    s.label(1050, 250, "auto-activation", fs=9, color=BL)
    # E2F -> EZH2 (module 2), E2F -> Skp2
    ezh = s.io(1280, 250, "EZH2 transcription", sub="Module 2 (with Cyclin E/A)", w=200)
    s.edge([e2f["r"], (1180, RY), (1180, 250), ezh["l"]], kind="fb", width=1.6)
    skp2 = s.node(1280, 470, 140, 52, "Skp2", sub="SCF ligase · E2F target", fs=13, fill="#ffffff", stroke=PAL["ink"])
    s.edge([(1180, RY), (1180, 470), skp2["l"]], kind="act", width=1.8, color=BL)
    s.label(1200, 400, "also degraded by\nAPC/C–Cdh1 in G0/G1", fs=9, anchor="start")

    # ---- brakes on CDK4/6 ----
    ink4 = s.node(320, 200, 150, 50, "p16 / p18 (INK4)", sub="competitive · raised in MB", fs=12, fill="#f4ecf7", stroke=PAL["ezh2_line"])
    s.edge([ink4["b"], cdk46["t"]], kind="inh", width=2.0)
    s.label(400, 270, "raise Cyclin D1 needed", fs=9.5, anchor="start")
    palbo = s.pill(140, 200, "palbociclib", w=110)
    s.edge([palbo["r"], (245, 200), cdk46["tl"]], kind="drug", width=2.0)
    s.label(140, 228, "Vmax block (not surmountable)", fs=9, color=PAL["drug"])

    # ---- p27 and Skp2 toggle ----
    p27 = s.node(320, 600, 170, 60, "p27 (CIP/KIP)", sub="reset HIGH at birth (MB higher)\n→ transient G0", fs=13, fill="#fbeee6", stroke="#b9770e", sw=2)
    s.edge([(p27["t"][0] - 12, p27["t"][1]), (cdk46["b"][0] - 12, cdk46["b"][1])], kind="inh", width=1.8)
    s.edge([(cdk46["b"][0] + 12, cdk46["b"][1]), (p27["t"][0] + 12, p27["t"][1])], kind="inh", width=1.8, color=PAL["cd_line"])
    s.label(340, 470, "mutual antagonism:\np27 brakes CDK4/6;\nCDK4/6 clears p27\n(growth-gated)", fs=9, anchor="start")
    s.edge([p27["r"], (1050, 600), ce["b"]], kind="inh", width=2.0)
    s.label(700, 592, "p27 binds and inhibits Cyclin E/A–CDK2", fs=10, bold=True)
    s.edge([skp2["b"], (1280, 680), (320, 680), p27["b"]], kind="inh", width=2.0)
    s.label(800, 672, "Skp2 (SCF) degrades p27 — E2F → Skp2 → ⊣ p27 → CDK2 → more E2F", fs=10, bold=True)
    s.callout(560, 700, 480, ["G0: E2F off → Skp2 low → p27 high → CDK2 off (self-sustaining, p27⁺)",
                              "Commit: mono-P leaks E2F → Skp2 up → p27 cleared → hyper-P → E2F on"],
              title="Skp2–p27 feedforward = bistable commitment switch", color=BL, fs=10)

    # ---- growth / Rb-dilution gates ----
    mass = s.node(560, 470, 170, 56, "cell growth (mass)", sub="grows, halves at division", fs=12.5, fill="#eceff1", stroke=PAL["growth_b"] if "growth_b" in PAL else "#607d8b")
    GC = "#607d8b"
    s.edge([mass["l"], (438, 470), (438, RY + 10)], kind="plain", width=1.4, dashed=True, color=GC)
    s.gate(438, RY + 26, color=GC)
    s.label(456, 425, "commit gate: cell must reach a\ncritical size ≡ Rb diluted by growth", fs=9, color=GC, anchor="start")
    sg = s.io(560, 560, "S-entry size gate → origin firing", sub="Module 4", w=250, h=40, fs=11)
    s.edge([mass["b"], sg["t"]], kind="plain", width=1.4, dashed=True, color=GC)

    s.legend(40, 690, [("act", "activation / conversion"), ("inh", "inhibition"), ("drug", "drug"),
                       ("fb", "hand-off to another module")])
    s.label(1200, 760, "Phase read-outs:  G0 = p27-high / hyper-P-Rb-negative  ·  G1 = hyper-P-Rb⁺, p27⁻, pre-S", fs=9.5, anchor="end")
    return s


if __name__ == "__main__":
    out = os.path.join(os.path.dirname(__file__), "module3_restriction_point.svg")
    build().save(out)
    print("wrote", out)
