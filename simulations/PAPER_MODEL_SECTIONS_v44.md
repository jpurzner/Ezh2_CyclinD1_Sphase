# Computational model — manuscript text for Chahin et al.

Manuscript-ready **Results**, **Methods**, **References**, and **Figure captions** for the
computational model. Written as a standalone description of the model (no reference to any earlier
model version). Figure numbering follows the manuscript: **Fig. 5** (main); **Fig. S8** (architecture
+ validation); **Fig. S9** (Ezh2 feedback, mitogen withdrawal, CDK4/6i). Detailed calibration data,
targets, and the CDK-inhibitor / RNA-seq tables live in `docs/PARAMETERIZATION.md`. Internal author
notes (what is deliberately framed or omitted) are at the end.

---

## 1. RESULTS

Having established that Ezh2 directly represses Cyclin D1 through H3K27me3 (Fig. 4L), we sought to
understand the quantitative consequences of this repression for cell-cycle control and drug response.
To do so, we constructed an ordinary differential equation (ODE) model integrating Hedgehog signaling
and its target gene Mycn, Ezh2, and the mammalian cell-cycle engine (Fig. S8A). Within the model
architecture, Cyclin D1 mRNA transcription serves as the central integration node, receiving three
inputs: basal transcription, Gli-driven transcription (via Shh/Ptch1/Smo/Gli), and Mycn-driven
transcription. All inputs are subject to Ezh2-mediated repression through H3K27me3 at the Cyclin D1
promoter. To capture Ezh2's dependence on S-phase duration, we based the cell-cycle core on the
real-time mammalian G1/S model of Heldt et al. (2018), which represents DNA replication explicitly so
that S-phase has a mechanistic, variable duration. To this core we added four components needed to
reproduce the relevant biology (Methods). A Cyclin B/CDK1 mitotic switch makes entry into and exit from
mitosis abrupt and all-or-none, through feedback between the Cdc25 activator and the Wee1 inhibitor
(Pomerening et al., 2003). An intra-S checkpoint lets ongoing DNA replication restrain mitotic entry,
enforcing the dependence of mitosis on completed replication (Hartwell and Weinert, 1989). A
growth-gated restriction point requires the cell to grow to a critical size before it licenses
commitment and S-phase entry (Ginzberg et al., 2015; Cadart et al., 2018). Finally, a Skp2–p27
feedforward switch, in which E2f-induced Skp2 degrades p27 (Carrano et al., 1999), makes the G0→G1
commitment a sharp, all-or-none decision (Yao et al., 2008). Two further features tie the model to the
tumour genetics: the Hedgehog module includes the canonical Gli→Ptch1 negative feedback, in which Ptch1
is itself a Gli target that re-inhibits Smo (Marigo and Tabin, 1996), intact in GNPs but broken in
Ptch+/− MB, where loss-of-function Ptch1 yields constitutive, mitogen-independent Gli activity (Goodrich
et al., 1997); and the Cyclin D1 commitment threshold is raised in MB by its elevated CDK-inhibitor
tone: the INK4 proteins p16 and p18, which competitively inhibit CDK4/6 (Serrano et al., 1993; Sherr and
Roberts, 1999), and the CIP/KIP proteins p21 and p27, which inhibit CDK2. The elevated p18 contribution
is consistent with the established tumor-suppressor role of Ink4c/p18 in medulloblastoma (Uziel et al.,
2005). Ezh2 transcription is driven by E2f during the S-phase window (Bracken et al., 2003; Pasini et
al., 2004), and because Ezh2 protein is stable and diluted only at division, Ezh2 integrates the
duration of S-phase, closing the feedback loop. The model reproduced 22 of 28 quantitative targets
across 13 conditions, including Shh-dependent GNP cycling, G0 arrest, Hedgehog inhibition (HHi) and Ezh2
inhibition (Ezh2i) responses, and MB drug combinations (Fig. S8B).

**Ezh2 feedback keeps Cyclin D1 coupled to mitogen input, so that withdrawing Hedgehog signaling
collapses the cell cycle.** To quantify the functional impact of Ezh2-mediated Cyclin D1 repression, we
compared simulations with and without the feedback by computationally removing the ability of Ezh2 to
repress Cyclin D1. Removing the feedback raised mean Cyclin D1 mRNA ~3-fold and converted it from a
cell-cycle-coupled oscillation into a constitutively high signal (Fig. S9C,D); this oscillation arises
because Ezh2 protein levels track pRb–E2f activity through the cell cycle, producing a phase-shifted
repression of Cyclin D1 that peaks during S/G2, in anti-phase with Cyclin D1 (Fig. S9E–H). Because the
feedback holds Cyclin D1 bounded and cell-cycle-coupled rather than constitutively high, it keeps
proliferation contingent on continued mitogen (Shh) input — the physiological driver of GNP cycling
(Wechsler-Reya and Scott, 1999) — and concentrates Ezh2's action at the G0/G1 commitment boundary, where
Cyclin D1–CDK4/6-dependent Rb phosphorylation is rate-limiting, consistent with our observation that Ezh2
overexpression specifically traps cells in G0 (Fig. 3F). When we withdrew mitogen by ramping Shh down,
Cyclin D1 fell and, because Ezh2 is cell-cycle-coupled, Ezh2 fell with it and transiently de-repressed
Cyclin D1, sustaining a few additional divisions before Cyclin D1 dropped below threshold and the cell
exited to a stable G0 (Fig. S9K,L). Vismodegib acts as a pharmacological form of this mitogen withdrawal,
removing the Gli-driven input to Cyclin D1. In Ptch+/− MB — where the Gli→Ptch1 feedback is broken and
Cyclin D1 is otherwise constitutive — vismodegib collapsed Cyclin D1 by ~86%, to roughly the level of a
cycling GNP (Fig. 5A); and because MB's elevated CDK-inhibitor tone (INK4 p16/p18 and CIP/KIP p21/p27)
holds the commitment threshold high, this drop crossed the threshold and arrested the cells, whereas a
GNP at the same Cyclin D1 — without the brake — continued to cycle (Fig. 5B).

**Ezh2 inhibition rescues vismodegib-arrested MB by overcoming the CDK-inhibitor threshold, whereas
CDK4/6 inhibition cannot be rescued.** Adding Ezh2 inhibition to vismodegib-treated Ptch+/− MB
de-repressed Cyclin D1 transcription (~2–3-fold), and because the INK4 brake is competitive — more
Cyclin D1 overcomes it — pushed cells back over the commitment threshold and restored proliferation;
across a heterogeneous tumour, the cycling fraction fell from 100% (MB) to ~20% (MB+HHi) and was restored
to ~74% by Ezh2 inhibition, reflecting cell-to-cell variation in the proliferation–quiescence decision
(Spencer et al., 2013; Cappell et al., 2016) (Fig. 5C). Critically, the model predicted that CDK4/6
inhibition — a non-competitive (Vmax) block of Cyclin D–CDK4/6 kinase activity — would arrest MB cells
in a manner that Ezh2 inhibition could not rescue, holding the cycling fraction at 0% with or without
Ezh2i (Fig. S9I,J), because raising Cyclin D1 transcript cannot bypass a kinase-level block. The rescue
of an upstream, Cyclin D1-reversible blockade (vismodegib) but not a downstream, kinase-level one
(CDK4/6i) indicates that Ezh2 inhibition operates specifically through transcriptional de-repression of
Cyclin D1, and motivated the experimental rescue experiments described below.

---

## 2. METHODS — Computational Modeling

### Model construction
We developed an ODE-based model of the Ezh2–Cyclin D1 negative feedback loop embedded within the
mammalian cell cycle, integrating four functional modules: (1) Hedgehog signaling, (2) Mycn regulation,
(3) Ezh2 regulation, and (4) the cell cycle engine. The model was implemented in Antimony format and
simulated using Tellurium with the CVODE integrator in Python 3 (69). The cell cycle engine runs in real
time (minutes), so that the durations of cell-cycle phases — in particular S-phase — emerge from the
underlying replication and growth kinetics rather than being imposed.

### Cell cycle engine
The core was built on the real-time mammalian G1/S model of Heldt et al. (2018) (Heldt 2018), which
represents DNA replication explicitly. An Rb–E2f restriction point governs the G0/G1→S transition:
Cyclin D/CDK4,6, Cyclin E/CDK2, and Cyclin A/CDK2 progressively phosphorylate Rb, and phosphorylated Rb
releases E2f, which transcribes Cyclin E and Cyclin A in a positive-feedback loop. Released E2f licenses
replication origins that fire and synthesize DNA (Dna: 0→1) at a defined fork speed, so that S-phase
duration is set mechanistically by replication kinetics. We extended this core with four features: (i) a
Cyclin B/CDK1 mitotic switch with Cdc25/Wee1 hysteresis (Pomerening 2003); (ii) an intra-S Chk1
checkpoint that holds mitotic entry until replication is complete; (iii) a cell-growth–gated restriction
point, in which cell mass grows continuously and halves at division, and size gates (M_commit, M_size)
couple G0→G1 commitment and S-entry to attainment of a critical size, so that G1 length reflects the time
to grow; and (iv) a Skp2–p27 feedforward switch, in which Skp2 is a dynamic E2f target degraded by
APC/C–Cdh1, and p27/p21 (reset high at each division) is cleared as Cyclin D–CDK4/6 activity rises,
making the G0→G1 commitment bistable (Yao 2008) and producing a transient, p27-high G0 dwell after each
division. The model tracks the Cyclin/CDK complexes (D, E, A, B) and their CDK-inhibitor complexes,
Rb/phospho-Rb, E2f, Skp2, Cdh1, Cdc20, PCNA, the replication machinery (licensed and firing complexes and
Dna), cell mass, and the upstream Hedgehog/Mycn/Ezh2/Cyclin D1 species. The overall period was calibrated
to ~22–23 h to match published GNP proliferation rates (70).

### Hedgehog signaling module
The Hedgehog pathway is an upstream driver of Cyclin D1 transcription. Shh binds Ptch1, relieving its
tonic inhibition of Smoothened (Smo); active Smo converts Gli repressor to Gli activator via a Hill
function (n = 2), and Gli activator drives Gli1 transcription. The module incorporates the canonical
Gli→Ptch1 negative feedback: Ptch1 is itself a Gli target, transcribed at a basal rate plus a
Gli-activator–induced term,

  d(Ptch1_mRNA)/dt = k_basal + k_Gli·Gli_act² / (K_Gli,Ptch² + Gli_act²) − k_deg·Ptch1_mRNA,

and functional Ptch1 protein re-inhibits Smo,

  Smo_active = k_act / (1 + Ptch1_free·f / K_Ptch,Smo) × (1 − HHi),

where *f* is the functional-Ptch1 fraction (a context input: 1.0 in wild-type GNP, reduced to ~0.1 in
Ptch+/− MB to represent loss-of-function Ptch1) and HHi = 0 (untreated) or 1 (Vismodegib). With *f* near
1 the loop is intact and Gli is mitogen-responsive; with *f* ≈ 0.1 the loop is broken and Gli activity is
constitutive and Shh-independent, as in MB (Marigo & Tabin 1996; Goodrich 1997).

### Ezh2 module
Ezh2 mRNA transcription is the sum of a basal term and an E2f-driven term gated on the S-phase window —
proportional to E2f and to the S-phase cyclins (Cyclin E + Cyclin A):

  d(Ezh2_mRNA)/dt = k_basal + k_E2f·E2f·(CyclinE + CyclinA) − k_deg·Ezh2_mRNA.

This makes Ezh2 a cell-cycle-dependent gene that is high in S/G2 (E2f active, replication ongoing) and
low in G0, consistent with the ~2-fold higher Ezh2 in S/G2M vs G0 in the scRNA-seq (Fig. 4A) (Bracken
2003; Pasini 2004). Ezh2 protein is translated from this mRNA, is stable (slow first-order degradation,
half-life ~35 h) and is diluted two-fold at each division; because turnover is slow relative to the
cycle, Ezh2 protein integrates the duration of S-phase. Ezh2 represses all Cyclin D1 transcription
through a repression factor,

  Ezh2_repression = K_Ezh2 / (K_Ezh2 + Ezh2·(1 − Ezh2i)),

with K_Ezh2 = 0.75 and Ezh2i = 0 (untreated) or 1 (full catalytic inhibition, modeling Tazemetostat),
producing ~2–3-fold Cyclin D1 de-repression upon Ezh2 inhibition, matching the experimental data
(Fig. 3C).

### Mycn module
Mycn synthesis has two components: an autonomous (basal) term scaled by an amplification factor (1.0 for
GNPs, 2.86 for MB) and a small Gli-dependent term that is not amplified. This makes MB Mycn predominantly
autonomous and Hedgehog-inhibitor–resistant. Mycn drives Cyclin D1 transcription through a cooperative
Hill function,

  Mycn_contribution = k_Mycn·Mycn^n / (K_Mycn^n + Mycn^n),  n = 3.7,

which supplies the HHi-resistant residual of Cyclin D1: at GNP Mycn levels the contribution is small,
whereas in MB it sustains the ~14% of Cyclin D1 that persists after Gli ablation by Vismodegib.

### Cyclin D1 mRNA integration
Cyclin D1 mRNA transcription integrates basal, Gli-driven, and Mycn-driven inputs, each multiplied by the
Ezh2 repression factor:

  d(CyclinD1_mRNA)/dt = (k_basal + k_Gli·f(Gli_act, Gli1, Gli_rep) + k_Mycn·g(Mycn))·Ezh2_repression
                        − k_deg·CyclinD1_mRNA,

where *f*() is Gli-activator-driven, Gli-repressor-inhibited transcription and *g*() is the Mycn Hill
term. Calibration to the RNA-seq (below) makes the Gli-driven term dominant — Gli1 and Cyclin D1 rise
~7-fold in parallel from GNP to MB, identifying Cyclin D1 as predominantly Gli-driven, with Mycn supplying
the HHi-resistant residual. Cyclin D1 protein is translated from this mRNA and feeds the cell cycle
through Cyclin D–CDK4/6-mediated Rb phosphorylation.

### Commitment threshold: INK4 and CIP/KIP CDK inhibitors
Cyclin D–CDK4/6 phosphorylation of Rb saturates with Cyclin D1 (CDK4/6 kinase activity is bounded) and is
competitively inhibited by the INK4 proteins:

  rate = kPhRbCd · Cd / (K_CdRb·(1 + p16 + p18) + Cd).

p16 (Cdkn2a) and p18 (Cdkn2c) are INK4-family CDK4/6 inhibitors entered as context inputs (Serrano 1993;
Sherr & Roberts 1999): p16 is silent in GNPs (0) and elevated in MB (0.88), while p18 is expressed in
GNPs (0.4) and further elevated in MB (1.2), consistent with the tumor-suppressor role of Ink4c/p18 in
medulloblastoma (Uziel 2005). Raising the INK4 tone increases the Cyclin D1 level required to
phosphorylate Rb and commit, but because the inhibition is competitive, additional Cyclin D1 overcomes
it. In parallel, the CIP/KIP inhibitors p21/p27 (the model's p27 species) inhibit Cyclin E/A–CDK2; their
synthesis rate is increased ~2-fold in MB relative to GNP, raising the threshold further and lengthening
the p27-high transient G0. Together the elevated INK4 + CIP/KIP tone in MB places the commitment
threshold above the Cyclin D1 level that remains after Vismodegib.

### Drug implementations
HHi (Vismodegib) blocks Smo activation (HHi = 1, above). Ezh2i (Tazemetostat) sets Ezh2i = 1, abolishing
Ezh2 repression of Cyclin D1. CDK4/6 inhibition (Palbociclib) sets the Cyclin D–CDK4/6-mediated
Rb-phosphorylation rate (kPhRbCd) to zero — a non-competitive (Vmax) block that, unlike the competitive
INK4 brake, cannot be bypassed by raising Cyclin D1, and that leaves Cyclin D1 transcription unaffected.
Hydroxyurea (HU) reduces replication fork speed (the DNA-synthesis rate is scaled by a Hill function of
HU dose), lengthening S-phase in a concentration-dependent manner while still permitting mitosis once
replication completes.

### Parameter calibration
Parameters were calibrated against quantitative RNA-seq from P7 wild-type GNPs and Ptch+/− MB, with and
without HHi (Vismodegib). The following experimental targets were used:

| Target | Experimental value | Source |
|---|---|---|
| Cyclin D1: GNP+HHi / GNP | 0.16 (84% reduction) | RNA-seq |
| Cyclin D1: MB+HHi / MB | 0.14 (86% reduction) | RNA-seq |
| Cyclin D1: MB / GNP | ~7× | RNA-seq |
| Gli1: MB / GNP | ~7× | RNA-seq |
| Gli1: GNP+HHi reduction | >99% | RNA-seq |
| Mycn: GNP+HHi / GNP | 0.78 (22% reduction) | RNA-seq |
| Mycn: MB+HHi / MB | 0.86 (14% reduction) | RNA-seq |
| Mycn: MB / GNP | ~2.8× | RNA-seq |
| Ezh2: G0 / cycling ratio | ~0.6 | scRNA-seq (Fig. 4A) |
| Ezh2i: Cyclin D1 fold change | ~2× | qPCR (Fig. 3C) |
| Cell cycle period (GNP) | ~22–23 h | Published GNP data (70) |
| Vismodegib MB proliferation (pRb⁺) | ~1/4 of DMSO; restored to ~3/4 by Ezh2i | Experimental data (Fig. 5) |

The Hedgehog/Mycn→Cyclin D1 and Ezh2 parameters were fit to the RNA-seq folds by a constrained parameter
search that rejected any parameter set failing to integrate numerically or breaking the figure-critical
cycling behaviors (GNP cycles only with Shh; GNP−Shh, GNP+HHi and serum-starved GNPs arrest; MB cycles;
MB+CDK4/6i arrests and is not rescued by Ezh2i), minimizing a weighted error on the expression ratios
among survivors. The CDK-inhibitor brake (p16, p18, and the p27 synthesis rate) and the cell-to-cell
Cyclin D1/p27 heterogeneity were then set to reproduce the experimental proliferation (pRb⁺) response —
Vismodegib reducing the MB cycling fraction to ~1/4 of control and Ezh2i restoring it to ~3/4. Key
calibrated values: Mycn Hill n = 3.7, k_Mycn = 35, K_Mycn = 1.7, Mycn amplification = 2.86; K_Ezh2 = 0.75;
K_CdRb = 0.357 with p16/p18 as above; MB p27 synthesis = 2× the GNP rate. (Full target compendium and the
CDK-inhibitor RNA-seq tables: `docs/PARAMETERIZATION.md`.)

### Simulation conditions
Simulations used Tellurium with the CVODE integrator (absolute tolerance 1e-9, relative tolerance 1e-6).
MB was modeled as reduced functional Ptch1 with Mycn amplification and Shh present (Ptch+/− loses
functional Ptch1 but the ligand-independent pathway is constitutively active), rather than as ligand
removal:

| Condition | SHH | Ptch1 (f) | HHi | Ezh2i | Mycn amp. | INK4 (p16/p18) | p27 synth. | CDK4/6i | HU |
|---|---|---|---|---|---|---|---|---|---|
| GNP + SHH | 0.5 | 1.0 | 0 | 0 | 1.0 | 0 / 0.4 | 1× | No | No |
| GNP − SHH | 0 | 1.0 | 0 | 0 | 1.0 | 0 / 0.4 | 1× | No | No |
| GNP + HHi | 0.5 | 1.0 | 1 | 0 | 1.0 | 0 / 0.4 | 1× | No | No |
| GNP + Ezh2i | 0.5 | 1.0 | 0 | 1 | 1.0 | 0 / 0.4 | 1× | No | No |
| GNP + HHi + Ezh2i | 0.5 | 1.0 | 1 | 1 | 1.0 | 0 / 0.4 | 1× | No | No |
| GNP serum-starved | 0.5 | 1.0 | 0 | 0 | 1.0 | 0 / 0.4 | 1× | No | No |
| MB (Ptch+/−) | 0.5 | 0.1 | 0 | 0 | 2.8 | 0.88 / 1.2 | 2× | No | No |
| MB + HHi | 0.5 | 0.1 | 1 | 0 | 2.8 | 0.88 / 1.2 | 2× | No | No |
| MB + Ezh2i | 0.5 | 0.1 | 0 | 1 | 2.8 | 0.88 / 1.2 | 2× | No | No |
| MB + HHi + Ezh2i | 0.5 | 0.1 | 1 | 1 | 2.8 | 0.88 / 1.2 | 2× | No | No |
| MB + CDK4/6i | 0.5 | 0.1 | 0 | 0 | 2.8 | 0.88 / 1.2 | 2× | Yes | No |
| MB + CDK4/6i + Ezh2i | 0.5 | 0.1 | 0 | 1 | 2.8 | 0.88 / 1.2 | 2× | Yes | No |
| MB + HU | 0.5 | 0.1 | 0 | 0 | 2.8 | 0.88 / 1.2 | 2× | No | Yes |

Validation simulations were run for ~200 h; feedback-impact and drug-rescue simulations for 168 h to
ensure stable oscillatory behavior. Cell divisions were counted as peaks in Cyclin B (Cyclin B/CDK1, MPF)
concentration (prominence > 0.15, minimum spacing 200 min) after an initial settling period; period was
the mean interval between consecutive peaks; steady-state species concentrations were the mean over the
settled window. Cells that failed to integrate at the commitment bifurcation under high-Cyclin D1
conditions were re-run with progressively looser tolerances and a shorter horizon, and excluded only if
integration failed at all settings.

### Ezh2 feedback analysis
The impact of Ezh2-mediated Cyclin D1 repression was assessed by comparing simulations with intact
feedback (K_Ezh2 = 0.75) versus computationally ablated feedback (K_Ezh2 = 1×10⁶, making the repression
factor ≈ 1.0 regardless of Ezh2 level), in both the GNP + SHH and MB contexts. Phase relationships between
Cyclin D1 mRNA, Ezh2 protein, and Cyclin B protein were visualized by normalizing each species to [0, 1]
over the last 70 h of the 168-h simulation. Mitogen withdrawal was modeled by ramping Shh down in a
cycling GNP and recording Cyclin D1, Ezh2, and divisions along the falling input.

### Population (cycling-fraction) analysis
The fraction of cycling cells under each MB drug condition (Fig. 5C) was computed as an ensemble of 120 MB
cells carrying log-normal cell-to-cell variation in two mother-set quantities — the Cyclin D1 setpoint
(translation rate) and the inherited p27 level (reset at division). Each cell was simulated independently
and scored as cycling if it produced ≥2 Cyclin B peaks after settling; the cycling fraction is the
percentage of integrable cells that cycled.

---

## 3. REFERENCES (model features; verify volume/page and renumber to the manuscript scheme)

1. **Heldt FS, Barr AR, Cooper S, Bakal C, Novák B.** A comprehensive model for the proliferation–quiescence decision in response to endogenous DNA damage in human cells. *Proc Natl Acad Sci USA*. 2018;115(10):2532–2537. — explicit-replication cell-cycle core.
2. **Bracken AP, Pasini D, Capra M, Prosperini E, Colli E, Helin K.** EZH2 is downstream of the pRB-E2F pathway, essential for proliferation and amplified in cancer. *EMBO J*. 2003;22(20):5323–5335. — Ezh2 is an E2F target.
3. **Pasini D, Bracken AP, Helin K.** Polycomb group proteins in cell cycle progression and cancer. *Cell Cycle*. 2004;3(4):396–400.
4. **Marigo V, Tabin CJ.** Regulation of *patched* by sonic hedgehog in the developing neural tube. *Proc Natl Acad Sci USA*. 1996;93(18):9346–9351. — Ptch1 is a Gli-induced feedback target.
5. **Goodrich LV, Milenković L, Higgins KM, Scott MP.** Altered neural cell fates and medulloblastoma in mouse *patched* mutants. *Science*. 1997;277(5329):1109–1113.
6. **Serrano M, Hannon GJ, Beach D.** A new regulatory motif in cell-cycle control causing specific inhibition of cyclin D/CDK4. *Nature*. 1993;366(6456):704–707. — INK4/p16 competitive CDK4/6 inhibition.
7. **Sherr CJ, Roberts JM.** CDK inhibitors: positive and negative regulators of G1-phase progression. *Genes Dev*. 1999;13(12):1501–1512. — INK4 (competitive) vs CIP/KIP (stoichiometric).
8. **Uziel T, Zindy F, Xie S, et al.** The tumor suppressors Ink4c and p53 collaborate independently with Patched to suppress medulloblastoma formation. *Genes Dev*. 2005;19(22):2656–2667. — p18^Ink4c^ in medulloblastoma.
9. **Pomerening JR, Sontag ED, Ferrell JE Jr.** Building a cell cycle oscillator: hysteresis and bistability in the activation of Cdc2. *Nat Cell Biol*. 2003;5(4):346–351. — Cyclin B/CDK1 mitotic switch.
10. **Yao G, Lee TJ, Mori S, Nevins JR, You L.** A bistable Rb-E2F switch underlies a Ki67 cell cycle entry threshold. *Nat Cell Biol*. 2008;10(4):476–482. — Rb–E2F restriction-point bistability.
11. **Wechsler-Reya RJ, Scott MP.** Control of neuronal precursor proliferation in the cerebellum by Sonic Hedgehog. *Neuron*. 1999;22(1):103–114. — Shh is the GNP mitogen.
12. **Spencer SL, Cappell SD, Tsai FC, Overton KW, Wang CL, Meyer T.** The proliferation-quiescence decision is controlled by a bifurcation in CDK2 activity at mitotic exit. *Cell*. 2013;155(2):369–383.
13. **Cappell SD, Chung M, Jaimovich A, Spencer SL, Meyer T.** Irreversible APC^Cdh1^ inactivation underlies the point of no return for cell-cycle entry. *Cell*. 2016;166(1):167–180.
14. **Hartwell LH, Weinert TA.** Checkpoints: controls that ensure the order of cell cycle events. *Science*. 1989;246(4930):629–634. — dependence of mitosis on completed replication (intra-S checkpoint).
15. **Ginzberg MB, Kafri R, Kirschner M.** On being the right (cell) size. *Science*. 2015;348(6236):1245075. — growth-gated cell-size control of cell-cycle progression.
16. **Cadart C, Monnier S, Grilli J, et al.** Size control in mammalian cells involves modulation of both growth rate and cell cycle duration. *Nat Commun*. 2018;9(1):3275. — growth-coupled timing of cell-cycle phases.
17. **Carrano AC, Eytan E, Hershko A, Pagano M.** SKP2 is required for ubiquitin-mediated degradation of the CDK inhibitor p27. *Nat Cell Biol*. 1999;1(4):193–199. — E2f-induced Skp2 degrades p27 (Skp2–p27 feedforward switch).

Retained: **Gérard C, Goldbeter A.** *Proc Natl Acad Sci USA*. 2009;106(51):21643–21648 (ref 30 — cyclin/CDK-oscillator lineage). Tellurium (ref 69) and the GNP-period source (ref 70) as in the manuscript.

---

## 4. FIGURE CAPTIONS

**Figure 5. The MB Cyclin D1 axis and the Ezh2-inhibitor rescue of vismodegib (model).**
(A) Calibrated expression: Gli1 and Cyclin D1 rise ~7-fold in lockstep from P7 GNP to Ptch+/− MB, and
vismodegib (HHi) collapses MB Cyclin D1 by ~86%, to roughly the level of a cycling GNP (model vs RNA-seq).
(B) Predicted Cyclin B/CDK1 trajectories: Ptch+/− MB cycles; MB+HHi arrests because the vismodegib-reduced
Cyclin D1 falls below the INK4/CIP-KIP-raised commitment threshold, whereas a GNP at the same Cyclin D1
(no brake) continues to cycle. (C) Cycling fraction of a heterogeneous MB ensemble (N = 120 cells, with
cell-to-cell Cyclin D1/p27 variation): 100% (MB) → ~20% (MB+HHi) → ~74% (MB+HHi+Ezh2i); CDK4/6 inhibition
holds the fraction at 0% with or without Ezh2i.

**Supplementary Figure 8. Model architecture and validation.**
(A) Architecture of the ODE model: Hedgehog (with the Gli→Ptch1 negative feedback) and Mycn → Cyclin D1;
the explicit-replication cell-cycle engine; the growth-gated restriction point and Skp2–p27 feedforward;
the INK4 (p16/p18) competitive and CIP/KIP (p21/p27) CDK-inhibitor commitment brakes; the Cyclin B/CDK1
mitotic switch and intra-S Chk1 checkpoint; and the Ezh2 layer that integrates S-phase duration and
represses Cyclin D1. (B) Validation: simulated vs measured ratios and behaviors across the 13 conditions
(22 of 28 targets passed) — Cyclin D1, Mycn, Gli1 and Ezh2 fold changes, cell-cycle period, GNP cycling
and arrest, and the MB drug combinations.

**Supplementary Figure 9. Ezh2 feedback, mitogen withdrawal, and the CDK4/6i contrast.**
(A,B) Cyclin B traces with and without the Ezh2→Cyclin D1 feedback, in GNP+SHH and Ptch+/− MB. (C,D)
Cyclin D1 mRNA with and without the feedback: removing it raises mean Cyclin D1 ~3-fold and abolishes its
oscillation, in GNP and MB. (E–H) Ezh2 protein and the normalized phase relationship of Cyclin D1 mRNA,
Ezh2 protein and Cyclin B over the last 70 h (Ezh2 peaks in S/G2, in anti-phase with Cyclin D1), in GNP
and MB. (I,J) MB + CDK4/6i and MB + CDK4/6i + Ezh2i: the non-competitive (Vmax) kinase block arrests MB
and cannot be rescued by Ezh2 inhibition. (K,L) Mitogen withdrawal: ramping Shh down in a cycling GNP —
Cyclin D1 and Ezh2 fall together, sustaining a few additional divisions before the cell exits to a stable
G0.

---

## 5. Internal notes for the authors (not for the manuscript)

- The cell-cycle core is the explicit-replication Heldt (2018) model run in real time (minutes); this is
  why S-phase duration is mechanistic and Ezh2 (an S-phase-window E2F target) can integrate it.
- The Ezh2→Cyclin D1 feedback sets Cyclin D1 *amplitude/oscillation* and keeps the cell mitogen-gated; it
  does **not** appreciably change the cell-cycle period (phase durations are set by replication and
  growth), so no period-lengthening claim is made and the period panels are not cited in the text.
- MB is modeled as **reduced functional Ptch1 (f = 0.1) + Mycn amplification with Shh present** (broken
  Gli→Ptch1 feedback → constitutive Gli), not as ligand removal.
- The vismodegib-treated MB Cyclin D1 is **~0.14 of MB** (86% drop, ≈ a cycling GNP), per the RNA-seq; the
  arrest comes from the MB CDK-inhibitor threshold, and the Ezh2i rescue is a **fractional/population**
  effect (100 → ~20 → ~74%).
- Validation 22/28: the 6 misses are the absolute Ezh2 MB/GNP ratio and the S/G1 phase-proportion split,
  documented soft/architectural limits (`docs/PARAMETERIZATION.md`).
- New S9K,L (mitogen withdrawal) uses the `mitogen_rampdown` / `differentiation_collapse` simulations;
  S9A,B (Cyclin B ± feedback) are present but uncited — keep or drop per layout.
- A senescence-vs-transient-G0 scRNA-seq analysis (spec: `docs/scRNAseq_senescence_vs_transientG0_spec.md`)
  found no senescence signature in the cycling compartment and no p16/p21 association of senescent cells —
  supporting the "reversible quiescence, not senescence" framing. Held for reviewer rebuttal, not in text.
