# Computational model — Results and Methods (with model-design rationale)

Self-contained **Results** and **Methods** for the manuscript, written standalone (no reference to any
earlier model version). This version makes explicit **why** each component was added on top of the core
cell-cycle model — the rationale that was implicit in the prior drafts. Figure numbering follows the
manuscript (Fig. 5; Fig. S8 architecture/validation; Fig. S9 feedback / mitogen withdrawal / CDK4/6i).
Calibration data and CDK-inhibitor / RNA-seq tables: `docs/PARAMETERIZATION.md`.

---

## RESULTS

Having established that Ezh2 directly represses Cyclin D1 through H3K27me3 (Fig. 4L), we sought to
understand the quantitative consequences of this repression for cell-cycle control and drug response.
To do so, we constructed an ordinary differential equation (ODE) model integrating Hedgehog signaling,
Mycn, Ezh2, and the mammalian cell cycle engine (Fig. S8A). Within the model architecture, Cyclin D1
mRNA transcription serves as the central integration node, receiving three inputs: basal transcription,
Gli-driven transcription (via Shh/Ptch1/Smo/Gli), and Mycn-driven transcription. All inputs are subject
to Ezh2-mediated repression through H3K27me3 at the Cyclin D1 promoter. Because our central observation
is that Ezh2 accumulates in proportion to the time a cell spends in S-phase (Fig. 4G) — a quantity that
is fixed in a phenomenological oscillator — we built the cell cycle core on the real-time mammalian G1/S
model of Heldt et al. (2018), which represents DNA replication explicitly, so that S-phase has a
mechanistic, concentration-dependent duration; we then extended this core so that it resolves mitosis,
gates mitosis on completed replication, times cell-cycle phases by cell growth, and makes the G0/G1
commitment switch-like (Methods). Two further features tie the model to the tumour genetics: the
Hedgehog module includes the canonical Gli→Ptch1 negative feedback, in which Ptch1 is itself a Gli
target that re-inhibits Smo (Marigo and Tabin, 1996), intact in GNPs but broken in Ptch+/− MB, where
loss-of-function Ptch1 yields constitutive, mitogen-independent Gli activity (Goodrich et al., 1997); and
the Cyclin D1 commitment threshold is raised in MB by its elevated CDK-inhibitor tone — the INK4 proteins
p16 and p18, which competitively inhibit CDK4/6 (Serrano et al., 1993; Sherr and Roberts, 1999), and the
CIP/KIP proteins p21 and p27, which inhibit CDK2 — consistent with the established tumor-suppressor role
of Ink4c/p18 in medulloblastoma (Uziel et al., 2005). Ezh2 transcription is driven by E2f during the
S-phase window (Bracken et al., 2003; Pasini et al., 2004), and because Ezh2 protein is stable and
diluted only at division, Ezh2 integrates the duration of S-phase, closing the feedback loop. The model
reproduced 22 of 28 quantitative targets across 13 conditions, including Shh-dependent GNP cycling, G0
arrest, Hedgehog inhibition (HHi) and Ezh2 inhibition (Ezh2i) responses, and MB drug combinations
(Fig. S8B).

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

## METHODS — Computational Modeling

### Model construction
We developed an ODE-based model of the Ezh2–Cyclin D1 negative feedback loop embedded within the
mammalian cell cycle, integrating four functional modules: (1) Hedgehog signaling, (2) Mycn regulation,
(3) Ezh2 regulation, and (4) the cell cycle engine. The model was implemented in Antimony format and
simulated using Tellurium with the CVODE integrator in Python 3 (69). The cell cycle engine runs in real
time (minutes), so that the durations of cell-cycle phases — in particular S-phase — emerge from the
underlying replication and growth kinetics rather than being imposed.

### Cell cycle engine — core and the rationale for each extension
The central hypothesis we set out to model is that Ezh2 accumulates in proportion to the time a cell
spends in S-phase (the hydroxyurea result, Fig. 4G). This cannot be represented in a phenomenological
relaxation-oscillator, in which S-phase duration is a fixed structural property; we therefore required a
core in which S-phase duration is a genuine, variable quantity. We built on the real-time mammalian G1/S
model of Heldt et al. (2018) (Heldt 2018), which represents DNA replication explicitly: an Rb–E2f
restriction point governs the G0/G1→S transition, with Cyclin D/CDK4,6, Cyclin E/CDK2 and Cyclin A/CDK2
progressively phosphorylating Rb; phosphorylated Rb releases E2f, which transcribes Cyclin E and Cyclin A
in a positive-feedback loop and licenses replication origins that fire and synthesize DNA (Dna: 0→1) at a
defined fork speed — so that the duration of S-phase is set mechanistically by replication kinetics.

The Heldt G1/S core does not, however, resolve mitosis, gate mitosis on completed replication, set phase
durations by cell growth, or make the commitment decision switch-like — each of which our data require —
so we extended it with four components, each added for a specific reason:

- **Cyclin B/CDK1 mitotic switch (Cdc25/Wee1 hysteresis; Pomerening 2003).** *Why:* the G1/S core does not
  produce mitosis, yet our experimental readout of proliferation is the rate of division. We added an
  explicit Cyclin B/CDK1 mitotic oscillator so that (i) the model completes a full division cycle and
  divisions can be scored as Cyclin B peaks — the basis for the period and division-count metrics — and
  (ii) the discrete events of division (cell-mass halving and resetting p27 to a high level) can occur,
  which in turn generate the post-division transient G0.

- **Intra-S Chk1 checkpoint.** *Why:* to model hydroxyurea, which slows replication forks. The checkpoint
  holds mitotic entry until DNA replication completes, so that lengthening S-phase lengthens the cell
  cycle (and increases Ezh2 accumulation, Fig. 4G) rather than allowing premature mitosis with
  unreplicated DNA. Without this coupling, slowing replication would not translate into a longer effective
  cycle or more Ezh2.

- **Cell-growth–gated restriction point.** *Why:* experimentally, cell-cycle phase durations are governed
  largely by cell growth and size, not by Cyclin D1 abundance. We let cell mass grow continuously, halve
  at division, and gate both G0→G1 commitment (M_commit) and S-entry (M_size) on attainment of a critical
  size. This makes Cyclin D1 act at the commitment *decision* (a yes/no threshold) rather than set cycle
  *speed*, which is required to reproduce both that Ezh2 (through Cyclin D1) acts specifically at the
  G0/G1 boundary where it traps cells in G0 (Fig. 3F), and that the cell-cycle period is largely
  insensitive to Cyclin D1 level.

- **Skp2–p27 feedforward switch (Yao 2008).** *Why:* the proliferation–quiescence decision in these cells
  is switch-like and produces a transient, p27-high G0 dwell after each division. Heldt's Skp2 is
  constant; we made it a dynamic E2f target that is degraded by APC/C–Cdh1, creating a bistable
  Cyclin D1/p27-ratio commitment switch. This bistability is what produces the sharp threshold that a
  vismodegib-reduced Cyclin D1 must cross to arrest MB (and that Ezh2 inhibition restores), and it
  generates the transient G0 observed in GNP and MB cells; with a graded commitment, the threshold-based
  arrest and rescue would not emerge.

The model tracks the Cyclin/CDK complexes (D, E, A, B) and their CDK-inhibitor complexes, Rb/phospho-Rb,
E2f, Skp2, Cdh1, Cdc20, PCNA, the replication machinery (licensed and firing complexes and Dna), cell
mass, and the upstream Hedgehog/Mycn/Ezh2/Cyclin D1 species. The overall period was calibrated to
~22–23 h to match published GNP proliferation rates (70).

### Hedgehog signaling module
The Hedgehog pathway is an upstream driver of Cyclin D1 transcription. Shh binds Ptch1, relieving its
tonic inhibition of Smoothened (Smo); active Smo converts Gli repressor to Gli activator via a Hill
function (n = 2), and Gli activator drives Gli1 transcription. *Why the Gli→Ptch1 negative feedback was
included:* the GNP–MB contrast turns on whether the Cyclin D1 drive is mitogen-dependent, and the
canonical mechanism for this is the Hedgehog negative feedback in which Ptch1 is itself a Gli-induced
target that re-inhibits Smo (Marigo & Tabin 1996). We modeled Ptch1 transcription as a basal rate plus a
Gli-activator–induced term,

  d(Ptch1_mRNA)/dt = k_basal + k_Gli·Gli_act² / (K_Gli,Ptch² + Gli_act²) − k_deg·Ptch1_mRNA,

and let functional Ptch1 protein re-inhibit Smo,

  Smo_active = k_act / (1 + Ptch1_free·f / K_Ptch,Smo) × (1 − HHi),

where *f* is the functional-Ptch1 fraction (a context input: 1.0 in wild-type GNP, reduced to ~0.1 in
Ptch+/− MB to represent loss-of-function Ptch1) and HHi = 0 (untreated) or 1 (Vismodegib). With *f* near
1 the loop is intact and Gli is mitogen-responsive; with *f* ≈ 0.1 the loop is broken and Gli activity is
constitutive and Shh-independent, reproducing the elevated, ligand-independent Gli1 of MB as a
consequence of feedback loss (Goodrich 1997).

### Ezh2 module
Ezh2 mRNA transcription is the sum of a basal term and an E2f-driven term gated on the S-phase window —
proportional to E2f and to the S-phase cyclins (Cyclin E + Cyclin A):

  d(Ezh2_mRNA)/dt = k_basal + k_E2f·E2f·(CyclinE + CyclinA) − k_deg·Ezh2_mRNA.

*Why this form:* it makes Ezh2 a cell-cycle-dependent gene, high in S/G2 (E2f active, replication
ongoing) and low in G0, matching the ~2-fold higher Ezh2 in S/G2M vs G0 in the scRNA-seq (Fig. 4A)
(Bracken 2003; Pasini 2004), and — combined with a stable, slowly-degraded protein (half-life ~35 h)
that is diluted two-fold at division — it is what lets Ezh2 protein *integrate the duration of S-phase*,
the property that the explicit-replication core was chosen to support. Ezh2 represses all Cyclin D1
transcription through a repression factor,

  Ezh2_repression = K_Ezh2 / (K_Ezh2 + Ezh2·(1 − Ezh2i)),

with K_Ezh2 = 0.75 and Ezh2i = 0 (untreated) or 1 (full catalytic inhibition, modeling Tazemetostat),
producing ~2–3-fold Cyclin D1 de-repression upon Ezh2 inhibition, matching the experimental data
(Fig. 3C).

### Mycn module
*Why included:* a fraction of MB Cyclin D1 persists after Gli ablation by Vismodegib, and Mycn is the
HHi-resistant, autonomous driver that accounts for it. Mycn synthesis has an autonomous (basal) term
scaled by an amplification factor (1.0 for GNPs, 2.86 for MB) and a small Gli-dependent term that is not
amplified, making MB Mycn predominantly autonomous and Hedgehog-inhibitor–resistant. Mycn drives
Cyclin D1 transcription through a cooperative Hill function,

  Mycn_contribution = k_Mycn·Mycn^n / (K_Mycn^n + Mycn^n),  n = 3.7,

so that at GNP Mycn levels the contribution is small, whereas in MB it sustains the ~14% of Cyclin D1
that remains under Vismodegib.

### Cyclin D1 mRNA integration
Cyclin D1 mRNA transcription integrates basal, Gli-driven, and Mycn-driven inputs, each multiplied by the
Ezh2 repression factor:

  d(CyclinD1_mRNA)/dt = (k_basal + k_Gli·f(Gli_act, Gli1, Gli_rep) + k_Mycn·g(Mycn))·Ezh2_repression
                        − k_deg·CyclinD1_mRNA.

Calibration to the RNA-seq makes the Gli-driven term dominant — Gli1 and Cyclin D1 rise ~7-fold in
parallel from GNP to MB, identifying Cyclin D1 as predominantly Gli-driven, with Mycn supplying the
HHi-resistant residual. Cyclin D1 protein is translated from this mRNA and feeds the cell cycle through
Cyclin D–CDK4/6-mediated Rb phosphorylation.

### Commitment threshold: INK4 and CIP/KIP CDK inhibitors
*Why this module was added:* Vismodegib lowers MB Cyclin D1 only to roughly the level of a cycling GNP,
so the MB arrest cannot be explained by loss of Cyclin D1 alone — a GNP cycles at that level. The RNA-seq
shows that MB acquires a high CDK-inhibitor tone across both inhibitor families, which raises the
Cyclin D1 level required to commit. We modeled Cyclin D–CDK4/6 phosphorylation of Rb as saturating in
Cyclin D1 (CDK4/6 kinase activity is bounded) and competitively inhibited by the INK4 proteins:

  rate = kPhRbCd · Cd / (K_CdRb·(1 + p16 + p18) + Cd).

p16 (Cdkn2a) and p18 (Cdkn2c) are INK4-family CDK4/6 inhibitors entered as context inputs (Serrano 1993;
Sherr & Roberts 1999): p16 is silent in GNPs (0) and elevated in MB (0.88), while p18 is expressed in
GNPs (0.4) and further elevated in MB (1.2), consistent with the tumor-suppressor role of Ink4c/p18 in
medulloblastoma (Uziel 2005). Raising the INK4 tone increases the Cyclin D1 threshold for Rb
phosphorylation and commitment, but because the inhibition is competitive, additional Cyclin D1 overcomes
it. In parallel, the CIP/KIP inhibitors p21/p27 (the model's p27 species) inhibit Cyclin E/A–CDK2; their
synthesis rate is increased ~2-fold in MB relative to GNP, raising the threshold further and lengthening
the p27-high transient G0. *Why the competitive form matters:* representing the INK4 brake as competitive
(surmountable by more Cyclin D1) while modeling the CDK4/6 inhibitor as a non-competitive (Vmax) kinase
block is precisely what produces the rescuable-vs-non-rescuable distinction between the two perturbations.

### Drug implementations
HHi (Vismodegib) blocks Smo activation (HHi = 1, above). Ezh2i (Tazemetostat) sets Ezh2i = 1, abolishing
Ezh2 repression of Cyclin D1. CDK4/6 inhibition (Palbociclib) sets the Cyclin D–CDK4/6-mediated
Rb-phosphorylation rate (kPhRbCd) to zero — a non-competitive (Vmax) block that, unlike the competitive
INK4 brake, cannot be bypassed by raising Cyclin D1, and that leaves Cyclin D1 transcription unaffected.
Hydroxyurea (HU) reduces replication fork speed (the DNA-synthesis rate is scaled by a Hill function of
HU dose), lengthening S-phase concentration-dependently while still permitting mitosis once replication
completes (the function of the intra-S Chk1 checkpoint above).

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

## Model equations

All species are concentrations in arbitrary units; time is in minutes; `Cell` is the (constant)
compartment volume. Inputs/switches: **SHH** (ligand, 0–0.5), **GDC0449** = HHi (0/1), **EZH2i** (0/1),
**HU** (0/1), **p16, p18** (INK4 tone), **Ptch1_copy_number = f** (functional Ptch1 fraction),
**MYCN_amplification = A**, **CDK4/6i** (sets kPhRbCd→0), serum-starve (sets k_Cd_translation→0).
Equations are grouped by module; the four cell-cycle components we added on top of the core are marked
**[added]**. The Heldt (2018) G1/S core supplies the remaining standard reactions (Emi1/APC-Cdh1, PCNA
cycling, p53/DNA-damage); these are reproduced in the model file (`src/build_model_v44_heldt.py`,
`models_external/heldt2018.ant`) and the published model (BioModels BIOMD0000000700).

### Symbol glossary
`Cd_mRNA, Cd` = Cyclin D1 mRNA / protein · `Gli_act, Gli_rep` = Gli activator/repressor · `Gli1` ·
`MYCN` · `Ptch1_mRNA, Ptch1_free, SHH_Ptch` · `Smo_active` · `EZH2m, EZH2` = Ezh2 mRNA/protein ·
`Rb, pRb, E2f, RbE2f` = Rb / phospho-Rb / E2f / Rb·E2f complex · `Ce, Ca` = Cyclin E– and Cyclin A–CDK2 ·
`P21` = p21/p27 (CIP/KIP) · `CeP21, CaP21` = CDK2·p21 complexes · `Skp2` · `C1, pC1` = APC/C-Cdh1
(active/phospho) · `Rc, pRc, aRc, iRc` = replication complexes (licensed/primed/active/inactive) ·
`Dna` (0→1) · `aPcna, iPcna` = PCNA · `MPF, preMPF` = Cyclin B–CDK1 (active/Tyr15-phos) · `Cdc20` ·
`mass` · `P53, Dam` = p53 / DNA damage.

### Module 1 — Hedgehog → Gli, with the Gli→Ptch1 negative feedback
```
d[Ptch1_mRNA]/dt = k_Ptch1_basal + k_Ptch1_Gli·Gli_act² / (K_Gli_Ptch² + Gli_act²)        ← Gli-induced
                   − k_Ptch1_mRNA_deg·Ptch1_mRNA
d[Ptch1_free]/dt = k_Ptch1_translation·Ptch1_mRNA − k_Ptch1_deg·Ptch1_free
                   − k_SHH_Ptch_bind·SHH·Ptch1_free + k_SHH_Ptch_release·SHH_Ptch
d[SHH_Ptch]/dt   = k_SHH_Ptch_bind·SHH·Ptch1_free − (k_SHH_Ptch_release + k_SHH_Ptch_deg)·SHH_Ptch
d[Smo_active]/dt = k_Smo_act / (1 + Ptch1_free·f / K_Ptch_Smo) · (1 − GDC0449) − k_Smo_inact·Smo_active
                                          ↑ functional Ptch1 (f) re-inhibits Smo; HHi blocks Smo
   let  S2 = Smo_active² / (K_Smo_Gli_switch² + Smo_active²)            (Smo→Gli switch, Hill n=2)
d[Gli_act]/dt   =  k_Gli_rep_to_act·S2·Gli_rep − k_Gli_act_to_rep·(1 − S2)·Gli_act
d[Gli_rep]/dt   = −(above)                                              (Gli_act + Gli_rep conserved)
d[Gli1_mRNA]/dt =  Vmax_Gli1_tx·Gli_act²/(K_Gli_act_Gli1² + Gli_act²)·(1 − Gli_rep²/(K_Gli_rep_Gli1² + Gli_rep²))
                   − k_Gli1_mRNA_deg·Gli1_mRNA
d[Gli1]/dt      =  k_Gli1_translation·Gli1_mRNA − k_Gli1_deg·Gli1
```
*f = Ptch1_copy_number*: 1.0 in GNP (loop intact, Gli mitogen-responsive) → ~0.1 in MB (loop broken,
Gli constitutive).

### Module 2 — Mycn
```
d[MYCN]/dt = k_MYCN_synth_basal·A + k_MYCN_synth_Gli·(Gli_act + Gli1)/(K_Gli_MYCN + Gli_act + Gli1)
             − k_MYCN_deg·MYCN                                          (A = MYCN_amplification)
```

### Module 3 — Cyclin D1 (the integration node) and the Ezh2 feedback
```
Ezh2 repression factor:   R_EZH2 = K_EZH2_repression / (K_EZH2_repression + EZH2·(1 − EZH2i))

d[Cd_mRNA]/dt = [ k_Cd_tx_basal
                + k_Cd_tx_Gli_max · (Gli_act+Gli1)²/(K_Gli_act_CycD² + (Gli_act+Gli1)²)
                                  · K_Gli_rep_CycD²/(K_Gli_rep_CycD² + Gli_rep²)            ← Gli-driven
                + k_Cd_tx_MYCN · MYCN^n_MYCN / (K_MYCN_Cd^n_MYCN + MYCN^n_MYCN) ] · R_EZH2   ← Mycn-driven
                − k_Cd_mRNA_deg·Cd_mRNA                                   (n_MYCN = 3.66)
d[Cd]/dt      = k_Cd_translation·Cd_mRNA − k_Cd_deg·Cd

Ezh2 (E2f-driven, gated on the S-phase cyclins → integrates time-in-S):
d[EZH2m]/dt = kEZbas + kEZE2f · E2f/(K_E2f_EZ + E2f) · [ wCe·Ce/(K_Ce_EZ + Ce) + (1−wCe)·Ca/(K_Ca_EZ + Ca) ]
              − kDeEZm·EZH2m
d[EZH2]/dt  = kTlEZ·EZH2m − kDeEZ·EZH2          (EZH2 and EZH2m additionally halved at each division)
```

### Module 4 — Cell-cycle engine

**Restriction point (Rb–E2f) with the INK4 commitment brake [added: saturating + competitive form]**
```
Rb-kinase activity per Rb:
   v_Rb = kPhRbCd·Cd / (K_CdRb·(1 + p16 + p18) + Cd) · commit_gate   +  kPhRbCe·Ce + kPhRbCa·Ca
                         ↑ Cyclin D–CDK4/6, saturating in Cd and competitively raised by INK4 (p16+p18)
   Rb → pRb              at  v_Rb · Rb
   RbE2f → pRb + E2f     at  v_Rb · RbE2f          (phospho-Rb releases E2f)
   pRb → Rb             at  kDpRb · pRb
d[E2f]/dt  = kSyE2f + kSyE2fE2f·E2f/(jSyE2f + E2f) − kDeE2f·E2f   (+ release from RbE2f; Rb+E2f ⇌ RbE2f)
d[Ce]/dt   = kSyCe·E2f − (kDeCe + kDeCeCa·Ca)·Ce  (± p21 binding)   ;  d[Ca]/dt = kSyCa·E2f − … (± p21)
```

**p21 / p27 (CIP/KIP) with Skp2-dependent degradation**
```
d[P21]/dt = (kSyP21 + kSyP21P53·P53)                                    ← kSyP21 is the MB-elevated knob
            − (kDeP21 + kDeP21Cy·Skp2·(Ce+Ca) + kDeP21aRc·Cdt2·aRc)·P21  (+ exchange with CeP21/CaP21/PCNA)
```

**Skp2–p27 feedforward [added: Skp2 made a dynamic E2f target]**
```
d[Skp2]/dt = (kSySkp2bas + kSySkp2·E2f) − kDeSkp2C1·C1·Skp2 − kDeSkp2bas·Skp2
```
(E2f raises Skp2 → Skp2 degrades p27 → CDK2 active → more E2f: the bistable commitment feedforward;
APC/C-Cdh1 (C1) keeps Skp2 low in G0.)

**Cell growth and size gates [added]**
```
d[mass]/dt  = mu·mass                                  (mu = 0.0005 min⁻¹; mass → mass/2 at division)
size_gate   = mass^n / (M_size^n + mass^n)             (M_size = 2.5; gates S-entry / origin firing)
commit_gate = mass^n / (M_commit^n + mass^n)           (M_commit = 1.3; gates G0→G1 commitment) ; n = 6
```

**Explicit DNA replication (Heldt core; HU enters here)**
```
Rc → pRc      at  kPhRc·(Ce+Ca)^n/(jCy^n + (Ce+Ca)^n) · size_gate · Rc      (origin licensing, size-gated)
pRc (+PCNA) → aRc                                                            (firing)
d[Dna]/dt    =  kSyDna · vfork · aRc                                         (DNA synthesis, 0 → 1)
vfork        =  vmin_fork + (1 − vmin_fork)·KmHU_fork^h / (KmHU_fork^h + HU^h)   ← HU slows fork speed
```

**Cyclin B/CDK1 mitotic switch + intra-S Chk1 checkpoint [added]**
```
g2gate     = Dna^nG2 / (KG2^nG2 + Dna^nG2)             (≈1 only when replication ≈ complete; KG2=0.85, nG2=6)
Chk1       = aRc / (jChk + aRc)                        (active forks = unfinished S → inhibits mitosis)
Cdc25a     = (a25 + (1−a25)·MPF^n/(KmMpf^n + MPF^n)) / (1 + wChk·Chk1)        (activating, hysteretic)
Wee1a      = (aWee + (1−aWee)·KmMpf^n/(KmMpf^n + MPF^n)) · (1 + wChkW·Chk1)    (inactivating)
d[preMPF]/dt = kSyCb·g2gate + kWee·Wee1a·MPF − k25·Cdc25a·preMPF − (kDeCbBas + kDeCb·Cdc20)·preMPF
d[MPF]/dt    = k25·Cdc25a·preMPF − kWee·Wee1a·MPF      − (kDeCbBas + kDeCb·Cdc20)·MPF
d[Cdc20]/dt  = kaCdc20·MPF^nCdc20/(KmCdc20^nCdc20 + MPF^nCdc20)·(1 − Cdc20) − kiCdc20·Cdc20
              (Cdc20 then drives mitotic Cyclin B and Cyclin A destruction → mitotic exit)
```

**Division event** (when MPF crosses MPF_div):
```
Dna→0; Rc→1, pRc=aRc=iRc=0; Rb→Rb+pRb, pRb→0; P21→P21_div (high); Skp2→0.05;
Ce→Ce_div, Ca→Ca_div; MPF=preMPF=0; Cdc20→0; mass→mass/2; EZH2→EZH2/2, EZH2m→EZH2m/2.
```
(The p27 reset-high + halved mass create the growth-timed, p27-high transient G0; halving EZH2 is the
dilution that makes EZH2 a time-in-S integrator.)

**Drug switches** (recap): GDC0449 in `Smo_active`; EZH2i in `R_EZH2`; HU in `vfork`; CDK4/6i sets
`kPhRbCd = 0` (a Vmax block of `v_Rb`'s Cyclin D term — note this cannot be overcome by raising Cd,
unlike the competitive p16/p18 term); serum-starve sets `k_Cd_translation = 0`.

---

## References
Heldt et al. 2018 (PNAS 115:2532); Bracken et al. 2003 (EMBO J 22:5323); Pasini et al. 2004 (Cell Cycle
3:396); Marigo & Tabin 1996 (PNAS 93:9346); Goodrich et al. 1997 (Science 277:1109); Serrano et al. 1993
(Nature 366:704); Sherr & Roberts 1999 (Genes Dev 13:1501); Uziel et al. 2005 (Genes Dev 19:2656);
Pomerening et al. 2003 (Nat Cell Biol 5:346); Yao et al. 2008 (Nat Cell Biol 10:476); Wechsler-Reya &
Scott 1999 (Neuron 22:103); Spencer et al. 2013 (Cell 155:369); Cappell et al. 2016 (Cell 166:167).
Retained from the manuscript: Gérard & Goldbeter 2009 (ref 30); Tellurium (ref 69); GNP period (ref 70).
Verify volume/page and renumber to the manuscript scheme.
