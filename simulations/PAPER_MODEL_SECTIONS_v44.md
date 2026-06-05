# Computational model — replacement text for Chahin et al. (v44 model)

This document replaces the v42 model text in the manuscript. It is organized as:
1. **Results** paragraph (replaces "A quantitative model of the Ezh2-Cyclin D1 feedback loop…")
2. **Methods → Computational Modeling** (replaces the v42 Computational Modeling subsection)
3. **Figure captions** (Fig. 5A–C, Supp. Fig. 7, Supp. Fig. 8) updated for the v44 model
4. **Summary of what changed from the published v42 text** (for the authors)

Numbers below are from `simulations/validate_v44.py` (the calibrated v44 model). Where the new
model gives a *different qualitative result* than the v42 model, this is flagged explicitly — the
most important being that, in v44, the EZH2→Cyclin D1 feedback sets Cyclin D1 *level/amplitude*
but does **not** lengthen the cell-cycle period (the 4.3 h lengthening in v42 was a consequence of
the global "ε" clock, which v44 removes).

---

## 1. RESULTS

### A quantitative model of the Ezh2–Cyclin D1 feedback loop predicts context-dependent drug responses in GNPs and medulloblastoma

To formalize the negative feedback loop identified in Fig. 4 and explore its quantitative
consequences, we constructed an ordinary differential equation (ODE) model integrating Hedgehog
signaling, Mycn, Ezh2, and the mammalian cell-cycle engine (Fig. S7A). As in the experimental
system, Cyclin D1 mRNA transcription is the central integration node, receiving basal,
Gli-driven (Shh→Ptch1→Smo→Gli), and Mycn-driven inputs, all subject to Ezh2-mediated repression
through H3K27me3 at the Cyclin D1 promoter. Unlike a phenomenological oscillator, the cell-cycle
core models **DNA replication explicitly**: it is built on the real-time mammalian G1/S model of
Heldt et al. (2018), in which an Rb–E2f restriction point licenses origins that fire and
synthesize DNA at a defined fork speed, so that **S-phase has a mechanistic duration** rather than
a fixed, clock-imposed one. We extended this core with a Cyclin B/CDK1 mitotic switch
(Cdc25/Wee1 hysteresis), an intra-S CHK1 checkpoint that holds mitosis until replication
completes, a cell-growth–gated restriction point (so G1 length is set by the time to reach a
critical size, decoupled from Cyclin D1 level), and a Skp2–p27 feedforward switch that makes the
G0→G1 commitment bistable. Ezh2 transcription is driven by E2f during the S-phase window, and
because Ezh2 protein is stable and diluted only at division, **Ezh2 integrates the duration of
S-phase**, closing the feedback loop: cycling cells with longer S accumulate more Ezh2, which
represses Cyclin D1.

Model parameters were calibrated against RNA-seq data from P7 GNPs and Ptch+/− medulloblastoma
± Hedgehog pathway inhibitor (vismodegib), including Cyclin D1, Mycn, Gli1, and Ezh2 transcript
ratios across conditions, the within-cell-cycle Ezh2 gradient, and the hydroxyurea-induced Ezh2
increase (see Methods). The model reproduced the Hedgehog-dependence of proliferation (GNPs cycle
only with Shh; GNP−Shh, GNP+HHi and serum-starved GNPs arrest, whereas Ptch+/− MB cycles
Shh-independently), the near-complete loss of Cyclin D1 in GNP+HHi versus the partial retention in
MB+HHi (Mycn-buffered), the cell-cycle Ezh2 transcript gradient (≈2-fold higher in S/G2 than G0),
and the ~1.3-fold increase in S-phase Ezh2 protein under hydroxyurea, while producing a cell-cycle
period of ~23 h (Fig. S7B–M).

**Ezh2 feedback constrains Cyclin D1 and makes it oscillate, but does not lengthen the cell
cycle.** To quantify the impact of Ezh2-mediated Cyclin D1 repression, we compared simulations
with and without the feedback (by computationally removing Ezh2's ability to repress Cyclin D1).
In both Shh-stimulated GNPs and Ptch+/− MB, removing the feedback raised mean Cyclin D1 mRNA
~3-fold and converted it from a cell-cycle-coupled oscillation into a constitutively high signal
(Fig. S8C,D). This oscillation arises because Ezh2 protein tracks E2f activity through the cycle
and represses Cyclin D1 most strongly in S/G2 (Fig. S8G,H). In contrast to a model with a global
cell-cycle clock, removing the feedback in the explicit-replication model changed the cell-cycle
period only marginally (Fig. S8A,B), because S-phase and G1 durations are set by replication and
cell growth rather than by Cyclin D1 abundance — consistent with reports that mammalian
cell-cycle phase durations are largely uncoupled. The model therefore predicts that Ezh2 acts
chiefly at the G0/G1 commitment boundary, where Cyclin D1–dependent Rb phosphorylation is
rate-limiting, consistent with the experimental observation that Ezh2 overexpression specifically
traps cells in G0 (Fig. 3F).

**CyclinD1 is predominantly Gli-driven, and the Ezh2i rescue of HHi is not a simple CyclinD1
threshold effect.** When the Hedgehog module is calibrated to the bulk RNA-seq (Gli1 and Cyclin D1
both ~7-fold higher in Ptch+/− MB than in P7 GNPs, rising in lockstep along the Ptch+/− series; see
Methods), the model reproduces the Gli-driven elevation of Cyclin D1 in MB. A direct consequence is
that vismodegib-treated MB retains a substantial Cyclin D1 level — the Mycn-driven, HHi-resistant
fraction is ~3-fold that of a cycling GNP — which remains **above** the commitment threshold. The
model therefore predicts that the reduction in MB proliferation by vismodegib, and its restoration
by Ezh2 inhibition (Fig. 5D,E), are **not** explained by Cyclin D1 falling below the proliferative
threshold; the residual Cyclin D1 is too high for that. This indicates that the HHi/Ezh2i effect on
proliferation acts through additional Hedgehog outputs beyond Cyclin D1 (e.g. other Gli targets or
survival pathways), consistent with the partial, population-level changes seen experimentally rather
than an all-or-nothing arrest. **The one prediction the Cyclin D1 axis does make cleanly is the
CDK4/6-inhibitor contrast:** CDK4/6 inhibition blocks Cyclin D–CDK4/6 kinase activity downstream of
Cyclin D1 transcription, so it arrests cells in a way Ezh2 inhibition cannot rescue (de-repressing
Cyclin D1 transcript cannot bypass a kinase-level block; Fig. S8I,J) — a robust, mechanism-specific
prediction that survives the data-matched calibration.

[NOTE TO AUTHORS: This replaces the earlier "Ezh2i de-represses Cyclin D1 above the threshold and
rescues cycling" framing. That clean single-cell rescue only held when the model under-estimated MB
Cyclin D1 (MB+HHi below the GNP level). With the measured Cyclin D1 (MB+HHi ~3x a cycling GNP), the
rescue is not a Cyclin D1-threshold effect — a genuine finding from fitting the RNA-seq. Decide how
much to foreground this vs the CDK4/6i contrast in the main text.]

---

## 2. METHODS — Computational Modeling

### Model construction

We developed an ODE-based model (v44) of the Ezh2–Cyclin D1 negative feedback loop embedded
within the mammalian cell cycle, integrating four functional modules: (1) Hedgehog signaling,
(2) Mycn regulation, (3) Ezh2 regulation, and (4) the cell-cycle engine. The model was implemented
in Antimony and simulated using Tellurium with the CVODE integrator in Python 3. In contrast to
the relaxation-oscillator core used previously, whose S-phase duration is a structural invariant
set only by a global time-scaling factor, the v44 cell-cycle engine is a real-time model with
explicit DNA replication, so that the duration of S-phase is determined mechanistically by
replication kinetics. This was required to model the dependence of Ezh2 accumulation on time spent
in S-phase.

**Cell-cycle engine.** The core was built on the real-time mammalian G1/S model of Heldt, Barr,
Cooper, Bakal & Novák (2018) (PNAS; BioModels BIOMD0000000700). This core contains an Rb–E2f
restriction point (Rb is phosphorylated by Cyclin D/CDK4,6, Cyclin E/CDK2 and Cyclin A/CDK2;
phospho-Rb releases E2f, which drives Cyclin E and Cyclin A transcription in a positive-feedback
loop), the CDK inhibitor p21/p27, APC/C-Cdh1, PCNA, and **explicit DNA replication** in which
licensed origins (Rc) are fired by CDK2 activity to active replication complexes (aRc) that
synthesize DNA (Dna: 0→1) at fork speed kSyDna. Because the Heldt model is a one-shot
proliferation–quiescence model with no mitosis, we added: (i) a **Cyclin B/CDK1 (MPF) mitotic
switch** with Cdc25/Wee1 hysteresis, in which Cyclin B/CDK1 is synthesized (as Tyr15-phosphorylated
pre-MPF) only as replication completes (a Dna-dependent gate), so that a genuine multi-hour G2
elapses before sharp mitotic entry; (ii) an **intra-S CHK1 checkpoint** gated on active
replication forks (Chk1 ∝ aRc/(jChk+aRc)), holding MPF inactive until replication finishes; (iii)
a mitotic **APC/Cdc20** that destroys Cyclin A/B at mitotic exit; and (iv) a **division-reset
event** that, at mitotic entry, resets DNA, re-licenses origins, dephosphorylates Rb, sets p27
high and Skp2 low (a newborn G0 state), halves cell mass, and dilutes Ezh2 two-fold. The model is
real-time (minutes); no global time-scaling factor is used.

**Cell growth and the restriction point.** A cell-mass variable grows exponentially and halves at
division (size homeostasis). Two size gates couple growth to the cycle: commitment (the Cyclin
D–driven step of Rb phosphorylation) and S-entry (origin firing) each require the cell to exceed a
critical size. Consequently the duration of G1 is set by the time required to grow to commitment
size and is decoupled from the level of Cyclin D1 — this allows MB cells with high Cyclin D1 to
retain a substantial G1, as observed. Commitment remains mitogen- (Cyclin D1-) dependent, so
low-Cyclin D1 cells (GNP without Shh, or +HHi) never commit and remain quiescent, preserving the
Hedgehog-dependence of proliferation.

**Skp2–p27 feedforward switch.** Skp2 was made a dynamic E2f target that is degraded by APC/C-Cdh1.
In G0, E2f is off and Cdh1 is active, so Skp2 is low and p27 is not degraded (p27-high,
phospho-Rb-negative quiescence). Cyclin D1–driven partial Rb phosphorylation releases some E2f,
raising Skp2, which degrades p27, activating Cyclin E/CDK2, which inactivates Cdh1 and further
stabilizes Skp2 — a bistable feedforward that constitutes the G0→G1 commitment switch. p27 is reset
high at division. (Note: the cyclin D1/p27 balance setting this switch is consistent with the
ultrasensitive Cyclin D1/p27-ratio threshold of Fan & Meyer 2021. An explicit two-step Rb —
Cyclin D–CDK4,6 mono- then Cyclin E/A–CDK2 hyper-phosphorylation, per Narasimha et al. 2014 and
Sanidas et al. 2019 — is available in the model as an optional variant but is not used here, because
it does not sustain the observed transient G0 stably; G0 is classified by the p27 marker, for which
the default model is accurate.)

**Hedgehog signaling module.** The Hedgehog pathway drives Cyclin D1 transcription. Shh binds
Ptch1, relieving inhibition of Smo; active Smo converts Gli repressor to Gli activator (Hill
n = 2), and Gli activator drives Gli1 transcription in a positive-feedback loop. Ptch1 gene dosage
is a copy-number input (1.0 = wild-type GNP; reduced in MB, reflecting Ptch1 loss and constitutive
Hedgehog activity). HHi (vismodegib/GDC0449) is modeled as a block on Smo activation
(Smo_active = k_act/(1+Ptch1_free/K)×(1−GDC0449)).

**Mycn module.** Mycn synthesis has an autonomous (basal) component scaled by an amplification
factor (1.0 for GNP, 2.8 for MB) and a small Gli-dependent component. Mycn drives a minority,
HHi-resistant component of Cyclin D1 transcription through a cooperative Hill function. This
explains why Ptch+/− MB retains a substantial Cyclin D1 fraction under HHi (only ~40% reduction)
whereas GNP loses nearly all Cyclin D1 (~86% reduction). Consistent with the bulk RNA-seq, the
dominant driver of the MB/GNP Cyclin D1 elevation is Gli (see Cyclin D1 integration), not Mycn.

**Ezh2 module.** Ezh2 mRNA transcription is an E2f-driven term gated on the S-phase window (the
rate is proportional to E2f×(Cyclin E + Cyclin A)), plus a small basal term; Ezh2 protein is
translated from this mRNA, is stable (slow first-order degradation), and is diluted two-fold at
each division. Because protein turnover is slow relative to the cycle, Ezh2 protein integrates the
duration of the S-phase window — a longer S-phase yields more Ezh2. Ezh2 represses all Cyclin D1
transcription through a repression factor K_Ezh2/(K_Ezh2 + Ezh2×(1−Ezh2i)), where Ezh2i = 0
(untreated) or 1 (full catalytic inhibition, modeling tazemetostat). Ezh2 inhibition therefore
de-represses Cyclin D1.

**Cyclin D1 integration and drug implementations.** Cyclin D1 mRNA transcription integrates a basal
term, a dominant Gli-driven term (calibrated to the ~7-fold Gli1/Cyclin D1 lockstep rise GNP→MB),
and a minority Mycn-driven term, each multiplied by the Ezh2 repression factor. Cyclin D1 protein is
translated from this mRNA and feeds the cell cycle through Cyclin D/CDK4,6-mediated Rb
phosphorylation, implemented as a *saturating* drive, kPhRbCd·Cd/(K_CdRb+Cd), so that the
data-matched (up to ~7×) transcript level maps to a bounded kinase activity (CDK4/6 saturates) and
the integrator stays stable. CDK4/6 inhibition (palbociclib) was implemented by setting this Cyclin
D/CDK4,6-mediated Rb-phosphorylation rate to zero, leaving Cyclin D1 transcription unaffected.
Hydroxyurea was implemented as a reduction of replication fork speed (vfork scales kSyDna with a
Hill dependence on the HU dose), so HU lengthens S-phase concentration-dependently while leaving
G1 unchanged and allowing mitosis once replication completes (no hard arrest).

### Parameter calibration

Parameters were calibrated against quantitative RNA-seq and immunofluorescence data from P7
wild-type GNPs and Ptch+/− medulloblastoma, ± vismodegib (HHi). Calibration targets and the
calibrated model values were:

| Target | Experimental | Model (v44) | Source |
|--------|--------------|-------------|--------|
| Cyclin D1: GNP+HHi / GNP | 0.14 (86% reduction) | **0.14** | RNA-seq |
| Cyclin D1: MB+HHi / MB | 0.40 (60% reduction) | **0.43** | RNA-seq |
| Cyclin D1: MB / GNP | 7.58× (bulk) / 5.07× (IF) | **7.05×** | RNA-seq, Fig. 4I |
| Gli1: MB / GNP | 6.90× | **6.30×** | bulk RNA-seq |
| Mycn: GNP+HHi / GNP | 0.78 | **0.79** | RNA-seq |
| Mycn: MB+HHi / MB | 0.86 | **0.91** | RNA-seq |
| Mycn: MB / GNP | ~2.8× | **2.7×** | RNA-seq |
| Gli1: GNP+HHi reduction | >99% | **>99%** | RNA-seq |
| Ezh2: G0 / cycling | ~0.6 | **0.7** | scRNA-seq (Fig. 4A) |
| Ezh2: MB / GNP | 2.05× | **see note** | RNA-seq, Fig. 4J |
| Ezh2i: Cyclin D1 fold change | ~2× | **2–3×** | qPCR (Fig. 3C) |
| Ezh2 transcript: S/G0 | 1.8–2.5× | **~2.0–2.9×** | scRNA-seq (Fig. 4A) |
| Ezh2 protein: G2/G0 | 1.48× | **~1.5–1.8×** | IF (Supp. 6D) |
| HU: Ezh2-in-S / DMSO | 1.31× | **~1.3×** | IF (Fig. 4G) |
| Cell-cycle period (GNP) | ~22–23 h | **~23 h** | published (ref. Nakashima) |

*Hedgehog de-saturation and the saturating Cyclin D1→Rb drive.* The bulk RNA-seq shows Gli1 and
Cyclin D1 rising ~7-fold in lockstep from P7 GNP to Ptch+/− MB (and along the intermediate Ptch+/−
series), identifying Cyclin D1 as predominantly Gli-driven. The original Hedgehog module was
saturated (Gli1 MB/GNP only ~1.1×), so we re-fit the Ptch→Smo→Gli and Gli→Cyclin D1 Hill terms to
de-saturate the pathway; the calibrated module now reproduces Gli1 MB/GNP 6.30× and Cyclin D1 MB/GNP
7.05×, with Mycn supplying the HHi-resistant residual (GNP+HHi 0.14, MB+HHi 0.43). Because Cyclin D1
protein in the original engine drove Rb phosphorylation *linearly* (rate ∝ kPhRbCd·Cd), a 7× Cyclin
D1 level produced a 7× cell-cycle drive that was both biologically implausible (CDK4/6 activity
saturates) and numerically stiff. We therefore made the Cyclin D1→Rb drive saturating,
kPhRbCd·Cd/(K_CdRb+Cd), which decouples the (now data-matched) transcript level from a bounded
cell-cycle drive. The Ezh2 MB/GNP ratio remains cell-cycle-coupled and is treated as a soft target.

*Consequence for the rescue interpretation.* With Cyclin D1 matched to the data, MB+HHi retains
~3-fold a cycling GNP's Cyclin D1 (above the commitment threshold), so the model no longer produces
an all-or-nothing Cyclin D1-threshold arrest/rescue. We report this directly in the Results: the
vismodegib/Ezh2i proliferation effect is not a simple Cyclin D1-threshold switch and requires
additional Hedgehog biology, while the CDK4/6i-cannot-be-rescued contrast (a downstream kinase
block) is preserved.

A constrained parameter search over the Hedgehog/Mycn→Cyclin D1 and Ezh2-repression parameters was
performed with a hard guard that rejected any parameter set that (i) failed to integrate
numerically on any condition or (ii) broke the figure-critical cycling behavior (GNP+Shh and MB
must cycle; GNP−Shh / GNP+HHi must arrest; MB+CDK4/6i must arrest and not be rescued by Ezh2i),
minimizing a weighted error on the expression ratios (Gli1 and Cyclin D1 MB/GNP, the HHi reductions,
and the Mycn folds) among survivors.

### Phase classification

Cell-cycle phase proportions were computed as time fractions of a single cycling trajectory
(equivalent to the proportions of an asynchronous population). Consistent with the experimental
immunofluorescence assay (phospho-Rb Ser807/811 + DAPI/PCNA), cells were classified as **G0**
(pre-replication, p27-high), **G1** (pre-replication, p27-low/committed), **S** (origins fired and
replicating, 0 < Dna < 1), or **G2/M** (Dna ≈ 1 until division). G0 was delineated by the **p27**
marker, consistent with the manuscript G0 definition (phospho-Rb-negative **or** p27-positive); by
this criterion the model reproduces a transient G0 of ~19% in MB and ~32% in GNP that collapses
under Ezh2 inhibition. The hydroxyurea-induced redistribution (S up, G2 down) is a prediction.

### Simulation conditions

All simulations used Tellurium with the CVODE integrator (absolute tolerance 1×10⁻⁹, relative
tolerance 1×10⁻⁶). Cell-type and drug context were set as follows (Shh; Ptch1 copy number; Mycn
amplification; GDC0449/HHi; Ezh2i; HU; CDK4/6i):

| Condition | Shh | Ptch1 | Mycn amp | HHi | Ezh2i | HU | CDK4/6i |
|-----------|-----|-------|----------|-----|-------|----|---------|
| GNP + SHH | 0.5 | 1.0 | 1.0 | 0 | 0 | 0 | No |
| GNP − SHH | 0.0 | 1.0 | 1.0 | 0 | 0 | 0 | No |
| GNP + HHi | 0.5 | 1.0 | 1.0 | 1 | 0 | 0 | No |
| GNP + EZH2i | 0.5 | 1.0 | 1.0 | 0 | 1 | 0 | No |
| GNP + HHi + EZH2i | 0.5 | 1.0 | 1.0 | 1 | 1 | 0 | No |
| GNP Serum-starved | 0.5 | 1.0 | 1.0 | 0 | 0 | 0 | No (Cyclin D1 translation = 0) |
| MB (Ptch+/−) | 0.5 | 0.1 | 2.8 | 0 | 0 | 0 | No |
| MB + HHi | 0.5 | 0.1 | 2.8 | 1 | 0 | 0 | No |
| MB + EZH2i | 0.5 | 0.1 | 2.8 | 0 | 1 | 0 | No |
| MB + HHi + EZH2i | 0.5 | 0.1 | 2.8 | 1 | 1 | 0 | No |
| MB + CDK4/6i | 0.5 | 0.1 | 2.8 | 0 | 0 | 0 | Yes |
| MB + CDK4/6i + EZH2i | 0.5 | 0.1 | 2.8 | 0 | 1 | 0 | Yes |

(In v44, MB is modeled by reduced Ptch1 copy number — constitutive Hedgehog — plus Mycn
amplification, with Shh present, rather than by zeroing Ptch1 and Shh as in the previous model.)

Validation simulations were run for ~200 h; feedback-impact and rescue simulations for 168 h. Cell
divisions were counted as peaks in Cyclin B/CDK1 (MPF) after an initial transient; the cell-cycle
period was the mean interval between consecutive peaks. Steady/mean species levels were the mean
over the settled portion of the trajectory.

### Ezh2 feedback analysis

The impact of Ezh2-mediated Cyclin D1 repression was assessed by comparing simulations with intact
feedback versus computationally ablated feedback (K_Ezh2 set very large, making the repression
factor ≈ 1 regardless of Ezh2). This was done in both GNP+Shh and MB contexts. Phase relationships
between Cyclin D1 mRNA, Ezh2 protein, and Cyclin B were visualized by normalizing each species to
[0, 1] over the last 70 h of a 168-h simulation.

### Software and code availability

The model was implemented in Antimony (Smith et al., 2009) and simulated with Tellurium (Choi et
al., 2018) in Python 3; the cell-cycle core derives from the SBML model of Heldt et al. (2018)
(BioModels BIOMD0000000700). SciPy was used for peak detection (scipy.signal.find_peaks) and
Matplotlib for figure generation. Model source, the validation suite (`simulations/validate_v44.py`),
and the figure scripts are available at https://github.com/jpurzner/Ezh2_CyclinD1_Sphase.

### New references to add

- Heldt FS, Barr AR, Cooper S, Bakal C, Novák B. A comprehensive model for the proliferation–
  quiescence decision in response to endogenous DNA damage in human cells. *PNAS* 2018;115:2532–2537.
  (BioModels BIOMD0000000700)
- Narasimha AM, Kaulich M, Shapiro GS, et al. Cyclin D activates the Rb tumour suppressor by
  mono-phosphorylation. *eLife* 2014;3:e02872.
- Sanidas I, Morris R, Fella KA, et al. A code of mono-phosphorylation modulates the function of RB.
  *Mol Cell* 2019;73:985–1000.
- Fan Y, Meyer T. Molecular control of cell density-mediated exit to quiescence. *Cell Reports*
  2021;36:109436. (the Cyclin D1/p27 ratio threshold)

---

## 3. FIGURE CAPTIONS (updated for v44)

**Figure 5.** v44 model of the MB Cyclin D1 axis, calibrated to the bulk RNA-seq. (A) Calibrated
expression: Gli1 and Cyclin D1 rise ~6–7-fold in lockstep from P7 GNP to Ptch+/− MB, with Mycn the
HHi-resistant residual (model vs data). (B) Predicted Cyclin B/CDK1 oscillation in cycling Ptch+/−
MB. (C) Population (N = 120 cells, cyclin D1 / p27 heterogeneity): the cycling fraction falls modestly
under HHi and is restored by adding Ezh2i — a fractional, not all-or-nothing, effect, matching the
fractional Fig. 5D/E proliferation bars. (Mechanistic note: because the data place MB+HHi Cyclin D1
~3-fold above a cycling GNP, this fractional rescue is *not* a Cyclin D1-threshold switch — see
Results.) (D) Contrast: CDK4/6 inhibition arrests MB and Ezh2i cannot rescue it (downstream kinase
block), the prediction the Cyclin D1 axis makes cleanly.

**Supplementary Figure 7. Model architecture and validation.** (A) Architecture of the v44 ODE
model: Hedgehog/Mycn→Cyclin D1 input, the Heldt-based cell-cycle engine with explicit DNA
replication, the Skp2–p27 growth-gated restriction point, the Cyclin B/CDK1 mitotic switch and
intra-S CHK1 checkpoint, and the Ezh2 layer that integrates S-phase duration and represses Cyclin
D1. (B) Cyclin B traces for GNPs (+SHH, −SHH, +HHi, +Ezh2i, +HHi+Ezh2i). (C) Simulated vs measured
Cyclin D1 mRNA fold change across GNP conditions. (D) Simulated vs measured Mycn in GNP and
GNP+HHi. (E) Simulated Ezh2 in cycling, serum-starved (G0) and HHi-treated GNPs. (F) Cyclin B
traces for Ptch+/− MB, +HHi, +Ezh2i, +HHi+Ezh2i. (G) Simulated vs measured Cyclin D1 mRNA across
MB conditions. (H) Simulated vs measured Mycn in MB and MB+HHi. (I) Simulated vs measured Gli1 in
MB conditions. (J) Cell-division counts across all conditions. (K) Cyclin D1 mRNA and (L) Ezh2
protein time courses in MB conditions (the Ezh2 sawtooth reflects accumulation across S-phase and
two-fold dilution at division). (M) Mycn protein in GNP+SHH, GNP+HHi, MB and MB+HHi.

**Supplementary Figure 8. Ezh2-mediated Cyclin D1 repression constrains Cyclin D1 and drives its
oscillation.** Cyclin B traces with and without Ezh2 feedback in (A) GNP+SHH and (B) Ptch+/− MB
(period is essentially unchanged — see note). Cyclin D1 mRNA with and without Ezh2 repression in
(C) GNP+SHH and (D) Ptch+/− MB (removing the feedback raises Cyclin D1 ~3-fold and abolishes its
oscillation). Ezh2 protein with and without feedback in (E) GNP+SHH and (F) Ptch+/− MB. Normalized
Cyclin D1 mRNA, Ezh2 protein and Cyclin B over the last 70 h in (G) GNP+SHH and (H) Ptch+/− MB
(Ezh2 peaks in S/G2, anti-phase with Cyclin D1). Predicted cell-cycle oscillations in (I) MB with
CDK4/6 inhibitor — arrest — and (J) MB with CDK4/6i + Ezh2i — no rescue, because the CDK4/6 block
is downstream of Cyclin D1 transcription.

---

## 4. SUMMARY OF CHANGES FROM THE PUBLISHED v42 TEXT (for the authors)

1. **Cell-cycle core replaced.** The Gérard–Goldbeter relaxation oscillator (38 species, global
   ε = 150 time-scaling) is replaced by the real-time Heldt 2018 core with **explicit DNA
   replication**. S-phase now has a mechanistic, concentration-dependent duration. New modules:
   Cyclin B/CDK1 mitotic switch, intra-S CHK1 checkpoint, growth-gated restriction point, Skp2–p27
   feedforward. (Motivation for the rebuild: in the relaxation oscillator S-phase duration is a
   structural invariant, so Ezh2's dependence on time-in-S could not be modeled.)
2. **Ezh2 now integrates S-phase duration.** Ezh2 transcription is an E2f target gated on the
   S-window; stable Ezh2 protein, diluted at division, accumulates with S-phase length. This is the
   mechanistic basis for the hydroxyurea result (Fig. 4G).
3. **MB is modeled as reduced Ptch1 copy number + Mycn amplification with Shh present** (rather than
   Ptch1 = 0 / Shh = 0). The Simulation-conditions table changes accordingly.
4. **Period is no longer lengthened by the Ezh2 feedback.** In v42 the feedback lengthened the GNP
   period by 4.3 h; that was a consequence of the global clock. In v44, removing the feedback raises
   Cyclin D1 ~3-fold and abolishes its oscillation (this part is preserved) but changes the period
   only marginally, because phase durations are set by replication and growth, not by Cyclin D1
   level. **The corresponding sentence in the Results and the Supp. Fig. 8 caption have been
   rewritten** (the "lengthened the cell cycle period by 4.3 hours" and "0.9 h lengthening, 18.1 h
   vs. 17.3 h" claims are removed). This is the only qualitative change to a stated result.
5. **Hedgehog de-saturated; Cyclin D1 now matches the bulk RNA-seq.** The Ptch→Smo→Gli and
   Gli→Cyclin D1 terms were re-fit so Gli1 (6.30×) and Cyclin D1 (7.05×) rise in lockstep GNP→MB,
   identifying Cyclin D1 as predominantly Gli-driven. The Cyclin D1→Rb drive was made saturating
   (kPhRbCd·Cd/(K_CdRb+Cd)) so the ~7× transcript maps to a bounded, stable cell-cycle drive.
6. **The HHi rescue is reframed; the CDK4/6i contrast is preserved.** With data-matched Cyclin D1,
   MB+HHi sits ~3× above a cycling GNP (above threshold), so the vismodegib/Ezh2i proliferation
   effect is *not* a Cyclin D1-threshold switch — it requires additional Hedgehog biology, and we
   report this directly. The model still cleanly predicts that CDK4/6i arrests MB and Ezh2i cannot
   rescue it (downstream kinase block). The Ezh2 MB/GNP ratio remains a soft (cell-cycle-coupled)
   target. **The corresponding Results sentence has been rewritten** (the earlier "Ezh2i de-represses
   Cyclin D1 above the threshold and rescues cycling" claim is removed).
7. **Citation added:** Heldt et al. 2018 (BioModels BIOMD0000000700). The Gérard–Goldbeter 2009
   citation for the core engine is replaced by Heldt 2018 (Gérard–Goldbeter may be retained as the
   lineage/precursor if desired).
