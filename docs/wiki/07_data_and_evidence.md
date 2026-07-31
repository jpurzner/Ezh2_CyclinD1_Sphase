# Data & Evidence

*The experimental measurements and published literature that ground the model. This page is the "what is real" companion to [Calibration & Validation](06_calibration_and_validation.md) (which covers how those numbers are scored) and the [Model Architecture](02_model_architecture.md) (which covers what they constrain). Its organizing principle is the distinction the project fights hardest to preserve: **measured** vs **inherited** vs **assumed**, and **P7 GNP** vs **SHH-MB** cell type.*

**Primary source docs** (all under `/Users/jpurzner/Dropbox/Q_research/py_projects/Ezh2_CyclinD1_Sphase/`): `docs/DATA_PROVENANCE.md`, `docs/Ezh2_CcnD1_model_targets.md`, `docs/validation_targets.md`, `data/validation_targets.json`, `docs/H3K27me3_CyclinD1_literature_review.md`, `docs/GNP_developmental_cell_cycle.md`, `docs/hh-signalling-timing.md`, `docs/model_gene_transcript_map.md` (+`.csv`), `docs/cdk46i_resistance_research.md`, `docs/research_concentration_dependent_checkpoint.md`, `docs/scRNAseq_senescence_vs_transientG0_spec.md`, `docs/EVIDENCE_FRAMEWORK_PROPOSAL.md`, `data/AvgKdegs_genes_v1.1.csv`, `data/skp2_timecourse.json`.

---

## 1. The guiding principle: ratios, not absolutes

The model species are in **arbitrary units**. It is never fit to absolute concentrations, only to **experimental ratios and qualitative behaviors** (fold-changes, phase fractions, direction of response). The manuscript-facing one-liner (`DATA_PROVENANCE.md` §5):

> The signalling→CyclinD1 cascade (Hedgehog/Ptch1/Gli, MYCN, EZH2 induction and EZH2⊣CyclinD1 repression) and the CDK-inhibitor commitment threshold are **calibrated to Chahin et al. measurements**; the real-time **cell-cycle engine is adopted from Heldt et al. 2018**; the **H3K27me3 replicative-dilution kinetics from the chromatin literature**; only the mark's steady-state repression depth is tied to our CyclinD1 data.

### Provenance legend (used throughout `DATA_PROVENANCE.md` and `validation_targets.json`)

| Tag | Meaning |
|---|---|
| **MEASURED** | Constrained by Purzner-lab data (Chahin et al.): RNA-seq, qPCR, scRNA-seq, IF, microscopy, CKI panel |
| **INHERITED** | From external published models/data (Heldt 2018 engine; H3K27me3 chromatin literature). Not fit to our data |
| **SOFT** | A real target, loosely enforced / known miss (approximate literature value) |
| **FREE** | Not directly measured; set by constrained search or numerical stability (direction, not value, matters) |
| Cell-type tags | **P7 GNP** (postnatal-day-7 cerebellar granule-neuron progenitor) · **MB** (SHH-medulloblastoma culture) · **MB↔GNP** (cross-type ratio) |
| Extra prov. tags (JSON) | **manuscript** · **ext-lit** · **unpub** · **prediction** · **assum** |

---

## 2. The primary datasets (Chahin et al., Purzner lab)

| Data type | What it measures | Where in the model |
|---|---|---|
| **Bulk RNA-seq** | Gli1, Ezh2, Ccnd1, Ccnd2, Mycn, CKI panel across GNP / Ptch⁺ᐟ⁻ / MB ± vismodegib | steady-state fold targets |
| **qPCR** | CyclinD1 after Ezh2 cKO / Tazemetostat / Ezh2-OE / Palbociclib (Fig 3C,E; 4F) | EZH2⊣CyclinD1 repression depth |
| **scRNA-seq** | Ezh2 transcript by cell-cycle phase (Sepp, Vladoiu, Ocasio, Luo) (Fig 4A,B) | EZH2 = E2f-gated writer |
| **Immunofluorescence** | Ezh2 protein by phase; under HU / RO3306 / palbociclib / rShh (Fig 4G,H; S6D) | EZH2 S-phase integration |
| **Cell-culture microscopy** (NOT flow — a correction, see §7.4) | G0/G1/S/G2 proportions in MB ± arrest drugs (§F) | MB phase-fraction targets |
| **CDK-inhibitor panel** (DESeq2) | INK4 (p16/p18) + CIP/KIP (p21/p27) counts, GNP vs MB | commitment threshold |

The **GNP RNA-seq series** is a developmental Hedgehog-dose escalation: WT P7 GNP → Ptch⁺ᐟ⁻ P7 GNP → Ptch⁺ᐟ⁻ P28 GNP → Ptch⁺ᐟ⁻ MB.

---

## 3. Bulk RNA-seq expression ratios (the load-bearing folds)

From `Ezh2_CcnD1_model_targets.md` §A and `validation_targets.md`. These are the calibration backbone.

| Target | Value | Source | Status | Model |
|---|---|---|---|---|
| **CyclinD1 MB/GNP** | **5.07×** (p<0.05) | Fig 4I | MEASURED — *load-bearing* | 4.61 ✓ |
| CyclinD1 GNP+HHi / GNP | 0.14–0.157 (86% ↓) | §A | MEASURED | ✓ |
| CyclinD1 MB+HHi / MB | **0.144** (harness/provenance) vs **0.40** (compendium §A) | §A | MEASURED — *unreconciled conflict* | 0.11 ✓ |
| Gli1 MB/GNP | 6.9× | RNA-seq (no fig#) | MEASURED (pins Ptch1 copy-number f≈0.1) | 7.15 ✓ |
| Gli1 GNP+HHi | >99% ↓ | §A | MEASURED | ✓ |
| EZH2 MB/GNP | 2.05× (p<0.001) | Fig 4J | MEASURED / SOFT (model ~1.9–2.03) | 2.03 ✓ |
| EZH2 GNP+HHi | 25–47% ↓ | §A | MEASURED (EZH2 is Hh-responsive) | ✓ |
| **MYCN MB/GNP** | **2.80–2.86×** | §A | MEASURED | 2.71 ✓ |
| MYCN GNP+HHi / GNP | 0.78 (22% ↓) | §A | MEASURED (model ~0.91, SOFT miss) | ✓ |
| MYCN MB+HHi / MB | 0.86 (14% ↓) | §A | MEASURED | ✓ |
| **Skp2 MB/GNP** | **2.3×** | unpub (scRNA + MB table) | MEASURED (unpublished) | 2.34 ✓ |

**Data-derived design conclusion (`DATA_PROVENANCE.md` §2.1):** CyclinD1 tracks Gli1 ~1:1 along the Ptch escalation series → **CyclinD1 is Gli-driven, not MYCN-driven** (Gli supplies ~86% of MB CyclinD1). This is a structural choice forced by the data, not a free parameter.

### 3.1 MYCN — elevated *expression*, not amplification (correction, 2026-07-27)

Memory `mycn-elevated-expression-not-amplified`. **SHH-MB is NOT MYCN-amplified** (that is Group-3 MB); MYCN is high because it is a Hedgehog/Gli target. The JP-provided normalized counts verify the folds exactly: MB/GNP 8808/3077 = 2.86; GNP+vismo 2474/3077 = 0.80; MB+vismo 7862/8808 = 0.89. Developmental escalation 3077 (WT) → 4251 (Ptch⁺ᐟ⁻) → 8808 (MB). **Key nuance:** the acute vismodegib response *shrinks* as baseline Hh rises (GNP −20%, MB −11%) — a developmental **threshold-lock** into a vismodegib-resistant MYCN state. The acute self-activating latch to realize this proved **unmodelable** (buffer strength and GNP/MB discrimination are the same knob) → resolved by an additive `MYCN_expr` elevated-expression term replacing the old `MYCN_amplification` multiplier (algebraically neutral). See [Mechanisms](04_mechanisms_explored.md), [Abandoned](05_things_tried_and_abandoned.md).

### 3.2 CyclinD2 — the dominant, Hedgehog-buffered D-cyclin (two-cyclin model, baked 2026-07-28)

Memory `cyclind1-d2-two-cyclin-model`. The Chahin RNA-seq (normalized counts) shows **CCND2, not CCND1, is the dominant D-cyclin** and is far more Hh-buffered:

| | P7 GNP | MB | GNP+vismo | MB+vismo |
|---|---|---|---|---|
| Ccnd1 | 3639 | 27590 | 572 (**0.157**) | 3971 (**0.144**) |
| Ccnd2 | 18174 | 43389 | 12083 (**0.66**) | 25802 (**0.59**) |
| D1+D2 pool | 21813 | 70979 | 12655 (0.58) | 29773 (0.42) |

CCND2 = 83% of the GNP D-cyclin pool / 61% of MB, and drops only ~35% under vismodegib vs CCND1's ~85%. **Implication:** total D-cyclin driving CDK4/6 does **not** collapse under Hh-block (a D2 floor persists) → GNP cycle exit is not total-D-cyclin starvation, reinforcing EZH2/CyclinD1 as a threshold **governor** rather than an on/off switch. A separate `Cd2` species was added (default-neutral, then baked to **30/32** with D2 folds MB/GNP 2.39, GNP+HHi 0.66, MB+HHi 0.77). JP's belief that CCND2 is also H3K27me3-marked was tested and found **data-incompatible for a strong dynamic mark** (equal-marking flips the vismo direction: D2 rises under HHi; 27/32) — reconciled as D2's H3K27me3 being a static/developmental mark. See [Mechanisms](04_mechanisms_explored.md).

---

## 4. The EZH2 evidence (writer coupling + repression of CyclinD1)

### 4.1 EZH2 is a cell-cycle-coupled, mitogen-dosed E2F target

- **scRNA-seq, EZH2 transcript by phase** (Fig 4A GNP; 4B MB; 5 datasets): cycling/G0 ≈ **2.0×** → G0/cycling ≈ 0.5. Per-dataset G1/S-vs-G0: Sepp 1.87, Vladoiu 2.07, Ocasio(GNP) 2.09, Luo(MB) 1.79, Ocasio(MB) 2.45. Grounds the E2f×(CycE+CycA) transcription gate.
- **EZH2 protein by phase / arrest** (IF, Fig 4G, S6D): G2/G0 = **1.48×** (p=2.6e-5); HU S-phase/DMSO = **1.31×** (p=0.005; *alt-text says 1.22×* — a documented text conflict); RO3306 G2/DMSO 1.16× (ns); palbociclib **abolishes** the cycle-dependent EZH2 elevation → grounds the long EZH2 half-life + S/G2 accumulation + CDK4/6→Rb-E2F→EZH2 arm.
- **EZH2 integrates S-phase duration:** HU-arrested S EZH2 = 1.31× cycling S — the anchor for the "EZH2 = S-phase-duration sensor" mechanism.
- **Directional/qualitative** (load-bearing but non-scored): rShh dose → EZH2 ↑ (Fig 4H); E2f1-OE → EZH2 ↑ (Fig 4D, p<2.2e-16); palbociclib → EZH2 mRNA ↓ 55.6% + E2f1 ↓ 56.7% (Fig 4F qPCR).

### 4.2 EZH2 ⊣ CyclinD1 repression depth

| Target | Value | Source | Status |
|---|---|---|---|
| CyclinD1, Ezh2 cKO GNP / WT | log2FC 1.301 (**~2.46×**), −log10 p=16.0 | Fig 3A RNA-seq | MEASURED (genetic) |
| CyclinD1, GNP + Tazemetostat 1µM / DMSO | **2.2×** (p<0.05) | Fig 3C qPCR | MEASURED — *model 1.63, real shortfall* |
| CyclinD1, MB + Ezh2-OE / control | significant ↓ (fold n.r.) | Fig 3E qPCR | MEASURED (direction of repression) |
| CyclinD1, MB + Palbociclib / DMSO | +57% (1.57×, p<0.01) | Fig 4F qPCR | MEASURED (orthogonal de-repression) |

The recombination-corrected per-cell de-repression is ~2.8–3×; the model's baked EZH2i CyclinD1 fold ≈ 1.63 is a **genuine, confirmed shortfall** (the genetic cKO 2.46× argues it is real, not an artifact). Phase-fraction perturbations: GNP+Taze G1 2.03× / G2 1.44× (EZH2i *promotes* proliferation, Fig 3B); MB+Ezh2-OE G0 **2.84×** / G1 0.51× (EZH2 GOF → arrest, Fig 3F — the only up-direction test).

---

## 5. The CDK-inhibitor panel & commitment threshold

From `DATA_PROVENANCE.md` §2.6, DESeq2 CKI panel. Modeled as **abundance-proportional constants** (there is no CKI time-course to plot against — a stated limitation, `model_gene_transcript_map.md` #8):

| CKI | GNP → MB | Notes |
|---|---|---|
| p16 (Cdkn2a) | 0 → high (~**139×**, small absolute) | H3K27me3-silenced in GNP |
| p18 (Cdkn2c) | 0.464 → higher | constitutive INK4, **dominant** |
| p21 (Cdkn1a) | — (kSyP21 0.002 GNP / 0.004 MB) | 3.1× MB/GNP |
| p27 (Cdkn1b) | high at birth | 1.4× MB/GNP; dominant CIP/KIP |

**Key data-derived behavior:** vismodegib barely moves the CKIs (RNA-seq) → arrest is **CyclinD1 dropping across a static brake, not CKI induction**. The mapping caveat (`model_gene_transcript_map.md` #1): Heldt's `P21` species is literally *Cdkn1a*, but the v44 comments reinterpret the same pool as p27/*Cdkn1b* — one pool cannot map cleanly to one transcript.

---

## 6. Perturbation data: the vismodegib / EZH2i / CDK4/6i contrast

This is arguably the central phenomenon the model exists to explain (`EVIDENCE_FRAMEWORK_PROPOSAL.md` §1.3).

- **The rescuable arrest (Fig 5C, MB):** pRb⁺ cycling 100% → **~25%** under vismodegib → **~75%** restored by EZH2i (model 100/20/74). EZH2i lowers the mitogen bar → vismo-rescue. This is a manuscript figure but the raw N/cell-line is *unrecorded in the repo*.
- **The NON-rescuable contrast (CDK4/6i, palbociclib):** JP culture data — pRb⁺ cycling drops to **16–18%** (a partial arrest with a real residual, NOT 0), and is **reversible**: rebounds to **~50%** after 72 h + washout. MB+CDK4/6i+EZH2i "no rescue" is a **model prediction** (Fig S9I,J), *not* an experiment (a corrected mislabel). The harness reduces this to a coarse binary "0-div" proxy.
- **CDK4/6i mechanism (`cdk46i_resistance_research.md`, 2026-07-15 literature commission):** the pRb⁺ residual is **CDK2 (cyclin E/A–CDK2) re-phosphorylating Rb**, not incomplete drug engagement or alternative pocket proteins. pRb⁺ residual = diagnostic of CDK2/cyclin-E bypass, not RB1 loss (Wander 2020; Herrera-Abreu 2016; Asghar 2017). JP's pRb stain = phospho-Rb Ser807/811 clone D20B12 = a **hyperphosphorylation/committed-cycling** marker (maps to the model's `pRb` hyper species). Ser807/811 is phosphorylated by *both* CDK4/6 and CDK2 (only Ser780 is CDK4/6-exclusive) → the residual signature is CDK2 bypass. Direct SHH-MB validation: **Cook Sangar 2017** (PDX, PMID 28637687) — palbo → ~97% G1 arrest, ~63% tumor reduction, but 15/20 tumors recurred within 60 days of withdrawal = reversible cytostatic arrest + quiescent reserve. A true pRb⁺-dividing CCNE clone proved **unbuildable** in v44 (see [Abandoned](05_things_tried_and_abandoned.md)); the defensible story is reversible arrest + reserve.

---

## 7. Cell-cycle length & phase timing — the biggest data-hygiene story

### 7.1 GNP period: settled at ~16 h (the "22 h" is a miscitation)

The docs long carried "~22–23 h Nakashima" as the GNP period target. The 2026-07-25 audit (`validation_targets.md`, memory `validation-targets-formal-spec`) found this is a **miscitation**: the two real sources **agree on ~16 h**:

- **Nakashima 2015** (live-imaging lineage tracking): mean **15.9 h**, SD 3.71, CV 0.23 (windowed) / 0.32 (full).
- **Contestabile 2009** (Brain Pathol, PMID 18482164), P2 mouse oEGL control, Table 2 phase durations:

| Phase | Duration (h) | Fraction |
|---|---|---|
| **Tc** | **16.25 ± 1.28** | 1.00 |
| TG1 | 7.63 ± 1.45 | 0.470 |
| Ts | 7.11 ± 0.66 | 0.437 |
| TG2 | 1.18 ± 0.12 | 0.073 |
| TM | 0.33 ± 0.04 | 0.020 |

The model's **22.83 h is ~40% too slow** vs both; 22 h is literature-unsupported — it is purely the growth default `mu=0.0005` (`Tc = ln2/mu`). Adopting ~16 h *fails* the current model, so the harness change is deferred to a coordinated CKI-split recalibration (see the [16 h attempt](03_evolution_timeline.md) / memory `recalibration-16h-cycle-attempt`). Cross-validation of the control set: Legué 2015 oEGL S-phase index ~40% ±4% (P7 EdU); Fujita 1967 same ~40%; Fujita/Mares & Lodin Tc 15–29 h P1–P10.

### 7.2 The growth-fraction warning (a central methodological result)

`GNP_developmental_cell_cycle.md` §4: every developmental study reporting *graded G1 lengthening* derives G1 **by subtraction** from an extrapolated cumulative-BrdU Tc under an **unmeasured, assumed-1 growth fraction**. The algebra shows only Tc is inflated (by 1/GF) and the entire artifact lands in G1 (S, G2, M are GF-unbiased). A 20% non-cycling fraction makes Ts65Dn G1 identical to control. The one GF-unbiased measurement (Legué's clone S-phase index) points the *opposite* way: index **rises** (46.9% vs 34.6%) as clones approach differentiation. **Modeling consequence: do NOT build buffering on graded G1 lengthening** (a route that was in fact attempted and abandoned — memory `withdrawal-graded-g1-lengthening`, `gnp-g1-lengthening-buffering`; see [Abandoned](05_things_tried_and_abandoned.md)). ~8 divisions/lineage, ~130 h (5.4 d) cycling window (Espinosa & Luo 2008; Legué 2015).

### 7.3 MB period: no formal target yet

MB is longer (JP: ↑CKI + transient G0s), but there is **no formal MB cycle-length target** — an open item. The model currently sets MB/GNP ≈ 0.99 (the missing feature). Candidate mechanism: `k_mu_cki` couples Tc to INK4 tone (GNP 16.1 / MB 22.0); JP's approved reframe is MB = same ~16 h core + transient-G0 excursions, not a slower engine (memory `mb-celltype-transient-g0-parameterization`).

### 7.4 MB phase fractions (`Ezh2_CcnD1_model_targets.md` §F)

**These are cell-culture microscopy count-fractions, NOT flow cytometry** (a documented mislabel in two older docs). 24 h treatment, N=12 (Palbo N=11). DMSO mean ± SE: G0 21.7±1.59, G1 37.9±2.79, S 13.7±1.14, G2 14.1±1.42. Renormalized to 100% (÷0.874) → the validate targets 24.8/43.4/15.7/16.1. Critical caveats: compare to the model's **age-weighted count-fraction** (not duration-fraction); fit the resolvable **2N pool (G0+G1 = 68.2%)** (G0/G1 split needs a p27⁻/Ki67⁻ marker, treat as soft); S-duration anchor is cumulative-BrdU (~3 h), not the 15.7% gate; G2+M ≈ 2.5 h from direct methods (Fujita), not the 4N gate. Table F2 fold-vs-DMSO (SD, not SE — a caught error): HU S 1.359, HU G2 **0.233** (the perennial model miss ~0.36); Palbo S 0.295, Palbo G0 1.957; RO3306 S 1.584, G2 1.256.

---

## 8. CCND1 mRNA half-life — atlas data vs the model's reservoir

Memory `ccnd1-halflife-data-vs-reservoir` (2026-07-24). JP downloaded the **RNAdecayCafe atlas** (Vock; `data/AvgKdegs_genes_v1.1.csv`, ~122k gene×line rows, 11 human lines; nucleotide-recoding SLAM/TimeLapse-seq, EZbakR pipeline, no transcription inhibitor). Measured median half-lives (donorm / raw):

| Gene | t½ (donorm/raw) | Percentile |
|---|---|---|
| MYC / FOS / MYCN | 0.30 / 0.31 / 0.39 h | labile (~0th) |
| **CCND1** | **1.52 / 2.52 h** | **28th (moderately labile)** |
| GLI1 | 1.70 / 4.31 h | — |
| CCND2 / CCND3 / EZH2 | 3.56 / 2.74 / 2.53 h | — |
| ACTB / GAPDH (controls) | 6.42 / 12.95 h | 86th / 98th |

The controls land where expected → the atlas is trustworthy. The model's baked `k_Cd_mRNA_deg=0.001` (t½ ≈ 11 h reservoir) is **~5–7× too slow**; it predicts a 50% Ccnd1 drop ~12 h after Hh-withdrawal vs the data's ~2 h. **Crucial result:** even with a *fast* transcript (t½ 1.5–2.5 h), the model still coasts ~1 division after Hh-off — so the Ho-2020 ~11–12 h memory is the **committed Rb-E2F state**, not the transcript half-life. A fast transcript is validation-neutral (28/29). This **challenges** the "memory = slow CyclinD1 transcript reservoir" bake (memory `v44-cyclind1-hh-memory-carrier`) — decision pending. Caveats: transformed human lines (not GNP); gene-level isoform blend (full-length ~30 min vs 3′UTR-truncated hours); asynchronous. (Independent Hh-review value: Wiestner 2007 full-length CCND1 mRNA t½ ~30 min.)

---

## 9. Literature that grounds the inherited layers

### 9.1 The Heldt 2018 cell-cycle engine (INHERITED)

The real-time oscillator core — Rb–E2F bistable restriction point, cyclin/CDK network, Skp2–p27 commitment, explicit DNA replication, CyclinB/CDK1 mitotic switch, APC/C + division reset, growth/size gate — is **Heldt, Barr, Cooper, Bakal & Novák 2018 (PNAS), BIOMD0000000700**. Literature parameters, not fit to our data. See [Origins](01_origins_and_foundations.md). The `research_concentration_dependent_checkpoint.md` digest evaluated adding a Gérard-Goldbeter 2009 ATR/Chk1 concentration-dependent intra-S checkpoint (to lengthen S under HU without breaking oscillation) and flagged the Chao 2018 caveat that basal cell-cycle phases are statistically *uncoupled* — so the "longer S → more EZH2 → longer G1" coupling should be sought under HU stress, not basally.

### 9.2 The H3K27me3 / PRC2 chromatin literature (INHERITED)

`H3K27me3_CyclinD1_literature_review.md` is the full evidence base for the epigenetic module. Six executive conclusions: (1) the **complex, not the mark**, is the proximal repressor (readout should modulate initiation/burst frequency with a leaky floor); (2) repression on three timescales — mark removal (days) ≫ PRC2 (~12 h) > PRC1 (~hours), so **mark loss is rate-limiting**; (3) bidirectional coupling (nascent RNA evicts PRC2); (4) PRC2 is dynamic (~20% bound, seconds residence); (5) restoration slow/density-dependent (τ_restore ≈ 10–18 h); (6) absolute output anchorable (~2 mRNA/h average). Key kinetic anchors: EZH2i H3K27me3 t½ ~1 day, 90% by 3–4 d (Knutson 2012); dilution-driven derepression takes multiple divisions, proliferation-rate-dependent (Jadhav 2020); read-write allostery (Margueron 2009); EZH2 is an E2F target (Bracken 2003). *Ccnd1* locus facts: bivalent, non-CpG-methylated, large (read-write sufficient), poised Ser5P Pol II. **Open gap:** no direct measurement of H3K27me3 restoration kinetics at a large bivalent locus like *Ccnd1* — τ_restore is inferred from genome-wide/other-locus data. The ChIP direction (H3K27me3 on CyclinD1 is **~half in MB vs GNP**, i.e. MB<GNP) is JP's own ChIP observation, the sole source of the MB<GNP mark and the reason for the Gli-driven Jmjd3/Kdm6b eraser (memory `v44-gli-jmjd3-eraser-reconception`, `v44-four-mark-mechanisms`).

### 9.3 Hedgehog signaling timing (INHERITED / structural)

`hh-signalling-timing.md` (modeling-ready review). Bottom line: model Shh as a **slow-integrating, low-pass, cell-cycle-gated (cilium duty-cycle) mitogen**, NOT a frequency-encoded pulse train (that is ERK). Best-constrained: Gli1 protein t½ ~40 min; Smo t½ ~2 h; adaptation ~20 h (neural tube). GNP-specific: **~94% of outer-EGL progenitors are ciliated** (Ott 2024), cilia resorb pre-mitotically → a phase-gated input. **Central gap:** adaptation has never been quantified in GNPs (they proliferate for days under sustained Shh → may resemble the non-adapting regime). Downstream: Shh→Gli induces Mycn, Ccnd1, Ccnd2 (Kenney & Rowitch 2000); Shh withdrawal → irreversible arrest + differentiation (Wechsler-Reya & Scott 1999); a cell-intrinsic ~1-week timer overrides sustained mitogen. This review informs the "realistic mitogen input" work (memory `v44-realistic-mitogen-input`).

### 9.4 Skp2 developmental trajectory

`data/skp2_timecourse.json`: Skp2 log2FC 1.73, direction down across pseudo-time (E15→P56), t50 at P12; absolute table GNP_E15 2369 / P1 2524 / P7 1584 / MB 3605 — grounds the Skp2 MB/GNP 2.3× target (memory `v44-skp2-mitogen-dose`). Contestabile's molecular result (`GNP_developmental_cell_cycle.md` §8): in Ts65Dn only Skp2 (0.69×) and CcnB1 (0.85×) changed; **Ccnd1/Ccnd2 transcripts were unchanged despite a Shh-responsiveness deficit** — i.e. "whatever CyclinD1 is doing, it is not reported by its own steady-state level," predicting a changed exit fraction with unchanged CyclinD1 (a constraint that shaped the governor framing).

---

## 10. Distribution features & the senescence discriminator

**Distribution data** (`validation_targets.md`, mostly non-formal / proposed):
- GNP CyclinD1 G0/1 **CV = 0.70** (JP flow) — grounds the model σ=0.68 heterogeneity scale; predicted EZH2i→0.82 (EZH2 buffers variance) is a model *prediction* (memory `distribution-model-thread`).
- Transient-G0 fraction ~12% (MB; GNP outer-EGL ≈ 0), proposed.
- GNP cycle-duration is **lineage-correlated**: siblings σ 1.94 h (CV 0.12) ≪ non-siblings σ 5.11 h (CV 0.32) → most cell-to-cell variance is heritable/between-lineage, a real constraint on the population-of-models.

**Senescence vs transient-G0 (`scRNAseq_senescence_vs_transientG0_spec.md`):** a scRNA-seq analysis spec to test whether the model's *reversible transient-G0* is being mislabeled senescence. The critical rule: **pull all CDKN genes out of the senescence signature** (they are the model's brake — leaving them in makes the correlation circular). Prediction: high-CKI MB cells are still cycling and senescence-LOW; any true senescence is a small, non-cycling, SASP⁺, EZH2-low minority the model does not describe.

---

## 11. Measured vs assumed — the honest summary

**Data-constrained (MEASURED):** the entire signalling→CyclinD1 cascade folds (Gli1 6.9×, CyclinD1 5.07×, MYCN 2.86×, EZH2 2.05×, Skp2 2.3×, D2 folds), all ± vismodegib responses, the EZH2⊣CyclinD1 repression depth (cKO 2.46× / Taze 2.2×), EZH2 phase-coupling (scRNA 2.0×, IF 1.48×, HU-S 1.31×), the CKI panel set-points, MB phase fractions (2N pool 68.2%), and the vismo/EZH2i/CDK4/6i behavioral contrast.

**Inherited (external literature, NOT our data):** the whole Heldt engine + rate constants; all H3K27me3 kinetic constants (t½, τ_restore, dilution); Hh timing constants; GNP period (~16 h, ext-lit); G2+M ~2.5 h (Fujita).

**Assumed / free (direction matters, not value):** `k_Cd_translation` (CyclinD1 protein scale, numerical); CyclinD1 heterogeneity σ=0.68 (magnitude of the rescue); `M_commit` (G0/G1 split); the PRC2 eviction-arm strength (kept weak per the 5.07× fold); the 22 h period (growth-default artifact, literature-unsupported).

**Standing tensions / known misses** (`DATA_PROVENANCE.md` §4, `validation_targets.md`): EZH2 transcript S/G0 gradient too steep; MB HU S/G2 folds (structural anti-correlation, memory `v44-hu-s-g2-structural-limit`); MYCN GNP+HHi 0.91 vs 0.78; MB G2+M 3.6 h vs 2.5 h; EZH2 MB/GNP ~1.9 vs 2.05; **EZH2i fold 1.63 vs 2.2 (real shortfall)**.

---

## 12. The provenance problem (why this page matters)

`EVIDENCE_FRAMEWORK_PROPOSAL.md` (2026-07-25, proposal — nothing implemented) diagnoses six leaks in how evidence enters the repo. The most serious: **primary lab measurements enter through conversation** — the load-bearing 16–18% pRb⁺ palbo residual, the 100/25/75 rescue, CyclinD1 G0/1 CV 0.70, and Skp2 2.3× exist only as sentences in derived docs and chat transcripts, with no cell line, N, replicate structure, date, or raw-file pointer. Other leaks: the harness (`validate_v44.py`) scores literals while `validation_targets.json` claims to be the source of truth but *no code reads it* (drift already visible — period 22 vs 16, CyclinD1 MB+HHi 0.144 vs 0.40); ~40 author-year citation strings with no `.bib`/DOI/PMID; 45 undifferentiated `docs/` files mixing durable references with dated notes; results JSONs untraceable to a git SHA. The proposed fix (evidence records cited by ID, harness reads the target file, a cell-type-coherence linter) is the durable answer to the GNP↔MB conflation this project is most prone to. See [Reviewer & Open Questions](10_reviewer_and_open_questions.md).

---

*Cross-references: [Origins & Foundations](01_origins_and_foundations.md) · [Model Architecture](02_model_architecture.md) · [Evolution Timeline](03_evolution_timeline.md) · [Mechanisms Explored](04_mechanisms_explored.md) · [Things Tried & Abandoned](05_things_tried_and_abandoned.md) · [Calibration & Validation](06_calibration_and_validation.md) · [Figures Catalog](08_figures_catalog.md) · [Repos, Directories & Reproduction](09_repos_directories_reproduction.md) · [Reviewer & Open Questions](10_reviewer_and_open_questions.md)*
