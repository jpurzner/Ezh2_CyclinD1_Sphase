# Developmental cell-cycle kinetics: does subthreshold mitogen lengthen G1 in vivo?

**Companion to:** `spencer_transient_G0_review.md` (Spencer/Meyer/Cappell transient-G0 literature, cell culture)
**Question driving this document:** the Ezh2/CyclinD1 model requires G1 lengthening to produce buffering under mitogen withdrawal, because Cyclin D1 responds rapidly to mitogen while baseline Ezh2 expression does not, so Ezh2 falls only if S-phase expression is delayed. Spencer-lab data are entirely cell culture. Do developmental systems, and GNPs specifically, show graded G1 lengthening under subthreshold mitogen?
**Short answer:** the phenomenology is reported, but every measurement that reports it derives G1 by subtraction from an extrapolated total cycle length under an unmeasured growth-fraction assumption. The reported G1 effect is exactly what a growing non-cycling fraction would counterfeit. See section 4.
**Papers read in full:** Contestabile 2009 (Brain Pathol), Legué 2015 (Development, pp. 1662-1668 only).
**Date:** 2026-07-20

---

## 1. The cell-cycle length hypothesis (cortex)

The governing framework in developmental neurobiology, and the reason "G1 lengthening causes differentiation" is treated as settled.

- **Calegari & Huttner 2003** (J Cell Sci 116:4947-4955). Olomoucine at 80 uM in whole mouse embryo culture (E9.5-E10.5) significantly lengthens but does not arrest the neuroepithelial cell cycle, and induces premature up-regulation of TIS21 and premature neurogenesis in telencephalon. Establishes G1 lengthening as a cause, not a consequence.
- **Calegari, Haubensak, Haffner, Huttner 2005** (J Neurosci 25:6533-6538). Cumulative BrdU in Tis21-GFP knock-in embryos. Progenitors undergoing neuron-generating divisions have a significantly longer cell cycle than those undergoing proliferative divisions, at both onset and advanced stages of telencephalic neurogenesis. Basal/SVZ progenitors have a longer G2 than apical ones.
- **Lange, Huttner, Calegari 2009** (Cell Stem Cell 5:320-331). Cdk4/cyclinD1 overexpression prevents G1 lengthening without affecting cell growth, cleavage plane or IKNM, delays neurogenesis, and expands basal progenitors. Conclusion: G1 lengthening is necessary and sufficient for the proliferation-to-neurogenesis switch.
- **Arai, Pulvers, Haffner, Schilling, Nusslein, Calegari, Huttner 2011** (Nat Commun 2:154). G1 lengthening tracks the apical-to-basal (fate-restricted) progenitor transition. Counterintuitively, expanding apical and basal progenitors have a substantially **longer** S phase than progenitors committed to neuron production.

**Note on numbers.** Absolute cortical phase durations across developmental stages (Takahashi et al. 1995 and successors) were not verified in this session and are deliberately omitted rather than quoted from memory. Retrieve before citing.

**Note on method.** All four use cumulative thymidine-analogue labelling or percent-labelled-mitoses, so all four inherit the growth-fraction problem in section 4. The two causal manipulations (olomoucine, Cdk4/cyclinD1) are pharmacological/genetic G1 manipulations, not mitogen manipulations, and therefore do not answer the subthreshold-mitogen question directly.

---

## 2. GNP baseline cell-cycle parameters (the usable numbers)

### Contestabile 2009, Table 2, P2 mouse, oEGL, control (euploid)

| Phase | Duration (h) | Fraction of Tc |
|---|---|---|
| **Tc** | **16.25 +/- 1.28** | 1.00 |
| **TG1** | **7.63 +/- 1.45** | 0.470 |
| **Ts** | **7.11 +/- 0.66** | 0.437 |
| **TG2** | **1.18 +/- 0.12** | 0.073 |
| **TM** | **0.33 +/- 0.04** | 0.020 |

n = 13 control mice. Values are mean +/- SEM.

### Independent cross-validation of the control set

- **Legué 2015, Fig 2B:** S-phase index of randomly sampled oEGL cells (30 min EdU pulse, P7) = **~40%, s.d. +/- 4%**, n = 18 clones, 412-5720 cells analysed per clone. Contestabile's Ts/Tc = 7.11/16.25 = **0.437**. Agreement to within one s.d.
- **Fujita 1967** (J Cell Biol 32:277-287) reports the same ~40% S-phase index; Legué cites it explicitly.
- **Fujita 1967; Mares & Lodin 1970** (Brain Res 23:343-352): Tc = **15-29 h from P1 to P10** by 3H-thymidine autoradiography. Contestabile's 16.25 h at P2 sits at the fast end, as expected for early postnatal.

**Conclusion: the control column is sound and can be used as the GNP parameter set.** The Ts65Dn column cannot (section 4).

### Shape of the GNP cycle vs. MCF10A

| | GNP (P2, control) | MCF10A (CDK2inc) |
|---|---|---|
| Tc / IMT | 16.25 h | 16-20 h |
| G1 (or M-to-S) | 7.63 h | 5-7 h (CV 0.26) |
| S | 7.11 h (43.7% of cycle) | shorter fraction |
| G2 | 1.18 h | longer |

Total cycle length is remarkably similar. The **distribution across phases is not**: GNPs spend 43.7% of the cycle in S with a very short G2 (1.18 h). For a replication-coupled dilution model this matters, because the dilution event is not brief relative to the cycle and the post-S restoration window (G2 + M = 1.51 h) is very short. Essentially all PRC2 restoration time in a GNP is G1 time.

---

## 3. Number of divisions and lineage timing

| Quantity | Value | Source |
|---|---|---|
| Descendants per GCP clone (median) | **250** | Espinosa & Luo 2008 |
| Average divisions per lineage | **~8** | Espinosa & Luo 2008 |
| E17.5 -> adult clone size (MADM) | 480 cells | Espinosa & Luo 2008 |
| E17.5 -> adult clone size (Atoh1-CreERT2) | **446.3** (n = 32) | Legué 2015, Fig 1E |
| E15.5 -> P7 clone size | 410.8 (n = 34) | Legué 2015, Fig 1E |
| E17.5 -> P7 clone size | 151.1 (n = 52) | Legué 2015, Fig 1E |
| P1 -> adult clone size | 118.6 (n = 44) | Legué 2015, Fig 1E |
| Spontaneous (no Tm) clone size | 71.3 (n = 11) | Legué 2015, Fig 1E |

Implied divisions, log2(clone size): E17.5->adult **8.8**; E15.5->P7 **8.7**; E17.5->P7 **7.2**; P1->adult **6.9**. Consistent with Espinosa & Luo's ~8 and with 2^8 = 256 vs. their median of 250.

**Total cycling time per lineage:** 8 divisions x 16.25 h = **~130 h (5.4 days)**. This is the window over which cumulative H3K27me3 dilution must be computed. Note that Tc is not constant across it, since average GCP proliferation slows progressively through development (Espinosa & Luo 2008).

Other lineage timing:
- gc transit time oEGL -> IGL: **~40 h** (Komuro & Yacubova 2003, via Legué 2015)
- Clone differentiates **en masse over a ~2-day window**, at a time chosen seemingly randomly between **E18.5 and P16.5** (Espinosa & Luo 2008)
- Clone size correlates with growth duration: Spearman **r = 0.83, P = 0.003** (Legué 2015, supp Fig S2D)
- Clone subtype homogeneity at P7: **<10%** of clones contained all three gc subtypes; **~80%** contained one subtype only or one representing **>75%** of cells (Legué 2015, Table S3)

---

## 4. The growth-fraction problem (the central methodological result)

### The measurement

Contestabile 2009 used an abridged cumulative BrdU protocol: 1, 2 or 3 BrdU injections (100 ug/g) spaced 2 h, sacrifice 2 h after the last, giving labelling times of **2, 4 and 6 h**. n = 6, 4, 3 per group per timepoint (13 control and 13 Ts65Dn total). A least-squares line y = a*t + b is fitted to three points, and the extrapolated x at y = 1 gives Tc - Ts, with the y-intercept giving Ts/Tc.

Their estimators:
- Tc = 1/a
- Ts = b/a
- TG2 = G2I x Tc, where G2I is the pHH3-punctate index
- TM = MI x Tc
- **TG1 = Tc - (Ts + TG2 + TM)**, by subtraction

Growth fraction is **assumed to be 1** ("100% for GCP located in the oEGL"). It is never measured, because the labelling curve never reaches plateau within 6 h. The extrapolation runs from 6 h out to 16.25 h (control) or 23.51 h (Ts65Dn), so 3 to 4 times beyond the data range.

### The algebra

With a true growth fraction GF < 1, the cumulative labelling index is

    LI(t) = GF * (Ts + t) / Tc

so the fitted parameters are a = GF/Tc and b = GF*Ts/Tc. Substituting into their estimators:

| Estimator | Result | Bias |
|---|---|---|
| Tc = 1/a | **Tc_true / GF** | **inflated by 1/GF** |
| Ts = b/a | Ts_true | **unbiased** |
| TG2 = G2I x Tc, where G2I = GF*TG2_true/Tc_true | TG2_true | **unbiased** (GF cancels) |
| TM = MI x Tc | TM_true | **unbiased** (GF cancels) |
| TG1 = Tc - (Ts+TG2+TM) | Tc_true/GF - (Ts+TG2+TM)_true | **absorbs the entire inflation** |

**Every phase estimated directly is unbiased. Only Tc is inflated, and because G1 is obtained by subtraction, the entire artifact lands in G1.**

### Sensitivity analysis on the Ts65Dn data

Ts65Dn measured: Tc 23.51, Ts 8.28, TG2 2.42, TM 0.44, so S+G2+M = **11.14 h**.

| Assumed GF (Ts65Dn) | Tc_true (h) | TG1_true (h) | vs. control TG1 = 7.63 h |
|---|---|---|---|
| 1.00 (their assumption) | 23.51 | 12.37 | +62% |
| 0.90 | 21.16 | 10.02 | +31% |
| **0.80** | **18.77** | **7.63** | **0% (identical)** |
| 0.69 | 16.22 | 5.08 | -33% |

**A non-cycling fraction of 20% in the Ts65Dn oEGL makes G1 identical between genotypes.** A non-cycling fraction of 31% removes the entire +45% Tc effect. Neither figure is implausible.

### Supporting evidence within the same paper for a non-cycling fraction

- oEGL BrdU labelling index (2 h, single pulse, Fig 2G): **~57% control vs. ~44% Ts65Dn**, ratio 0.77
- iEGL labelling index (Fig 2G): **~12.5% control vs. ~1.4% Ts65Dn**
- Total EGL BrdU+ cells (Fig 2H): **~440,000 control vs. ~190,000 Ts65Dn**
- Total pyknotic cells, P2 EGL (Table 1): **321.3 +/- 12.5 control vs. 1632 +/- 98.0 Ts65Dn**, a ~5-fold increase
- Total GCP number (Table 1): **2.203 +/- 0.496 x10^6 control vs. 1.246 +/- 0.251 x10^6 Ts65Dn** (-45%)

### Second, independent problem: steady state

Cumulative labelling assumes a closed population at steady state. The oEGL is not closed: cells continuously leave for the iEGL as they exit the cycle. See **arXiv:1905.12368**, "Cumulative labelling with thymidine analogues when the steady-state assumption is violated," for the formal treatment. Read before citing any cumulative-labelling Tc from EGL.

### What survives in Contestabile

**The G2 result is real.** It is supported by two assays that do not depend on Tc:
- At 2 h after a single BrdU pulse, **>60% of mitotic figures were BrdU+ in controls, and none in Ts65Dn** (Fig 6A,B). So TG2 + 1/2 TM is <2 h in control and >2 h in Ts65Dn.
- BrdU/pHH3 double label: all pHH3+ cells were BrdU+ in controls; in Ts65Dn only a few G2 cells and **no M-phase cells** were BrdU+ (Fig 6C,D).
- G2 index by pHH3 (Fig 5D): **10% of total GCP in Ts65Dn vs. 7% in control**; M index ~2% in both, no difference.

**The G1 result has no independent support anywhere in the paper.**

---

## 5. The CMV paper: the same artifact, independently generated

*JCI Insight* 2024, article 175525 (published version of bioRxiv 2022.10.02.510512). MCMV infection of P8 mouse GCPs.

- Tc increased by **+8.06 h**
- Ts increased by **+2.02 h**
- TG2+M determined separately by single-pulse BrdU with harvest at 1, 1.5 and 2 h, scoring BrdU+pHH3+ / total pHH3+. Maximum mitotic LI reached **within 2 h in both** groups.
- G1 derived by subtraction, reported as **+6.04 h** in one passage and **+3.5 h** in another. The discrepancy is in the published text and I have not resolved it; check the table before using either number.
- IdU max in iEGL at **32 h control vs. 40 h MCMV**, indicating delayed oEGL -> iEGL movement.
- The paper states that MCMV **increased cell-cycle exit**.

**This is the growth-fraction violation asserted in the same paper that assumes GF = 1.** Two of two studies reporting GNP G1 lengthening have this structure.

---

## 6. Legué 2015: clone-level heterogeneity and the S-phase index

Method: *Atoh1-CreERT2/+*;*R26R-LSL-YFP/+*, tamoxifen titrated to <25 clones per cerebellum. Spontaneous (Tm-independent) clones 1.8 per Cb, <13% of clones analysed. EdU 30 min before harvest.

### The key measurement (Fig 2)

| Population | S-phase index | Variation | n |
|---|---|---|---|
| Random oEGL cells surrounding clones | **~40%** | s.d. **+/- 4%** | 18 clones, 412-5720 cells each |
| **Pure oEGL clones** | **range 9.1% to 64.1%** | **s.d. +/- 15.2%** | 11 |
| Differentiating EGL clones | mean **46.9%** | +/- 5.1% | 7 |
| Pure oEGL clones (mean) | **34.6%** | | 11 |

- EdU incorporation did **not** correlate with clone size: **r = 0.41, P = 0.2** (Fig 2C).
- Clone-level s.d. is **~4x** the s.d. of randomly sampled cells from the same tissue.

### Legué's interpretation and an alternative

They reject strong cell-cycle-phase synchronisation within clones, reasoning that synchronised clones with S occupying a modest fraction of the cycle should mostly show a *low* S-phase index, whereas the observed indices were broadly distributed. They favour each clone having a distinct proliferation rate that changes over time.

**Alternative not considered in the paper:** the S-phase index of an asynchronous population is (Ts/Tc) x GF. A clone at 9.1% against a tissue background of 40% is either cycling ~4x more slowly **or** has roughly three-quarters of its cells out of cycle at that instant. Given that clones exit en masse over a ~2-day window at a stochastically chosen time, a clone caught partway through its exit window is precisely a clone with a large and rising non-cycling fraction. This is the best existing GNP evidence for the population structure the Spencer framework predicts, and it is read in the paper as a rate difference rather than a fraction difference.

### The result that argues directly against progressive G1 lengthening

Clones already containing differentiating iEGL cells had a **higher** S-phase index (**46.9%**) than pure oEGL clones (**34.6%**), with greater cell number. Legué interprets this, consistent with Espinosa & Luo at P9, as gcps **increasing** their proliferation rate during the final rounds of division.

**Progressive G1 lengthening as the mechanism of GNP cycle exit predicts the S-phase index falls as a clone approaches differentiation. It rises.** An abrupt binary exit from a cycle running at or above full speed does not have this problem. Note also that the S-phase index is unbiased with respect to growth fraction in a way that cumulative-labelling Tc is not, so this measurement is more trustworthy than the Contestabile G1 number it contradicts.

### Ancillary numbers (foliation, not load-bearing)

- Cb midline AP circumference increases **17.6-fold** between E17.5 and P14; at P14 AP is **7.8x** the ML axis, having been **half** the ML length at E17.5
- **87.9%** of clones (n = 149) more elongated along AP than ML
- AP/ML ratio distribution: <0.5 **3.4%**; 0.5-1 **8.7%**; 1-2 **31.5%**; >2 **56.4%**
- AP/ML >2: **80.8%** of EGL clones, ~50% of migrating clones, **29.7%** of IGL clones (anisotropy generated during proliferation, lost during differentiation)
- Vertical/horizontal anaphase ratio higher in sagittal than coronal sections, **P = 0.0019** (Fig 3F): division plane preferentially oriented along AP
- Nearest-neighbour dispersion: EGL clones **26.2 um**; migrating **32.4 um** (P = 0.11); IGL P7 **64.5 um**; IGL adult **77.2 um**
- Long vs. short folia: **178.8 vs. 101.4 cells** (P = 0.007); AP/ML **2.8 vs. 1.5** (P < 0.0001); AP length **1755 vs. 787 um**
- *En1+/-;En2-/-*: **197 vs. 446.3 cells** (P = 0.01); AP length **1430 vs. 2766 um**

---

## 7. Consolidated GNP parameter table

| Quantity | Value | Source | Confidence |
|---|---|---|---|
| Tc, P2 control | 16.25 +/- 1.28 h | Contestabile T2 | Good |
| TG1, P2 control | 7.63 +/- 1.45 h | Contestabile T2 | Moderate (by subtraction) |
| Ts, P2 control | 7.11 +/- 0.66 h | Contestabile T2 | Good (GF-unbiased) |
| TG2, P2 control | 1.18 +/- 0.12 h | Contestabile T2 | Good (GF-unbiased) |
| TM, P2 control | 0.33 +/- 0.04 h | Contestabile T2 | Good (GF-unbiased) |
| Ts/Tc | 0.437 | derived | Good |
| S-phase index, oEGL, P7 | ~40%, s.d. +/-4% | Legué F2B; Fujita 1967 | Good |
| S-phase index, clone-level spread | 9.1-64.1%, s.d. +/-15.2% | Legué F2 | Good |
| S-phase index, differentiating clones | 46.9% (vs 34.6% pure oEGL) | Legué F2 | Good |
| Tc range, P1-P10 | 15-29 h | Fujita 1967; Mares & Lodin 1970 | Moderate (3H-thymidine) |
| Divisions per lineage | ~8 (log2 of clone size: 6.9-8.8) | Espinosa & Luo; Legué F1E | Good |
| Descendants per GCP | median 250 (MADM 480 at E17.5->adult) | Espinosa & Luo 2008 | Good |
| Total cycling time per lineage | ~130 h (5.4 d) | derived | Moderate (Tc not constant) |
| oEGL -> IGL transit | ~40 h | via Legué | Moderate |
| Clone differentiation window | ~2 d, onset random E18.5-P16.5 | Espinosa & Luo 2008 | Good |
| **TG1, Ts65Dn** | **12.37 +/- 0.83 h** | Contestabile T2 | **DO NOT USE as fitting target** (section 4) |
| TG2, Ts65Dn | 2.42 +/- 0.12 h (+100%) | Contestabile T2, F6 | Good (independent assays) |
| Ts, Ts65Dn | 8.28 +/- 0.30 h (n.s.) | Contestabile T2 | Good |
| G2 index, control / Ts65Dn | 7% / 10% of total GCP | Contestabile F5D | Good |
| M index, control / Ts65Dn | ~2% / ~2% (n.s.) | Contestabile F5D | Good |

---

## 8. Contestabile's molecular result: the Overton loop, rediscovered

Cell-cycle microarray (112 genes, 81 expressed) plus RT-qPCR, P2 cerebellum. **Only two genes changed.**

| Gene | Array (Ts65Dn/ctrl) | RT-qPCR | Protein (western) |
|---|---|---|---|
| **Skp2** | **0.691 +/- 0.131*** | **0.725 +/- 0.125*** | **-41%** |
| **CcnB1** | **0.852 +/- 0.054*** | **0.780 +/- 0.083*** | **-35%** |
| Cdkn1a (p21) | 1.026 +/- 0.099 | 0.997 +/- 0.171 | not measured |
| Cdkn1b (p27) | 0.820 +/- 0.582 | 0.923 +/- 0.208 | **not measured** |
| Ccna2 | 0.922 +/- 0.172 | 0.945 +/- 0.219 | |
| **Ccnd1** | **1.193 +/- 0.381** | **0.934 +/- 0.120** | |
| **Ccnd2** | **1.040 +/- 0.172** | **0.947 +/- 0.148** | |
| App (trisomic control) | N/A | 1.512 +/- 0.254* | +46% |

GEO accession GSE9216.

### Two implications

**(a) Cyclin D1 and D2 transcripts are unchanged despite a documented Shh-responsiveness deficit** in this exact model (Roper 2006 shows a reduced but dose-dependent GCP mitogenic response to Shh, cell-autonomous). This is a second, independent instance of the tension already logged in `spencer_transient_G0_review.md` section 7:

- Spencer 2013 F4: cyclin D1 **protein** does not distinguish CDK2inc from CDK2low at the bifurcation
- Min 2020: mother-G2 cyclin D **does** control daughter fate, with translation rate as the integrator
- Contestabile T3: cyclin D1/D2 **transcript** unchanged under reduced Shh output in GNPs

Whatever cyclin D1 is doing, it is not reported by its own steady-state level, at either the transcript or the protein layer. **Modelling consequence: if Ezh2 loss perturbs cyclin D1, do not expect a changed cyclin D1 level as the readout. Expect a changed exit fraction with unchanged cyclin D1.**

**(b) Their proposed mechanism (Fig 8) is Skp2 -> p27 ubiquitination -> degradation, so less Skp2 means stabilised p27. This is the Overton 2014 double-negative loop, gene for gene:** CDK2 -> SCF/Skp2 -> degradation of the KIP pool -> de-repression of CDK2. Overton could not fit p21-null MCF10A with p21 alone and had to add p27 and p57 with identical kinetic constants; Contestabile arrives at the same module from the developmental side, in the tissue where p27/p57 are the dominant CKIs.

Note that **p27 mRNA is unchanged (0.923)**, which is what that model predicts, since the regulation is post-translational. **The missing experiment is a p27 protein western**, never run despite p27 accumulation being the proposed mechanism. Combined with Moser 2018 (p21 rise begins 5-10 h before anaphase, in mother G2) and Contestabile's G2 doubling, the specific prediction is: **p27 protein elevated in Ts65Dn GNP G2, with p27 mRNA flat.**

---

## 9. What this changes for the model

1. **Use the Contestabile control column as the GNP parameter set.** Tc 16.25 h, G1 7.63 h, S 7.11 h, G2 1.18 h, M 0.33 h. Do not use the Ts65Dn G1 value as a fitting target.

2. **Do not build buffering on graded G1 lengthening.** No GNP measurement supports it once the growth-fraction assumption is examined, and the one GF-unbiased measurement (Legué's clone S-phase index) points the other way: the S-phase index **rises** as a clone approaches differentiation (46.9% vs 34.6%).

3. **Build buffering on a variable-duration pre-G1 interval instead.** The CDK2low state satisfies the model requirement without stretching G1: its duration is variable and mitogen-sensitive, E2F is off during it (so Ezh2, an E2F target, is not transcribed, which is the "delayed S-phase expression" condition), and the population response is graded through the mixing fraction even though each cell is binary (Overton 2014: WT graded across 4-7 ng/mL EGF, p21-null switch-like).

4. **Almost all PRC2 restoration time in a GNP is G1 time.** G2 + M is only 1.51 h, and S is 43.7% of the cycle. Any restoration term needs to be evaluated against a 7.63 h window, not against a nominal "interphase."

5. **Cumulative dilution runs over ~8 divisions and ~130 h.** With Tc slowing across that window (Espinosa & Luo 2008), a fixed-Tc simulation will misestimate cumulative loss.

6. **New heterogeneity source: clone-synchronised exit with unsynchronised cell-cycle phase.** Legué shows clones differentiate en masse over ~2 days while being poorly synchronised in cycle phase (s.d. +/- 15.2% vs +/- 4% for random cells). If a differentiation trigger arrives across a clone within ~1 day while cells are scattered across a 16.25 h cycle, each cell is at a different point in the dilution-restoration cycle at exit. **H3K27me3 at cycle exit should be heterogeneous within a clone even with identical PRC2 kinetics per cell**, with the spread set by (restoration half-time / Tc). This predicts **within-clone** variance, distinct from the **between-clone** variance predicted by the G0-duration mechanism. Both are measurable, and they are separable.

---

## 10. Discriminating experiments

**Fixed tissue, feasible now:**

1. **Measure the plateau.** Extend cumulative EdU past the expected Tc (sampling to at least 36 h at P2, ideally 48 h) so the growth fraction is measured rather than assumed. If the plateau is below 1, and lower in the perturbed condition, the G1 claim collapses and the size of the non-cycling pool falls out directly. This is the single highest-value experiment in this document.
2. **Phospho-Rb (S807/811) IF on P7 EGL sections.** Per Moser 2018 the hyper-pRb fraction estimates the committed fraction directly from fixed tissue, without live imaging. Combined with Ki67 (Miller 2018: continuous decay in G0/G1 regardless of exit cause, so it clocks time since division), this gives branch fraction and residence time in situ.
3. **p27 protein western / IF in Ts65Dn or Shh-manipulated GNP G2.** Tests Contestabile's own proposed mechanism, which they left unmeasured. Predicted: p27 protein up, p27 mRNA flat.
4. **H3K27me3 IF on Atoh1-CreERT2 clones at P7**, partitioning variance within vs. between clones. Distinguishes the two heterogeneity mechanisms in item 6 above.

**Live imaging:**

5. **DHB CDK2 sensor in purified GNPs under graded Shh**, on the Nakashima 2015 dissociated-culture platform. Yields the CDK2low fraction and residence-time distribution as a function of Shh dose, which is the exact input the dilution module needs. Caveat: if an ERK KTR is used alongside, apply the CDK2-crosstalk correction (Spencer lab, Cell Systems 2025).

**Modelling, before any of the above:**

6. **Run the dilution model twice**, matched on mean inter-S interval: (a) uniformly stretched G1; (b) fixed G1 preceded by an exponentially distributed G0 with a mitogen-dependent hazard rate. If mean H3K27me3 matches but variance and single-cell distribution shape differ, the architectures are experimentally separable and experiment 2 above is the readout.

---

## 11. References

**Verified PMIDs**

| Paper | PMID | DOI |
|---|---|---|
| Contestabile et al. Brain Pathol 2009;19:224-237 | 18482164 | 10.1111/j.1750-3639.2008.00168.x |
| Legué, Riedel, Joyner. Development 2015;142:1661-1671 | 25834018 | 10.1242/dev.120287 (PMC4419279) |
| Roper et al. PNAS 2006;103:1452-1456 | 16432181 | 10.1073/pnas.0510750103 |
| Calegari & Huttner. J Cell Sci 2003;116:4947-4955 | 14625388 | 10.1242/jcs.00825 |
| Calegari, Haubensak, Haffner, Huttner. J Neurosci 2005;25:6533-6538 | 16014714 | 10.1523/JNEUROSCI.0778-05.2005 |
| Lange, Huttner, Calegari. Cell Stem Cell 2009;5:320-331 | 19733543 | 10.1016/j.stem.2009.05.026 |
| Arai et al. Nat Commun 2011;2:154 | 21224845 | 10.1038/ncomms1155 |

**DOI confirmed, PMID not verified**

| Paper | DOI |
|---|---|
| Espinosa JS, Luo L. J Neurosci 2008;28(10):2301-2312 | 10.1523/JNEUROSCI.5157-07.2008 |
| Nakashima K, Umeshima H, Kengaku M. Dev Dyn 2015;244(6):748-758 | 10.1002/dvdy.24276 |
| Mathematical model of granule cell generation, Bull Math Biol 2016 | 10.1007/s11538-016-0163-3 |

*(I previously wrote "PMID 18322077" for Espinosa & Luo from memory and did not confirm it. Do not use without checking.)*

**No verified identifier yet**

- Fujita S. J Cell Biol 1967;32:277-287. (Open on PMC; search title there for PDF and PMID.)
- Mares V, Lodin Z. Brain Res 1970;23:343-352.
- Koppel H, Lewis PD. Neuropathol Appl Neurobiol 1983;9:207-214.
- JCI Insight 2024, article 175525 (CMV lengthens GCP cell cycle). Published version of bioRxiv 2022.10.02.510512.
- arXiv:1905.12368, cumulative labelling when the steady-state assumption is violated. Preprint, no PMID by design.

**Worth pulling next**

- Locker M et al. Genes Dev 2006. Hh signalling shortens both G1 and G2 in developing Xenopus retina. Cleanest positive control for a morphogen tuning G1 length directly.
- Otsuki L, Brand AH. Science 2018;360:99-102. G2 quiescence in Drosophila neural stem cells. Counterexample to assuming developmental quiescence is always G1. Cited by Moser 2018.
- Yang R et al. Stem Cell Reports 2015;5:816-828. Cell division mode change in GNP neurogenesis controlled by Shh signalling.

---

## 12. Caveats and open items

- **Legué 2015 was read from pp. 1662-1668 only.** The abstract, opening Introduction, and pp. 1669-1671 (end of Discussion, Materials and Methods) were not available. EdU dose/timing details and clone-inclusion criteria sit partly in the missing section.
- Legué's EdU measurement is a **single 30-min pulse at one timepoint per clone**, so "each clone has a distinct proliferation rate that changes over time" is inferred cross-sectionally, not observed longitudinally.
- Legué's key comparison has **small n** (11 pure oEGL clones, 7 differentiating clones).
- Legué is **wild-type Swiss Webster with no mitogen perturbation**. It characterises baseline heterogeneity, not the response to subthreshold Shh.
- Contestabile's Ts65Dn is a **chronic, developmental, genetic reduction in Shh responsiveness**, not acute mitogen withdrawal. The causal chain "subthreshold Shh -> longer G1" is assembled across two labs (Roper measured the Shh response, Contestabile measured cycle length) and was never measured in the same cells.
- No study in the GNP literature separates a stretched G1 from a fixed G1 plus a variable G0, because cumulative labelling, EdU index, FUCCI (Cdt1 is high in G0 and G1 alike) and Ki67 all conflate the two.
- Cortical absolute phase durations (Takahashi 1995 and successors) were **not verified** in this session and are omitted above.
- The CMV paper reports **two different G1 values** (+6.04 h and +3.5 h) in different passages of the published text. Unresolved.
