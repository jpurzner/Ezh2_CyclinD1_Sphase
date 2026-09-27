# Inferring burst kinetics for GNP/MB from ChIP-seq + RNA-seq: feasibility

Proposal: take cell lines that have **both** transcriptional-kinetics data **and** ChIP-seq, learn
`bursting ~ f(ChIP, RNA-seq)`, then apply that to our own ChIP + RNA-seq to infer burst
characteristics for the model. Assessed below with the data in hand.

**Verdict: the premise checks out, and the problem is smaller than it looks.** Two results change
the design substantially — and one of them means a cheap route may get most of the way without
any ChIP model at all.

---

## 1. The premise is correct: burst frequency is set by cell state, burst size by the gene

Larsson Table S6 pairs 1609 orthologues with burst kinetics in **mouse and human fibroblasts** —
different species, *same* cell type. Against the earlier mouse fibroblast vs mouse ESC comparison
— same species, *different* cell type:

| | across **species**, same cell type | across **cell type**, same species |
|---|---|---|
| burst **frequency** | **+0.36** | +0.13 |
| burst **size** | +0.27 | **+0.50** |

The two go opposite ways, and it is exactly the split the mechanism predicts:

* **Burst size is gene-intrinsic.** It survives a change of cell type (r = 0.50) and degrades
  across a species boundary (r = 0.27) — consistent with it being partly promoter-sequence-encoded
  (Larsson's TATA effect, `protein_noise_bursting_chromatin.md` §3).
* **Burst frequency is cell-state-determined.** It survives a change of *species* (r = 0.36) but
  not a change of *cell type* (r = 0.13) — i.e. what sets it travels with the cell, not the gene.

**So the thing worth predicting from chromatin is burst frequency, and only burst frequency.**
Predicting burst size from ChIP is solving the easy half the hard way; borrow it instead.

## 2. The identity constraint halves the problem

For the telegraph model, mean = `kon·ksyn/(kon+koff)`. Checked against Larsson's own reported
means on 6180 ESC genes: **log–log r = 0.9999, median ratio 1.000.** The constraint is exact, not
approximate. (The cruder bursty limit `f·B` gives r = 0.988 and overestimates the mean by 13%;
97% of genes have `koff > kon`, 58% have `koff > 5·kon`.)

**Consequence: given mean expression from RNA-seq, you only need to predict *one* of {f, B}. The
other follows.** That removes a whole degree of freedom, and it means our own RNA-seq is not just
a covariate — it is half the answer.

## 3. What that implies — including a cheap route that may not need ChIP at all

Two ways to close the system, both using measured GNP/MB mean expression:

**Route A — borrow B, derive f.** Take burst size from mouse fibroblast or ESC (the
best-transferring quantity available, r = 0.50 across cell types), combine with measured mean, and
solve for frequency. **No ChIP model, no training, available today.**

**Route B — predict f from chromatin, derive B.** The proposed approach.

Route A is the honest baseline, and Route B has to beat it. On current evidence that is not
obvious: borrowing B carries r ≈ 0.50, whereas a ChIP model for f would have to do better than
that, and the only guide to how much variance ChIP explains is weak (single-antigen marginal
correlations of ~0.2, `protein_noise_bursting_chromatin.md` §2b–2c).

**But Route B has one thing Route A structurally cannot: condition-dependence.** A borrowed burst
size is a single static number. If the question is how bursting changes under **±vismodegib and
±EZH2i** — which is the actual biology of this model — only the chromatin route can answer it,
because only the chromatin changes between conditions. That, not genome-wide accuracy, is the case
for building it.

## 4. Training data that exists

| cell line | burst kinetics | ChIP depth | notes |
|---|---|---|---|
| **mouse ESC** | Larsson S3 (7309 genes); Ochiai 2020 | **very deep** — 469 distinct antigens in ChIP-Atlas alone (see below) | best-characterised pair |
| **GM12878** | Jin, Feng & Bush, bioRxiv 2022.09.21.508743 | **deepest that exists** — ENCODE Tier 1, hundreds of antigens | strongest candidate; human LCL |
| mouse fibroblast (primary CAST×C57) | Larsson S1 (7870 genes) | **thin** — NIH/3T3 is only a proxy and has few experiments | weakest link |
| human fibroblast | Larsson S6 (1609 orthologues) | moderate | gives the species contrast |

**mESC + GM12878 is the pair to train on** — both have deep ChIP, and they are about as different
as two cell types get, which is what makes the transfer test meaningful.

### 4b. Feature richness is far better than a single-antigen extrapolation suggests

An earlier pessimistic estimate here was extrapolated from **one** antigen (Suz12, marginal
r ≈ 0.2). That was the wrong basis. ChIP-Atlas alone holds **469 distinct antigens for mm10
`ES cells`**, including the marks that matter, with deep coverage:

| | experiments | | experiments |
|---|---|---|---|
| H3K27ac | 3182 | ATAC-seq | 3969 |
| **H3K27me3** | **2137** | RNA Pol II | 709 |
| H3K4me3 | 1099 | Ctcf | 420 |
| H3K9me3 | 595 | **Suz12** | 333 |
| H3K4me1 | 488 | Rnf2 (Ring1B) | 260 |
| H3K36me3 | 284 | **Ezh2** | 133 |
| H2AK119ub | 247 | Jarid2 | 102 |

GM12878 has ENCODE's full complement on top of this. A genome-wide model with 20–50 features is a
completely different proposition from a single-antigen correlation, and the earlier pessimism
should not be carried forward.

**Correction:** this document previously stated that H3K27me3 "cannot be tested" because
ChIP-Atlas excludes histone marks. That applies only to their precomputed *Target Genes*
convenience tables (`antigen class = "TFs and others"`). The **per-experiment peak BEDs and
bigWigs are available for histone marks**, so H3K27me3 — the mark v44 actually uses — can be used;
it just needs the promoter/gene-body aggregation step doing directly rather than downloaded
ready-made.

### 4c. This does not appear to have been published

Searching for a genome-wide predictive chromatin→burst-kinetics model turns up **inference**
methods (DeepTX, BurstLink — recovering bursting *from* scRNA-seq) and Ochiai's tail-selected
correlations, but not a trained, cross-cell-type-validated predictor. That cuts both ways: there
is **no published R² to set expectations against**, and the work would be novel rather than a
reimplementation.

## 5. The design — and the validation that actually matters

The point of "a couple of cell lines" is not more data. It is that **cross-cell-type validation is
the only test of whether the approach works at all.**

The question is not "can ChIP predict f?" — within one cell type it certainly can to some degree.
The question is whether the **mapping** `f = g(chromatin)` is conserved across cell types even
though `f` itself is not (r = 0.13). If the mapping is conserved, feeding GNP chromatin gives GNP
bursting. If the mapping is itself cell-type-specific, the method fails and no amount of training
data fixes it.

So:

1. **Features**, per gene: ChIP RPM at the **promoter (−1000/+100) and gene body separately** —
   Ochiai's split is the key methodological point, the two carry different information — plus mean
   expression, gene length, and TATA presence.
2. **Target**: log burst frequency (§1). Derive B from the identity (§2).
3. **Train on cell type A, test on cell type B**, and report the held-out-cell-type R². Within-cell-
   type cross-validation is not informative here and will flatter the model badly.
4. **Baseline to beat: Route A.** Report the improvement over borrowing B, not the raw R².
5. **Apply** to GNP/MB ChIP + RNA-seq, ideally per condition.

### Things that will bite

* **The expression confound.** ChIP signal, burst frequency and mean expression are all mutually
  correlated; mean expression must be in the model as a covariate or the ChIP coefficients will
  mostly be reading expression back (`protein_noise_bursting_chromatin.md` §2b).
* **Antigen availability must intersect.** The model can only *apply* marks present in both the
  training lines and GNP/MB. Two mitigations rather than restricting the model up front: train the
  rich version and **ablate** to our available marks to measure what is lost; and consider **ChIP
  imputation** (ChromImpute / Avocado / Epitome) to fill missing marks for GNP from the ones we
  have plus sequence — a mature method that exists precisely for this situation.
* **Train genome-wide; *evaluate* on the genes that matter.** These are separate concerns and an
  earlier version of this document conflated them, suggesting the training set be weighted toward
  genes resembling *Ccnd1* etc. **That was wrong** — it discards most of the ~6000-gene sample for
  no gain. Train on everything; then report held-out performance *stratified* by the regime our
  genes occupy (expression decile, PRC2-target status), and per-gene prediction intervals rather
  than a single global R². The stratification is the real concern: a model can be accurate overall
  and useless for lowly-expressed PRC2 targets, which is where several of our genes sit.
* **Larsson's counts are capture-deflated** (UMI capture p ≈ 0.21). Keep `B` and `⟨m⟩` in the same
  unit system throughout — see the unit hazard in `protein_noise_HANDOFF.md` §3.3.

## 6. The comparison worth making first

For a model with ~10 species, direct measurement may beat imputation on both cost and reliability.
**smFISH on *Ccnd1* in GNP** gives burst size and frequency for the single most important gene
directly, in the right cell type, with no transfer assumption. Against a genome-wide imputation
whose held-out-cell-type performance is unknown and whose feature set may be one histone mark,
that is a strong competitor.

The imputation earns its place if the goal is (a) covering the whole Hedgehog arm at once, or
(b) predicting how bursting **shifts between conditions** — which smFISH would require repeating
per condition per gene.

## 7. Recommended order

1. **Route A now.** Borrow B, take mean from our RNA-seq, derive f. Costs nothing, gives the model
   a defensible first parameterisation today, and is the baseline everything else is judged against.
2. **Test the transfer assumption** on public data alone: train `f ~ g(ChIP, expression)` on mESC,
   test on GM12878. If held-out-cell-type R² is near zero, stop — the method does not transfer and
   our own ChIP will not rescue it.
3. **Only if step 2 works**, apply to GNP/MB ChIP, prioritising the ±vismodegib / ±EZH2i contrast,
   which is the part Route A cannot do.

Step 2 is the decision point and it needs no new lab data — the burst kinetics and the ENCODE ChIP
are all public.
