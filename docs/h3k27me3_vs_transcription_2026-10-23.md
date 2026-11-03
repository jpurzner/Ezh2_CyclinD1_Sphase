# H3K27me3 pattern vs transcription — first analysis of the lattice dataset
2026-10-23. Data: `hpc5412@frontenac:/global/home/hpc5412/promoter_segmentation/lattice_dataset`
(47 cell types, TSS +/- 20 kb, 200 x 200 bp strand-oriented bins, paired `log2tpm`).
Analysis run server-side; this is mouse **GNP** (mm39), 77,455 genes, 16,080
me3-marked (20.8%).

---

## 1. The core is a SWITCH; the flank is a RHEOSTAT

me3-marked genes, by expression decile (expressed genes only):

| decile | log2tpm | **core me3** (±1 kb) | **flank me3** (2–5 kb) | ac core | k4 core |
|---|---|---|---|---|---|
| 1 | 0.00 | 1.026 | 0.679 | −0.106 | 0.072 |
| 2 | 0.02 | 0.804 | 0.486 | −0.045 | 0.129 |
| 3 | 0.05 | 0.854 | 0.458 | −0.012 | 0.245 |
| 4 | 0.12 | 0.826 | 0.422 | 0.021 | 0.354 |
| 5 | 0.25 | 0.900 | 0.429 | 0.059 | 0.479 |
| 6 | 0.51 | 0.948 | 0.420 | 0.085 | 0.555 |
| 7 | 0.97 | 1.045 | 0.412 | 0.107 | 0.707 |
| 8 | 1.74 | 0.920 | 0.348 | 0.211 | 0.941 |
| **9** | **2.98** | **0.458** | 0.275 | 0.440 | 1.333 |
| **10** | **4.49** | **−0.153** | 0.141 | 0.708 | 1.628 |

**Core me3 is flat (0.80–1.05) across deciles 1–8 — no relationship to expression
at all — then collapses in the top two deciles**, going below input in decile 10.
The threshold sits near **log2tpm ≈ 3 (TPM ≈ 8)**, the top ~20% of expressed genes.

**Flank me3 declines steadily and monotonically** (0.679 → 0.141) across the whole
range. So the two compartments encode different things: the flank tracks
expression continuously, the core switches.

## 2. It is NOT a nucleosome-depletion artefact

The obvious confound is that highly expressed promoters have an NDR, so *all*
histone signal falls there. Ruled out by the internal control: across the same
deciles **ac rises −0.106 → 0.708 and k4 rises 0.072 → 1.628** while me3 falls.
Active marks go up at exactly the sites where me3 goes down. The clearance is
specific.

Confirmed independently on the `state` layer, which is immune to signal scaling:

| core state | decile 1 → 10 |
|---|---|
| Polycomb | 0.26 → **0.05** |
| Active | 0.03 → **0.54** |
| Bivalent | 0.05 → 0.16 (peak at deciles 7–8) → 0.06 |
| Low | 0.66 → 0.35 |

**Bivalency peaks at intermediate expression and falls off at both ends** — a
transitional state, not an endpoint. That is Sneppen & Ringrose's prediction
(bivalent chromatin as an unstable background population) seen directly.

## 3. Correlation coefficients are blind to all of this

| feature | Pearson r with log2tpm | Spearman |
|---|---|---|
| total me3 (±10 kb) | −0.150 | −0.125 |
| flank me3 | −0.123 | −0.082 |
| core me3 | −0.166 | −0.073 |
| depletion (flank − core) | +0.123 | **−0.001** |

Every |r| <= 0.17, and the depletion Spearman is **exactly zero** — yet the decile
table shows a dramatic, highly structured relationship. The effect is
non-monotonic and threshold-like, so linear and rank correlation both miss it.
AUC for discriminating the top expression quintile: depletion 0.592, core 0.374,
total 0.389 — all weak, for the same reason.

**Any analysis of this dataset that reports correlation coefficients will
conclude there is no relationship. There is a very strong one.**

## 4. This contradicts our model in a specific, fixable way

Below the threshold the core is **me3-ENRICHED relative to the flanks** — the
core/flank ratio rises to **2.46** at moderate expression before inverting in the
top decile. Our particle lattice cannot produce that: its K4me3 block is a static
core mask, so the core is carved at every gene regardless of transcription.

The data say the core is where PRC2 nucleates (CpG island) and is therefore the
**most** marked region, until transcription crosses a threshold and clears it.

This is exactly the fix already flagged as highest priority in
`lattice_model_results_2026-10-22.md` §4 — **make K4me3 dynamic and tie it to
transcriptional state** rather than to a geometric mask. The dataset now demands
it independently, and supplies the target behaviour to fit: flat core me3 below
log2tpm ~3, sharp clearance above.

## 5. Next

- Repeat across the other 46 cell types: is the threshold at log2tpm ~3 universal
  or cell-type specific? The 20 GBM tumours and 20 ENCODE mouse tissues make this
  a real test of transferability.
- Fit the switch: what is the width of the transition in expression units? A sharp
  switch implies cooperativity, which is a direct constraint on the model.
- Use `me3`/`ac`/`k4` continuous layers to extract per-gene boundary geometry
  (core width, edge position) — the measurement the peak-call data could not give.
- GNP has extra marks available (H2AK119ub, H3K4me1, H3K36me3, Pol2, Gli1) if the
  three-mark lattice underfits.

## Caveats
- One cell type so far (GNP). Nothing here is yet shown to generalise.
- `log2tpm` is zero-inflated (26% of me3-marked genes are exactly 0); silent genes
  are analysed as their own bin rather than forced into a quantile.
- Direction of causation is not addressed: core clearance could cause expression
  or follow it. The EZH2ko samples would speak to this.


---

# ALL 47 CELL TYPES (added overnight 2026-10-23)

Figure: `simulations/fig_lattice_threshold_47cells.png`.
Per-cell curves: `lattice_allcells3.json` (on Frontenac, and copied to the
session scratchpad).

## Two QC problems had to be solved first

1. **The `state` layer is not comparable across samples.** Genome-wide Active
   fraction ranges **0.001 to 0.577** — where ac/k4 binarisation failed, every
   me3+ site collapses into "Polycomb" and nothing is Bivalent. The index file
   shows the cause: `cut_ac` ~0.005 and `cut_k4 = inf` in several samples.
2. **`me3_bin` is SATURATED for the ENCODE samples.** Core me3_bin plateau =
   **1.000** for 11 of 20 ENCODE mouse tissues — their me3 cutoffs are near or
   below zero (index `cut_me3` down to −2.035), so within me3-marked genes
   everything is called positive and the binary layer carries no information.

Both are avoided by using **continuous me3**, normalised within each sample as
`(core − min) / (plateau − min)`, which is invariant to both offset and scale and
so survives the SES-vs-log-ratio difference the README warns about. Samples are
kept only if the amplitude (plateau − min) exceeds 0.30; below that the fit is
noise.

## RESULT 1 — the threshold is universal as a RANK, not as an expression level

27 of 47 samples pass amplitude QC:

| group | n | pass QC | threshold percentile | transition width |
|---|---|---|---|---|
| GNP | 1 | 1 | 89.6 | 7.3 |
| ENCODE mouse | 20 | 4 | 85.5 (80.8–86.7) | 14.4 |
| ENCODE human | 6 | 2 | 64.7 (47.6–81.8) | 11.3 |
| **GBM** | 20 | **20** | **79.5 (73.2–89.6)** | **33.3** |

Pooled over all 27: **threshold = 80.3rd expression percentile, IQR 76.8–82.5,
range 47.6–89.6.**

The same samples in absolute units: **log2tpm median 1.75, range 0.04–3.99.**

    coefficient of variation:  percentile 0.10   vs   log2tpm 0.48

So the threshold is **~5x more consistent as a rank than as an absolute
expression level**, across two species, primary tissue, immortalised lines and
primary tumours. The GNP figure of log2tpm ~3 from the first pass is a GNP-specific
number; the transferable statement is **"the top ~20% of expressed genes"**.

## RESULT 2 — GBM transitions are significantly broader

| | n | median width (expression percentile points) |
|---|---|---|
| GBM | 20 | **33.3** |
| normal mouse tissue + GNP | 5 | **8.2** |

Mann–Whitney, GBM broader: **p = 0.0021.**

In normal tissue the switch completes within ~8 percentile points of expression;
in tumours it is spread over ~33. Two readings, and they are distinguishable:

- **Degraded switching per cell** — the boundary machinery is less sharp in tumour cells.
- **Population mixing** — a sharp per-cell switch measured over a heterogeneous
  bulk tumour necessarily looks broad.

GBM is notoriously heterogeneous, so the second is the more conservative reading,
and it is the same confound flagged throughout this project: bulk measurement
cannot distinguish a blurred boundary from a sharp boundary at varying positions.
Single-cell or single-molecule data separates them.

**Caveat:** GBM uses SES input-subtraction while the others use log-ratio or RPGC.
Width is measured on a normalised curve in expression-percentile units, so it is
largely protected, but the comparison is not between identical pipelines.

## RESULT 3 — amplitude differs greatly, and may be technical

| group | core me3 amplitude (plateau − min) |
|---|---|
| GNP | 1.092 |
| GBM | 0.762 (0.398–1.681) |
| ENCODE mouse | 0.162 (0.032–0.388) |
| ENCODE human | 0.179 (0.063–0.436) |

Only **4/20 ENCODE mouse and 2/6 ENCODE human** samples show enough core
clearance to fit at all, against **20/20 GBM**. Whether ENCODE tissues genuinely
lack the depleted-core architecture, or whether their shallower/differently
processed data cannot resolve it, is **not settled here** — the saturated
`me3_bin` layer in the same samples argues for a data-quality contribution.

## What this does to the model

The universality of the rank threshold is the useful constraint. A model whose
core clearance is driven by an absolute transcription rate would have to be
retuned per cell type; one driven by **relative** transcriptional rank transfers
unchanged. That favours a mechanism in which clearance depends on competition for
a limiting resource — which is what the particle model already assumes, since
PRC2 per domain is set by total abundance divided by the number of targets.

It also sharpens the dynamic-K4me3 fix: the K4me3 / transcription coupling should
saturate, producing a switch ~8 percentile points wide in normal tissue, and the
model should reproduce a broader apparent switch when averaged over a
heterogeneous population.
