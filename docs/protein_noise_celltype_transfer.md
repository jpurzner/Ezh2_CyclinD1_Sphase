# How cell-type-specific are these rates? Measured, not guessed

The concern: almost every prior in `turnover_burst_priors.csv` is borrowed from NIH3T3, mESC,
RPE1, H1299 or primary immune cells — none of them a cerebellar granule-neuron progenitor.

Rather than hunt for more datasets, this note **measures the transfer**, using the fact that
several datasets we already hold cover **multiple cell types with a single method**. That
separates cell-type variance from measurement variance, which is the question that actually
matters.

**Headline: the worry is real but misdirected.** Protein and mRNA half-lives transfer at
r ≈ 0.6 *no matter how similar the cell types are* — matching the cell type more carefully buys
almost nothing. The rate that genuinely does not transfer is **burst frequency** (r = 0.13), and
that one should not be borrowed at all.

---

## 1. The measurement

Every comparison below is log–log Pearson r on per-gene protein half-life.

| comparison | what varies | n | **r** | median fold-difference |
|---|---|---|---|---|
| **Dörrbaum** neuron-enriched vs glia-enriched | cell type only — same prep, method, species, lineage | 2738 | **0.63** | 1.8× |
| **McShane** mouse NIH3T3 vs human RPE | cell type **+ species**, same method | 1751 | **0.62** | 1.7× |
| **Mathieson** across 5 primary cell types | cell type only, same method | ~3000/pair | **0.57 – 0.79** | 3.9× |
| **Schwanhäusser vs McShane**, both NIH3T3 | **method only** — cell type held fixed | 2550 | **0.54** | 2.0× |

Read the table vertically. Neurons versus glia — same culture, same lab, same day, same lineage,
arguably the most similar pair obtainable — give **r = 0.63**. Mouse fibroblast versus human
epithelial, across a species boundary, give **r = 0.62**. And the same cell line measured by two
methods gives **r = 0.54**.

**The correlation is flat at ≈ 0.6 regardless of what you vary.** There is a reproducibility floor
that you hit immediately and cannot beat by matching cell type more carefully. Chasing a
better-matched dataset is not a route to a better number.

## 2. But the *level* does shift systematically

Rank transfers poorly-but-consistently; absolute level shifts in a structured way.

| | median protein t½ |
|---|---|
| Dörrbaum glia-enriched | 2.9 days |
| Dörrbaum neuron-enriched | 5.3 days |
| Mathieson mouse neurons | 46.6 h |
| Mathieson hepatocytes | 64.9 h |
| Mathieson monocytes | 66.1 h |
| Mathieson B cells | 78.3 h |
| Mathieson NK cells | 169.9 h |

A 3.6× spread across Mathieson's medians, and glia turn over ~2× faster than neurons in the same
culture. So the practical split is:

* **Borrow the ordering.** It transfers at r ≈ 0.6 and will not get better.
* **Do not borrow the level.** Rescale it, or rely on the fact that v44 is calibrated to *ratios*
  and not absolute concentrations (`docs/DATA_PROVENANCE.md`), which already absorbs a global
  scale factor.

## 3. Proliferation status may matter more than lineage

Glia (dividing) and neurons (post-mitotic) in the same culture differ ~2× in median half-life —
about as much as any cross-lineage comparison here. GNPs are **proliferating**.

That creates a tension worth being explicit about: the *lineage*-matched references (Mathieson
mouse neurons, Dörrbaum) are post-mitotic, while the *cycle-state*-matched references (NIH3T3,
RPE1, mESC) are the wrong lineage. **For a cycling progenitor the proliferating references are
probably the safer borrow**, and it is not obvious the neural ones are closer in any useful sense.

This also re-frames the dilution audit: post-mitotic references are dilution-free by construction
(the class-(b) anchors in `protein_noise_data_exploration.md` §1.5), which is a real advantage —
but it comes bundled with a cell state that differs from GNP in exactly the way that matters.

## 4. By quantity: what to trust

| quantity | cross-cell-type transfer | verdict |
|---|---|---|
| **protein half-life** | r ≈ 0.54–0.79, flat | borrow the rank, rescale the level. As uncertain as the measurement itself. |
| **mRNA half-life** | r ≈ 0.53 (Herzog mESC vs Schwanhäusser NIH3T3 — confounded with method) | same, with the caveat that no two mRNA datasets share a cell type, so this cannot yet be decomposed |
| **burst size `B`** | **r = 0.50** (Larsson fibroblast vs ESC) | borrow with caution; partly promoter-encoded (TATA), so there is a mechanistic reason it travels |
| **burst frequency `f`** | **r = 0.13** | **do not borrow.** Enhancer/chromatin-driven, no sequence determinant (R² = 0.008), does not survive a change of cell type |
| **cycle-phasing of degradation** | unknown — only Battich, one cell type | assume phase is conserved, fold-range is not |

**The single actionable conclusion: burst frequency is the parameter to worry about, and it is the
one your own ChIP-seq could actually constrain** (`protein_noise_bursting_chromatin.md` §6). The
half-lives are fine to borrow at the precision they are measured to.

## 5. What would genuinely help, and what would not

**Would not help much:** more protein half-life datasets in more cell types. The r ≈ 0.6 floor is
already reached and another cell type will not raise it.

**Would help:**

1. **Any proliferating neural progenitor turnover dataset.** Nothing found covers this — every
   neural dataset is post-mitotic (Mathieson neurons, Dörrbaum, iPSC-derived neurons, dopaminergic
   neurons). The gap is cycle state, not lineage.
2. **Burst kinetics in a second mouse cell type beyond fibroblast/ESC**, to test whether
   r = 0.13 is a floor or an artefact of that particular pair. Candidate:
   *Allelic transcription factor binding shapes transcriptional kinetics in human cell lines*
   (bioRxiv 2022.09.21.508743) — multiple human lines, allele-resolved.
3. **Two mRNA half-life datasets sharing a cell type**, so the method/cell-type decomposition can
   be done for mRNA as it has been for protein. We currently cannot tell whether the Herzog vs
   Schwanhäusser r = 0.53 is cell type or platform.

None of these is blocking. The honest position is that borrowed half-lives carry roughly 2-fold
per-gene uncertainty whatever you do, and that should be propagated as a sensitivity band rather
than chased.

---

## Sources added for this note

* **Dörrbaum et al. 2018** *eLife* 7:e34202, Suppl. Table 1 (PMC6008053, OA) — dynamic SILAC
  half-lives in rat primary hippocampal culture, split into mixed / neuron-enriched /
  glia-enriched. `data/priors/raw/dorrbaum2018/`.
