# How the lattice gets parameterised from ChIP data (2026-10-23)

Not "training" in the ML sense: the lattice is stochastic and
non-differentiable, and one 24-cell ensemble costs ~8 s. There is no gradient
to descend and no fitting 21,600 genes in a loop. This is simulation-based
parameter inference via an emulator, structured around what the identifiability
analysis says is actually recoverable.

---

## Stage 0 — QC and gene selection

Filter samples on FRiP / enrichment first. Of four H3K27me3 samples testable
via peak calls, two were unusable (one over-called at 97% window coverage, one
with 260 targets). Then select genes with a RESOLVABLE boundary: a profile with
no edge carries no geometry, and fitting it returns noise dressed as a
parameter. Expect this to keep a few thousand PcG targets, not all 21,600.

## Stage 1 — Reduce profiles to features, not raw bins

Per gene per sample, extract the 7 features in `src/lattice/geometry_fit.py`
(`core_up/dn`, `rise_up/dn`, `plateau_up/dn`, `dip`). Cheap, vectorised, no
simulation. This turns ~150 bins into 7 numbers whose noise behaviour we can
estimate from replicates.

Fitting features rather than raw bins is deliberate: the bins are
heteroscedastic and strongly autocorrelated, so a least-squares fit over bins
would silently weight the flat plateau (many bins, little information) far above
the boundary (few bins, all the information).

## Stage 2 — Build an emulator of the lattice in feature space

The step that makes this tractable. Rather than calling the lattice inside an
optimiser:

1. Latin-hypercube design over the identifiable parameters only — from the SVD
   that is ~3-4 directions (kappa = a_rw/del_mk, k_h33, kdm_core_width x
   k_kdm6b), plus the per-gene recruitment terms.
2. Run the lattice once per design point (~300-500 runs, roughly an hour).
3. Fit a cheap regressor (GP or quadratic) mapping parameters -> the 7 features,
   carrying the ensemble noise as the observation model.

Then inversion is milliseconds per gene instead of hours, and the same emulator
serves every gene and every sample.

## Stage 3 — Hierarchical structure: share the chemistry, fit the recruitment

The design decision that makes per-gene inference possible at all.

**Kinetic constants are properties of the ENZYMES, not of individual
promoters.** They should be shared across all genes and samples:

  k_me1, k_me2, k_me3, del_mk, k_h33, a_rw_prc2, p_local, k_ac, g_ac

**What genuinely varies per gene is RECRUITMENT** — the spatial profiles:

  nucleation amplitude and width (CpG island), psi_KDM6B amplitude,
  transcription level, core K4me3 extent

This collapses the parameter count enormously and is biologically right. A
flat per-gene fit of all parameters would be both unidentifiable and wrong.

Partial pooling also handles the long tail of weakly-constrained genes
gracefully: they shrink toward the global posterior instead of returning
nonsense.

## Stage 4 — Cross-tissue data is leverage, not just volume

The same gene across tissues shares its sequence-determined geometry (the CpG
island does not move) but differs in Gli1/KDM6B, transcription and EZH2. That is
repeated measures with shared structure — exactly what identifies a hierarchical
model, and the reason the cross-tissue pull is worth more than the same number
of extra replicates in one tissue.

It also supplies the within-gene H3K27me3 -> TPM dose-response, which is the
step-2 gate identification that a cross-gene scatter at one condition cannot
give.

## Stage 5 — Validation

- held-out GENES (does the global chemistry transfer?)
- held-out TISSUES (does it transfer across contexts?)
- held-out PERTURBATION: the EZH2ko / EZH2i samples. This is the real test, and
  the surrogate work already shows the shape of the risk — a fitted mean-field
  model predicted held-out EZH2i at R2 0.925 but with a SYSTEMATIC bias,
  over-predicting mark retention by 58% at partial inhibition.

---

## Honest limits

1. **Normalisation quality is a hard ceiling.** Three of the seven features are
   levels (`plateau_up`, `plateau_dn`, `dip`) and they carry a large share of the
   sensitivity. Levels are only comparable across samples if normalisation is.
   This makes the `chipseq.norm` work a prerequisite for the cross-tissue step,
   not a parallel project. The length features (`core_*`, `rise_*`) are
   scale-free and will survive poor normalisation; the level features will not.

2. **ChIP occupancy is not me3 fraction.** The model's plateau is an absolute
   occupancy in [0,1]; ChIP gives relative enrichment. Without spike-in
   normalisation or matched mass-spec there is an unknown monotone transform
   between them, so absolute rate constants are not recoverable — only ratios
   (which is fine, since kappa is a ratio and is the best-determined direction).

3. **`plateau_up`'s noise floor came out at 0.001**, which is almost certainly
   saturation (it sits at 0.91, near ceiling) rather than real precision. Any
   SNR driven by that feature is inflated and should be discounted until the
   reference point is moved off the ceiling.

4. **Two parameters are not recoverable from me3 profiles at all** —
   `p_local` (sv 5.7) and `nucleation_gain` (degenerate, sv 2.9). Nucleation
   needs the EZH2/SUZ12 track. Kernel locality probably cannot be measured this
   way and should be fixed by assumption, with sensitivity reported.

5. **Biological vs ensemble noise.** The noise floor used in the identifiability
   analysis is simulation ensemble noise (28 cells). Real feature noise must be
   estimated from biological replicates and will be larger; the identifiable
   directions should be recomputed with that number before trusting the
   condition number of 50.7.
