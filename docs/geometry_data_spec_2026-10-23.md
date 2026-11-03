# Data spec: what the lattice needs to be parameterised per gene (2026-10-23)

For the server agent / whoever produces the profile-resolution ChIP. Written
against the identifiability analysis in
`simulations/geometry_identifiability.py` — i.e. this asks only for what
measurably constrains the model, and section 3 says what it will and will not buy.

---

## 1. The core request

**Binned H3K27me3 coverage around TSSs, per gene, per sample.**

| field | requirement | why |
|---|---|---|
| bin size | **100 bp preferred, 200 bp maximum** | the boundary rise width is 600-1500 bp; we need >=5 bins across it to fit `rise_up`/`rise_dn`, which carry the kernel and eraser information |
| window | **-6 kb to +9 kb** relative to TSS, asymmetric | the downstream domain runs to ~5 kb; the plateau must be measured BEYOND the boundary or it is not a plateau. A symmetric +/-2 kb window cannot see the edge |
| orientation | **strand-flipped** (minus-strand genes reversed) | up/down asymmetry is one of the strongest signals; mixing strands averages it away |
| granularity | **per gene**, not only a metagene | per-gene geometry is the point; a metagene is a useful fallback but cannot give per-gene N_core |
| normalisation | input-normalised, or depth+input consistent ACROSS samples | `plateau` is an identifiable feature, so its scale must be comparable between cell lines/tissues |
| format | genes x bins matrix (rows = gene, cols = bin offset), or bigWig + a TSS bed | either is fine; a matrix is less work for me |

The existing `alternateTSS` quantification does **not** satisfy this: all 21,599
regions are a single 500 bp window, so it has no spatial information at all.
A single 500 bp window is one bin — it cannot resolve a boundary.

## 2. Strongly desirable, same windows, same samples

| track | what it buys |
|---|---|
| **matched RNA-seq TPM** | step 2 of the parameterisation (the gate); without it this is geometry only |
| H3K4me3 | the core's active-mark profile; currently a static assumption in the model |
| H3K27ac | the p300 arm; the model currently predicts p300 contributes almost nothing to core depletion, which this would test |
| H3K36me2/3 | needed only if the downstream flank underfits (deferred Option C) |
| **KDM6B and EZH2/SUZ12 ChIP** | these become psi_KDM6B(x) and the nucleation profile DIRECTLY. `GNP_P7_KDM6B` and `GNP_P7_EZH2ko_Gli1-Flag` already exist in the ChIP repo — getting them binned on the same windows is probably the single highest-value addition |

## 3. What the data will and will not determine

From the sensitivity Jacobian (28-cell ensembles, 25% perturbations, scaled by
the ensemble noise floor). Condition number **50.7** — well posed, not degenerate.

**Best determined (will be pinned by a profile):**

1. **sv 146 — `a_rw_prc2` vs `del_mk` (weights -0.72 / +0.66).** This is the
   propagation-to-turnover ratio, i.e. Hodges & Crabtree's kappa. The profile
   determines it better than anything else, which is the quantitative version of
   the claim in the options doc that profile shape alone fixes the spatial
   parameters.
2. **sv 55 — transcription-coupled exchange `k_h33`**, separately identifiable.
3. **sv 37 — `kdm_core_width` x `k_kdm6b`**, the KDM6B extent and strength as a
   combination (not separable from each other, but jointly well constrained).

**Poorly determined (a profile will NOT give these):**

- **`p_local`, the kernel locality (sv 5.7).** Within the geometric-kernel
  family, an average profile barely constrains how local spreading is. Note this
  is different from the gross question — a 1/|i-j| contact kernel *is*
  distinguishable, because it destroys the boundary entirely. But fine-tuning
  locality from profile shape will not work.
- **`nucleation_gain` vs `del_mk` vs `p_local` (sv 2.9)** — mutually degenerate.
  Nucleation strength needs the EZH2/SUZ12 ChIP track, not the me3 profile.

Individually every parameter moves some feature by >15 sd, so detectability is
not the issue; correlation between parameters is.

## 4. A prediction the data can check immediately

With transcription-coupled H3.3 exchange acting only over the gene body, the
model predicts the **upstream flank carries MORE H3K27me3 than the downstream
flank** (plateau 0.91 vs 0.51 at the reference parameters), even though the
downstream domain is longer. A strand-flipped metagene shows this on the first
plot, before any fitting. If the real profiles are symmetric in plateau height,
the body-exchange term is wrong or much weaker than assumed.

## 5. Sample-quality warning from the peak-based attempt

Working from epic2 peak calls (the only spatial information currently
available), only 2 of 4 H3K27me3 samples were usable:

| sample | n PcG targets | coverage at TSS | verdict |
|---|---|---|---|
| EL_GNP1_H3K27me3_1 | 10069 | 0.213 | usable |
| EL_Tumor1_H3K27me3_1 | 11537 | 0.104 | usable |
| Scott_MB_H3K27me3_rep1 | 10348 | 0.969 | **over-called** — 97% of the window is "peak" |
| P7_H3K27me3_rep5 | 260 | 0.050 | **under-called** — only 260 targets |

Whatever pipeline produces the binned coverage, please include per-sample FRiP
or a similar enrichment metric so these can be filtered before fitting rather
than after. Geometry fitted on an over-called sample will report no boundary at
all.
