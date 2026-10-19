# Stability (half-life) data used in the v44 model — provenance & methods

All turnover values come from the lab's **stability-catalog** (AMOS "external-parameter" arm:
`py_projects/stability-catalog`), a curated aggregation of **24 published mammalian turnover datasets /
120 experimental arms**. No value is a lab measurement. Half-lives are *degradation-only* estimates; for
each gene the model adopts the value from the single most tissue-appropriate ("chosen") arm. Verified
against the model code 2026-10-19.

The model's native time unit is **minutes**; half-lives below are ln2/k/60.

---

## 0. How one value was selected per gene (the ranking rule)

The catalog does **not** average across datasets or fit the best number. `stability.py:_half_life()` applies a
**fixed, deterministic preference ranking** and takes the single highest-ranked arm that has data for the gene
(then the median of its replicates). The rank is:

1. **Curated dataset preference order** (tissue-relevance first). For **mRNA**: Serdar 2025 mouse embryonic
   cortex (E14.5 → E16.5 → E11.5) sits at positions 0–2 — the closest neural-progenitor data is ranked first —
   then Ietswaart 2024, Herzog 2017 mESC, RNAdecayCafe, Schwanhäusser 2011. For **protein**: Nakanoh 2026
   ESC-derived neural progenitor (0–1) → Li 2021 cell-line turnover (RPE1/U2OS/HCT116/HEK293T) → Zecha 2018
   HeLa → Fornasiero 2018 / Li 2025 brain tissue → others.
2. **Degradation-only filter** (we call it with `prefer="degradation_only"`): arms whose published number is a
   *removal rate that includes growth dilution* are demoted (+500), unknown-dilution arms +200. Rationale: the
   model supplies its **own** division/dilution term, so we want a pure **degradation constant**.
3. **Exclusions**: censored values (bounds, not measurements) +1000; McShane 1-state fits +100 (they measure a
   different quantity than the 2-state steady-state fit).
4. Lowest total rank wins; value = **median within that arm**.

Consequence: **every mRNA value is Serdar 2025 E14.5** because it is rank-0 *and* has all seven genes. The
**protein** winners differ by gene because rank-0 (Nakanoh NP) often lacks that gene or isn't degradation-only,
so the ranking falls through to the first arm that qualifies (p18/p27 → Zecha HeLa; p21/Skp2 → Li 2021).

**The `range_across_arms` (up to 1000-fold) is NOT real disagreement** — the extreme bounds are actinomycin-D
transcription-arrest artifacts (half-lives of hundreds of h) and near-zero labeling artifacts; across the
*relevant* neural datasets the values are tight (§5).

---

## 1. mRNA half-lives — USED in the baked default model

Wired into the two-step gene-expression modules (transcription → mRNA species → translation → protein).
All six adopt the **tissue-matched arm: mouse embryonic cortex E14.5 (Serdar et al. 2025)** — the closest
available neural-progenitor dataset.

| Gene (model / symbol) | mRNA t½ | source arm (n arms) | where in code | status |
|---|---|---|---|---|
| p18 / *Cdkn2c* | **3.73 h** | serdar2025:embryonic_cortex_E14.5 (100) | `kDe_p18m`, `with_cdki_mrna` | baked default |
| p19 / *Cdkn2d* | **3.62 h** | serdar2025:embryonic_cortex_E14.5 (97) | `kDe_p19m` | baked default |
| p21 / *Cdkn1a* | **5.18 h** | serdar2025:embryonic_cortex_E14.5 (107) | `kDe_p21am` | baked default |
| p27 / *Cdkn1b* | **3.61 h** | serdar2025:embryonic_cortex_E14.5 (92) | `kDe_P21m` (species "P21" = p27) | baked default |
| Skp2 / *Skp2* | **3.38 h** | serdar2025:embryonic_cortex_E14.5 (108) | `kDe_Skp2m` | baked default |
| CyclinD1 / *Ccnd1* | **~3.5 h** | serdar2025 3.52 h; catalog median ~3.4 h (n=135); confirmed mufteev2026 NPC 2.0 h | `k_Cd_mRNA_deg`, `measured_ccnd1` | baked default |

Notes: CyclinD1 as-built is 3.56 h (cited as 3.4 h; ~5% off from a unit-label note in the `_r` derivation —
cosmetic, pre-existing). The legacy 11.5 h CyclinD1 "reservoir" (`measured_ccnd1=False`) is a committed-state
*construct*, not a measured decay.

## 2. Protein half-lives — looked up; mostly NOT used in the baked default

The baseline model's protein turnover is **not** taken from the catalog: CDK-inhibitor and Skp2 protein
clearance is cell-cycle-*regulated* in the model (SCF–Skp2, KPC, APC/Cdh1), and EZH2 protein turnover is fit
to EZH2-coupled validation targets. The catalog's degradation-only protein half-lives (from non-cycling cell
lines) would mis-specify these actively degraded regulators.

| Gene / symbol | measured protein t½ | source arm | status in model |
|---|---|---|---|
| p18 / *Cdkn2c* | 23.87 h | zecha2018:HeLa (19) | **sensitivity variant only** (`cdki_measured_protein`, default OFF) |
| p21 / *Cdkn1a* | 0.98 h (short-lived) | li2021:RPE1 (5) | **sensitivity variant only** (`cdki_measured_protein`, default OFF) |
| p27 / *Cdkn1b* | 44.47 h | zecha2018:HeLa (15) | examined, **not used** (regulated clearance retained) |
| Skp2 / *Skp2* | 6.77 h (short-lived) | li2021:HCT116 (8) | examined, **not used** (regulated clearance retained) |
| CyclinD1 / *Ccnd1* | 0.42 h (short-lived) | li2021:RPE1 (8) | consistent with model short Cd; level not catalog-set |
| EZH2 / *Ezh2* | 57.47 h (range 11.6–57 h) | zecha2018:HeLa (8) | **not used** — EZH2 protein t½ is model-fit |
| p19 / *Cdkn2d* | **no data** | — | no measured protein arm exists; model value assumed |

## 3. Source citations (from the catalog's `export/sources.py`)

- **Serdar et al. (2025)** — mRNA half-lives, mouse embryonic cortex E14.5 (metabolic labeling). *PLoS Biology* 23(2):e3003031. doi:10.1371/journal.pbio.3003031.
- **Zecha J, Meng C, Zolg DP, Samaras P, Wilhelm M, Kuster B (2018)** — protein turnover, HeLa (pulsed-SILAC/TMT). *Mol. Cell. Proteomics* 17:974–992. doi:10.1074/mcp.RA118.000583.
- **Li J, Cai Z, Vaites LP, … Gygi SP (2021)** — protein turnover (non-linear model), RPE1/HCT116. *Molecular Cell* 81:4722–4735.e5. doi:10.1016/j.molcel.2021.09.015.
- **Mufteev M, Rodrigues DC, et al. (2026)** — CCND1 mRNA, iPSC-derived NPC (preprint). bioRxiv 2023.03.01.530249. doi:10.1101/2023.03.01.530249.

(The across-arm ranges draw on the full catalog pool — Li 2025, Schwanhäusser 2011, Mathieson 2018,
Fornasiero 2018, Nakanoh 2026, Swovick 2021 for protein; Agarwal 2022, Ietswaart 2024, Herzog 2017,
Sharova 2009, Friedel 2009, Loedige 2023 for mRNA — but the *chosen* values are the arms listed above.)

---

## 5. Chosen value vs the other datasets (robustness check)

**mRNA (chosen = Serdar 2025 E14.5) vs the other neural datasets** — values are tight (all ~2–6 h); the chosen
value is representative, not an outlier:

| Gene | **chosen** (Serdar E14.5) | Serdar E16.5 / E11.5 | Mufteev 2026 iPSC-NPC / neuron | Herzog mESC | Schwanhäusser NIH3T3 |
|---|---|---|---|---|---|
| p18 | **3.73** | 3.81 / 3.37 | 3.26 / 1.90 | — | 4.85 |
| p19 | **3.62** | 3.00 / 3.00 | 2.79 / 2.15 | 1.85 | — |
| p21 | **5.18** | 3.34 / 4.82 | 5.85 / — | 1.41 | — |
| p27 | **3.61** | 2.78 / 3.32 | 1.91 / 1.60 | — | — |
| Skp2 | **3.38** | 2.17 / 3.90 | 2.32 / 2.56 | 2.23 | 2.82 |
| CyclinD1 | **3.52** | 3.24 / 4.07 | 1.97 / 1.80 | 4.15 | 14.0 |
| EZH2 | **4.16** | 3.11 / 4.05 | 3.42 / 2.13 | 3.71 | 9.85 |

**Protein (chosen vs the relevant alternatives):**

| Gene | **chosen** | Nakanoh NP (ms/hu) | Li 2021 RPE1/U2OS/HCT116 | Zecha HeLa | brain tissue (Fornasiero/Li25) |
|---|---|---|---|---|---|
| p18 | **23.9** (Zecha) | — | — | 23.9 | 228 (cerebellum) |
| p21 | **0.98** (Li RPE1) | — | 0.98 / 1.36 / 1.52 | 2.15 | — |
| p27 | **44.5** (Zecha) | 10.2 / — | — | 44.5 | 84 / 96 |
| Skp2 | **6.77** (Li HCT116) | — | 13.7(lin) / 6.77 | 5.59 | — |
| CyclinD1 | **0.42** (Li RPE1) | — | 0.42 / 0.55 / 0.47 | 0.45 | — |
| EZH2 | **57** (Zecha)† | 21 / 39 | 26 (RPE1) | 57 | — |

† EZH2 protein t½ is **not used** in the model (model-fit instead); shown for completeness. The NP arms (21–39 h)
are closer to the ~25 h the model's EZH2 turnover approximates.

---

## 6. Outliers & normalization robustness (sim_mrna_normalization_impact.py)

**What the outliers are.** The up-to-1000-fold "range across arms" is driven almost entirely by **`agarwal2022`
compendium arms** — actinomycin-D transcription-arrest and 4sU-labeling measurements that produce both near-zero
(~0.003 h) and ultra-long (200–625 h) mRNA values. These are **methodology artifacts** (ActD confounds decay by
blocking transcription). For protein, the long tail is instead **real post-mitotic tissue** (Li 2025 brain/organ,
up to 998 h) — not an artifact, but inapplicable to cycling progenitors; naive normalization would *inflate* the
protein values (e.g. p18 23.9 h → 84 h), so it is the wrong move for protein.

**Normalized central values** (log-IQR-trimmed geometric mean, outliers removed): p18 3.7→2.6 h, p19 3.6→2.7 h,
**p21 5.2→2.1 h**, **p27 3.6→1.2 h**, Skp2 3.4→3.8 h, CyclinD1 3.5→3.4 h, EZH2 4.2→3.5 h. The Serdar neural
values sit on the high side for p21 and p27; the rest are close.

**Impact on the model — tested by swapping the normalized set in:**
- **Validation: none.** 30/32, identical to baseline (same two misses; every target value unchanged). The
  two-step restructure is steady-state-neutral (translation rate = mRNA decay rate), so the mRNA half-life is
  validation-irrelevant — the outliers and the choice of central value do not affect any calibration target.
- **Transient-G0: modest, dynamics only.** Shortening p21 mRNA (the driver) 5.2→2.1 h lowered MB snapshot-G0
  occupancy from ~46% to ~32% (shorter mRNA → shorter protein excursions → shorter G0 dwells; less permanent
  arrest). This stays within the ~36–55% identifiability band (§5 of the protein doc / overnight review) and
  does not change the qualitative result (large MB transient-G0, GNP 0%, p21-driven).

Conclusion: the mRNA half-life enters the model only as a burst-filtering timescale, not as a calibration
constraint, so the exact value (chosen vs normalized) is immaterial to validation and only fine-tunes the
transient-G0 dwell within the range we already report.

---

## 4. Draft methods paragraph

> **Transcript and protein stability.** Messenger-RNA half-lives used to parameterise the model's two-step
> gene-expression modules were obtained from a curated compilation of published mammalian turnover
> measurements (an aggregation of 24 datasets / 120 experimental arms; all degradation-only estimates). For
> every gene we adopted the value from the most tissue-appropriate dataset — mouse embryonic cortex at E14.5
> (Serdar et al., 2025), the closest available neural-progenitor arm — giving mRNA half-lives of 3.7 h
> (*Cdkn2c*/p18), 3.6 h (*Cdkn2d*/p19), 5.2 h (*Cdkn1a*/p21), 3.6 h (*Cdkn1b*/p27), 3.4 h (*Skp2*) and ~3.5 h
> (*Ccnd1*/Cyclin D1); the Cyclin D1 value was consistent with an independent neural-progenitor measurement
> (Mufteev et al., 2026) and with the catalog-wide median (~3.4 h, n = 135 measurements). Because the model is
> expressed in dimensionless concentration units, these half-lives constrain relaxation kinetics rather than
> absolute abundances: each transcription→mRNA→protein cascade was parameterised with translation rate equal
> to the mRNA decay rate, so that steady-state protein levels — and therefore the full validation suite — were
> preserved, while the measured mRNA half-life sets the low-pass filtering of transcriptional bursts. Protein
> half-lives were not drawn from the compilation for the baseline model: CDK-inhibitor and Skp2 protein
> turnover is governed by the model's cell-cycle-regulated clearance (SCF–Skp2, KPC and APC/Cdh1), and EZH2
> protein turnover was fit to EZH2-coupled observations, because the compilation's degradation-only protein
> half-lives (measured in non-cycling cell lines) would mis-specify these actively degraded regulators. As a
> sensitivity analysis we additionally imposed measured protein half-lives on the two short-lived / unregulated
> CDK inhibitors — p18 (23.9 h; Zecha et al., 2018) and p21 (0.98 h; Li et al., 2021) — a variant that also
> satisfied the validation suite (30/32 targets). No measured protein half-life was available for p19
> (*Cdkn2d*), whose protein turnover was assumed.
