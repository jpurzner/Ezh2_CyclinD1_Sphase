# Ezh2–Cyclin D1 model: experimental target compendium

All quantitative experimental values from the manuscript, assembled as calibration/validation targets for the model rebuild. Folds are relative to the stated reference. "n.r." = not reported in text.

---

## Read-before-fitting notes

1. **Table F2 "SE" column is standard deviation, not SE.** SD = SE × √N. Confirmed: DMSO G0 5.514/√12 = 1.59 = Table F1 SE. Use Table F1 for true SE; N = 12 for all groups except Palbociclib (N = 11).
2. **Phase proportions do not sum to 100%** (DMSO 87.4, HU 91.7, Palbo 81.3, RO3306 90.5). An unclassified/dead fraction of ~8–19% is excluded per group. Decide explicitly whether to renormalize to 100% before fitting.
3. **Ezh2-protein arrest values conflict between two text versions** (Section E). Reconcile before use.
4. **Ocasio et al. appears twice** in Section D with different folds (GNP panel Fig.4A vs MB panel Fig.4B). Likely different cell subsets; verify which population each refers to.
5. **Two distinct cell-cycle proportion datasets**: Section A (Fig.3B) is GNP + Tazemetostat; Section F is MB + arrest drugs. Do not pool.
6. Model-derived numbers (Section H) are outputs of the *current* model, not experimental targets. Listed for reference only.

---

## A. RNA-seq steady-state expression ratios (bulk)

| Target | Experimental value | Source |
|---|---|---|
| Cyclin D1: GNP+HHi / GNP | 0.14 (86% reduction) | RNA-seq |
| Cyclin D1: MB+HHi / MB | 0.40 (60% reduction) | RNA-seq |
| Cyclin D1: MB / GNP | 5.07x (p<0.05) | RNA-seq, Fig.4I |
| Mycn: GNP+HHi / GNP | 0.78 (22% reduction) | RNA-seq |
| Mycn: MB+HHi / MB | 0.86 (14% reduction) | RNA-seq |
| Mycn: MB / GNP | ~2.8x | RNA-seq |
| Gli1: GNP+HHi | >99% reduction | RNA-seq |
| Ezh2: MB / GNP | 2.05x (p<0.001) | RNA-seq, Fig.4J |
| Ezh2: GNP+HHi | 25–47% reduction | RNA-seq |

GNP series = WT P7 GNP, Ptch+/− P7 GNP, Ptch+/− P28 GNP, Ptch+/− MB.

---

## B. Ezh2 perturbation of Cyclin D1 transcript

| Target | Experimental value | Source |
|---|---|---|
| Cyclin D1: Ezh2 cKO GNP / WT | log2FC 1.301 (~2.46x); adj −log10 p = 16.0 | RNA-seq, Fig.3A |
| Cyclin D1: GNP + Tazemetostat 1µM / DMSO | 2.2x (p<0.05) | qPCR, Fig.3C |
| Cyclin D1: MB + Ezh2 OE (TMP 10µM) / control | significant decrease (p<0.05; fold n.r.) | qPCR, Fig.3E |
| Cyclin D1: MB + Palbociclib 1µM / DMSO | +57% (1.57x, p<0.01) | qPCR, Fig.4F |

---

## C. Palbociclib (CDK4/6i) transcript effects in MB (Fig.4F)

| Transcript | Fold vs DMSO | p |
|---|---|---|
| Ezh2 | 0.444x (−55.6%) | <0.001 |
| Cyclin D1 | 1.57x (+57%) | <0.01 |
| E2f1 | 0.433x (−56.7%) | <0.001 |

---

## D. Ezh2 transcript by cell-cycle phase (scRNA-seq, fold vs G0)

GNP datasets (Fig.4A):

| Dataset | G1/S vs G0 | p | G2/M vs G0 | p |
|---|---|---|---|---|
| Sepp et al. | 1.87x | 2.7e-8 | 2.17x | 1.03e-8 |
| Vladoiu et al. | 2.07x | 4.5e-4 | 2.11x | 1.56e-5 |
| Ocasio et al. | 2.09x | 5.79e-4 | 2.17x | 2.07e-4 |

MB datasets (Fig.4B):

| Dataset | G1/S vs G0 | p | G2/M vs G0 | p |
|---|---|---|---|---|
| Luo et al. (SHH MB) | 1.79x | 4.95e-3 | 1.89x | 3.16e-3 |
| Ocasio et al. | 2.45x | 2.48e-3 | 2.46x | 7.2e-3 |

Summary target: cycling/G0 ≈ 2.0, i.e. **G0/cycling ≈ 0.5**.

---

## E. Ezh2 protein by cell-cycle phase and arrests (IF; Fig.4G, Supp.Fig.6D)

| Target | Experimental value | Source | Note |
|---|---|---|---|
| DMSO: Ezh2 G2 / G0 (peak in G2) | 1.48x (p=2.6e-5) | Supp.6D | |
| HU 10µM: Ezh2 S-phase / DMSO | 1.31x (p=0.005) | Fig.4G | **alt text: 1.22x, p<0.001** |
| RO3306 5µM: Ezh2 G2 / DMSO | 1.16x (p=0.137, ns) | Fig.4G | **alt text: 1.20x in S, p<0.01** |
| Palbociclib: cycle-dependent Ezh2 elevation | abolished (no phase > G0) | Supp.6D | |
| E2f1 OE (TMP 10µM): Ezh2 protein / control | significant increase (p<2.2e-16; fold n.r.) | Fig.4D | |
| rShh dose escalation: Ezh2 | dose-dependent increase (qualitative) | Fig.4H | |

---

## F. Cell-cycle phase proportions — MB arrest experiments

### F1. Mean proportion (%) ± SE — true SE (24 h treatment, MB cells)

| Phase | DMSO | HU 10µM | Palbo 1µM | RO3306 5µM |
|---|---|---|---|---|
| G0 | 21.7 ± 1.59 | 25.7 ± 1.09 | 42.4 ± 4.45 | 23.7 ± 0.973 |
| G1 | 37.9 ± 2.79 | 44.1 ± 1.78 | 25.3 ± 2.32 | 27.4 ± 0.53 |
| S  | 13.7 ± 1.14 | 18.6 ± 1.87 | 4.04 ± 0.378 | 21.7 ± 1.50 |
| G2 | 14.1 ± 1.42 | 3.28 ± 0.513 | 9.56 ± 2.22 | 17.7 ± 0.866 |

> **PROVENANCE & how to compare these to the model (added after a count-vs-duration error was caught).**
> These are **flow-cytometry COUNT-fractions** — the fraction of cells *sitting* in each phase in an
> asynchronous, exponentially growing population. The DMSO targets used in `validate_v44.py`
> (24.8/43.4/15.7/16.1) are the above DMSO column renormalized to 100% (÷0.874).
> 1. **Compare to the model's age-weighted COUNT-fraction, NOT its duration-fraction** (phase time /
>    period). In exp. growth there are ~2× more just-divided than about-to-divide cells, so duration%
>    over-counts late phases. Apply the age density n(a)=2λe^(−λa), λ=ln2/Tc (the λ/Tpot correction);
>    `classify()` now returns `count_frac`. Comparing duration% to these flow numbers spuriously
>    "passed" S and G2 — a real bug.
> 2. **G0 vs G1 cannot come from DNA content** (both 2N); the G0 column requires a separate marker
>    (p27⁻/Ki67⁻), i.e. it is a **quiescent-fraction**, not a phase duration. For MB (growth fraction
>    < 1) the 2N pool also contains truly-quiescent cells a single-cycle classifier misses — needs a
>    growth-fraction term (not modeled). So fit the **resolvable 2N pool (G0+G1 = 68.2%)**; treat the
>    G0/G1 split as soft/model-internal.
> 3. **S duration:** the reliable anchor is cumulative-BrdU (control GCP Tc−Ts ≈ 20 h → Ts ≈ 3 h at
>    Tc ≈ 23 h), not the 15.7% gate; S tracks the assumed period (true Tc ≈ 28 h → Ts ≈ 8 h).
> 4. **G2/M duration:** the 16.1% 4N gate is an **unreliable duration proxy** (late-S near-4N content,
>    doublets, tetraploidy). Direct methods (Fujita G2 ≈ 2 h + M ≈ 0.5 h; pHH3/BrdU mitotic index
>    peaking < 2 h) → **G2+M ≈ 2.5 h**; anchor to that direct duration, not the 4N count.
> 5. **M ≈ 0** in the model (a brief MPF spike); report 0.5–1 h in tables and note it displaces ~2%.

### F2. Fold change vs DMSO (SD shown; N=12, Palbo N=11)

| Phase | Treatment | Mean % | SD | Fold | p | Sig |
|---|---|---|---|---|---|---|
| G0 | HU 10µM | 25.71 | 3.77 | 1.187 | 0.0485 | * |
| G0 | Palbo 1µM | 42.39 | 14.76 | 1.957 | 7.99e-4 | *** |
| G0 | RO3306 5µM | 23.69 | 3.37 | 1.094 | 0.290 | ns |
| G1 | HU 10µM | 44.06 | 6.17 | 1.163 | 0.0774 | ns |
| G1 | Palbo 1µM | 25.30 | 7.68 | 0.668 | 2.34e-3 | ** |
| G1 | RO3306 5µM | 27.43 | 1.84 | 0.724 | 3.24e-3 | ** |
| S  | HU 10µM | 18.62 | 6.49 | 1.359 | 0.0376 | * |
| S  | Palbo 1µM | 4.04 | 1.26 | 0.295 | 1.78e-6 | *** |
| S  | RO3306 5µM | 21.70 | 5.19 | 1.584 | 3.74e-4 | *** |
| G2 | HU 10µM | 3.28 | 1.78 | 0.233 | 5.00e-6 | *** |
| G2 | Palbo 1µM | 9.56 | 7.36 | 0.678 | 0.103 | ns |
| G2 | RO3306 5µM | 17.70 | 3.00 | 1.256 | 0.0433 | * |

DMSO is the reference (fold = 1.000) for every phase.

---

## A2. GNP cell-cycle proportions — Tazemetostat (Fig.3B; fold vs control)

| Phase | Fold | p |
|---|---|---|
| G0 | 0.919x | <0.0001 |
| G1 | 2.034x | <0.0001 |
| G2 | 1.439x | <0.001 |

## A3. MB cell-cycle proportions — Ezh2 overexpression (Fig.3F; fold vs control)

| Phase | Fold | p |
|---|---|---|
| G0 | 2.84x | <0.0001 |
| G1 | 0.51x | <0.0001 |
| G2 | 0.95x | ns |

---

## G. Dynamics / period

| Target | Value | Source |
|---|---|---|
| Cell cycle period, GNP | ~22–23 h | Published (ref 70, Nakashima) — **soft target** |

**The cell cycle period is a literature estimate, not measured in our system, so treat it as approximate.** It sets the overall timescale (the current model uses it only to fix the epsilon scaling factor) but is not a tight constraint. If the rebuilt model lands somewhat off ~22–23 h, that is acceptable and should not drive parameter choices at the expense of the measured ratios and proportions above.

Model-derived (current model, reference only — not fitting targets): GNP Ezh2-feedback lengthens period by 4.3 h; MB period 18.1 vs 17.3 h (±feedback); removing feedback doubles mean Cyclin D1 mRNA (2.4 → 4.7 a.u.).
