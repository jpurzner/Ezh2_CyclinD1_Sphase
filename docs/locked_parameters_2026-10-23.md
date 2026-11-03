# Locked parameter table for the particle-based rewrite
2026-10-23. Every row carries its source and a confidence grade. Nothing here
is from memory; anything I could not verify is in section D as an open gap.

**Grades:** A = multiple independent measurements agree · B = one good primary
measurement · C = single source, second-hand, or wide spread · X = not found.

---

## A. Locked (grade A–B)

| quantity | value | source | grade |
|---|---|---|---|
| nucleosomes / diploid human cell | **3 x 10^7** (6.4 Gb / ~200 bp repeat) | standard | A |
| PRC2 chromatin-bound fraction | **~20%** | [Youmans 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6049511/), endogenous HaloTag, U2OS | B |
| PRC2 D_free | **2.09 ± 0.12 µm²/s** (EZH2), 1.99 ± 0.13 (SUZ12) | Youmans 2018 | B |
| PRC2 residence time | **seconds** — scored at 3 s; Cbx7 (same study family) 4.5–11 s stable, 0.7–1.0 s transient | Youmans 2018; [Zhen 2016 eLife](https://elifesciences.org/articles/17667) | B |
| p300 bound fraction | **0.54** (0.55 ± 0.02, 0.52 ± 0.03, two clones) | [Ferrie 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10245728/), endogenous knock-in | B |
| p300 residence time | **~26 s** (FRAP — soft, FRAP/SMT differ 2–4x) | Ferrie 2023 | C |
| PRC2 affinity | **mid-nanomolar**, DNA-dominated not me3-dominated | [Choi 2017](https://www.nature.com/articles/nsmb.3488); Wang 2017 | A |
| PRC2 array saturation | bridged states saturate **5–100 nM**; majority bound at 500 nM | [Leicher 2020 PNAS](https://pmc.ncbi.nlm.nih.gov/articles/PMC7720148/) | B |
| allosteric stimulation | **~3-fold** on nucleosomal substrate | reviews; up to 3x | B |
| catalytic ordering | me1 fastest, **me2→me3 slowest** | multiple | A |
| H3K27me3 lifetime | **multiple cell divisions** (vs K4me3 6.8 h, K36me3 57 h) | [Zheng/Garcia 2014](https://ncbi.nlm.nih.gov/pubmed/24826939) | B |
| H3K27me3 domain formation | **~36 h** to steady state after EED induction | Oksuz et al. | B |
| cohesin (calibration anchor) | **250,000** nuclear complexes G1 = 1 per 120 nucleosomes | [Holzmann 2019](https://elifesciences.org/articles/46269v1) | A |

## B. Abundance — PaxDb v6.1, histone-normalised, cohesin-validated

Method: normalise each protein's ppm against total H3 in the SAME dataset
(nucleosomes = H3/2), which cancels the total-protein assumption. Validated
against cohesin, the one chromatin complex measured absolutely.

| dataset | cohesin per nucleosome | verdict |
|---|---|---|
| measured (Holzmann) | 1 per 120 | reference |
| human WHOLE_ORGANISM | 1 per 217 | usable (~1.8x) |
| human HeLa Geiger | 1 per 168 | usable (~1.4x) |
| mouse WHOLE_ORGANISM | 1 per 249 | usable (~2x) |
| mouse BRAIN | 1 per 272 | usable (~2x) |
| **mouse CEREBELLUM** | **1 per 19** | **FAILS QC — 6x wrong, histones under-detected. Do not use.** |

That the cerebellum set fails is worth noting: it is the most relevant tissue
and it is the one that cannot be trusted. Numbers below use human WO/HeLa and
mouse WO/BRAIN only.

| protein | human ppm | mouse ppm | nucleosomes per molecule | copies/cell (cohesin-anchored) | grade |
|---|---|---|---|---|---|
| **EZH2** | 1.81 / 2.39 | 3.49 | **830 – 4,200** | **15,000 – 36,000** | B |
| EZH1 | 0.295 | 0.616 | 5,100 – 23,800 | 2,600 – 5,900 | C |
| SUZ12 | 3.33 / 10.9 | 26.0 | 450 – 560 | 67,000 – 110,000 | B |
| EED | 3.48 / 9.95 | 8.19 | 430 – 1,790 | 35,000 – 70,000 | B |
| **EP300** | 0.735 / 0.892 | 9.52 | **1,540 – 5,500** | **5,500 – 40,000** | C (6x spread) |
| CREBBP | 0.616 / 2.17 | 13.9 | 1,055 – 2,435 | 12,000 – 59,000 | C |
| KDM6A | 0.409 | 0.458 | 3,700 – 32,000 | 1,900 – 8,200 | C |
| **KDM6B** | 0.241 | 0.056 / 0.001 / 4.96 | **6,200 – 11,800,000** | **240 – 4,800** | **X — see D** |
| HDAC1 | 34.8 / 289 | 106 | 17 – 138 | 450,000 – 1,800,000 | B |
| HDAC2 | 33.5 / 458 | 61.1 | 11 – 240 | 259,000 – 2,800,000 | B |

**EZH2 is the limiting PRC2 subunit** — SUZ12 and EED run 2–7x higher, so most
is not in a complete complex. PRC2 holocomplex count ≈ EZH2 + EZH1.

## C. Derived numbers for the model

Per **50-nucleosome promoter domain**, using measured bound fractions and
assuming concentration at ~3,000 PcG target elements:

| | copies | bound | per target domain |
|---|---|---|---|
| PRC2 | 15,000–36,000 | 20% → 3,000–7,200 | **1.0 – 2.4** |
| p300 | 5,500–40,000 | 54% → 3,000–21,600 | 1.0 – 7.2 |
| KDM6B | 240–4,800 | ~20%? → 50–960 | 0.02 – 0.32 |
| HDAC1+2 | ~10^6 | — | ~1 per 11–45 **nucleosomes** |

Spread uniformly over the genome instead, all writers/erasers fall to ~0.01
per domain.

**Visit statistics.** With 1–2 bound PRC2 at ~5 s residence, a domain
experiences one unbinding every ~2.5–5 s, i.e. **~700–1,400 binding events per
hour**. Domain formation takes ~36 h and deposits of order 30 me3 marks, so
roughly **1 in 10^3 visits is catalytically productive**. Binding and catalysis
must be modelled as separate events.

## D. NOT locked — open gaps

| gap | status |
|---|---|
| **PRC2 k_cat on nucleosomes** | **X — not found** after four targeted searches. Patents give the method, the 2026 S21A preprint gives kcat/Km on *peptides*. This is the single most important missing number: it decides whether a ~5 s visit can deposit anything. |
| **KDM6B abundance** | **X — unlockable from proteomics.** Spans 0.001–4.96 ppm across datasets, a ~5,000x range. KDM6B is an inducible immediate-early gene, so a tissue average is probably meaningless; it needs direct measurement in GNP/MB. **The lattice currently makes KDM6B a dominant antagonist on no abundance and no kinetic data at all.** |
| **KDM6B / UTX k_cat** | **X.** Only KDM6A ~**0.29 min⁻¹**, and that is second-hand via a QM/MM paper citing an inaccessible kinetic study; substrate and conditions unverified. Grade C at best. |
| PRC2 residence, exact tau | reported ~7.3 s in K27M and ~1.6x lower in WT (≈4.6 s) in the DIPG study, but the values live in figures I could not extract. Treat "a few seconds" as the usable statement. |
| GNP / MB specific abundances | the one relevant tissue dataset fails QC (above). Needs its own proteomics. |
| HDAC residence / k_cat | not searched. |
| p300 k_cat | not searched. |

## E. What this changes in the model, concretely

1. `prc2_activity()` computes a continuous field over 51 nucleosomes. Replace
   with **1–2 discrete PRC2 particles per domain**.
2. **Separate binding from catalysis** — ~10^3 visits per productive methylation.
3. **Read-write should modulate catalysis once bound, not binding affinity.**
   Youmans: blocking the EED cage for 30 min does not change residence; only
   72 h, after the mark is gone, does. Choi/Wang: DNA dominates affinity.
   `a_rw_prc2` currently does the opposite.
4. **Drop the allosteric boost from 5–10x to ~3x** and re-test the bistability.
5. **HDACs are the one enzyme that IS a field** (1 per 11–45 nucleosomes) — the
   existing uniform-constant treatment is correct.
6. **The KDM6B arm is currently unfalsifiable.** Either measure it in GNP/MB or
   demote it to a sensitivity parameter and stop drawing conclusions from it.
