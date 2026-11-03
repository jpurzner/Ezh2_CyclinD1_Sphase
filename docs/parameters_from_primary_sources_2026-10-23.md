# Parameters read from the primary papers
2026-10-23. Supersedes `locked_parameters_2026-10-23.md`, which was built from
abstracts, reviews and search snippets. Everything here comes from Methods and
Results text I have actually read. **Corrections to my earlier numbers are
flagged in bold.**

Full texts pulled: Holzmann 2019 (PMC6606026), Youmans 2018 (PMC6049511),
Tatavosian 2018 (PMC5970213), Wisniewski 2014 (PMC4256500), Ferrie 2023
(PMC10245728).

---

## 1. CORRECTION — the cohesin anchor was the wrong number

My whole copy-number calibration used "250,000 cohesin". Reading
[Holzmann 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6606026/), that paper
reports **four different numbers** and I picked the wrong one:

| quantity | value | method |
|---|---|---|
| **total cellular SCC1/RAD21** | **~350,000 copies / G1 cell** | SRM with isotope-labelled reference peptides (EtEP), 5500 QTRAP, 2 biological x 2 technical |
| total cellular SMC3 | ~417,000 | same (excess over SCC1 = SMC1-SMC3 dimers not bound to SCC1) |
| chromatin-bound, **biochemical** | **60,000–70,000** complexes (SMC1/SMC3/SCC1 61–69k; STAG1 15k + STAG2 45k = 60k) | MS on chromatin fraction, differential centrifugation |
| nuclear SCC1-GFP | ~250,000 | **FCS**, in vivo |
| chromatin-bound, **in vivo** | **160,000 ± 10,000** (64 ± 4% of nuclear) | **iFRAP** |

Two things matter here.

**(a) PaxDb measures total cellular protein, so the correct anchor is 350,000,
not 250,000.** My numbers were ~1.4x low on that account alone.

**(b) The same paper reports chromatin-bound cohesin as 60–70k by biochemical
fractionation and 160k by in-vivo iFRAP — a 2.5x discrepancy, unresolved.**
No abstract or review says this. It is a standing caution for every
"chromatin-bound fraction" in this field, including the PRC2 and p300 numbers
below.

Also worth having: Holzmann concludes from ChIP-seq comparison that **"some
genomic cohesin and CTCF enrichment sites are unoccupied in single cells at any
one time."** That is the same conclusion I am drawing for PRC2, already
established for a better-measured complex.

## 2. CORRECTION — the two "bound fractions" are not comparable

I wrote that "p300 is ~2.7x more chromatin-bound than PRC2 (0.54 vs 0.20)".
**Retract that.** The two numbers come from different algorithms with different
definitions:

| | PRC2 (Youmans 2018) | p300 (Ferrie 2023) |
|---|---|---|
| method | **Spot-On, two-state kinetic model** | **State Array, Bayesian diffusion spectrum** |
| definition of bound | fitted 2nd state, `D_bound` constrained to **[0.0005, 0.08] µm²/s** | fractional likelihood of **D < 0.1 µm²/s** |
| value | ~20% | 0.55 ± 0.02 and 0.52 ± 0.03 (2 clones) |
| acquisition | 97.5 Hz, 10.3 ms lag, dZ 0.7 µm, LocError 0.035 µm | 223 Hz (4.48 ms), focal depth 0.7 µm |
| cells | U2OS, endogenous HaloTag knock-in | U2OS, endogenous knock-in |
| n | 3 biological reps, >12 cells each | bootstrap mean ± SD, 2 clones |

A two-state fit forces every molecule into one of two states. Integrating a
continuous diffusion spectrum below a threshold additionally captures slow-but-
unbound molecules. These will not agree even on identical data. The cell line
and tag strategy do match, so a like-for-like reanalysis would be worth doing,
but as published the ratio is not interpretable.

## 3. CORRECTION — there is no fitted PRC2 residence time

Youmans measure residence by **empirical survival probability**, not a fitted
tau: intermittent acquisition at **2 Hz, 500 ms lag**, tracking only particles
with expected D_max 0.1 µm²/s, then an empirical CDF. The reported result is
*"a significant decrease in EZH2 that remains chromatin-bound after 3 sec"*
under 72 h A-395.

So: no tau exists in that paper, and 500 ms lag with no stated photobleaching
correction will truncate long events. **The defensible statement is "PRC2
chromatin binding is resolved on a seconds timescale", not any specific value.**
Tatavosian 2018 reports Ezh2 residence and search times for H3.3K27M vs WT, but
those values live in figures I could not extract from the text.

Their other finding is clean and verified from the abstract text: **H3.3K27M
prolongs Ezh2 residence and search time but does NOT change its bound
fraction**, whereas for Cbx7 it does the reverse.

## 4. VERIFIED — the proteomic ruler, and a correction to how I used it

[Wisniewski 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4256500/) is the
method I improvised. Their actual approach and numbers:

- DNA mass of a diploid human cell = **6.5 pg**; histone mass ≈ DNA mass.
- In their four-cell-line dataset the **histone MS signal is 2.07–4.03% of
  total MS signal**.
- Equating that with 6.5 pg recovers cellular protein mass within
  **1.24 ± 0.29x** of cell counting.
- They use **total histone**, not H3 alone.

**Correction to my method:** I normalised to H3 only (nucleosomes = H3/2). On a
total-histone basis (nucleosomes = sum of core histones / 8) the answer differs
by **2.1x** in the human whole-organism set. Total histone is the right basis
because it averages out variant-specific MS bias — and that bias is large: in
molar ppm, H3:H4:H2A:H2B should be 1:1:1:1 and comes out **5,435 : 11,624 :
11,194 : 16,639**.

**Independent confirmation of the nucleosome count:** 6.5 pg DNA / 13.8 kDa mean
histone mass = 2.84 x 10^8 histone molecules, / 8 = **3.55 x 10^7 nucleosomes**.
This agrees with the 6.4 Gb / 200 bp route (3 x 10^7) by a different argument.

**Dataset QC from this:** human WHOLE_ORGANISM gives total core histone = 4.5%
of molecules (close to Wisniewski's 2–4%). HeLa Geiger gives **8.0%**, out of
range, and its H2A:H2B ratio is 1:11 when it should be 1:1 — so **HeLa Geiger's
histone quantification is unreliable** and the whole-organism set should be
preferred.

## 5. REVISED abundance table

Anchored on the verified RAD21 = 350,000 copies/G1 cell, total-histone basis,
human WHOLE_ORGANISM (the dataset that passes histone QC).

PaxDb under-reports RAD21 relative to histone by **5.2x** (whole-organism) or
**2.8x** (HeLa Geiger) against the SRM truth. That factor is applied below, and
it is the dominant systematic uncertainty.

| protein | ppm | nucleosomes / molecule (corrected) | copies/cell | previously reported | 
|---|---|---|---|---|
| **EZH2** | 1.81 | **~594** | **~50,700** | was 15,000–36,000 → **1.4–3x low** |
| EZH1 | 0.295 | 3,644 | 8,260 | — |
| SUZ12 | 3.33 | 323 | 93,240 | was 67,000 |
| EED | 3.48 | 309 | 97,440 | was 61,000–70,000 |
| **EP300** | 0.735 | **1,463** | **~20,600** | was 5,500–15,000 → **low** |
| CREBBP | 0.616 | 1,745 | 17,250 | — |
| KDM6A | 0.409 | 2,628 | 11,450 | was 8,200 |
| **KDM6B** | 0.241 | 4,461 | **~6,750** | was 240–4,800 |
| HDAC1 | 34.8 | 31 | 974,000 | — |
| HDAC2 | 33.5 | 32 | 938,000 | — |

EZH2 remains the limiting PRC2 subunit (SUZ12/EED ~2x higher).

## 6. REVISED model-relevant numbers

PRC2 holocomplexes ≈ EZH2 + EZH1 ≈ **59,000**. Applying the (caveated, two-state)
20% bound fraction gives ~12,000 bound PRC2.

| targets assumed | bound PRC2 per 50-nucleosome domain |
|---|---|
| 3,000 PcG targets | **3.9** |
| 10,000 | 1.2 |
| 30,000 | 0.4 |
| uniform over genome | 0.02 |

**So: ~1–4 PRC2 complexes per target promoter domain.** The qualitative
conclusion from before survives — PRC2 is discrete and scarce, not a continuous
field — but the number is ~1.6x higher than I reported, and the bound fraction
underlying it is model-dependent.

HDACs remain ~1 per 31 nucleosomes, i.e. genuinely field-like.

## 7. Still NOT verified

| gap | status |
|---|---|
| PRC2 k_cat on nucleosomes | **X** — still not found. The productive-visit estimate (~1 in 10^3) remains an inference from domain-formation time, not a measurement. |
| PRC2 residence tau | **X** — no fitted value exists in Youmans; Tatavosian's are in figures. |
| KDM6B abundance in GNP/MB | **X** — PaxDb whole-organism gives 6,750 copies, but mouse datasets span 0.001–4.96 ppm and the one cerebellum dataset fails histone QC (cohesin 1 per 19). |
| KDM6B/UTX k_cat | **X** — only KDM6A ~0.29 min⁻¹, second-hand via a QM/MM paper citing an inaccessible study. |
| Choi 2017 Kd / residence values | not retrieved — NSMB, not open access, not on PMC. |
| Oksuz 36 h figure | not verified from primary text (PMC7700016 is not OA). |
| Zheng H3K27me3 half-life | not verified from primary text (PMC4184988 is not OA); the 6.8 h / 57 h figures for K4me3 / K36me3 are from the abstract only. |
