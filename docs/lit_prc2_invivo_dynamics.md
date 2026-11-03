# Live-cell PRC2 / PcG dynamics: a method-annotated census

Compiled 2026-10-23. **Every number below was read in the Methods, Results, or a
supplementary table of the cited paper unless the Grade column says otherwise.**
Numbers that exist only inside a plotted figure are marked `FIGURE-ONLY`.

## Grading key

| Grade | Meaning |
|---|---|
| **A** | Numeric value printed in Results/Methods text **or** a supplementary table, with the full method stack recoverable (model, bound-definition, frame rate, dZ, bleach correction, tagging). |
| **B** | Numeric value printed in text, but one or more method parameters are missing or must be inferred. |
| **C** | Value reported only as a fold-change, a secondary citation, or a qualitative statement; or value is a plotted bar with no printed number. |
| **X** | Abstract-only, or in-vitro / reconstituted rather than live-cell (kept for contrast, **never comparable to live-cell numbers**). |

---

## 1. The core problem, stated concretely

Four published "bound fractions" of **the same protein (EZH2) in live cells**:

- **0.20** (Youmans 2018, U2OS, Spot-On **2-state**, D_bound forced into [0.0005, 0.08] µm²/s)
- **0.23** (Tatavosian 2018, mESC, Spot-On **3-state**, D_bound forced into [0.0001, 0.01] µm²/s)
- **0.046** (Tatavosian 2018, mESC, *specific* long-lived fraction F1sb from a 2-exponential dwell fit)
- **0.08** (Youmans 2018, U2OS, EZH2 with SUZ12-VEFS only — and the authors say this is indistinguishable from the ~0.10 background of free HaloTag)

These are four different quantities, not four estimates of one quantity. Section 5
lays out which pairs can legitimately be compared.

Three published residence times of **the same protein (Cbx7) in the same cell type
(mESC), from the same lab, with nominally the same protocol**: **7.3 s**
(Zhen 2016), **9.4 s** (Tatavosian 2018), **6.9 s** (Kent 2020). This ±30% spread
across identical nominal methods is the practical floor on how precisely any of
these numbers should be read.

---

## 2. Master table

Abbreviations: `τsb` = long-lived/"specific" residence time; `τtb` = short-lived/
non-specific residence time; `F1` = total chromatin-bound fraction of all molecules;
`F1sb` = long-lived bound fraction of all molecules; `f1sb` = long-lived fraction
*of bound molecules*; `τsearch` = target-search time; `KI` = endogenous CRISPR knock-in.

### 2.1 PRC2 subunits — EZH2, EED, SUZ12

| Protein | Value | Quantity | METHOD / MODEL | Bound definition | Cell line | Tagging | Source | Grade |
|---|---|---|---|---|---|---|---|---|
| EZH2 | 0.20 ± 0.03 | F1 (fraction bound) | Spot-On, **2-state**, CDF fit, 3 iterations, 97.5 Hz (10.3 ms lag), 7 timepoints, dZ 0.7 µm, LocError 0.035 µm fixed, SLIMfast+evalSPT localisation | D_bound **constrained** to [0.0005, 0.08] µm²/s; D_free [0.15, 25] | U2OS (lacO-array parental, clone F42B8) | **endogenous KI**, 3xFlag-HaloTag-EZH2, JF646 | Youmans 2018, Genes Dev, PMC6049511 | **A** |
| SUZ12 | 0.25 ± 0.03 | F1 | same as above | same | U2OS | **endogenous KI**, 3xFlag-HaloTag-SUZ12 | Youmans 2018, PMC6049511 | **A** |
| EZH2 | 2.09 ± 0.12 µm²/s | D_free | Spot-On 2-state | — | U2OS | endogenous KI | Youmans 2018, PMC6049511 | **A** |
| SUZ12 | 1.99 ± 0.13 µm²/s | D_free | Spot-On 2-state | — | U2OS | endogenous KI | Youmans 2018, PMC6049511 | **A** |
| EZH2 | survival probability curve only; **no fitted τ anywhere in the paper** | residence | **empirical survival probability (MATLAB ecdf)** on tracked particles; 2 Hz, 500 ms lag; expected D_max 0.1 µm²/s; track length >3 frames; **photobleaching correction not stated** | "static" = tracked with D_max 0.1 µm²/s | U2OS | endogenous KI | Youmans 2018, PMC6049511, Fig 4C | **B** (quantity exists, no τ) |
| EZH2 | ~0.40 | immobile fraction at H3K27me3-rich lacO array | **FRAP**, confocal, bleach 2 s, 5 min recovery, fit y=A(1−e^−Bt), immobile = C − A | immobile fraction of FRAP fit | U2OS F42B8 | endogenous KI | Youmans 2018, PMC6049511 (value stated in Discussion; Fig 3C is the plot) | **C** (FIGURE-ONLY for the per-site values; 40% is from the Discussion prose) |
| Ezh2 | **(23 ± 1)%** | F1 | **Spot-On, 3-state**, CDF fit (ModelFit=2), FitIterations=5, lag Δτ = 30 ms (33.3 Hz), TimePoints=2, GapsAllowed=2, JumpsToConsider=4, dZ **0.700 µm**, **LocError 0.035 µm fixed** (FitLocError=0), U-track localisation, HILO, JF549 | D_bound **constrained** to **[0.0001, 0.01] µm²/s**; D_free2 [0.01, 1]; D_free1 [1, 5] | **mouse ES cells** (PGK12.1 WT) | tet-inducible lentiviral HaloTag-Ezh2, imaged **uninduced** (basal) | Tatavosian 2018, Nat Commun, PMC5970213, **Suppl. Table 1** | **A** |
| Ezh2 | (22 ± 1)% | F1 | as above | as above | mESC **Ezh2−/−** | HaloTag-Ezh2 rescue | Tatavosian 2018, Suppl. Table 1 | **A** |
| Eed | (20 ± 1)% / (22 ± 1)% | F1 (WT / Eed−/−) | as above | as above | mESC | HaloTag-Eed | Tatavosian 2018, Suppl. Table 1 | **A** |
| Ezh2 | F2 = (23 ± 4)% at D2 = 0.54 ± 0.09 µm²/s; F3 = (54 ± 5)% at D3 = 4.2 ± 0.4 µm²/s | confined / free populations | Spot-On 3-state | — | mESC | as above | Tatavosian 2018, Suppl. Table 1 | **A** |
| Eed | F2 = (29 ± 4)% at 0.52 ± 0.07 µm²/s; F3 = (51 ± 3)% at 4.4 ± 1.2 µm²/s | confined / free | Spot-On 3-state | — | mESC | as above | Tatavosian 2018, Suppl. Table 1 | **A** |
| Ezh2 | **τtb = 0.9 ± 0.1 s; τsb = 10.0 ± 0.2 s** | residence (2-exponential) | **separate slow-acquisition experiment**: 30 ms integration + 170 ms dark (200 ms cycle, 5 Hz), 400 frames, low laser (3.4 mW/cm² vs 10.6 for diffusion); survival curves **photobleaching-normalised** against a measured JF549/H2A-HaloTag bleach curve (2-component); fit chosen by F-test | **D_m ≤ 0.032 µm²/s** per-track maximum-likelihood D (note: *different* threshold from the 0.01 µm²/s Spot-On constraint in the same paper) | mESC WT | HaloTag-Ezh2 | Tatavosian 2018, **Suppl. Table 2** | **A** |
| Ezh2 | F1tb = (18 ± 1)%, **F1sb = (4.6 ± 0.3)%** | transient / stable bound fractions of all molecules | F1tb = F1×(1−f1sb), F1sb = F1×f1sb | as above | mESC WT | HaloTag-Ezh2 | Tatavosian 2018, Suppl. Table 2 | **A** |
| Ezh2 | **τsearch = 210 ± 19 s**; Ntrial = 5.0 ± 0.4; τ3D = 41.6 ± 5.2 s | target search | Loffreda-type 2-binding-state kinetic model: τsearch = Ntrial·τ3D + (Ntrial−1)·τtb, Ntrial = 1/f1sb | derived from τtb/τsb/F1sb above | mESC WT | HaloTag-Ezh2 | Tatavosian 2018, Suppl. Table 2 | **A** |
| Eed | τtb = 0.9 ± 0.2 s; τsb = 9.5 ± 0.1 s; F1tb = (16 ± 3)%; F1sb = (4.3 ± 0.1)%; τsearch = 214 ± 6 s | residence + search | as Ezh2 | D_m ≤ 0.032 µm²/s | mESC WT | HaloTag-Eed | Tatavosian 2018, Suppl. Table 2 | **A** |
| Eed | τsb = 9.2 ± 0.2 s; F1sb = (4.1 ± 0.1)%; τsearch = 220 ± 11 s | residence + search | as above | as above | mESC **Eed−/−** | HaloTag-Eed rescue | Tatavosian 2018, Suppl. Table 2 | **A** |
| Ezh2 (catalytic lobe only) | F1 = (17 ± 3)%; τtb 0.8 ± 0.1 s; **τsb 4.7 ± 0.1 s**; F1sb (2.8 ± 0.5)%; τsearch 162 ± 25 s | all | as above | as above | HEK293T + H3.3-FLAG | HaloTag-Ezh2^Catalytic, stable overexpression | Tatavosian 2018, Suppl. Tables 7–8 | **A** |
| EZH2 | 0.55 ± 0.02 and 0.52 ± 0.03 | f_bound | **State Array** (Bayesian diffusion spectrum), 223 Hz | fractional likelihood of **D < 0.1 µm²/s** | U2OS | — | *(p300, not PRC2 — supplied by the requester as the method contrast; Ferrie 2023, PMC10245728)* | — |

### 2.2 PRC2: Tatavosian 2018 DIPG / H3.3K27M — the values the requester asked for

All from **Supplementary Tables 5 and 6** of PMC5970213 (the figures plot these).
Method stack is identical to §2.1 Tatavosian rows: **Spot-On 3-state, D_bound
[0.0001, 0.01] µm²/s, 30 ms lag, dZ 0.7 µm, LocError 0.035 µm** for fractions;
**30 ms + 170 ms dark, D_m ≤ 0.032 µm²/s, photobleach-normalised 2-exponential**
for residence; Loffreda model for search. Tagging is **stable lentiviral
overexpression of HaloTag fusions**, not knock-in. Grade **A** throughout.

**HEK293T, H3.3-FLAG vs H3.3K27M-FLAG (30 ms exposure rows):**

| Protein / condition | F1 | F1tb | F1sb | τtb (s) | **τsb (s)** | Ntrial | τ3D (s) | **τsearch (s)** |
|---|---|---|---|---|---|---|---|---|
| HT-Ezh2 / H3.3 | (23 ± 4)% | (19 ± 4)% | (4.8 ± 1.0)% | 1.0 ± 0.1 | **10.2 ± 0.2** | 5.0 ± 0.4 | 40.0 ± 8.6 | **204 ± 40** |
| HT-Ezh2 / **H3.3K27M** | (23 ± 3)% | (18 ± 3)% | (5.2 ± 0.7)% | 1.0 ± 0.1 | **15.3 ± 0.2** | 4.4 ± 0.4 | 62.6 ± 8.5 | **278 ± 30** |
| HT-Cbx7 / H3.3 | (46 ± 1)% | (35 ± 1)% | (11.8 ± 0.3)% | 1.0 ± 0.1 | 9.7 ± 0.1 | 4.0 ± 0.3 | 17.7 ± 1.5 | **73 ± 3** |
| HT-Cbx7 / **H3.3K27M** | (32 ± 2)% | (25 ± 1)% | (7.0 ± 0.4)% | 1.1 ± 0.1 | 9.9 ± 0.1 | 4.6 ± 0.4 | 27.8 ± 3.0 | **133 ± 10** |

**Patient-derived DIPG line SF8628 (heterozygous H3.3K27M) vs control brain-tumour line 9427:**

| Protein / line | F1 (30 ms) | F1tb | F1sb | τtb (s) | **τsb (s)** | Ntrial | τ3D (s) | **τsearch (s)** |
|---|---|---|---|---|---|---|---|---|
| HT-Ezh2 / **9427** | (21 ± 2)% | (17 ± 2)% | (4.4 ± 0.4)% | 0.9 ± 0.1 | **8.7 ± 0.1** | 4.9 ± 0.3 | 37.8 ± 3.8 | **188 ± 12** |
| HT-Ezh2 / **SF8628** | (22 ± 2)% | (17 ± 1)% | (4.3 ± 0.3)% | 1.0 ± 0.1 | **11.8 ± 0.2** | 3.4 ± 0.4 | 76.2 ± 10.9 | **265 ± 24** |
| HT-Cbx7 / **9427** | (42 ± 1)% | (32 ± 1)% | (10.1 ± 0.3)% | 1.0 ± 0.1 | 8.5 ± 0.1 | 4.2 ± 0.3 | 17.5 ± 1.6 | **76 ± 4** |
| HT-Cbx7 / **SF8628** | (30 ± 3)% | (23 ± 2)% | (7.1 ± 0.6)% | 1.1 ± 0.1 | 8.5 ± 0.1 | 3.9 ± 0.3 | 28.1 ± 4.1 | **119 ± 15** |

**Ezh2 catalytic lobe, HEK293T (Suppl. Tables 7–8):** H3.3 τsb 4.7 ± 0.1 s,
τsearch 162 ± 25 s; **H3.3K27M τsb 7.3 ± 0.2 s, τsearch 277 ± 24 s**; F1 17–18%
in both; F1sb 2.8% vs 2.6%.

**Exposure-time series (motion-blur control), F1 for HT-Ezh2 in HEK293T/H3.3:**
5 ms 18%, 10 ms 18%, 20 ms 21%, 30 ms 23%. For HT-Cbx7/H3.3: 5 ms 39%, 10 ms 40%,
20 ms 44%, 30 ms 46%. **The bound fraction moves by 5–7 percentage points purely
from exposure time at fixed lag.** The authors state explicitly that motion blur
still inflates the bound fraction even below 5 ms exposure.

**In-vitro contrast (same paper, Fig 6h–i, grade X for live-cell comparison):**
FLAG-HaloTag-Ezh2-PRC2 purified from mESC, on tethered mononucleosomes, TIRF,
30 ms integration + 470 ms dark, **1-component exponential, photobleaching NOT
corrected (authors say the value is a lower limit)**: **τ = ~257 s on
H3.3-nucleosome, ~404 s on H3.3K27M-nucleosome.** This is 25–40× the in-vivo τsb
measured in the same paper.

### 2.3 PRC1 — CBX proteins

| Protein | Value | Quantity | METHOD / MODEL | Bound definition | Cell line | Tagging | Source | Grade |
|---|---|---|---|---|---|---|---|---|
| Cbx7 | **F1 = (29 ± 1)%** (F2 49%, F3 22%) | bound fraction | **NOT Spot-On.** 3-component **Gaussian mixture fit to the per-track log D_m histogram**, 30 ms interval, U-track; **no dZ / defocalisation correction, no localisation-error term** | **Bound component centre pinned at D_m1 = 0.032 µm²/s** (taken from H2A-HaloTag) | mESC WT | lentiviral HaloTag-Cbx7 | Zhen 2016, eLife, PMC5056789 | **A** |
| Cbx7 | F1 = (30 ± 1)% | bound fraction | as above | as above | mESC **Cbx7−/−** | HaloTag-Cbx7 rescue | Zhen 2016 | **A** |
| Cbx2 / Cbx4 / Cbx6 / Cbx8 | F1 = 30% / 10% / 17% / 26% (±1) | bound fraction | as above | as above | mESC WT | lentiviral HaloTag | Zhen 2016 | **A** |
| Cbx7 | **τtb = 0.79 ± 0.01 s (F1tb 23.4 ± 1.1%); τsb = 7.3 ± 0.1 s (F1sb 5.3 ± 0.4%)** | residence | 30 ms + 170 ms dark, 2-exponential on cumulative dwell distribution, **photobleaching corrected** against a measured 9-curve H2A-HaloTag/JF549 bleach decay (2-component) | **D_m < 0.10 µm²/s** (97% of H2A-HaloTag falls below) — *a 3× looser threshold than the 0.032 used by the same lab in 2018/2020* | mESC WT | HaloTag-Cbx7 | Zhen 2016 | **A** |
| Cbx7 variants | CD^Cbx7 τsb 4.7 ± 0.1 s; Cbx7^F11A 5.8 ± 0.1 s; Cbx7^ΔCD 4.7 ± 0.1 s; Cbx7^ATLm 6.0 ± 0.2 s; Cbx7^ΔATL 5.8 ± 0.1 s; Cbx7^ΔCD-ATL 4.6 ± 0.1 s | residence | as above | as above | mESC WT | HaloTag fusions | Zhen 2016 | **A** |
| Cbx7 | bound fraction collapses so far in **Eed−/−** that the 3-component fit fails to converge (only F2 65 ± 3%, F3 35 ± 6% reported) | H3K27me3 dependence | as above | as above | mESC Eed−/− | HaloTag-Cbx7 | Zhen 2016 | **A** |
| Cbx8 | F1 falls 26% → **(8 ± 2)%** in Eed−/− | H3K27me3 dependence | as above | as above | mESC Eed−/− | HaloTag-Cbx8 | Zhen 2016 | **A** |
| Cbx2 / Cbx4 / Cbx6 | F1 30→26%, 10→14%, 17→13% in Eed−/− | H3K27me3 dependence | as above | as above | mESC Eed−/− | HaloTag | Zhen 2016 | **A** |
| Cbx7 | **F1 = (40 ± 2)%** (F2 40%, F3 20%, D2 0.60 ± 0.10, D3 4.1 ± 0.1 µm²/s) | bound fraction | **Spot-On 3-state** (see §2.1 Tatavosian stack) | D_bound [0.0001, 0.01] µm²/s | mESC WT | lentiviral HaloTag-Cbx7 | Tatavosian 2018, Suppl. Table 1 | **A** |
| Cbx7 | τtb 0.8 ± 0.1 s; **τsb 9.4 ± 0.1 s**; F1tb (27 ± 1)%; **F1sb (12.9 ± 0.6)%**; Ntrial 3.1 ± 0.2; τ3D 19.9 ± 2.1 s; **τsearch 64 ± 5 s** | residence + search | 30 ms + 170 ms dark; D_m ≤ 0.032 µm²/s; bleach-corrected 2-exponential | — | mESC WT | HaloTag-Cbx7 | Tatavosian 2018, Suppl. Table 2 | **A** |
| CBX2 | **τsb = 8.9 s** (Cbx2−/− rescue: 9.1 s) | residence | 30 ms + 170 ms dark, 400 frames, 2.25 mW; 2-exponential, **photobleach-corrected**; U-track | **D_m ≤ 0.032 µm²/s** | mESC (PGK12.1) | tet-inducible lentiviral HT-CBX2 | Kent 2020, Cell Rep, PMC7593837 | **A** |
| RING1B / MEL18 / PHC1 | **τsb = 5.5 / 5.7 / 5.6 s** | residence | as above | D_m ≤ 0.032 µm²/s | mESC | HT fusions | Kent 2020, PMC7593837 | **A** |
| CBX6 / CBX7 | **τsb = 2.7 s / 6.9 s** | residence | as above | D_m ≤ 0.032 µm²/s | mESC | HT fusions | Kent 2020, PMC7593837 | **A** |
| CBX2 | **F1sb ≈ 11%**; **τsearch ≈ 72 s** (Ntrial ≈ 4, Ntrial·τ3D ≈ 69 s, Ntrial·τtb ≈ 3 s) | stable bound fraction + search | Spot-On 3-state for F1, then the Loffreda 2-binding-state model | D_m ≤ 0.032 µm²/s | mESC WT | HT-CBX2 | Kent 2020, PMC7593837 | **A** |
| CBX2 | τsb = **9.7 s** in Eed−/− (vs 8.9 s WT, i.e. unchanged) but **F1sb falls to 7.9%** (≈70% of WT) and Ntrial, τ3D, τsearch all rise | H3K27me3 dependence | as above | as above | mESC Eed−/− | HT-CBX2 | Kent 2020, PMC7593837 | **A** |
| CBX2 | τsb = 7.8 s in Ring1a−/−/Ring1b−/− and in Mel18−/−/Bmi1−/−; search parameters unchanged | PRC1 dependence | as above | as above | mESC | HT-CBX2 | Kent 2020, PMC7593837 | **A** |
| CBX2 | CBX2^CD alone: **τsb = 0 s** (no detectable tight binding); CBX2^65–532 τsb 8.9 s; CBX2^1–498 τsb 9.4 s | domain mapping | as above | as above | mESC | HT fusions | Kent 2020, PMC7593837 | **A** |
| CBX2 | NPC bound fraction > mESC; **stable bound fraction ~2× larger in NPCs**; **τl similar in both**; **τsearch ~2× shorter in NPCs**; ~1400 bound CBX2 molecules in mESC vs ~4900 in NPC | cell-type dependence | Spot-On 3-state at 30 ms/0 ms dark for F1; 30 ms + 170 ms dark for τ; D_m < 0.032 µm²/s; photobleach-corrected | D_m < 0.032 µm²/s | **mESC and mouse NPCs** | **endogenous CBX2-HaloTag knock-in** (Cbx2^fl/fl-HT) | Yang/Ren 2026, Mol Cell, PMC12981346, Fig 2L–O | **C** (FIGURE-ONLY: fold-changes in text, absolute values only in the plots) |
| CBX2 | ~3 molecules per condensate (~8 nM) | condensate stoichiometry | single-molecule photobleaching step counting, Chung-Kennedy filter, 30 ms exposure | — | mESC | endogenous KI | Yang/Ren 2026, PMC12981346 | **B** |
| H2A-HaloTag (bound control) | 79% (2026) / **82 ± 2%** (2018) / **72 ± 1%** (2016) | F1 | same lab, three pipelines | varies with the pipeline above | mESC | HaloTag | Yang 2026 / Tatavosian 2018 / Zhen 2016 | **A** |
| HaloTag-NLS (free control) | ~10% (2026) / **8 ± 1%** (2018) / F1 not resolvable (2016) | F1 | same | same | mESC | HaloTag-NLS | as above | **A** |

### 2.4 PRC1 — Klose lab, endogenous knock-ins, mESC

Method stack for every Huseyin & Klose 2021 row (Nat Commun, PMC7873255):
**Spot-On, 3-state** (a 2-state model "fit the data poorly"), **67 Hz / 15 ms
exposure**, PA-JF549 photoactivatable dye, HILO, bin width 0.01 µm, **8 timepoints,
4 jumps, max jump 5.05 µm, localisation error 40 nm, dZ = 0.7 µm, CDF fitting,
3 iterations**; 25 nm localisation precision, 768 nm linking radius; tracks fitted
**per movie** and statistics taken over movies. **Bound = the immobile state of
the 3-state fit; no explicit D cut-off is published.** Tagging = **homozygous
endogenous CRISPR HaloTag knock-in** (RING1B, PCGF1/2/3/6, CBX7, RYBP), except the
RING1B^NBM and PCGF1-fragment mutants, which are stable transgenes.

| Protein | Value | Quantity | Notes | Grade |
|---|---|---|---|---|
| RING1B | **20%** | bound fraction | n > 50 movies, 2 experiments | **A** |
| H2B (control) | **59%** | bound fraction | same pipeline | **A** |
| HaloTag-3xNLS (control) | **11%** | bound fraction | same pipeline | **A** |
| PCGF2 | **18%** | bound fraction | n > 50 movies | **A** |
| CBX7 | **18%** | bound fraction | n > 50 movies | **A** |
| PCGF1 | **26%** | bound fraction | highest of all PCGFs | **A** |
| PCGF3 | **13%** | bound fraction | | **A** |
| PCGF6 | **17%** | bound fraction | | **A** |
| RYBP | **13%** | bound fraction | | **A** |
| RING1B^NBM (nucleosome-binding mutant) | ≈ wild type | bound fraction | imaged with endogenous RING1B present | **C** (FIGURE-ONLY) |
| RING1B catalytic domain only (RPCD) | **8%** | bound fraction | indistinguishable from HaloTag-NLS | **A** |

Residence / stable binding for the same paper uses a **completely different
acquisition**: **2 Hz, 0.5 s exposure, 600 frames**, reduced laser, tracking radius
192 nm; or **0.033 Hz, 1 s exposure, 30 frames**, radius 672 nm.
Dwell distributions fitted with a **double exponential**; **photobleaching corrected
by the H2B-HaloTag survival time** (t_bleach ≈ **50 s**, remeasured each experiment)
via t_bound = t_dwell·t_bleach/(t_bleach − t_dwell).

| Protein | Value | Quantity | Grade |
|---|---|---|---|
| RING1B | **38% of binding events are long-lived**; the long-lived τ itself is **beyond the photobleaching limit** (dwell ≈ t_bleach); some molecules visible ≥100 s, and at 0.033 Hz some survive hundreds of seconds | stable fraction *of bound events* | **A** for the 38%, **B/X** for the time (explicitly not estimable) |
| PCGF2 | **44%** of binding events long-lived; stable time also beyond photobleaching | stable fraction of bound events | **A** |
| CBX7 | **50%** of binding events long-lived; stable time also beyond photobleaching | stable fraction of bound events | **A** |
| RYBP, PCGF3, PCGF6 | stable binding **between 30 and 50 s** (i.e. measurable, unlike canonical PRC1) | stable residence | **B** (range printed; per-protein values are FIGURE-ONLY except PCGF3) |
| PCGF3 | **t_bound = 39 s** for the rare stable state (<4% of molecules), **short binding events 1.5 s**, **t_search > 1000 s** | residence + search | **A** |
| PCGF1 | occupancy estimate **0.01–0.02 molecules per kb** across PRC1 target sites | occupancy | **A** |
| RING1B, PCGF2 | **H3K27me3 depletion** (endogenous dTAG-SUZ12, dTAG-13 96 h): bound fraction **modestly but significantly reduced**; **stable-binding fraction unchanged** | H3K27me3 dependence | **C** (direction + significance in text; magnitudes FIGURE-ONLY) |
| PCGF6 | bound fraction unaffected by SUZ12 depletion | H3K27me3 dependence | **C** (FIGURE-ONLY) |
| RING1B | FRAP in Polycomb bodies: **slower and less complete recovery** than nucleoplasm (n = 20 cells each) | FRAP | **C** (curves only, no τ or immobile fraction printed) |

### 2.5 FRAP-era and non-mammalian measurements (for contrast only)

| Protein | Value | Quantity | METHOD / MODEL | Cell type | Tagging | Source | Grade |
|---|---|---|---|---|---|---|---|
| PC::GFP | **residence 0.46 s** (neuroblast interphase), **1.39 s** (SOP interphase) | residence = 1/koff | **FRAP with full reaction–diffusion model fitting** (Mueller 2008 MatLab), extracting D_f, k*on, koff simultaneously | *Drosophila* larval neuroblasts and SOP lineage, in vivo | GFP knock-in-style transgene under own promoter | Steffen/Fonseca 2012, Genes Dev, PMC3337459 | **A** |
| PH::GFP | **2–9× longer residence than PC::GFP** in all cell types; ~2× longer in SOP than neuroblast | residence | same | same | same | Fonseca 2012 | **B** |
| PC::GFP | **bound fraction 57%** (SOP) → **66–68%** (SOP daughters pIIa/pIIb), computed as k*on/(koff + k*on) | bound fraction | same | same | same | Fonseca 2012 | **A** |
| PC::GFP / PH::GFP | **0.4–2% bound** on metaphase chromatin, with residence up to **>5 min** (>300-fold longer than interphase; PC neuroblast 0.46 → 3.35 s, 7.3-fold) | mitotic binding | same | same | same | Fonseca 2012 | **A** |
| PC::GFP, PH::GFP | **2–6 min** residence, same fusions, different tissue and different FRAP model | residence | FRAP, earlier model | *Drosophila* | same fusions | Ficz 2005, as quoted by Fonseca 2012 | **C** (secondary citation) |
| CBX2/4/6/7/8 | FRAP t½ of the major mobile fraction; **t½ ~2× longer in neuronal progenitors than in undifferentiated ESC for CBX2 and CBX6**; fast-recovery fraction nearly eliminated in neuronal progenitors; t½ grows ~40% by day 6 of differentiation | FRAP mobility | confocal FRAP, triple-fraction decomposition (fast / major / non-recovering), single-exponential on the major fraction | **mouse ES cells and ESC-derived neuronal and adipocyte precursors** | fluorescent-protein CBX fusions, stable lines | Ren, Vincenz & Kerppola 2008, MCB, PMC2293085 | **C** (all absolute t½ values are in figures/supplementary tables, not in the prose) |
| BMI1, Pc, Ph | exchange half-times **20–250 s** in Polycomb bodies / polytene chromosomes | FRAP t½ | FRAP | mammalian + *Drosophila* | GFP fusions | cited in Ren 2008, PMC2293085 | **C** (secondary citation) |
| PRC1 | **~40% chromatin-bound at enriched Polycomb domains** | FRAP bound fraction | FRAP | U2OS | — | Hernandez-Munoz 2005, as cited by Youmans 2018 | **C** (secondary citation) |
| CBX2 | **>90% stably bound to mitotic chromosomes without dissociation** | mitotic binding | FRAP/live imaging | mESC | HaloTag/FP | Zhen 2014, PMC4230780, as summarised in the Ren-lab review PMC8266577 | **C** |
| PRC2 (+/− PHF1) | "multisecond residence time, midnanomolar affinity" on free DNA; PHF1 extends it | in-vitro residence | TIRF single molecule on recombinant chromatin / DNA | **in vitro** | recombinant | Choi 2017 | **X** |
| PRC1^CBX8 | FRAP **T½ = 71 ± 8 s** (GFP-CBX8) inside reconstituted PRC1–chromatin condensates | in-vitro FRAP | confocal FRAP on droplets, exponential-association fit | **in vitro** | GFP-CBX8 | Baumgart 2025, NSMB, PMC11919719 | **X** |

### 2.6 Context: PRC1 depletion changes the dynamics of the *transcription* machinery

Szczurek 2024 (Nat Cell Biol, PMC11469961) is the most recent Klose-lab SPT study.
It does **not** re-measure PRC1/PRC2 dynamics; it measures PIC components. Same
pipeline as §2.4 (Spot-On for bound fraction, low-frame-rate 1-CDF for stable
binding time, ≈100 cells/replicate for fractions, ≈20 for binding times, ≥3
replicates, endogenous HaloTag knock-ins in mESC). Direction of effect on PRC1/
H2AK119ub1 depletion (IAA): TBP bound fraction **+~50%** and binding time up;
TAF11 bound fraction and binding time up; TAF1 binding time up, fraction
unchanged; NC2β binding time sharply down; MED14 both up; TFIIB binding time up;
CDK7 unaffected; CDK9 binding time up; NELF-B binding time slightly down.
All magnitudes other than the TBP "+~50%" are **FIGURE-ONLY** (grade **C**).

---

## 3. Answers to the specific questions

**Q2 — Tatavosian 2018 Ezh2 residence and search, WT vs H3.3K27M.** Fully
recovered from Supplementary Tables 5–8 of the Nature Communications supplementary
PDF (`41467_2018_4455_MOESM1_ESM.pdf`), not from the figures. See §2.2.
Headline: **τsb 10.2 → 15.3 s and τsearch 204 → 278 s in HEK293T**;
**τsb 8.7 → 11.8 s and τsearch 188 → 265 s in 9427 → SF8628**; **F1 and F1sb
unchanged in both systems**. For Cbx7 the pattern inverts: **τsb unchanged
(9.7→9.9 s; 8.5→8.5 s)** while **F1sb falls 11.8→7.0% and 10.1→7.1%** and
**τsearch rises 73→133 s and 76→119 s**.

**Q3 — three-state analyses of PRC2.** There is **no published three-state
*binding* model for PRC2** (i.e. no vbSPT/State-Array decomposition into
free / transient-bound / stable-bound in one fit). What exists is two different
three-state things that are easy to confuse:
- Tatavosian 2018 and Kent 2020 use a **3-state *diffusion* model in Spot-On**
  (bound / confined / free). The "confined" state (D ≈ 0.5 µm²/s, 23–29% of
  molecules for Ezh2/Eed) is explicitly *not* counted as bound, and the authors
  warn it "can also be molecules that non-specifically associate with chromatin."
- The transient-vs-stable split comes from a **separate experiment** (5 Hz
  stroboscopic) fitted with a **2-component exponential**. So the "three states"
  of PRC2 (free / transient / stable) are assembled from **two incompatible
  acquisitions**, and the transient/stable fractions are rescaled onto the total
  population by multiplying by F1 from the fast acquisition.
- Youmans 2018 used a **2-state** Spot-On model and reports **no state
  decomposition of binding at all**.
- No PRC2 paper found uses **vbSPT**, **State Array / saSPT**, or any Bayesian
  diffusion-spectrum method. Huseyin & Klose explicitly tested 2-state vs 3-state
  Spot-On and rejected 2-state for PRC1.

**Q4 — diffusion coefficients and search times.** D_free: EZH2 **2.09 ± 0.12**,
SUZ12 **1.99 ± 0.13** µm²/s (U2OS, Spot-On 2-state) vs Ezh2 **4.2 ± 0.4**, Eed
**4.4 ± 1.2** µm²/s (mESC, Spot-On 3-state fastest component). **These differ by
~2×, and essentially all of that is the model, not the biology** — in a 3-state
fit the "free" component excludes a 23–29% confined population at ~0.5 µm²/s that
a 2-state fit folds into D_free and drags downward. Search times: Ezh2 **210 ± 19 s**
(mESC), 204 ± 40 s (HEK293T), 188 ± 12 s (9427); Eed 214 ± 6 s; Cbx7 **64 ± 5 s**
(mESC), 73 ± 3 s (HEK293T); CBX2 ~72 s (mESC); PCGF3 **>1000 s** (Klose, different
formalism). **Tatavosian's own framing: τsearch for Ezh2 is ~20× its residence
time** — PRC2 spends ~95% of its life searching.

**Q5 — perturbations.**
- *H3K27me3 loss (EED-side):* Cbx7 bound fraction collapses below fit convergence
  and Cbx8 falls 26→8% in Eed−/− mESC (Zhen 2016); CBX2 τsb unchanged (8.9→9.7 s)
  but F1sb falls to 70% and search lengthens (Kent 2020); RING1B/PCGF2 bound
  fraction "modestly but significantly" reduced with stable fraction unchanged
  after SUZ12 degron (Huseyin 2021). **EZH2's own bound fraction is unchanged by
  H3.3K27M-driven ~7-fold H3K27me3 loss** (Tatavosian 2018).
- *EED inhibitor A-395 (U2OS, Youmans 2018):* 30 min — **no change** in bound
  fraction, **no change** in residence; 72 h — residence **decreases** (significant
  drop in the fraction still bound at 3 s, P < 0.05) which the authors attribute
  to possible secondary effects. **No fitted τ is reported for either condition.**
- *Accessory-subunit mutants:* SUZ12(PHF1+ABHmut) and SUZ12(VEFS) each cut the
  EZH2 bound fraction significantly (P < 0.006), down to **0.08**, which the
  authors say is at the HaloTag-alone background of ~0.10 (Youmans 2018).
  Ezh2^ΔSANT2 (cannot bind Suz12) in mESC: F1 unchanged at 22%, but **F1sb halves
  4.6→1.9%, τsb halves 10.0→5.6 s, τsearch rises 210→287 s** (Tatavosian 2018,
  Suppl. Table 4). Eed−/−: **τsb 10.0→6.9 s, τsearch 210→138 s, F1sb unchanged**
  (the two changes cancel) — a clean demonstration that bound fraction alone is
  blind to real kinetic change.
- *Cbx7 overexpression in DIPG SF8628 (dox titration):* **τsb rises 8.5 → 13.2 →
  16.2 s** and **τsearch 119 → 283 → 369 s** across 0.1 / 1.0 / 2.0 µg/ml dox
  (Tatavosian 2018, Suppl. Table 10). **Expression level by itself moves the
  residence time 2-fold.** This is a serious caveat for every
  overexpression-based number in this document.

**Q6 — mESC and neural measurements.** mESC: Zhen 2016, Tatavosian 2018,
Kent 2020, Huseyin & Klose 2021, Szczurek 2024, Yang 2026 — all mESC.
**Neural:** only three sources, and all three are weak numerically —
(i) Yang/Ren 2026 CBX2 in mouse NPCs (endogenous knock-in, fold-changes only,
absolute values figure-only); (ii) Ren 2008 FRAP of CBX2/CBX6 in ESC-derived
neuronal progenitors (t½ ~2× longer than ESC, absolute values in figures);
(iii) Tatavosian 2018 in DIPG line SF8628 and control brain-tumour line 9427
(full Supplementary Table 6, but transformed glioma lines, not normal neural
cells, and overexpressed fusions). **There is no published EZH2/SUZ12/EED
residence time or bound fraction in normal neural progenitors or neurons.**

**Q7 — FRAP vs SMT disagreements.** See §4.

---

## 4. Where FRAP and SMT disagree, and by how much

1. **Polycomb residence, FRAP vs FRAP, same fusion proteins:** Fonseca 2012
   measured PC::GFP and PH::GFP interphase residence at **0.5–10 s**; Ficz 2005
   measured the **same fusions at 2–6 min**. That is a **~30–100-fold
   disagreement** between two FRAP studies, which Fonseca attributes to cell type
   *and* to the FRAP model used. This alone sets the credibility ceiling on any
   single FRAP-derived residence time.
2. **Polycomb residence, FRAP vs SMT:** Fonseca's FRAP τ of 0.46–1.4 s for PC::GFP
   is **5–20× shorter** than the SMT τsb of 5.5–9.4 s for mammalian CBX/RING1B/PHC1
   (Zhen 2016, Tatavosian 2018, Kent 2020). Different organism, but the direction
   is the systematic one: **reaction–diffusion FRAP fits return short residence
   times, SMT dwell-time fits return long ones**, because SMT's D threshold
   pre-selects the slowest molecules and bleach correction inflates the tail.
3. **PRC2 bound fraction, FRAP vs SMT, same paper, same cells:** Youmans 2018 gets
   **0.20 genome-wide by SMT** but **~0.40 immobile by FRAP at the H3K27me3-rich
   lacO array**. The authors read this as a real locus effect, but note the two
   numbers are not the same quantity: FRAP "immobile fraction" is whatever has not
   recovered within the 5-minute window at one spot, SMT "bound fraction" is a
   displacement-model state over the whole nucleus. **A factor of 2 that is partly
   locus and partly definition, with no way to separate them from the published data.**
4. **PRC1 bound fraction, FRAP vs SMT:** the older FRAP literature put PRC1 at
   **~40% bound** at Polycomb domains (Hernandez-Munoz 2005); Huseyin & Klose get
   **20% for RING1B, 18% for CBX7 and PCGF2** by SMT genome-wide, and their own
   FRAP shows slower, incomplete recovery in Polycomb bodies without a quantified
   fraction. **~2-fold disagreement**, again confounded by genome-wide vs
   locus-restricted.
5. **PRC2 residence, in vivo SMT vs in vitro TIRF, same paper:** Tatavosian 2018
   gets **τsb ≈ 10 s in cells** and **τ ≈ 257 s on reconstituted nucleosomes**
   (uncorrected for bleaching, so a lower bound). **25-fold.** This is the single
   largest discrepancy in the PRC2 literature, and it cannot be dismissed as a
   fitting artefact in either direction.
6. **Same lab, same cells, same nominal protocol, Cbx7 τsb:** 7.3 s (2016),
   9.4 s (2018), 6.9 s (2020). **±30%.** The 2016 paper used a bound threshold of
   D_m < 0.10 µm²/s; 2018 and 2020 used D_m ≤ 0.032 µm²/s. **Same lab, same cells,
   threshold changed 3-fold mid-stream, and this is never reconciled in print.**
7. **Same lab, same cells, Cbx7 bound fraction:** 29–30% (Zhen 2016, Gaussian
   mixture with the bound peak pinned at 0.032 µm²/s) vs **40%** (Tatavosian 2018,
   Spot-On 3-state with D_bound ≤ 0.01 µm²/s). **A 1.4× difference produced purely
   by swapping the fitting framework.** Klose's endogenous knock-in CBX7 in the
   same cell type gives **18%**, i.e. **2.2× lower than the Ren-lab transgene
   value** — overexpression vs knock-in is the obvious suspect, and Tatavosian's
   own dox titration (τsb 8.5→16.2 s across a 20-fold expression range) shows the
   effect is real.

---

## 5. Which numbers are mutually comparable

### Comparable sets (compare freely within a set)

**Set A — Spot-On 3-state, D_bound ≤ 0.01 µm²/s, 30 ms lag, dZ 0.7, LocError 0.035
(Ren lab fast-acquisition).** Tatavosian 2018 all rows (Ezh2 23%, Eed 20%,
Cbx7 40%, Ezh2^Catalytic 17%, all HEK293T/9427/SF8628 rows) + Kent 2020 F1 values
+ Yang 2026 F1 values. Within-set comparisons are valid **provided the exposure
time matches** — the 5/10/20/30 ms series moves F1 by 5–7 points.

**Set B — Spot-On 3-state, localisation error 40 nm, 15 ms exposure, 67 Hz,
dZ 0.7, endogenous knock-ins (Klose lab).** RING1B 20%, PCGF2 18%, CBX7 18%,
PCGF1 26%, PCGF3 13%, PCGF6 17%, RYBP 13%, RPCD 8%, H2B 59%, NLS 11%, plus the
dTAG-SUZ12 perturbation. **Internally the cleanest dataset in the field**: one
pipeline, one cell line, endogenous expression, shared controls.

**Set C — Ren-lab 5 Hz residence pipeline, D_m ≤ 0.032 µm²/s, bleach-corrected
2-exponential.** Tatavosian 2018 (all τtb/τsb/F1sb/τsearch), Kent 2020
(CBX2 8.9 s, RING1B 5.5, MEL18 5.7, PHC1 5.6, CBX6 2.7, CBX7 6.9), Yang 2026.
Ezh2 10.0 s, Eed 9.5 s, Cbx7 9.4 s, CBX2 8.9 s, RING1B 5.5 s are a genuinely
comparable ladder.

**Set D — Youmans 2018 internal.** EZH2 0.20 vs SUZ12 0.25 vs the SUZ12-mutant
0.08; and D_free 2.09 vs 1.99. Valid internally; **none of it transfers to Sets A–C.**

**Set E — Fonseca 2012 internal.** PC vs PH, neuroblast vs SOP vs pIIa/pIIb,
interphase vs metaphase. One model, one organism, one lab.

### Non-comparable pairs (do **not** put these in the same plot)

| Pair | Why not |
|---|---|
| Youmans EZH2 0.20 vs Tatavosian Ezh2 0.23 | 2-state vs 3-state; D_bound ceiling 0.08 vs 0.01 µm²/s (8-fold); 97.5 Hz vs 33 Hz; U2OS vs mESC; knock-in vs overexpression. **The near-agreement is a coincidence of four offsetting differences, not a replication.** |
| Youmans EZH2 0.20 vs Klose RING1B 0.20 | different complexes, different labs, different exposure (10.3 vs 15 ms), different LocError (35 vs 40 nm), 2-state vs 3-state. Equal numbers, unequal quantities. |
| Zhen Cbx7 29% vs Tatavosian Cbx7 40% vs Klose CBX7 18% | three fitting frameworks, two expression regimes. **Spread is 2.2-fold and is entirely methodological.** |
| Any Spot-On F1 vs any F1sb | F1 includes transient non-specific binding; F1sb is the long-lived subset. For Ezh2 these are **23% vs 4.6% — a 5-fold difference in the same cells in the same paper.** Papers that report only one of the two are not comparable to papers that report only the other. |
| Youmans "survival probability at 3 s" vs any τsb | Youmans reports **no fitted τ and no bleach correction**. Reading a residence time off that curve is not supported by the paper. |
| Klose "38% stable" vs Tatavosian "F1sb 4.6%" | Klose's 38% is a fraction **of bound events**; Tatavosian's 4.6% is a fraction **of all molecules**. Klose's comparable quantity is 0.38 × 0.20 ≈ **7.6% of all RING1B molecules**; Tatavosian's comparable quantity for Cbx7 is f1sb = 12.9/40 ≈ **32% of bound molecules**. |
| Klose stable residence vs Ren-lab τsb | Klose's canonical-PRC1 stable binding is **beyond the photobleaching limit** (t_bleach ≈ 50 s) and explicitly not estimated; Ren-lab τsb of 5–10 s comes from a bleach-corrected fit at 5 Hz. The Klose data say PRC1 stable binding is **longer than 50 s**; the Ren data say 5.5 s. These are in open conflict and the papers do not address it. |
| Ferrie f_bound 0.55 (State Array, D < 0.1 µm²/s, 223 Hz) vs anything above | State Array is a Bayesian spectrum with a **post-hoc** D cut-off of 0.1 µm²/s — 10× looser than Tatavosian's Spot-On constraint and ~1.25× looser than Youmans'. A looser bound-D ceiling mechanically raises f_bound. |
| Any live-cell τ vs Tatavosian's 257 s in vitro or Baumgart's 71 s droplet FRAP | different system entirely. |

### The three parameters that drive most of the spread

1. **The D ceiling that defines "bound"**: 0.01 (Tatavosian Spot-On) < 0.032
   (Ren residence) < 0.08 (Youmans) < 0.10 (Zhen residence; Ferrie) µm²/s. **A
   10-fold range.** A looser ceiling pulls confined molecules into "bound" and
   inflates both fraction and dwell time.
2. **Number of states**: a 2-state fit has no home for the ~25% confined population
   at D ≈ 0.5 µm²/s, so it redistributes it — raising the bound fraction, lowering
   D_free, or both.
3. **Tagging and expression**: endogenous knock-in (Youmans, Klose, Yang 2026) vs
   tet-inducible lentiviral overexpression (Zhen, Tatavosian, Kent). Tatavosian's
   own dox titration shows **τsb doubling across a 20-fold expression range**.

### Practical recommendation for parameterising a model

For PRC2 in a mouse neural/stem context the defensible numbers are
**Tatavosian 2018 Supplementary Table 2** (mESC, grade A, complete method stack):
**τ_specific ≈ 10 s, τ_nonspecific ≈ 0.9 s, specifically-bound fraction ≈ 4.6%,
total bound ≈ 23%, τ_search ≈ 210 s, D_free ≈ 4.2 µm²/s.** Carry the caveats that
these are overexpressed fusions, that the bound/confined boundary is set by an
imposed D_bound ≤ 0.01 µm²/s, and that the same lab's own exposure-time series
shows F1 is inflated 5–7 points by motion blur even at 5 ms. If the model needs
K27M-perturbed values, use the SF8628/9427 pair (τsb 8.7 → 11.8 s, τsearch
188 → 265 s) rather than the HEK293T pair, since those are at least brain-tumour
lines with an endogenous heterozygous mutation.

---

## 6. Sources

| Short name | Full reference | Access |
|---|---|---|
| Youmans 2018 | Youmans DT, Schmidt JC, Cech TR. Genes Dev 32:794–805 | PMC6049511 |
| Tatavosian 2018 | Tatavosian R, Duc HN, Huynh TN, et al. (Ren X). Nat Commun 9:2080 | PMC5970213; suppl. PDF `41467_2018_4455_MOESM1_ESM.pdf` (Supplementary Tables 1–12) |
| Zhen 2016 | Zhen CY, Tatavosian R, Huynh TN, et al. (Ren X). eLife 5:e17667 | PMC5056789 |
| Kent 2020 | Kent S, Brown K, Yang CH, et al. (Ren X). Cell Rep 33:108248 | PMC7593837 |
| Huseyin & Klose 2021 | Huseyin MK, Klose RJ. Nat Commun 12:887 | PMC7873255 |
| Szczurek 2024 | Szczurek A, et al. (Klose RJ). Nat Cell Biol | PMC11469961 |
| Yang/Ren 2026 | Self-clustering of three CBX2 molecules drives PRC2… Mol Cell | PMC12981346 |
| Ren 2008 | Ren X, Vincenz C, Kerppola TK. Mol Cell Biol | PMC2293085 |
| Fonseca 2012 | Fonseca JP, Steffen PA, et al. (Ringrose L). Genes Dev 26:857–871 | PMC3337459 |
| Ren-lab review 2021 | Single-molecule imaging of epigenetic complexes in living cells. Nucleic Acids Res | PMC8266577 |
| Zhen 2014 | Cbx2 stably associates with mitotic chromosomes… Mol Biol Cell | PMC4230780 |
| Baumgart 2025 | Dynamic PRC1-CBX8 stabilizes a porous structure of chromatin condensates. NSMB | PMC11919719 (in vitro) |
| Ferrie 2023 | p300 State Array analysis (supplied by requester) | PMC10245728 |
