# Enzyme occupancy and binding kinetics — what the biochemistry actually says
2026-10-23. For replacing the lattice's continuous "PRC2 activity field" with
discrete, countable complexes.

Values are from the primary sources where I could reach them; **everything I
could not verify is marked as such** rather than filled in from memory.

---

## 1. PRC2 in live cells — VERIFIED

[Youmans, Schmidt & Cech 2018, Genes Dev 32:794](https://pmc.ncbi.nlm.nih.gov/articles/PMC6049511/)
(HaloTag knocked into the **endogenous** EZH2 and SUZ12 loci, U2OS, single-particle tracking):

| quantity | value |
|---|---|
| chromatin-bound fraction | **~20%** (80% freely diffusing) |
| D_free, EZH2 | **2.09 ± 0.12 µm²/s** |
| D_free, SUZ12 | **1.99 ± 0.13 µm²/s** |
| residence time | **seconds** — scored as the fraction still bound after 3 s; no single tau reported |
| EED–K27me3 inhibitor, 30 min | **no change** in residence time |
| EED–K27me3 inhibitor, 72 h | residence time **decreases** |
| at an H3K27me3-rich locus (lacO array) | larger immobile fraction, slower recovery than bulk nuclear sites |

The 30-min-vs-72-h split matters: blocking the EED cage does not immediately
change how long PRC2 stays on chromatin. Only after the MARK is depleted does
residence drop. So the read-write feedback acts on *where PRC2 ends up over
hours*, not on its instantaneous dwell time.

## 2. PRC2 in vitro — VERIFIED

- Binding is **DNA-dominated, not H3K27me3-dominated**
  ([Choi 2017 NSMB 24:1039](https://www.nature.com/articles/nsmb.3488); Wang 2017).
  PRC2 binds free DNA with **multisecond residence and mid-nanomolar affinity**.
  Adding H3K27me3 to nucleosomes has little effect on affinity; increasing
  linker length does.
- **PHF1 prolongs residence ~2–3x** (Choi 2017) — accessory subunits, not the
  mark, are the main affinity determinant. Consistent with Youmans: removing
  AEBP2/PCL binding strongly depletes the bound fraction.
- Occupancy titration ([Leicher 2020 PNAS](https://pmc.ncbi.nlm.nih.gov/articles/PMC7720148/),
  optical tweezers on 12-mer arrays): PRC2-bridged states saturate between
  **5 and 100 nM**; at **500 nM the majority of nucleosomes are PRC2-bound**.
  Bridged (distal) complexes rupture at **21.4 ± 1.0 pN** vs 16.6 ± 0.3 pN for
  ordinary mononucleosome unwrapping — the bridge is the mechanically stable mode.
- **Allosteric stimulation is only up to ~3-fold** on nucleosomal substrate.
  Note our v44 module and the Jost model both assume **5–10x**. That is a
  discrepancy worth resolving.
- me2→me3 is the slowest catalytic step (agrees with the lattice's `k_me3`).
- Domain formation in vivo takes **~36 h to reach steady state** after EED
  induction (Oksuz et al.), spreading outward from nucleation sites.

## 3. p300 — VERIFIED

[Ferrie et al. 2023, Mol Cell](https://pmc.ncbi.nlm.nih.gov/articles/PMC10245728/)
(endogenous knock-in, U2OS):

| quantity | value |
|---|---|
| bound fraction | **0.54** (0.55 ± 0.02 and 0.52 ± 0.03 in two clones) |
| residence time | **~26 s** (FRAP, not SMT) |

p300 is **~2.7x more chromatin-bound than PRC2** (0.54 vs 0.20) and stays
~10x longer. Caveat: FRAP and single-molecule estimates of residence can differ
by 2–4x (see the live-cell benchmark literature), so 26 s is soft.

## 4. Abundance — NOW RESOLVED (2026-10-23)

Source: **PaxDb v6.1 bulk files** (`pax-db.org/downloads/latest/datasets/9606/`),
whole-organism integrated (19,483 proteins) and HeLa Geiger 2012. The web UI is
JS-rendered and unusable; the flat files are not.

### Method, and why it is trustworthy

PaxDb reports **ppm, not copies**. Converting via total protein per cell is
unreliable — it gave cohesin ~37,500 against a measured 250,000. So instead
each enzyme is normalised **against histone H3 in the same dataset**, which
cancels the total-protein assumption (nucleosomes = H3/2).

**Validation.** Cohesin is the one chromatin complex measured absolutely
([Holzmann 2019 eLife](https://elifesciences.org/articles/46269v1): ~250,000
nuclear complexes in G1), i.e. 1 per 120 nucleosomes.

| | cohesin per nucleosome |
|---|---|
| measured (Holzmann) | 1 per **120** |
| PaxDb ratio, whole-organism | 1 per **217** |
| PaxDb ratio, HeLa Geiger | 1 per **168** |

Agreement within ~1.5x. The ratio method is sound; the earlier shortfall was the
total-protein assumption, not PaxDb. Copy numbers below are anchored on RAD21 =
250,000 cohesin.

### The numbers

| protein | ppm (WO) | ppm (HeLa) | copies/cell | nucleosomes per molecule |
|---|---|---|---|---|
| **EZH2** | 1.81 | 2.39 | **15,000–36,000** | **830–2,000** |
| EZH1 | 0.295 | – | 5,900 | 5,100 |
| SUZ12 | 3.33 | 10.9 | 67,000 | 450 |
| EED | 3.48 | 9.95 | 61,000–70,000 | 430–490 |
| **EP300** | 0.735 | 0.892 | **5,500–15,000** | **2,000–5,500** |
| CREBBP | 0.616 | 2.17 | 12,000–13,000 | 2,300–2,400 |
| KDM6A | 0.409 | – | 8,200 | 3,700 |
| **KDM6B** | 0.241 | – | **4,800** | **6,200** |
| HDAC1 | 34.8 | 289 | 700,000–1,800,000 | 17–43 |
| HDAC2 | 33.5 | 458 | 670,000–2,800,000 | 11–45 |

**EZH2 is the limiting PRC2 subunit** — SUZ12 and EED are present at ~2–4x
more, so most SUZ12/EED is not in a catalytically complete PRC2. The number of
PRC2 holocomplexes is set by EZH2(+EZH1): **~20,000–40,000 per cell**.

## 5. The number that matters for the model

Combining abundance with the **measured** bound fractions:

| | copies | bound fraction | bound molecules | per 50-nucleosome domain if concentrated at 3,000 targets |
|---|---|---|---|---|
| PRC2 (EZH2) | 36,000 | 0.20 (Youmans) | 7,240 | **2.4** |
| p300 | 14,700 | 0.54 (Ferrie) | 7,938 | **2.7** |
| KDM6B | 4,800 | unknown | ~1,000 if 0.2 | **~0.3** |

Spread uniformly over the genome instead, every one of these falls to
**~0.01 per domain**. At 10,000 target elements rather than 3,000 they fall to
~0.7.

**So: of order ONE enzyme molecule on a target promoter at any instant — and
for KDM6B, less than one.** Chromatin state is maintained by rare, brief visits,
not by resident enzymes. HDACs are the exception, at 1 per 11–45 nucleosomes,
which is three orders of magnitude denser than the writers and erasers.

This is decisive for the rewrite: `prc2_activity()` currently evaluates a
continuous PRC2 field at all 51 nucleosomes simultaneously. The biochemistry
says the correct object is **one or two discrete particles per domain**, binding
for seconds, mostly doing nothing.

## 6. What is still missing

| needed | status |
|---|---|
| EZH2 / SUZ12 / EED copy number | **RESOLVED** — PaxDb v6.1, histone-normalised, cohesin-validated (section 4) |
| EP300, KDM6B, KDM6A copy number | **RESOLVED** (section 4) |
| per-cell-type abundance (GNP vs MB) | not done — PaxDb whole-organism is a cross-tissue average |
| histone MS normalisation | H3:H4:H2B come out 1:2:3 in ppm when they should be 1:1:1, so the nucleosome denominator carries ~2-3x uncertainty |
| PRC2 k_cat on nucleosomes | **not found.** Pivotal: decides whether a 3 s visit can deposit a methyl at all. Qualitatively PRC2 is described as inefficient, and me2->me3 slowest, but no turnover number surfaced |
| KDM6B / UTX binding kinetics | **nothing found** — the eraser arm has no kinetic data at all, which is awkward since the lattice says KDM6B is a dominant antagonist |
| HDAC residence / rate | not searched |

## 7. Design implications for the rewrite

1. **Model enzymes as particles, not fields.** A small integer number of PRC2
   complexes per domain, each binding, sliding/hopping, occasionally catalysing,
   and unbinding. This is Owen/Mirny's enzyme-limitation regime and it is where
   the biochemistry actually sits.
2. **Separate binding from catalysis.** The measured residence (seconds) and the
   measured domain formation time (~36 h) differ by ~10^4, so most visits must be
   unproductive. The current model conflates "PRC2 active here" with "PRC2 writes
   here" — they are different events with very different rates.
3. **Recruitment, not the mark, sets occupancy.** Both Youmans (accessory-subunit
   mutants) and Choi/Wang (DNA dominates) say the same thing. The read-write
   feedback should act on catalytic efficiency once bound, NOT on binding
   affinity — which is the opposite of how `a_rw_prc2` enters now.
4. **Reduce the allosteric boost from 5–10x to ~3x** unless a better source
   turns up, and re-check whether the bistability survives — it may not, which
   would be a result in itself.
5. **p300 is the more chromatin-resident enzyme** (0.54 bound, ~26 s vs 0.20,
   ~3 s). The model currently treats p300 as a weak static profile. If anything
   the occupancy asymmetry runs the other way.
