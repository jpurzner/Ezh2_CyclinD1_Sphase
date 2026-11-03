# Consolidated parameter table — census in progress
2026-10-23. **Supersedes `parameters_from_primary_sources_2026-10-23.md` and
`locked_parameters_2026-10-23.md`, both of which contain numbers now known to be
wrong.** Sources: six parallel literature censuses (`docs/lit_*.md`) plus my own
primary reading.

**STATUS: 3 of 6 censuses reported (in-vivo dynamics, mark turnover, abundance).
Catalysis, in-vitro binding and antagonist enzymes still outstanding.** k_cat in
particular is still unknown, so section 6 is a prediction awaiting test.

Grades: A = direct measurement in primary text · B = primary but model-derived or
label-free · C = secondary/second-hand · X = does not exist.

---

## 1. Corrections to what I previously told you

| claim | status |
|---|---|
| "cohesin anchor = 250,000" | **WRONG** — that is nuclear SCC1-GFP by FCS. Total cellular RAD21 = 350,000 (SRM). |
| "PaxDb under-reports RAD21 5.2x, so scale EZH2 up by that" | **WRONG METHOD** — Itzhak (same histone-ruler method) over-reports RAD21 1.43x. The 5.2x was a PaxDb integration artefact and I propagated it onto EZH2. |
| "EZH2 ~50,700 copies/cell" | **WRONG, 2.8x low** — Itzhak direct read gives 142,635. |
| "p300 is 2.7x more chromatin-bound than PRC2" | **RETRACTED** — different algorithms (Spot-On 2-state vs State Array spectrum), not comparable. |
| "PRC2 residence is 'a few seconds', no value exists" | **SUPERSEDED** — tau_stable = 10.0 ± 0.2 s (mESC). But see the endogenous caveat below. |
| "H3K27me3 persists multiple cell divisions" | **PARTLY WRONG** — true for me3 (4.7 cycles), false for me1/me2 (~1.5 cycles). |
| "domain formation reaches steady state at ~36 h" | **OVERSTATED** — 36 h is Oksuz's last time point, not a fitted plateau. No time constant published. |
| "nucleosomes = 3.0e7" | **human only. MOUSE = 2.7e7** (5.4 Gb genome). Our system is mouse. |
| "the proteomic ruler independently confirms the nucleosome count" | **WRONG — it is circular**; it assumes histone mass = DNA mass. |

## 2. Abundance (grade B — no grade A exists for any PRC2 subunit)

**There is no SRM/PRM or recombinant-standard quantification of any PRC2 subunit
in the literature.** PRC2 has no Holzmann-2019 equivalent. Everything is
label-free. Source below: Itzhak 2016 eLife 5:e16950, HeLa, proteomic ruler,
median of 6 replicate organellar maps.

| protein | copies/cell | per nucleosome (mouse, 2.7e7) |
|---|---|---|
| **EZH2** | 142,635 | 1 per 189 |
| EZH1 | 9,699 | 1 per 2,784 |
| SUZ12 | 194,816 | 1 per 139 |
| EED | 351,948 | 1 per 77 |
| EP300 | 27,961 | 1 per 966 |
| CREBBP | 54,919 | 1 per 492 |
| KDM6A | 4,183 | 1 per 6,455 |
| **KDM6B** | **not detected** | — |
| HDAC1 | 461,170 | 1 per 59 |
| HDAC2 | 1,542,148 | 1 per 18 |

Assembled core is strictly 1:1:1 (Smits 2013). **Use EZH1+EZH2 to set the number
of competent complexes** = **152,334**; SUZ12/EED excess is partly iBAQ bias.
Uncertainty band: ~3x. EPOP can be on at most ~4.5% of PRC2 (it is 22x less
abundant than EZH2) — the AP-MS figure of 0.47 per core is off ~10x.

## 3. In vivo dynamics (grade A for fractions, B for tau)

| quantity | value | source + method |
|---|---|---|
| EZH2 chromatin-bound | **20% ± 3%** | Youmans 2018, **endogenous knock-in**, Spot-On 2-state, D_bound ∈ [0.0005, 0.08] |
| SUZ12 chromatin-bound | 25% ± 3% | same |
| EZH2 total bound | 23% | Tatavosian 2018, mESC, Spot-On **3-state**, D_bound ∈ [0.0001, **0.01**] |
| EZH2 **stably** bound | **4.6% ± 0.3%** | Tatavosian Suppl. Table 2 |
| tau_stable | **10.0 ± 0.2 s** | Tatavosian, 2-exponential, bleach-corrected |
| tau_transient | 0.9 s | same |
| tau_search | 210 ± 19 s | same |
| D_free | 4.2 µm²/s (mESC) ; 2.09 (U2OS) | Tatavosian ; Youmans |

**Internal consistency check that passes:** a molecule bound 10 s and searching
210 s should be bound 10/220 = 4.5% of the time. Measured stable fraction 4.6%.
Independent agreement.

**Caveats that matter.** Youmans' 20% and Tatavosian's 23% look like replication
but the D_bound constraints differ 8x and the models differ (2- vs 3-state).
Tatavosian's values are **overexpressed fusions**, and the same paper shows
residence doubles with expression level (Cbx7 8.5 → 16.2 s under dox titration).
**No endogenous-expression PRC2 residence time exists.** No PRC2 bound fraction
or residence time has ever been measured in normal neural progenitors.

## 4. Mark turnover (grade A)

Zheng 2014 Table 1, HeLa, dual-SILAC chase on old histones, dilution normalised out:

| mark | t1/2 | loss rate |
|---|---|---|
| **K27me3** | **158 ± 23 h** | 0.0044 /h |
| K27me2 | 57 ± 20 h | 0.0122 /h |
| K27me1 | 52 ± 16 h | 0.0133 /h |
| K4me3 | 6.8 h | 0.102 /h |
| K36me3 | 57 h | 0.0122 /h |

**In vivo loss is essentially all replicative dilution**, by three independent
routes: Reverón-Gómez 2018 (parental landscape unchanged 24 h across mitosis
under EZH2i; only 3.5% of regions lose moderately), Jadhav 2020 (17 d of decay
in Eed-/- fits a pure 50%-per-division model), Alabert 2020 (the winning domain
model contains **no demethylation terms**; a demethylation model is rejected by BIC).

Zheng's **forward** constants are bulk averages over all nucleosomes, most of
which never see PRC2, so they do NOT transfer to a PcG target domain. The loss
rate does.

**H3K27 state stoichiometry:** me0 10–20%, me1 10–15%, **me2 50–70%**, me3 5–15%.
me2 is the dominant state, not a transit state.

## 5. Consequence already tested in the lattice

Our `del_mk` 0.06 /h + exchange 0.06 /h = 0.12 /h is **27x faster** than the
measured 0.0044 /h. Correcting it:

| | core | plateau | depletion | 20–80% width |
|---|---|---|---|---|
| del_mk 0.060 (ours) | 0.049 | 0.457 | 0.408 | 1900 bp |
| **del_mk 0.004 (measured)** | 0.121 | **0.835** | **0.714** | **950 bp** |

The boundary becomes **2x sharper and 1.7x deeper**, with a more plausible
plateau. The mean-field **bistability disappears entirely** at realistic turnover
— it was an artefact of the inflated rate. The model is better with the right number.

Uncertainty now concentrates in transcription-coupled exchange, which becomes the
only significant turnover term: raising k_h33 0.15 → 0.4 collapses the domain
(plateau 0.835 → 0.024). **No absolute mammalian nucleosome exchange rate exists**
— only relative indices.

## 6. PREDICTION awaiting the catalysis census

From measured abundance + bound fraction + stoichiometry + half-life:

- 66–230 marked nucleosomes per chromatin-engaged PRC2
- maintenance demand 50,000–177,000 me3/h genome-wide (dilution : turnover ≈ 5:1)
- = 1.7–5.8 me3 per engaged PRC2 per hour
- a PRC2 is stably bound only 166 s/h
- **=> required effective k_cat while bound: 0.6–2.1 per minute**

If the measured k_cat is far below this, then either PRC2 works mainly in the
transient mode, or the stable fraction is underestimated, or the maintenance
picture is wrong. This is the test to run when the catalysis census lands.

## 7. Known gaps (3 censuses still out)

- PRC2 k_cat on nucleosomes — **pending**
- KDM6B abundance (not detected by Itzhak) and any KDM6B kinetics — **pending**
- No kb/h spreading velocity for H3K27me3 exists; only H3K9me3/HP1a at 0.036 kb/h
- No absolute mammalian nucleosome exchange rate
- No CGNP / medulloblastoma copy numbers
- No endogenous-expression PRC2 residence time
- Retrievable and worth chasing: Spruijt 2013 supplementary (mouse, neural, nuclear,
  absolute Eed/Ezh2/Suz12), Javasky 2018 Fig 2B (cell-cycle-resolved absolute
  H3K27me1/2/3), Leroy 2013 (24 cancer lines, K27/K36 states)

---

# CENSUS COMPLETE — all 6 reported (2026-10-23)

Per-topic detail in `lit_prc2_catalysis.md`, `lit_prc2_binding_invitro.md`,
`lit_prc2_invivo_dynamics.md`, `lit_h3k27_turnover.md`,
`lit_antagonist_enzymes.md`, `lit_abundance_stoichiometry.md`.

## 8. PRC2 k_cat on nucleosomes — EXISTS, in 5 independent datasets

My four searches missed it; the census found it. All from primary Tables.

| source | complex / substrate | kcat (min⁻¹) | K_M |
|---|---|---|---|
| Ge/Muir JACS 2019 | core 4-mer, 12-mer 601 array | 0.00096 | 120 nM |
| Jani PNAS 2019 | 5-mer +AEBP2, 12-mer array | 0.0057 | 89 nM/site |
| Sneeringer PNAS 2010 | 5-mer, oligonucleosome | 0.0071 | 141 nM |
| Nat Commun 2021 | core 4-mer, ATOH1 chromatin, 35 mM KCl | 0.034 | 111 nM |
| EMBO J 2025 | 6-mer AEBP2ᴸ-JARID2, 100 mM KCl | 0.071 | 38 nM |
| McCabe PNAS 2012 | 5-mer, HeLa mono/dinucleosome | 0.11 | 80 nM |
| **Nat Commun 2021** | **+PALI1 PIR, 100 mM KCl** | **0.474** | 392 nM |

**0.06–28 methyl transfers per PRC2 per hour; ~480-fold spread.** Do not average —
match subunit composition and salt to the biology. Most are formally kcat,app
(sub-saturating SAM), so they are underestimates.

## 9. THE PREDICTION WAS TESTED AND FAILED — productively

I predicted 0.6–2.1 /min required. **Even the fastest measured value (0.474) is
1.3x short; typical values are 10–600x short.** The assumption that broke was
mine: that only the *stably* bound 4.6% catalyses.

| catalysing population | bound time | required kcat |
|---|---|---|
| stable only (4.6%) | 166 s/h | 0.62 – 2.10 /min — **exceeds all measurements** |
| stable + transient (23%) | 828 s/h | **0.12 – 0.42 /min — inside the measured range** |

**Falsifiable claim: PRC2 must deposit marks during sub-second transient
encounters, not only during 10 s stable binding.** At kcat 0.4 /min a 0.9 s visit
is productive with p ≈ 0.006 (1 in ~170). A model in which only stably-bound
PRC2 writes cannot maintain H3K27me3 at measured abundance and turnover.

## 10. Corrections to my own earlier statements

| claim | status |
|---|---|
| "allosteric stimulation ~3-fold" | **WRONG, too low.** Margueron 2009 Fig 4A: **~7x** on oligonucleosomes, effect on kcat with K_M unchanged. Lee/Armache 2018: **12x for EZH2, 2x for EZH1**. The 3x I had is probably Ge/Muir's *distal in-cis* value. |
| "read-write should modulate catalysis not binding; a_rw does the opposite" | **I WAS WRONG ABOUT OUR OWN CODE.** `a_rw_prc2` feeds `P`, which scales k_me1/2/3 — i.e. catalysis. That is exactly right, and is confirmed: H3K27me3 affects catalysis, **not binding** (null result across Wang 2017, Cookis 2025, Finogenova 2020 — no paired K_d shift exists anywhere). |
| "Leicher: majority of nucleosomes bound at 500 nM" | **Not a measurement** — their justification for a working concentration, cited to Wang 2017. Measured: bridging saturates ~100 nM, at 200 mM KCl. |
| "no KDM6B kcat exists" | **WRONG** — Walport 2016 Nat Commun Table 1: 0.70 /min, K_M 6.7 µM. Missed because the paper is indexed as an *arginine* demethylation study. |

## 11. Spreading range — the kernel is measurable after all

Ge/Muir JACS 2019, H3K27me3 *in cis* on 12-mer arrays: **~8.8x proximal, ~3.2x
distal, spreading restricted to 1–4 nucleosomes, no stimulation in trans.**

This matters because the geometry analysis showed kernel locality (`p_local`) is
**not identifiable from ChIP profiles** (sv 5.7, near-degenerate). In vitro
supplies exactly what the profile cannot. Our geometric kernel at p_local = 0.75
has a mean range ~1.3 nucleosomes — at the short end of the measured 1–4.

## 12. The KDM6B problem, and a resolution

The census makes the KDM6B arm untenable as currently used:
- kcat 0.70 /min is **peptide-only, catalytic domain only**. **No KDM6 enzyme has
  ever been assayed on a nucleosome.** Calibration from HDACs (10–200x slower on
  nucleosomes than free H3) suggests 0.70 is an upper bound by 1–2 orders.
- **KDM6B is undetectable** in the Itzhak proteome; KDM6A is 4,183 copies
  (EZH2:KDM6A ≈ 34:1). No protein-level induction fold exists anywhere.
- De Santa 2009 (grade A): **88% of Jmjd3 target genes had no H3K27me3 before
  stimulation**; of Jmjd3 peaks near an H3K27me3 peak, 16.3% lost >2-fold mark
  **against a 17.6% genome background** — no enrichment. Authors: activity
  "largely independent of H3K27me3 demethylation." Supported by enzyme-dead UTX
  knock-in rescuing mesoderm (Wang 2012) and UTX catalysis being dispensable for
  viability (Shpargel 2012).

**But removing it breaks the model.** At measured turnover (del_mk 0.004):

| | core | plateau | depletion |
|---|---|---|---|
| full model | 0.121 | 0.835 | **+0.714** |
| no KDM6B | 0.799 | 0.841 | +0.042 |
| no KDM6B, no p300 | 0.824 | 0.841 | +0.017 |
| exchange only | 0.953 | 0.840 | **−0.113** (inverted) |

At the measured low turnover the mark is so stable that *blocking writing* does
nothing — you need active removal. So three things cannot all be true.

**Partial resolution — replication carves it.** With division ON and no
demethylase at all:

| cycle | core | plateau | depletion | edge |
|---|---|---|---|---|
| no division | 0.799 | 0.841 | +0.042 | none |
| 48 h | 0.585 | 0.702 | +0.117 | 1160 bp |
| **24 h** | 0.398 | 0.541 | **+0.143** | **1540 bp** |
| 12 h | 0.209 | 0.251 | +0.042 | 1540 bp |

S phase halves the mark everywhere; the core fails to recover because PRC2
catalysis is blocked there, while the flanks are rewritten. **This needs no
demethylase** — which is what Alabert (BIC rejects demethylation terms) and
De Santa both require. It is non-monotonic in cycle length, peaking near 24 h.

It is also **5x weaker** than the KDM6B-driven boundary (0.143 vs 0.714), so it
is a mechanism, not yet a full explanation. Removing the K4me3 block inverts the
profile (−0.308), so K4me3 remains essential.

**Testable consequence: the depleted core should be weaker or absent in
non-dividing cells.** Our 0 h-division run gives depletion 0.042 vs 0.143 at a
24 h cycle.

## 13. In vitro ↔ in vivo cross-validation (not previously made in the literature)

| | fast / transient | slow / stable |
|---|---|---|
| Choi 2017 in vitro, 12-mer | sub-second | 2.28 ± 0.44 s |
| Choi + PHF1, 12-mer | sub-second | 6.48 ± 0.52 s |
| Tatavosian 2018 in vivo, mESC | 0.9 s | 10.0 ± 0.2 s |

Same biexponential structure, same order of magnitude, in vivo modestly slower —
consistent with accessory subunits (PHF1 alone 2.3 → 6.5 s) plus crowding.

**But k_on does not transfer.** Choi's k_on = 1.4–4.4 × 10⁷ M⁻¹s⁻¹ predicts a
67 ms search time at 3,000 targets; measured in vivo search is **210 s** —
a ~3,000x discrepancy. Parameterise search from the in vivo number, never from
solution kinetics.

## 14. Verified Kd values (all grade A, primary text)

| | Kd |
|---|---|
| 30 bp PRE DNA, PRC2 core | 64 ± 4 nM |
| + PHF1 | 30 ± 3 nM |
| 215 bp mononucleosome, core | 163 ± 13 nM |
| + PHF1 | 66 ± 4 nM (**P = 0.057, not significant**) |
| 2.5 kb DNA, AFM (Heenan 2020) | 150 ± 12 nM |
| DNA looping (Heenan) | 900 ± 400 nM |

PHF1 acts **entirely through k_off**; k_on is statistically identical across
complexes. Heenan also finds **no positional preference for an embedded CpG
island** at 2 nm accuracy — so CGI nucleation must come from accessory proteins
(PCL winged-helix), not the PRC2 core. Our `nucleation_gain` term should be
attributed to PCL, not to PRC2 itself.

---

# 15. LATE CORRECTION — abundance has a 3–5x error bar, not a value

The antagonist census returned a post-hoc correction after a background download
completed. **KDM6B is NOT undetectable** — 859 copies/cell (Bekker-Jensen 2017
Cell Syst, PMC5493283, Table S3, iBAQ calibrated against HeLa standards,
14,237 protein groups, grade A). It was missing from Itzhak only because 859
copies sits at the depth limit of an 8,710-protein experiment.

More importantly, the two grade-A datasets disagree systematically:

| protein | Itzhak 2016 | Bekker-Jensen 2017 | ratio |
|---|---|---|---|
| **EZH2** | 142,635 | **30,794** | **4.6x** |
| SUZ12 | 194,816 | 58,010 | 3.4x |
| EED | 351,948 | 88,382 | 4.0x |
| EP300 | 27,961 | 37,256 | 1.3x |
| CREBBP | 54,919 | 13,901 | 4.0x |
| KDM6A | 4,183 | 974 | 4.3x |
| KDM6B | not detected | 859 | — |
| HDAC1 | 461,170 | 850,903 | 1.8x |
| HDAC2 | 1,542,148 | 514,720 | 3.0x |

**Both datasets get PRC2 stoichiometry wrong** (EED > SUZ12 > EZH2 for a strictly
1:1:1 complex), which calibrates how far to trust either. **Any single-dataset
abundance parameter should carry ±3–5x, not the quoted precision.**
Section 2's table should be read with that band.

**KDM6A and KDM6B are equimolar** (974 vs 859). There is therefore **no abundance
basis for privileging KDM6B over KDM6A** in the model. Whole-cell ratios:
HDAC1 : HDAC2 : EP300 : EZH2 : CREBBP : KDM6A : KDM6B ≈ 874 : 529 : 38 : 32 : 14 : 1.0 : 0.9.
KDM6B is ~36x rarer than EZH2 and ~1,000x rarer than HDAC1; even a 10-fold
immediate-early induction leaves it ~100x below HDAC1.

## The maintenance constraint, recomputed across the abundance range

| abundance source | PRC2 | catalysing population | required kcat | vs measured max 0.474 /min |
|---|---|---|---|---|
| Itzhak (high) | 152,334 | stable only | 0.60 – 2.10 | exceeds all |
| Itzhak (high) | 152,334 | **stable + transient** | **0.12 – 0.42** | **within range** |
| Bekker-Jensen (low) | 32,794 | stable only | 2.79 – 9.75 | exceeds all |
| Bekker-Jensen (low) | 32,794 | stable + transient | 0.56 – 1.95 | exceeds all |

So the conclusion is **abundance-dependent**:
- At the high estimate, transient catalysis closes the gap and the maintenance
  picture is self-consistent.
- **At the low estimate it does not close.** Even with every bound molecule
  catalysing at the fastest measured nucleosomal rate, PRC2 falls 1.2–4x short.

If the low abundance estimate is right, one of these must give: fewer than
2–7 x 10^6 marked nucleosomes; an in-cell kcat above any in vitro value
(plausible — most published values are kcat,app at sub-saturating SAM); or
maintenance that is not PRC2-autonomous. **This is now the sharpest open
question the census produced, and it is answerable by measuring EZH2 copy number
properly — which, per section 2, has never been done by a gold-standard method.**
