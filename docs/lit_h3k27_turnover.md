# H3K27 methylation turnover and domain kinetics — primary-literature census

Compiled 2026-10-23. Every number below was read in the **Methods, Results, or a table of the primary
paper** (full text retrieved via EuropePMC `fullTextXML`, the built-in browser on PMC, or a local PDF).
Nothing here is taken from an abstract unless explicitly flagged `UNVERIFIED — abstract only`.

---

## 0. The central confound, stated up front

In a dividing cell an apparent "half-life" of a histone mark can contain three physically distinct
processes:

| Term | Process | What makes it go away |
|---|---|---|
| **(a) Enzymatic** | demethylation in situ (KDM6A/B, KDM7, etc.) | catalysis on a nucleosome that stays put |
| **(b) Exchange** | replication-independent nucleosome eviction/replacement (H3.3 deposition, Pol II passage, remodellers) | the marked histone leaves chromatin and is degraded |
| **(c) Dilution** | replication-coupled deposition of naive new H3.1/H3.2 | 50% of nucleosomes at each locus are replaced by unmarked ones every S phase |

**Every row of the table below carries an explicit "includes" column.** The short version of what the
literature actually supports:

- Methods that normalise *within the old-histone pool* (Zheng 2014) remove **(c)** and uniform **(b)**,
  leaving enzymatic removal plus *site-selective* exchange.
- Methods that follow a single replicated DNA segment with the writer inhibited (Reverón-Gómez 2018)
  isolate **(a)+(b)** on parental nucleosomes and find them **near zero over one cell cycle**.
- Methods that follow bulk mark loss after writer loss in vivo (Jadhav 2020) find the decay is
  **quantitatively accounted for by (c) alone**.
- Methods that measure *formation* kinetics (Zee 2010, Zee 2012, Alabert 2015, Oksuz 2018, Højfeldt
  2018) measure a composite of synthesis + deposition + de novo methylation **racing against ongoing
  dilution**, and almost none of them separate the two.

A recurring misquotation trap: **Zee et al. 2010 `t½` values are half-maximal times of *formation*,
not decay half-lives.** They are routinely cited as if they were decay. See row Z10 and the caveat.

Confidence grades: **A** = numeric value in a table/Results with explicit n and an explicit treatment of
the dilution confound. **B** = numeric value in Results, confound handled only by argument or partially.
**C** = read in the primary paper but qualitative, figure-only, or derived by me from read values.
**X** = abstract only / not verifiable.

---

## 1. Master table

| # | Quantity | Value (as printed) | Method | Cell type | Cycling? | How dilution was separated | **Half-life/rate includes** | n | Source (location) | Grade |
|---|---|---|---|---|---|---|---|---|---|---|
| **Z14-1** | H3.1 **K27me1** decay t½ on old histones | **52 ± 16 h** | Dual SILAC pulse–chase ([13C6]Arg + [13C1,D3]Met, 5 d label → light chase, 0/4/8/24/48 h), SRM-MS on HPLC-fractionated H3.1 | HeLa | **Proliferating**, asynchronous | All species normalised to **total old (heavy-Arg) histone signal**, which divides out bulk replicative dilution and uniform histone loss | **(a) enzymatic + (b) site-selective exchange.** Excludes (c) and excludes uniform (b) | 2 biol × 3 tech | Zheng 2014, PMC4184988, **Table 1A** | **A** |
| **Z14-2** | H3.1 **K27me2** decay t½ on old histones | **57 ± 20 h** | ″ | HeLa | Proliferating | ″ | ″ | 2 biol × 3 tech | Zheng 2014, Table 1A | **A** |
| **Z14-3** | H3.1 **K27me3** decay t½ on old histones | **158 ± 23 h** | ″ | HeLa | Proliferating | ″ | ″ | 2 biol × 3 tech | Zheng 2014, Table 1A | **A** |
| **Z14-4** | H3.1 **K27 overall methyl** (me1+me2+me3) t½ | **89 ± 37 h** | ″ | HeLa | Proliferating | ″ | ″ | 2 biol × 3 tech | Zheng 2014, Table 1A | **A** |
| **Z14-5** | Reference: **H3.1 protein** t½ (YQKSTELLIR) | **33.6 ± 0.5 h** — the paper's own cell-cycle length estimate | ″ | HeLa | Proliferating | This *is* the dilution rate | **(c) dilution + bulk (b)** | 2 biol | Zheng 2014, Table 1C | **A** |
| **Z14-6** | Comparators in same experiment | K4me3 **6.8 ± 0.6 h**; K4me2 24 ± 2.1; K4me1 418 ± 6.8; K9me1 23 ± 3.1; K9me2 28.5 ± 0.04; K9me3 **221 ± 105**; K36me1 44 ± 11; K36me2 61 ± 7.4; K36me3 **57 ± 3.2**; K79me1 110 ± 33; K79me2 176 ± 106 | ″ | HeLa | Proliferating | ″ | (a)+(b, selective) | 2 biol × 3 tech | Zheng 2014, Table 1A | **A** |
| **Z14-7** | H3.1 **K27 effective rate constants** (×10⁻² h⁻¹) | methylation: me0→me1 **3.58 ± 0.02**, me1→me2 **0.41 ± 0.08**, me2→me3 **0.54 ± 0.08**; demethylation: me1→me0 **1.4 ± 0.5**, me2→me1 **1.1 ± 0.4**, me3→me2 **0.4 ± 0.2** | ODE fit (M4K) to 10 SILAC species, ODE45/MATLAB | HeLa | Proliferating | Fit to old-histone-normalised relative levels | Demethylation terms are **(a)+(b, selective)** lumped | 2 biol | Zheng 2014, **Table 1B** | **A** |
| **Z14-8** | Derived: K27me3→me2 implied half-life | 0.4×10⁻² h⁻¹ ⇒ **≈173 h** | derived by me from Z14-7 | — | — | — | (a)+(b, selective) | — | derived | **C** |
| **Z10-1** | **K27me1/2/3 half-maximal time of FORMATION** | K27me1:1 **0.470 ± 0.005 d (11.3 h)**; K27me2:2 **1.145 ± 0.001 d (27.5 h)**; K27me3:3 **3.128 ± 0.032 d (75.1 h)** | Heavy-methyl SILAC (13CD3-Met), daily sampling for 7 d, LTQ-Orbitrap | HeLa (suspension, Joklik) | **Proliferating**, asynchronous | **Not separated.** Reference unmodified peptides give overall H3 t½ = **1.298 ± 0.007 d** (≈ the ~24 h doubling time), which is the dilution benchmark | **Formation**, i.e. new-histone deposition + de novo methylation, measured against a background of (c). *Not a decay half-life.* | 2 tech replicates | Zee 2010, PMC2823435, **Table 1** | **B** |
| **Z10-2** | K27 kinetic parameters (day⁻¹) | α 0.3848 ± 0.054; k0→1 **7.0139 ± 0.072**; k1→2 **4.1676 ± 0.007**; k2→3 **0.2329 ± 0.0014**; k1→0 0.0016 ± 0.0084; k2→1 **0.0635 ± 0.0022**; k3→2 **0.1680 ± 0.0101**; μ0 1.8727; μ1 0.0047; μ2 0.0000; μ3 0.0003 | ODE fit (MATLAB, 1000–10000 iterations, RMSD ≤125% of min) | HeLa | Proliferating | Table legend: "All four μ parameters **already account for dilution** and thus represent degradation **above** the dilution rate". ⚠️ Methods state μ were constrained ≥ ln(2) = 0.693/d, yet μ1–μ3 are reported far below that — internal inconsistency, treat μ with caution | k3→2 is **(a) enzymatic only** (by model construction); μ is **(b) above (c)** | 2 tech | Zee 2010, **Table 2** | **B** |
| **Z10-3** | Derived: K27me3 enzymatic-only half-life | k3→2 = 0.168 d⁻¹ ⇒ **≈4.1 d ≈ 99 h** | derived by me from Z10-2 | — | — | — | (a) only | — | derived | **C** |
| **Z12-1** | Timing of K27me3 re-formation across the cycle | **Begins in S phase, completed only in G1 of the next cell cycle.** De novo me3 from the unmodified state (me3:3, during S and G2/M) reaches **~50% of newly synthesised H3**; the remaining ~50% comes from conversion of **pre-existing me2** (me3:1) **during G1** | Heavy-methyl SILAC (13CD3-Met) in double-thymidine- and thymidine/nocodazole-synchronised cells; time-resolved MS | HeLa S3 | **Synchronised, cycling** | Synchronisation + explicit labelled/unlabelled methyl bookkeeping; dilution is the denominator, not removed | Formation only; **(c)** is the thing being out-run | 3 biol replicates | Zee 2012, PMC3434498, Results + Fig 5C/6C/7 | **A** |
| **Z12-2** | Excess dimethylation flux | New K27me2 events exceed newly synthesised H3 by **>25% before M phase** (1.25 × 40% = 50%) and by **~50% in the subsequent G1** (1.5 × 40% = 60%). Total K27me2 ≈ **2×** total K27me3 abundance | ″ | HeLa S3 | Cycling | ″ | Formation | 3 biol | Zee 2012, Results | **A** |
| **Z12-3** | K27 dimethylation mechanism | **Processive** (me2:2 exceeds me2:1 a few hours *before* me1:1 exceeds me1:0) | ″ | HeLa S3 | Cycling | ″ | — | 3 biol | Zee 2012, Fig 6A/B | **B** |
| **SC09-1** | New-histone K27 state after deposition | Transient **monomethyl** peak at/immediately after deposition; at **6 h** after release new H3 is mostly me0/me1 with **no appreciable K27me3** (me3 is the dominant form on old H3); patterns on old vs new still differ at **20 h**. Conclusion: "the histone code needs **one full cell cycle** in order to become fully re-established" | Pulsed SILAC (R0/R4/R10), MALDI-TOF on H3 peptide 27–40 | HeLa | Double-thymidine synchronised, cycling | Not separated; old (R0) vs new (R4) compared directly | Formation on new histones | 3 biol replicates | Scharf 2009, PMC2731903, Results/Fig 3D, Fig 5 | **B** (MALDI-TOF cannot resolve me3 from acetylation on this peptide, and does not resolve K27 from K36 on peptide 27–40 — the authors say so) |
| **AL15-1** | **Dilution factor per S phase** | Di- and trimethyl marks are **diluted exactly 2-fold immediately after replication** (t = 0 h) | Nascent chromatin capture (NCC; 15 min b-dUTP pulse) + double/triple pulsed SILAC, Orbitrap MS | Cell line not restated in Methods; NCC protocol is from Alabert 2014 (**HeLa S3**) | **Thymidine-synchronised, cycling** | NCC isolates only newly replicated DNA; parental-histone PTM level is the pre-replication reference | This *is* (c), measured directly | 3 ≤ n ≤ 9 | Alabert 2015, PMC4378191, Fig 3A | **A** (cell line provenance = B) |
| **AL15-2** | Restoration of total PTM level on replicated DNA | **Fully restored within one cell cycle (24 h)**; across marks the establishment kinetics on new histones span **2–24 h** | ″ | ″ | Cycling | ″ | Formation against (c) | 3 ≤ n ≤ 9 | Alabert 2015, Fig 3A–C | **A** |
| **AL15-3** | **H3K27me3 lag on new histones** | H3K9me3 and H3K27me3 **remain underrepresented on new histones at 24 h**; the pre-replication level is restored by continued methylation of **both new and old** histones, and addition of these trimethyl marks **continues for up to three cell generations** after histone incorporation (so me3 accumulates with histone age) | ″ | ″ | Cycling | ″ | Formation; the "3 generations" figure is the key number for mark–histone age coupling | 3 ≤ n ≤ 9 | Alabert 2015, Fig 3D/3E, Suppl S6B | **A** |
| **AL15-4** | Mitotic requirement | Nocodazole block **impedes** H3K9me3 and H3K27me3 acquisition on new histones (10, 16, 18.5 h after labelling; mean of three time points) | NCC ± nocodazole | ″ | **Arrested in mitosis** | — | Isolates the cell-cycle dependence of formation | n ≥ 3, p < 0.0001 / p < 0.05 | Alabert 2015, Fig 4A | **A** |
| **AL15-5** | **Proliferating vs G0** | Contact-inhibited G0 cells show a **significant increase** in H3K27me2/me3 (also H3K79me1/2, H4K20me2/3) vs exponentially growing. H3K9 methylation did **not** further accumulate in G0 | Histone MS on exponential vs contact-inhibited cells | **Primary human TIG-3 fibroblasts** | Both states compared | Direct proliferating-vs-arrested contrast | Demonstrates that steady-state K27me3 is **set by (c)**, not by (a) or (b) | n = 3, p < 0.01 | Alabert 2015, Fig 4B | **A** |
| **RG18-1** | Positional fidelity of recycling | H3K27me3 occupancy patterns reproduced on newly replicated DNA with no blurring; mean nascent-vs-parental peak displacement (measured for H3K4me3) **≈170 bp**, below the 250 bp fragment resolution | ChOR-seq (10–20 min EdU pulse → Click biotin → streptavidin → ChIP-seq), *Drosophila* spike-in | HeLa S3 | Thymidine-synchronised, cycling | Only EdU-labelled (once-replicated) DNA is analysed | — | not stated per-figure | Reverón-Gómez 2018, PMC6202308, Fig 2, S2D/E | **A** |
| **RG18-2** | **H3K27me3 restoration kinetics, locus-resolved** | Time points 0, 4 h (G2), 10 h (early G1), 24 h (arrested at G1/S to prevent re-replication). Signal accumulates at **all** time points, with the **major increase after mitosis**. **~80% of H3K27me3 regions are "R24"** (take the full 24 h); a substantial subset R10; **almost none restored before mitosis (R0/R4)**. Faster-restoring regions have higher H3K27me3 and higher EZH2 (ENCODE) occupancy | Spike-in-normalised qChOR-seq, 2 kb windows | HeLa S3 | Cycling, then G1/S-arrested at 24 h | Re-replication explicitly prevented by the G1/S arrest at the 24 h point | **Formation**, with (c) fixed at exactly one round | not stated per-figure | Reverón-Gómez 2018, Fig 4A–H | **A** |
| **RG18-3** | **Domain width is not restored — it is already correct** | Restoration increases H3K27me3 **within** domains without changing their width; domain borders faithfully demarcated at all time points | ″ | HeLa S3 | Cycling | ″ | — | — | Reverón-Gómez 2018, Fig 4C/4E | **A** |
| **RG18-4** | ⭐ **Active removal over one cell cycle ≈ 0** | With **EZH2 inhibitor** added before EdU labelling: the nascent (T0) H3K27me3 landscape is unchanged (so it is **entirely parental recycling**, with "little, if any, contribution from de novo methylation"); and at **T24** that landscape **persists**. **Only 3.5% of regions show a moderate reduction**; the large majority are unchanged. Regions of moderate loss do **not** correlate with H3.3 occupancy | qChOR-seq ± EZH2i | HeLa S3 | Cycling through one division | Follows a once-replicated DNA segment, writer off | This isolates **(a) + (b)** with (c) held to exactly one round — and finds them **negligible over ~24 h** | not stated per-figure | Reverón-Gómez 2018, Fig 5A–H | **A** |
| **RG18-5** | H3K4me3 comparator | **~50% of H3K4me3 regions restored within 1 h** of replication, the remaining ~50% within 6 h; complete by G2, no further gain at 12 h | qChOR-seq | HeLa S3 | Cycling | ″ | Formation | — | Reverón-Gómez 2018, Fig 3C–E | **A** |
| **JA20-1** | ⭐ **Dilution vs active removal, in vivo** | Eed⁻/⁻ Lgr5⁺ intestinal stem cells: H3K27me3 ChIP-seq (TSS ± 1 kb) at days **0, 5, 9, 12, 14, 17**. A stochastic crypt model that **halves H3K27me3 at every division and contains nothing else** fits the data best at a **mean ISC division interval of ~4 days**, which matches (or is slightly slower than) the reported 3–3.5 d ISC cycle. Authors: "this agreement … suggests that **other factors, such as enzymatic demethylation, have a minor or no role**" | Conditional KO in vivo + quantitative ChIP-seq + graph-based stochastic model (iGraph, ~2,000 trajectories) | **Mouse small-intestinal Lgr5⁺ ISCs, in vivo** | Proliferating in vivo | The model *is* the separation: a pure-dilution model exhausts the data | Decay here is **(c) only**; (a) and (b) are bounded to be small | ChIP-seq replicates at 6 time points | Jadhav 2020, PMC7376365, Fig 4C–E, S2E | **A** |
| **JA20-2** | Residual mark after writer loss | Ezh2⁻/⁻ villus cells: **25%–62% reduction** of all methylated H3K27 forms by immunoblot; **~40% uniform residual H3K27me3** across all PRC2 target sites, and that residual keeps targets inactive | Immunoblot quantitation + ChIP-seq | Mouse intestinal villus epithelium, in vivo | Proliferating | — | Steady-state level, not a rate | — | Jadhav 2020, Fig 1A/1B, S1C | **B** |
| **JA20-3** | **Arrest abolishes the loss** | SU-DHL-10 (EZH2-mutant lymphoma) + UNC1999: **geometric decline of H3K27me3 after successive cell divisions**, and **cell-cycle arrest attenuates this loss**. Palbociclib co-treatment **negates** target-gene activation. KARPAS-422 (cycles ~2.5× slower) shows correspondingly delayed activation. Most targets activate after **~4 cell divisions** (reproduced with the EED inhibitor MAK-683) | ChIP-qPCR + CellTrace division-resolved RNA-seq | Human SU-DHL-10, KARPAS-422 | Proliferating vs palbociclib-arrested | **Direct arrest control** | Confirms loss is **(c)** | biological replicates, n not printed per panel | Jadhav 2020, Fig 6A–C, S4C–E | **A** |
| **OK18-1** | **De novo domain formation time course** | EED-KO mESC + 4-OHT-inducible EED rescue. Time points are **0, 12, 24, 36 h only**. At **12 h**: H3K27me3 + EED + SUZ12 confined to **nucleation sites** (237 strong, 1,389 weak); H3K27me2 at nucleation sites with some already at spreading sites. At **36 h**: H3K27me3 present at both nucleation and **spreading sites** (36,618 total), and "the steady-state WT distribution of H3K27me3 was **approximated**"; H3K27me2 has moved off nucleation sites onto spreading sites | Cre-lox inducible rescue + spike-in-normalised ChIP-seq, immunofluorescence, immuno-FISH | **mESC** (E14-TG2a / C57BL/6) | **Proliferating** (mESC, ~12 h cycle) throughout | **Not separated at all.** De novo deposition is racing against ~3 rounds of dilution over the 36 h window; the paper does not address this | **Formation**, confounded by (c) | ChIP-seq, not replicated per time point in the text | Oksuz 2018, PMC7700016, Fig 2B–F | **B** |
| **OK18-2** | Global stoichiometry | H3K27me3 = **20% of H3K27** in WT mESC; **~1%** in EED aromatic-cage mutants (F97A, Y365A); **0%** in EED KO | Quantitative histone MS on acid-extracted histones | mESC | Proliferating | — | Steady state | — | Oksuz 2018, Fig 1B | **A** |
| **OK18-3** | Spatial scale of spreading | 73% of the ~200 cage-mutant-retained (nucleation) sites lie within **±5 kb of a TSS**; **71% of H3K27me3 spreads ±5 kb outside TSSs**; "proximal" spreading sites are defined as within **±10 kb** of a nucleation site, "distal" beyond that (including far-cis intra-chromosomal contacts) | ChIP-seq region classification | mESC | Proliferating | — | — | — | Oksuz 2018, Fig 1E/2C | **B** |
| **OK18-4** | Derived spreading velocity | Spreading sites are populated between 12 h and 36 h. If the proximal class (±10 kb) is filled in that 24 h window, the implied front velocity is **~0.4 kb/h ≈ 2 nucleosomes/h** — **an order-of-magnitude bound, not a measurement**. The paper reports **no kb/h figure and fits no time constant** | derived by me from OK18-1/OK18-3 | mESC | Proliferating | — | would include (c) | — | derived | **C** |
| **HJ18-1** | **De novo re-establishment after total erasure** | E14 mESC, EPZ6438 10 µM × 7 d (**>10 divisions; parental histone pool diluted >1000-fold**), then washout and 4 d of observation. **Both global H3K27me1/2/3 levels and genome-wide patterns are restored within 4 days.** H3K27me3 never appears transiently in regions where it is not normally found; me1 and me2 appear only as intermediates | Inhibitor washout + spike-in-normalised ChIP-seq + ChIP-qPCR, western | mESC E14 (2i/LIF) | **Proliferating** | Not separated — restoration occurs while cells keep dividing | Formation against (c) | Western ≥3×; **ChIP-qPCR kinetics n = 1 biological with technical duplicates** (one biological replicate in Suppl. Fig 4c) | Højfeldt 2018, PMC5842896, Fig 3a–c | **B** (the kinetic panel is n=1) |
| **HJ18-2** | Rank order of regional rates | "The methylation rate in different genomic regions **directly correlates with the steady-state methylation level** of the region: H3K27me3-positive regions methylate faster than H3K27me2 regions, and **H3K27me1-dominant regions are slowest**." Global me1 recovers first, but that early me1 sits in regions destined to become me2/me3 | ChIP-seq heatmaps + ChIP-qPCR (units in figure: **% of input / hr washout**) | mESC E14 | Proliferating | — | Formation | n = 1 biol + replicate | Højfeldt 2018, Fig 3b/3c | **C** (values are figure-only; no table) |
| **HJ18-3** | Self-propagation is not required | Ezh1/Ezh2 dKO mESC (all H3K27me patterns erased, no trace detectable) + EZH2 re-expression → **accurate restoration of the whole genome-wide pattern**; SUZ12 binding to CGIs is independent of H3K27 methylation and of EED | Genetic KO + rescue, spike-in ChIP-seq | mESC | Proliferating | — | — | ≥3× for westerns | Højfeldt 2018, Fig 2a–c, Fig 4a | **A** |
| **AL20-1** | Model structure (rate constants) | Domain model (8 or 15 domains, **38 parameters: 23 methylation rate constants k_ij→kl, the dilution/cell-cycle rate c, and domain weights w — and NO demethylation terms**) beats the global model (which does include demethylation) by BIC. A model with demethylation **cannot** explain the stability of K27me3 under EZH2 inhibition | Triple-SILAC generational MS + ODE model selection (AMICI/PESTO, MATLAB) | mESC E14 (2i/LIF) | Proliferating | `c` is an explicit fitted dilution parameter in the ODE | Explicitly **rejects global (a)**; stability attributed to domain structure | 6 ≤ n ≤ 9 (Fig 1); n = 3 (Figs 2D/E, 3, 4A–C); n = 2 (Fig 4) | Alabert 2020, Cell Rep 30:1223, Results + STAR Methods | **A** |
| **AL20-2** | Steady-state fluxes (the only absolute rates printed) | Fluxes are in **h⁻¹** (flux = state abundance × rate constant). The K27me0→me3 path at K36me0 carries **57.4% of all flux**. Worked examples: **2.5 h⁻¹** (further K27 methylation of K27me1K36me0) vs **0.21 h⁻¹** (K36 methylation of the same state) = 11.9-fold; **0.67 h⁻¹** vs **0.15 h⁻¹** = 4.5-fold. The antagonism rule is stated for "all cases with fluxes over **0.5 h⁻¹**" | ″ | mESC E14 | Proliferating | ″ | Forward methylation flux only | n = 3 | Alabert 2020, Fig 2F + STAR Methods | **A** |
| **AL20-3** | **Parental K27me3 feeds forward onto de novo rate** | After 7 d EZH2i (K27me3 drops **below 3%**, me2 strongly reduced, me1 accumulates via EZH1), establishment of K27me2/3 on new (generation-2) histones is **delayed**. Model selection: the rate constant **k(K27me0K36me0 → K27me1K36me0) is roughly halved (×0.47)** in cells lacking K27me3; k(K27me2K36me0 → K27me3K36me0) is also reduced; K36 mono-methylation on K27me2K36me0 **increases ~10-fold** | ″ | mESC E14 | Proliferating | ″ | De novo methylation rate | n = 3 | Alabert 2020, Fig 3B–F | **A** |
| **AL20-4** | Cycle length sets steady-state K27me3 | Model prediction: **shortening the cell cycle, all other parameters unchanged, reduces the steady-state level of K27me3**. Cited as consistent with mESCs having higher K27me3 in 2i than in serum/LIF (faster cycle) | ODE simulation | mESC (in silico) | — | — | Pure (c) effect on steady state | — | Alabert 2020, Fig S2I | **B** (prediction, not measurement) |
| **DE10-1** | ⭐ **Nucleosome turnover genome-wide (CATCH-IT)** | **Mean nucleosome lifetimes, one full turnover cycle: ~1 h at the peak just downstream of the TSS of active genes; ~1 h at GAF (trxG) sites; ~1.5 h at EZ+PSC (PcG) sites** (Table S1). Authors call these **conservative** (upper bounds): any delay in Aha incorporation inflates them. "Nucleosomes within active genes and at epigenetic regulatory elements turn over multiple times during each **~15-h cell cycle**" | CATCH-IT: azidohomoalanine metabolic labelling → click-biotin → streptavidin → urea/salt wash to (H3/H4)₂ → tiling microarray. Kinetics from the **20 min** datasets | ***Drosophila* S2 cells** | **Proliferating**, unsynchronised, late log | Argued, not measured: replication through a nucleosome takes 1–2 s ≈ 1/50,000 of the cell cycle, so in an unsynchronised population the replication contribution is genome-wide homogeneous and only compresses dynamic range | **(b) replication-independent exchange only** | 2 independent experiments for the 20′ rate calculation | Deal 2010, PMC2879085, Results + Table S1 | **A** (for flies) |
| **DE10-2** | Interpretation the authors draw | "The fact that EZ is responsible for di- and tri-methylation of H3K27, but the nucleosomes that it modifies **turn over faster than a cell cycle**, argues against proposals that histone modifications required for cellular memory themselves transmit epigenetic information" | ″ | *Drosophila* S2 | Proliferating | ″ | — | — | Deal 2010, Discussion | **A** — note this **directly contradicts** RG18-4 in human cells; see §8 |
| **KR13-1** | **H3.3 exchange: promoters vs gene bodies vs heterochromatin** | Three categories. **Rapid**: sharp peaks appear **within 2 h**; **5′UTR and promoters have the highest turnover index, enhancers second**. **Intermediate**: broad peaks in transcribed gene bodies, apparent at **12 h**, maximal at **72 h**. **Slow**: telomeres, detectable only at **24–48 h**. 3′UTR/TES lowest. Turnover is **negatively correlated with H3K27me3 at regulatory regions** and with H3K36me3 at gene bodies | TET-ON HA/FLAG-H3.3, DOX time course 0,1,2,3,4,5,6,12,18,24,48 h (to 72 h), HA ChIP-seq, turnover index by linear regression | **Mouse NIH/3T3 MEFs** | **Replication-blocked**: aphidicolin 1 µg/ml × 18 h before DOX | **Aphidicolin arrest removes (c) by design** | **(b) replication-independent exchange only** | 2 independent experiments (TI reproducibility) | Kraushaar 2013, PMC3983652, Figs 3–4 | **B** — TI is a relative, unitless slope, **not** an absolute half-life |
| **DT16-1** | **H3.3 turnover, mammalian, locus-resolved** | Turnover index from 0/3/6/12 h. **Highest at active enhancers, with super-enhancers highest of all**; high at active TSSs (particularly **upstream** of the TSS) and at the 3′ ends of expressed genes; **low at silent TSSs**. **Polycomb-bound regions carry H3.3 but show LOW turnover** compared with active regions. Heterochromatin has no H3.3 at all. H3.1 also shows high turnover at these regions. ESC cycle **12 h**, NSC cycle **18 h** | time-ChIP: SNAP-tagged H3.3 (DOX ×4 doublings), 40 min CP-biotin pulse, chase, streptavidin pulldown + sequencing, spike-in controlled | **Mouse ESC → NSC** | **Proliferating** | Not separated; the TI is a relative loss slope over 12 h, i.e. about one ESC cycle | **(b) + (c) mixed**, expressed as a relative index | 3 ESC time-ChIP replicates | Deaton 2016, PMC4965263, Figs 2–4 | **B** — relative index only; ~1 h labelling means the very fastest turnover is **underestimated** |
| **HA12-1** | ⭐ **The only published spreading velocity — but it is H3K9me3, not H3K27me3** | "H3K9me3 propagated **symmetrically and continuously at rates of ~0.18 nucleosomes/hr** to produce domains of **up to 10 kb**." (≈ **0.036 kb/h** at a 200 bp repeat length — my conversion.) csHP1α recruitment significant by 6 h and near-saturated by 24 h; H3K9me3 first detectable at **18 h**; complete reporter silencing by **5 days** | Chromatin in vivo Assay (CiA): rapamycin-induced tethering of csHP1α to a modified *Oct4* allele; ChIP-qPCR along the locus + flow cytometry; kinetic lattice model (nucleation k₊ at target, propagation k₊ to adjacent marked sites, turnover k₋ everywhere; κ = k₊/k₋ sets the boundary) | **Mouse ESC and MEFs** | Both; **also tested serum-starved, non-dividing MEFs** | **Serum-starved non-dividing CiA MEFs still silenced GFP and built the full heterochromatic domain within 5 days** → replication-dependent histone exchange is **not required** for domain formation | Propagation rate is (a)-inverse; the fitted **k₋ lumps demethylation + nucleosome turnover** and "approximates previously measured rates of global histone turnover in HeLa S3 (Zee et al., 2010)" | — | Hathaway 2012, PMC3422694, Abstract(=Results claim), Figs 2, 4, 6, S5 | **B** — wrong mark for your purpose, and the rate is model-derived |
| **HA12-2** | Generality of the model | The κ framework fits the spatial profiles of **99.2% of non-centromeric H3K9me3 domains** in the mouse genome | ″ | mouse | — | — | — | — | Hathaway 2012, Results | **B** |
| **VI19-1** | Supporting (same CiA system) | Silencing velocity is **unaffected by promoter CpG content**; after short-term csHP1α recruitment both CpG-full and CpG-depleted reporters lose the domain and re-express GFP within **4–6 days**; H3K9me3 reduced but not ablated 4 d after washout | CiA + flow cytometry + ChIP | Mouse embryo fibroblasts | Proliferating | — | Domain decay includes (a)+(b)+(c) | — | Vignaux 2019, PMC6608945 | **C** |
| **LE18-1** | In vitro enzymology underpinning the feed-forward loop | PRC2-EZH2 activity is increased **~12-fold** by an H3K27me3 peptide; PRC2-EZH1 only **~2-fold**. PRC2-EZH2 activity on 40 bp-linker dinucleosomes is **5–6× higher** than on 30 bp-linker dinucleosomes; activity rises with linker length to 40 bp then declines | ³H-SAM HMT assays, 60 min, recombinant human PRC2 | recombinant / mESC validation | in vitro | n/a | n/a | n = 3 per data point, mean ± SD | Lee & Holder 2018, PMC5949877, Figs 1D/1E, 2A | **A** (in vitro) |

---

## 2. Answers to the seven questions

### Q1 — Half-lives of H3K27me1, me2 and me3 separately

The only paper that reports all three separately, with dilution explicitly removed, is **Zheng 2014**
(rows Z14-1→4). In proliferating HeLa, with the old-histone pool as the denominator:

| | t½ (h) | in cell cycles (33.6 h) | includes |
|---|---|---|---|
| K27me1 | 52 ± 16 | 1.5 | enzymatic + selective exchange |
| K27me2 | 57 ± 20 | 1.7 | enzymatic + selective exchange |
| K27me3 | **158 ± 23** | **4.7** | enzymatic + selective exchange |
| K27 overall | 89 ± 37 | 2.6 | enzymatic + selective exchange |

**Correction to your abstract-level note.** The abstract's "half-lives of multiple cell divisions"
applies to **K27me3 and K9me3 only**. K27me1 (52 h) and K27me2 (57 h) are *not* multi-division marks —
they turn over on roughly 1.5 cell cycles, i.e. faster than H3K36me3 (57 h) is slower, and comparable
to it. Your two quoted comparators are exact: **K4me3 = 6.8 ± 0.6 h, K36me3 = 57 ± 3.2 h.**

Also note the paper has a **units typo in Results §3.3**: it writes "the rate constants of K4me3, K9me1
and K9me2 demethylation are high (8, 5 and 2 h⁻¹)". Table 1B makes clear these are **×10⁻² h⁻¹**.
Do not propagate the text version.

Second, independent estimate of the *enzymatic-only* K27me3 term: **Zee 2010** fits
k(me3→me2) = **0.168 day⁻¹**, i.e. a demethylation-only half-life of **≈4.1 d ≈ 99 h** (row Z10-3).
Zheng's lumped me3→me2 constant implies ≈173 h (row Z14-8). These bracket **~100–175 h** for whatever
removes K27me3 from a nucleosome that is not being diluted. Both are longer than one cell cycle in
every mammalian line measured.

**Do not use Zee 2010 Table 1 as decay half-lives.** K27me1 0.47 d, K27me2 1.15 d, K27me3 3.13 d are
**half-maximal times of formation** of the fully heavy-labelled species. Against an overall-H3 formation
t½ of 1.30 d (≈ the doubling time), the interpretable statement is: **K27me3 formation is ~2.4× slower
than bulk histone replacement**, K27me2 is at parity, and K27me1 is ~2.8× faster.

### Q2 — De novo H3K27me3 on new histones, and the lag

This is the best-constrained question in the whole census, and four independent methods agree.

1. **New histones arrive essentially unmethylated at K27.** Transient K27me1 appears at or immediately
   after deposition; at 6 h post-release new H3 is still me0/me1 with no appreciable me3 (Scharf 2009,
   SC09-1).
2. **The mark is diluted exactly 2-fold at the fork** (Alabert 2015, AL15-1).
3. **Positional information is nonetheless perfect** — parental H3K27me3 is recycled onto the daughter
   strand in place, with ≤170 bp displacement, and this happens with EZH2 inhibited, so it is pure
   recycling (Reverón-Gómez 2018, RG18-1/RG18-4).
4. **Level restoration is slow and post-mitotic.** Almost nothing is restored before mitosis (R0/R4
   ≈ none); **~80% of H3K27me3 regions take the full 24 h** (RG18-2). Bulk MS agrees: me3 formation
   begins in S but completes only in the next G1, half from de novo methylation of new histones and
   half from conversion of pre-existing me2 (Zee 2012, Z12-1). Mitotic passage is required — nocodazole
   blocks acquisition (AL15-4).
5. **The lag runs past one cell cycle.** H3K27me3 (and H3K9me3) are still underrepresented on new
   histones at 24 h, and addition **continues for up to three cell generations** after a histone is
   incorporated, so K27me3 accumulates with histone age (Alabert 2015, AL15-3).

**Locus dependence:** high-PRC2/high-EZH2 regions restore faster; domain *borders* restore slowest
(RG18-2). Højfeldt agrees from the opposite direction: after total erasure, me3 regions re-methylate
fastest and me1-dominant regions slowest (HJ18-2).

**From zero:** after 7 d of EZH2i (>10 divisions, parental histones diluted >1000-fold), global levels
*and* the genome-wide pattern are rebuilt in **4 days** in mESC (HJ18-1) — and the rebuild is accurate,
with no transient ectopic me3.

**Feed-forward term:** parental K27me3 is not required for accuracy but it is required for speed —
without it, the K27me0→me1 rate constant is **halved** and the me2→me3 constant is also reduced
(Alabert 2020, AL20-3). This is the in vivo counterpart of the ~12-fold allosteric activation of
PRC2-EZH2 by an H3K27me3 peptide in vitro (LE18-1).

### Q3 — Spreading velocity of an H3K27me3 domain, in kb/h

**No primary paper reports a kb/h spreading velocity for H3K27me3.** This is the clearest gap.

- **Oksuz 2018 does not provide one.** Its ChIP-seq time course has only two informative points, 12 h
  and 36 h, and 36 h is simply the last time point, not a fitted plateau. The paper's claim is the
  qualitative one that at 36 h the WT distribution is "approximated". No time constant is fitted
  anywhere in the paper, and the word "kb/h" does not appear. (This corrects your abstract-level note:
  "reaches steady state ~36 h" should be read as "**has approached the WT distribution by the last time
  point, 36 h**", with no evidence about what happens after 36 h.)
- A crude bound can be derived from Oksuz's own region definitions: proximal spreading sites are
  ±10 kb from a nucleation site and are populated between 12 h and 36 h, giving **~0.4 kb/h ≈ 2
  nucleosomes/h** (row OK18-4). This is my derivation, grade C, and it is confounded: mESCs divide
  roughly every 12 h, so ~3 rounds of dilution occur inside that window and the paper does not
  correct for them.
- **The only measured propagation rate in the literature is for the wrong mark**: Hathaway 2012 reports
  **~0.18 nucleosomes/h (≈0.036 kb/h) for H3K9me3**, producing domains up to 10 kb, in the CiA:Oct4
  system (HA12-1). That is **~10× slower** than the Oksuz-derived K27 bound. Two useful features of
  that study for your purposes: (i) serum-starved, non-dividing MEFs still built the full domain in
  5 days, so replication is not required for spreading; (ii) the fitted turnover term k₋ was found to
  approximate global histone turnover rates, i.e. domain boundaries are set by a race between
  propagation and nucleosome turnover.

### Q4 — Fraction of per-S-phase H3K27me3 loss attributable to dilution vs active removal

**Essentially all of it is dilution.** Three independent lines:

- **Direct, one cell cycle, human cells.** With EZH2 inhibited, the parental H3K27me3 landscape laid
  down at the fork **persists unchanged 24 h later**, across mitosis. Only **3.5%** of regions show even
  a moderate reduction, and those do not correlate with H3.3 occupancy. Reverón-Gómez's own words:
  "domains laid down at the time of replication are remarkably stable … without challenges by
  demethylation or histone exchange" (RG18-4).
- **In vivo, many cell cycles.** H3K27me3 decay in Eed⁻/⁻ mouse Lgr5⁺ ISCs over 17 days is fully
  accounted for by a model containing **only** 50%-per-division halving, at a division interval
  (~4 d) matching the independently measured ISC cycle. The authors conclude enzymatic demethylation
  has "a minor or no role" (JA20-1).
- **Arrest control.** Blocking the cycle with palbociclib **attenuates** H3K27me3 loss and **prevents**
  target-gene derepression in EZH2-inhibited lymphoma cells (JA20-3).

So the working number for your model is: **per S phase, loss ≈ 50%, of which ≈100% is dilution and
≤~3.5% of loci show any additional removal.** The ~100–175 h enzymatic/selective-exchange half-life
from Q1 is consistent: over a 24–34 h cycle that term removes only ~12–20% of the mark, and most of
that is at the minority of loci with low PRC2 occupancy (domain borders).

**One caveat on the Oksuz/Alabert-2020 side.** Alabert 2020's model selection formally **rejects**
global K27 demethylation (AL20-1) — a demethylation-containing model cannot reproduce the observed
stability of K27me3 under EZH2i. That is the same conclusion from an entirely different data type.

**The dissenting datum is Drosophila.** Deal 2010 measures a ~1.5 h mean nucleosome lifetime at
PcG (EZ+PSC) sites in S2 cells — i.e. ten turnovers per 15 h cell cycle — and explicitly argues this
"challenges models for epigenetic inheritance that rely on stability of histone marks" (DE10-2). See §8.

### Q5 — H3.3 / H3.1 exchange at promoters and gene bodies

No mammalian study reports absolute exchange half-lives per genomic compartment; all three mammalian
measurements give **relative turnover indices**. The consistent picture:

| Compartment | Kraushaar 2013 (MEF, aphidicolin-arrested, DOX H3.3) | Deaton 2016 (mESC, SNAP H3.3 time-ChIP) |
|---|---|---|
| Promoters / 5′UTR / TSS-upstream | **Highest** TI; sharp peaks within **2 h** | High; highest **upstream** of active TSS |
| Enhancers | Second highest | **Highest overall; super-enhancers highest of all** |
| Gene bodies (active) | Intermediate; apparent at **12 h**, maximal at **72 h** | Elevated; high also at 3′ ends of expressed genes |
| Silent promoters | — | Low |
| **Polycomb-bound regions** | Turnover **negatively correlated with H3K27me3** at regulatory regions | H3.3 present but turnover **low** |
| Heterochromatin / telomeres | Slowest; telomeres detectable only at **24–48 h** | No H3.3 at all |

Kraushaar is the cleaner measurement of the pure transcription-coupled exchange term because the cells
were **aphidicolin-arrested**, so replication is excluded by design. Deaton's 0–12 h window in a 12 h
ESC cycle mixes exchange with dilution. Deaton also measured H3.1 and found it **also** turns over fast
at enhancers and super-enhancers in ESC, which is worth noting because it means "H3.1 = replicative
only" is not safe in stem cells.

The only absolute numbers in this space are Deal 2010's, and they are from flies (Q7).

### Q6 — How long a single nucleosome retains H3K27me3

Nobody has measured this on single nucleosomes. The two closest bounds:

- **Upper-bound-ish, human:** on a once-replicated DNA segment with EZH2 off, H3K27me3 occupancy is
  unchanged after 24 h at >96% of loci (RG18-4). So the retention time of an individual marked
  nucleosome is **≫ 24 h**.
- **Population half-life:** 158 ± 23 h for the me3 species on old histones (Z14-3), and ~100–175 h for
  the enzymatic/selective-exchange component (Q1). This is the best available proxy: a nucleosome that
  is never diluted retains K27me3 with a **~4–7 day half-life**.
- **Marked histones keep gaining methylation with age** — K27me3 addition continues for up to three
  cell generations after a histone is deposited (AL15-3), so the "retention" and "acquisition"
  timescales overlap heavily and a single-nucleosome retention time is not cleanly separable from them.
- **Counter-bound, Drosophila:** ~1.5 h mean nucleosome lifetime at PcG sites (DE10-1).

### Q7 — Nucleosome turnover genome-wide (metabolic labelling)

The canonical numbers, all from **Deal, Henikoff & Henikoff 2010, *Drosophila* S2 cells, ~15 h cell
cycle**, CATCH-IT with Aha, kinetics from 20-minute labelling, two independent experiments:

| Site | Mean nucleosome lifetime (one full turnover cycle) |
|---|---|
| Peak just downstream of the TSS of active genes | **~1 h** |
| GAF (trithorax-group) binding sites | **~1 h** |
| **EZ + PSC (Polycomb-group) binding sites** | **~1.5 h** |

The authors call these **conservative** (any delay in Aha incorporation inflates the estimate), and
state that nucleosomes in active genes and at regulatory elements "turn over multiple times during each
~15-hr cell cycle." Turnover across gene bodies scales strongly with expression quintile. These are
**replication-independent exchange only**: the authors argue the replication contribution is negligible
in an unsynchronised population because passage of the fork through a nucleosome takes 1–2 s,
~1/50,000 of the cell cycle.

Deal also replicates the yeast result they cite: <30 min mean lifetime at active yeast promoters and
boundary elements; by contrast yeast gene *bodies* turn over slowly, which they attribute to a genuine
biological difference from flies.

There is **no mammalian CATCH-IT-equivalent with absolute lifetimes** — see §9.

---

## 3. Parameter set I would hand to a model

| Parameter | Value | Source | Caveat |
|---|---|---|---|
| Dilution per S phase | **exactly 2-fold** (50% naive nucleosomes) | AL15-1, JA20-1 | measured and modelled independently |
| Parental recycling fidelity | positional, ≤170 bp displacement, EZH2-independent | RG18-1, RG18-4 | human, one cycle |
| k(me3 → me2), enzymatic | **0.004 h⁻¹** (Zheng) to **0.007 h⁻¹** (Zee) | Z14-7, Z10-2 | HeLa; lumped with selective exchange in Zheng |
| K27me3 species half-life, dilution removed | **158 ± 23 h** | Z14-3 | HeLa, proliferating |
| K27me1 / me2 half-lives, dilution removed | **52 ± 16 h / 57 ± 20 h** | Z14-1/2 | HeLa |
| k(me0→me1) : k(me1→me2) : k(me2→me3) | **3.58 : 0.41 : 0.54** ×10⁻² h⁻¹ (Zheng, HeLa) | Z14-7 | note Zee's HeLa ratio differs: **7.01 : 4.17 : 0.23** d⁻¹ — me1→me2 is 10× faster in Zee |
| Fraction of new me3 from pre-existing me2 vs from me0 | **~50 : 50** | Z12-1 | HeLa, synchronised |
| Lag before new histones match parental K27me3 | **> 24 h; continues ~3 cell generations** | AL15-3 | HeLa |
| Fraction of loci fully restored only by 24 h | **~80%** | RG18-2 | HeLa |
| Feed-forward gain from parental me3 on k(me0→me1) | **×2.1** (rate halves without it) | AL20-3 | mESC |
| Nucleosome exchange at PcG sites | **~1.5 h** mean lifetime | DE10-1 | ***Drosophila*** — do not transplant to mammals without flagging |
| Nucleosome exchange at active TSS | **~1 h** mean lifetime | DE10-1 | ***Drosophila*** |
| Spreading velocity | **~0.036 kb/h** (H3K9me3, measured) to **~0.4 kb/h** (H3K27me3, derived upper bound) | HA12-1, OK18-4 | 10× spread; neither is a clean K27 measurement |
| Steady-state K27me3 stoichiometry | **20% of all H3K27** in WT mESC | OK18-2 | quantitative MS |
| Effect of faster cycling | lowers steady-state K27me3 | AL20-4, AL15-5, JA20-3 | predicted + observed in three systems |

---

## 4. Corrections to the three priors you supplied

**Zheng 2014 (PMC4184988) — verified, with one correction and one caveat.**
K4me3 = 6.8 h and K36me3 = 57 h are exact. K27me3 = **158 ± 23 h** (4.7 cell cycles), K9me3 =
**221 ± 105 h** — so "multiple cell divisions" is right for the trimethyl marks. But the abstract's
framing does **not** extend to K27me1 (52 h) and K27me2 (57 h), which are ~1.5 cell cycles. Caveat: the
Results text misprints the demethylation rate constants as h⁻¹ when Table 1B gives ×10⁻² h⁻¹.

**Oksuz 2018 (PMC7700016) — partially corrected.**
Nucleation → spreading is confirmed: at 12 h the mark is confined to 237 strong + 1,389 weak nucleation
sites (CGI-dense, GA-rich and GCN-tandem-repeat motifs, 78% of strong sites carry one); by 36 h it has
reached 36,618 spreading sites, with H3K27me2 depleted from nucleation sites and enriched at spreading
sites. **But "reaches steady state ~36 h" overstates the data.** 36 h is the final time point; the
paper says the WT distribution was "approximated", fits no time constant, reports no velocity, and runs
the whole experiment in cells dividing roughly every 12 h with no correction for the ~3 rounds of
dilution that occur during the window.

**Alabert 2020 (Cell Rep 30:1223) — rate constants extracted as requested.**
The paper prints **no table of rate constants**; the 23 fitted k's live in the Zenodo code
(zenodo.org/record/3353481), not in the text. What is printed, in h⁻¹, are **fluxes** (abundance × rate
constant): 2.5, 0.67, 0.21, 0.15 h⁻¹ for the worked comparisons, with "all cases with fluxes over
0.5 h⁻¹" as the scope of the antagonism rule, and 57.4% of total flux on the K27me0→me3 / K36me0 path.
Two further model results you should have: (i) **the winning domain model contains no demethylation
terms at all** — a demethylation model is rejected by BIC because it cannot reproduce K27me3 stability
under EZH2i; (ii) **lack of parental K27me3 halves the k(K27me0K36me0→K27me1K36me0) rate constant
(×0.47)** and also reduces k(me2→me3), the quantitative in vivo measure of the allosteric feed-forward
loop. The fitted parameter vector also contains an explicit dilution/cell-cycle rate `c`.

---

## 5. The one genuine contradiction in the field

**Deal 2010 (flies) vs Reverón-Gómez 2018 (human).**

- Deal: PcG-bound nucleosomes in *Drosophila* S2 cells have a **~1.5 h mean lifetime** — ten exchanges
  per cell cycle — and therefore "nucleosome turnover … erases histone modifications", which "argues
  against proposals that histone modifications required for cellular memory themselves transmit
  epigenetic information."
- Reverón-Gómez: in HeLa with EZH2 off, **>96% of H3K27me3 regions are unchanged after a full cell
  cycle**, and the regions that do lose signal are **not** the H3.3-rich ones. They address Deal
  directly: "In *Drosophila*, H3K27me3-marked regions, such as Polycomb group response elements, have
  been identified as fragile high-turnover sites (Deal et al., 2010). However, in our inhibitor-treated
  HeLa cells, regions of moderate H3K27me3 loss did not correlate with higher occupancy of … H3.3."
- Deaton 2016 lands with Reverón-Gómez for mammals: PcG regions in mESC carry H3.3 but show **low**
  turnover; Kraushaar 2013 likewise finds turnover **anti-correlated with H3K27me3** at regulatory
  regions in MEFs.

**Resolution for modelling purposes:** treat the ~1.5 h PcG nucleosome lifetime as a *Drosophila*-
specific parameter. In mammals the transcription-coupled exchange term at Polycomb domains is small
enough that one cell cycle of EZH2 inhibition does not visibly erode the landscape.

---

## 6. What is still missing

1. **A kb/h spreading velocity for H3K27me3.** Nothing in the literature. Oksuz has two usable time
   points in a system with three cell divisions inside the window. The only measured propagation rate
   (0.18 nucleosomes/h) is for H3K9me3/HP1α at a single engineered locus. **This is the single biggest
   gap for a lattice/domain model**, and the two available anchors differ by ~10×.
2. **Separate turnover numbers for H3K27me1 and me2 in a non-cancer, non-HeLa line.** Every me1/me2
   half-life in the census is HeLa. Given that HeLa's H3K27me2 stoichiometry and cycle length are both
   atypical, these numbers are weakly transferable.
3. **Absolute mammalian nucleosome lifetimes.** All three mammalian H3.3 studies (Kraushaar, Deaton,
   Ha) report **relative, unitless turnover indices**. There is no mammalian CATCH-IT giving "mean
   lifetime = X h at promoters". You cannot currently write a mammalian exchange rate into a model
   from a primary measurement — only a fly one.
4. **H3K27me1 and me2 at single-locus resolution during the cell cycle.** ChOR-seq has only been run
   for me3 (and K4me3/K36me3/K79me3). Højfeldt has me1/me2/me3 kinetics but as ChIP-qPCR at seven loci
   with **n = 1 biological replicate**, and the values are figure-only.
5. **A direct measurement of K27me3 turnover in an arrested mammalian cell.** Alabert 2015's G0
   comparison (TIG-3) is a steady-state level comparison, not a decay rate. Hathaway's serum-starved
   MEF experiment is the closest to a non-dividing kinetic measurement and it is for H3K9me3. Nobody
   has run a heavy-methyl SILAC chase in a cleanly G0-arrested mammalian cell for K27.
6. **Single-nucleosome retention of H3K27me3.** No direct measurement; the best bound is ">24 h at
   >96% of loci" plus a 158 h population half-life.
7. **Nucleation site number and strength in a non-ESC mammalian cell.** Oksuz's 237 strong / 1,389 weak
   / 36,618 spreading partition is mESC-specific; nothing equivalent exists for a differentiated or
   tumour cell.
8. **Error bars on the Oksuz time course.** ChIP-seq at 12 and 36 h appears unreplicated in the text;
   the quantitative heatmap scales (0–2) are the only quantitation.
9. **Alabert 2020's individual rate constants** are only in the Zenodo deposit
   (zenodo.org/record/3353481, MATLAB/AMICI/PESTO). If you want the 23 k's as numbers, they have to be
   pulled from that code — the paper prints fluxes only.
10. **Replicate depth on the de novo establishment kinetics.** The two cleanest "from zero" experiments
    (Højfeldt washout kinetics, Oksuz induction) are n = 1 and n unstated respectively.

---

## 7. Sources, with access route

| Short | Full | PMCID | Access used |
|---|---|---|---|
| Zheng 2014 | Zheng, Tipton, Thomas, Kelleher, Sweet. *Proteomics* 14:2190–2199 | PMC4184988 | built-in browser (non-OA on PMC) |
| Zee 2010 | Zee, Levin, Xu, LeRoy, Wingreen, Garcia. *J Biol Chem* 285:3341 | PMC2823435 | built-in browser |
| Zee 2012 | Zee, Britton, Wolle, Haberman, Garcia. *Mol Cell Biol* 32:2503 | PMC3434498 | built-in browser |
| Scharf 2009 | Scharf, Barth, Imhof. *Nucleic Acids Res* 37:5032 | PMC2731903 | EuropePMC fullTextXML |
| Alabert 2015 | Alabert, Barth, Reverón-Gómez, Sidoli, … Imhof, Groth. *Genes Dev* 29:585 | PMC4378191 | EuropePMC fullTextXML |
| Reverón-Gómez 2018 | Reverón-Gómez, González-Aguilera, Stewart-Morgan, … Groth. *Mol Cell* 72:239 | PMC6202308 | EuropePMC fullTextXML |
| Oksuz 2018 | Oksuz, Narendra, Lee, Descostes, … Reinberg. *Mol Cell* 70:1149–1162 | PMC7700016 | built-in browser (non-OA) |
| Højfeldt 2018 | Højfeldt, Laugesen, Willumsen, … Helin. *Nat Struct Mol Biol* 25:225 | PMC5842896 | EuropePMC fullTextXML |
| Alabert 2020 | Alabert, Loos, Voelker-Albert, … Hasenauer, Imhof, Groth. *Cell Rep* 30:1223–1234 | (no PMCID) | local PDF → pdftotext |
| Jadhav 2020 | Jadhav, Manieri, Nalapareddy, … Shivdasani. *Mol Cell* 78:141–151 | PMC7376365 | built-in browser |
| Deal 2010 | Deal, Henikoff JG, Henikoff S. *Science* 328:1161–1164 | PMC2879085 | built-in browser |
| Kraushaar 2013 | Kraushaar et al. *Genome Biol* 14:R121 | PMC3983652 | EuropePMC fullTextXML |
| Deaton 2016 | Deaton, Gómez-Rodríguez, … Jansen, Kingston. *eLife* 5:e15316 | PMC4965263 | EuropePMC fullTextXML |
| Hathaway 2012 | Hathaway, Bell, Hodges, … Crabtree. *Cell* 149:1447 | PMC3422694 | built-in browser |
| Vignaux 2019 | Vignaux, Bregio, Hathaway. *PLoS ONE* 14:e0217699 | PMC6608945 | EuropePMC fullTextXML |
| Lee & Holder 2018 | Lee, Holder, Grau, … Armache. *Mol Cell* 70:435–448 | PMC5949877 | built-in browser |
