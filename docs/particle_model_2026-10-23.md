# Particle-based PRC2 lattice — build and first result (2026-10-23)

Code: `src/lattice/particle_lattice.py`. Figure:
`simulations/fig_particle_abundance_sweep.png`. Parameters from the six-part
census (`docs/parameters_consolidated_2026-10-23.md`).

---

## 1. What changed from the field model

| | field model | particle model |
|---|---|---|
| PRC2 | continuous activity at all 59 nucleosomes | **discrete molecules**, bind / dwell / catalyse / unbind |
| KDM6, p300 | rate fields | **discrete molecules** (both sub-stoichiometric) |
| HDAC | field | field — **correctly**, at 1 per 18–59 nucleosomes |
| binding modes | none | stable 10.0 s + transient 0.9 s, both catalysing |
| search | none | 210 s unbound between events |
| read-write | on catalysis | on catalysis (**confirmed correct** — H3K27me3 does not affect binding) |
| nucleation | PRC2 property | **PCL property** (Heenan: PRC2 core has no CGI preference at 2 nm) |
| basal turnover | 0.06 /h | **0.0044 /h** (measured) |

**Free consistency check:** tau_stable/(tau_stable + tau_search) = 10/220 = 4.5%
against the independently measured 4.6% stable bound fraction.

## 2. Three errors found during the build

1. **Exchange geometry (inherited from the field model, so it affects every
   earlier result).** Transcription-coupled exchange was applied across an 8 kb
   gene body, demanding **18 methylation steps/h from a domain that can supply
   ~1.4**. Histone turnover is promoter-local (Kraushaar 2013, Deaton 2016:
   promoters highest, PcG domains lowest). Moved to an NDR-centred peak.
2. **Cascade bug.** One random draw was reused across all three demethylation
   transitions in forward order, so a single draw took me3 -> me2 -> me1 -> u in
   one step, tripling the effective loss rate. Fixed by reversing the order.
3. **kcat normalisation.** I normalised to the fastest methylation step, which
   double-penalises the rate-limiting one ~9x. Published nucleosomal kcat is an
   aggregate turnover number and must refer to the limiting step. Before the fix
   the chain stalled: me1 accumulated to 0.22 while me2 sat at 0.03.

## 3. The capacity arithmetic reproduces the census crisis mechanistically

For a 40-nucleosome flank held at 50% me3, at n_prc2 = 2 and kcat 0.1 /min:

| | steps/h |
|---|---|
| capacity, unstimulated | 0.15 |
| capacity, with 9x allosteric boost | **1.38** |
| demand: basal demethylation | 0.09 |
| demand: replication (24 h, x3 steps) | 1.25 |
| **total demand** | **1.34** |

Marginal — and sustainable **only** because of the allosteric boost. So the
marked state sustains itself while the unmarked state cannot bootstrap, which is
the asymmetry Polycomb needs. The field model obtained that asymmetry from an
inflated turnover rate instead of from the enzymology.

## 4. RESULT: boundary depth peaks at INTERMEDIATE PRC2 abundance

24 h cycle, 8 cells per point, 1000 h, measured kcat range.

| PRC2/domain | kcat 0.47 core | plateau | **depletion** | edge (bp) | corresponds to |
|---|---|---|---|---|---|
| 3 | 0.000 | 0.076 | +0.076 | 1750 | low abundance, 10k targets |
| 11 | 0.009 | 0.717 | **+0.708** | 1370 | low abundance, 3k targets |
| **16** | 0.037 | 0.812 | **+0.775** | 1370 | **high abundance, 10k targets** |
| 30 | 0.203 | 0.885 | +0.682 | 1180 | — |
| 51 | **0.434** | 0.930 | +0.496 | 1180 | high abundance, 3k targets |

Three things fall out:

- **Too little PRC2 gives no domain** (n = 3: plateau 0.076).
- **Too much PRC2 fills the core** (n = 51: core 0.434, half-marked) and the
  boundary degrades even though the flank plateau is at its highest (0.930).
- **Depletion is maximal at n ≈ 16**, and the viable window is roughly
  n = 11–30 with kcat >= 0.3 /min.

The boundary also **moves inward** as abundance rises: edge 1750 -> 1370 -> 1180 bp.

The four literature-mapped regimes bracket the optimum, and only the lowest
(low abundance x 10k targets) fails outright.

## 5. Prediction

**There is an optimal PRC2 dosage for maintaining a promoter boundary, and
EZH2 overexpression should ERODE the H3K27me3-depleted core at active promoters
rather than strengthen repression.** In the model, going from n = 16 to n = 51
raises the flank plateau 0.81 -> 0.93 while the core rises 0.04 -> 0.43 and
depletion falls 0.78 -> 0.50.

This is directly testable in the EZH2-overexpression and MB systems, and it is
measurable with the metagene analysis already built
(`src/lattice/geometry_fit.py`): the prediction is a *narrowing and shallowing*
of the promoter dip with rising EZH2, not a uniform gain of H3K27me3.

## 6. Caveats

- kcat 0.47 /min is the single fastest measured nucleosomal value (PALI1 PIR,
  100 mM KCl). The viable window requires kcat >= 0.3, i.e. the top ~quarter of
  the measured 480-fold range. If the true in-cell kcat is nearer the median,
  the domain is not maintainable at any abundance in this model.
- n_prc2 assumes PRC2 is partitioned evenly across targets. Real targets differ
  in affinity, so a strong target would see more and a weak one less.
- KDM6 is present at n_kdm6 = 0.3 with an assumed PRC2-like residence; no
  residence measurement exists for any KDM6 enzyme.
- Single geometry, single transcription level, no PRC1, no H2Aub.
- 8 cells/point: depletion values carry roughly +/-0.03.

---

## 7. TEST OF THE OVEREXPRESSION PREDICTION — not possible with current data

### What was predicted

EZH2 CPM (Scott bulk RNA-seq, 55 samples): P7 GNP wt **187**, MB Ptch-het **320**
= 1.7x. Measured target counts: GNP 10,069, MB 11,537 = 1.15x. So

    n_prc2(MB)/n_prc2(GNP) = 1.71 / 1.15 = 1.49  ->  15.2 -> ~22.7

Interpolating the sweep at kcat 0.47: core 0.037 -> ~0.11, plateau 0.81 -> ~0.85,
depletion 0.775 -> ~0.73, edge 1370 -> ~1280 bp. **Prediction: MB has a
shallower, slightly narrower promoter dip than GNP.**

### What the measurement can actually resolve

Coverage-normalised metagene dip from epic2 peaks, all 20 H3K27me3 samples:

| | |
|---|---|
| samples passing QC | 9 of 20 (7 over-called at plateau >0.9, 4 too sparse) |
| predicted effect (relative dip) | **~0.09** |
| within-condition replicate sd | **0.043** |
| between-condition spread | 0.091 |

Worse than the variance ratio suggests: **width estimates are bimodal.** A sample
either resolves a dip (~1,200–2,200 bp) or returns the window edge
(10,400–12,000 bp = no dip detected). Within P7 H3K27me3, rep2 returns 10,400
and rep3 returns 2,200 — **the same condition both shows and fails to show a
boundary depending on replicate.** Among MB samples the relative dip spans
0.019–0.319, a 17-fold range.

Most conditions have n = 1. The test has no usable power.

### Why, and what would fix it

Peak calling **thresholds a continuous signal**, so a coverage dip registers as a
peak-coverage dip only if it crosses the caller's threshold. The measurement is a
nonlinear, threshold-dependent transform of the quantity the model predicts. This
is the same failure mode that made 11 of 20 samples unusable.

This vindicates `docs/geometry_data_spec_2026-10-23.md` precisely: the test needs
**coverage-level metagenes at <=200 bp bins over -6/+9 kb**, not peak calls.

### A second problem the data cannot solve either

Even with a perfect measurement, GNP -> MB does not isolate EZH2:

| | GNP P7 | MB | fold | model's predicted effect on the core |
|---|---|---|---|---|
| Ezh2 | 187 | 320 | 1.7x | **fills** the core |
| Gli1 | 108 | 563 | 5.2x | **widens/deepens** the core (KDM6B arm) |
| cycle length | longer | shorter | — | lowers the plateau |

Three variables move together and two of them predict the opposite of EZH2.
Whatever the dip does between GNP and MB, it cannot be attributed.

**The clean test is the EZH2 conditional knockout with H3K27me3 ChIP** — a
single-variable perturbation. The mice and the RNA-seq exist
(`P7_GNP_EZH2ko`, `MB_EZH2ko`); the H3K27me3 ChIP does not. That is one
experiment, and it would test the prediction directly in the falsifiable
direction (less EZH2 -> domain fails to form, per the n = 3.3 row of the sweep).
