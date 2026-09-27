# Protein-noise work: handoff to the modeller

**Read this first.** It is the entry point for the six detail documents listed at the end. Written
for whoever implements the noise modules (brief tasks T1–T8). Nothing in `src/` has been changed;
`simulations/validate_v44.py` is **28/28** and untouched.

---

## 0. TL;DR — the six things that change what you build

1. **The Hedgehog arm is missing a `/60`.** Its constants were ported from v42 (hours) into v44's
   minute-frame core. Fixing it is **free**: 28/28 → 28/28, zero flips. **Do this before anything
   else** — bursts injected into a module that equilibrates in 4 minutes are averaged to nothing.
2. **Almost nothing fluctuates fast.** 99% of proteins have τ > 8 h, 0% < 2 h. The brief treats
   frozen heterogeneity and fast fluctuation as equally plausible; for most of the proteome it is
   frozen heterogeneity. Default to **long** drift correlation times.
3. **`Skp2` is the burst candidate the brief missed** — the only v44 species with genuinely fast
   noise (τ = 4.0 h), and it sits on the Skp2–p27 commitment switch.
4. **PRC2 lowers burst frequency** (~1.6× at matched expression, two independent datasets). It may
   also raise burst size, but that arm is weak — see §4.1, which was **revised after review**.
5. **mRNA degradation is cycle-modulated** for 90% of genes (median 3.3×). v44 treats it as
   constant. This may be the missing mechanism behind the model's worst standing target.
6. **Do NOT "fix" `kDeEZm`** to its measured value — it breaks 7 targets. Tested. See §2.

---

## 1. What exists now

| artefact | what it is |
|---|---|
| `data/priors/turnover_burst_priors.csv` | **924 rows**, 313 numeric slots filled, every value sourced, 512 explicit `NA`. Never imputed. |
| `data/priors/SOURCES.md` | URL + date + sha256 + **dilution convention** per file |
| `simulations/build_turnover_priors.py` | regenerates the CSV from `data/priors/raw/` |
| `simulations/analyze_noise_budget.py` | per-species CV and correlation time vs the model's own timescales |
| `models_external/biomodels/` + `SOURCES.md` | 11 reference SBML, sha256-verified, load-tested |
| `.claude/skills/biomodels/` | skill + `bmfetch.py` for the BioModels REST API |

Nine datasets ingested: Schwanhäusser 2011, Mathieson 2018, Larsson 2019 (S1/S2/S3), Sigal 2006,
Eden 2011, McShane 2016, Ochiai 2020, Herzog 2017, Battich 2020, Dörrbaum 2018. **No literature is
outstanding.** Remaining gaps are lab data only (§6).

Raw files are gitignored (repo policy); `SOURCES.md` has checksums so they are reproducible.

---

## 2. Tested changes — evidence, not opinion

All run against `validate_v44.py` via its `VALIDATE_PARAMS` override. **Baseline 28/28.**

| # | change | result |
|---|---|---|
| **F1** | divide the 28 HH-arm rate constants by 60 | **28/28 → 28/28**, no flips, max move 3.7% |
| **F3** | set `kDeEZm` to the measured 9.85 h | **28/28 → 21/28** — do not do this |

**F1 detail.** Constants ported from `build_model_v42_mycn.py`, whose time unit is hours
(`t_end=96`/`168`), into v44's minute-frame core. `k_Cd_deg = 1.0` meant 1.0/h (t½ 42 min) and is
being integrated as 1.0/min (t½ 42 s). Measured on the built model: `Cd` relaxes in **3.8 min**.

It is free because scaling every rate in a module by one factor leaves steady states
algebraically invariant, and that arm's calibration is entirely steady-state. The two targets that
moved most both moved **toward** target, including `EZH2 transcript S/G0` (1.563 → 1.621).
One real behavioural change: within-cycle CyclinD1 amplitude drops ~2.5× — the intended effect,
and what makes the noise work possible. Constant list and repro: `protein_noise_fixes.md` §F1.

**Open item on F1:** under `MB+HHi+EZH2i` a side harness found 2 division peaks baseline vs 6
rescaled. That is the condition `validate_v44.run()` flags as stiff and the suite passes either
way — likely an integration artefact, but confirm before merging.

**F3 detail.** All 7 broken targets are EZH2 cell-cycle-phase gradients, and all collapse to ~1.0:
an mRNA living 9.85 h cannot resolve 2–6 h phases. v44's short `kDeEZm` is load-bearing. **But see
§4.2 — this may be a missing mechanism rather than a data conflict.**

---

## 3. Parameters to use

Best current values. Full provenance per row in the CSV.

### 3.1 Half-lives — prefer these sources in this order

**mRNA: Herzog 2017 (SLAM-seq, mouse ESC) over Schwanhäusser.** Metabolic labelling, no
transcription block. Schwanhäusser runs **2.31× long** on 3071 shared genes.

| gene | Herzog t½ (h) | Schwanhäusser | note |
|---|---|---|---|
| `Ccnd1` | **3.37** | 14.0 | was the project's biggest lever; now ~4× not 16× |
| `Ezh2` | **3.71** | 9.85 | model uses 1.18 h — see §4.2 |
| `Cdkn1a` | 1.41 | — | |
| `Gli1` | 3.05 | — | **only source** |
| `Gli2` | 1.90 | — | **only source** |
| `Mycn` | 2.00 | — | **only source** |
| `Ccne1` / `Ccnb1` / `E2f1` / `Trp53` | 4.01 / 3.09 / 3.38 / 5.99 | — | only source |
| `Skp2` | 2.23 | 2.82 | |

Still missing everywhere: `Gli3`, `Ptch1`, `Smo`, `Rb1`, `Cdkn1b`.

**Protein: match the pool, not the gene.** McShane fits a 2-state model (fast unassembled pool
`kA` → stable assembled pool `kB`). For a species representing a **complexed** form, use `kB`.
Schwanhäusser runs **1.74× long** vs McShane in the same cell line, and the bias is concentrated
in non-exponential proteins (**3.41× for NED** vs 1.24× for ED) — a 1-exponential fit to 2-pool
kinetics. **Rule: Schwanhäusser is fine for ED proteins, unusable for NED.**

EZH2 worked example — v44's `kDeEZ` (τ = 14.3 h) against the candidates:

| candidate | t½ | τ with dilution | model/bio |
|---|---|---|---|
| Schwanhäusser 1-exp | 22.0 h | 13.4 h | 1.1× |
| **McShane assembled pool `kB`** | **12.2 h** | **10.0 h** | **1.4×** |
| McShane 1-state blend | 5.0 h | 5.5 h | 2.6× |
| McShane fast pool `kA` | 0.5 h | 0.7 h | 21.5× |

EZH2 in the model is catalytic PRC2, so `kB` is correct. **The constant is well-placed — leave it.**

### 3.2 Dilution classes (brief T4)

| source | measures | class | add `-µ·x`? |
|---|---|---|---|
| Mathieson (post-mitotic) | degradation only | (b) | **yes** |
| Schwanhäusser | degradation only — **established empirically, see below** | (b) | **yes** |
| Eden (bleach-chase) | degradation **+ dilution** (SOM Box S1) | (a) | **no** |
| McShane | cycling cells, convention unverified | (c) | unknown |
| Heldt-inherited constants | fitted in cycling RPE1 | (a) | **no** |

Schwanhäusser's class was settled from data, not the Methods: on the 413 most stable proteins, the
with-dilution hypothesis predicts a median t½ of **23.1 h** against **121.3 h** observed (and
116.9 h predicted without). The gap vanishes exactly where dilution would be largest.

**`Rb1` = 18.9 h** (Mathieson, post-mitotic). At T_cc = 16 h dilution (0.043/h) and degradation
(0.037/h) are **the same size** — the Zatulovskiy premise with a number on it.

### 3.3 Noise budget — current best estimates

`python3 simulations/analyze_noise_budget.py`. Rows without a measured protein half-life are
flagged `[!]`; their τ is a dilution-only assumption, not a measurement.

| gene | CV_p | **τ_c (h)** | reading |
|---|---|---|---|
| **`Skp2`** | 0.94 | **4.0** | **fast + loud — the burst candidate** |
| `Ccnd1` | 0.77 | 4.9 | intermediate |
| `Ccna2` | 1.69 | 6.9 | loud |
| `Chek1` | 1.31 | 7.5 | |
| `Cdk4` | 0.45 | 11.6 | |
| `Suz12` | 0.76 | 13.2 | |
| **`Ezh2`** | **0.66** | **13.4** | **slow — frozen heterogeneity** |
| `Eed` | 0.52 | 16.1 | |
| `Cdk2` | 0.54 | 17.2 | |

Genome-wide over 3298 complete genes: **99% have τ > 8 h, 74% > 16 h, 0% < 2 h.**

**Unit hazard — do not repeat this mistake.** Larsson counts *detected* molecules (UMI capture
p ≈ 0.21, calibrated over 3080 shared genes). Thinning scales `B` and `⟨m⟩` together, so
`(p + B)/⟨m⟩` is capture-invariant **within one dataset** — but mixing Schwanhäusser's absolute
`⟨m⟩` with Larsson's deflated `B` understates noise by ~1/p (measured 4.4–15.6× on our genes).
Take `⟨m⟩` and `B` from the same source. This bug was in an earlier version of the analysis script
and is fixed.

### 3.4 Extrinsic / global OU

Sigal 2006, 40 CD-tagged clones: **τm = 1.75 generations** [0.8–4.3], **CV = 0.20** [0.12–0.28].
Both confirm the brief's §3.1 assumption.

Joining Eden's half-lives to Sigal's mixing times on the **same clones**: observed memory is
**2.7× longer** than the removal rate alone predicts (n=6, range 2.0–4.1, none near 1.0). So the
extrinsic component is real and dominant, not a bookkeeping artefact of dilution.

**Calibration target:** with the model's own removal rates in place, the added drift should extend
each species' autocorrelation half-time by **≈2.7×**. n = 6 — treat as order 2–3.

### 3.5 Burst kinetics in real units

Larsson reports frequency per mRNA-degradation-rate. Herzog supplies the denominator **in the same
cell type**, so `f[1/h] = f[per δ] · ln2 / t½(mRNA)`:

| gene | B | **h per burst (per allele)** |
|---|---|---|
| `Cdkn1a` | 9.8 | **0.8** |
| `Ezh2` | 3.1 | 3.0 |
| `Ccnd1` | 5.7 | 7.2 |
| `Mycn` | 18.5 | 10.6 |
| **`Skp2`** | **12.8** | **21.9 — about once per cell cycle** |

A diploid cell fires twice as often. GNP cycle ≈ 22 h.

---

## 4. Structural changes worth making

### 4.1 PRC2 lowers burst frequency; the burst-size arm is weak  *(revised after review)*

Ochiai 2020, genome-wide in mESC with ChIP-seq at promoter and gene body, state that PRC2 promoter
occupancy "correlated inversely with burst frequency, while they correlated positively with burst
size and normalized intrinsic noise". **An earlier version of this handoff recommended coupling
H3K27me3 to both `k_on` and `B` as a default structural change. That was overstated.**

Ochiai give this as one sentence pointing at a heatmap — no coefficients, no statistics — and
their causal test is four genes, one of which (*Peg3*) goes the opposite way. An independent
rebuild from ChIP-Atlas `ES_cells` binding × Larsson ESC bursting
(`simulations/test_prc2_bursting.py`) gives:

| | marginal r | partial r given expression | bound/unbound at matched expression |
|---|---|---|---|
| burst **frequency** | **−0.23** | −0.02 | **0.61×** (Suz12), 0.96× (Ezh2, n=16) |
| burst **size** | **−0.06** | +0.09 | 1.15× (Suz12), 1.01× (Ezh2) |

The confounder matters: PRC2 marks lowly expressed genes and both f and B correlate positively
with expression, so marginal correlations are negative for both and only expression-conditioned
estimates are interpretable.

**Recommendation:** couple PRC2/H3K27me3 to **`k_on` (down)** — supported, ~1.6×. Treat the
**burst-size arm as a sensitivity**, not a default: it is worth <1% of variance and marginally the
sign reverses. The earlier "frequency-only coupling understates the noise" claim is technically
true but worth only **~6% on CV**, not a qualitative correction.

Still worth keeping: they perturbed the **complex**, not catalysis. Since the correlation is with
PRC2 *occupancy* and tazemetostat leaves PRC2 bound, this predicts **EZH2i de-represses the mean
while leaving PRC2-imposed burstiness largely intact** — testable, and it follows from the repo's
existing catalytic-inhibitor distinction.

### 4.2 Cycle-modulated mRNA degradation — and a possible resolution of F3

Battich 2020 (scEU-seq, **RPE1 — Heldt's line**): **90% of genes vary their mRNA degradation rate
>2-fold across the cycle** (median 3.3×). 35% "cooperative" — degradation drops when synthesis
rises, amplifying the swing. v44 uses constant mRNA degradation throughout, so it must generate the
whole cycle amplitude from synthesis alone and needs a **larger** synthesis swing than reality.

EZH2 is cooperative (cos = −0.554): degradation varies **2.96×** peaking early (pseudotime 0.20),
synthesis **4.02×** peaking later (0.39), anti-phased so they reinforce the S-phase rise.

**Hypothesis worth testing:** v44 needs an artificially short constant `kDeEZm` (1.18 h) because it
is compensating for missing cycle-modulation. With degradation free to drop ~3× during S, a longer
mean half-life could still give the observed gradient — which would make **F3 a missing mechanism,
not a data conflict**, and might also close the `EZH2 transcript S/G0` miss (1.563 vs ~2, the
README's worst).

Concrete test, not yet run:
```
kDeEZm_eff := kDeEZm * (1 + a_mk*(gS - <gS>))     // a_mk < 0 = cooperative
```
with `gS` an existing S-phase indicator (`Ce`, `Ca`, `aRc`), behind `with_cyclic_mrna_deg` (default
OFF). Ask whether S/G0 reaches ~2 **with `kDeEZm` set to Herzog's 3.71 h**. Prior: fold-range ≈ 3,
phase ≈ a quarter-cycle ahead of synthesis.

### 4.3 Non-exponential degradation

McShane: 11% of NIH3T3 and 10% of RPE proteins are NED, median **138×** split between fast
unassembled (t½ 99 min) and stable assembled (t½ 181 h) pools. Among v44 species, **`Cdk4`, `Pcna`
and `Ccnb1` are NED**. v44 degrades free and complexed pools at the same rate — structural, not a
parameter choice. Relevant to `CeP21`, `CaP21`, `RbE2f`, `E1C1`, `iPcna`, `iRc`.

---

## 5. Revised task order

1. **F1** — the `/60` rescale. Tested, free, unblocks everything.
2. **T4 dilution audit** — cheap now; the CSV carries `growth_corrected` per row and §3.2 gives the
   classes.
3. **§4.2 test** — cyclic `kDeEZm`. Could convert a standing validation miss and a data conflict
   into one structural fix.
4. **T5 drift.** Keep the exact OU discretisation and segmented integration — they are sound.
   Default to **long** correlation times (§0.2). Only four v44 species have every input for the
   full CV formula: `Ccnd1`, `Ccna2`, `Ezh2`, `Skp2`.
5. **T6 bursts.** Retarget: **`Cdkn1a`** (B = 19, 96th percentile) and **`Skp2`** (τ = 4.0 h, the
   only fast node, on the commitment switch) ahead of `Ccnd1` (τ ≈ 5 h, and only after F1).
   **Drop `Ezh2`** — slow, and its burst size is 45th percentile, i.e. unremarkable. `Mycn` has
   burst data in ESC only and no protein half-life anywhere — sensitivity only.
6. **T2 / T7 / T8** — as specified.

---

## 6. What is still blocked — all of it lab data

1. **HCI per-cell tables.** Every CV in §3.3 is a *prediction* with nothing to test it against.
   Need: mean intensity (concentration) **and** integrated, DNA content, EdU, nuclear area,
   condition, plate/well/batch. **Raw rows, not means ± SD.**
2. **An absolute scale (Ω).** `DATA_PROVENANCE.md` says species are in arbitrary units calibrated
   to ratios. T2 (Gillespie) and T6 (discrete bursts) need **counts**. Any one of: GNP nuclear
   volume, one absolute copy number, one quantitative western.
3. **Measured T_cc** for GNP and MB. Scales the dilution analysis and converts §3.4 from
   generations to hours.

Full spec: `protein_noise_data_requests.md` §A.

---

## 7. How much to trust a borrowed number

Measured, not assumed (`protein_noise_celltype_transfer.md`). Log–log r on protein half-life:

| comparison | r |
|---|---|
| neuron vs glia — same culture, method, species, lineage | 0.63 |
| mouse fibroblast vs human epithelial, one method | 0.62 |
| 5 primary cell types, one method | 0.57–0.79 |
| **same cell line, different method** | **0.54** |

**Flat at ≈0.6 regardless of what varies** — matching the cell type more carefully buys almost
nothing, and method disagreement is as large as cell-type disagreement. So:

* **borrow the rank; rescale the level** (Mathieson medians span 3.6×). v44 being ratio-calibrated
  already absorbs a global scale.
* **protein and mRNA half-lives**: fine to borrow, propagate ~2× per-gene uncertainty.
* **burst size** (r = 0.50): borrow with caution — partly TATA-encoded, so it travels for a reason.
* **burst frequency** (r = **0.13**): **do not borrow.** No sequence determinant (R² = 0.008) and
  it does not survive a change of cell type. Chromatin is the only route — which is the argument
  for using the lab's own ChIP-seq (`protein_noise_bursting_chromatin.md` §6).

Also: GNPs **cycle**, so the cycle-state-matched references (NIH3T3, RPE1, mESC) are arguably a
safer borrow than the lineage-matched but post-mitotic ones.

---

## 8. Corrections made during this work — do not re-derive the wrong versions

| claim | status |
|---|---|
| "The `/60` rescale needs a recalibration pass" | **wrong** — tested, it is a mechanical edit |
| "EZH2 is slow *and quiet*, CV ≈ 0.26" | **wrong** — dataset-mixing bug; CV ≈ 0.66–0.89. Slow, but not quiet. The recommendation stands for the τ reason, not the amplitude reason |
| "Set τ_ext to 1.75 generations would double-count dilution" | **too strong** — Eden×Sigal shows memory is 2.7× longer than removal alone; the extrinsic term is real and dominant |
| "`kDeEZm` is 8.4× too fast and should be corrected" | **wrong** — breaks 7 targets; and §4.2 suggests it is compensating for a missing mechanism |
| Battich rates read as min–max normalised | **wrong** — they are log2(rate ÷ cycle median) |
| "PRC2 modulates burst frequency and size oppositely; frequency-only understates the noise" | **overstated** — frequency arm holds (~1.6×), size arm is <1% of variance and marginally sign-reversed. See `protein_noise_bursting_chromatin.md` §2b |

---

## 9. Detail documents

| doc | contents |
|---|---|
| `protein_noise_fixes.md` | the tested change list F1–F5, with per-target diffs and repro commands |
| `protein_noise_data_exploration.md` | dilution audit, genome-wide noise budget, the Skp2 finding, Eden/McShane |
| `protein_noise_bursting_chromatin.md` | Ochiai PRC2 result, Larsson promoter regression, the mESC join |
| `protein_noise_cellcycle_rates.md` | Battich, cycle-modulated degradation, the §4.2 test |
| `protein_noise_celltype_transfer.md` | how much a borrowed rate is worth |
| `protein_noise_data_requests.md` | what is still needed from the lab |
| `protein_noise_plan.md` | original reconnaissance (superseded in parts by this document) |
