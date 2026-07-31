# Calibration & Validation

> How the v44 model is scored, how it has been fit over ~4 months, what the current best set achieves, and where it structurally cannot go. This is the "how do we know it's right (and where it isn't)" page. See also [Model Architecture](02_model_architecture.md), [Evolution Timeline](03_evolution_timeline.md), [Things Tried & Abandoned](05_things_tried_and_abandoned.md), and [Data & Evidence](07_data_and_evidence.md).

---

## 1. The single source of truth: `validate_v44.py`

All calibration for the active model funnels through one harness:
`Ezh2_CyclinD1_Sphase/simulations/validate_v44.py`. It is the direct successor to `validate_v42.py` (the old Gérard–Goldbeter model). Everything else — the optimizers, the figure scripts, the identifiability probes — reads its output or drives it.

**What it does, in order:**

1. **Builds the model once** (`build_model_v44` from `src/build_model_v44_heldt.py`) with the current default feature flags, then reuses the one roadrunner instance across all conditions via `reset()` + per-condition parameter sets. (This reuse is the source of the `reset()` gotcha — see §7.)
2. **Runs ~15 experimental conditions** (`CONDITIONS` dict): GNP ±SHH, GNP+HHi, GNP+EZH2i, GNP+HHi+EZH2i, GNP serum-starved, GNP+CDK4/6i, MB, MB+HHi, MB+EZH2i, MB+HHi+EZH2i, MB+CDK4/6i, MB+CDK4/6i+EZH2i, MB+HU. Cell type is encoded structurally, not by a single knob — MB = `Ptch1_copy_number=0.1` (Hh-constitutive), `MYCN_expr` elevated, `Cd2_expr=3.59`, and the INK4/CIP-KIP tones `p16/p18/p19/kSyP21`. See [Mechanisms Explored](04_mechanisms_explored.md).
3. **Classifies cell-cycle phase, counts divisions, measures settled levels**, then scores each target with a single `check(name, actual, target, tol)` function.
4. **Writes `validation_v44_results.json`** (path overridable via `VALIDATE_RESULTS`) — the machine-readable output every optimizer parses.

The model is **real time in minutes** (unlike v42's ε-scaled "hours"). Levels are `mean_settled` (mean over `t ≥ 4000 min`); the transcript readout is `Cd_mRNA` (matches bulk RNA-seq), not the fast protein `Cd`.

### How a target is scored — `check()`

```python
def check(name, actual, target, tol=0.20, unit='', zero_ok=False):
    if zero_ok:      ok = abs(actual) <= 0.001          # exact-zero (arrest imposed by construction)
    elif target==0:  ok = abs(actual) <= 1              # "0 divisions" = <1 division counted
    else:            ok = abs(actual - target)/abs(target) <= tol   # relative band
```

So a quantitative target passes if the model is within a **relative tolerance band** (`tol`, per-target). This is a **step function** — it is deliberately loose on noisy data (§6) but, as the identifiability review notes (§8), it is exactly zero inside the band and therefore a poor *cost surface*.

### How the readouts are computed

- **Divisions / period** — `count_divisions()`: peaks of `MPF` (mitotic CyclinB-Cdk1) after a 3000 min settle, `find_peaks(prominence=0.15, distance=200 min)`. Period = mean inter-peak time in hours. A "0-division arrest" is the absence of MPF peaks.
- **Phase classification** — `classify()`, anchored to the manuscript IF assay (phospho-Rb Ser807/811 + DAPI/PCNA):
  - **S** = *active synthesis flux* `vfork*aRc > SYN_THR (0.02)` **and** `Dna < 0.98`. Crucially, S is the synthesis flux, not mere replication-complex (`aRc`) presence — under HU, forks slow (`vfork≈0.1`) so stalled cells are BrdU-negative and fall out of S. This single design choice is what makes the two HU folds structurally coupled (§5).
  - **G2** = `Dna >= 0.98` (replication complete).
  - **G0** = pre-S with Rb **hypo**phosphorylated (`pRb < PRB_G0_THR = 1.5`) — the 2026-07-15 transient-G0 reframe made pRb-hypophosphorylation the functional arrest marker, replacing the old p27-level marker.
  - **G1** = pre-S with Rb phosphorylated.
- **Count-fraction vs duration-fraction** — flow cytometry measures the *count-fraction* (cells sitting in each phase), not the duration-fraction (phase-time / period). In exponential growth there are ~2× more just-divided than about-to-divide cells, so `classify()` applies the **age-density λ-correction** `n(a)=2λe^{−λa}`, `λ=ln2/Tc`, re-weighting each settled time point by the time since its last division. Comparing raw duration% to flow count% was a real earlier bug that spuriously passed S and G2. (If <2 MPF peaks are found, it silently falls back to duration-fractions — a known hazard flagged in the identifiability plan.)

---

## 2. The validation-target set

The canonical spec is `data/validation_targets.json` (machine-readable, edit first) with a readable view in `docs/validation_targets.md`. It was built (commit `6b4383d`, 2026-07-25) after a **3-source cross-audit** (harness vs the Chahin compendium `docs/Ezh2_CcnD1_model_targets.md` vs `docs/DATA_PROVENANCE.md`) that JP requested to stop MB↔GNP conflation and provenance drift. Every target now carries a **cell type** (P7 GNP / MB / MB↔GNP), a **provenance tag** (manuscript-fig / ext-lit / unpub / prediction / assum), a current `model_actual`, and pass/fail.

The harness scored **29 checks** at the time of the spec (**28/29**), grew to **32** when the two-cyclin (CyclinD2) checks and dynamic-CDKI module were baked (**30/32**, current). The count is conditional: the 3 CyclinD2 checks only score when `w_Cd2 > 0`.

### Full target table

Values below are as recorded in `docs/validation_targets.md` (default model, 28/29 snapshot, 2026-07-25) plus the three CyclinD2 checks from `validate_v44.py`. `tol` is the relative band.

**P7 GNP**

| Target | Quantity | Target | Tol | Model | ✓ | Provenance |
|---|---|---|---|---|---|---|
| CyclinD1 GNP+HHi / GNP | Ccnd1, vismo vs +Shh | 0.157 | 0.40 | ~0.10 | ✓ | manuscript §A |
| MYCN GNP+HHi / GNP | Mycn, HHi vs +Shh | 0.78 | 0.30 | ~0.91 | ✓* | manuscript §A (soft) |
| Gli1 GNP+HHi reduction | 1 − Gli1(HHi)/Gli1 | 0.99 | 0.05 | — | ✓ | manuscript §A |
| EZH2i CyclinD1 fold | Ccnd1, **Tazemetostat** vs +Shh | 2.2 | 0.42 | **1.63** | ✓* | manuscript Fig 3C |
| EZH2 G0/cycling | EZH2 protein, starved vs cycling | 0.6 | 0.40 | 0.61 | ✓ | manuscript §D |
| EZH2 Palbo mRNA drop | EZH2 mRNA, CDK4/6i vs +Shh | 0.44 | 0.42 | 0.59 | ✓ | Fig 4F ⚠ *measured in MB* |
| GNP+Shh proliferates | ≥1 division | yes | — | yes | ✓ | ext-lit (Wechsler-Reya '99) |
| GNP−Shh arrests | 0 div | 0 | — | 0 | ✓ | ext-lit |
| GNP+HHi arrests | 0 div | 0 | — | 0 | ✓ | inferred |
| GNP serum-starve arrests | 0 div | 0 | — | 0 | ✓ | assum (tautological: `k_Cd_translation=0`) |
| **Period GNP** ⚠ | cycle length | **~16 h** | 0.30 | 22.83 | ✗ | ext-lit (Nakashima 15.9h + Contestabile 16.25h) |

**MB**

| Target | Quantity | Target | Tol | Model | ✓ | Provenance |
|---|---|---|---|---|---|---|
| CyclinD1 MB+HHi / MB | Ccnd1, vismo vs MB | 0.144 | 0.45 | 0.11 | ✓ | manuscript §A (86% drop) |
| MYCN MB+HHi / MB | Mycn, HHi vs MB | 0.86 | 0.25 | 0.94 | ✓ | manuscript §A |
| MB proliferates | ≥1 div | yes | — | yes | ✓ | manuscript |
| MB+CDK4/6i (palbo) ⚠ | pRb⁺ cycling + reversibility | 16–18%, reversible | — | 0-div (proxy) | — | manuscript (harness-mismatch) |
| MB+CDK4/6i+EZH2i | EZH2i does *not* rescue | 0 | — | 0 | ✓ | **prediction** (Fig S9I,J) |
| MB 2N (G0+G1) count% | 2N fraction | 68.2% | 0.30 | 68.7 | ✓ | manuscript §F (microscopy) |
| MB S count% | S fraction | 15.7% | 0.42 | 21.9 | ✓ | manuscript §F |
| MB G2+M duration | pHH3/BrdU | 2.5 h | 0.55 | 3.39 | ✓ | ext-lit (Fujita), soft |
| MB HU S fold | S, HU/DMSO | 1.36 | 0.40 | 0.95 | ✓ | manuscript §F |
| **MB HU G2 fold** | G2, HU/DMSO | 0.23 | 0.40 | **0.36** | **✗** | manuscript §F — *the perennial miss* |
| EZH2 transcript S/G0 | within-cycle EZH2 mRNA | 2.0 | 0.45 | 1.53 | ✓ | Fig 4A/4B |
| EZH2 protein G2/G0 | within-cycle EZH2 protein | 1.48 | 0.45 | 1.86 | ✓ | Supp 6D (IF) |
| HU EZH2-in-S boost | EZH2, HU-S vs cycling-S | 1.31 | 0.35 | 1.25 | ✓ | Fig 4G |

**MB ↔ GNP cross-ratios**

| Target | Target | Tol | Model | ✓ | Provenance |
|---|---|---|---|---|---|
| CyclinD1 MB/GNP | 5.07 | **0.20** | 4.61 | ✓ | Fig 4I — **load-bearing, tightest band** |
| MYCN MB/GNP | 2.80 | 0.30 | 2.71 | ✓ | manuscript §A (~2.86) |
| Gli1 MB/GNP | 6.90 | 0.40 | 7.15 | ✓ | raw RNA-seq |
| EZH2 MB/GNP | 2.05 | 0.35 | 2.03 | ✓ | Fig 4J (soft) |
| Skp2 MB/GNP | 2.3 | 0.30 | 2.34 | ✓ | unpub (scRNA + MB table) |

**CyclinD2 (added 2026-07-29, only scored when `w_Cd2>0`)**

| Target | Target | Tol | ✓ | Provenance |
|---|---|---|---|---|
| CyclinD2 MB/GNP | 2.39 | 0.25 | — | Chahin RNA-seq |
| CyclinD2 GNP+HHi/GNP | 0.66 | 0.30 | — | RNA-seq (Hh-buffered) |
| CyclinD2 MB+HHi/MB | 0.59 | 0.45 | — | RNA-seq (trend, padj 0.076 → wide band) |

**Two flagged harness-mismatches** (⚠) — the check is a proxy that should change: **Period GNP** (harness still targets 22 h; literature says ~16 h, see §4) and **MB+CDK4/6i** (JP's culture data is a *partial, reversible* 16–18% residual pRb⁺ population; the harness scores a binary 0-division on one lineage). The **MB+CDK4/6i+EZH2i "no rescue"** is a *model prediction* (Fig S9I,J), not experimental — it was mislabeled and corrected in the audit.

---

## 3. Calibration methods, chronologically

The calibration approach itself evolved. From the git history (`simulations/`, 233-commit repo):

**Phase 0 — v43 lock (June 2 2026).** The precursor GG-core model reached 12/12 with the ε-in-S intra-S checkpoint. A verification study (`docs/v44_verification_findings.md`, `diag_v44_*`) proved that S-phase duration is a **structural invariant** of the Gérard–Goldbeter relaxation oscillator — no concentration (Cdc25/Cdk1/APC, even hard-gated) can lengthen S; only the global clock can. This *justified* moving to the Heldt-2018 explicit-replication core (v44), where S is a real rate-limited process (`vfork`, `aRc`, `Dna`).

**Phase 1 — between-condition calibration framework (June 2–5).** Introduced the `VALIDATE_PARAMS` env-var override hook so a driver can inject a parameter dict without editing source. First HH/MYCN ratio fits (`v44_optimize_hh.py`, `v44_recalibrate_gli.py`) hit the central discovery that **CyclinD1 *level* is linearly coupled to cell-cycle *drive*** (`Rb-P = kPhRbCd*Cd`), so a 7× transcript ratio → 7× drive → stiffness/instability. Fixed with the **saturating CyclinD1→Rb drive** (`with_cd_sat`, `kPhRbCd*Cd/(K_CdRb+Cd)`) that decouples transcript from bounded drive. The vismodegib-on-MB data (`MB_GDC0449`, CyclinD1 drops 86% to 0.144) recalibrated the whole HH module and made the population EZH2i-rescue emerge.

**Phase 2 — wide random search (July 8–9).** `overnight_search.py` then `v44_wide_search.py`: 17 parameters, hybrid local+wide search, against ~11 soft targets under **8 hard figure-critical guards** (GNP cycles/arrests, MB cycles, MB+HHi arrests, rescue emerges, CDK4/6i no-rescue). Seed 1 ran 3,334 evaluations → 8/8 guards + 11/11 soft. This recovered the raised commitment threshold (`kPhRbCd=0.35`) *with the rescue intact* — something hand-tuning couldn't do. Documented in `docs/v44_recalibration_and_search_results.md`. (Seed 2 died on a parallel-launch bug — run seeds sequentially.)

**Phase 3 — target-flexibility reform (July 10).** JP's directive (memory `target-flexibility-noisy-data`): **stop exact-matching noisy data.** The optimizers were chasing measurement noise and landing on knife-edges. Two changes: (a) widened `validate_v44` tolerances (RNA-seq folds ±30–40%, flow fractions ±30%, ChIP "half" → [0.4,0.6]); (b) switched optimizer loss from exact-distance to **band-hinge penalties**: `penalty = max(0, |x−t|/t − band)`, zero inside the noise band. Qualitative targets (arrest = 0 div, cycles > 0) stayed tight.

**Phase 4 — module-specific joint optimizers (July 10–14).** A family of focused subprocess optimizers, each freezing the validated machinery and searching one module's knobs. They share a design: sample params → write `VALIDATE_PARAMS` → run `validate_v44.py` in a subprocess → parse `validation_v44_results.json` → band-hinge loss + `10 × n_hardfail` → `ThreadPoolExecutor` over 8 workers → random phase then hill-climb (`perturb` with shrinking scale). Inventory:
  - `optimize_chain.py` — serial me1→me2→me3 methylation chain (promoted to default July 11, 28/28).
  - `optimize_prc2.py`, `optimize_cdh1_ezh2.py` — PRC2 occupancy / EZH2 stability.
  - `optimize_jmjd3_chip.py`, `optimize_no_chip.py` — Gli-Jmjd3 eraser vs ChIP direction (eraser turned OFF by default July 18).
  - `optimize_skp2.py` — Skp2 mitogen-dose fold.
  - `optimize_featureC.py` — commitment carryover.
  - **`optimize_twostep_g0.py`** — the load-bearing one (below).

**Phase 5 — the July 14 joint re-optimization to 28/28.** An adversarial critical review found the shipped 26/28 default was net-negative: the CyclinD1 MB/GNP fold was *really 7.1* (passing only on a wide 0.45 band, not 5.07), inflated by a **birth-p27 = 1.8 "crutch"** (3× the ~1.33 transcript grounding); there were **two real fails** (GNP+HHi 0.293; MB HU G2 0.685), not one; and the MB G0 rested on a saturated switch. `optimize_twostep_g0.py` (14 knobs, seeded `20260714`) genuinely landed the fold at **5.6** on a **tight 0.20 band**, recovered a real MB G0 dwell (~12%) at birth-p27 **1.39**, and fixed both standing fails (GNP+HHi 0.216, MB HU G2 0.316). Two mechanisms: a structural p27-toggle fix (CDK4/6↔p27 mutual antagonism) plus a mitogen-dose EZH2 re-fit (`Kez_cd 3.2→9.2`) that decouples the fold from birth-p27 (memory `v44-joint-reopt-28of28`). Merged as PR #2.

  A **key data-exclusion finding** came out of this: strong EZH2i de-repression (≥2.2) is *incompatible* with the honest fold + GNP+HHi + ChIP — every config recovering EZH2i≥2.2 re-inflates the fold to 6.2–6.9. The re-opt sits at EZH2i = 1.80. This *strengthens* the paper's weak-feedback story rather than being a mere shortfall.

**Phase 6 — parallel subprocess recal searches for the 16 h cycle (July 24).** `recal_16h.py` / `recal_split.py`: 320-candidate parallel searches to test dropping the 22 h period to Contestabile ~16 h, with a **hard Shh-dependence constraint** (only 5/176 candidates kept GNP arresting at SHH=0) and a **cKO post-filter** (does OFF/`f0_prc2=1` still require Shh?). See §4.

**Phase 7 — structural additions re-validated (July 25–30).** Cell-type-separated target spec (July 25); CyclinD2 as a second D-cyclin (July 29, default-neutral, 30/32); dynamic-CDKI species module baked as default (July 30, `CDKI_SPECIES=1`, **30/32** current). See [Evolution Timeline](03_evolution_timeline.md).

---

## 4. The 16 h period re-calibration (a large, deferred effort)

The single hard failure in the target table is **Period GNP**: model 22.83 h vs literature ~16 h (Nakashima 2015 live-imaging mean 15.9 h and Contestabile 2013 phase-resolved 16.25 h *agree* — the "~22–23 h Nakashima" in old docs is a **miscitation**). The 22 h is not this system's data; it is a **growth-timed default**: `Tc ≈ ln2/mu` with `mu=0.0005`, scaling exactly as `1/mu` and nearly invariant to mitogen (CyclinD1/Shh gate go/no-go at the R-point; growth sets the clock).

JP asked (July 24): abandon 22 h, let it fall to 16 h, retest. The search (`recal_16h.py`, 320 candidates) exposed a **fundamental tension**: a faster cycle → lower effective commitment threshold → three casualties simultaneously:
- **cKO (OFF, `f0_prc2=1`) Shh-independence** — OFF divides at SHH=0 (not basal-fixable; residual Gli drives un-repressed CyclinD1 over the lowered threshold).
- **Governor thresholds halve** (ON 0.18→0.09).
- **HU-S + HU-G2 folds both break** worse (§5).

Best global-`mu` result preserving ON-Shh-dependence = **25/29** (< 28/29 default). The resolution JP identified: **MB and GNP are different cell types** — the HU experiments are MB, which has a *longer* cycle than P7 GNP. The model wrongly set MB/GNP ≈ 0.99. The fix is a **CKI-coupled cell-cycle length** feature (`k_mu_cki`, default off): `mu_eff := mu/(1 + k_mu_cki*(p16+p18))`. Engaging `mu=0.000788, k_mu_cki=0.307` gives GNP Tc 16.1 h / MB Tc 22.0 h, and a focused split-search (`recal_split.py`) reached **28/29** at that split — **both HU folds passing**. But at the validation-optimal split the behavioral modifications (cKO Shh-dependence, governor) weaken; preserving them costs 1–2 checks. **Nothing is baked** — the `k_mu_cki` feature is in the builder default-off; the decision (adopt 16 h and accept a weaker governor, or keep 22 h if it reflects this system's data) is pending JP (memory `recalibration-16h-cycle-attempt`).

---

## 5. The one structural wall: the HU S-fold / G2-fold pair

The perennial lone miss is **MB HU G2 fold = 0.36 vs target ≤0.32**. Proven (2026-07-18, `_frontier_hu.py` 50-point Pareto sweep → **0/50** pass both, plus an independent adversarial-refuter agent) that **no parameter point passes both HU-S and HU-G2**.

**Why it's structural** (memory `v44-hu-s-g2-structural-limit`): in the model, the BrdU signal and the DNA-completion flux are the *same quantity* — `BrdU = vfork*aRc`, and completion rate `= kSyDna*vfork*aRc = kSyDna*BrdU`. So "counted as S (BrdU+)" ⟺ synthesis flux > 0 ⟺ the cell *is* completing replication and will reach G2. "G2 depleted" ⟺ cells persist at Dna<0.98 with flux≈0 ⟺ BrdU-negative ⟺ not counted as S. The two are **mutually exclusive by construction**. Every classifier-level decoupling relocates the failure (breaks DMSO-S, or overshoots HU-S). A true fix needs a distinct "replicating-but-stalled, BrdU+, Dna-capped" *state* — a structural change. The **decisive biology**: real HU-stalled cells incorporate *little* BrdU, which *favors* the flux definition already in the code. So forcing 28/28 here would be *less* faithful, not more. JP's ruling: accept 0.36 as a noise-consistent flow near-miss; do not hack the classifier.

---

## 6. Target philosophy: loose bands, not exact fits

A recurring theme (memory `target-flexibility-noisy-data`): the targets are noisy experimental data — bulk RNA-seq folds, flow phase fractions, a ChIP that's only "~half," Palbo "~56%." Demanding exact fits makes the model rigid and manufactures false target-vs-target tension. The bands encode measurement uncertainty deliberately. Note in particular the **CyclinD1 MB/GNP normalization ambiguity**: the target is 5.07 by the Fig 4I normalization but 7.58 in the raw RNA-seq table — a factor of 1.5, *wider than the 0.20 tolerance the model is scored against*. And the fold targets rest on small n (`MB_Ptch_het` n=2, `MB_GDC0449` n=2 vs GNP n=5), so the true error bar likely exceeds ±20%.

---

## 7. The `reset()` parameter gotcha

A calibration-infrastructure hazard that corrupted at least one result (memory `roadrunner-reset-param-gotcha`): **`roadrunner.reset()` resets floating species to initial values but does NOT reset parameters** set via `r['param'] = value`. In a sweep that calls `reset()` then sets only the swept parameter, every *other* parameter retains its value from the previous iteration — silent leakage. Symptom seen here: a sensitivity tornado (`sim_ramp_parameter_exploration.py`) whose H3K27me3 effect (0.02→0.037) contradicted the 2-D plane (up to 1.26) for the same values, because it inherited `mu` from the prior loop. **Fix:** re-apply a full DEFAULTS dict for every parameter any sweep touches, on every call. `validate_v44.run()` is safe because it sets a full baseline per condition; audit any ensemble/sweep script that assumes unset params are at default.

---

## 8. Current best, and the honest caveats

**Current state: 30/32** (dynamic-CDKI + two-cyclin default, commit 2026-07-30). The **2 misses are both HU folds** — the same structural cause (§5). The **Period GNP** target is not currently in the scored 32 in the way the identifiability plan counts it (harness still uses the 22 h check).

Two caveats surfaced by the reviewer-driven identifiability plan (`docs/PARAMETER_IDENTIFIABILITY_PLAN.md`, 2026-07-29, *plan only, nothing implemented*), important for how much weight "30/32" should carry:

1. **The pass count hides marginal passes.** Converting to band-fraction (`|model−target|/(target×tol)`), **7 of 25 continuous targets sit at ≥65% of their band** — the model is resting against ~a third of its constraints. Three of the five worst non-HU residuals (CyclinD1 MB/GNP 94%, GNP+HHi 91%, EZH2i fold 65%) are the load-bearing CyclinD1 numbers the manuscript's central claim rests on, all pulling low. The band-hinge optimizer loss is exactly zero inside the band, so the search has *no gradient* to move a target from 94% to 10% of band — some flat directions are manufactured by the objective, not the model.

2. **Verified exact parameter degeneracies.** The EZH2 → PRC2 → repression chain carries a **3-dimensional exactly-flat subspace** — PRC2-amplitude, EZH2-protein-scale, and EZH2-mRNA-scale rescalings each leave all 32 outputs unchanged to ~1e-7. Those degenerate parameters (`a0_prc2, a_rw_prc2, K_prc2, kme1–3, kTlEZ, kEZbas*`) are *precisely* the ones the module optimizers (`optimize_prc2/chain/cdh1_ezh2`) have been searching — so some overnight searches were wandering on a manifold where no evaluation can improve anything, and the landed values are one arbitrary point on a flat curve. The plan's recommendation: fix a gauge (pin `kTlEZ`, `a0_prc2`, `k_Cd_translation`), report identifiable *combinations*, run an FIM via central finite differences (h≈0.02–0.05), and reframe the deliverable as "the model's three testable predictions project onto the well-determined subspace" rather than fitted parameter values. Also flagged as inert (spurious zero eigenvalues if included): `MYCN_amplification`, `K_EZH2_repression`, `M_commit` (all dead in the default build).

This is the standing honest framing: **the model reproduces ~30 cell-type-separated behaviors within noise; the two misses are one structural wall; the period is a known unresolved gap; and individual parameter values are not the deliverable** — the calibration constrains ratios and prediction directions, not rate constants.

---

## Related pages
[Model Architecture](02_model_architecture.md) · [Evolution Timeline](03_evolution_timeline.md) · [Mechanisms Explored](04_mechanisms_explored.md) · [Things Tried & Abandoned](05_things_tried_and_abandoned.md) · [Data & Evidence](07_data_and_evidence.md) · [Reviewer & Open Questions](10_reviewer_and_open_questions.md)
