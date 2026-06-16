# v44 EZH2 recalibration + wide parameter search — results

_June 2026. Session summary written before context ran out. Current committed model = the EZH2
mitogen-dose recalibration (commit `41485f7`); the search's best params are a CANDIDATE, not yet baked._

---

## 1. What was changed and committed (`41485f7`, `f8926b2`)

**EZH2 synthesis gained a saturating mitogen-dose term** — `EZH2_tx ×= Cd/(Kez_cd+Cd)` (Kez_cd=2,
kEZE2f 0.010→0.026). Motivation (JP): mitogen-responsiveness statements need a real dynamic range, and
EZH2 increases dose-dependently over a **wide rShh range** in cultured GNPs (Fig 4H). Diagnosis: v44's
EZH2 (=E2f×cyclins) tracked only the binary commitment (E2f saturates at commitment) → EZH2 plateaued.

**Result (now in the committed model):**
- EZH2 climbs **1.07 → 2.30 over SHH 0.5→8** (wide, monotonic, no plateau) = **Fig 4H** ✓
- **EZH2 MB/GNP = 2.05** (Fig 4J; was 1.1, a hard fail before) ✓
- The stronger MB EZH2 represses MB CyclinD1 down to **MB/GNP ≈ 4.65 — reaching the measured 5.07**
  (Fig 4I), which the prior calibration had logged as "unreachable." `validate_v44` retargeted 7.58→5.07.
- **Validation 22/27, rescue intact** (MB+HHi+EZH2i cycles). New near-miss: MB+HHi/MB 0.197 vs 0.144 —
  a *correct* consequence (EZH2 declines on cell-cycle exit → vismo's CyclinD1 drop is slightly buffered).

**Honest limit found:** the *proliferation threshold* could not be raised by hand — `kPhRbCd` is global
(raising it breaks the MB+HHi+EZH2i rescue), and raising GNP `p18` violates the measured MB/GNP p18 ratio
(3.75×). So the dynamic range lives in the now-wide **CyclinD1/EZH2 dose-responses**, not the threshold.

**Figures:** all 20 `fig_v44_*` regenerated; `fig_v44_mitogen_sensitivity` is now a 2×2 (CyclinD1 row +
the EZH2 dose-response row; GNP driven by SHH, MB tuned by vismodegib, on a Gli1 axis).

---

## 2. The overnight wide search (`simulations/v44_wide_search.py`)

17 parameters searched against a focused objective (11 soft targets: CyclinD1/EZH2/MYCN/Gli folds +
the wide-rShh EZH2 dose-response) under 8 hard figure-critical guards (GNP cycle/arrest, GNP+HHi arrest,
MB cycle, MB+HHi arrest, the rescue, CDK4/6i no-rescue). Hybrid local+wide search. Seeded from the
current params (score 899 = 8 guards + 10/11 soft).

**Seed 1: 3,334 evaluations, best score 909.5 = 8/8 guards + 11/11 soft (the objective ceiling).**
(Seed 2 failed to launch — only seed 1 produced results. Best saved in `v44_wide_search_s1_best.json`,
all evals in `v44_wide_search_s1.jsonl`.)

**Key finding: the search recovered the threshold dynamic range that I could not get by hand.** Its best
set uses **kPhRbCd = 0.35 (raised commitment threshold) WITH the rescue intact** — it co-adjusted the CKI
tones and CyclinD1 drive so the higher threshold doesn't break MB+HHi+EZH2i. It also **fixes MB+HHi/MB
(0.143, target 0.144)** — the near-miss in the committed model.

### Best parameter set (seed 1)
| param | value | param | value |
|---|---|---|---|
| kPhRbCd | **0.35** (was 0.5; raises threshold) | Kez_cd | 2.94 |
| K_CdRb | 0.319 | kEZE2f | 0.0220 |
| p18 (GNP) | 0.464 | kEZbas | 0.00027 |
| k_Cd_tx_basal | 0.483 | kDeEZ | 0.00015 |
| k_Cd_tx_Gli_max | 59.2 | k_Cd_translation | 0.801 |
| k_Cd_tx_MYCN | 21.66 | P16_MB | 0.306 |
| K_EZH2_repression | 0.539 | P18_MB | 1.553 |
| wCe | 0.582 | KSYP21_MB | 0.002 |
| KmHU_fork | 0.2 | | |

---

## 3. Full `validate_v44` on the best params → **22/27** (same count, different profile)

The search optimized its 11-target subset to 11/11, but the FULL 27-target validation stays **22/27**,
because 4 targets the objective did NOT include still fail:

| target | best-params value | committed-model value | note |
|---|---|---|---|
| MYCN GNP+HHi/GNP | 0.908 (t 0.78) | 0.908 | pre-existing; cascade, not EZH2 |
| MB G2+M duration | 3.945 (t 2.5h) | ~3.4 | **regressed** — kPhRbCd=0.35 lengthened G2 |
| MB HU S fold | 4.60 (t 1.36) | 4.35 | pre-existing; HU coupling needs re-tune for short S |
| MB HU G2 fold | 0.847 (t 0.23) | 0.877 | pre-existing |
| EZH2 transcript S/G0 | 4.43 (t 1.8–2.5) | 4.18 | within-cycle EZH2 gradient too steep |

**Net:** the search **fixes MB+HHi/MB and achieves the raised-threshold dynamic range with the rescue
intact** (the headline win), but **trades MB G2+M duration**, and does **not** fix the pre-existing
HU / EZH2-S-gradient / MYCN-GNP+HHi fails (they weren't scored). So it is not a strict improvement on the
27-target count — it's a *better dynamic range* at the cost of one cell-cycle-duration target.

---

## 4. Recommended next steps (when resuming)

1. **Decide on the dynamic-range tradeoff.** The best params give the wider GNP sub-threshold range JP
   asked for (kPhRbCd=0.35) with the rescue preserved, but cost MB G2+M duration (~3.9 h vs 2.5 h). If
   the dynamic range matters more, bake the best params (`build_model_v44(params=…)` + MB-CKI constants
   in validate/figures) and re-tune kSyCb/the G2 timer to recover the ~2.5 h G2.
2. **Add the missing targets to the search objective** (MB G2+M duration, MB HU S/G2 folds, EZH2 S/G0,
   MYCN GNP+HHi) and re-run — then a high-score set is a *strict* 27-target improvement, not a subset win.
   Also re-enable seed 2 (the parallel launch failed; run seeds sequentially or fix the `&`/cwd issue).
3. **EZH2 transcript S/G0 (4.4 vs 1.8–2.5)** is the most persistent miss — the within-cycle EZH2 gradient
   is too steep; the mitogen-dose term slightly worsened it. Worth a targeted tune of the cyclE/cyclA gate
   (wCe, K_Ce_EZ, K_Ca_EZ) + kDeEZ.
4. The **HU folds** need the fork-coupling (KmHU_fork/vmin_fork) re-tuned for the current short S-phase
   (flagged open since the S-phase recalibration).

**Files:** best params `simulations/v44_wide_search_s1_best.json`; all evals `…_s1.jsonl`; search script
`simulations/v44_wide_search.py`; committed recalibration `41485f7`; figures `f8926b2`. The committed
model is the recalibration (without the search params) — the search best is a documented candidate.
