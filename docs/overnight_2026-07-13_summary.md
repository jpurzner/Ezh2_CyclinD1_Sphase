# Overnight run — 2026-07-13

Two tasks, per your choices (bake-but-keep-opt-in · prototype the mother-G2 integrator · local commits, no push).
**Nothing was pushed. The single-step chain remains the default working model at 28/28 (untouched).**

---

## 1. Two-step Rb — re-optimized, baked as a ready opt-in (27/28)

**Motivation (from your literature review + the Rb review):** the default single-step Rb lets CyclinD directly
release E2F, which is wrong at the R-point (Narasimha/Sanidas: CyclinD only *mono*-phosphorylates; CyclinE-CDK2
*hyper*-phosphorylates) and is the cause of the pRb-saturation artifact. The two-step path existed but was stale.

**What I did:**
- Set it up properly (commit 82cde4e): division reset now *inherits* mono/hyper Rb via `f_commit_carry` (committed
  daughters exit mitosis pRb-high, per Moser 2018) instead of full dephosphorylation; integrated the INK4/CIP brake
  (p16/p18/w_p27) so it composes with the MB conditions; fixed the `M_commit` no-op. Framed the CKI as **p27**.
- Re-optimized (commit 6f81592 optimizer; 500 random + 9 climb over 21 commitment/cell-cycle/EZH2/PRC2 params).
- **Baked branch-specifically** (commit 03707e5): applied only when `with_two_step_rb=True`, so the single-step
  default is unchanged. `build_model_v44(with_two_step_rb=True)` now validates without any extra params.

**Result — 27/28** (from 23/28 untuned): ChIP 0.58, CyclinD1 5.82, EZH2i 1.84, Palbo 0.52, MB S% 20.6, read-write
rw_frac 0.61. **Single-step default still 28/28.**
- **The payoff is confirmed:** pRb(hyper) is now a *valid* G0 marker — norm pRb **0.48 in G0 vs 0.82 committed**
  (single-step was 0.88/0.98, saturated). This matches Moser 2018 and removes the reason you had to use p27 as the
  marker.
- **The one remaining fail:** `MB HU S fold` 0.324 vs target 1.36 — the hydroxyurea replication-checkpoint
  redistribution interacts badly with the two-step commitment (HU *lowers* MB S-fraction instead of raising it).
  Specific and mechanistically interesting; everything else passes.

**To promote it to default** (if you decide it's worth it): flip `with_two_step_rb=True` in the `build_model_v44`
signature + `validate_v44`'s `TWO_STEP_RB` env default to `'1'`. Then it becomes the working model (with pRb restored
as a valid marker). I did **not** do this — you asked to keep it opt-in.

---

## 2. Mother-G2 p27 integrator — prototyped (opt-in, NOT promoted)

**Motivation:** the literature's priority-2 change — the daughter's birth p27 should be set by the *mother's G2
mitogen history* (Spencer R1 / Min 2020), making the fate a threshold on an inherited condition.

**What I did (commit ba0a595):** added `with_mother_g2` (default off). New species `Sg2` integrates a p27-inducing
signal during the mother's G2 (`Dna~1`) ∝ mitogen deficit; at division the daughter's birth p27 = `P21_div + g_moth·Sg2`,
and `Sg2` resets.

**Sanity check:** directionally correct — low mother-G2 mitogen → higher `Sg2` → more G0 (G0-fraction 0.59→0.76 as
SHH drops 1.2→0.3, monotone). This establishes the *mechanism*. It is **not calibrated** and a single deterministic
cell can't show a sharp bifurcation — reproducing the Overton graded→binary split + the Spencer 98% sister
concordance needs the heterogeneity ensemble. That's the natural next step.

---

## Decisions waiting for you

1. **Promote the two-step Rb to default?** It fixes the R-point logic and restores pRb as a valid marker, at 27/28.
   One-line flip (above). My lean: yes, it's a real structural improvement — but the MB HU S fold and the fact that
   it's a bigger change than the chain are worth your eye first.
2. **Does the MB HU S fold matter enough to chase?** It's a single perturbation-prediction target; I can dig into
   the two-step HU/checkpoint interaction, or you accept 27/28.
3. **Mother-G2 next step:** calibrate it + add the ensemble to reproduce the Overton bifurcation, and benchmark
   against your Section-5 anchors (CDK2low fraction, M-to-S CV) — as a *sanity check*, not a fit (GNP/MB ≠ those lines).

## How to use what's here
- Two-step model: `build_model_v44(with_two_step_rb=True)` or validate with `TWO_STEP_RB=1`.
- Mother-G2 integrator: `build_model_v44(with_mother_g2=True)`.
- Both are opt-in; default `build_model_v44()` is the single-step chain (28/28), unchanged.

---

## MORNING UPDATE — you said "proceed with all 3"; here is what actually happened

**1. Two-step Rb PROMOTED to default (done).** `build_model_v44()` is now the two-step model (27/28).
Single-step is opt-in via `with_two_step_rb=False` (still 28/28). The daughter inherits mono/hyper Rb
through `f_commit_carry`. Branch-specific two-step calibration is baked (`_ts_bake`), applied tolerantly
so every flag combination still builds. **Whole `fig_v44_*` suite (41 figures) audited and regenerated on
the two-step default** — 3 promotion-induced breakages fixed (bake tolerance for missing params;
`g0_mechanism` period-cap + tight-bbox canvas blow-up; `differentiation_collapse` CVODE stiffness at the
abrupt mid-cycle withdrawal → smaller step + coarser output), 4 lumped-mark studies pinned to single-step.

**2. MB HU S fold — chased to the structural ceiling (27/28 is the two-step max).** The two HU folds
(MB HU-S ↑ and MB HU-G2 ↓) trade off against each other in the two-step checkpoint; the joint optimizer
(`optimize_twostep.py`, band penalties on both) could not clear both simultaneously across two rounds.
I baked the calibration where the **central phenotype (MB HU-S fold, 1.34 vs target 1.36) passes** and the
MB HU-G2 fold is the one documented miss. Treat 27/28 as the two-step's ceiling, not a tuning gap.

**3. Mother-G2 ensemble — NEGATIVE result (honest).** Built the heterogeneity ensemble
(`fig_v44_mother_g2_ensemble.py`, N=60 cells, SHH sweep) on the two-step base. It does **NOT** reproduce
the Overton graded→binary bifurcation. With the integrator ON the quiescent fraction **saturates ≈1.0 at
all mitogen levels** (and is slightly inverted); OFF it is 0. Two design flaws the ensemble exposed:
  - **Birth-p27 floor:** `birth p27 = P21_div + g_moth·Sg2` can only *raise* p27 above the `P21_div=0.6`
    floor, so high-mitogen daughters can never fall into the immediate/CDK2-inc basin — there is no
    low-quiescent tail to make a graded curve.
  - **Frequency confound:** `Sg2` integrates over G2 and is confounded by cycle frequency, competing with
    the mitogen-deficit signal.
  - **Fix for a future session** (identified, not implemented): birth p27 should be *set by* (not added to)
    a **low-pass mitogen tracker** with a *low* floor — `dSg2/dt = k·(deficit − Sg2)`, birth p27 spanning
    immediate↔quiescent as mitogen varies. The ensemble scaffold + diagnostic figure are kept for that work.

**Net:** #1 landed as a real structural improvement; #2 hit a documented structural ceiling; #3 is a clean
negative result that tells us exactly what to redesign. Everything is committed locally (no push, per your call).
