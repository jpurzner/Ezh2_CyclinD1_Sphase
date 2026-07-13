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
