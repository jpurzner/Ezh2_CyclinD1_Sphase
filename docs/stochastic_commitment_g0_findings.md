# Stochastic commitment (p21/p18/p27): can the GNP-no-G0 / MB-transient-G0 dichotomy come from CKI-vs-CyclinD1? (2026-08-02)

**Question (JP):** rebuild the stochastic-commitment model with p21, p18, p27; test whether transient G0 can appear
in MB but NOT GNP — specifically, can **MB's increased CKI expression produce transient G0 even with MB's high
CyclinD1**? Which settings permit the dichotomy; if it doesn't work, why?

## TL;DR
- **NO — the dichotomy cannot come from the drive-brake stoichiometry, and the measured MB changes push it the
  WRONG way (MB commits *harder* → *less* G0 than GNP).** Verified across a 72,000-eval parameter hunt + 5 independent
  adversarial wirings (all confirm; max MB−GNP G0 gap = **+0.0 pp**).
- **Root cause:** the three things that differ GNP→MB (CyclinD1 ×4, CDK4/6 ×2, p18 ×3, from TMT) **all sit on the
  CDK4/6 arm**, where the drive (kCd·Cd, up 8×) dwarfs the INK4 brake (1+p18, up 2×). The CDK4/6 drive ratio
  MB/GNP is **≥ 2.0 for every possible parameter choice** (an un-removable `kCd×2` floor). The only arm that can
  bifurcate a cell into G0 is the **CDK2 arm (p21/p27)** — and the proteomics pins those **flat/shared** (p27 flat,
  p21 undetectable), so no lever there can *differentiate* MB from GNP (tightening it arrests both together).
- **What DOES work (the one route not excluded):** an **MB-specific, rare, high-birth-p21 subpopulation** (the
  Spencer/Barr mother-G2 replication-stress mechanism) — a *variance/tail* difference, not a mean-abundance one.
  It reproduces MB G0 with a bulk p21 mean that stays ~undetectable (consistent with the TMT), but it requires the
  tail to be MB-specific and the per-cell p21 high (~2.5–3.5) to beat MB's high drive.
- **Bottom line for the rebuild:** don't source MB G0 from the CyclinD1/CDK4-6/INK4 abundances. Source it from
  **inherited CDK2-arm birth heterogeneity** (an MB-specific birth-p21 tail), OR reframe the MB phenotype as
  **graded G1 lengthening** (which the baked model already does: MB G1 21h vs GNP 11h, G0-frac 0%), OR from the
  separate **p57/OLIG2 lineage-dormancy compartment** (wiki 11 Module B). p18/p27 do not contribute.

Model + all scripts: `simulations/stoch_commit_g0/` (reduced numpy ensemble; `model.py`, `analysis.py`,
`robustness.py`, `constructive.py`, + the adversarial `refute_*.py`).

## The model
Reduced stochastic-commitment ensemble (Yao Rb-E2F bistable switch + Overton p21/p27⊣CDK2 double-negative +
INK4/CIP brake on a saturating CyclinD-CDK4/6 drive):
```
Dd  = kCd*Cd / (Kd*(1 + w18*p18) + Cd)                 # CyclinD-CDK4/6 (mono-Rb); INK4(p18) brake
inh = 1/(1 + (p21+p27)/Kp)                              # CIP/KIP inhibition of the CyclinE-CDK2 (E2F self-act) arm
dE  = kb + km*(Dd + Erb*E^n/(Krb^n+E^n))*inh - kdE*E    # E2F switch: mono-Rb drive + CDK2 self-activation
dp21= s21 - (d0 + dS*E^md/(Kd2^md+E^md))*p21            # Skp2 (E2F target) degrades p21 ULTRASENSITIVELY in E
dp27= s27 - (d0 + dS*E^md/(Kd2^md+E^md))*p27            #   -> p21/p27 MAINTAINED at low E => genuine bistability
```
Each cell is born post-mitosis (E low) with p21/p27 drawn from a birth distribution (inherited from mother G2 =
the noise source). Fate after settling: E>thr = committed (CDK2-inc); E<thr = transient G0 (CDK2-low). The G0
**fraction** across the ensemble is the readout. Cell-type knobs (TMT-grounded, GNP=1): **MB Cd=4, kCd=2, p18=3**;
p27 flat, p21 low (shared birth). Genuinely bistable (birth-p sets the basin; verified via separatrix scan).

## Result 1 (negative): stoichiometry gives MB LESS G0
- **G0 fraction falls monotonically with CyclinD1:** 100% (Cd0.5) → 56% (Cd1) → 5.5% (Cd2) → **0.4% (Cd4)**.
  Higher CyclinD1 *suppresses* G0. MB's Cd↑ moves the wrong way.
- **Data-grounded dichotomy (shared flat birth-p, only Cd/kCd/p18 differ):** GNP ~6% G0, **MB 0%**. Inverted.
- **The separatrix is higher for MB:** GNP goes G0 above birth-p27≈0.3; MB commits at *all* birth-p27 up to 2.0.
- **p18 is on the wrong arm:** to get MB G0 at Cd4/kCd2 you need **p18 ≈ 18–22×**, not the measured 3× — because
  p18 brakes the CDK4/6 arm the 8× drive dominates.
- **Decomposition:** Cd×4 alone and kCd×2 alone each *remove* G0; only p18 adds it, but 3× is swamped in full MB.

## Adversarial verification (5 angles, all CONFIRM; 72,000 evals; max gap +0.0 pp)
1. **CdP21 buffer** (v44's CyclinD1↔p27 sequestration): *amplifies* commitment — high Cd sponges p27 off the CDK2
   arm → MB commits harder (MB 0% across 45 combos). The v44 feature most likely to help arrest works against it.
2. **CDK4/6-low exit** (Yang-2020, p21-independent): the only angle that cleared +5pp on paper, but it's a
   **noise artifact** — a multiplicative-CV + `max(0,·)` zero-floor pile-up in a knife-edge INK4 regime; it
   **reverses sign** under additive noise (−5 to −6pp) and reverses hard (−13 to −32pp) once coupled to the real
   E2F latch. MB's *median* CDK4/6 activity is still ≥ GNP's (commits harder; variance broadening around a higher mean).
3. **Free-parameter hunt** (72,000 evals over 12 knobs × birth-p × threshold): max MB−GNP gap **+0.0 pp**; the
   MB/GNP CDK4/6-drive ratio has a hard **minimum of 2.015** across the entire (Kd, w18) plane.
4. **Drive/brake topology** (INK4→CDK2 displacement, hard saturation, growth gate, non-monotonic Cd): no
   defensible rewire works — `kCd×2` is an un-removable drive floor; any brake big enough to bite arrests GNP first.
5. **v44 structural skeptic** (two-step Rb, CdP21-counts-in-drive, assembly-factor cap): all keep MB committing
   harder; the one theoretical escape (assembly cap making Cd inert) is foreclosed by the shared flat-p27 anchor
   *and* contradicts v44, where CyclinD1 is *the* governor.

## Result 2 (positive): an MB-specific birth-p21 tail reproduces the dichotomy
An MB-specific rare high-birth-p21 tail (Spencer/Barr; replication-stress → mother-G2 p21 → daughter CDK2-low):
- GNP (no tail) → ~0–6% G0; **MB (20–30% of cells at birth-p21 ≈ 2.5–3.5) → 10–26% G0.**
- Bulk p21 mean stays **~undetectable** even with the tail → consistent with the TMT ("p21 undetectable").
- **Requires the tail to be MB-specific** (the counterfactual: GNP *with* the same tail would have MORE G0 than MB,
  because GNP's lower drive needs less p21) — i.e. driven by MB's oncogenic replication stress, absent in GNP.
- Requires per-cell p21 ≈ 2.5–3.5 to beat MB's high drive (a stringent, testable requirement).

## Why (the deep structural reason)
The dichotomy needs a lever that **differs between GNP and MB** *and* acts on the **arm that can bifurcate** (CDK2).
The only measured differences (Cd, kCd, p18) act on the CDK4/6 arm and all deepen commitment in MB. The CDK2-arm
inhibitors (p21/p27) *can* bifurcate but are flat/shared → can't differentiate. So the only way to make MB-specific
G0 is to put an MB-specific difference *on the CDK2 arm* — which, given flat mean p27 and undetectable p21, must be
a **birth-state variance/tail** (inherited pulses), not a mean abundance. This is exactly what the earlier
deterministic work already implied (`cdki-baked-g0-g1-behavior`: "discrete G0 needs population heterogeneity, not a
clearance change").

## Recommendation for the rebuild
1. **Keep the drive-brake commitment module as the GNP engine + bulk-MB proliferation** — it correctly gives GNP
   commit-above-threshold/no-G0 and correctly makes bulk MB commit harder (matches the proteomics + biology).
2. **Add MB transient-G0 as a SEPARATE population layer**, sourced from **inherited CDK2-arm birth heterogeneity**:
   an MB-specific birth-p21 (and/or p27) tail whose *variance* is set by MB's oncogenic replication stress, not by
   mean abundance. This is the v45-style stochastic ensemble with birth-**p21** as the noise variable (not birth-p27).
3. **Specify the noise model explicitly** — the cdk46-low artifact shows results flip between multiplicative+floor
   and additive noise; validate any claimed MB G0 against the real E2F latch, not a static threshold.
4. **Consider the alternative framings**: MB's extra cycle length as **graded G1 lengthening** (baked model already
   does this, G0-frac 0%), and/or the therapy-resistant deep quiescence as the **p57/OLIG2 lineage-dormancy
   compartment** (wiki 11 Module B) — a rare, transcriptionally-imposed state, not a stoichiometric one.

## What would confirm / what's needed (JP)
- **Single-cell** CDK2-activity (sensor) or scRNA in MB vs GNP: does MB have an MB-specific high-p21 / CDK2-low
  subpopulation that GNP lacks? (Bulk cannot see it — that's the whole point.)
- Replication-stress markers (γH2AX, p53 activity) in MB vs GNP to ground the birth-p21 tail's MB-specificity.

## Caveats
- Reduced prototype (not the full v44/v45 engine) — but the adversarial v44-structural-skeptic angle confirms the
  conclusion carries to the full model's features (two-step Rb, CdP21 buffer).
- MB proteomics is n=1 (Cd/kCd/p18 folds single-sample) — but the conclusion holds across a wide fold range
  (the drive ratio ≥2 for any plausible values), so it's robust to the exact n=1 numbers.
- The refutation excludes *stoichiometric abundance shifts under a shared birth-p*; it does NOT exclude birth-state
  variance differences (the positive route) or orthogonal machinery (differentiation/EZH2 arm, p57/OLIG2).
