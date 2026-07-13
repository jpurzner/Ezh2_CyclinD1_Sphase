# Daughter-transfer mechanism for v44: CDK inhibitors + other key species

*Design plan, 2026-07-13. Produced by a judge-panel design workflow (5 independent proposals →
adversarial scoring → synthesis), then hand-audited. Builds directly on the mother-G2 negative
result ([memory: v44-mother-g2-negative-result]).*

## Why

Today every species is transferred at division by its own ad-hoc rule (p27 resets toward 0.6 with a
commitment carryover; Rb carries a committed fraction; EZH2 halves; the mark halves at S; INK4 is a
fixed constant; cyclins hard-reset). There is no unifying logic, and the one principled attempt — the
`with_mother_g2` prototype — **failed**: `birth p27 = P21_div + g·Sg2` can only *raise* p27 above the
0.6 floor, so high-mitogen daughters can never become "immediate," and the Overton graded split
collapsed to all-quiescent (plus a cycle-frequency confound from the reset-integral form).

This plan replaces the pile of per-species rules with **five transfer modes**, routes the cell's entire
inherited fate through **one species (p27)** whose birth value is *set* (not added) by a low-pass
mother-mitogen tracker with a **low floor**, and derives the Overton graded→binary split + Spencer ~98%
sister concordance from a clean two-layer heterogeneity source — all without touching the 27/28
validation or the hard 5.07× MB/GNP CyclinD1 fold.

## Convention: amounts vs concentrations (foundational — fixes the EZH2 anomaly)

Every species is an **amount** N (halves at symmetric division; set by synth−deg between divisions) or a
**concentration** C=N/V (**unchanged** at division; diluted by growth −μ·C between divisions). The taxonomy
below assumes ONE consistent convention:

- **Concentration** (all soluble regulators — cyclins, p27, Skp2, **EZH2/EZH2m**): **not halved at division**;
  reset only by their own kinetics/biological events (APC for cyclins, setpoint for p27, turnover for EZH2).
- **Amount / size-titration** (`Rb` family, `mass`): the *only* species that legitimately halve/conserve as
  amounts, with `/mass` readout — this IS the Zatulovskiy Rb-dilution timer.
- **Chromatin occupancy** (the mark `Mk, m1_me, m2_me`): a fraction in [0,1], halved at **S** (`Dna>0.05`) by
  DNA replication (parental modified histones spread over 2× DNA + unmodified new histones). Independent of cell
  volume; NOT a division-partition effect.

> **Correction to current code:** EZH2/EZH2m **and E1/Emi1** (`E1 = E1/2`) are currently halved at division
> ([L416](../src/build_model_v44_heldt.py) / [L80](../src/build_model_v44_heldt.py)) — the soluble proteins so
> treated (EZH2 the one that matters — E1 turns over fast). Halving is inconsistent here (halved like an amount, read like a concentration via
> fixed-K Hills with no `/mass`, and no growth-dilution term). It "works" only because EZH2 (t½≈9.9 h) is the sole
> soluble protein slow enough vs the ~22 h cycle for the halving to leave a mark; cyclins re-equilibrate in minutes
> so halving them is a no-op. Under the concentration convention EZH2 should **not** be halved: its ~10 h turnover
> (optionally + explicit `−μ·EZH2`) does the reset, and the S-duration integration lives in the synthesis gating.
> Requires a bounded re-fit of `kDeEZ`/`kTlEZ` (co-fit with the halving); makes EZH2 consistent AND gives a more
> realistic decaying mother→daughter EZH2 carryover.

## The five transfer modes

| Mode | Definition | Species |
|---|---|---|
| **1. setpoint-from-mother-history** (fast, p27) | Birth value **set** onto a low-floored range by a never-reset low-pass tracker of the mother's R-window mitogen, re-referenced by a *small* carryover `fc_p27` decoupled from the pRb-marker `f_commit_carry`. Sole carrier of inherited fate. | `P21` (=p27) |
| **2. setpoint-from-mother-stress** (medium, p21) | Minor additive `+w21·Str` term inside the p27 birth rule, set by a medium-timescale low-pass of the damage axis. **Zero at the unstressed baseline** (`Dam=0`) → validation-neutral. | p21/stress arm |
| **3. conserved-partition + commitment-carryover** (fast) | Pool total invariant; amounts redistribute among phospho/complex sub-forms with `f_commit_carry` retained; concentration falls only via `mass/2`. `tRb`, `tE2f` exactly conserved. | `Rb, pRb, Rbm, RbE2f, RbmE2f, E1` |
| **4a. size-titration** (amount) | Halves/conserved as an AMOUNT at mitosis; readout via `/mass`. The Zatulovskiy Rb-dilution timer. | `Rb` family, `mass` |
| **4b. chromatin-occupancy dilution** | Occupancy fraction halved at **S** (`Dna>0.05`) by DNA replication; re-established by weak read-write. Not a volume effect; **not** a fate carrier. | `Mk, m1_me, m2_me` |
| **5. hard-reset / mitogen-clamp / concentration-carried** | Set to a constant independent of mother state (APC-degraded machinery → 0; re-licensed origins; identity constants), **or** never assigned and left to re-equilibrate by kinetics — for fast species within minutes (`Cd`, cyclins), for slow **EZH2** a decaying mother→daughter carryover via its ~10 h turnover. **No division halving on any concentration.** | `MPF, preMPF, Cdc20, Ca, Ce, Skp2, CeP21, CaP21, Dna, Rc/pRc/aRc/iRc, p16, p18, Cd, EZH2, EZH2m` |

Overlaid on the CKI family is a **signal lens** (NOT a pure timescale clock): **p27** = fast, mitogen-gated
(per-cycle) — the sole CKI in GNP, where fate is p27-only; **p21** = medium, replication-stress-gated;
**p16/p18** = **MB-specific identity CKIs** — substantially present in MB, ~absent in GNP. They are cell-type
constants, **not a differentiation clock** (that framing was wrong). ⚠ The code currently has GNP `p18=0.464`
(non-zero, cited as a constitutive Shh-maintained Cdkn2c) — reconcile with "GNP is all p27" (open question).

## Species-by-species

| Species | Now | Proposed | Mode |
|---|---|---|---|
| **p27 (`P21`)** | `P21 = P21_div − fc·(P21_div − P21)`, floor 0.6 | `P21 = (1−fc_p27)·(P27set + w21·Str) + fc_p27·P21` | 1 |
| **p21 (stress)** | folded into pool | enters only as `+w21·Str`; `w21≈0.25` | 2 |
| **p16 / p18** | identity params (p16=0 GNP / 0.306 MB; p18=0.464 GNP / 1.553 MB) | **unchanged as constants** — MB-specific identity CKIs (present in MB, GNP fate is p27-only), **NOT a differentiation clock**. (Flag: is GNP p18=0.464 real or should it drop toward 0?) | 5 (identity constant) |
| **CyclinD1 (`Cd`)** | not in `E_div`, mitogen-clamped | **unchanged** — fate rides on p27, not daughter `Cd`; so `Cd` is fate-independent across sisters (Spencer) and the 5.07 fold is structurally protected | 5 |
| **Rb / pRb / Rbm / RbE2f / RbmE2f** | two-step carryover | **unchanged** (keeps Moser pRb-high 12–20 min post-anaphase) | 3 |
| **EZH2 / EZH2m** | halved at mitosis (L416) | **CHANGE**: drop the halving — EZH2 is a concentration, preserved at division, reset by its ~10 h turnover (± `−μ·EZH2`). Re-fit `kDeEZ`/`kTlEZ`. | 5 (concentration-carried) |
| **mark (Mk, m1_me, m2_me)** | halved at early-S `Dna>0.05` | **unchanged** — chromatin-occupancy dilution by replication; fate memory must NOT live in a self-reinforcing mark (data-excluded) | 4b |
| **Skp2** | `= 0.05` | **unchanged** — the Skp2-p27-E2F toggle IS the binarizer that turns graded birth-p27 into binary CDK2low/inc | 5 |
| **Ce / Ca / E1 / MPF / preMPF / Cdc20 / Dna / origins / CeP21 / CaP21 / mass** | resets/halving | **unchanged** | 3/4/5 |

**The only substantive change to `E_div` is the `P21` line.** Everything else is a taxonomy relabel of
rules that already exist.

## New machinery (antimony)

```
// ---- Mitogen-DEFICIT low-pass tracker (frequency-clean; NOT reset at division; shared by sisters) ----
species Mtr in Cell;  Mtr = 0.5;
k_M = 0.0025; K_m = 0.90; n_m = 3;          // 1/k_M ≈ 6.7 h ≈ mother R-window
mdef := K_m^n_m/(K_m^n_m + Cd^n_m);          // deficit: HIGH at low Cd (low mitogen), LOW at high Cd
Mtr_track: => Mtr;  Cell*k_M*mdef;
Mtr_relax: Mtr => ; Cell*k_M*Mtr;            // net dMtr/dt = k_M*(mdef − Mtr)  → level tracker, no frequency confound

// ---- Low-floored birth-p27 setpoint: SET, not ADD  (CORRECTED SIGN — see note) ----
P21_lo = 0.06; P21_hi = 1.10; K_S = 0.60; h_S = 4; fc_p27 = 0.10; w21 = 0.25;
P27set := P21_lo + (P21_hi − P21_lo)*Mtr^h_S/(K_S^h_S + Mtr^h_S);   // INCREASING in deficit Mtr

// ---- Medium-timescale replication-STRESS/p21 tracker (exactly 0 at HU=0 → validation-neutral) ----
// NB: v44 has NO damage/p53 axis, and Chk1=aRc/(jChk+aRc) saturates ≈1 in EVERY normal S (jChk=0.03),
// so neither Dam nor Chk1/aRc is a clean stress signal. STALLED FORKS aRc*(1-vfork) are: vfork≡1 at HU=0
// (so the term is identically 0 in every non-HU validation run), and >0 only under HU/replication stress.
species Str in Cell;  Str = 0.0;
k_str = 0.0012;                               // ~14 h relaxation (medium timescale)
sstress := aRc*(1 - vfork);                   // stalled forks: 0 at HU=0, >0 under HU → daughter p21
Str_track: => Str;  Cell*k_str*sstress;
Str_relax: Str => ; Cell*k_str*Str;           // dStr/dt = k_str*(sstress − Str)
```

Proposed `E_div` (only the `P21` assignment differs from current L80; two-step Rb block already default):

```
E_div: at (MPF > MPF_div): Dna = 0, Rc = 1, pRc = 0, aRc = 0, iRc = 0,
  Rb = Rb + (1 - f_commit_carry)*(pRb + Rbm), pRb = f_commit_carry*pRb, Rbm = f_commit_carry*Rbm,
  RbE2f = RbE2f + (1 - f_commit_carry)*RbmE2f, RbmE2f = f_commit_carry*RbmE2f,
  P21 = (1 - fc_p27)*(P27set + w21*Str) + fc_p27*P21,        // <-- the change (was: P21 = P21_div - fc*(P21_div - P21))
  CeP21 = 0, CaP21 = 0, Skp2 = 0.05, Ce = Ce_div + f_commit_carry*Ce, Ca = Ca_div,
  E1 = E1/2, MPF = 0, preMPF = 0, Cdc20 = 0 ;
// Mtr and Str are DELIBERATELY absent from E_div → carried across division (shared by both sisters, frequency-clean).
```

> **⚠ Sign correction (hand-audit).** The workflow's synthesis wrote `P27set` with `K_S^h/(K_S^h+Mtr^h)`
> (decreasing in `Mtr`). Since `Mtr` tracks the *deficit* (high at low mitogen), that is **inverted** —
> it gives *low* birth-p27 at *low* mitogen, the exact failure that sank the mother-G2 prototype.
> Verified numerically: as-written, low-mitogen (`Mtr=1.2`) → p27=0.12, high-mitogen (`Mtr=0.1`) → p27=1.10.
> The corrected `Mtr^h/(K_S^h+Mtr^h)` gives low-mitogen → 1.04, high-mitogen → 0.06. ✔ Use the corrected form
> above. (Equivalent alternative: track the mitogen *level* `mlev := Cd^n/(K_m^n+Cd^n)` and keep the
> decreasing sigmoid — pick whichever naming you find clearer; the corrected version keeps the deficit story.)

## Why the two failure modes are fixed

- **Birth-p27 floor** (killed the prototype): the setpoint is now **SET** onto `[P21_lo=0.06, P21_hi=1.1]`,
  and — critically — `fc_p27` (≈0.10) is **decoupled from** the pRb-marker `f_commit_carry` (≈0.375). With
  `fc_p27` small, the effective floor ≈ `P21_lo`, so high-mitogen daughters actually reach the immediate /
  CDK2-inc basin. (The Moser carryover constraint applies to *Rb*, nothing forces it on p27.) A judge flagged
  that the naive "0.08 floor" is *not* the effective floor unless this decoupling is done — it is.
- **Frequency confound**: `Mtr` is a **relaxation level-tracker** (steady state = the signal), not a
  reset-integral over G2, so its value doesn't scale with how often the cell divides.

## How Overton and Spencer both hold

Two decoupled layers:

- **Upstream (mother→mother, the Overton knob):** a per-lineage lognormal mitogen-sensitivity multiplier
  `s_mit ~ LogNormal(1, CV=0.70)` on the Gli/Shh→`Cd` drive — the **CV grounded in the measured CyclinD1
  G0/1 CV 0.70**, not a placeholder. This gives the birth-p27 distribution real width at fixed mean mitogen;
  the retained Skp2-p27-E2F toggle binarizes each daughter, and across N≈60 the fraction above the separatrix
  sweeps smoothly 0→1 as mean mitogen crosses the window. Because the spread is **wide** (CV 0.70), MB's lower
  mean birth-p27 still has an upper tail crossing threshold → MB ≈19% while GNP ≈32% **on the same map**
  (the graded curve is an ensemble property, never a deterministic mean — this resolves the "single monotonic
  map forces MB→0%" objection).
- **Sister→sister (Spencer ~98%):** both daughters inherit the **same** never-reset `Mtr` and the same
  `s_mit` (a mother property fixed pre-anaphase), so `P27set` is shared; the only sister difference is symmetric
  partition noise capped at `CV ≤ 2%`. Concordance = 1 − P(`P27set` within ~2·`eps_part` of the separatrix) —
  check this at **mid-window** mitogen (cells clustered near `K_S`), not just the extremes.

## Cell-type variability of the transient-G0 fraction (the central requirement)

Moser 2018 (Fig. S1) measured the CDK2low/transient-G0 fraction across cell lines: **1% (HCT116) → 31%
(RPE-hTERT)**, with MCF10A 26%, HLF 23%, MCF7 18%, U2OS 17%. Two constraints follow:

1. **The fraction is NOT ordered by cancer-vs-normal or by any single mitogen law** — HCT116 (1%) and MCF7
   (18%) are both cancer. The **decision mechanism is universal; only the fraction is cell-type-specific.**
2. **`pRb`-hyper fraction right after anaphase ≈ CDK2inc fraction**; the p21-high fraction ≈ CDK2low; and
   **both the hypo-Rb and high-p21 birth populations decay with time-since-anaphase** (the transient signature).

**Design consequence:** the transient-G0 fraction must be an **emergent, cell-type-specific** quantity — never
forced by a universal mitogen→fate map. It is set per cell type by the *combination* of {mean mitogen
(Gli/MYCN drive), CKI complement (p27 only in GNP; p27 + MB-specific p16/p18 in MB), commitment threshold
`M_commit`, and the birth-p27 setpoint `K_S`}. This is exactly what resolves the panel's "a single monotonic map
pins GNP and forces MB→0%" objection — with cell-type-specific parameters, GNP≈32% and MB≈19% coexist naturally
(both inside the observed 1–31% band), GNP > MB emerges from MB's higher mitogen (fewer CDK2low), and **no
"cancer → low CDK2low" rule is hard-coded.** Any of {`K_S`, mean mitogen, the CKI complement, the `s_mit` spread}
may be tuned per cell type to hit its measured fraction.

**Two new validation targets** (both directly checkable via the two-step `pRb`-hyper birth marker):
- birth `pRb`-hyper fraction **=** CDK2inc (immediate) fraction, per cell type;
- the hypo-Rb and high-p27 birth populations **decay monotonically with time-since-anaphase** (transience).

## Calibration & validation — this is a fixed-point problem, not a no-op

The most-cited flaw across the panel: making birth-p27 a function of `Cd` is **not** a validation no-op,
because there is a closed loop `p27 → CDK4/6+CDK2 → Rb → E2f → EZH2 (E2f target) → PRC2 → represses Cd_mRNA →
Cd → Mtr → p27`. So "GNP birth-p27 = 0.6" becomes a **fixed-point + stability** condition, not a preserved
constant. Plan:

1. **Single-lineage re-pin** (`validate_v44.py`): jointly fit `{P21_lo, P21_hi, K_S, h_S, K_m, n_m, k_M, fc_p27}`
   so the deterministic GNP settled `Mtr → P27set` self-consistently lands ≈0.6 *through the loop*, holding
   two-step Rb / `fc` / Skp2 / EZH2 / PRC2 / `Cd` params byte-identical. Target 27/28.
2. **Arrest guard** (two regimes — keep them distinct): **steady-state** arrest conditions (a cell *born* into
   GNP−SHH / GNP+HHi / serum-starve / MB+CDK4/6i) must still give **0 divisions** — this is hard and unchanged,
   because a low-mitogen cell settles to high `Mtr` → high birth-p27. **Acute withdrawal** of an already-cycling
   cell may coast **one** extra division while the ~6.7 h `Mtr` catches up — **JP: that's fine**, so `k_M` stays
   at 6.7 h. Only verify the steady-state targets stay at 0; don't over-tighten `k_M` to kill the acute coast.
3. **5.07-fold guard**: confirm effective sensitivity stays in the weak regime (`S ≈ 0.81`) and the fold stays
   in [4.75, 5.45] — protected because `Mtr` only **reads** `Cd` (no new feedback onto CyclinD1 transcription)
   and no chromatin param moves.
4. **Ensemble** (`sim_g0_bifurcation.py`, extend): sweep `s_mit ~ LogNormal(1, CV=0.70)`, N≈60, score the
   **per-daughter/birth-state** CDK2low fraction (not the period-averaged definition that inflates GNP~57%/MB~94%);
   target GNP ≈32% / MB ≈19% and a smooth graded curve. Fork sister pairs with independent `eps_part (CV ≤ 2%)`,
   confirm ≥98% same-fate **and** fate-independent daughter `Cd`.

## Phased rollout (each opt-in, shippable)

0. **Byte-identical refactor** — rebuild the `E_div` string from a taxonomy-labeled assembly (five explicit
   mode groups), no behavior change; add `with_mitogen_tracker` / `with_stress_p21` flags (default OFF).
   `validate_v44.py` must reproduce 27/28 exactly. Pure cleanup.
0b. **EZH2 convention fix** — ✅ **DONE (opt-in, `with_ezh2_conc`, default OFF).** Dropped `EZH2=EZH2/2, EZH2m=EZH2m/2`
   from the mitosis event; EZH2 reset now comes from its own turnover (no `−μ` term needed). Dropping the halving
   inflated EZH2 ~47% → re-fit `kTlEZ ×0.78` (`kDeEZ` unchanged) restores the level. **Result: 26/28.** All EZH2
   phase ratios + levels hold (G0/cycling 0.372, G2/G0 0.998, EZH2i fold 2.466, MB/GNP 4.865). The one new miss is
   the **CyclinD1 GNP+HHi/GNP** fold (0.29 vs a noisy 0.16 target) — it sat at its passing edge (0.22) even *with*
   the halving, and trades off against EZH2 G0/cycling in `(kTlEZ,kDeEZ)` (2-param structural ceiling). E1/Emi1
   halving KEPT (reclassified as real SCF-βTrCP mitotic destruction, not partitioning). Harness: `refit_ezh2_conc.py`.
   **PROMOTED to default (JP 2026-07-13):** `with_ezh2_conc=True` is now the builder default and the validator's
   default; `with_ezh2_conc=False` recovers the legacy halving (27/28). The lost target is accepted as a noisy HHi
   fold per the loosen-noisy-targets policy. Figure suite regenerated on the new default.
1. **`with_mitogen_tracker`** — ✅ **BUILT (opt-in, default OFF).** `Mtr` deficit tracker + `P27set` (corrected
   INCREASING-in-`Mtr` sign) + `fc_p27=0.10`; `E_div` P21 line swapped to `(1−fc_p27)·P27set + fc_p27·P21`.
   Key finding: GNP's *settled* Cd is ~1.9 (not the 0.70 init the plan assumed), so `K_m` was re-placed 0.90→**2.5**
   to center GNP on the sigmoid shoulder. Result: **GNP P27set≈0.6** (transient-G0) vs **MB=0.06** (proliferative) —
   the mitogen-dependent birth-p27 mechanism works, GNP > MB by construction, both cycling. **Validation holds 26/28**
   and the stubborn CyclinD1 GNP+HHi/GNP fold *improved* 0.29→0.242 (Phase 1 helps, doesn't hurt). Stays opt-in
   until the Phase-2 ensemble validates Overton + Spencer. Env hook `MITOGEN_TRACKER`.
2. **Ensemble enablement** — wire `s_mit` into `sim_g0_bifurcation.py`; validate the graded curve + per-daughter
   32%/19% + sister concordance ≥98%.
3. **`with_stress_p21`** (OFF) — add `Str` + reactions + `w21·Str`; verify `Dam=0` leaves all non-HU targets
   untouched.
4. **Docs-only** — record the p16/p18 slow/identity classification and the optional differentiation clock (OFF).
   Promote `with_mitogen_tracker` to default only after Phases 1–2 demonstrably hold 27/28 + 5.07 + the fractions.

## Open questions for JP

1. **Tracker driver**: `Cd` protein (default here, known scale) vs `Cd_mRNA` (phase-flatter) vs raw Gli/MYCN
   drive — which best represents the mother's late-G2 mitogen sense and is least loop-contaminated?
2. **What does CV 0.70 measure** — the mitogen *sensitivity* `s_mit` (mother→mother, sets the Overton spread) or
   only the instantaneous CyclinD1 pool? If sensitivity is narrower, MB 19% may need `K_S` re-placement instead
   of the wide-CV tail.
3. ~~First-daughter transient~~ **RESOLVED (JP):** one extra coasting division after acute withdrawal is fine →
   `k_M` stays at the ~6.7 h R-window value (no need to force it fast), which also relieves the Spencer squeeze.
4. **`fc_p27` = 0** (pure SET, cleanest floor, no p27 memory) **or ≈0.1** (a whisper of Spencer p27 carryover)?
5. Does the toggle separatrix (`M_commit` baked 2.4764) hold under the widened birth-p27 range, or budget a joint
   `M_commit` nudge in the Phase-1 re-fit?
6. ~~HU driver for `Str`~~ **RESOLVED:** no `Dam` axis exists, and `Chk1=aRc/(jChk+aRc)` saturates ≈1 in every
   normal S (`jChk=0.03`) so neither is a clean stress signal. Drive `Str` from **stalled forks `aRc·(1−vfork)`** —
   identically 0 at `HU=0` (validation-neutral), >0 only under HU/replication stress.
7. ~~p16/p18 differentiation clock~~ **RESOLVED (JP):** p16/p18 are MB-specific identity CKIs, not a differentiation
   clock; no clock. Remaining sub-question: reconcile GNP `p18=0.464` (cited constitutive Cdkn2c) with "GNP is all p27."

## Deepest residual risk

The **Spencer decorrelation squeeze** (flagged by the adversarial judge): `Mtr` is a lagged copy of the *same*
`Cd` that re-clamps daughter CyclinD1, so at fixed mitogen the fate it assigns is not fully independent of the
daughter's own `Cd`. The 5–10 h lag (mother pre-anaphase vs daughter post-anaphase mitogen) is the defense, but
whether it buys *true* fate-vs-daughter-`Cd` independence is an empirical question for the Phase-2 sister fork.
If it fails, the p21/`Str` arm (mitogen-orthogonal) is the pressure-relief valve for decorrelation.
