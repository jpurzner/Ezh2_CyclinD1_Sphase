# CDKI module — Phase A verification + option B (p27-inhibitory) build (2026-08-05, autonomous session)

**Context.** JP: "start with A and then move on." Implementing the dynamic-CDKI design spec v2.5 (§9a), with JP's
direction that **p27 drives G0** and the **CyclinD1/p27 ratio tunes GNP G1**. Standing rule respected: **nothing baked
into the default** — all new work is behind `with_p27_optionB` (default OFF); the shipped default is unchanged at 30/32.

## 1. Phase A — already done and baked (the spec's §5 "NOT YET APPLIED" is STALE)

The code as-shipped (dynamic-CDKI default, `ec1a141`-lineage, 30/32) already implements all of Phase A. Verified each
stop-check against the running model (`simulations/phaseA_verify.py`):

| Phase A item | spec status | actual code | verified |
|---|---|---|---|
| A1 KPC free-pool + magnitude | not applied | `kDeKPC` on free P21 only; CdP21 carries basal (`kDeP21`) | ✅ P21+CdP21 nonzero (GNP 0.73, MB 0.42) — **the P21≈0 bug is fixed** |
| A2 CdP21 stoichiometric+saturating | not applied | `Cd + P21 -> CdP21` (Cd consumed), bimolecular | ✅ CdP21 sublinear in the cycling regime |
| A3 P21 out of Rb denominator | not applied | denominator = `K_CdRb*(1 + p18_prot + p19_prot)` (INK4-only) | ✅ pRb monotonic in p27 (no non-monotonic dip) |
| A4 division mass conservation | not applied | `CdP21 => Cd` (no loss); CdP21 not reset at division | ✅ tP21 drift +3.1%, tp21a +0.5% (no sink) |
| A5 PIP-degron default | not applied | `with_p21_pip_degron=True` default | ✅ 30/32 |
| A6 release products (4 complexes) | verify | `Cep21a=>p21a`/`=>Ce` etc. present | ✅ |

**Note (kDeKPC magnitude):** baked `kDeKPC=0.02` gives free-p27 t½ ≈ 30 min, faster than the spec's 1–2 h target
(0.006) — the calibration pulled it up. Not a bug; worth revisiting if p27 stability becomes load-bearing.

**Key finding from A1 — why option B matters.** Under option A the p27 pool is ~99% **buffered** (CdP21), and the
buffered trimer is treated as **active** (it sits in the Rb-drive numerator). So free p27 is ~0.005–0.008 and
**p27 is nearly inert as a brake** — it only acts through the tiny free pool on CDK2. To "use p27 to drive G0 / tune
G1" (JP's goal), the buffered pool has to *do something* → **option B**.

## 2. Option B (Fan-Meyer p27-inhibitory CDK4/6) — implemented behind a flag

`build_model_v44(with_p27_optionB=True)` (default False). One structural change: the buffered p27·CyclinD-CDK4/6
complex (CdP21) is **removed from the Rb-drive** (it's inactive). The buffer reaction still **consumes Cd**, so raising
p27 now sequesters CyclinD1 out of the active pool → **p27 inhibits CDK4/6 via the CyclinD1/p27 ratio** (Fan-Meyer
2021), on top of free p27 → CDK2. Harness flag: `P27_OPTIONB=1`.

**Structural probe (`optionB_probe.py`):** on the *unchanged* option-A params, option B makes **GNP arrest (1 div)
while MB keeps cycling (6 div)** — the drive-asymmetry (my superset theorem): p27-as-a-brake hits the lower-drive GNP
harder. So option B does **not** hand MB-specific G0 for free (MB's high CyclinD1 overwhelms the brake). It needs a
re-tune of the drive/brake balance, and MB-specific G0 must come from the reversibility/population layer — consistent
with the earlier stochastic finding.

## 3. Re-calibration — option B reaches 31/32 with a trivial coarse tune (beats option A)

Coarse regime search (`optionB_recal_coarse.py`) found 9 feasible combos where GNP cycles ~17.5 h, MB ~24.4 h, and
Shh-dependence is preserved (GNP−SHH / +HHi / serum all arrest). Raising drive (`kPhRbCd` 0.14→0.32) offsets p27's new
brake; moderate `kSeqCd` (0.3–0.8).

**Full 32-target harness under option B (two independent combos):**
- `{kSeqCd 0.8, kSyP21 0.0006, kPhRbCd 0.32}` → **31/32** (CyclinD1 MB/GNP 4.63, period 17.6 h)
- `{kSeqCd 0.3, kSyP21 0.00116, kPhRbCd 0.32}` → **31/32** (CyclinD1 MB/GNP 4.51)

Both fail only the **same structural MB-HU-G2** target that option A also fails (`v44-hu-s-g2-structural-limit`). So
**option B matches-or-beats the baked default (30/32) with only a 3-parameter coarse adjustment, no full re-opt** —
and it makes p27 a genuine, Fan-Meyer-correct CDK4/6 brake.

**Refinement optimizer (`optionB_optimize.py`, ~150 candidates over 7 params): ceiling = 31/32**, CyclinD1 MB/GNP tops
out ~4.65 (passes on band; the 32nd is the unfixable HU-G2). So 31/32 is robust and is the option-B ceiling — same
structural limit as option A. Best option-B params (`optionB_optimize_best.json`):
`kSeqCd 1.23, kRelCd 0.0141, kDeKPC 0.00687, w_ink4 6.0, kSyP21 0.000707, kPhRbCd 0.339, K_CdRb 0.417`.

**Bonus — option B is more biologically consistent on p27 half-life.** The optimizer's best `kDeKPC=0.0069` gives free-
p27 t½ ≈ 1.5 h — exactly the spec's biologically-motivated target (§5.1) — whereas option A needed `kDeKPC=0.02`
(t½ ≈ 30 min, ~3–4× too fast) to fit. Option B accommodates the correct p27 stability where option A had to distort it.

**⚠ Caveat the 32 targets miss — partial-Shh dose-response.** The validation only tests SHH=0 (full withdrawal). The
fine SHH sweep shows this option-B combo (kPhRbCd raised to 0.32) makes GNP **cycle down to SHH≈0.10** where option A
arrests at SHH≈0.42 — i.e. option B, as coarse-calibrated, is **less Shh-sensitive on partial withdrawal** (it still
arrests at SHH=0, so validation passes). Neither option's partial-Shh behavior is currently a validation target;
before promoting option B, **add a partial-Shh dose-response target** (graded proliferation drop as Shh falls) so the
re-cal doesn't buy 31/32 at the cost of Shh dose-sensitivity.

## 4. Mechanism characterization — option B is necessary-not-sufficient (the honest result)

**p27 does NOT tune GNP G1 in either option — GNP's cycle is CyclinD1-governed.**
- kSyP21 sweep (`characterize_p27_v2.py`): option A period flat 17.6 h to 2× p27 then abrupt arrest at 3× (p27 inert
  then cliff); **option B period flat 17.6 h across 0.5–8× p27** (no lengthening).
- Fine SHH sweep (`shh_g1_sweep.py`): option A is a sharp switch (17.5 h to SHH≈0.42 → arrest below 0.4); **option B
  cycles at ~17.6 h down to SHH=0.10** (the raised kPhRbCd pushed the arrest threshold near SHH=0). **Decisive control:
  option B with p27 removed (kSyP21≈0) is IDENTICAL to option B at every SHH** → p27 contributes nothing to the GNP
  cycle. GNP operates with CyclinD1/p27 ≫ 1 (Cd≈3, free p27 tiny), i.e. **far above the ultrasensitive threshold,
  where the switch is insensitive to p27** (exactly the Fan-Meyer "G1-knob only near the switch point" caveat).
- **Why:** to keep GNP cycling at ~16 h you need enough CyclinD-CDK4/6 drive, which pushes the ratio *away* from
  threshold — the same drive-dominance tension as the stochastic-commitment finding. There is no option-B calibration
  where GNP both cycles at 16 h *and* sits near the p27 threshold.

**Reversibility asymmetry is absent (needs Phase B5).** CDK4/6i washout (`characterize_p27_v2.py`): option B restarts
*both* GNP (29.7 h) and MB (33.7 h) — MB is if anything *slower*, no MB-reversible/GNP-terminal asymmetry. The current
model has no differentiation/DREAM arm, so it cannot distinguish GNP terminal exit from MB reversible G0. That
asymmetry is spec **Phase B5 (DREAM/MMB + Neurod1 differentiation)**.

**Implications for JP's two goals:**
- *"p27 tunes GNP G1"* — not delivered by option B alone. Graded GNP G1-lengthening in this model comes from the
  latent `K_g1len` growth-slowdown (already present, set below GNP CyclinD1), or from the DREAM shallow-switch (B5), or
  requires operating GNP near the CyclinD1/p27 threshold (a calibration trade-off vs the 16 h cycle). p27 tunes G1
  **only near the exit transition** (partial mitogen withdrawal), consistent with Uziel (p27 comes on as GNPs exit).
- *"p27 drives MB G0"* — option B makes p27 a real CDK4/6 brake, but MB's high CyclinD1 keeps MB above threshold too;
  MB transient-G0 needs the **birth-p27 population layer** (rare/high-p27 subpopulation) and/or **B5 reversibility**,
  not the mean-field p27 brake — again consistent with the stochastic-commitment result.

## 5. Recommendation (for JP's review — NOT baked)

Option B is a real, low-cost mechanistic upgrade (31/32, Fan-Meyer-correct p27 brake, fixes the "p27 inert" defect),
but the characterization shows it is **necessary-not-sufficient** for the two goals. Sequenced plan:

1. **Add a partial-Shh dose-response target to the harness FIRST**, then re-run the option-B optimizer against it. This
   guards against the dose-sensitivity regression (§3 caveat) and is the prerequisite for promoting option B. Only then
   **promote option B to default** (pending JP sign-off; the standing no-bake rule is respected — it's flag-off now).
2. **p27 tunes GNP G1 requires operating near the CyclinD1/p27 threshold** — which conflicts with the 16 h GNP cycle at
   high drive. Realistic routes: (a) the latent `K_g1len` growth-slowdown (already present) for graded G1; (b) B5's
   DREAM shallow-switch; (c) confine p27's G1 role to the exit transition (partial withdrawal), matching Uziel.
3. **MB transient-G0 needs the population layer** (rare/high birth-p27 subpopulation) — the mean-field p27 brake keeps
   MB above threshold (high CyclinD1). This matches the stochastic-commitment finding: G0 is a variance/tail + exit
   phenomenon, not a mean-abundance one.
4. **Next structural module = B5 (DREAM/MMB + Neurod1 differentiation)** — delivers the MB-reversible / GNP-terminal
   reversibility asymmetry that option B alone cannot (both cell types restart after CDK4/6i washout today).
5. Deferred: B2 (explicit INK4 complexes), B3 (p53 carrier — Cd_mRNA reservoir already exists), B6 (p19/E2F).

**Bottom line for JP:** Phase A was already done (spec stale). Option B (your endorsed Fan-Meyer p27-inhibitory
mechanism) is implemented, reversible (flag-off), and reaches 31/32 — promote it after adding a Shh dose-response
target. But the model confirms what the stochastic work predicted: **p27 as a mean-field brake can't source the
GNP-G1-tuning or the MB-G0 by itself** (GNP/MB both sit far above the CyclinD1/p27 threshold); those need the
near-threshold/exit-transition regime, the population layer, and the DREAM module (B5).

## Files (this session)
`simulations/phaseA_verify.py`, `optionB_probe.py`, `optionB_recal_coarse.py`, `optionB_optimize.py`,
`characterize_p27.py`, `characterize_p27_v2.py`; builder flag `with_p27_optionB` in `src/build_model_v44_heldt.py`;
harness flag `P27_OPTIONB` in `simulations/validate_v44.py`.
