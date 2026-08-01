# Mechanistic p21–PCNA / CRL4^Cdt2 arm — the "un-map" (2026-08-01)

Status: **exploratory, behind `with_p21_pip_degron=False` (default OFF); nothing baked.**
Harness: `P21_PIP_DEGRON=1` env flag on `validate_v44.py`. Builder: `src/build_model_v44_heldt.py`
(un-map block at the end of the `with_cdki_species` section). Probes: `scratchpad/p21arm/test_unmap*.py`.

## The provenance bug this fixes

Heldt-2018's species named **`P21`** is literally **p21 (Cdkn1a)** — it binds PCNA (`iPcna` = p21·PCNA),
binds loaded replication complexes (`iRc` = p21·loaded-PCNA, the "RCi"), and is cleared by CRL4^Cdt2 on
those complexes. v44 **remapped `P21` to p27 (Cdkn1b)**. Consequence in the shipping model: **p27 binds
PCNA and inactivates replication forks (`aRc + P21 -> iRc`)** — biologically impossible (p27 has no PIP box).
The PIP-degron machinery is *inherited, not missing*; it is attached to the wrong protein. My first design
(a parallel `p21Pc` subsystem) reinvented `iPcna`/`iRc` while leaving the mislabeled originals running.

## What the un-map does (`with_p21_pip_degron=True`)

Structural only; quantitative knobs left at runtime. Repoints the PCNA/Rc binding **p27→p21a**
(`aPcna/aRc + p21a -> iPcna/iRc`, and `iPcna => p21a` export); restores the **CRL4^Cdt2 degron
(`kDeP21aRc·Cdt2·aRc`) on the p21a-bound `iPcna`/`iRc`**; drops the free-pool `aRc` proxy; and fixes the
moieties (`tP21 = P21+CeP21+CaP21+CdP21`; new `tp21a = p21a+Cep21a+Cap21a+iPcna+iRc`). p27 keeps
CDK2 (CeP21/CaP21) + CyclinD1 buffer (CdP21) + basal/KPC/Skp2; **no PCNA**. p57 gets no Cdt2 arm
(no PIP box; Thr310→Skp2/FBL12).

## THE key distinction: sequestration vs degradation

At a cycling operating point, free p21 ≈ 1e-4 while `iPcna` ≈ 0.19–0.50 — **essentially all p21 is parked
on *soluble* PCNA (`iPcna`)**, which is correctly *not* a CRL4^Cdt2 substrate (Havens & Walter 2009:
CRL4^Cdt2 needs DNA-loaded PCNA). So in this model p21 is held low by **sequestration onto soluble PCNA**,
**not** by the degron. The degron acts only on the loaded pool (`iRc`), which is small unless p21 is abundant
*during S*. This is the single most likely thing to be re-derived wrongly later, hence this note.

## Findings (probes test_unmap1–4)

1. **Positivity / moiety — solid, and the round-3 diagnosis was right (JP).**
   Moiety conserved to machine epsilon (all p21 sources+sinks off, `drift = 3.9e-16`). Degradation-side
   positivity robust: `p21a_min = +7e-8` even at `Cdt2 = 50` (fast degron, zero supply). Round-3's `−0.4`
   was a **broken moiety** (near-zero pool × `kAsPcP21=100` stiffness), not the un-map itself. "Make every
   sink first-order" was a weak claim (positivity in exact arithmetic only); the real invariant is the moiety,
   and it holds.

2. **Steady-synthesis sweep — p21-at-abundance arrests in G1; binding affinity is not the bottleneck.**
   As `kSyp21a` rises past ~0.003, the cell arrests with `aRc = 0` (never enters S) — a CDK2/G1 block, not
   fork poisoning. A 2D sweep `kSyp21a × kAsPcP21` (100→3000, i.e. giving p21 its true tight PIP binding)
   does **not** open a cycling-with-meaningful-`iRc` window. So "p21 protein stays low in cycling cells" is
   enforced by the commitment architecture — consistent with observed biology (p21-GFP undetectable in S).
   *This is a sanity check that the model isn't wrong, not a new claim.*
   **Caveat (important):** a steady sweep tests a *threshold*; the Barr phenotype is *hysteresis* and cannot
   be shown this way. See (3).

3. **Step experiment (spike p21 mid/early-S) — the degron IS load-bearing, in the Barr direction.**
   Commit at cycling params, spike `kSyp21a` at S-entry, compare degron ON (`Cdt2=1`) vs OFF (`Cdt2=0.01`,
   Barr's ~1% depletion):
   - `iRc` (p21 on loaded PCNA): **0.059 OFF vs 0.010 ON** — the degron clears loaded-PCNA p21.
   - **S-phase duration: 4.0 h (no spike) → 4.2 h ON → 4.8 h OFF** — Cdt2 depletion lengthens S via p21
     accumulating on forks. **Correct Barr direction, degron-dependent.** Modest magnitude (~+14%).
   - **But cell *fate* is degron-independent here:** both ON and OFF complete the current S and both arrest in
     the *next* G1 (via free p21 ≈ 0.9). There is **no feedback from `iRc` clearance back onto G1 commitment**,
     so the degron modulates S-progression (Barr's "irreversibility" role) but not the commitment decision.
   This overturns my earlier "the degron acts too late / has no vote" conclusion — that was the wrong
   experiment (threshold, not hysteresis) and the wrong call.

4. **The un-map is NOT validation-neutral: 27/32 (vs 30/32).**
   New fails: **CyclinD1 MB/GNP 4.61→3.79** (real), CyclinD1 GNP+HHi 0.094 and EZH2 MB/GNP 2.77 (both at
   band edges). Mechanism: removing p27's (wrong) fork inhibition changes the replication/S-phase dynamics,
   which shifts EZH2-in-S accumulation (EZH2 MB/GNP rises → represses CyclinD1 more → fold drops). **The
   current 30/32 calibration silently depends on p27's biologically-impossible fork inhibition.** So the
   un-map is a correctness fix that *requires re-calibration*; it cannot simply be promoted to default.

## Caveats on interpretation

- **kSy vs transcript:** Cdkn1a 890→2511 (MB/GNP 2.8×) is **transcript**, not protein. Pinning `kSyp21a`
  at 2.8× assumes translational efficiency is unchanged — probably false here: the **miR-17/92 cluster**
  (amplified in SHH-MB; a Shh/N-myc target; Dicer1-loss raises p21, miR-17-5p mimic lowers p21 protein in
  GNPs) is a post-transcriptional brake that scales *up* with the same axis driving the transcript. The
  honest pin is `kSyp21a × (transcript fold) × (translational-efficiency ratio ≤ 1)`, ratio **unconstrained**,
  not silently 1. "mRNA overstates protein" may be explained by miR-17/92 without invoking the degron.
- **The real MB question is heterogeneity, not on/off.** The sharp deterministic cliff at `kSyp21a ≈ 0.003`
  is the deterministic shadow of Barr's bistability (~21% arrest, set by inherited damage). What a 2.8×
  transcript shift speaks to is *the fraction of cells either side of the cliff* — a population-layer question.

## Decision (open; nothing baked, flag default-off)

- **A′ — keep the un-map as the correct structure, re-calibrate to recover 30/32.** Two of three new fails are
  band-edge; CyclinD1 MB/GNP (3.79 vs ≥4.06) needs ~7% recovery, likely reachable by re-tuning the
  EZH2-in-S / fold knobs. Payoff: fixes the p27-binds-PCNA absurdity and gains a modest but real Barr-direction
  degron effect on S-duration. Recommended if the p21/Cdt2 axis is worth carrying correctly.
- **B′ — add the degron→commitment feedback (Barr-style bistable G1/S).** The only way to make the degron
  *fate*-relevant (full Cdt2-depletion phenotype) and to address the heterogeneity/fraction question. A new
  commitment module coupling `iRc`/p21 clearance into CDK2 activation; needs the population layer; re-validates
  the whole 30/32 (commitment is load-bearing for the cKO/governor results). Scope as a separate project.

Not promoted to default. `with_p21_pip_degron=False` keeps the shipping model byte-identical (30/32).

## A′ re-calibration outcome (2026-08-01): a Pareto wall — best 29/32

JP authorized A′ (re-cal the un-map to recover 30/32). Result across ~20 full-harness evals
(`P21_PIP_DEGRON=1`, `VALIDATE_PARAMS` sweeps; scratchpad/p21arm/recalA*, rA/rB/rC/rD*):

- **EZH2 MB/GNP is recoverable.** Raising `Kez_cd` 9.18→11 (NB the naive direction is inverted by the
  Cd→EZH2 feedback: *lower* Kez_cd *raises* the fold) restores EZH2 MB/GNP to band and keeps Palbo + S/G0.
  → 28/32.
- **CyclinD1 MB/GNP is the lone non-HU casualty, and it sits on a Pareto wall vs GNP+HHi.** Under vismo
  Gli→0, so raising the Gli drive to lift MB/GNP *lowers* GNP+HHi, and any K_prc2/basal comp that fixes
  GNP+HHi drags MB/GNP back down. Best MB/GNP-in-band config (Gli×1.7, t1) = **4.008** but GNP+HHi 0.084 ✗
  (28/32); best all-else-pass config (MYCN×1.5+basal, rD_mycn) = **29/32** with MB/GNP 3.59 ✗. No config of
  {Kez_cd, k_Cd_tx_Gli_max, k_Cd_tx_MYCN, k_Cd_tx_basal, K_prc2, kEZbas_Cd} threaded both. **30/32 not reached.**

**Two things this means.**
1. **The correctness fix costs ~1 validation target, and it's the single most normalization-ambiguous one:**
   CyclinD1 MB/GNP is 5.07 (Fig 4I) vs 7.58 (raw RNA-seq) — a 1.5× fork *wider than its own 0.20 band*. So
   the un-map's cost (MB/GNP ~3.6–4.0 vs a [4.06,6.08] band drawn on a disputed normalization) is arguably
   within the target's real uncertainty; t1's 4.008 is essentially at the band floor.
2. **The baseline 30/32 was leaning on p27's biologically-impossible fork inhibition to hit CyclinD1 MB/GNP.**
   Removing that wrong edge (the whole point of the un-map) is *why* the fold drops — the calibration had
   silently absorbed the error. This is a real finding about the model, not just a re-fit failure.

Decision (JP): (a) accept **29/32** and promote the un-map (recommended — correct biology + real
Barr-direction degron effect, cost is the most-ambiguous target); (b) launch a full multi-hour joint
optimizer (more knobs — kDeEZ, a_rw_prc2, cell-type-specific EZH2 — might squeak 30/32, but the Pareto
evidence says tight); (c) revert and keep 30/32 with the known-wrong p27-PCNA wiring. **Nothing baked;
flag default-off; shipping model still 30/32; all re-cal params are runtime-only (not written to `_ts_bake`).**

## A′ RESOLVED → **30/32 BAKED** (2026-08-01) — the "Pareto wall" was a species-conflation artifact

The 29/32 wall above was measured against a **mis-specified** target set: **EZH2 MB/GNP was scored on the
model's EZH2 _protein_ while the data (Fig 4J) is RNA-seq _transcript_** (JP confirmed), and Gli1 folds likewise
scored protein vs RNA-seq. Fixing the harness to score the assay-matching species (`EZH2m`, `Gli1_mRNA`),
Palbo in MB (not GNP), and Period at 16 h (already met — the GNP model runs ~17.8 h) changed the landscape:

- **Species-correct DEFAULT (flag-off) = clean 30/32** — EZH2 MB/GNP 2.025 on `EZH2m` (≈ the protein value, so
  the fix didn't move it), Period 17.6 h passes 16 h with no mu re-tune.
- **Un-map, species-correct = 28/32** — the EZH2 MB/GNP "conflict" (protein 2.77) **vanished** on transcript
  (2.31, in band). The only real residual was CyclinD1 MB/GNP (3.79 vs a 4.06 floor the default itself sits at)
  and GNP+HHi at its floor. The un-map's damage is **MB-specific** (removing p27 fork-inhibition raises MB
  EZH2-in-S more than GNP — MB has more p27), so it was fixed with an **MB-specific CyclinD1 driver that spares
  GNP+HHi**: the MYCN→CyclinD1 coupling (`k_Cd_tx_MYCN` ×2) + a small Gli bump (`k_Cd_tx_Gli_max` ×1.15).

**BAKED (verified 30/32, no env overrides):** `with_p21_pip_degron=True` (builder default) + harness
`P21_PIP_DEGRON` default `1`; `_ts_bake` `k_Cd_tx_Gli_max 0.152118→0.17494`, `k_Cd_tx_MYCN 0.031157→0.062`.
Fails only the 2 structural HU folds — identical to the prior default, now with correct p21/PCNA biology and
species-correct scoring. `Kez_cd` unchanged (the earlier Kez_cd sweeps were chasing the mis-scored protein fold).

**Lesson:** the whole "un-map costs a target / Pareto wall" conclusion was an artifact of a conflated target.
JP's insistence on annotating measured species (protein vs mRNA) before re-optimizing was exactly right and
directly unlocked the 30/32. **Still open (orthogonal):** CyclinD1 MB/GNP target 5.07 (Fig 4I) vs 7.58 (raw
DESeq2) — the model does ~4.09 here either way; 7.58 would fail the default too, so 5.07 is the operative value.
