# JP's H3K27me3 amplifier hypothesis — findings (2026-08-09)

**JP (2026-08-09):** "there is likely substantial variance in CDKi that is generating this population, note
the longer they are out of cycle presumably the more ezh2 mediated H3K27me3 so we should see over a couple
of divisions this impact will accumulate on CDK6 and Cyclin D1."

Tested in the integrated governor (with_p27_optionB + with_cdk6_gli, integrated_optimize_best.json, 31/32).
Both CDK6 and CyclinD1 are mark-repressed IFFL targets, so the loop JP describes is structurally present:
slow/out-of-cycle → less H3K27me3 dilution → mark banks → represses CDK6+CyclinD1 → drive erodes.

## What the model actually does

**1. Premise HOLDS in slow-cycling cells.** As CDKI rises and the cycle lengthens (24→34→43 h at 1×/2×/3×),
EZH2 and the mark bank up:

| state | E2f | EZH2 | Mk | PRC2_rep | cyc (h) |
|---|---|---|---|---|---|
| cycling 1× | 0.88 | 5.25 | 0.935 | 0.024 | 24 |
| slow 3× | 0.90 | **7.28** | **0.956** | **0.034** | 43 |

**2. But the governor IFFL is SATURATED**, so the banked mark barely represses its targets: cdk6 flat
(5.084→5.082), Cd −3% across 1×→3×. No compounding — each CDKI level settles to a fixed slower limit
cycle and plateaus.

**3. The governor CANNOT be de-saturated to activate it.** De-saturating the CDK6/CyclinD1 IFFLs (raising
K_prc2 / K_prc2_cdk6, baseline auto-matched) does switch the amplifier on (cdk6 −20%, Cd −28% in slow
cells) and drops the MB arrest threshold from 5× to a tunable 2–4× — but it BREAKS the governor folds,
because a single Hill cannot be both saturated (big static fold) and high-slope (big dynamic response) at
one operating point. Full validation of 6 de-saturated candidates: **all 27–29/32**, failing consistently:
`EZH2i CycD1 fold (GNP) → 1.03 (target 2.2)` and `CyclinD1 MB/GNP → 3.3 (target 5.07)`. The EZH2i
de-repression fold and the MB/GNP fold *need* the IFFL saturated; the amplifier *needs* it at midpoint.

**4. Resolution attempt — a separate, excursion-gated amplifier term (`with_mark_amplifier`, flag, default
OFF).** A second Hill `mark_amp = f_mamp + (1-f_mamp)/(1+(PRC2_rep/K_mamp)^n_mamp)` multiplied onto CDK6
synthesis + CyclinD1 transcription, with K_mamp set *above* the cycling mark so it is ≈1 during normal
cycling (governor + all validation conditions preserved) and engages only on a mark excursion. This
preserves baseline (cdk6 5.08, Cd 7.8) but **does not lower the arrest threshold** (still 5×): PRC2_rep
*oscillates within each cycle* (diluted at each division), so even a sharp high-set Hill only fires briefly
near the mark peak (~15% time-averaged repression) — not enough against MB's large drive headroom.

**5. The WALL — the amplifier cannot deepen a full arrest, and the premise inverts there.** In a fully
arrested cell, E2f collapses → **EZH2 collapses** (EZH2 is a proliferation-coupled Rb-E2f target — a baked,
validated feature: EZH2 MB/GNP 2.05×, cell-cycle-phase ratios) → the writer stops → the mark decays:

| state | E2f | EZH2 | Mk | PRC2_rep |
|---|---|---|---|---|
| arrested 5× | 0.03 | 2.06 | 0.847 | **0.009** |
| arrested 8× | 0.03 | 2.02 | 0.844 | 0.008 |

So "longer out of cycle → more H3K27me3" runs **opposite** to the model: an arrested cell has *less* EZH2
and *less* mark than a cycling one. Hysteresis test confirms it — arrest at 5×, drop CDKI to 1×: re-entry
is identical (81 h) with the amplifier on or off, because PRC2_rep in arrest (0.009) is far below K_mamp.

Moreover, repression is EZH2-**occupancy** based (`PRC2_rep = EZH2 × mark-terms`), so even a persistent mark
`Mk` would not hold repression once EZH2 falls. JP's amplifier needs **writer-independent, mark-based
repression that persists through G0** — i.e. a bistable chromatin-memory latch.

## Bottom line

JP's amplifier operates (modestly) in **slow-cycling** cells but cannot create or deepen a **transient-G0
dwell**, because it collides with two established, validated model features:
- the **governor** needs the IFFL *saturated* for the EZH2i (2.2×) and MB/GNP (5×) folds — can't also be a
  high-slope dynamic sensor; and
- **EZH2 is proliferation-coupled** (Rb-E2f target) + repression is EZH2-occupancy, so out-of-cycle cells
  *lose* the mark and the repression — the opposite of the premise.

The amplifier is, in effect, a **chromatin-memory / bistable-latch** mechanism (mark sustains repression
through non-cycling), which prior work already found **data-excluded** by the vismo (0.144) and MB/GNP
(5.07) walls (~450 candidates, zero feasible; see memory `v44-chromatin-memory-validation-bound`,
`prc2-buffering-memory-topologies`, `v44-latch-regime-consequences`). This is the same wall, reached from a
new direction.

## Design fork for JP (nothing baked; `with_mark_amplifier` is a flag, default OFF, 31/32 preserved)

1. **Keep EZH2 proliferation-coupled** (validated) → the transient-G0 population comes from **CDKI variance
   alone** (the rare high-CDKI tail, memory `mb-atoh1-g0-rare-cdki-tail`); the amplifier stays a modest
   slow-cycling modulator, not the G0 driver. *(Recommended — preserves all validated biology.)*
2. **Overturn EZH2=Rb-E2f-coupling** (let the writer stay active, or make repression a persistent
   writer-independent `Mk` mark) so the amplifier can deepen G0 — but this re-opens the chromatin-memory
   walls (vismo + MB/GNP + ChIP MB<GNP) and needs a full re-derivation. A biological decision for JP:
   does H3K27me3 on CDK6/CyclinD1 **persist** in Atoh1+ G0 cells (memory), or **track** proliferation?
3. **Two-pool mark** (`h3k27me3-kinetics-vs-local-repression`): broad proliferation-coupled EZH2 for
   kinetics + a *local, stable* H3K27me3 on the CDK6/CyclinD1 loci that persists into G0 and does the
   repression. This is the most faithful to JP's "broad=kinetics, local=repression" framing but is a
   substantial new module.

Recommendation: **option 1** for the current deliverable (governor 31/32 + variance-driven G0 tail), and
option 3 as the principled next module if the transient-G0 memory needs to be mechanistic. The empirical
question that decides 2 vs 3 vs 1 is JP's: **is H3K27me3 on the CDK6/CyclinD1 loci higher in Atoh1+ G0
cells than in cycling cells?** If yes → local persistent mark (option 3); the model currently predicts the
opposite (proliferation-coupled), so this is a testable discriminator (e.g. CUT&RUN/ChIP on sorted
G0 vs cycling Atoh1+ cells).
