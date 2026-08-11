# Population layer for MB Atoh1+ transient-G0 (2026-08-09)

Closes item #3 of the 2026-08-07 refocus ("some form of transient G0 with high p27 periodically in
cells"), now tested in the **full integrated governor** rather than the reduced stochastic model.

## The model under test

Integrated v44 governor, all three second-arm mechanisms ON (flags, not baked):
- `with_p27_optionB` — buffered p27 (CdP21) is INACTIVE and brakes CDK4/6 (Fan-Meyer).
- `with_cdk6_gli` — CDK6 as a Gli-driven, EZH2-marked second governor arm (mark-memory IFFL) + an
  unbound-CDK6 INK4 sink.
- INK4 = p18 (dominant, stable) + p19 (short) sinked by unbound CDK6; p21 dropped (transcript ~0).

Re-calibrated to **31/32** (`integrated_optimize_best.json`); CyclinD1 MB/GNP = 5.07 exact; GNP still
arrests only at SHH=0; MB+CDK4/6i arrests. Nothing here is baked into the default model.

## Three analyses (`simulations/population_layer.py`)

**(1) Bulk CDKI heterogeneity (lognormal, CV=0.6 on p18/p19/p27) → 0% G0** in every condition
(MB, INK4-KO, p27-KO, CDK6×0.5). No cell in a realistic-noise population arrests.

**(2) Threshold sweep — how high must a cell's CDKI be to arrest?**

| CDKI fold-of-mean | 1× | 2× | 3× | 5× | 8× | 12× |
|---|---|---|---|---|---|---|
| divisions (default CDK6) | 6 | 4 | 3 | **0** | 0 | 0 |
| divisions (CDK6 ×0.3) | 6 | 4 | 3 | **0** | 0 | — |

MB arrests (G0) only at **CDKI ≥ ~5×** the mean, with a **graded slowdown** (6→4→3 divisions) before
the switch. Cutting CDK6 to 0.3× does *not* lower the threshold — the drive, not CDK6 titration, sets it.

**(3) Bimodal high-CDKI tail → discrete G0 fraction.** A distinct subpopulation at ~5–6× CDKI on top of
normal noise:

| tail | G0 fraction |
|---|---|
| 15% at 5× | ~10% |
| 20% at 5× | ~5% |
| 25% at 6× | ~12% |

The Atoh1+ G0 fraction tracks the tail prevalence (± N=40 sampling noise), and it is **CDKI-gated**
(the high-CDKI cells are the arrested ones).

## Conclusion

The MB Atoh1+ transient-G0 subpopulation is a **rare high-CDKI tail** — a distinct stress/niche state at
~5× the mean p18/p19/p27 — **not bulk expression noise**. This is the drive-dominance / superset result
recurring in the full model: MB's calibrated drive is high enough that only an extreme-CDKI cell can
arrest, so G0 lives in the **variance (a bimodal tail), not the mean**. It reproduces, in the full v44
governor, the reduced-model stochastic-commitment finding (`stoch_commit_g0/`): MB G0 works only via an
MB-specific high-CKI tail.

Biologically consistent with JP's framing: MB has *some* Atoh1+ G0 cells (the high-CDKI tail), *largely
overwhelmed by very high CDK6* (the bulk cycles); the excess unbound CDK6 sinks INK4 and keeps the bulk
committed, while the rare tail escapes the sink and arrests.

## Design fork for JP (nothing baked)

1. **Ship the tail explicitly** — add a population layer where a fraction `f_g0` of MB cells occupy a
   high-CDKI (~5×) state → an Atoh1+ G0 fraction = `f_g0`. Controllable, matches biology, but the tail
   prevalence is a free parameter (no direct measurement yet).
2. **Accept the integrated governor (31/32) as the deliverable** — the discrete G0 subpopulation is a
   rare-cell/stress phenomenon that the deterministic core correctly *cannot* produce from bulk
   heterogeneity; document it as a tail requirement rather than baking a population layer.
3. **Reduced population model (v45)** — keep the transient-G0 population dynamics in the dedicated
   stochastic successor and leave v44 as the deterministic governor.

Recommendation: option 2 for the paper's core (the governor is the thesis and it's clean at 31/32), with
the tail (option 1) as a clearly-labelled population addendum if a measured Atoh1+ G0 fraction becomes
available to pin `f_g0`.
