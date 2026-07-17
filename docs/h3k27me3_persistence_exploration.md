# H3K27me3 persistence: why the mark decays too fast, and how to lengthen it

*2026-07-16 (overnight, JP request: "H3K27me3 decays too fast (~2h), reduce the Gli eraser, EZH2 is higher
expressed, and transcription is sporadic — punctate modeling may make the memory easier"). Exploration in the
chain default model + reduced 1-locus prototypes. Sims: simulations/sim_mark_persistence_sweep.py, sim_punctate_meanfield_proto.py, sim_punctate_bistable_proto.py, fig_h3k27me3_persistence.py.*

## TL;DR
- **The fast decay is the Gli ERASER, and reducing it works (your lever).** The me3→me2 demethylation rate is
  `del_mk + k_jmjd3_gli·Gli1` = **0.044 in MB → intrinsic me3 half-life 0.2h**. Reducing the eraser (`k_jmjd3_gli`
  0.221→0) *and* turnover (`del_mk` 0.0023→0.0006) → **18.4h** (~90× longer, a genuine ~1-cycle memory).
- **Dilution is NOT the limiter:** EZH2i *slows* the MB cycle (32h→95h in-model), so replicative dilution is
  infrequent after EZH2i — the decay is demethylation-(eraser)-dominated. Good: the eraser is the right knob.
- **The catch: MB<GNP (ChIP) is currently produced ONLY by the eraser.** Zeroing it sends MB/GNP 0.56→1.02 (the mark
  saturates equally). Transcription-eviction + dilution alone do NOT separate MB from GNP.
- **No continuous-model regime does all three** (long mark + MB<GNP ~0.5 + correct 5.07 fold): strong eviction
  (g_prc2=0.9) restores MB<GNP but over-de-represses CyclinD1 to ~10× (vs 5.07 target).
- **Punctate transcription, mean-field: marginal** (~20% longer half-life, identical MB/GNP) — soft bursts of a
  linear-in-M eviction ≈ continuous.
- **The real lever is BISTABILITY:** cooperative read-write (M^n, n≥2) + suppressed nucleation (a0→0) makes the mark
  LATCH; MB's higher transcription tips it OFF while GNP stays ON → **MB<GNP via the switch, no eraser, with
  persistence** (prototype: GNP latched ON=0.67, MB OFF=0.01). BUT prior work found strong bistable memory is
  **DATA-EXCLUDED** (conflicts with the 5.07 MB/GNP fold + soft vismodegib data — see the
  chromatin-memory-validation-bound analysis), so this needs a joint re-optimization to check.

## The core tension (why it isn't a one-liner)
The Gli→Jmjd3/Kdm6b eraser was introduced to reproduce JP's ChIP (MB has ~half GNP H3K27me3 at Ccnd1). It does two
things at once: (i) makes **MB<GNP** (MB high-Gli → strong erasure), and (ii) makes the mark **decay fast** (the
same erasure). These are the same term, so you cannot slow the decay by weakening the eraser without also losing
the ChIP separation — unless the ChIP comes from somewhere else.

| eraser `k_jmjd3_gli` | intrinsic me3 half-life (MB) | MB/GNP mark ratio (ChIP) |
|---|---|---|
| 0.221 (current) | **0.2h** | **0.56** ✓ |
| 0.08 | 0.6h | 0.79 |
| 0.03 | 1.4h | 0.92 |
| 0.00 + del 0.0006 | **18.4h** ✓ | **1.02** ✗ (MB≥GNP) |

## What each proposed lever buys
- **Reduce the Gli eraser (your main suggestion): YES for persistence** — the dominant fix (0.2h→up to 18h). Cost:
  the ChIP. → must be paired with a replacement MB<GNP mechanism.
- **Reduce turnover `del_mk`:** extends the mark further once the eraser is down (del 0.0023→0.0006 roughly 4×).
- **Higher EZH2 (you noted EZH2 is higher expressed — confirmed, MB EZH2 3.6 vs GNP 2.0):** strengthens the *writing*
  (read-write), which helps MAINTAIN the mark against dilution during cycling (the memory), and supports a weaker
  eraser. Does not change the after-EZH2i *decay* (writing is off then).
- **Punctate transcription:** marginal in the mean-field. Its genuine value is only in the BISTABLE/stochastic
  regime, where strong brief bursts stochastically switch a latched locus OFF — which is a real (larger) model
  change and the same regime flagged as data-tension.

## Recommendation (for your review — not applied to the default)
1. **Reduce the eraser + turnover** to lengthen the mark (target ~5–18h, a genuine ~1-cycle memory).
2. **Move MB<GNP off the eraser onto a MILD cooperative-read-write bistability** (make the read-write `a_rw·Mk^n`,
   n≈2, with a small nucleation floor a0) + transcription-eviction: MB's higher transcription tips the switch, GNP
   stays latched — this gives BOTH the ChIP separation AND the persistence, without a strong eraser.
3. **Jointly re-optimize** the chain (eraser, turnover, a_rw, n, a0, g_prc2) against the FULL target set, and
   critically re-check the two data walls (the 5.07 MB/GNP CyclinD1 fold and the soft vismodegib MB+HHi/MB), because
   prior work excluded *strong* bistable memory. The open empirical question is how long the mark *actually* persists
   — a genuine chromatin timescale (ChIP-seq over an EZH2i time-course) would pin the target and tell us how much
   bistability the data allows.
4. **Punctate transcription** is worth a stochastic prototype only if (2)–(3) land in a bistable regime; in the
   mean-field it is not the lever.
