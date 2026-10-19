# Overnight run — CDKI two-step mRNA/protein + transient-G0

Started 2026-10-18 21:02:21.019220 · finished 2026-10-18 21:10:52.125929

## Phase A — locked ~34% reference (two-step, model-tuned protein)
- **MB**: snapshot-G0 **36%**, transient 72%, permanent 2% (n=182/300, mortality 39%)
  - config: A_ON=2.9, duty=0.28, tau=700 min; p21 mRNA CV 0.95 -> protein CV 1.21
  - EZH2 quiescent 2.15 vs cycling 3.69 (fold 0.58)
- **GNP**: snapshot-G0 **0%**, pure-cycling 100% (n=300/300, mortality 0%)

## Phase B — measured-protein INK4 (p18 23.9h, p21 0.98h)
- validation: VALIDATION SUMMARY: 30/32 targets passed
- 2x2 attribution (MB snapshot-G0):
  - model-tuned protein:  p21+p18 = 33%, p21-only = 38%
  - measured protein:     p21+p18 = 57%, p21-only = 50%
  - GNP controls: tuned 0%, measured 0%
  - read: if measured p21+p18 ~= measured p21-only, stable p18 has dropped out as a driver -> p21 alone carries it.

## Phase C — bake confirmation
- baked default (two-step ON): VALIDATION SUMMARY: 30/32 targets passed
- legacy single-step (CDKI_MRNA=0): VALIDATION SUMMARY: 30/32 targets passed

## Files
- overnight/results.json (all metrics) · overnight/progress.log (timeline)
