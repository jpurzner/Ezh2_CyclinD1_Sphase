# Resolving the PRC2 / transcription paradox (and why fast H3K27me3 is wrong)

*2026-07-17 (JP: "transcription is antagonistic with PRC2, yet we see very high PRC2 + high expression, peaking at
maximal division; fast H3K27me3 doesn't make biological sense — resolve this paradox"). Reduced-locus + allele-
population prototypes; simulations/fig_prc2_transcription_paradox.py.*

## The paradox
If transcription evicts PRC2, high-transcription cells should have low PRC2/H3K27me3. Yet MB (high CyclinD1
transcription) has **high EZH2/PRC2, peaking at maximal division** — the opposite. And the working model makes
H3K27me3 decay in ~0.2h, which is not a memory mark.

## Resolution, in three moves

**1. Global PRC2 ABUNDANCE ≠ local PRC2 OCCUPANCY (this dissolves most of it).** "High PRC2" (IF/WB) is *global
EZH2 protein*. EZH2 is an **Rb-E2F target**, so it is HIGH in proliferating cells and **peaks at division** — it is
a *proliferation readout*, not a CyclinD1 repressor. At the CyclinD1 *locus*, the story is local: the high
transcription evicts PRC2 → **low local H3K27me3** (your ChIP MB<GNP). So high global PRC2 and high CyclinD1
expression coexist trivially — both are downstream of mitogen/E2F, and the peak-PRC2-at-division is just E2F driving
EZH2 (and the G2/M restoration re-methylating the S-phase-diluted mark). *(Figure A.)*

**2. The mark is STABLE; MB<GNP comes from the transcription ANTAGONISM, not a Gli eraser.** The model's ~0.2h decay
is an artifact of the Gli→Jmjd3/Kdm6b eraser, a phenomenological device inserted only to make MB<GNP. Biologically
H3K27me3 is a memory mark (slow turnover). The *right* source of MB<GNP is the eviction JP names: MB's higher
transcription locally excludes PRC2 → less mark. A prototype with **no eraser** and a slow (10h+) turnover confirms
the mark stays stable while mitogen-tuned eviction pushes it down.

## The deep obstacle we found (and why the eraser existed)
The transcription↔PRC2 antagonism is a **double-negative = positive feedback → a BISTABLE SWITCH**. In mean-field
this has NO useful graded intermediate: weak eviction → the mark stays saturated (MB/GNP ≈ 0.92, no ChIP); strong
eviction → the switch flips all-or-nothing (MB mark → 0, CyclinD1 jumps ~9×, **losing the graded dose-response**).
The scan found essentially nothing between 0.92 and 0.00. So a *single deterministic locus* cannot give a moderate
MB<GNP (~0.5) AND graded CyclinD1 AND a stable mark — which is exactly why the model reached for the eraser (it
decouples the mark level from the switch, at the cost of a fast mark).

## The full resolution — and it validates the punctate-transcription idea
Move the switch to the **single-ALLELE** level and let transcription be **PUNCTATE (stochastic)**:
- each allele is **bistable with a STABLE mark** — latched ON (high mark, silent-ish) or OFF (low mark, active);
- **punctate transcription bursts stochastically flip an allele ON→OFF**, at a rate that rises with mitogen; PRC2
  re-methylation flips OFF→ON at a steady rate;
- the **POPULATION AVERAGE** (over alleles/cells) is then a **smooth, graded** readout: fraction-OFF rises with
  mitogen → average H3K27me3 falls smoothly (**MB/GNP ≈ 0.40, matching the ~2× ChIP**) and CyclinD1 rises
  (dose-graded, fold ~3) — with **every allele's mark stable/latched, and NO fast eraser**. *(Figure B.)*

So the mean-field couldn't see the answer, and that is the point: **JP's "transcription is sporadic / punctate" is
mechanistically load-bearing** — the stochastic flipping of bistable alleles, averaged over the population, is what
converts an all-or-nothing switch into the graded MB<GNP + dose-graded CyclinD1 we measure, while keeping the mark a
genuine (stable, latched) memory. Peak PRC2 at division stays E2F-driven and global.

## Model implications (a real restructuring, for review)
1. **Remove the Gli eraser** — it was standing in for the ChIP-level control; the fast decay was its side-effect.
2. **Make the CyclinD1 locus an allele-level bistable switch** (cooperative read-write + suppressed nucleation),
   with a STABLE mark, and drive OFF-flipping by **stochastic/punctate transcription** (rate ∝ mitogen).
3. **Read out at the population level** (the model already is a population-of-models) — the ChIP fold and CyclinD1
   fold are population averages of ON/OFF alleles; the mark half-life is the (long) per-allele latch time.
4. This is a **stochastic** module (Gillespie/tau-leap on the allele state), a larger change than a parameter
   re-fit — but it is the one that reproduces the biology without a fast eraser, and it turns the 5.07 fold + the
   ChIP + the stable mark from a three-way conflict into three population averages of the same switch.
5. **Empirical test that would pin it:** single-molecule/allele H3K27me3 + nascent-transcription (or an EZH2i ChIP
   time-course) — the model predicts *bimodal* per-allele mark (latched ON/OFF), not a graded per-allele level, and
   a long per-allele persistence with a mitogen-graded OFF-fraction.
