# Transcriptional bursting versus H3K27me3 encroachment

2026-10-25. SLURM 12385903 (sweep), 12385983 (rasters). Ccnd1, 60 cells per
condition, 400 h, statistics from the second half. Code:
`src/lattice/bursting.py`, initiation recording in `polii_lattice.py`.

## Setup

Pol II initiation events are the observable. PRC2 abundance is the knob (the
EZH2-dose experiment), and for each cell we record how far me3 has come into
the central ±1 kb together with its full initiation train. The result is the
**per-cell cloud of 480 cells**; the condition means only summarise it.

Burst scoring is validated against synthetic trains: a Poisson train returns
CV² = 1.03 and Fano = 1.16, a train with planted bursts of 8 returns CV² = 14.3
and mean burst size 8.4.

## The result

**Expression falls monotonically as me3 encroaches, and firing turns bursty.**

| | cleared core (me3 < 0.02) | encroached (me3 > 0.40) | |
|---|---|---|---|
| cells | 202 | 69 | |
| initiation rate | 1.08 /h | 0.14 /h | **7.7× lower** |
| CV² of intervals | 0.98 | 2.84 | **2.9× higher** |
| bursts per hour | 0.84 | 0.08 | 10× lower |

Per-cell correlations with encroachment, all 480 cells:

| quantity | Pearson r | Spearman ρ |
|---|---|---|
| initiation rate | **−0.961** | −0.901 |
| bursts per hour | **−0.953** | −0.898 |
| fraction of events in bursts | +0.811 | +0.699 |
| initiations per burst | +0.790 | +0.732 |
| Fano factor | +0.644 | +0.660 |
| CV² of intervals | +0.636 | +0.781 |

So the predicted direction holds, and tightly: a cleared core fires in a
near-Poisson, effectively constitutive way (CV² = 0.98), and encroaching me3
converts it to sparse bursty firing.

## It is frequency modulation, not size modulation

Across the encroachment range, bursts per hour fall **0.86 → 0.03** (a factor
of about 29) while initiations per burst rise only **1.27 → 1.9**. PRC2
encroachment therefore controls *how often* the promoter fires, not how much
it produces once it does. That is the signature usually reported for
repressive chromatin, and it is a prediction that single-molecule imaging
could test directly.

## Dose–response

| core me3 | cells | rate /h | CV² | events per cell |
|---|---|---|---|---|
| 0.00–0.10 | 272 | 1.05 | 1.02 | 211 |
| 0.10–0.20 | 59 | 0.73 | 1.65 | 146 |
| 0.20–0.30 | 45 | 0.51 | 2.52 | 102 |
| 0.30–0.40 | 33 | 0.35 | 2.89 | 69 |
| 0.40–0.50 | 28 | 0.22 | 3.26 | 44 |
| 0.50–0.60 | 19 | 0.13 | 2.65 | 26 |
| 0.60–0.70 | 17 | 0.07 | 2.53 | 14 |

The rate decays smoothly with no threshold — consistent with the monostable
structure established earlier. There is no switch; encroachment is a dial.

**CV² peaks around 0.4 and then declines.** That is a measurement artefact, not
biology: by 0.6 each cell fires only ~14 times in 200 h, so the interval
statistic is estimated from a handful of gaps. Read CV² only up to ~0.5.

## Caveats

- One gene (Ccnd1), one parameter set for the Pol II arm. The second parameter
  set in the sweep (`k_load` 1e-3, `p_evict` 0.1) reproduces the same pattern
  at higher absolute rates, so the shape is not a single-parameter artefact,
  but `k_load`, `p_evict`, `K_m` and `h` remain unconstrained by measurement.
- Encroachment here is driven by PRC2 dose. Driving it another way — moving
  the nucleation site, narrowing the active zone — should be checked to
  confirm the relationship is to encroachment itself rather than to PRC2
  abundance.
- Absolute rates are not calibrated to any measured transcription rate; the
  comparison that means something is the *ratio* across encroachment.
