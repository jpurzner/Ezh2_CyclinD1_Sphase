# Per-gene lattice model: Ccnd1 and Ccnd2

2026-10-25. SLURM 12382646 (nucleation x abundance) and 12382673 (antagonism
sweep). Code: `src/lattice/gene_lattice.py`. 120 cells per condition, 300 h,
replication every 24 h, profiles time-averaged over the second half.

## Setup

The lattice is 100 nucleosome sites at 200 bp over TSS +/- 10 kb, oriented
5'->3', built directly from each gene's 25 bp ChIP table.

- **Input (measured):** H3K4me3 and H3K27ac give the "active" occupancy that
  antagonises PRC2; the NDR sets where histone exchange peaks.
- **Output (predicted):** the H3K27me3 profile.
- **Population:** cells differ through the stochastic dynamics plus a
  lognormal spread in PRC2 abundance (CV 0.4) and in how strongly the active
  machinery is engaged (CV 0.25).

So me3 is not fitted to itself. The question is whether the measured active
profile plus literature PRC2 kinetics puts the me3 boundary in the right place
with the right sharpness.

## 1. The window contains two promoters, and that matters

Every one of these genes has another TSS within a few kb upstream:

| gene | neighbour | offset |
|---|---|---|
| Ccnd1 | Gm71874 / Gm57215 | −1,809 / −2,703 bp |
| Ccnd2 | Gm56941 | **−3,726 bp** |
| Cdk6 | Cdk6os1 | −1,678 bp |

Ccnd2's upstream me3 block peaks at −3.9 kb — on top of Gm56941. Seeding PRC2
at the neighbouring promoters rather than uniformly is the single biggest
determinant of fit:

| nucleation | Ccnd2 r | Ccnd1 r |
|---|---|---|
| uniform | 0.226 | −0.778 |
| antisite (away from active) | 0.284 | −0.777 |
| **at neighbouring promoters** | **0.876** | 0.320 |

## 2. Ccnd2 is reproduced well; Ccnd1 is not

Best fit for both is nucleation at neighbours with **PRC2 = 8 molecules per
domain** (8 beats 15 beats 30 everywhere — too much PRC2 fills the window,
the same result the earlier abundance sweep gave).

- **Ccnd2: r = 0.876**, RMSE 0.097. The simulated domain sits over Gm56941 and
  clears across the Ccnd2 core, matching the measured block at −5 to −2.5 kb.
  The simulation runs slightly hot downstream (0.1–0.2 against a measured
  ~0.02).
- **Ccnd1: r = 0.320**. It reproduces the −4 to −1.5 kb block but completely
  misses the measured me3 peak at 0 to +1.5 kb — the bivalent core, where the
  real promoter carries me3 15.1 *and* K4me3 10.0 simultaneously.

## 3. The cleared core is made by histone exchange, not by the K4me3 block

I expected Ccnd1's failure to be the antagonism term forcing me3 to zero
wherever the active mark is high, so I swept `w_inh_active` from 0.90 to 0.00.
**It changes nothing.** The core error holds at −0.63 across the whole range:

| w_inh | 0.00 | 0.30 | 0.60 | 0.90 |
|---|---|---|---|---|
| Ccnd1 core error | −0.633 | −0.635 | −0.639 | −0.643 |
| Ccnd1 r | 0.296 | 0.317 | 0.325 | 0.320 |

With catalytic antagonism switched off entirely, PRC2 is free to methylate the
core and still no me3 accumulates there. The mechanism actually clearing the
core is **nucleosome turnover at the NDR** (`k_exch` peaked at +150 bp), which
strips marks irrespective of what PRC2 is doing.

That reframes the model: the depleted core is a turnover phenomenon, and the
K4me3/ac antagonism does something different — it sharpens the *boundary*.
For Ccnd2, more antagonism is steadily better (r 0.721 → 0.876 as w_inh goes
0.00 → 0.90) even though it does not touch the core.

## 4. The population is NOT a mixture

The expectation was that overlapping marks in bulk would require a variety of
single-cell patterns. At these parameters they do not. In the simulated
ensembles, only 2–3 of 100 sites have both a substantial high population
(>15% of cells above 0.6) and a substantial low one. Inside the domain ~85% of
cells are high and essentially none are low; outside, ~80% are low. The bulk
profile is reproduced by cells that all look alike.

The one place cells genuinely differ is a ~1 kb zone at each domain edge,
where the boundary position jitters. That jitter is what smooths the bulk
edge — but it is boundary noise, not two populations.

**Caveat:** this is a statement about this parameter regime, not about the
system. Whether sites are bistable depends on the read-write gain against
turnover; a bistable regime would give genuine mixtures. Testing that means
sweeping `allo_max` and `k_exch` together, which has not been done.

## What this says to do next

1. **Ccnd1's bivalent core is the open problem.** It needs either much lower
   exchange at that NDR than the model assumes, or an me3 species that
   survives turnover, or genuine per-cell alternation (me3 in some cells,
   K4me3 in others) that the current parameters cannot produce.
2. **Check the exchange rate against data.** It is currently the dominant
   core-clearing term and is set from a single literature value
   (`k_exch = 0.60/h` peak H3.3 exchange at the NDR). The whole core result
   rests on it.
3. **Sweep for bistability** before concluding anything about population
   heterogeneity.
4. Cdk6 has not been run; its geometry is already built.
