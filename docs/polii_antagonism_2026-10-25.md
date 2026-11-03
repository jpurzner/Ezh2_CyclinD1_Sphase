# Adding RNA Pol II and PRC2 antagonism

2026-10-25. SLURM 12383613, 12383619, 12383630, 12383696, 12383751, 12383813.
Code: `src/lattice/polii_lattice.py`.

## What was added

Previously transcription entered only as a constant `tx` scaling a fixed,
NDR-shaped exchange field — turnover was an input, not a consequence, and the
antagonism ran one way. Pol II is now explicit:

- **Loading** Poisson at the promoter, repressed by promoter me3:
  `k_load_eff = k_load / (1 + (m_prom/K_m)^h)`
- **Elongation** at `v_elong` sites/s
- **Eviction** each site a polymerase crosses loses its mark with probability
  `p_evict` — so turnover is now produced *by* transcription

That closes the loop: me3 ⊣ Pol II ⊣ eviction ⊣ me3, a double negative, which
together with the existing active ⊣ PRC2 term is a mutual-antagonism toggle.

## The antagonism is real in the data

GNP Pol II S5P against core H3K27me3, measured:

| gene | core me3 | core Pol II |
|---|---|---|
| Ccnd1 | 15.1 | **4.4** |
| Cdk6 | 6.9 | 6.2 |
| Ccnd2 | 1.0 | **15.5** |

Perfect rank anticorrelation across all three. Three genes is three genes, but
it is the ordering the model assumes, and it was not used to fit anything.

## Mean field: the motif can be bistable

With Pol II-driven loss the loss curve is no longer a straight line — it is
steep while Pol II loads and flattens once me3 shuts Pol II off — so the
writing parabola can cross it three times:

    dm/dt = beta(1-m)[1+(A-1)m] - [gamma_0 + c k_load/(1+(m/K_m)^h)] m

At `c·k_load = 24`, `K_m = 0.25`, `h = 3` the fixed points are
**0.057 / 0.336 / 0.914** — a repressed state, a threshold, and an active one.
Requirements: **h ≥ 2** (h = 1 is monostable at every coupling tested) and
Pol II-driven loss roughly 10–30× the basal writing rate.

This is the bistability that read-write alone could never give, because there
`W(0) = beta > 0` forces a single crossing.

## In the lattice it does not engage

53 hysteresis pairs across `k_load` 2e-5 to 2e-3, `K_m` 0.04 to 0.40, PRC2 8
to 30, on both Ccnd2 and Ccnd1. **Largest gap 0.033** — noise.

The reason is traceable, and it is different for the two genes.

**Ccnd2.** Promoter me3 sits at 0.001 in every condition, so `(m/K_m)^h ≈ 0`
and Pol II is never repressed — the loop is open. Scaling the measured active
field to zero only lifts it to 0.015–0.034. The real cause is geometric: the
Polycomb domain sits over the neighbour Gm56941 at −3.7 kb, so the Ccnd2
promoter is *outside* the PRC2 recruitment zone and cannot hold me3 whatever
Pol II does. At Ccnd2 a promoter toggle is simply not the right model.

**Ccnd1**, where PRC2 nucleates at its own promoter, gets promoter me3 up to
0.09–0.15, and with more PRC2 up to 0.66. But the response to PRC2 dose is
smooth — 8 → 0.14, 16 → 0.42, 30 → 0.66 — a graded dial again, with no gap
from either initial condition.

So the lattice is more dissipative than the mean field: stochastic turnover,
50% dilution every 24 h, and finite PRC2 numbers smear out a switch the
deterministic system supports.

## What did improve

**Ccnd1's profile fit roughly doubled, r = 0.320 → 0.637**, from two changes:
nucleating PRC2 at its own promoter (rather than only at the neighbours) and
letting Pol II produce the turnover instead of a static NDR field. Promoter
me3 is non-zero for the first time (0.14–0.16 at the best fit), so the
bivalent core is partly captured — the failure flagged in
`gene_lattice_2026-10-25.md` is reduced, not yet solved.

## Caveats

- The model's Pol II occupancy trace is not yet trustworthy: the accumulator
  credited only the site where a polymerase ended each step, which biases it
  by step phase. Fixed in the module, but every run above predates the fix, so
  the *eviction dynamics* are sound while the reported occupancy profile is
  not. Re-run before using it.
- `h`, `K_m`, `p_evict` and `k_load` have no measurements behind them here.
  The mean-field window is a statement about the motif, not a prediction with
  numbers attached.
- Three genes, one cell type, for the anticorrelation.

## Next

The honest reading is that mutual antagonism is the right structure — it is in
the data and it improves the fit — but a switch needs the lattice to be less
dissipative than it currently is. The testable moves are a longer cell cycle
(or no dilution), higher local PRC2, or cooperative rather than linear
read-write, each of which the mean field says should widen the window.
