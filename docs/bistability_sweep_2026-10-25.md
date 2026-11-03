# Is the lattice bistable? No — across everything tested

2026-10-25. SLURM 12382824, 12382996, 12383172, 12383267. Ccnd2 geometry,
PRC2 = 8, 60 cells per condition, 300 h. Bistability tested by **hysteresis**:
each parameter pair is run from a fully marked domain and from a bare one. If
they settle apart, two states coexist.

## The answer

Every condition gives the same level from both starts. The largest gap
anywhere is 0.043, and signs flip at random, so it is run-to-run noise.

| sweep | conditions | largest |gap| |
|---|---|---|
| allo 1–25 x k_exch 0.15–2.4/h, NDR-local turnover, +replication | 35 pairs | 0.019 |
| allo 1–25 x k_exch 0.05–2.5/h, **uniform** turnover, +replication | 35 pairs | 0.011 |
| allo 6–50 x recruitment range 0–100% non-local, uniform turnover, +replication | 16 pairs | 0.043 |
| same grid, **replication off** | 16 pairs | 0.043 |

Bimodality coefficients of the per-cell domain average run 0.19–0.46 against
a 0.555 threshold; the distributions are unimodal. Response to allosteric gain
is smooth and monotone (0.25 → 0.76 as allo goes 1 → 25), which is a graded
dial, not a switch.

**So genuine cell-to-cell mixtures cannot arise from bistability in this
model.** The intermediate levels in bulk ChIP are true intermediates.

## A flaw in the first sweep, and the fix

The first sweep varied `k_exch` while exchange was NDR-local, so the domain
never experienced it: the exchange weight at Ccnd2's me3 block (−5.5 to −2.0
kb) is 2.4e−07 against 0.994 at the NDR. Turnover rows were therefore flat by
construction (allo=9 gave 0.567 → 0.574 across a 16x range of k_exch) and the
sweep said nothing about turnover. Repeating it with uniform exchange fixed
that, and showed turnover is destructive rather than bistability-generating:
at k_exch >= 0.4/h the domain is wiped (0.004–0.011), and only below ~0.15/h
does anything survive.

That also means the committed `k_exch = 0.60/h` only works *because* it is
NDR-local. Applied uniformly it would erase every domain.

## Why it is monostable — a structural reason, not a parameter one

Two things in `particle_lattice.py`:

```python
def _allo(self):
    return 1.0 + (self.p.allo_max - 1.0) * (self.K @ (self.state == ME3))
...
self.w_prc2 = w / w.sum()      # set once in __init__, never updated
```

1. **The read-write feedback is linear** in the local me3 fraction. Linear
   positive feedback on a bounded variable gives a graded response. A switch
   needs supra-linear (cooperative) feedback — several me3 neighbours required
   to engage EED, not a proportional boost.
2. **me3 boosts catalysis but does not recruit PRC2.** `w_prc2` is fixed by
   nucleation at initialisation and never responds to the mark. Read-write in
   the Dodd/Sneppen sense means existing me3 *brings in more enzyme*; here the
   enzyme lands wherever nucleation says regardless of what is already there,
   so the loop that would close the feedback is absent.

Pushing allosteric gain to 50x (4–7x above the measured 7–12x) and making
recruitment fully domain-wide does broaden the per-cell spread — at allo=50,
f=1 about 20% of cells sit above 0.6 and 27% below 0.2 — but it stays a wide
unimodal distribution, not two states.

## What to change if bistability is wanted

Make `w_prc2` depend on the current me3 state, so that me3 recruits PRC2 as
well as stimulating it, and/or make `_allo` cooperative (a Hill term in local
me3 rather than a linear one). Either closes the feedback loop properly. Both
are departures from what the in-cis measurements support, so they should be
presented as hypotheses to test rather than corrections.

Until then the model's position is: **the H3K27me3 domain is a graded,
monostable structure, and the bulk ChIP profile reflects what every cell is
doing rather than a mixture of marked and unmarked cells.**
