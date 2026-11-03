# Gli1, KDM6B and the peripheral H3K27me3 reservoir

2026-10-25. SLURM 12386333, 12386344, 12386366, 12386385.
Code: `polii_lattice.py` (Gli1-targeted KDM6B, `withdraw_gli1()`).

## The mechanism being tested

JP's reading: Gli1 binds the proximal promoter and the enhancers and, via its
co-IP partner KDM6B, keeps H3K27me3 off the proximal region while allowing it
to persist further out. The gene is therefore highly expressed *despite*
carrying the mark. Remove the mitogen, Gli1 leaves, and the peripheral
reservoir re-encroaches quickly — giving **fast off-kinetics**. A promoter
whose me3 had been fully purged has nothing to draw on and should shut down
**slowly**.

This is the lattice-level version of the reservoir in
`two_compartment_prc2_2026-08-09.md`.

## Implementation

`GeneGeom` now carries the measured Gli1 profile; at Ccnd1 it peaks at
**−1.5, −1.3, +2.1 and +2.3 kb**. KDM6B is recruited to those sites rather
than spread over the core, and `withdraw_gli1()` removes the targeting.

## Two corrections made along the way

**The first withdrawal run was not a valid test.** Pre-withdrawal rates
differed 7–25× between arms (1.06/h purged against 0.04–0.16/h with a
reservoir), so the reservoir arms were already repressed and never occupied
the marked-but-expressed state the mechanism requires.

**The periphery was mis-defined.** Scoring the "reservoir" as 1.5–6 kb
includes the Gli1 sites at ±1.3–2.3 kb, so Gli1-recruited KDM6B was erasing
the region being measured. The reservoir has to be scored *beyond* the Gli1
sites, at 3–8 kb.

## The marked-but-expressed state exists, but needs far more KDM6B

Sweeping KDM6B against PRC2 (regions scored correctly):

| KDM6B | PRC2 | rate /h | core <1 kb | Gli1 1–3 kb | reservoir 3–8 kb |
|---|---|---|---|---|---|
| 1 | 24 | 0.05 | 0.558 | 0.428 | 0.305 |
| 3 | 24 | 0.12 | 0.394 | 0.351 | 0.237 |
| 10 | 24 | 0.50 | 0.161 | 0.254 | 0.168 |
| **30** | **24** | **0.82** | **0.057** | 0.194 | **0.115** |
| 30 | 8 | 1.01 | 0.013 | 0.074 | 0.006 |

At KDM6B = 30 and PRC2 = 24 the gene runs at 82% of maximum rate with a clear
core and a marked periphery — exactly the configuration described. Note this
needs **100× the committed KDM6B default of 0.3 molecules per domain**. Either
that default is far too low, or local Gli1-templated recruitment concentrates
KDM6B far above its genome-average abundance.

## The result

Both arms run an identical phase 1 to that state (0.774 initiations/h). At
withdrawal both lose Gli1 and keep PRC2 at 24; they differ **only** in whether
me3 is erased at t = 0.

| | reservoir | purged |
|---|---|---|
| pre-withdrawal rate | 0.774 /h | 0.774 /h |
| peripheral me3 at t=0 | 0.095 | 0.000 |
| **t½ of shutdown** | **10 h** | **20 h** |
| transcripts in the first 50 h | 5.9 | 12.8 |

Bootstrap over cells: purged is **+9.7 h slower** (95% CI +5.0 to +10.0),
slower in 99.9% of resamples.

So the prediction holds: **the peripheral reservoir halves the shutdown time
and more than halves the transcript made after the mitogen is gone.** The
third panel shows the mechanism — proximal me3 is already at 0.135 by 5 h in
the reservoir arm, while the purged promoter starts at 0.078 and takes 15 h to
catch up. Both end at the same place; the reservoir simply gets there first.

This also makes the broader point concrete: **H3K27me3 does not preclude high
expression as long as it is held far enough from the proximal promoter.** The
marked-but-expressed state runs at 82% of the unmarked rate.

## Caveats

- One gene, one parameter set. `k_load`, `p_evict`, `K_m`, `h` and the KDM6B
  abundance are all unconstrained by measurement.
- "Purging" is implemented as erasing every me3 at the instant of withdrawal,
  which is cleaner than any real EZH2i pre-treatment would be.
- The absolute timescale is set by the model's turnover and dilution, so the
  **ratio** (2×) is the claim, not the 10 h.
