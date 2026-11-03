# The lattice model, mathematically

Reference for `src/lattice/particle_lattice.py` and `gene_lattice.py`.

## 1. State

N sites (nucleosomes) at 200 bp pitch. Each site i carries one state

    s_i in {U, ME1, ME2, ME3, AC}

H3K27 can hold one modification at a time, so acetylation and methylation are
mutually exclusive by construction. H3K4me3 is a separate static field
k4_i in [0,1] (measured, not simulated).

## 2. Time stepping

Fixed step dt = 30 s, not Gillespie. Each step applies discrete enzyme events
then per-site field processes. The step is small enough that per-site event
probabilities stay << 1.

## 3. Enzyme binding is a renewal process

Each enzyme molecule alternates search and bound:

    cycle = tau_search + tau_bound
    N_events in dt ~ Poisson( n_mol * dt / cycle )

Each event picks a site from a weight vector w (sum 1) and draws a dwell

    d ~ Exponential(tau_bound)

PRC2 has two binding modes (Tatavosian 2018), chosen per event:

    P(stable) = p_stable * (0.3 + 0.7 * nucl_i)
    d ~ Exp(tau_stable = 10 s)  or  Exp(tau_transient = 0.9 s)

so nucleation sites get longer dwells, not just more visits.

## 4. Catalysis within a dwell

First-order over the dwell:

    P(write) = 1 - exp(-k_i * d)

with

    k_i = (kcat/60) * step_rel[s] / min(step_rel) * allo_i * inhib_i

`step_rel = (8.7, 1.0, 1.3)` for u->me1 : me1->me2 : me2->me3 (Zheng 2014 in
vivo). Dividing by `min(step_rel)` makes kcat refer to the RATE-LIMITING step;
normalising to the fastest step instead penalises the limiting step ~9x twice.

## 5. The two modulating fields

A kernel K couples neighbours, exponential over `allo_range_nuc`, excluding
self, truncated at 6 nucleosomes, row-normalised:

    K_ij ∝ exp( -(|i-j| - 1) / 2.5 ),  K_ii = 0,  K_ij = 0 for |i-j| > 6

**Read-write (allosteric stimulation of EED):**

    allo_i = 1 + (A - 1) * (K m)_i ,    m_j = [s_j == ME3]

**Antagonism by the active marks:**

    a_j = clip(k4_j + [s_j == AC], 0, 1)
    inhib_i = clip(1 - w_inh * (K a)_i, 0, 1)

Both multiply the catalytic rate. Neither changes where PRC2 binds.

## 6. Field processes, per site per step

Each is an independent Bernoulli trial with probability rate*dt:

    AC -> U            at hdac_rate = 0.093/60 /s
    ME3->ME2->ME1->U   at del_mk = 0.0044/h  (iterated low->high so one random
                       draw cannot cascade a site down three levels in a step)
    anything -> U      at k_exch * tx * exch_i, with exch_i a Gaussian on the
                       NDR (centre +150 bp, sigma 450 bp)

The exchange geometry matters more than it looks: at Ccnd2's me3 block
(-5.5 to -2.0 kb) exch_i = 2.4e-07 against 0.994 at the NDR, so turnover acts
essentially only at the promoter core.

## 7. Replication

Every `cycle_h` (24 h), each site independently resets to U with probability
`dil_frac` = 0.5.

## 8. The ensemble

C cells run independently. Cell c draws

    n_prc2^c  ~ lognormal(mean n_prc2, CV 0.4)
    active^c  = active * N(1, 0.25)

Each cell's profile is time-averaged over the second half of the run, then
averaged across cells. That ensemble mean is what gets compared to bulk ChIP.

## 9. Why it is monostable — the mean-field argument

Reduce to one number m, the local me3 fraction. Writing happens on unmarked
sites at a rate stimulated by neighbours; loss is first-order:

    W(m) = beta (1 - m) [1 + (A-1) m]
    L(m) = gamma m
    dm/dt = W - L = -beta(A-1) m^2 + [beta(A-2) - gamma] m + beta

This is a downward parabola with

    f(0) = +beta > 0        f(1) = -gamma < 0

for any positive beta, gamma, A. A curve positive at 0 and negative at 1
crosses zero exactly once in (0,1). **One fixed point, for every parameter
choice.** No amount of allosteric gain produces a second state — which is
exactly what the 102-pair hysteresis sweep found empirically.

The geometric statement: the writing curve starts *above* zero at m = 0,
because unmarked chromatin still gets methylated at the basal rate. A straight
loss line through the origin can only cut such a curve once.

## 10. What would make it bistable

Let me3 *recruit* PRC2, not merely stimulate it. The writing flux then gains a
factor m:

    W(m) = beta m (1 - m) [1 + (A-1) m]
    dm/dt = m * { beta (1-m)[1+(A-1)m] - gamma }

Now m = 0 is always a fixed point (bare chromatin attracts no enzyme, so
nothing starts it), and the bracket can vanish twice. Three fixed points,
outer two stable:

| A | gamma/beta | m=0 | interior roots | |
|---|---|---|---|---|
| 9 | 0.5 | unstable | 0.941 | monostable high |
| 9 | 1.5 | stable | 0.079, 0.797 | **bistable** |
| 9 | 2.0 | stable | 0.180, 0.695 | **bistable** |
| 9 | 4.0 | stable | — | monostable off |
| 25 | 2.0 | stable | 0.046, 0.913 | **bistable** |
| 25 | 4.0 | stable | 0.156, 0.803 | **bistable** |
| 25 | 8.0 | stable | — | monostable off |

Bistability needs gamma > beta (loss outruns basal writing) and the parabola
still clearing zero. It is a window, not a guarantee.

In the code the missing ingredient is one line: `w_prc2` is computed once in
`__init__` from nucleation and never updated, so the enzyme lands in the same
places regardless of what is already marked.
