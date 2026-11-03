# Nucleosome-lattice promoter model — build + first results (2026-10-22)

Built outside Antimony as agreed, then fitted an Antimony-expressible surrogate
to it. Code:

- `src/lattice/promoter_lattice.py`  — the lattice (ground truth)
- `src/lattice/mean_field.py`        — single-site chain, fixed points, bifurcation
- `simulations/lattice_experiments.py` — the six experiments below
- `simulations/fit_lattice_surrogate.py` — the two-box surrogate fit

---

## 1. What the model is

Lattice of 51 nucleosomes at 190 bp pitch, asymmetric about the TSS
(-1.5 kb to +8 kb), each carrying one mutually exclusive H3K27 state from
{u, ac, me1, me2, me3}. PRC2 is the only enzyme with read-write feedback; all
antagonists are local and non-spreading:

| term | form | source |
|---|---|---|
| PRC2 write | `P_i = EZH2(1-EZH2i)(a0 + nucl_i + a_rw*(K@me3)_i^n_rw)(1 - w*active_i)` | EED read-write + K4me3/ac blocking the allosteric pocket |
| spreading kernel | geometric, nearest-neighbour dominant | Hodges & Crabtree: bounded domains need LOCAL propagation |
| KDM6B/JMJD3 | `me3->me2->me1` at `k_kdm6b*Gli1*psi_kdm(x)` | existing `k_jmjd3_gli*Gli1`, now spatially resolved |
| p300 | `u->ac` in the core, gated by promoter activity | Newar/Jost 5-state competition on the same residue |
| H3.3 exchange | `any -> u` over the body, ∝ nascent output | existing `k_h33` term |
| replication | each nucleosome independently -> u with p=1/2 at S | replaces `Mk -> Mk/2`, and gives partition noise free |

Time unit hours. me2->me3 is rate-limiting (Newar/Jost).

## 2. Core/flank geometry and the depleted core EMERGES

With CpG-island nucleation plus core antagonists the model produces a depleted
core and a bounded flanking domain:

```
core me3 0.048   flank plateau 0.428   depletion +0.380
boundary edge  2300 bp      20-80% width  1140 bp
```

The depleted zone (~2.3 kb) is **much wider than the 1.1 kb core region I
specified as input**. Core width is therefore an OUTPUT, set by the KDM6B
profile and the boundary dynamics, not an assumption. That is the single
clearest argument for going nucleosome-resolved.

## 3. Mean-field is bistable; the lattice is NOT

Mean field (one homogeneous `m`) gives a clean saddle-node at `a_rw ~ 2.7` and
a wide bistable window (low branch ~0.05, high branch 0.63-0.88).

The lattice, run at identical rates in a homogeneous flank, does not:

| a_rw | MF low | MF high | lattice from m=0 | lattice from m=1 |
|---|---|---|---|---|
| 2.00 | 0.044 | 0.044 | 0.142 | 0.161 |
| 2.75 | 0.045 | 0.634 | 0.159 | 0.216 |
| 3.50 | 0.047 | 0.768 | 0.540 | 0.379 |
| 4.50 | 0.050 | 0.836 | 0.679 | 0.683 |
| 6.00 | 0.057 | 0.884 | 0.824 | 0.830 |

**EQUILIBRATION CORRECTION (2026-10-23).** The 500 h runs in this table were
NOT equilibrated at a_rw = 2.75 and 3.5. Re-run to 3000 h:

| a_rw | 500 h (from 0 / from 1) | 1500 h | 3000 h |
|---|---|---|---|
| 3.50 | 0.593 / 0.363 (gap 0.230) | 0.424 / 0.468 (gap 0.043) | 0.459 / 0.437 (**gap 0.022**) |
| 4.50 | 0.674 / 0.690 (gap 0.016) | 0.656 / 0.684 | 0.691 / 0.705 (gap 0.013) |

So the apparent 0.540-vs-0.379 spread at a_rw=3.5 was a transient, not noise or
hysteresis. **The conclusion is unchanged and in fact better supported**: at
full equilibration both initial conditions converge to the same value, so there
is no bistability. The a_rw=2.75 and 3.5 rows of the table above should be read
as transients; the 4.50 and 6.00 rows were already converged.

Note also **critical slowing down**: equilibration is slow precisely near the
mean-field saddle-node (a_rw ~ 2.7-3.5) and fast away from it — the classic
signature of approaching a bifurcation, and a useful independent confirmation
that the mean-field bifurcation is a real feature of the parameter set even
though the lattice does not inherit the bistability.

Both initial conditions converge: no hysteresis. And the distribution is
**unimodal** (sd 0.15-0.19), so the lattice is not hopping between branches —
it genuinely occupies levels the mean field says cannot exist.

**Why:** spatial coexistence. Observed me3 run lengths are 1.8-2.2x longer than
independent sites at the same mean, so the lattice holds patches of marked and
unmarked nucleosomes (3-6 nucleosomes long) side by side. A single homogeneous
`m` has no variable for that. This is Anink-Groenen's mobile-domain-wall
picture, reproduced here.

Consequence for the current module: **its bistability may be an artefact of
being mean-field.** The permanent-lock latch we found in
`stochastic_amplifier_2026-10-22.md` is the same kind of object.

## 4. Model variants

| variant | core | plateau | depletion | edge (bp) |
|---|---|---|---|---|
| v1 local kernel | 0.048 | 0.428 | +0.380 | 2300 |
| **v2 3D contact 1/d** | **0.149** | **0.177** | **+0.028** | **20** |
| v2b 3D contact 1/d^2 | 0.076 | 0.447 | +0.371 | 2680 |
| v3 collapsed me chain | 0.150 | 0.763 | +0.612 | 1540 |
| v4 no read-write | 0.254 | 0.046 | **-0.208** | none |
| v5 no KDM6B | 0.292 | 0.507 | +0.214 | 1920 |
| v6 no p300 | 0.050 | 0.491 | +0.441 | 2490 |

- **v2 confirms the literature warning.** A 1/|i-j| polymer-contact kernel
  DESTROYS the boundary — the profile flattens to 0.149 vs 0.177 and the edge
  collapses to 20 bp. This is exactly why Option C (3D contacts) was deferred.
  At 1/d^2 the kernel is local enough and the boundary returns.
- **v4 inverts the profile** (depletion -0.208): without read-write the mark
  sits only where nucleation is, i.e. AT the TSS. Read-write is what creates
  the flanking domain at all.
- **v3**: collapsing me1/me2 changes core (0.048 -> 0.150) and plateau
  (0.428 -> 0.763) substantially. The methylation chain is a real kinetic
  delay, not decoration — a 3-state reduction is not adequate.
- **v6**: removing p300 barely moves the core (0.050 vs 0.048).

**CORRECTION (2026-10-23).** I originally read v6 as "KDM6B does essentially all
the core depletion and p300 almost none". The p300 half is right; the KDM6B half
is wrong, and v6 does not test what I said it tested. `w_inh_active` acts on
`k4 + is_AC`, and `k4` is a STATIC mask set to 1 across the core. Setting
`k_ac=0` therefore leaves the core's allosteric block fully intact. Separating
the three channels properly (1500 h, 14 cells):

| variant | core | depletion |
|---|---|---|
| full model | 0.049 | +0.407 |
| no p300 (`k_ac=0`) — the old v6 | 0.049 | +0.419 |
| **no static K4me3 block** | **0.446** | **+0.050** |
| no p300 AND no K4me3 block | 0.475 | -0.024 |
| no KDM6B | 0.341 | +0.126 |
| no KDM6B, no p300, no K4me3 | 0.855 | -0.340 |

**The static K4me3 allosteric block is the single largest contributor to core
depletion** — removing it alone raises the core ~9x (0.049 -> 0.446), more than
removing KDM6B (-> 0.341). p300 contributes nothing measurable.

This matters for the claim in section 2 that the depleted core is an OUTPUT.
Its WIDTH still is (2.3 kb emerged from a 1.1 kb input region), but its
existence and depth are largely IMPOSED by a hard-coded static mask. Until
K4me3 is dynamic — written by SET1/MLL at the CpG island and removed with
transcriptional state — the model should not be cited as showing that the
depleted core self-organises. That is the highest-priority fix to the lattice.

## 5. Gli1 (KDM6B) dose-response = the GNP -> MB axis

| Gli1 | core me3 | flank plateau | edge (bp) |
|---|---|---|---|
| 0.0 | 0.292 | 0.507 | 1920 |
| 0.3 | 0.048 | 0.428 | 2300 |
| 1.0 | 0.007 | 0.402 | 2870 |
| 2.2 | 0.001 | 0.430 | 3440 |
| 3.0 | 0.001 | 0.412 | 3440 |

**Gli1 collapses the core and pushes the boundary outward while leaving the
flank plateau essentially unchanged (0.40-0.51 throughout).** That is the
signature of a strictly local eraser, and it is directly measurable: a
wide-window H3K27me3 readout would barely move between GNP and MB while the
promoter-proximal signal collapses by >100x. Quantification window choice alone
could produce apparently contradictory results.

## 6. Replicative dilution vs cycle length

| cycle | no div | 96 h | 48 h | 36 h | 24 h | 18 h | 12 h |
|---|---|---|---|---|---|---|---|
| flank plateau | 0.428 | 0.326 | 0.287 | 0.210 | 0.143 | 0.102 | 0.064 |
| core | 0.047 | 0.045 | 0.042 | 0.032 | 0.021 | 0.022 | 0.017 |

Monotonic, and steep: a 12 h cycle nearly abolishes the domain. Reproduces
Reinig (memory scales with cycle length) and Alabert (shorter cycle -> less
K27me3, the 2i vs serum/LIF observation) without being fitted to either.

## 7. The Antimony surrogate: does a two-box model suffice?

Two scalars (core, flank) with cross-coupled read-write — i.e. exactly the shape
of the existing `with_proximal_distal` module — fitted to the lattice's Gli1 and
cell-cycle sweeps.

| fit | R2 | RMSE | n_rw |
|---|---|---|---|
| n pinned to the lattice's 2.0 | 0.987 | 0.0215 | 2.0 (fixed) |
| **n free** | **0.994** | **0.0142** | **1.28** |

```
free-n parameters:
  k_eff 0.0588   a0 0.0709   a_rw 2.4834   n_rw 1.2784
  w_cf  0.0779   w_fc 0.1362  psi_c 0.8067  psi_f 0.0000
  ex_c  0.0275   ex_f 0.0000
```

**So yes — "good enough" is achievable and the two-compartment shape survives.**
Four things about HOW it fits, which matter more than the R2:

1. **The effective cooperativity is ~1.3, not the microscopic 2.0.** The free
   fit softens the read-write Hill substantially, and that is what buys the
   improvement. Spatial coexistence averages away the sharp feedback, so a
   mean-field model needs a gentler nonlinearity than the nucleosomes actually
   have. Reading microscopic parameters straight into a mean-field equation is
   the mistake to avoid.
2. **Fitted, the surrogate is MONOSTABLE** — a single flank fixed point at
   0.431 for every Gli1 tested. So it does not invent the bistability that
   mean-field showed in section 3. The resolution: a two-box model is bistable
   when handed the lattice's MICROSCOPIC parameters, and monostable when FITTED
   to lattice BEHAVIOUR. The surrogate is safe provided it is fitted, never
   derived.
3. **KDM6B is core-only in the fit** (psi_c=0.81, psi_f=0.0000), confirming the
   spatial separation is real rather than imposed.
4. Worst residuals are at long cycle (48-96 h) and high Gli1, both within
   ~0.03 absolute.

### Cross-validation — the in-sample R2 is optimistic

| fit on | predict | in-sample RMSE | out-of-sample RMSE | out-of-sample R2 |
|---|---|---|---|---|
| Gli1 axis | cycle axis | 0.0163 | 0.0439 | **+0.883** |
| cycle axis | Gli1 axis | 0.0079 | 0.2114 | **+0.009** |

The second row is a total failure — no better than predicting the mean. It is
largely **structural rather than damning**: along the cycle axis Gli1 is fixed
at 0.3, so `k_kdm6b*Gli1*psi_c` collapses to one lumped constant and psi_c is
simply unidentifiable from that data. But the practical lesson stands:
**this surrogate cannot be calibrated on one perturbation axis and trusted on
another.** Any calibration must span every axis it will be used to predict.

The Gli1 -> cycle direction (R2 = 0.883) is the more meaningful test and passes.

### Held-out perturbation: fit on Gli1 + cycle, predict EZH2i

The decision-relevant test, because EZH2i experiments are already in the 32-target
validation set. The surrogate never saw any EZH2i data.

| EZH2i | core (lattice) | core (surrogate) | plateau (lattice) | plateau (surrogate) |
|---|---|---|---|---|
| 0.00 | 0.047 | 0.043 | 0.427 | 0.444 |
| 0.20 | 0.019 | 0.026 | **0.163** | **0.257** |
| 0.40 | 0.011 | 0.015 | **0.040** | **0.091** |
| 0.60 | 0.004 | 0.009 | 0.011 | 0.040 |
| 0.80 | 0.001 | 0.004 | 0.001 | 0.016 |
| 0.95 | 0.000 | 0.001 | 0.000 | 0.004 |

**HELD-OUT RMSE 0.0326, R2 = +0.925.** The surrogate generalises to an unseen
perturbation — so a fitted two-box module CAN stand in for the lattice inside
v44.

But the error is **systematic, not random**: the surrogate over-predicts mark
retention at every intermediate EZH2i, by up to +0.094 absolute at EZH2i=0.2
(0.257 vs 0.163, a 58% relative overestimate) and +0.051 at 0.4. The lattice's
domain collapses FASTER under partial EZH2 inhibition than any mean-field
version of it, because losing PRC2 breaks read-write locally — patches collapse
and cannot be reseeded by their neighbours — whereas a mean-field domain
declines smoothly.

**Practical consequence:** the surrogate is sound for the qualitative shape of
an EZH2i response but will systematically UNDER-estimate drug effect in the
partial-inhibition regime, which is exactly where a dose-response is measured.
Fitting EZH2i dose-response data with the surrogate would infer the wrong
potency. Use the lattice for anything quantitative about partial EZH2
inhibition.

## 8. Geometry from the existing data — what is and is not possible

- The `chipseq_norm` / `alternateTSS` quantification is **a single 500 bp window
  per TSS** (all 21,599 regions exactly 500 bp). It cannot give a spatial
  profile. No BAMs or bigWigs for these samples are in the Dropbox tree.
- The **epic2 peak calls can**, via a binary-mask metagene. Results:

| sample | n PcG | coverage at TSS | plateau | dip | core width | up | dn |
|---|---|---|---|---|---|---|---|
| EL_GNP1_H3K27me3_1 | 10069 | 0.213 | 0.242 | 0.029 | 1400 | 500 | 900 |
| EL_Tumor1_H3K27me3_1 | 11537 | 0.104 | 0.153 | 0.049 | 2200 | 700 | 1500 |
| Scott_MB_H3K27me3_rep1 | 10348 | 0.969 | 0.978 | 0.009 | 1200 | 300 | 900 |
| P7_H3K27me3_rep5 | 260 | 0.050 | 0.062 | 0.012 | 10600 | 9900 | 700 |

Only the two EL samples are usable: Scott_MB is over-called (97% of the window
is "peak") and P7_rep5 has only 260 targets. For those two:
**core width 1.4-2.2 kb, asymmetric with the core extending further downstream
(900-1500 bp) than upstream (500-700 bp)** — the same direction as the
transcriptional-model kernel fit (lambda_up 248 bp vs lambda_dn 4.7 kb).

**Suggestive but NOT evidence:** EL_Tumor1 (MB, high Gli) has a wider core
(2200 vs 1400 bp) and deeper dip than EL_GNP1, which is the direction the
lattice predicts for raised Gli1/KDM6B. But n=2, same lab only, and the MB
sample also has lower overall coverage (0.104 vs 0.213), which would widen an
apparent core mechanically. It is a hypothesis to test on the cross-tissue data,
not a result.

## 9. Verdict and recommended next steps

1. **The two-box shape survives the redesign** — fitted, it reaches R2 0.994
   in-sample and 0.925 on a held-out perturbation. The interface to v44 does not
   need to change.
2. **It must be FITTED to lattice behaviour, never derived from lattice
   parameters.** Derived, it is bistable and wrong; fitted, it is monostable and
   right. Effective cooperativity is ~1.3, not the microscopic 2.0.
3. **Calibrate across every axis you will predict.** Fitting on the cycle axis
   alone leaves psi_c unidentifiable and predicts the Gli1 axis at R2 = 0.009.
4. **Keep the lattice for partial EZH2 inhibition** and for anything where the
   boundary POSITION (not just the core level) is the observable.
5. Next: close the loop (needs the roadrunner-RHS wrapper from the options doc,
   not a full Antimony port), and fit geometry per gene once profile-resolution
   ChIP exists.

## 10. Equilibration

Added 2026-10-23 after the kymograph showed the flank still rising at 300 h.

Reference promoter geometry re-run to 3000 h from both initial conditions:

| hours | init | core | plateau | edge (bp) |
|---|---|---|---|---|
| 500 | 0.8 | 0.060 | 0.449 | 2110 |
| 500 | 0.0 | 0.045 | 0.484 | 2490 |
| 3000 | 0.8 | 0.043 | 0.495 | 2490 |
| 3000 | 0.0 | 0.040 | 0.497 | 2490 |

**The qualitative results hold.** Core stays depleted, the domain stays bounded,
the boundary stays at ~2.3-2.5 kb, and the two initial conditions agree at
3000 h. But the 500 h numbers used throughout sections 2, 4, 5 and 6 carry a
systematic drift of roughly **+10% on plateau and +18% on edge position** versus
full equilibration. Within-table comparisons are safe (identical protocol
throughout); absolute values are not, and the surrogate fit targets would shift
correspondingly.

Anything quantitative from here should use >=1500 h, or start from a marked
state (init_me3=0.8) which equilibrates considerably faster.

## 11. Caveats

- Rate constants are literature-scaled (Newar/Jost, Hodges & Crabtree) but not
  fitted to these promoters. Absolute me3 levels are not calibrated.
- The working point `a_rw=3.5` sits inside the mean-field bistable window by
  construction, to make the MF-vs-lattice comparison informative. Whether real
  promoters sit there is unknown.
- `kappa_param` here is not numerically comparable to Hodges & Crabtree's: they
  use a 2-state model, this is a 3-step methylation chain, so the transit time
  differs. The qualitative constraint (bounded domains need local propagation)
  is what transfers, and v2 confirms it.
- Ensemble sizes are 16-30 cells; the variant and Gli1 tables carry roughly
  +/-0.02-0.04 on each occupancy.
- One-way coupling only. No feedback from the lattice to the cell cycle yet.
