# Plan: parameter identifiability, sensitivity and robustness analysis for v44

**Status:** plan, nothing implemented. Written 2026-07-29, in response to an external reviewer's
six-stage proposal (FIM spectrum → sloppy directions → continuation → profiles).

The reviewer's framing is right and the proposal is worth doing. This document adapts it to what
the model actually is today, corrects the parts that are addressed to a model we no longer have,
and reports four measurements made while writing it that change the priority order.

Companion docs: [`DATA_PROVENANCE.md`](DATA_PROVENANCE.md) (the measured/inherited/soft/free tagging
the plan leans on), [`validation_targets.md`](validation_targets.md) (the target spec),
[`EVIDENCE_FRAMEWORK_PROPOSAL.md`](EVIDENCE_FRAMEWORK_PROPOSAL.md) (the provenance-plumbing proposal
this one depends on for σ).

---

## 0. What changed since the proposal was written

The proposal describes v44 as *"22/27 validation, 5 standing failures, ~130 parameters, an AUM mark
module with `f0_mk`/`K_mk`/`del_mk`"*. That was the model as of roughly the `9802d8a` bake. The
current tree (`transient-g0-reframe`, HEAD `eb388d0`) differs enough that several stages need
re-aiming:

| | proposal assumes | actually true now |
|---|---|---|
| validation | 22/27, **5** standing failures | **30/32**, **2** misses — both HU folds, one shared cause |
| parameters | ~130 | **213 global parameters**, 49 floating species |
| mark module | lumped `Mk`, params `k_w_mk, k0_mk, K_mk, n_mk, f0_mk` | **serial me1→me2→me3 chain + PRC2 occupancy**; the named params exist only in the `with_prc2=False` legacy branch and are **absent from the default build** |
| D-cyclin | one (`Cd`) | **two** (`Cd`, `Cd2`) — and the check count is 29 vs 32 depending on whether `w_Cd2 > 0` |
| cell types | one operating point | **P7 GNP / SHH-MB split baked as default** (`k_mu_cki` growth slowdown) |
| CyclinD1 | fast protein | **slow mRNA reservoir** (τ≈16 h) driving a fast protein |
| MYCN | `MYCN_amplification` multiplier | reframed to an additive `MYCN_expr`; **`MYCN_amplification` is now dead code, referenced by no rate law** |

Two consequences worth stating plainly:

- **Stage 6 as written cannot run.** It proposes profiling `f0_mk`, `K_mk`, `del_mk`. Two of those
  three are not in the model. The live analogues are `f0_prc2`, `K_prc2`, `del_mk`.
- **Stage 4 as written is nearly moot.** "Are the 5 standing failures simultaneously satisfiable?"
  now has only 2 failures, both from the same HU fork-coupling, and `simulations/_frontier_hu.py`
  already probes exactly that Pareto front. What replaced it as the interesting question is
  different, and worse — see §1.1.

---

## 1. What to fix before any eigenvalue is meaningful

The reviewer named three prerequisites. Two are real and under-stated; the third — numerical noise
from event handling — turns out to be narrower than feared once measured. A fourth, which the
reviewer had no way to see from outside the repo, is added.

### 1.1 The pass/fail score is hiding the model's actual state — measured

The reviewer is right that `check()` is a step function and unusable as a cost surface. What that
conceals is larger than "22/27 is a coarse number". Converting the current run to the harness's own
band-fraction scale — `|model − target| / (target × tol)`, so ≥100% is a fail:

| target | model | target | band used | harness |
|---|---|---|---|---|
| MB HU S fold | 0.672 | 1.360 | **127%** | FAIL |
| MB HU G2 fold | 0.340 | 0.230 | **119%** | FAIL |
| MB S count% | 22.171 | 15.700 | **98%** | pass |
| CyclinD1 MB/GNP | 4.116 | 5.070 | **94%** | pass |
| CyclinD1 GNP+HHi/GNP | 0.100 | 0.157 | **91%** | pass |
| MB G2+M duration | 3.500 | 2.500 | 73% | pass |
| Period GNP | 17.205 | 22.000 | 73% | pass |
| CyclinD2 MB+HHi/MB | 0.772 | 0.590 | 69% | pass |
| EZH2i CyclinD1 fold | 1.600 | 2.200 | 65% | pass |
| *(16 further targets)* | | | ≤59% | pass |

**7 of 25 continuous targets sit at ≥65% of their tolerance band.** "30/32" reads as a model that
is nearly right with two known structural gaps. The residual picture says something else: the model
is resting against roughly a third of its constraints simultaneously, and the pass count cannot
distinguish "comfortably satisfied" from "pinned to the band edge". Three of the five worst
non-HU residuals — CyclinD1 MB/GNP, CyclinD1 GNP+HHi, EZH2i CyclinD1 fold — are the load-bearing
CyclinD1 numbers the manuscript's central claim rests on.

This also explains a pathology in the search infrastructure. Every optimizer in the repo
(`optimize_*.py`, `v44_wide_search.py`) uses a **band-hinge** loss:
`max(0, |x−t|/t − tol)`, plus `10 × n_fail`. That loss is **exactly zero everywhere inside the
band**. So the search has no gradient over most of the feasible region and no reason to move a
target from 94% of band to 10% of band. The flat directions the reviewer wants to find are partly
being *manufactured by the objective*, not just inherited from the model.

**Action.** Define the residual vector once, in one module, and have every downstream consumer
(optimizers, FIM, profiles) use it:

```
r_i = log(model_i / target_i) / σ_i
```

with the two binary bifurcation checks and the five arrest checks excluded (§1.3). Note the harness
band is *asymmetric in log space* — `|x−t|/t ≤ tol` maps to `[log(1−tol), log(1+tol)]`, which for
tol=0.4 is `[−0.51, +0.34]`. Don't silently symmetrize it; either use `σ_i = (log(1+tol) −
log(1−tol))/2` and say so, or (better) get real σ (§1.2).

### 1.2 σ should not come from the tolerances — and the data to set it is not in the repo

The tolerances are calibration artifacts, not measurement uncertainties, and using them as σ puts
the whole FIM on an arbitrary footing. Two specific problems:

- **The tightest band is on the most ambiguous target.** `CyclinD1 MB/GNP` carries `tol=0.20`, the
  tightest continuous band in the harness. Its target is **5.07 by the Fig 4I normalization but
  7.58 in the raw RNA-seq table** ([`PARAMETERIZATION.md` §A2](PARAMETERIZATION.md)). The
  normalization ambiguity is a factor of 1.5 — **2.5× wider than the tolerance the model is being
  scored against**. Every optimizer run to date has been pulled hard by a number whose own
  definition is looser than its weight.
- **The fold targets rest on n=2.** `MB_Ptch_het` has **n=2** and `MB_GDC0449` **n=2**, against
  `P7_GNP_wt` n=5 and `P7_GNP_wt_GDC0449` n=4. A fold between an n=2 and an n=5 group has a wide
  confidence interval; the ±20% band is almost certainly narrower than the data's own error bar.

**The repo contains only group means** — `data/` holds no counts matrix, and no per-group SD or
DESeq2 `lfcSE` appears in any doc. So σ cannot currently be computed from the repo.

**Action, and it is the highest-leverage input to the entire analysis:** recover, for the ~8
RNA-seq fold targets, the DESeq2 `log2FoldChange` **and `lfcSE`** (or the normalized per-replicate
counts). Everything downstream — which directions are stiff, which residuals matter, what a profile
confidence interval means — is a linear function of σ. Doing the FIM before this is doing it on
invented weights. Where a real σ genuinely cannot be recovered (the IF and microscopy targets),
record the assumed σ explicitly in the target spec rather than inheriting a tolerance.

### 1.3 Drop the checks that carry no independent information

Of the 32 checks, four are not evidence about the model:

| check | why it must be excluded |
|---|---|
| `GNP serum-starve arrest` | **tautological** — the condition sets `k_Cd_translation = 0`, so the arrest is imposed by construction |
| `MB+CDK4/6i+EZH2i no rescue` | a **model prediction** (Fig S9I,J), not a measurement — fitting to it is circular |
| `MB+CDK4/6i arrest (0)` | the data is a **16–18% residual pRb⁺ population that is reversible on washout**; the binary 0-division proxy on one lineage discards both facts |
| `EZH2 Palbo mRNA drop` | measured in **MB**, scored on **GNP** |

Plus the two genuine bifurcation constraints — `GNP+SHH cycles` and `MB cycles` — are binary and
must be excluded from the FIM entirely (the reviewer is right on this) and handled as feasibility
constraints in §2.5.

That leaves **~23 informative continuous residuals** against 213 declared parameters. State that
ratio in the manuscript; it is the honest framing and it is not a weakness, because most of the 213
are inherited and frozen (§2.1).

### 1.4 The discontinuity worry is real but narrower than feared — measured

The reviewer predicted that naive finite differences would produce noise swamping the small
eigenvalues, because of event-time sensitivity (`∂t_event/∂θ` through the state jump). `E_div: at
(MPF > MPF_div)` does reset 21 state variables at a strongly parameter-dependent crossing time, so
the concern is well-founded in principle.

**Measured, it does not happen** (§5.4). Central-difference elasticities over six representative
parameters at relative steps from 0.20 down to 0.005 plateau cleanly by h ≈ 0.05 and are stable
thereafter; successive-step cosine similarity is ≥ 0.99 for five of the six. The reason is that the
validation outputs are not trajectory values — they are settled means, age-reweighted phase
fractions and peak-derived periods over a limit cycle. Those functionals average over the event, and
the event-time term enters them smoothly.

**The one exception is the exception that matters:** `kPhRbCd` at h = 0.20 gives an elasticity of
−100.9 on `Skp2 MB/GNP` and a cosine similarity of **−0.011** against h = 0.10 — a complete sign
flip. That is a bifurcation crossing: a ±20% perturbation of the CyclinD→Rb phosphorylation rate
pushes a condition across the arrest/cycling boundary. At h ≤ 0.10 it is clean.

Two consequences:

- **Use h ≈ 0.02–0.05, central.** Do not use the ±25% one-sided step the recovered
  `sensitivity.py` used (§5.2) — it is larger than the step at which `kPhRbCd` already breaks, so
  that matrix is contaminated for exactly the commitment-adjacent parameters.
- **The Hessian picture really does fail at the saddle-node**, exactly as the reviewer said in
  Stage 5 — and now with a measured demonstration rather than an assertion.

Separately, three *harness-level* discontinuities remain live hazards that the probe above did not
happen to trip:

1. **The retry ladder.** `validate_v44.run()` wraps every simulation in a 3×3×5 retry over
   `max_step ∈ {1e9, 20, 5}` × horizon `∈ {12000, 8000, 6000} min` × `atol ∈ {1e-9 … 1e-5}`. A
   perturbation that trips stiffness **silently shortens the integration horizon from 12000 to 6000
   minutes**, which changes every settled average discontinuously. This is invisible in the output.
2. **`classify()` falls back.** The count-fraction reweighting needs ≥2 MPF peaks; below that it
   silently returns duration-fractions instead — a different quantity.
3. **`find_peaks(prominence=0.15, distance=200)`** — division count changes by ±1 as a peak crosses
   the prominence threshold, stepping period and every count-fraction.

Each of these can still produce the symptom the reviewer predicted, in regions the probe did not
visit. **Fix before measuring:** have `run()` *record* which retry rung it used and
*reject* any evaluation that did not run at the nominal rung; make the `classify()` fallback raise
instead of substituting; and record the peak count alongside every derived quantity so a ±1 step is
detectable rather than absorbed.

There is a fourth, structural, discontinuity that cannot be fixed in the harness: the inherited
Heldt replication-complex disassembly uses **`piecewise(0, Dna < 1, piecewise(1*Rc, Dna > 1,
0.5*Rc))`** for four species. The right-hand side is not C⁰ at `Dna = 1`, so the variational
equations are undefined at a crossing that happens **once per cell cycle**. This is why the AMICI
recommendation does not buy what it promises here (§2.3).

---

## 2. The revised stage plan

### 2.1 Stage 1 — Block the parameters, and freeze most of them

Correct in the proposal; the blocking is already done and just needs to be made machine-readable.
`DATA_PROVENANCE.md` tags every parameter MEASURED / INHERITED / SOFT / FREE, and a
`param_meta.json` with `{module, role, provenance, range, celltype}` for 37 load-bearing parameters
already exists (recovered — see §5).

Of 213 declared globals, the analysis set is far smaller. Remove:

- **Structurally inert parameters — zero sensitivity by construction.** `MYCN_amplification` (no
  rate law references it), `K_EZH2_repression` (textually replaced by the PRC2 occupancy term in
  the default build — the declaration survives orphaned), `k_MYCN_auto`/`K_MYCN_auto`/
  `n_MYCN_auto`/`n_MYCN_Gli`, `k_ptch1_slow_*`/`w_ptch1_slow` (gated off), `w_cilium`/`K_cil`/
  `n_cil` (gated off), `dil_frac` (=0, making the `Mk_replicative_dilution` event a numerical
  no-op that still fires), and the entire legacy `k_w_mk`/`k0_mk`/`K_mk`/`n_mk`/`f0_mk`/`g_mk`/
  `w_ezdir` set (absent from the default build). **Including any of these in a FIM produces a
  spurious exact zero eigenvalue that has nothing to do with the science.**
- **`M_commit` is inert in the default build — verified.** Since `decouple_commit=True` became the
  default (commit `5e4f63e`), `_apply_two_step_rb` drops `*commit_gate` from the p27-clearance term,
  and `commit_gate := mass^n_size/(M_commit^n_size + mass^n_size)` is then **defined as an
  assignment rule but referenced by no rate law**. Doubling `M_commit` from 2.476 to 4.95 changes
  **zero of the 32 validation outputs** (max relative change 0.0). Yet it is baked with a
  15-significant-figure value and is actively searched by `optimize_twostep.py` and
  `optimize_featureC.py`. This is the same pathology as §2.2 in its simplest form: evaluations spent
  on a knob that cannot move anything. Audit for others by checking every assignment rule for
  downstream references before building the FIM.
- **Experimental design inputs, not parameters:** `SHH`, `HHi`, `EZH2i`, `HU`, `Ptch1_copy_number`,
  `MYCN_expr`, `Cd2_expr`, `p16`, `p18`, and the runtime overrides `kPhRbCd=0` (CDK4/6i) and
  `k_Cd_translation=0` (serum starve).
- Note `KSYP21_MB = 0.002` and `P21_DIV_MB = 0.6` are set in every MB condition dict but **equal
  the GNP defaults** — they are currently no-ops. The genuinely cell-type-differentiating knobs are
  only `Ptch1_copy_number`, `MYCN_expr`, `Cd2_expr`, `p16`, `p18`.

Run the analysis twice, as the reviewer suggests — **engine frozen** (the working assumption) and
**engine free** — and report the difference as the inherited share of the uncertainty.

### 2.2 Stage 2 — Analytic degeneracies: the reviewer's list is wrong, and the truth is worse

This is the stage that pays for itself, and it is worth doing *first*, not second.

The proposal offers three predicted null vectors as a pipeline self-test. **All three are wrong as
stated**, for reasons that matter:

| claimed | reality |
|---|---|
| Gli1 rescale absorbed by (`k_Gli1_translation`, `K_Gli_act_CycD`, `K_Gli_act_Gli1`, `K_Gli_MYCN`) | **No such degeneracy exists.** `K_Gli_act_Gli1` is the half-max on `Gli_act`, not on `Gli1`. And `Gli1` enters downstream only as the **sum** `(Gli_act + Gli1)`, where `Gli_act + Gli_rep` is a **conserved moiety pinned at 0.6** by its initial conditions — so `Gli_act` cannot be rescaled with it. At baseline `Gli_act ≈ 0.5` vs `Gli1 ≈ 1.4`, so the degeneracy is not even approximately satisfied. |
| Cd rescale absorbed by 5 named params | Direction right, **membership wrong**. `K_EZH2_repression` is dead code; `K_EZH2_Cd` exists only in the `with_hh=False` placeholder. The **correct 7-member group** is `k_Cd_translation`, `Cd(0)`, `K_CdRb`, `Kez_cd`, **`w_Cd2`**, **`K_Skp2_Cd`**, **`K_g1len`** — the last three were missed, and `w_Cd2` matters because `Cd` enters as `(Cd + w_Cd2·Cd2)`. |
| `k_w_mk`/`k0_mk` against `K_mk` | **Not in the default build.** And even in the legacy branch this is a genuine shape trade-off, not a reparametrization: the `(1 − Mk)` saturation factor pins the `Mk` scale. In the default chain, `Mk` is an **occupancy fraction in [0,1]**, pinned by the literal `1` in `(1 − m1_me − m2_me − Mk)`. It is not an AU species and cannot be rescaled at all. |

The degeneracies that **do** exist are in the sector the model's central claim depends on:

- **PRC2 amplitude:** `{a0_prc2, a_rw_prc2} × λ`, `K_prc2 × λ`, `{kme1, kme2, kme3} ÷ λ` — exactly
  invariant. `PRC2` enters the methylation chain only as a linear prefactor, and `PRC2_rep` only
  through `(PRC2_rep/K_prc2)^n_prc2`.
- **EZH2 protein scale:** `kTlEZ × λ`, `{a0_prc2, a_rw_prc2} ÷ λ`.
- **EZH2 mRNA scale:** `{kEZbas, kEZbas_Cd, kEZE2f} × λ`, `kTlEZ ÷ λ`.
- **CyclinD1 mRNA scale:** `{k_Cd_tx_basal, k_Cd_tx_Gli_max, k_Cd_tx_MYCN} × λ`, `Cd_mRNA(0) × λ`,
  `k_Cd_translation ÷ λ`, `K_tx_mk × λ`.
- **Whole-Gli-sector rescale** (the real version of the reviewer's first claim): `Gli_rep(0)`,
  `Gli_act(0)`, `Gli1(0)`, `K_Gli_Ptch`, `K_Gli_rep_Gli1`, `K_Gli_act_Gli1`, `K_Gli_act_CycD`,
  `K_Gli_rep_CycD`, `K_Gli_MYCN`, `K_Cd2_Gli`, `Vmax_Gli1_tx`, `k_Gli1_auto` all `× λ`;
  `k_jmjd3_gli ÷ λ`.

**This is the finding that most changes what to do next.** The EZH2 → PRC2 → repression chain
carries a **verified 3-dimensional exactly flat subspace**, and it contains
`a0_prc2`, `a_rw_prc2`, `K_prc2`, `kme1–3`, `kTlEZ`, `kEZbas`, `kEZbas_Cd` and `kEZE2f` —
*precisely* the parameters that `optimize_prc2.py`, `optimize_chain.py` and `optimize_cdh1_ezh2.py`
have been searching, and that the `_ts_bake` block reports as freshly optimized.

`optimize_chain.py`'s `SPACE` varies `a0_prc2`, `a_rw_prc2`, `K_prc2`, `kme1`, `kme2`, `kme3`
simultaneously and log-scaled — **all six members of the verified null direction are jointly free
inside its search box.** So a curve along which the objective is provably constant lies entirely
within the region being searched. Those overnight runs have been partly wandering on a manifold
where no evaluation can improve anything. The landed values are one arbitrary point on that curve,
and reporting them as fitted values is not defensible.

**Three of these were verified numerically** (§5.3): the PRC2-amplitude, EZH2-protein-scale and
EZH2-mRNA-scale degeneracies each leave **all 32 validation outputs unchanged to 6×10⁻⁷ or better**
across a factor-4 range of λ, while a control (scaling `a0_prc2` alone) moves outputs by 38–45%. The
CyclinD1-mRNA-scale degeneracy could not be tested through the harness — `_apply_overrides` uses the
regex `name = value;` and **cannot reach species declared in comma-separated `species ... in Cell;`
lists**, which is how `Cd_mRNA(0)` is declared. Fixing that regex is a one-line prerequisite.

**Actions.**
1. **Fix a gauge.** Pin one parameter per exact scaling group before any search or profile — the
   obvious choices are `kTlEZ`, `a0_prc2`, and `k_Cd_translation`. Without this, every profile runs
   away along a unit direction and every optimizer wastes evaluations.
2. Re-express the sector in the identifiable combinations: `kTlEZ·a0_prc2/K_prc2`,
   `a_rw_prc2/a0_prc2`, `kme_i·a0_prc2`. Those are the quantities that can be reported.
3. Use the verified null vectors as the pipeline self-test the reviewer proposed — but with the
   corrected list. If the FIM does not return these at machine zero, the sensitivities are wrong.
4. **Do not run SIAN/StructuralIdentifiability.jl.** The proposal suggests it on the signalling
   submodule. The exact degeneracies above were found by reading rate laws in an afternoon; the
   remaining question is *practical*, not structural, identifiability, which those tools do not
   answer. Skip the Julia dependency.

### 2.3 Stage 3 — FIM spectrum: stay with finite differences, and here is why

The proposal argues for SBML export → AMICI, on the grounds that AMICI handles event sensitivities
correctly. Three reasons that is the wrong call here:

1. **The environment blocks it.** The venv is **Python 3.9**; current AMICI needs ≥3.10. This is a
   new-interpreter, new-toolchain project, not a `pip install`.
2. **The `piecewise` discontinuity defeats it anyway.** AMICI corrects sensitivities across
   *declared events*. The four RC-disassembly `piecewise` terms at `Dna = 1` are **not events** —
   they are discontinuous right-hand sides, and no forward-sensitivity engine handles that without
   reformulating them as explicit events first. If you are going to do that reformulation, you can
   do it in the Antimony source and keep roadrunner.
3. **The outputs are not trajectory values.** Every validation target is a *functional* of a limit
   cycle — a settled mean, a phase count-fraction with age-density reweighting, a peak-derived
   period. Trajectory sensitivities on a limit cycle also grow secularly with time through phase
   drift. `∂x(t)/∂θ` is not the object you want; `∂(functional)/∂θ` is, and finite differences on
   the functional is the direct route to it.

**Recommendation: central finite differences on the existing `VALIDATE_PARAMS` subprocess harness**,
at relative step **h = 0.02–0.05**, after the §1.4 fixes. This is now an evidence-backed
recommendation, not a fallback: the step-size probe (§5.4) shows the functionals are smooth and the
elasticities converge. The cost is small — a measured **44 s** per validation run, so a
central-difference Jacobian over 40 parameters is 80 runs ≈ **8 minutes at 8 workers**.

Re-run the convergence check at three step sizes for any parameter whose elasticity you intend to
interpret, and flag any parameter (like `kPhRbCd`) whose Jacobian column is step-size dependent —
those sit near a bifurcation and belong in §2.5, not in the FIM.

Then: log-parameters, ~23 weighted residuals, `FIM = JᵀJ`, eigendecomposition, spectrum plot, and
read the sloppy eigenvectors for opposite-sign loadings.

**Note that half this infrastructure already exists and was nearly lost.** `simulations/
fig_param_biplot.py` already computes the **SVD of the parameters × targets elasticity matrix** —
a de-facto sloppy-model analysis — and `fig_param_settings.py` renders a per-parameter sensitivity
strip. Both read `sensitivity.json` + `param_meta.json` from a hard-coded path in a **deleted
session scratchpad**, and no committed script generates them. The generator and both JSONs have
been recovered (§5); **commit them** as `simulations/compute_sensitivity.py` before doing anything
else, or this happens again.

### 2.4 Stage 4 — Reachability of the standing misses (much cheaper than proposed)

With 5 failures this was the expensive centrepiece. With **2 failures sharing one cause** it is a
small, well-posed question: **is `(MB HU S fold = 1.36, MB HU G2 fold = 0.23)` jointly reachable?**

`simulations/_frontier_hu.py` already grids `KmHU_fire × vmin_fork × kSyDna` with a
`ProcessPoolExecutor` and probes the S-vs-G2 anti-correlation Pareto front — it is the closest
thing in the repo to a structural-identifiability test. Extend it rather than starting over: add
`KmHU_fork` and `hHU_fork`, and plot the attainable set in the `(S fold, G2 fold)` plane with the
target point marked. If the target lies outside the attainable set, that is a **quantitative
structural claim** for the manuscript — "no parameterization of the current fork-coupling reaches
the measured pair" — replacing "we tried hard".

The reviewer's actual technique (project residual gradients into the stiff subspace, check for
antiparallel corrections) should instead be applied to the **five marginal-but-passing** targets in
§1.1, which is where the unexamined tension now lives — specifically whether `CyclinD1 MB/GNP`
(94% of band, low) and `CyclinD1 GNP+HHi` (91%, low) and `EZH2i CyclinD1 fold` (65%, low) are
pulling in mutually opposed directions. All three are low, all three are CyclinD1, and all three
are load-bearing for the manuscript's central claim. That is the tension worth resolving.

### 2.5 Stage 5 — Continuation: no toolchain, so do reachability maps instead

The reviewer is right that commitment is a saddle-node and that the Hessian picture fails at the
arrest/cycling boundary. But: **PyDSTool and AUTO are not installed and will not install cleanly on
Python 3.9** alongside this stack.

Two realistic options, in order of cost:

- **Cheap and sufficient:** 2D reachability maps. `fig_v44_mitogen_paramspace.py` already draws
  `SHH × f0_mk` and `SHH × p16` planes coloured by division count with the cycle/arrest boundary
  traced. Extend the same idiom to the planes the reviewer names — `(k_Cd_translation, M_commit)`
  and `(K_Gli_act_CycD, n_Gli_act)`. This gives the compensation curve as a boundary in a grid,
  which is what is actually publishable, without a continuation package.
- **If a genuine fold locus is wanted:** `struct/B_timescales.py` already computes the full species
  Jacobian spectrum and does automatic fast/slow splitting, and proposes a slow core. A
  pseudo-arclength tracker with `scipy.optimize.fsolve` on that reduced core is a contained piece of
  work. Do this only if a reviewer asks for the bifurcation diagram specifically.

The `struct/A_loops.py` docstring already flags the relevant caveat: **event edges (`E_div`, mark
dilution) are outside the Jacobian**, so any continuation on the reduced core describes the
between-division dynamics only.

### 2.6 Stage 6 — Profiles, on the claims only (corrected parameter list)

Agreed: do not profile 213, or even 40, parameters. Profile the ones carrying manuscript claims.
Correcting the proposal's list to parameters that exist in the default build:

| claim | profile | note |
|---|---|---|
| repression depth → the 2.2× EZH2i CyclinD1 fold | **`f0_prc2`**, **`K_prc2`** | *not* `f0_mk`/`K_mk`; and `K_prc2` must be gauge-fixed first (§2.2) or the profile runs along the exact degeneracy |
| the ~16 h de-repression half-life | **`del_mk`** | **do not profile yet** — measured max elasticity 0.016 across all targets (§5.4). It is unconstrained by the current target set, so its profile would be flat for a trivial reason. Add a dynamic target first. |
| the Gli-vs-MYCN split of the HHi-resistant CyclinD1 floor | **`k_Cd_tx_MYCN`**, `K_MYCN_Cd`, `n_MYCN_Cd` vs `k_Cd_tx_Gli_max`, `K_Gli_act_CycD` | this is the ~78%/86% decomposition a reviewer will ask about |

Plot co-varying parameters along each profile path — that trace *is* the compensation relationship,
and it is the figure that answers the robustness question directly.

One caution the proposal does not mention: **profile likelihood needs a smooth cost**, so §1.1 and
§1.4 are hard prerequisites. Profiling against the band-hinge loss would produce flat-bottomed
profiles that are artifacts of the objective.

### 2.7 Stage 7 — The reframing (adopt it; it is the strongest part of the proposal)

> Given AU species and ratio-only calibration, individual parameter values were never the
> deliverable. The defensible claim is that the *predictions* project onto stiff directions.

This is right, it is cheap once §2.3 exists, and it directly answers the reviewer question the
current Methods text invites. The three predictions to project:

1. **vismodegib arrests MB** (CyclinD1 crosses the CKI-raised commitment threshold)
2. **EZH2i rescues that arrest** (competitive INK4 brake, overcome by more CyclinD1)
3. **CDK4/6i arrest is not rescuable** (Vmax block cannot be bypassed by transcript)

For each, compute the gradient of the prediction with respect to log-parameters and decompose it in
the FIM eigenbasis. A prediction whose gradient lies in the stiff subspace is well-determined by the
data *regardless of how poorly the individual parameters are determined*. That is a far stronger
robustness statement than any parameter table, and it is the honest one for a model in arbitrary
units.

Note prediction 3 is currently also a **scored check** fit to the model's own output (§1.3). It must
be removed from the residual vector before it can be honestly analysed as a prediction.

---

---

### 2.8 Recommended order of work

The stages are not equally urgent, and the proposal's numbering is not the right order. Sequenced by
value per hour:

| # | do this | why first | cost |
|---|---|---|---|
| 1 | **Fix the gauge** (§2.2): pin `kTlEZ`, `a0_prc2`, `k_Cd_translation` | Stops every future search from wasting evaluations on a provably flat manifold. Nothing else is meaningful until this is done. | minutes |
| 2 | **Fix `_apply_overrides`** to reach comma-declared species (§5.3), and the three harness discontinuity recorders (§1.4) | Prerequisites; both are small | ~1 h |
| 3 | **Commit `compute_sensitivity.py`** (§5.2) with central h=0.02 differences | Restores lost infrastructure, feeds everything downstream | ~1 h + 8 min run |
| 4 | **Recover σ from the RNA-seq replicates** (§1.2) | Every weight downstream depends on it; it is a data request, so start it early | external |
| 5 | **Smooth residual module** (§1.1), retro-fit the optimizers to use it | Removes the band-hinge flat spot that has been shaping every search | ~2 h |
| 6 | **FIM spectrum + prediction projection** (§2.3, §2.7) | The manuscript deliverable | ~1 h |
| 7 | **HU reachability** (§2.4) | Settles the only two standing failures, one way or the other | ~3 h |
| 8 | Profiles (§2.6), 2D maps (§2.5) | Figures, once the above is trustworthy | overnight |

Steps 1–3 are worth doing even if the rest of the programme is dropped.

---

## 3. Compute budget

Everything here is small. The unit cost is one validation run.

| item | runs | wall (8 workers) |
|---|---|---|
| measured baseline validation | 1 | **44 s** |
| degeneracy verification (§2.2) | ~9 | ~2 min |
| FD step-size probe (§2.3) | ~72 | ~12 min |
| central-difference Jacobian, 40 params | 80 | ~8 min |
| Jacobian at 3 step sizes (convergence check) | 240 | ~25 min |
| HU reachability grid, 5 knobs × 8 levels (§2.4) | ~2000 | ~3 h |
| 2D reachability maps, 3 planes × 40×40 (§2.5) | ~4800 | ~7 h |
| profiles, 5 params × 25 points × re-opt (§2.6) | ~5000 | ~8 h |

The whole programme is roughly **one overnight run**, against the many overnight *searches* already
spent. Worst case per evaluation is ~600 s if a parameter set climbs the retry ladder, so budget
generously and — per §1.4 — **reject rather than accept** ladder-climbing evaluations.

---

## 4. What this changes in the manuscript

The current Methods "Parameter calibration" section is also **stale**: it states 22 of 28 targets,
`CyclinD1 MB/GNP ~7×`, a `MYCN amplification = 2.86` multiplier that no longer exists in the model,
and a "~22–23 h" GNP period that `validation_targets.md` documents as a **miscitation** (Nakashima
is 15.9 h). That section needs rewriting regardless of this analysis.

The replacement framing, in one paragraph:

> Species are in arbitrary units and the model is calibrated to ratios, so individual rate constants
> are not identifiable and are not reported as fitted values. We instead characterize which
> *combinations* the data constrain (FIM eigen-spectrum, §), verify the exact scaling degeneracies
> analytically and fix a gauge, and show that the model's three testable predictions project onto
> the well-determined subspace.

---

## 5. Measurements made while writing this plan

Four things were measured rather than assumed. **Scripts and results are in
[`simulations/identifiability/`](../simulations/identifiability/)** — placed in the repo rather than
a scratchpad precisely because §5.2 is a story about losing exactly this kind of artifact.

### 5.1 Baseline validation, timed
30/32, **44 s** wall. Both misses are the HU folds. Residual table in §1.1.

### 5.2 The lost sensitivity generator — recovered
`sensitivity.py`, `sensitivity.json` (37 params × 29 targets) and `param_meta.json` were found in a
prior session's scratchpad and copied to `simulations/identifiability/`. The generator used a **one-sided
+25% finite difference**, which is far too coarse to trust for small eigenvalues — but its SVD is a
useful preview and its `param_meta.json` provenance blocking is directly reusable. **Commit it as
`simulations/compute_sensitivity.py`.**

Preview SVD over the 22 smooth targets: 4 singular values at ~1e-14 (numerical zero, consistent
with the exact degeneracies of §2.2), effective rank 16 of 22 at a 1e-3 relative cut. Stiffest
directions are dominated by `KmHU_fire`, `kSyDna`, `mu`, `k_mu_cki` (the HU/timing block) and by
`kEZbas`, `K_prc2`, `kDeEZ_C1`, `a0_prc2` (the EZH2/PRC2 block). Treat as indicative only.

### 5.3 Exact-degeneracy verification — three confirmed

Each candidate null direction was applied at λ = 2.0 and λ = 0.5 and the full 32-check output
compared to baseline. `simulations/identifiability/degen_test{,2}.py`.

| transformation | max relative change, all 32 outputs | verdict |
|---|---|---|
| **PRC2 amplitude** — `{a0_prc2, a_rw_prc2, K_prc2} × λ`, `{kme1, kme2, kme3} ÷ λ` | 2.4×10⁻⁷ (λ=2), 3.4×10⁻⁷ (λ=0.5) | **exact** |
| **EZH2 protein scale** — `kTlEZ × λ`, `EZH2(0) × λ`, `{a0_prc2, a_rw_prc2} ÷ λ` | **0.0** (both λ) | **exact** |
| **EZH2 mRNA scale** — `{kEZbas, kEZbas_Cd, kEZE2f} × λ`, `EZH2m(0) × λ`, `kTlEZ ÷ λ` | 6.3×10⁻⁷, 3.9×10⁻⁷ | **exact** |
| **control** — `a0_prc2 × λ` alone | 0.45, 0.38 (17 targets moved) | not degenerate — the test has power |
| **CyclinD1 mRNA scale** | 1.2×10⁻², 6.4×10⁻² | **untested** — see below |

Two methodological lessons worth carrying into the implementation:

- **Initial conditions are part of the null vector.** A first pass omitted `EZH2(0)` and `EZH2m(0)`
  and got ~1% residuals, which reads exactly like "degeneracy is only approximate". Adding the
  initial conditions took both to machine zero. Any AU-species rescaling must scale the species'
  initial value too.
- **The CyclinD1-mRNA case is unresolved for a mechanical reason, not a scientific one.**
  `_apply_overrides` matches `name = value;` and `Cd_mRNA` is declared as
  `species MYCN = 0.4, Cd_mRNA = 2.6, Cd in Cell;` — comma-terminated, so the override raises. The
  residual left by omitting `Cd_mRNA(0)` shows the same signature the omitted-IC runs showed for
  EZH2 (worst targets are the two near-threshold arrest conditions, `CyclinD1 GNP+HHi/GNP` and
  `CyclinD1 MB+HHi/MB`, where a transient overshoot perturbs the bistable commitment decision). Fix
  the regex, then re-test.

This is also the argument for the reviewer's pipeline self-test, sharpened: **do not assert null
vectors, verify them.** Of the three the proposal offered, none survived contact with the rate laws;
of the four derived by reading the rate laws, three verified exactly and one is blocked on a
one-line harness bug.

### 5.4 Finite-difference step-size convergence — the outputs are smooth

Six parameters (`kSyDna`, `kPhRbCd`, `k_Cd_tx_Gli_max`, `kDeEZ`, `del_mk`, `K_prc2`), central
differences at relative steps 0.20 / 0.10 / 0.05 / 0.02 / 0.01 / 0.005, 72 validation runs.
`simulations/identifiability/fd_probe.py`, `fd_analyze.py`.

Cosine similarity between the elasticity vector at successive step sizes:

| parameter | 0.20→0.10 | 0.10→0.05 | 0.05→0.02 | 0.02→0.01 | 0.01→0.005 |
|---|---|---|---|---|---|
| `kSyDna` | 0.986 | 1.000 | 1.000 | 1.000 | 0.998 |
| **`kPhRbCd`** | **−0.011** | 0.971 | 0.943 | 0.945 | 0.989 |
| `k_Cd_tx_Gli_max` | 1.000 | 1.000 | 1.000 | 0.999 | 1.000 |
| `kDeEZ` | 1.000 | 1.000 | 0.999 | 0.998 | 0.999 |
| `del_mk` | 0.995 | 0.999 | 1.000 | 0.998 | 0.992 |
| `K_prc2` | 1.000 | 1.000 | 0.999 | 0.999 | 0.999 |

Conclusion in §1.4. The `kPhRbCd` sign flip at h=0.20 is a bifurcation crossing, not numerical noise.

**One further result from the same data, which changes Stage 6.** `del_mk` — the mark turnover rate
carrying the manuscript's ~16 h de-repression half-life — has a **maximum elasticity of 0.016**
across every validation target. Its total influence ranks near the bottom of the recovered
sensitivity matrix (0.05, against 7.5 for `mu`), alongside `k_dil_S` (0.07), `k_h33` (0.06),
`g_prc2` (0.06) and `k_jmjd3_gli` (0.01).

That is not a degeneracy — it is a **gap in the target set**. Every scored target is a settled ratio,
a phase fraction or a period; **none is a time constant**. So the mark module's kinetics are
unconstrained by construction: they come entirely from the chromatin literature (correctly tagged
INHERITED in `DATA_PROVENANCE.md`) and the fit neither confirms nor refutes them. Profiling
`del_mk` against the current targets would produce a flat profile that says nothing.

**The fix is not a better profile, it is a dynamic target.** The de-repression time course after
EZH2i, and the mark-restoration τ after washout, are the measurements that would constrain this
block — and `sim_cdk46i_washout.py` and `fig_ezh2_prc2_kinetics.py` already simulate exactly those
kinetics without scoring them. Add them to the target set before claiming the 16 h half-life is
supported by the model rather than assumed by it.
