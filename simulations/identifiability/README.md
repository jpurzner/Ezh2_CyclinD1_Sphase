# Identifiability probes — evidence behind `docs/PARAMETER_IDENTIFIABILITY_PLAN.md`

Scripts and results from the 2026-07-29 scoping pass. Nothing here is part of the figure pipeline;
these exist so the measurements in the plan are reproducible and so the artifacts are not lost to a
deleted scratchpad a second time (see plan §5.2).

Each script drives `simulations/validate_v44.py` as a subprocess via `VALIDATE_PARAMS`, the same
interface every optimizer uses. One validation run ≈ 44 s.

| file | what it measures | result |
|---|---|---|
| `degen_test.py` → `degen_test.json` | Are the claimed PRC2/EZH2 scaling degeneracies exact? λ = 2.0 and 0.5, all 32 outputs vs baseline, plus a non-degenerate control | PRC2 amplitude **exact** (2.4e-7); control moves outputs 38–45% |
| `degen_test2.py` → `degen_test2.json` | Same, with species **initial conditions** included in the transformation | EZH2 protein scale **exact (0.0)**; EZH2 mRNA scale **exact (6.3e-7)** |
| `degen_test3.py` | CyclinD1 mRNA scale, without `Cd_mRNA(0)` | 1.2e-2 / 6.4e-2 — **inconclusive**; `_apply_overrides` cannot reach comma-declared species, so the initial condition could not be scaled |
| `fd_probe.py` → `fd_probe.json` | Central-difference elasticities at relative steps 0.20 … 0.005 for 6 parameters (72 runs) | Converged by h ≈ 0.05; `kPhRbCd` sign-flips at h = 0.20 (bifurcation crossing) |
| `fd_analyze.py` | Reads `fd_probe.json`, prints convergence tables and successive-step cosine similarity | — |
| `val_base.json` | Baseline validation, 30/32, for diffing | — |

## Recovered from a deleted scratchpad

`compute_sensitivity_RECOVERED.py`, `sensitivity.json`, `param_meta.json`, `param_cats.json`.

`simulations/fig_param_settings.py` and `simulations/fig_param_biplot.py` both read
`sensitivity.json` and `param_meta.json` from a hard-coded path in a session scratchpad that no
longer exists, and **no committed script generated them**. They were found in a stale scratchpad and
copied here. Both figures are currently unreproducible without them.

Caveat on `sensitivity.json`: it was computed with a **one-sided +25% finite difference** over 37
parameters × 29 targets. Per `fd_probe`, 25% is larger than the step at which commitment-adjacent
parameters cross a bifurcation, so the matrix is contaminated for exactly those columns. Treat it as
indicative, and regenerate with central h = 0.02 differences (plan §2.3, step 3).

## `M_commit` inertness check

`val_mcommit.json` — validation output with `M_commit` doubled (2.476 → 4.95). Diff against
`val_base.json` is **exactly zero on all 32 checks**. With `decouple_commit=True` (default since
`5e4f63e`), `commit_gate` is defined as an assignment rule but referenced by no rate law, so
`M_commit` cannot affect anything. It is nonetheless searched by `optimize_twostep.py` and
`optimize_featureC.py`. Reproduce with:

```bash
VALIDATE_PARAMS='{"M_commit":4.95}' VALIDATE_RESULTS=/tmp/m.json TWO_STEP_RB=1 H3K27_CHAIN=1 EZH2_CONC=0 ./venv/bin/python simulations/validate_v44.py
```
