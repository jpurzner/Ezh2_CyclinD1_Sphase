# Proposed modifications before the noise work starts

Companion to `docs/protein_noise_plan.md`. That document reports what the priors say; this one
is the **list of concrete changes to make, each with the test that was actually run**. Written to
be passed along as-is.

Baseline for everything below: `simulations/validate_v44.py` on current `main`, **28/28 targets**,
17 s runtime. Every test below re-ran that suite via its `VALIDATE_PARAMS` env-var override, so
these are measurements on the real acceptance criterion, not on a side harness.

| # | change | verdict | evidence |
|---|---|---|---|
| **F1** | divide every HH-arm rate constant by 60 | **do it — free** | 28/28 → **28/28**, no flips, max move 3.7% |
| **F2** | *Ccnd1* mRNA lifetime | **largely solved by F1**; sweep the remainder | t½ becomes 52 min, the ActD-consistent end |
| **F3** | "correct" `kDeEZm` to the measured 9.85 h | **do NOT** | 28/28 → **21/28**; it is a real data conflict |
| **F4** | dilution audit classification | do it — cheap now | Mathieson supplies the class-(b) anchors |
| **F5** | declare the HH arm's parameters assumed | do it — documentation only | zero external coverage, §1.4 of the plan |

---

## F1 — The Hedgehog arm is missing a `/60`. Fix it; it costs nothing.

### What is wrong

`HH_MYCN_BLOCK` (`src/build_model_v44_heldt.py:203`) was ported from `build_model_v42_mycn.py`
with its numeric constants intact. **v42's time unit is hours** — its driver scripts run
`t_end=96` and `rescue_t_end=168`, i.e. 4 days and 1 week. **v44's Heldt core is minutes**
(`validate_v44.run()` says "Real-time minutes"; `kDeEZ = 0.0011685/min` gives the documented
EZH2 t½ ≈ 10 h).

So a constant that was authored to mean *1.0 per hour* is now being integrated as *1.0 per
minute*. The Hedgehog arm — and only the Hedgehog arm — runs 60× fast.

The original v42 values are biologically sensible **as hourly rates**. Nobody chose bad numbers;
a unit conversion was skipped. The code comment at line 199 already flags it:

> parameter MAGNITUDES are inherited from v42 and not yet rescaled to the Heldt minute/Cd~0.65
> frame -- deferred to the whole-data calibration pass

| constant | current (read as 1/min) | implied t½ | after `/60` | literature |
|---|---|---|---|---|
| `k_Cd_deg` | 1.0 | 42 s | **41.6 min** | 25 min (classic) / 1.3 h (Schwanhäusser) |
| `k_Cd_mRNA_deg` | 0.8 | 52 s | **52 min** | ~1 h (ActD chase) |
| `k_MYCN_deg` | 1.0 | 42 s | **41.6 min** | ~30 min |
| `k_Gli1_deg` | 0.8 | 52 s | **52 min** | hours |
| `k_Ptch1_deg` | 0.5 | 83 s | **83 min** | hours |

Measured directly on the built model (step `SHH`→0, time the 1/e decay): `Cd` 3.8 min,
`Cd_mRNA` 2.7 min, `Gli1` 3.8 min. The whole cascade equilibrates inside five minutes.

### Why it blocks the noise work specifically

A fluctuation with a 1-minute correlation time, fed into an hour-scale cell-cycle engine, is
averaged to nothing — the engine integrates over hundreds of independent excursions per cycle.
Bursts injected into *Ccnd1* today would give a visibly noisy `Cd_mRNA` trace and **no change in
cycle behaviour at all**, and that null result would be an artefact of the time unit. This is the
one change that has to precede T5/T6 at the CyclinD1 node.

### The fix

Divide the **28 rate constants** by 60. Leave every `K_*` (concentration), `n_*` (Hill exponent)
and dimensionless quantity (`Ptch1_copy_number`, `MYCN_amplification`, `g_smo_Gli1_broken`,
`K_EZH2_repression`) **untouched**.

```
k_Ptch1_basal        k_Ptch1_Gli          k_Ptch1_mRNA_deg     k_Ptch1_translation
k_Ptch1_deg          k_SHH_Ptch_bind      k_SHH_Ptch_release   k_SHH_Ptch_deg
k_Smo_act            k_Smo_inact          k_Gli_rep_to_act     k_Gli_act_to_rep
Vmax_Gli1_tx         k_Gli1_mRNA_deg      k_Gli1_translation   k_Gli1_deg
k_Gli1_auto          k_epi_on             k_epi_off            k_Cd_tx_basal
k_Cd_tx_Gli_max      k_Cd_mRNA_deg        k_MYCN_synth_basal   k_MYCN_synth_Gli
k_MYCN_deg           k_Cd_tx_MYCN         k_Cd_translation     k_Cd_deg
```

Behind a flag `with_realtime_signalling`, default **OFF**, per the brief's opt-in pattern.

### Why it is free, and the test that shows it

Scaling *every* rate in a module by the same factor leaves its steady states **algebraically
invariant** — a steady state is a ratio of rates (`mRNA* = k_tx/k_deg`), and the factor cancels.
Only the transient timescales change. The Hedgehog arm's calibration is entirely steady-state
(dose-responses and fold-changes), so the calibration survives by construction.

That is the argument. Here is the measurement — the full 28-target suite, baseline vs `/60`:

```
target                                       base       resc   ratio
CyclinD1 MB/GNP                            5.6120     5.5060   0.981
CyclinD1 MB+HHi/MB                         0.1300     0.1340   1.031
CyclinD1 GNP+HHi/GNP                       0.2160     0.2160   1.000
Gli1 MB/GNP                                7.1530     7.0580   0.987
MYCN MB/GNP                                2.7110     2.7110   1.000
EZH2 MB/GNP                                1.7210     1.7280   1.004
EZH2i CycD1 fold (GNP)                     1.7160     1.7180   1.001
EZH2 transcript S/G0 (1.8-2.5)             1.5630     1.6210   1.037
Period GNP (~22h)                         22.6890    22.8260   1.006
MB 2N (G0+G1) count% (flow)               68.7060    68.7290   1.000
MB HU S fold                               1.3650     1.3640   0.999
...
VALIDATION SUMMARY: 28/28  ->  28/28      (zero pass/fail flips; largest move 3.7%)
```

All five arrest/cycling calls are preserved exactly (GNP−SHH, GNP+HHi, serum-starve,
MB+CDK4/6i, MB+CDK4/6i+EZH2i all still arrest; GNP+SHH and MB still cycle).

**Two targets moved most, and both moved toward their targets** — `EZH2 transcript S/G0`
1.563 → 1.621 (target 2; this is the worst standing fail in the README) and
`CyclinD1 MB+HHi/MB` 0.130 → 0.134.

> **Correction to `protein_noise_plan.md` §2 R1.** That document said to "budget this as a
> calibration pass, not an edit". That was written before the test and is wrong: it is a
> mechanical edit with no recalibration cost.

### The one real behavioural change

Within-cycle CyclinD1 swing (`Cd_amp`, peak-to-trough over the mean) drops:

| condition | baseline | after `/60` |
|---|---|---|
| GNP | 0.55 | 0.24 |
| MB | 0.72 | 0.27 |
| MB+HHi | 0.40 | 0.19 |

This is the intended effect, not a regression. With a 42-minute half-life CyclinD1 **low-pass
filters** the cell cycle instead of tracking it instantaneously. It is also the mechanism that
makes the noise work possible at all. It should still be eyeballed against any figure that
displays the `Cd` trace over a cycle.

**One thing to check that this pass did not settle:** under `MB+HHi+EZH2i` the baseline harness
found 2 division peaks and the rescaled one found 6. That is the condition `validate_v44.run()`
names as stiff, and the suite passes it either way, so it is most likely an integration artefact
in the baseline rather than a real change — but it is worth confirming before merging.

---

## F2 — *Ccnd1* mRNA lifetime: F1 does most of the work

`docs/protein_noise_plan.md` §1.2 flags this as the single biggest lever: CyclinD1's noise
*amplitude* is robust (CV 0.51–0.91) but its *correlation time* swings 1.7 h → 20 h over the
plausible mRNA half-life, which is the difference between "fluctuates several times per G1" and
"frozen for more than a cycle".

After F1, `k_Cd_mRNA_deg` implies **t½ = 52 min**, which sits at the actinomycin-D-chase end of
the literature and against Schwanhäusser's 14.0 h. So F1 implicitly picks the fast end.

**Recommendation: take 52 min as the default and treat 14 h as the documented upper sensitivity
bound, rather than silently inheriting either.** Concretely, report escape rates as a curve over
t½(mRNA) ∈ [0.5, 14] h. `simulations/analyze_noise_budget.py` already prints the CV/τ band.

If one measurement could be added to this project, this is it: *Ccnd1* mRNA half-life in GNP or a
comparable neural progenitor. Everything downstream is a function of it.

---

## F3 — Do **not** "correct" `kDeEZm`. It is a real conflict, and flagging it is the deliverable.

`protein_noise_plan.md` §1.1 notes `kDeEZm = 0.0098/min` (t½ **1.18 h**) against Schwanhäusser's
measured *Ezh2* mRNA t½ of **9.85 h** — 8.4× apart — and implied this was a discrepancy to close.
**It was tested and it must not be closed.**

Setting `kDeEZm = ln2/(9.85·60) = 0.001173/min`:

| | baseline | kDeEZm→9.85 h | +F1 as well |
|---|---|---|---|
| **suite** | **28/28** | **21/28** | **22/28** |
| EZH2 transcript S/G0 (target 2) | 1.563 ✅ | **0.990** ❌ | **0.988** ❌ |
| EZH2 Palbo mRNA drop (target 0.44) | 0.519 ✅ | **0.995** ❌ | **0.980** ❌ |
| EZH2 G0/cycling (target 0.6) | 0.454 ✅ | **0.966** ❌ | **0.929** ❌ |

The failure is mechanistically transparent: **all** the broken targets are cell-cycle-phase
gradients of EZH2, and they all collapse to ~1.0. An mRNA with a 9.85 h lifetime cannot resolve
phases that last 2–6 h — it averages over them. v44's short `kDeEZm` is load-bearing: it is what
lets EZH2 mRNA track the cycle at all.

So the two cannot both hold for the same system:

* the phase data this repo fits (EZH2 transcript S/G0 ≈ 2, Palbo mRNA drop to 0.44) **require**
  an *Ezh2* mRNA lifetime of roughly an hour;
* Schwanhäusser measures 9.85 h in NIH3T3.

Plausible resolutions, in rough order of likelihood: the GNP/MB transcript measurement is nascent
or intronic (so it reports transcription rate, not mRNA pool, and needs no short lifetime);
*Ezh2* mRNA really is shorter-lived in GNPs than in NIH3T3; or Schwanhäusser's half-life is long
for the same systematic reason its *Ccnd1* value is. **This is a question for the bench, not a
parameter to tune.**

**Good news: the noise conclusion does not depend on which is right.** EZH2's correlation time is
`max(1/g_m, 1/g_p_eff)` and is pinned by the 22 h *protein*, so it is ~13.4 h either way. The CV
only gets smaller under the short-mRNA option (0.26 → 0.12). Either way EZH2 is **slow and
quiet**, and the plan's recommendation — frozen heterogeneity for EZH2, not bursts — stands.

---

## F4 — Dilution audit (brief's T4): now cheap, one item still open

`data/priors/turnover_burst_priors.csv` already carries `growth_corrected` per row. Proposed
classification:

* **(a) effective — already includes dilution; do NOT add `-µ·x`.** Everything inherited from
  Heldt (`kDeE2f`, `kDeP21`, `kDeCe`, `kDeCa`, `kDeE1`, …): fitted in cycling hTert-RPE1.
* **(b) pure degradation — `-µ·x` is a real, missing term.** The Mathieson mouse-neuron values,
  which are post-mitotic and therefore dilution-free by construction:

  | gene | neuron t½ (h) |
  |---|---|
  | `Rb1` | **18.9** |
  | `Cdkn1b` (p27) | 15.6 |
  | `Cdk4` | 28.4 |
  | `Suz12` | 31.2 |
  | `Eed` | 60.7 |
  | `Pcna` | 64.6 |
  | `Cdk1` | 127.1 |

  `Rb1` is the one that matters: at T_cc = 16 h, dilution (µ = 0.043/h) and degradation
  (0.037/h) are **the same size**. That is the Zatulovskiy 2020 premise with a number on it, and
  it is what the size-independent-synthesis option for Rb should be calibrated against.
* **(c) unknown.** Every Schwanhäusser-derived value. Their Methods are said to describe a growth
  correction; **that has not been re-read and verified, and the CSV says `unknown` rather than
  guessing.** Verifying it is the one remaining piece of work in T4 — it is a Methods read, not a
  simulation.

Note that v44's `mass` variable and division event mean dilution is partly already implicit;
the audit must check for double-counting against that, not only against the rate constants.

---

## F5 — State up front that the Hedgehog arm runs on assumed values

Not a code change. From `protein_noise_plan.md` §1.4: `Gli1`, `Gli3`, `Ptch1`, `Smo`, `Ccnb1`,
`Fbxo5` and `Cdkn2d` are absent from **all three** turnover/burst datasets, and BioModels has
**zero** curated entries annotating EZH2, MYCN, GLI2 or GLI3.

The risk is discovering this one parameter at a time and quietly filling gaps with defaults. The
plan should say in advance that every noise parameter on the HH arm is **assumed**, and that the
arm therefore supports sensitivity analysis only — not any claim resting on a parameter value.

---

## Suggested order

1. **F1** — tested, free, unblocks everything else. Do it first.
2. **F4** — cheap, and the Methods check can run in parallel with anything.
3. **F3** — no code change; write the conflict into the model documentation so it is not
   "fixed" by someone later. (Worth a comment at `kDeEZm` in the builder saying *why* it is short.)
4. **F2** — set the default, define the sweep range, then start T5.
5. **F5** — one paragraph in the plan.

## Reproducing the tests

```bash
# baseline
./venv/bin/python simulations/validate_v44.py

# F1: the /60 rescale (writes the override JSON, then runs the real suite)
./venv/bin/python - <<'PY' > /tmp/f1.json
import sys, json, warnings; warnings.filterwarnings("ignore"); sys.path.insert(0,"src")
import tellurium as te
from build_model_v44_heldt import build_model_v44
RATES = "k_Ptch1_basal k_Ptch1_Gli k_Ptch1_mRNA_deg k_Ptch1_translation k_Ptch1_deg \
k_SHH_Ptch_bind k_SHH_Ptch_release k_SHH_Ptch_deg k_Smo_act k_Smo_inact k_Gli_rep_to_act \
k_Gli_act_to_rep Vmax_Gli1_tx k_Gli1_mRNA_deg k_Gli1_translation k_Gli1_deg k_Gli1_auto \
k_epi_on k_epi_off k_Cd_tx_basal k_Cd_tx_Gli_max k_Cd_mRNA_deg k_MYCN_synth_basal \
k_MYCN_synth_Gli k_MYCN_deg k_Cd_tx_MYCN k_Cd_translation k_Cd_deg".split()
r = te.loada(build_model_v44())
print(json.dumps({p: r[p]/60.0 for p in RATES}))
PY
VALIDATE_PARAMS="$(cat /tmp/f1.json)" ./venv/bin/python simulations/validate_v44.py

# F3: the EZH2-mRNA change that must NOT be made
VALIDATE_PARAMS='{"kDeEZm": 0.001173}' ./venv/bin/python simulations/validate_v44.py
```
