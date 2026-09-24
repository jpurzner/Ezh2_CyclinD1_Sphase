# Protein-level noise in v44: reconnaissance results and a revised plan

Status: **reconnaissance complete, no model code changed.** This document reports what the
external models and datasets actually say, and proposes how the task list in the hand-off brief
should change in light of it. Nothing in `src/` has been touched; `validate_v44.py` is unaffected.

**Concrete, tested modifications to act on are in [`protein_noise_fixes.md`](protein_noise_fixes.md).**
This document is the evidence; that one is the change list.

Companion artefacts produced by this pass:

| artefact | what it is |
|---|---|
| `.claude/skills/biomodels/` | reusable skill + `bmfetch.py` CLI for the BioModels REST API |
| `models_external/biomodels/` + `SOURCES.md` | 11 reference SBML models, sha256-verified, load-tested |
| `models_external/repos/earm` | EARM (Sorger lab), the ensemble-ODE reference implementation |
| `data/priors/raw/` + `data/priors/SOURCES.md` | Schwanhäusser 2011, Larsson 2019, Mathieson 2018 supplementary tables |
| `data/priors/turnover_burst_priors.csv` | T3 prior table, 528 rows, every value sourced, 230 explicit NA |
| `simulations/build_turnover_priors.py` | regenerates the CSV from the raw files (never imputes) |
| `simulations/analyze_noise_budget.py` | predicted CV and correlation time per species vs the model's own timescales |

---

## 1. Four findings that change the plan

### 1.1 The signalling arm is 100–1000× faster than the biology, so it cannot carry noise

`HH_MYCN_BLOCK` in `src/build_model_v44_heldt.py:195` says so itself:

> parameter MAGNITUDES are inherited from v42 and not yet rescaled to the Heldt minute/Cd~0.65
> frame -- deferred to the whole-data calibration pass

The Heldt core runs in **real-time minutes** (`validate_v44.run()` is explicit about this, and
`kDeEZ = 0.0011685/min` gives the documented EZH2 t½ ≈ 10 h). The bolted-on v42 constants are
interpreted in the same minutes:

| constant | value (1/min) | implied t½ | literature | model relaxation vs biology |
|---|---|---|---|---|
| `k_Cd_deg` (CyclinD1 protein) | 1.0 | **42 s** | 1.3 h (Schwanhäusser, mouse) | **105× too fast** |
| `k_Cd_mRNA_deg` (*Ccnd1* mRNA) | 0.8 | **52 s** | 14 h (Schwanhäusser) / ~1 h (ActD chase) | **970× / 70× too fast** |
| `k_MYCN_deg` | 1.0 | 42 s | ~30 min (MYCN is famously short-lived) | ~40× too fast |
| `k_Gli1_deg`, `k_Gli1_mRNA_deg` | 0.8 | 52 s | hours | ~100× too fast |
| `k_Ptch1_deg` | 0.5 | 83 s | hours | ~50× too fast |

Measured directly rather than inferred: stepping `SHH` to 0 in the built model and timing the
1/e decay gives `Cd` **3.8 min**, `Cd_mRNA` **2.7 min**, `Gli1` **3.8 min**. The whole signalling
cascade equilibrates in under five minutes.

This is a defensible modelling choice as long as the module is a quasi-static mitogen relay —
which is exactly how the code describes it. **It stops being defensible the moment the module is
asked to carry noise.** Fluctuations with a 1-minute correlation time, fed into a cell-cycle
engine with hour-scale dynamics, are averaged to nothing: the engine integrates over hundreds of
independent excursions per cell cycle. Injecting bursts into `Ccnd1` in the model as it stands
would produce a visually noisy `Cd_mRNA` trace and **no change whatsoever in cycle behaviour**,
and the null result would be an artefact of the timescale, not biology.

**Consequence for the brief: T6 (explicit bursts on *Ccnd1*) cannot be done first.** It depends on
a prior rescaling of the HH/MYCN arm to real minutes, which is a structural change requiring
recalibration — the "whole-data calibration pass" the code comment already defers to.

### 1.2 For CyclinD1 the noise *amplitude* is well determined and the *timescale* is not

From `analyze_noise_budget.py`, using
`CV_p² ≈ 1/⟨p⟩ + (1+B)/⟨m⟩ · g_p_eff/(g_p_eff+g_m)` with dilution at T_cc = 16 h:

```
Ccnd1   ⟨p⟩ = 89,632   ⟨m⟩ = 8.6   B = 6.76   t½ protein = 1.3 h
        t½ mRNA:   0.5 h    1 h     2 h     4 h     9 h     14 h
        CV_p        0.51    0.64    0.75    0.83    0.89    0.91
        tau_c/h     1.7     1.7     2.9     5.8     13.0    20.2
```

⟨p⟩ and t½(protein) are Schwanhäusser; `B` and ⟨m⟩ are Larsson (mouse fibroblast, both alleles).

The amplitude is robust: **CV ≈ 0.5–0.9 across the entire plausible range** — CyclinD1 is a
genuinely noisy protein however you resolve the mRNA lifetime. The correlation time is not: it
swings by a factor of twelve, from "fluctuates several times per G1" to "essentially frozen for
longer than a cell cycle". Those are opposite predictions for the one question the noise work
exists to answer — whether fluctuations drive transitions across the G0 basin.

The brief's §2.2 table assumed a *Ccnd1* mRNA half-life of "~1 h (placeholder)". The only direct
measurement we now hold says **14.0 h** (replicates 10.5 h and 18.7 h). Schwanhäusser's mRNA
half-lives are known to run long, so this should not simply be adopted — but it should not be
ignored either, and it is a 12× lever on the headline conclusion.

**Consequence: the mRNA half-life of *Ccnd1* is the highest-value single measurement in this
project.** Until it is pinned, every escape rate is a function of an unmeasured parameter, and
the honest output is a curve over t½(mRNA), not a number.

### 1.3 EZH2 is the mirror image — slow and quiet — and the model already has it right

```
Ezh2    ⟨p⟩ = 2,156   ⟨m⟩ = 30.0   B = 2.77   t½ protein = 22.0 h   t½ mRNA = 9.85 h
        CV_p = 0.26      tau_c = 14.2 h
        t½ mRNA 0.5→9.85 h:  CV 0.08 → 0.25,  tau_c pinned at 13.4 h throughout
```

Because the EZH2 protein is long-lived, its correlation time is set by the protein, not the mRNA,
so it is **insensitive to the one parameter that wrecks the CyclinD1 estimate**. EZH2 fluctuations
are slow (≈ 13 h ≈ one cell cycle) and small (CV ≤ 0.28).

And a genuine independent validation: v44's `kDeEZ` implies a relaxation time of **14.26 h**; the
Schwanhäusser-derived biological value is **14.2 h**, a ratio of 0.9×. That constant was fitted in
this repo to HU-arrest and G0-withdrawal data with no reference to proteomics, and it landed on
the measured value. (The same is not true of `kDeEZm`: model 1.7 h vs biology 14.1 h, 8.4× fast.)

**Consequence: EZH2 does not need a burst module.** For a species with τ ≈ one generation, a cell
essentially keeps its EZH2 rank for the whole cycle — which is frozen heterogeneity, i.e. the
brief's **option 2** (per-cell sampled parameters), not option 3 or 4. This is worth stating
positively, because it means the epigenetic brake is a *stable* per-cell property while the
mitogenic drive is the noisy one. That asymmetry is a model prediction with a clear experimental
reading, and it falls out of the priors rather than being assumed.

### 1.4 Data coverage, not method, is the binding constraint — and it is worst where the biology is novel

| dataset | genes matched | of |
|---|---|---|
| Schwanhäusser 2011 (copy number, half-lives, rates) | 13 | 39 |
| Larsson 2019 (burst size + frequency) | 27 | 39 |
| Mathieson 2018 (post-mitotic protein half-life) | 8 | 33 |

Absent from **all three**: `Gli1`, `Gli3`, `Ptch1`, `Smo`, `Ccnb1`, `Fbxo5`, `Cdkn2d`.

The BioModels UniProt probe says the same thing from the other direction. Of 29 proteins probed:

* well covered — `CCNE1` 22 models, `E2F1` 17, `CDK2` 17, `RB1` 14, `CCNB1` 12, `CDK1` 12;
* **`EZH2` 0 models. `MYCN` 0. `GLI2` 0. `GLI3` 0. `GLI1` 1. `SMO` 1. `PTCH1` 2.**

The models with the broadest overlap with v44 are Gerard2009 (`BIOMD0000000730`, 16/29),
Iwamoto2010 (`…939`, 13/29) and Csikász-Nagy2006 (`…1044`, 11/29) — all cell-cycle cores.

So: **the cell-cycle half of v44 sits in a dense, well-parameterised literature; the EZH2/Hedgehog
half — the part that is actually novel in this paper — has essentially no external support at
all.** Every noise parameter on that arm will be "assumed", and the plan should say so in advance
rather than discovering it per-parameter. It also means the noise work cannot be used to make
claims about the HH arm that rest on parameter values; it can only be used for sensitivity.

### 1.5 A bonus: the dilution audit now has a clean class-(b) anchor

Mathieson's `Mouse Neurons` columns are **primary post-mitotic neurons**. They do not divide, so
there is no dilution term and the fitted half-life is unambiguously pure degradation. These are
the only certain **class (b)** values available for T4:

| gene | neuron t½ (h) | note |
|---|---|---|
| `Rb1` | **18.9** | at T_cc = 16 h, dilution (µ = 0.043/h) ≈ degradation (0.037/h) |
| `Cdkn1b` (p27) | 15.6 | |
| `Cdk4` | 28.4 | |
| `Suz12` | 31.2 | |
| `Eed` | 60.7 | |
| `Pcna` | 64.6 | |
| `Cdk1` | 127.1 | |

The `Rb1` number is the useful one: it puts a measured value on the Zatulovskiy 2020 premise that
Rb is diluted by growth. Dilution and degradation are the same size for Rb in a 16 h cycle, so
`-µ·x` on Rb is a real term, not a rounding correction.

---

## 2. Proposed revision to the task list

The brief's ordering (T1 → T8) front-loads infrastructure and defers the units question to T4.
The findings above invert that: **the units problem gates the two tasks the project exists for.**

### Do first

**R1 — Rescale the HH/MYCN/CyclinD1 arm to real minutes.** *(new; blocks T5 and T6 at the CyclinD1
node)* Divide the 28 HH-block rate constants by 60, behind a flag `with_realtime_signalling`
defaulting **OFF**.

**This has since been tested and it is free.** The root cause is narrower than first described:
the constants were ported from `build_model_v42_mycn.py`, whose time unit is hours (its drivers
run `t_end=96`/`168`), into v44's minute-frame core — so they are simply missing a `/60`, and the
original values are biologically sensible as hourly rates. Because scaling every rate in a module
by one factor leaves its steady states algebraically invariant, and the arm's calibration is
entirely steady-state, `validate_v44.py` goes **28/28 → 28/28** with no pass/fail flips and a
largest per-target move of 3.7%. See `docs/protein_noise_fixes.md` §F1 for the full per-target
diff and the one real behavioural change (within-cycle CyclinD1 amplitude drops ~2.5×, which is
the intended effect).

*An earlier revision of this section said to "budget this as a calibration pass, not an edit".
That predated the test and is wrong — it is a mechanical edit.*

**R2 — Settle the *Ccnd1* mRNA half-life, or commit to reporting a curve over it.** *(new)* This is
the 12× lever. Either find a direct measurement in GNP/MB or a comparable neural progenitor, or
make every downstream escape-rate result explicitly a function of it.

**T4 — dilution audit** — promote to third. It is cheap now: `data/priors/turnover_burst_priors.csv`
already carries `growth_corrected` per row, Mathieson supplies the class-(b) anchors, and the
Heldt-inherited constants are class (a) by construction (fitted in cycling RPE1). The one open
item is verifying Schwanhäusser's Methods, currently recorded as `unknown` rather than guessed.

### Then, with the timescales fixed

**T5 — drift module (option 3)** — unchanged; keep the exact OU discretisation and the
segmented-integration scheme as specified, they are sound. Note what can actually be
parameterised, though: the full `CV_p²` formula needs ⟨p⟩, ⟨m⟩, `B` **and both half-lives**, and
only four v44 species have all of them — **`Ccnd1`, `Ccna2`, `Ezh2`, `Skp2`** (plus `Cdk4`, which
is a parameter not a species). `Cdkn1a`, `Ccne1`, `Trp53`, `E2f1`, `Cdc20`, `Cdk1` and `Dtl` have
Larsson burst data but no matching half-life, so for those the drift amplitude is constrained and
the correlation time is not. Everything else on the HH arm is assumed.

**T6 — burst module (option 4)** — **narrow the target list.** The brief nominates *Ccnd1*,
*Cdkn1a* and *Mycn*. On the priors:

* ***Cdkn1a* is the strongest candidate, ahead of *Ccnd1***: burst size **B = 20.0** at mean
  expression 26 (fibroblast) — by far the burstiest gene in our set, and Barr 2017
  (`BIOMD0000000660`, downloaded) already models p21 heterogeneity as discrete stochastic events.
* ***Ccnd1*** (B = 6.8) is second, and **only after R1**.
* ***Mycn*** has burst data in ES cells only (B = 18.5) and **no protein half-life in any of the
  three datasets**. It cannot be parameterised; treat it as assumed-value sensitivity only.
* **Drop *Ezh2* from consideration** — §1.3.

**T2 — stochastic Heldt reproduction** — keep, but note it is a toolchain check, not a result.
The `hasOnlySubstanceUnits` / compartment-volume warning in the brief is well aimed: the same
class of bug (`Mk` as a compartment species rather than a parameter with a rate rule) already bit
this repo once.

### Keep as specified

**T1** — done (this pass). **T7** (division), **T8** (ensemble runner) — unchanged.

### Reconsider

**T3 as originally scoped is not fully achievable.** "Every value has provenance" is satisfied;
"for each gene, extract where available" yields 199 of 429 numeric slots filled, and the gaps are
concentrated in the HH arm. The CSV records them as explicit `NA`, per the brief's instruction
never to impute silently. The plan should state up front that the HH arm will run on assumed
values, rather than treating each missing cell as a lookup still to do.

---

## 3. Two things the brief listed that resolved differently

* **`MODEL1703030000` is not a separate model.** It redirects to `BIOMD0000000700`, and the
  downloaded file is byte-identical (md5 `bd4e342d…`) to `models_external/heldt2018_BIOMD700.xml`
  already in this repo. v44's core is the published Heldt model, unmodified — confirmed, not
  assumed.
* **Geva-Zatorsky Models 3 and 6 cannot be run.** `BIOMD0000000154` and `BIOMD0000000155` are
  delay differential equation models; libRoadRunner refuses them. Models 1, 2, 4 and 5 integrate.
  The brief wanted these as precedent for low-frequency noise on production rates — Models 1/2/4/5
  serve that purpose, so nothing is lost, but the DDE pair should not be chased.

Also resolved: **HARISSA** is `harissa` 3.1.2 on PyPI, github.com/ulysseherbach/harissa, by Ulysse
Herbach as the brief guessed. Not installed; it is a design reference for the PDMP module.

Two models **not** in the brief were added because the UniProt probe surfaced them and both model
the transition the noise work is meant to drive: **Pandey2018** (`BIOMD0000000954`, reversible
quiescence↔proliferation) and **Yan2012** (`BIOMD0000000720`, Rb–E2F bistability).

---

## 4. Reproducing this

```bash
python3 simulations/build_turnover_priors.py     # regenerate the prior CSV from data/priors/raw/
python3 simulations/analyze_noise_budget.py      # the CV / correlation-time tables above
python3 simulations/analyze_noise_budget.py --tcc 24   # same, for a slower cycle
```

Both are stdlib + `xlrd`/`openpyxl` only and do not import the model, so they run without the venv
and cannot perturb `validate_v44.py`.
