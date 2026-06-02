# Research digest — moving from an `eps`-clock S-checkpoint to a concentration-dependent one

*Deep-research run wf_053c0fc0 (24 sources, 25 claims verified, 21 confirmed). Citations at end.*

## The headline (and it's a good one)

**Our model already descends from a published model that contains exactly the concentration-dependent
intra-S checkpoint we want — and the `eps` clock we've been fighting is not part of it.**

Two verified facts change the strategy:

1. **Gérard & Goldbeter 2009 already specifies a mechanistic, concentration-dependent ATR/Chk1 DNA-replication
   checkpoint.** It is *activated by cyclin E/Cdk2 at the initiation of DNA replication* and *inhibits the Cdc25
   phosphatases (Pa, Pb, Pe) that activate Cdk2 and Cdk1*, "slow[ing] down the dynamics of the cell cycle
   **without altering its oscillatory nature**." It extends the 39-variable core to **44 variables**. The kinetic
   equations are in the **2009 PNAS SI, section 3, scheme Fig. S5**. *This is the ATR-CHK1 → Cdc25-inactivation →
   sustained CDK-Tyr15-inhibition → delayed S/G2 route, verbatim.* [GG2009]

2. **The `eps` global clock is NOT fundamental and is removable.** The 2009 PNAS *main text uses rate constants
   in h⁻¹* (Table S2) and produces a **~19 h period natively, with no dimensionless `eps` in the main text**.
   `eps` was added downstream in the CellML/SBML encodings (`eps=17` in BioModels BIOMD0000000730; we inherited
   `eps=150`). Lowering `eps` just rescales time uniformly — it is a time-rescaling device, not real-time
   parameterization. **The native `eps`-free h⁻¹ form is the real-time alternative.** [GG2009, BIOMD730, CellML]

So the cleanest move is **not** to switch model families — it's to *complete our model back toward the full GG 2009
spec we had simplified away*: drop `eps` (use the native h⁻¹ rates) and add GG's own ATR/Chk1 module, gated on
replication state, with HU as the input.

---

## Why a concentration checkpoint *can* lengthen S even though amplitude-scaling can't

This is the crux, and it resolves the tension with our earlier finding ("period is robust to cyclin amplitude;
scaling a CDK activation flux by a constant does nothing").

- A **constant multiplier** on a flux only changes the *height* of the bistable branch → the relaxation oscillator
  is robust to it → period unchanged. (What we proved in `diag_phi_levers.py`.)
- A **checkpoint** is a *state-dependent gate*, not a constant. CHK1 → Cdc25 inhibition **raises the activation
  threshold** the slow variable must cross, and — crucially — **holds the system at that threshold (the "knee")
  until a replication signal clears.** Holding at a saddle-node/threshold makes the dwell time *diverge*
  controllably: the closer you sit to the bifurcation, the longer the wait. That is a *slow-leg / threshold*
  modification, which the period is **not** robust to. This is precisely why GG report their checkpoint slows the
  cycle while preserving oscillation. The duration of the prolongation = how long the checkpoint stays active =
  how long replication takes → **concentration-dependent, exactly what we want.** [GG2009, Lemmens2018]

The experimental rationale is the **"brake model"** (Lemmens et al. 2018, *Mol Cell*): ongoing replication, via
CHK1 and p38-MK2 phosphorylating CDK regulators, **continuously raises the CDK1/PLK1 activation threshold and
restricts mitotic entry until S completes.** Blocking origin licensing/firing caused *premature* CDK1/PLK1
activation. This is an ODE-amenable, concentration-driven checkpoint→CDK coupling. [Lemmens2018]

---

## Q1 — Models without a global `eps` (real-time-parameterized)

| Model | Structure | Timescale | Code | Verdict for us |
|---|---|---|---|---|
| **Heldt, Barr, Cooper, Bakal, Novak 2018** (PNAS) | Two bistable switches (restriction point + G1/S) around CycE:Cdk2, CycA:Cdk2, E2F, Rb, **p21, PCNA, p53** | **Real-time, minutes; NO `eps`.** 56 explicit rate constants; reproduces measured G1/S lengths | **SBML BIOMD0000000700**, Dataset S1 @ cellcycle.org.uk | **Cleanest eps-free scaffold — BUT no Cdk1/Wee1/Cdc25/Tyr15 machinery.** Its checkpoint is p53/p21 stoichiometric Cdk2 inhibition + p21–PCNA blocking replication. An ATR-CHK1-CDC25 intra-S module would have to be **added**, not reused. [Heldt2018] |
| **Csikász-Nagy et al. 2006** (Biophys J) generic eukaryotic | Generic CDK-cyclin network | **Real-time, min⁻¹; no `eps`** | **SBML BIOMD0000001044** | Real-time backbone but **yeast-generic, no intra-S replication-stress module**. Less mammalian-specific. [CsikaszNagy2006] |
| **Gérard-Goldbeter 2009** (our lineage) | 39-var core (4 Cdk modules) → 44-var with ATR/Chk1 | **Native h⁻¹ (~19 h); `eps` added only in encodings** | BIOMD0000000730 (eps=17; checkpoint presence in the XML unconfirmed) | **Best path (see below)** — has the checkpoint we want; `eps` is removable. [GG2009] |
| **Jung et al. 2021** (npj Syst Biol) p53-Plk1 | DNA-damage G2/M | Uses global scaler **τ=1.65** (like `eps`) → *not* real-time | DOI 10.1038/s41540-021-00203-8 | Not eps-free, and G2/M not intra-S. But has a **ready Cdc25-inactivation-by-ATM/ATR rate law** to graft: `Cdc25_synth ∝ KA2/(KA2+[ATM/ATR])`. [Jung2021] |

---

## Q2 — Concentration-dependent intra-S checkpoint modules

| Source | Mechanism | Ties checkpoint conc. → CDK inhibition? | S responds to replication stress? |
|---|---|---|---|
| **GG 2009 ATR/Chk1** | cyclin E/Cdk2 activates ATR/Chk1 at replication init → **Chk1 inhibits Cdc25 (Pa/Pb/Pe)** → sustained Tyr15-P on Cdk2/Cdk1 | **Yes — directly** | Conceptually yes; **HU not modeled** (would add as input) |
| **Jung 2021** | ATM/ATR ↓ Cdc25 synthesis → sustained inhibitory phospho on CDK1 in CycB:CDK1 (MPF) | Yes (G2/M) | Damage signal, not HU/dNTP |
| **Heldt 2018** | p53→p21 stoichiometric Cdk2 inhibition + **p21–PCNA blocks replication directly** | Yes, but **no Tyr15 / Cdc25 / Wee1** | DNA-damage halt; no ATR-CHK1 cascade |
| **Lemmens 2018** (experimental) | Replication → CHK1/p38-MK2 raise CDK1/PLK1 activation threshold | Yes (the "brake") | Yes — direct experimental basis |

**Gap found:** no ready-to-import SBML/ODE model surfaced that explicitly computes *"% genome replicated"*
(origin licensing/firing, fork progression, dNTP pools) to gate mitotic entry. That layer, if we want it, we'd
build ourselves (Lemmens gives the qualitative target).

---

## ⚠️ Caveat that bears directly on the core hypothesis

**Chao, Spencer et al. 2018** (*Mol Syst Biol*) found that in human cells (RPE, U2OS, H9) **G1, S, and G2 durations
are statistically independent in single cells** under basal conditions (no pair with P<0.01 *and* R²>0.1) — titled
"Evidence that the human cell cycle is a series of uncoupled, memoryless phases." **This tensions our "longer S →
more EZH2 → longer subsequent G1" inter-phase-memory premise** *under basal conditions.* The escape hatch (and it's
real): they note heritable factors create sibling correlations, and the claim "uncoupled *even under perturbation*"
was **refuted 0-3** in verification. So the coupling we hypothesize most plausibly manifests **under replication
stress (HU)**, not in unperturbed cycling — which is exactly our experimental regime. Worth designing the test to
show the coupling *appears under HU* rather than expecting it basally. [Chao2018]

---

## Recommended path: v44 = "complete the model, drop the clock"

**Stay on the GG/MYCN/EZH2 core; make two structural changes.**

### Change 1 — Remove `eps`, reparameterize in real time (h⁻¹)
`eps` is not in GG's main text. Either (a) set `eps=1` and re-fit the handful of rate constants so the period stays
~19–20 h (native GG values are published, Table S2), or (b) keep `eps` only as a units bookkeeping constant. This
removes the phenomenological clock entirely and makes every duration emerge from kinetics.

### Change 2 — Replace the phase-gated `eps` brake with a CHK1→Cdc25 checkpoint
Add a minimal concentration-dependent module (GG's own, simplified to fit our string-surgery approach):

```
# Replication progress R: integrates while in S (CycA/Me high), at fork speed v_fork(HU)
R' = v_fork(HU) * sphase_onset - k_reset * (mitosis)      # R rises during S, resets at division
v_fork(HU) = v0 * Kd^h/(Kd^h + HU^h)                      # HU slows forks (dNTP depletion)

# CHK1 active while replication INCOMPLETE (R < 1):
CHK1' = k_on * Me * (1 - R/(K_R+R)) - k_off * CHK1        # cyclin E/Cdk2 arms it; clears when R complete

# CHK1 inhibits Cdc25 (Pa, Pb) activation  →  raises CDK Tyr15-P  →  holds S/G2 transition
Pa_activation: ... * (Ki_CHK1/(Ki_CHK1 + CHK1))           # concentration-dependent gate, NOT a constant
Pb_activation: ... * (Ki_CHK1/(Ki_CHK1 + CHK1))
```

- **At HU=0:** forks fast → R reaches 1 quickly → CHK1 clears early → no delay → identical baseline.
- **At HU>0:** forks slow → R lingers below 1 → CHK1 stays high → Cdc25 held down → CDK activation threshold
  raised → **S/G2 transition waits until replication completes → S prolonged concentration-dependently.**
- **Key:** because the gate is *dynamic and state-dependent* (clears when R completes), it sits the system at the
  transition threshold and lengthens the dwell — it does NOT merely shrink a cyclin amplitude, so it escapes the
  relaxation-robustness trap.

### What to verify before committing (the open questions)
1. **Does CHK1→Cdc25 inhibition actually lengthen S in *our* parameterization?** GG say it slows their cycle; we
   must confirm it isn't washed out by relaxation robustness in our variant. Quick test: scale `Pa/Pb` activation
   *dynamically* (gated on a replication signal) and check S-dwell vs a *constant* scaling.
2. **Is the deposited BIOMD0000000730 the 39-var core or the 44-var checkpoint extension?** Unconfirmed — the
   checkpoint equations live in the 2009 SI (section 3 / Fig S5), which we should pull from the PNAS SI PDF.
3. **Calibrate `v_fork(HU)`** so the S-prolongation magnitude matches the data, and confirm the emergent EZH2 boost
   still lands ~1.31× at the 10 µM anchor.

### The fallback if GG's Cdc25 route proves too robust in our hands
Keep the **replication-completion variable R** but have it gate **mitotic entry directly** (a "licensing" gate on
the CycB/Mb transition that only opens when R≈1), per Lemmens' brake model — still concentration-dependent, still
HU-driven, and dynamically holds the transition rather than scaling amplitude.

---

## Citations
- **[GG2009]** Gérard C, Goldbeter A. "Temporal self-organization of the cyclin/Cdk network driving the mammalian
  cell cycle." *PNAS* 2009;106(51):21643-21648. DOI 10.1073/pnas.0903827106. (PMC2799800; checkpoint in SI §3,
  Fig S5; rate constants Table S2, h⁻¹; ~19 h period; no `eps` in main text.)
- **[GG2010]** Gérard C, Goldbeter A. (model reduction paper). PMC3262247. (States the detailed model "goes up to
  44 [variables] when the ATR/Chk1 checkpoint … is incorporated.")
- **[BIOMD730/CellML]** BioModels BIOMD0000000730; CellML gerard_2009 exposure (shows `eps=17` added in encoding).
- **[Heldt2018]** Heldt FS, Barr AR, Cooper S, Bakal C, Novák B. "A comprehensive model for the proliferation–
  quiescence decision in response to endogenous DNA damage in human cells." *PNAS* 2018;115(10):2532-2537.
  DOI 10.1073/pnas.1715345115. SBML **BIOMD0000000700** (minutes; no `eps`; p53/p21/PCNA, no Cdc25/Wee1/Cdk1).
- **[CsikaszNagy2006]** Csikász-Nagy A, et al. "Analysis of a generic model of eukaryotic cell-cycle regulation."
  *Biophys J* 2006;90(12):4361-4379. DOI 10.1529/biophysj.106.081240. SBML **BIOMD0000001044** (min⁻¹).
- **[Jung2021]** Jung J, et al. p53-Plk1 cell-cycle/DNA-damage model. *npj Syst Biol Appl* 2021.
  DOI 10.1038/s41540-021-00203-8 (PMC8660825). (Global scaler τ; ATM/ATR ↓ Cdc25 → CDK1 inhibition.)
- **[Lemmens2018]** Lemmens B, et al. "DNA replication determines timing of mitosis by restricting CDK1 and PLK1
  activation." *Mol Cell* 2018;71(1):117-128. DOI 10.1016/j.molcel.2018.05.026 (PMC6039720). (The "brake model.")
- **[Chao2018]** Chao HX, et al. "Evidence that the human cell cycle is a series of uncoupled, memoryless phases."
  *Mol Syst Biol* 2018;14(3):e8604. DOI 10.15252/msb.20188604 (PMC6423720). (⚠️ phase durations uncoupled basally.)
