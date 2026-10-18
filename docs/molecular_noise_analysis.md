# Units, molecule numbers, and how to inject molecular noise over time

Written 2026-09-29 in response to JP's question, to read on waking. Verified against the current
`transient-g0-reframe` model (`build_model_v44()`), not from memory.

---

## 1. What units are the equations in? — Concentration, dimensionless, one well-mixed compartment

- The model has **one compartment, `Cell`, with size = 1.0** (verified). Every species is declared
  `species X in Cell`.
- Every reaction is written `Cell * (mass-action rate in concentrations)`, e.g.
  `Synthesis_of_p21: => P21; Cell*(kSyP21 + kSyP21P53*P53)` and
  `Association…: Rb + E2f -> RbE2f; Cell*(kAsRbE2f*Rb*E2f - kDsRbE2f*RbE2f)`.
- Because `Cell = 1`, the **amount and the concentration are numerically identical**, and there is no
  physical volume, no Avogadro factor, and **no molecule count anywhere**. The species values are
  **dimensionless / arbitrary concentration units** (this is inherited from Heldt 2018 and the
  Novak–Tyson lineage it comes from). `Cd = 8.1`, `Cd_mRNA = 10.1`, `EZH2 = 5.4`, `P21 (p27) = 0.003`
  are all "a.u.", not molecules.

**Consequence that matters for your question:** the model is a *deterministic, continuous,
mean-field* ODE. It represents the *average* trajectory of an infinitely large, perfectly-stirred
pool. It therefore has **exactly zero intrinsic noise** — run it twice from the same state and you
get the identical curve. All the cell-to-cell and over-time variation we've been adding is *bolted
on from outside*. That is exactly why the recent drift looked arbitrary: nothing in the equations
tells us how big it should be. To fix that we have to re-introduce the molecule scale we threw away.

---

## 2. Molecule numbers — for our own edification

To talk about molecules we have to pick a **system size Ω** = molecules per unit of a.u.
concentration (Ω folds in the cell volume and the a.u.→molar scale). Once Ω is chosen, molecules =
Ω × (a.u. value). The model itself cannot tell us Ω; we anchor it to measured absolute abundances.

**A good anchor (Schwanhäusser 2011, NIH3T3, genome-wide):** median **mRNA ≈ 17 copies/cell**,
median **protein ≈ 50,000 copies/cell**, i.e. a protein:mRNA ratio ~ few thousand, protein t½ median
~46 h, mRNA t½ median ~9 h. Cell-cycle regulators are lower-copy and more bursty than the median.

If we anchor the model so that a moderately-expressed protein sits near 50,000 molecules, then in
round numbers:

| species (model a.u.) | what it is | plausible copies/cell | noise floor CV≈1/√N |
|---|---|---|---|
| CyclinD1 protein (`Cd` ≈ 8) | protein | ~10^4–10^5 | 0.003–0.01 (tiny on its own) |
| CyclinD1 mRNA (`Cd_mRNA` ≈ 10) | mRNA | **~10–50** | **0.14–0.32** |
| EZH2 protein (`EZH2` ≈ 5.4) | protein | ~10^4–10^5 | ~0.005 |
| EZH2 mRNA (`EZH2m` ≈ 0.26) | mRNA | **~5–20** | **0.22–0.45** |
| p27 (`P21` ≈ 0.003) | protein (held low) | low here by design | — |

**This table is the whole point.** A protein taken *by itself* is high-copy and essentially
noiseless (CV ~0.01). **The mRNA is low-copy (tens), and that is where the variance is born.** Your
intuition — 10 mRNA → ~400 protein, 11 mRNA → ~440, "but probably there is variance there" — is
exactly the physics: a ±1 mRNA fluctuation is a ~10% swing, and because each mRNA is translated into
a *burst* of many proteins, that swing is stamped onto the protein pool. Whether any individual
molecule is made or destroyed "depends on being next to the complex" — that is Poisson shot noise on
the reaction events, and it is loudest where the counts are smallest (the mRNA and the enzyme
encounters), not on the abundant protein.

**The governing result (Thattai & van Oudenaarden 2001; Paulsson 2004).** For a gene → mRNA →
protein cascade with burst size `b = k_translation / γ_mRNA` (proteins made per mRNA over its
lifetime) and `η = γ_protein/γ_mRNA = τ_mRNA/τ_protein`:

```
⟨protein⟩ = b · ⟨mRNA⟩
CV²_protein = 1/⟨protein⟩ + (1/⟨mRNA⟩) · 1/(1+η)
            ≈ 1/⟨mRNA⟩         (for a stable protein, τ_protein ≫ τ_mRNA, and b ≫ 1)
```

**Protein noise is set by the mRNA copy number, not the protein copy number.** A measured protein
CV of 0.70 (our CyclinD1 number) would, *if it were purely intrinsic bursting*, imply ⟨mRNA⟩ ≈
1/0.70² ≈ **2** — implausibly low. So the measured 0.70 is **not all intrinsic**: it is intrinsic
(a few–tens of mRNA, fast) **plus** extrinsic (slow, shared) — see §4.

**Model-specific caution.** On this branch the CyclinD1 mRNA was retuned into a slow reservoir
(`k_Cd_mRNA_deg = 0.001`, t½ ≈ 11.5 h) with `k_Cd_translation = 0.801`, giving a burst size
`b = k_tl/γ_m ≈ 800`. EZH2's is `b ≈ 1.5` (normal). A burst size of 800 is molecularly extreme and
means the CyclinD1-mRNA "concentration" (a.u. 10) is decoupled from any real copy number. When we
add molecular noise we will have to decide the *real* CyclinD1 mRNA lifetime — which is the same
open question flagged as F2 in `protein_noise_fixes.md`. **The noise work and the mRNA-lifetime
question are the same problem.**

---

## 3. Two kinds of noise, and why the distinction is everything for G0

Elowitz et al. 2002 split single-cell variability into two components; both are real and they have
opposite time signatures:

- **Intrinsic noise** — the shot noise of discrete birth/death reactions on *this* gene. Magnitude
  ~1/√(copies), so it lives in the low-copy mRNA and the enzyme-encounter events. **Correlation time
  ≈ the species' own lifetime** (mRNA lifetime, low-pass-filtered by the protein lifetime). It is
  *mean-reverting*: a downward blip is pulled back toward the deterministic steady state on the
  timescale τ_protein. **It does not wander.**
- **Extrinsic noise** — slow, cell-wide factors shared across many genes (cell size, ribosome/
  polymerase capacity, upstream signaling tone, metabolic state). Correlated across proteins,
  **correlation time ≈ one cell cycle or longer**. This is the component that makes a cell "have more
  protein and keep it for hours," and it is what sibling/cousin correlations (Sigal–Alon memory,
  Nakashima sibling heritability) actually measure.

Your objection — "I can't imagine all the cell-to-cell variance is a cell just having more protein
and staying that way" — is precisely correct, and it names the failure of the two things we tried:

- the **static per-cell draw** (`sim_population_g0`, the bimodal tail) = *all extrinsic, frozen
  forever*. No dynamics, so a cell can never fall into G0 on its own.
- the **AR(1) redraw-every-division** (the recent figures) = one oversized walk with **no permanent
  identity and a made-up amplitude** (CV ~1.0–1.2, 2–3× the data). It moved too much and for the
  wrong reason.

The right model is **both components at once**: a small–moderate *intrinsic* part that fluctuates on
the hours timescale (this is what lets a cell dip below the mitogen threshold and fall into transient
G0, then recover) **plus** a slow *extrinsic* part that sets the persistent cell-to-cell spread and
the sibling correlations. Neither alone is right.

---

## 4. How to inject it — from the molecules, with the magnitude and the timescale *derived*, not tuned

The key win: once we commit to a molecule scale Ω, **both the size and the speed of the fluctuations
come out of the chemistry that is already in the model.** We stop hand-setting σ and ρ.

### Option A (rigorous) — Chemical Langevin Equation (CLE)
For each reaction j with propensity aⱼ(x) and stoichiometry νⱼ, replace the deterministic ODE
`dx = Σⱼ νⱼ aⱼ dt` with:

```
dx = Σⱼ νⱼ aⱼ dt  +  Σⱼ νⱼ √(aⱼ/Ω) dWⱼ
```

- The `√(aⱼ/Ω) dW` term is the discreteness/shot noise of individual reaction firings — literally
  "whether this molecule got degraded depends on the encounter." It is **automatically** large where
  counts/fluxes are small (mRNA) and negligible where they're large (protein). 1/√N falls out.
- **Ω is the single physical knob.** Bigger Ω = more molecules = quieter. We set Ω *once* by matching
  the measured single-cell CV of a reference (CyclinD1, IF), and every other species' noise is then a
  prediction, cross-checkable against its own measured CV.
- Correlation times are not chosen — they are the relaxation times of the network (mRNA lifetime,
  protein lifetime, the R-point relaxation), so the fluctuations mean-revert on their own timescales.
- Cost: needs a fixed-step SDE integrator (Euler–Maruyama) instead of CVODE; the stiff mitotic
  switch wants a small step. Feasible but heavier; best applied to a subset of nodes (below).

### Option B (recommended, interpretable) — explicit bursty transcription on the low-copy nodes
Keep the deterministic core; make **only the mRNAs that gate the decision** (CyclinD1, p27/`P21`
synthesis, EZH2) genuinely stochastic and let deterministic translation amplify them. Two equivalent
implementations:
- *Telegraph / bursting:* the gene switches on/off; while on it fires mRNA in bursts of mean size
  `b`; mRNA decays at γ_m. This reproduces the measured super-Poissonian protein noise with two
  interpretable numbers (burst frequency, burst size) that are literally measurable by smFISH.
- *mRNA-only CLE:* run just the mRNA species with the √(a/Ω) term, protein deterministic. Cheaper.

This targets the dominant physiological noise source, matches your 10→11-mRNA picture directly, and
the protein fluctuation timescale = the protein lifetime (hours) with no extra parameter.

### Option C (reference) — full Gillespie / τ-leaping
Scale all species to molecule counts (×Ω) and run exact stochastic simulation. Correct by
construction but the stiff switch + 47 species makes it slow and fragile; use it only to *validate*
Option A/B on a reduced module, not for production runs.

### The extrinsic term (add alongside A or B, small and slow)
A single slowly-varying multiplier on synthesis capacity, shared across the noisy genes, as an
Ornstein–Uhlenbeck process with **correlation time ≈ 1 cell cycle** and an amplitude set so that
`CV_total² = CV_intrinsic² + CV_extrinsic²` matches the measured total (e.g. CyclinD1 0.70). This is
the *right* version of what the AR(1) was groping at: slow, shared, and sized to the data — not
σ=0.9 re-rolled every division.

**Suggested split for CyclinD1** to hit total CV ≈ 0.70: intrinsic CV ≈ 0.25 (fast, hours, from
~10–15 effective mRNA/bursts) + extrinsic CV ≈ 0.65 (slow, ~1 cycle, shared). The intrinsic part is
what produces transient G0; the extrinsic part sets which cells sit near the edge.

---

## 5. Why this gives G0 the way you want

The EZH2/H3K27me3 governor sets a **mitogen threshold**: a cell commits only if CyclinD1 clears the
mark-set bar. In the deterministic model a cell is forever on one side of that bar. Add the intrinsic
noise and: a downward CyclinD1 fluctuation — a stretch with fewer mRNA/bursts, lasting ~τ_mRNA and
smoothed over ~τ_protein — transiently drops CyclinD1 **below the threshold**, the cell misses the
restriction point that cycle and sits in G0, then climbs back out when transcription resumes. The
dwell length is **set by the CyclinD1 lifetime**, not by a hand-picked number, and it is **reversible
and mean-reverting by construction** — the cell is pulled back to its deterministic setpoint, it does
not random-walk away. Cells nearer the threshold (via the slow extrinsic component) dip below it more
often — so the G0 fraction emerges from measured variances rather than from a tuned tail.

This is the mechanism you asked for: variation over time, born in the molecules, that lets cells fall
into and out of G0 — without proteins wandering randomly and without freezing heterogeneity for days.

---

## 6. What we'd need to parameterize it (and its provenance)

1. **Ω (molecule scale)** — set by matching one measured single-cell CV (CyclinD1 IF; already have
   ~0.70). One number, physically interpretable, cross-checkable against absolute proteomics.
2. **Real CyclinD1 mRNA lifetime** — currently 11.5 h on this branch (b≈800, implausible); needs a
   measured value (smFISH / ActD chase in GNP). This is the F2 question — noise + lifetime are one job.
3. **Intrinsic vs extrinsic split per gene** — from paired single-cell measurements or the
   sibling/cousin correlation structure (Sigal–Alon; Nakashima siblings CV 0.12 ≪ non-sib 0.32 gives
   the extrinsic fraction directly).
4. **Burst size/frequency** for CyclinD1, p27, EZH2 — smFISH would pin these; until then, inferred
   from the target CVs.

The one experiment that would anchor most of this: **smFISH (or scRNA) for CyclinD1 in GNP/MB** — it
gives ⟨mRNA⟩, the burst statistics, and the intrinsic CV in one shot.

---

## 7. Concrete first step (for when we pick this up)

1. Choose Ω by matching CyclinD1 IF CV, on `main` (where the CyclinD1 mRNA is not the 800-burst
   reservoir) or after resolving the mRNA lifetime.
2. Implement **Option B** (bursty CyclinD1 mRNA) as a hybrid over short chunks (chunk ≈ τ_mRNA), the
   rest deterministic — a small, testable change, and it reuses the chunked-lineage harness we
   already have.
3. Add the slow extrinsic OU multiplier (τ ≈ 1 cycle) sized so total CV matches data.
4. Read out the G0 fraction and the dwell distribution and check both against measured values — no
   hand-tuning of amplitude or timescale; only Ω and the intrinsic/extrinsic split, both anchored.

Nothing here has been coded or changed — this is analysis only.
