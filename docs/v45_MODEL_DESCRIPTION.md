# v45 model — complete description

**The stochastic-commitment model: a mechanistic bistable CDK2–p27 commitment toggle + Rb-dilution
G1 timer + mitogen-gated CyclinD1 bootstrap + mass-capped quiescence, with EZH2 as a dynamic
cell-cycle species that represses CyclinD1.** This is the authoritative reference for the v45 model.

Model + harness: `simulations/v45_stochastic_commitment.py` (a standalone Tellurium/Antimony kernel
plus a Python ensemble harness — NOT a string-surgery builder like v44).
Validation: `simulations/validate_v45.py` → **20/20** against the v44 target compendium
(`docs/Ezh2_CcnD1_model_targets.md`).

> **v45 is a PARALLEL successor prototype. v44 remains the committed working model** for the paper
> figures (`src/build_model_v44_heldt.py`, `docs/v44_MODEL_DESCRIPTION.md`). v45 demonstrates that the
> same biology and perturbation responses **emerge from a mechanistic commitment kernel** with no
> phenomenological size gate.

---

## 1. Why this model exists (lineage)

- **v44** (working model): Heldt-2018 G1/S core + a **growth-gated restriction point** (`M_commit`,
  `M_size`). That gate is *phenomenological*: the period is a growth clock, the CKIs set a threshold
  not a delay, and commitment fires off **inherited** Cyclin E/A rather than a CyclinD1-driven switch.
  We also found a **count-fraction over-constraint**: a homogeneous model with period ~22 h and the
  measured short S/G2 *durations* forces the 2N count-fraction to ~80% (vs the flow value 68%).
- **v45** (this model): replaces the size gate with an **emergent bistable CDK2–p27 commitment
  toggle** (Spencer/Cappell/Yao) plus a **Rb-dilution G1 timer** (Zatulovskiy 2020) and **noisy
  birth partitioning**. The G0-vs-fast-G1 split is then a *bifurcation outcome*, not a phase every
  cell traverses. A **mixture** of fast committed cyclers + a G0-dwell tail resolves the
  over-constraint (count-fractions are cycler-dominated; the tail lifts the mean period).

Central claims v45 makes mechanistic: (i) commitment is a CyclinD1/p27-ratio bistable switch;
(ii) mitogen (CyclinD1) withdrawal causes real **arrest**, and the **MYCN floor protects MB** from
it; (iii) **EZH2 represses CyclinD1** as a within-cycle feedback, and EZH2 inhibition rescues
commitment.

---

## 2. Architecture (the kernel)

Kernel species (minutes): `mass, Rb, CDK2, p27, Cdh1, G2p, Dna, EZH2`. CyclinD1 (`Cd`) is an
algebraic assignment (below). One cell is simulated from a noisy birth state to its first division;
the population is an ensemble of such cells (Section 5).

### A. Bistable CDK2–p27 commitment toggle (the restriction point)
- `CDK2act := CDK2/(1 + p27/Ki)` — p27 buffers/inhibits CDK2.
- E2F release: `drive := d0 + wE·CDK2act`; `RbP := drive/(K_CdRb·RbC + drive)`;
  `E2F := RbP^nE/(Km^nE + RbP^nE)` (ultrasensitive). CDK2 synthesis ∝ E2F (`ksE·E2F + ksE0`).
- p27: synthesis `ksp := ksp0/(1 + Cd/Kp)` (**mitogen sets the brake** — more CyclinD1 → less p27);
  cleared by CDK2 (Skp2 analog) `(kdp·CDK2act + kdp0)·p27`. The CDK2⊣p27⊣CDK2 loop is **bistable**
  for birth `RbC≈5–11` (verified by `fixed_points()`, the nullcline analyzer).
- **Cd enters commitment via the p27-clearance dynamics, NOT the separatrix** (the separatrix is on
  the Cd-independent E2F nullcline). Low CyclinD1 → high p27 synthesis → slower escape from the low
  basin → longer/permanent G0.

### B. Rb-dilution G1 timer (Zatulovskiy 2020)
- Rb is synthesized **size-independently** (`RbSyn: => Rb; ksRb`, constant), so `RbC := Rb/mass`
  **dilutes as the cell grows**. Birth `RbC≈9` (= `ksRb/kdRb`) sits inside the bistable window.
- Origins fire only once Rb concentration dilutes below `KfireRb` (=3.5): the `fire` term gates
  S-entry on **committed (CDK2act high) AND Rb-diluted (RbC low)**. G1 length = Rb-dilution time
  (~17 h), which sets the committed-cycler period without a CyclinD1-dependent collapse.

### C. Mitogen-gated CyclinD1 bootstrap → real arrest
- `d0 := dmax·Cd^nd0/(Kd0^nd0 + Cd^nd0)` (saturating Hill-2). This is the **CyclinD–CDK4/6
  mono-phosphorylation** of Rb: at Cd→0 there is **no basal Rb-P**, so dilution alone cannot start
  E2F → the cell cannot commit. Cd enters ONLY this basal term, never the `wE·CDK2act` toggle gain
  (gain-coupling collapses the MB period). Saturating + Hill-2 keeps MB (Cd≈7) from over-committing
  and sharpens the low-Cd cutoff (Cd≈0.25 arrests; Cd≈0.78 still cycles).

### D. Mass-capped quiescence (permanent arrest)
- `Growth: => mass; mu·mass·(1 − (mass/Mmax)^nM)` (Mmax=4.2). Cycling cells divide at mass ~3.5
  (below the cap); a cell stuck in the low basin **grows to the cap and stops**, flooring `RbC` so it
  stays bistable-low → **permanent arrest**. Without the cap, unbounded growth dilutes RbC→0 and
  forces commitment at *any* Cd>0 (low-Cd cells would only slow, never arrest). Together C+D give a
  sharp **commit-or-arrest bifurcation** at Cd≈0.2–0.3.

### E. Clean S phase (latch) + real G2 (timer)
- **S**: `DnaRep: => Dna; kFire·(fire + latch − fire·latch)·(1−Dna)`. `fire` (gated on commitment +
  Rb dilution = the G1 timer) only *initiates*; once replication starts, `latch := Dna^8/(0.15^8+…)`
  drives it to completion at full rate. This **decouples S duration from the slow RbC-firing ramp** →
  a clean ~2.5–3 h S. Threshold 0.15 sits above the pre-commit fire-leak so the latch never bypasses
  the G1 timer.
- **G2**: `G2acc: => G2p; kG2·(1−Cdh1)·g2gate` — a clock that runs only when committed (Cdh1 off) AND
  replicated (`g2gate := Dna^nG2/…`); division fires at `G2p > 1` → G2 duration ~2.5 h. `Cdh1`
  (APC/C) is inactivated by CDK2act (`koff·CDK2act·Cdh1`) = the point of no return.

### F. EZH2 → CyclinD1 feedback (the paper's mechanism)
- **EZH2 is a dynamic cell-cycle species**: synthesis `kEZbas + kEZsyn·CDK2act` **tracks CDK2act →
  peaks in S/G2** (Fig 4A/B; model G2/G0 ratio ~1.7×). It is **STABLE** (small `kDeEZ`) so it
  **integrates** over the cycle, and is **halved at division** (dilution).
- **EZH2 represses CyclinD1**: `Cd := Cd_drive·K_EZH2_Cd/(K_EZH2_Cd + EZH2·(1−EZH2i))`. `Cd_drive`
  is the EZH2-*independent* Gli/MYCN drive (Section H); `EZH2i=1` removes the brake → `Cd → Cd_drive`
  (de-repression). So the feedback is **internal to the v45 cycle** and the EZH2i rescue is
  **emergent** (GNP Cd 0.94→2.56, ≈ Fig 3C qPCR 2.2×).
- *Protein-vs-transcript note:* the repressing EZH2 is only ~1.34× in MB, not the 2.05× RNA-seq
  transcript fold (Fig 4J). Forcing 2× (a steeper `CDK2act^2` synthesis) over-represses MB CyclinD1
  (Cd 3.06 vs the required ~7×) — the same reason v44's repressing EZH2 *protein* is ~condition-
  independent (1.03×) while the transcript is 2.05×. v45 validates the cycle-dependence, MB-direction,
  and de-repression, not the bulk transcript magnitude.

### G. Birth noise & division (in the Python harness, not the kernel)
- `run_cell` resets a daughter at birth: `mass=1`, `Rb=ksRb/kdRb·mass`, `Cdh1=1`, `G2p=Dna=0`, and
  **inherited** `p27_0`, `CDK2_0`, `EZH2_0`. Division = `G2p>1`.
- The ensemble draws each from a lognormal: `p27_0 = P21_div·LN(0,sig_p)`, `CDK2_0 = phi·LN(0,sig_c)`,
  `EZH2_0 = E_eq·LN(0,sig_ez)`. **Noisy p27/CDK2 partitioning across the separatrix** gives the
  born-committed fraction + a G0 tail — a right-skewed period distribution **from birth noise alone**
  (no SDE needed). `equilibrate_ezh2()` finds the self-consistent birth EZH2 (the halving fixed point).

### H. CyclinD1 drive from the v44 cascade (wiring)
- `condition_drive(name)` runs the v44 HH/MYCN→CyclinD1 cascade (`build_model_v44(with_hh=True)`) to
  steady state with **EZH2i=1** (repression removed) and GNP-normalizes → `Cd_drive` (GNP 2.56,
  MB 19.0, etc.). The cascade is quasi-static vs the ~hours cycle, so a steady drive is principled.
  v45's own EZH2 then re-applies the repression. `condition_cd()` (repressed cascade Cd) is kept as a
  reference.

---

## 3. Inputs / conditions

Conditions map to the v44 cascade boundary inputs (mirrors `validate_v44.py`); `CONDITIONS` dict in
the module. Each yields a `(Cd_drive, EZH2i)` via `condition_drive`.

| Input | Meaning | GNP | MB |
|---|---|---|---|
| `SHH` | Hedgehog ligand | 0.5 | 0.5 |
| `Ptch1_copy_number` | Ptch1 dosage (1=WT, 0.1=loss/MB) | 1.0 | 0.1 |
| `MYCN_amplification` | MYCN level | 1.0 | 2.8 |
| `GDC0449` | SMO inhibitor / HHi-vismo (0–1) | 0 | 0 |
| `EZH2i` | EZH2-inhibitor toggle (1 = repression OFF) | 0 | 0 |

Normalized `Cd_drive` (cascade, EZH2i=1, GNP-normalized): GNP 2.56, GNP−SHH 0.52, GNP+HHi 0.30,
MB 19.0, MB+HHi 2.44. Cycle-mean **effective** CyclinD1 (after v45 EZH2 repression): GNP 0.93, MB 5.66.

Kernel parameters of note (defaults are the validated set): `mu=0.0008`, `Mmax=4.2 nM=16`,
`ksRb=0.0135 kdRb=0.0015` (Rb_ss=9), `KfireRb=3.5`, `wE=2.6 dmax=0.10 Kd0=0.82 nd0=2 K_CdRb=0.45`,
`ksE=0.040 ksE0=0.0015 kdE=0.013 Ki=0.12`, `ksp0=0.040 Kp=0.4 kdp=2.5`, `kG2=0.0067`, `kFire=0.022`,
`kEZbas=0.0002 kEZsyn=0.0009 kDeEZ=0.0001 K_EZH2_Cd=1.0`. Birth noise: `P21_div=0.42 sig_p=0.55
phi=0.55 sig_c=0.5 sig_ez=0.25`.

---

## 4. Phase markers (how to classify G0/G1/S/G2M)

`phase_at()` per time point, on DNA content + commitment:
- **S** = `0.05 ≤ Dna < 0.95` (replicating).
- **G2M** = `Dna ≥ 0.95` (until division at `G2p>1`).
- **2N** = not S, not G2M; split into **G0** (`CDK2act < 0.30`) vs **G1** (`CDK2act ≥ 0.30`).
  Flow cannot resolve G0/G1 (both 2N) — fit the **2N pool** and treat the split as model-internal.

Population read-outs (`count_fractions`, age-weighted by `n(a)=2λe^{−λa}`, λ from Euler-Lotka):
**flow COUNT-fractions**, not duration-fractions. `arrest_frac` = newborns that permanently exit
(mass-capped, low basin) = the **pRb⁻/Ki67⁻ analog**. `ezh2_stats` = cycle-mean EZH2, G2/G0 ratio,
effective Cd.

---

## 5. Status — validated (`validate_v45.py` → 20/20)

**CyclinD1 folds** (cycle-mean Cd_eff = Gli/MYCN drive repressed by v45 EZH2): MB/GNP 6.1 (t 7.58);
MB+HHi/MB, GNP+HHi/GNP (vismo crash); **GNP+EZH2i/GNP 2.74** (EZH2i de-repression, Fig 3C ~2.2×).
**EZH2 dynamics**: G2/G0 1.68 (Fig 4A/B ~2×); elevated in MB vs GNP; EZH2i de-represses CyclinD1.
**MB DMSO count-fractions**: 2N 70 (t 68.2), S 19 (t 15.7), G2+M duration 2.5 h (direct anchor),
quiescent 25% (data 21.7%). **Period**: GNP median 22.8 h (lit ~22–23 h). **Viability**: MB/GNP
cycle; GNP−SHH strongly quiescent; GNP+HHi permanent arrest.

**Emergent biology (the payoff):**
- **MYCN floor protects MB from Hh-withdrawal arrest**: GNP+HHi (Cd 0.16, below threshold) → 74%
  permanent arrest, but MB+HHi (Cd 0.78, MYCN-floored *above* threshold) → 0% permanent arrest, only
  transient G0 — the paper's claim, emergent not imposed.
- **EZH2i rescue**: MB+HHi quiescent 45→32%.
- **Mixture resolves the over-constraint**: fast committed cyclers (clean ~15 h mode) + a G0 tail →
  2N count ~68–70% with a right-skewed period (median ~22 h), from birth noise alone.

**Known limitation / out of scope:** the bulk EZH2 transcript fold (2.05×, Fig 4J) is not reproduced
by the repressing pool (see §2F). The within-cycle EZH2 feedback does not lengthen the period (matches
v44; the steady EZH2 is ~condition-independent).

---

## 6. Build & run

```bash
./venv/bin/python simulations/v45_stochastic_commitment.py    # single-cell + ensemble + rescue panel
./venv/bin/python simulations/validate_v45.py                 # full 20/20 validation suite
```
```python
import sys; sys.path.insert(0, "simulations")
import v45_stochastic_commitment as v45
drive, ezh2i = v45.condition_drive("MB+HHi")          # Cd_drive + EZH2i from the v44 cascade
cells, st = v45.ensemble(drive, EZH2i=ezh2i, N=300, return_stats=True)
cf, lam = v45.count_fractions(cells)                  # flow count-fractions
print(st["arrest_frac"], cf, v45.ezh2_stats(cells))
v45.fixed_points(RbC=9, Cd=7.0)                       # nullcline analyzer (bistability check)
```

Harness API: `run_cell(p27_0, CDK2_0, Cd_drive, EZH2_0, EZH2i)` → `(traj, CDK2_div, EZH2_div, status)`
with status `divided|arrested|failed`; `ensemble(...)`; `equilibrate_ezh2(...)`; `count_fractions`;
`population_growth_rate` (Euler-Lotka); `fixed_points` (nullcline analyzer); `condition_drive`,
`condition_cd`, `rescue_panel`, `ezh2_stats`.

---

## 7. Provenance

- Spencer, Cappell, et al. (2013) *Cell* 155:369; Cappell et al. (2016/2018); Yao et al. (2008) —
  CDK2/p27 bistable commitment and the proliferation–quiescence bifurcation.
- Zatulovskiy et al. (2020) *Science* 369:466 — Rb dilution by growth as the G1 size/timer.
- Fan & Meyer; Guiley et al. (2019) — CyclinD1/p27 ratio commitment, p27→CDK4/6.
- v44 working model: `src/build_model_v44_heldt.py`, `docs/v44_MODEL_DESCRIPTION.md`; targets:
  `docs/Ezh2_CcnD1_model_targets.md`.
