# v44 model — complete description

**The MYCN–EZH2–CyclinD1–Hedgehog cell-cycle model, rebuilt on the Heldt 2018 core with explicit
DNA replication, a real mitotic switch, a growth-gated restriction point, and a Skp2–p27
feedforward.** This is the authoritative reference for the v44 model.

Builder: `src/build_model_v44_heldt.py` → `build_model_v44(hu, with_ezh2, with_hh, with_growth,
with_skp2, params)` returns an Antimony string (load with `tellurium.loada`).
Base model (frozen): `models_external/heldt2018.ant` (Heldt et al. 2018 PNAS, BioModels
BIOMD0000000700), converted to Antimony; the builder injects everything below by string surgery.

---

## 1. Why this model exists (lineage)

- **v42** (published): MYCN–EZH2–CyclinD1–HH on the **Gérard–Goldbeter** relaxation-oscillator
  cell cycle. Frozen baseline (`src/build_model_v42_mycn.py`).
- **v43**: intra-S checkpoint via a **phase-gated global clock** (`eps` slowed during S). Worked,
  but `eps` is phenomenological, and we proved the GG oscillator's **S-phase duration is a
  structural invariant** — no molecular concentration can change it.
- **v44** (this model): rebuilt on **Heldt 2018** — real-time (minutes, **no `eps`**) with
  **explicit DNA replication**, so S-phase duration is mechanistic and concentration-dependent.
  Then extended with mitosis, growth control, and the Skp2–p27 switch.

Central hypothesis being tested: **longer S-phase → more EZH2 → stronger CyclinD1 repression →
longer G0/G1** (in SHH-driven cerebellar granule-neuron progenitors, GNPs, and SHH medulloblastoma,
MB).

---

## 2. Architecture (modules)

The model is the Heldt G1/S core + 6 injected modules. Heldt is a one-shot proliferation–quiescence
model (CyclinE/A, Rb–E2F, p21, PCNA, explicit replication) with **no mitosis**; v44 makes it cycle
and adds the regulatory layers.

### A. Heldt G1/S core (inherited)
- **Restriction point**: `Rb` (hypophosphorylated) ⇄ `pRb` (phosphorylated) ⇄ `RbE2f` (Rb·E2F).
  Rb is phosphorylated by `(kPhRbCd·Cd + kPhRbCe·Ce + kPhRbCa·Ca)`; phospho-Rb releases `E2f`.
- **E2F positive feedback**: `E2f` synthesis = `kSyE2f + kSyE2fE2f·E2f/(jSyE2f+E2f)`.
- **CyclinE / CyclinA (CDK2)**: `Ce`, `Ca` synthesized ∝ `E2f`; bound/inhibited by p21 (`CeP21`,
  `CaP21`); degraded (CycA by Cdh1 `C1`).
- **p21/p27 (CKI)**: `P21` — inhibits CDK2; degraded by `(kDeP21 + kDeP21Cy·Skp2·(Ce+Ca) +
  kDeP21aRc·Cdt2·aRc)`. (Heldt calls it p21; we use it as the **p27** marker.)
- **APC/C-Cdh1**: `C1` (active) ⇄ `pC1` (phospho, inactive); inactivated by CyclinE/A
  (`kPhC1Ce·Ce + kPhC1Ca·Ca`); `E1`=Emi1. Active in G0/G1, off in S/G2.
- **DNA replication**: origins `Rc` → fired `pRc` (by CDK2: `kPhRc·(Ce+Ca)^n/(jCy^n+…)`) → loaded
  `aRc` (with PCNA `aPcna`) → **`Dna` synthesized at `kSyDna·aRc`** until `Dna`=1. p21 blocks via
  PCNA (`iPcna`, `iRc`). DNA damage `Dam`→`P53`→p21.

### B. Explicit DNA replication + HU coupling
- `Synthesis_of_DNA: aRc => aRc + Dna; kSyDna·vfork·aRc` — **`kSyDna` = fork speed**.
- **HU (hydroxyurea)**: `vfork := vmin_fork + (1−vmin_fork)·KmHU^h/(KmHU^h+HU^h)`. HU depletes
  dNTPs → slower forks → **S lengthens concentration-dependently** (G1 unaffected; mitosis waits).

### C. CyclinB/CDK1 (MPF) mitotic switch + intra-S checkpoint + real G2 (**redesign #1**)
- `MPF` (active CyclinB-CDK1), `preMPF` (Tyr15-P inactive), `Cb := MPF+preMPF`.
- **Cdc25/Wee1 hysteresis**: `Cdc25a` (MPF-activated, Chk1-inhibited), `Wee1a` (MPF-inhibited,
  Chk1-activated), basal `a25/aWee` ignite the switch.
- **CHK1 intra-S checkpoint**: `Chk1 := aRc/(jChk+aRc)` — active replication forks hold MPF off
  until `Dna`→1 (forks disassemble, `aRc`→0).
- **Real G2**: CyclinB synthesized (as `preMPF`) only as replication completes —
  `SynCycB: => preMPF; kSyCb·g2gate`, `g2gate := Dna^nG2/(KG2^nG2+Dna^nG2)`. So CyclinB-CDK1
  **builds during G2** (a genuine multi-hour gap) instead of being pre-stocked and flipping instantly.
- **APC/Cdc20**: `Cdc20` activated by MPF, degrades CyclinA/B → mitotic exit.

### D. Division reset (event)
`E_div: at (MPF > MPF_div)` resets to a daughter: `Dna=0`, re-license origins (`Rc=1`), `Rb=Rb+pRb,
pRb=0` (hypophosphorylated), **`P21=P21_div` (high p27)**, `CeP21=CaP21=0`, **`Skp2=0.05` (low)**,
`Ce=Ce_div, Ca=Ca_div` (low — mitotic APC degradation), `E1/=2`, `MPF=preMPF=Cdc20=0`,
`mass/=2`, `EZH2/=2, EZH2m/=2`.

### E. Growth-gated restriction point (**redesign #2**)
- `mass` grows exponentially (`Growth: => mass; mu·mass`), halves at division (size homeostasis).
- **Size gate on S-entry**: origin firing × `size_gate := mass^n/(M_size^n+mass^n)`. Committed cells
  **wait for size** before replicating → **G1 length = growth time, DECOUPLED from CyclinD level**
  (high-CyclinD MB no longer has a collapsed G1). Sets period ~22 h.

### F. Skp2–p27 feedforward restriction-point switch (**redesign #3**)
- **Skp2 is now a dynamic E2F target degraded by APC/C-Cdh1 (`C1`)**:
  `Skp2_synthesis: => Skp2; kSySkp2bas + kSySkp2·E2f`; `Skp2_degradation_Cdh1: Skp2 => ; kDeSkp2C1·C1·Skp2`.
- **G0**: E2F off → low Skp2 synthesis + active Cdh1 degrades Skp2 → Skp2 low → **p27 not degraded
  → p27 high (constitutive)** → CDK2 inhibited → Rb hypo-P → E2F off. Self-reinforcing, p27⁺.
- **Commit**: CyclinD partial Rb-P → some E2F → Skp2 rises → p27 degraded → CyclinE/CDK2 active →
  Cdh1 phosphorylated OFF → Skp2 stabilized → full Rb-P → more E2F. **This IS the bistable R-point.**
- `kPhRbCd` raised 0.2→0.5 so the threshold sits between GNP−SHH (Cd~0.28, quiescent) and GNP+SHH
  (Cd~0.5, cycling).

### G. EZH2 epigenetic layer
- `EZH2m` (transcript), `EZH2` (protein). Transcription gated on **`E2f·(wCe·Ce + (1−wCe)·Ca)`**
  (Heldt analog of v43's E2F×(CycE+CycA) S-window gate).
- **EZH2 protein is stable (slow `kDeEZ`) and diluted 2× at division → it INTEGRATES S-phase
  duration**: a longer S → more EZH2. This is the mechanistic core of the hypothesis.
- **EZH2 represses CyclinD1**: `Cd` transcription × `K_EZH2/(K_EZH2 + EZH2·(1−EZH2i))`. `EZH2i`
  toggles the feedback (1 = OFF / EZH2 inhibitor).

### H. Hedgehog + MYCN → CyclinD1
- SHH → Ptch1 → Smo → Gli (rep/act) → Gli1 → **CyclinD1 (`Cd`)** transcription. **MYCN** drives
  CyclinD1 GDC-independently (Hill). EZH2 represses. `Cd_mRNA → Cd` (protein).
- Quasi-static (equilibrates fast vs the cell cycle). Inputs below set the cell-type/drug context.

---

## 3. Inputs / controls

| Input | Meaning | GNP | MB |
|---|---|---|---|
| `SHH` | Hedgehog ligand | 0.5 | 0.5 |
| `Ptch1_copy_number` | Ptch1 gene dosage (1=WT, ~0.1=loss/MB) | 1.0 | 0.1 |
| `MYCN_amplification` | MYCN level (1=WT, 2.8=MB) | 1.0 | 2.8 |
| `GDC0449` | SMO inhibitor / HHi (0–1) | 0 | 0 |
| `HU` | hydroxyurea (0=none, 1≈10 µM) | 0 | 0 |
| `EZH2i` | EZH2-inhibitor toggle (1 = feedback OFF) | 0 | 0 |

Set at runtime: `rr['SHH']=…`. Const kinetic params changed via `build_model_v44(params={...})`.

---

## 4. Phase markers (how to classify G0/G1/S/G2)

Per the experimental assays:
- **G0 = phospho-Rb(Ser807/811)-NEGATIVE OR p27-POSITIVE** (union). In the model: `pRb` low (`<~2`)
  **OR** `P21` high (`>~0.1`). A *transient* G0 occurs in all cells (newborn high p27, low Skp2).
- **G1** = phospho-Rb⁺ AND p27⁻, pre-S (no origins fired).
- **S** = replicating: `aRc` high (origins fired) and `Dna < ~0.98`.
- **G2/M** = `Dna ≥ ~0.98` until division.

Key readouts: `pRb`, `P21`(=p27), `Skp2`, `Ce`,`Ca`, `E2f`, `aRc`, `Dna`, `MPF`/`Cb`, `mass`,
`Cd`(CyclinD1), `EZH2`/`EZH2m`, `Gli1`, `MYCN`.

---

## 5. Status (what works / what's pending)

**Working & validated (committed on `main`):**
- Concentration-dependent S-phase (HU slows forks → S lengthens, G1 flat, mitosis waits, no arrest).
- EZH2 integrates S-duration; within-cycle EZH2 **transcript** gradient G0<G1<S<G2, S/G0≈2.0
  (Section D target 1.8–2.5).
- Real G2 phase (redesign #1).
- Growth-gated, CyclinD-decoupled G1, ~22–23 h period (redesign #2).
- **GNPs proliferate only under Hedgehog**; GNP−SHH and GNP+HHi quiesce; MB cycles SHH-independently;
  MB partially GDC-resistant (MYCN buffers).
- Skp2–p27 feedforward bistable R-point (redesign #3); proliferation–quiescence correct.
- HH/MYCN ratios: MYCN MB/GNP, Gli1+HHi >99%, CyclinD1 GNP+HHi all matched; CyclinD1 MB/GNP fittable
  to ~4× via the optimizer (`simulations/v44_optimize_hh.py`; params recorded, not yet baked in).

**Pending re-calibration (on the new ~23 h Skp2-p27 structure):**
- **Transient G0** — DONE: commitment is coupled to cell size (`commit_gate` at `M_commit`), so the
  cell grows in a p27-high G0 until the Skp2-p27 feedforward fires. MB G0 ~21% (target ~22%), GNP
  ~37%, both cycling, proliferation-quiescence intact. (G0 currently carried by p27; making
  phospho-Rb also negative in G0 would need lower inherited CyclinA — minor refinement.)
- **HU→EZH2-in-S boost** — DONE: kDeEZ=5e-5, KmHU_fork=0.35 -> boost **1.37×** (target 1.31×).
  MB G0~26% (target ~22%), proliferation/gradient intact. **Remaining: S high / G2 low (S/G2 proportions).**
- **S/G2 phase proportions** — re-tune fork speed (`kSyDna`) + G2 (`kSyCb`,`a25`) toward MB cycling
  G1 58 / S 21 / G2 21.
- **Absolute CyclinD1 ratios** — apply the recorded HH-optimizer params now that G1 is Cd-decoupled.

Experimental targets compendium: `docs/Ezh2_CcnD1_model_targets.md`. Full history/decisions:
`docs/v44_status.md`.

---

## 6. Build & run

```python
from src.build_model_v44_heldt import build_model_v44
import tellurium as te
rr = te.loada(build_model_v44())          # full model (all modules on)
rr['MYCN_amplification'] = 2.8; rr['Ptch1_copy_number'] = 0.1   # MB context
rr['HU'] = 1.0                            # 10 µM HU
res = rr.simulate(0, 20000, 80000, selections=["time","Dna","pRb","P21","Skp2","EZH2","Cd"])
```
Toggles: `with_ezh2`, `with_hh`, `with_growth`, `with_skp2` (all default True). Override const params:
`build_model_v44(params={"kSyDna": 0.012, "kDeEZ": 0.0003, ...})`. Builder self-check + a cycling
sanity run: `python -m src.build_model_v44_heldt`.

**Key simulation scripts** (`simulations/`): `v44_fork_doseresponse.py` (HU→S), `v44_phase_metrics.py`
(phase proportions), `v44_calibrate_ezh2.py` (EZH2 gradient + HU boost), `v44_calibrate_hh.py` /
`v44_optimize_hh.py` (HH/MYCN ratios), `v44_calibrate_phases.py` (phase balance).

---

## 7. Provenance

- Heldt, Barr, Cooper, Bakal, Novák (2018) *PNAS* 115:2532 — real-time core (BIOMD0000000700).
- Gérard & Goldbeter (2009) *PNAS* 106:21643 — the v42/v43 engine (now superseded).
- v42 (frozen): `src/build_model_v42_mycn.py` — Chahin et al. EZH2/CyclinD1/MYCN/HH integration.
- Retreat tags: `v42-published-baseline`, `v43-milestone`.
