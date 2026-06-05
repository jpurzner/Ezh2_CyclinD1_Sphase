# Data used to parameterize / validate the v44 model

Every quantitative constraint the model was fit to or checked against, with its source. These are
encoded as the target checks in `simulations/validate_v44.py` (run it to score the current model).
The model is calibrated against **ratios and behaviors**, not absolute concentrations (species are in
arbitrary units). Full target compendium with caveats: `docs/Ezh2_CcnD1_model_targets.md`.

## A. Bulk RNA‑seq — steady‑state expression ratios (cell‑type / drug context)

| Measurement | Experimental value | Constrains | Source |
|---|---|---|---|
| Cyclin D1: GNP+HHi / GNP | 0.14 (86% ↓) | Gli‑driven fraction of CyclinD1; HH pathway gain | RNA‑seq |
| Cyclin D1: MB+HHi / MB | 0.40 (60% ↓) | MYCN‑driven (HHi‑resistant) fraction of CyclinD1 | RNA‑seq |
| Cyclin D1: MB / GNP | 5.07× (p<0.05) | overall MB CyclinD1 drive (soft target — see note 1) | RNA‑seq, Fig. 4I |
| Mycn: GNP+HHi / GNP | 0.78 (22% ↓) | Gli‑dependent vs autonomous Mycn synthesis | RNA‑seq |
| Mycn: MB+HHi / MB | 0.86 (14% ↓) | autonomous (amplified) Mycn fraction | RNA‑seq |
| Mycn: MB / GNP | ~2.8× | `MYCN_amplification` (MB context) | RNA‑seq |
| Gli1: GNP+HHi | >99% ↓ | HHi (vismodegib) block on Smo→Gli | RNA‑seq |
| Ezh2: MB / GNP | 2.05× (p<0.001) | cell‑cycle‑coupled EZH2 (soft — see note 1) | RNA‑seq, Fig. 4J |
| Ezh2: GNP+HHi | 25–47% ↓ | EZH2 falls when cycling stops | RNA‑seq |

## B. EZH2 perturbation → Cyclin D1 (the feedback strength, `K_EZH2_repression`)

| Measurement | Experimental value | Constrains | Source |
|---|---|---|---|
| Cyclin D1: **Ezh2 cKO** GNP / WT | log2FC 1.301 → **2.46×** | EZH2 ⊣ CyclinD1 repression strength | RNA‑seq, Fig. 3A |
| Cyclin D1: GNP + Tazemetostat / DMSO | 2.2× | same (catalytic inhibition) | qPCR, Fig. 3C |
| Cyclin D1: MB + Ezh2‑OE / control | significant ↓ | direction of the repression | qPCR, Fig. 3E |
| Cyclin D1: MB + Palbociclib / DMSO | +57% (1.57×) | CDK4/6i ↑ CyclinD1 (via ↓EZH2) | qPCR, Fig. 4F |
| Ezh2, E2f1: MB + Palbociclib | −55.6%, −56.7% | CDK4/6→Rb‑E2F→EZH2 arm | qPCR, Fig. 4F |

> **Note 1 — the cKO 2.46× / Tazemetostat 2.2× are LOWER BOUNDS on the true feedback.** The model
> gives ~3× de‑repression (`K_EZH2_repression = 0.5`). The bulk cKO under‑reads the per‑cell effect
> because (i) Math1‑Cre is incomplete, so WT escapers dilute the signal — at 70–80% recombination the
> bulk 2.46× implies a **per‑cell 2.8–3.1×**, matching the model — and (ii) the Ezh2‑flox deletes the
> SET (methyltransferase) domain but leaves EZH2 protein (residual non‑catalytic repression), and
> H3K27me3 dilutes only with division, so even the per‑cell fold is a floor. The model's 3× is the
> recombination‑corrected value and is plausibly conservative (`fig_v44_ckoo_feedback_strength.py`).

## C. scRNA‑seq — Ezh2 transcript by cell‑cycle phase (the EZH2 cell‑cycle gate)

| Measurement | Experimental value | Constrains | Source |
|---|---|---|---|
| Ezh2 cycling / G0 (transcript) | ~2.0× (S/G0 ≈ 1.8–2.5 across Sepp, Vladoiu, Ocasio, Luo) | EZH2 transcription gated on E2f×(CycE+CycA) | scRNA‑seq, Fig. 4A,B |

## D. Immunofluorescence — Ezh2 protein by phase / under arrest drugs

| Measurement | Experimental value | Constrains | Source |
|---|---|---|---|
| Ezh2 protein G2 / G0 (DMSO) | 1.48× | EZH2 protein accumulation across S/G2 | IF, Supp. 6D |
| Ezh2 in S: HU 10 µM / DMSO | 1.31× (alt 1.22×) | EZH2 integrates S‑phase **duration** (HU→fork speed) | IF, Fig. 4G |
| rShh dose → Ezh2 | dose‑dependent ↑ | Hh→CyclinD1→cycle→EZH2 (whole loop) | IF, Fig. 4H |

## E. Cell‑cycle phase proportions (G0/G1/S/G2)

| Dataset | Values | Constrains | Source |
|---|---|---|---|
| MB, DMSO (renormalized) | G0 24.8 / G1 43.4 / S 15.7 / G2 16.1 % | phase‑duration balance; G0 by p27 marker | Fig. 4G / Supp. 6, Sec. F |
| MB, HU 10 µM (fold vs DMSO) | S ×1.36 ↑, G2 ×0.23 ↓ | HU lengthens S (replication‑fork model) | Sec. F |
| GNP + Tazemetostat (fold) | G0 ×0.92, G1 ×2.03, G2 ×1.44 | EZH2i shrinks G0 / drives cycling | Fig. 3B |
| MB + Ezh2‑OE (fold) | G0 ×2.84, G1 ×0.51 | EZH2‑OE traps in G0 | Fig. 3F |

## F. Dynamics

| Measurement | Value | Constrains | Source |
|---|---|---|---|
| Cell‑cycle period, GNP | ~22–23 h (soft) | overall timescale (growth rate `mu`, fork speed) | Published (Nakashima 2015, ref. 70) |

---

## Optimized / fit parameters

A constrained search (`simulations/v44_optimize_hh.py`, `v44_calibrate_*.py`) over the
HH/MYCN→CyclinD1 and EZH2 parameters, with a hard guard requiring the figure‑critical behaviors
(GNP cycles only with SHH; GNP+HHi & MB+HHi arrest; MB+HHi+Ezh2i rescues; CDK4/6i not rescued).
Key fit values (full list in `src/build_model_v44_heldt.py`):

- MYCN→CyclinD1 cooperative Hill (`n_MYCN_Cd = 3`, `k_Cd_tx_MYCN`, `K_MYCN_Cd`) — creates the
  threshold that makes MB Mycn HHi‑resistant while GNP Mycn is negligible.
- `MYCN_amplification = 2.8` (MB) — from Mycn MB/GNP ratio.
- `K_EZH2_repression = 0.5` — EZH2 ⊣ CyclinD1 feedback strength (~3× de‑repression; see note 1).
- EZH2 gate (`kEZbas`, `kEZE2f`, `K_Ce_EZ`, `K_Ca_EZ`) and turnover (`kDeEZ`) — the ~2× transcript
  gradient and the 1.3× HU‑in‑S protein boost.
- `KmHU_fork`, `vmin_fork` — HU→replication‑fork‑speed (S‑phase lengthening).

**Soft / not‑enforced targets** (documented limits of a single‑cell threshold model fit to bulk
population means): CyclinD1 MB/GNP 5.07× and Ezh2 MB/GNP 2.05× (forcing them breaks the MB+HHi arrest
the rescue depends on); absolute S/G2 proportions (S is biochemically long in the real‑time core).
