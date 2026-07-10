# EZH2 cycle-gating: the Palbociclib mRNA-drop vs within-cycle-gradient tension

**Status (2026-07-10):** open issue, documented here because it is subtle. The model is left at the
clean PRC2 calibration (25/28); the `EZH2 Palbo mRNA drop` target is intentionally left **failing** as a
marker of the under-gating, to be resolved by **Feature C** (commitment dynamics), NOT by an EZH2-module hack.

## The observation

The paper (Chahin Fig 4) reports that **CDK4/6 inhibition (Palbociclib) drops EZH2 mRNA ~56%** — because EZH2
is a direct Rb–E2f target, so blocking CDK4/6 → Rb hypophosphorylated → E2f down → EZH2 transcription down.

The current v44 model (PRC2 mechanism) drops EZH2 mRNA only **~27%** under `cdk46i` (kPhRbCd=0). The writer is
**under-cycle-gated** — it stays too high in arrest. This matters for JP's Point 2 (the writer should shut off
in a CKI-arrested MB cell so the Gli/Jmjd3 eraser wins and the CyclinD1 mark is stripped).

## Why you cannot just crank up the E2f-gating

EZH2 transcription is `kEZbas + (kEZbas_Cd + kEZE2f·E2f_gate·CycEA_gate)·Cd/(Kez_cd+Cd)`. The `kEZbas_Cd` term
is a **mitogen-dose-scaled but CYCLE-FLAT baseline** (E2f-independent) — it was added to flatten the within-cycle
EZH2 mRNA gradient to the measured value. It is also what keeps EZH2 up in arrest.

Shifting the balance toward the E2f-gated part (raise `kEZE2f`, lower `kEZbas_Cd`) to reach the 56% Palbo drop
**over-steepens the within-cycle gradient** past its data ceilings:

| config | Palbo mRNA drop | EZH2 transcript S/G0 (≤2.5) | EZH2 protein G2/G0 (≤2.0) |
|---|---|---|---|
| current (clean PRC2) | 0.73 (27% drop) ✗ | 1.77 ✓ | 1.81 ✓ |
| more E2f-gated (config A) | 0.41 (59% drop) ✓ | **3.41 ✗** | **2.00 ✗** |
| optimizer best compromise | 0.56 (44% drop) ✓ (edge) | 2.57 ✓ (edge) | 1.94 ✓ (edge) |

The full 13-parameter optimizer (`optimize_prc2.py`, EZH2 + PRC2 params jointly, 28 targets) found a compromise
at **26/28 but fragile** — Palbo only reaches 44% (not 56%) and transcript/protein/CyclinD1 are all pinned to
their pass-ceilings. Not baked.

**Root cause:** the within-cycle EZH2 gradient (S/G0) and the sustained-arrest drop (Palbo) are driven by the
**same E2f×CycE/A gate**, so they move together — and, critically, in the current model a cycling cell's **G1
has about as little E2f as a Palbo-arrested cell**, so "transient within-cycle G0" ≈ "sustained arrest" as far as
EZH2 is concerned. You cannot deepen one without steepening the other.

## Model structure is correct — the entities are NOT merged (verified)

JP flagged a worry that CyclinD1 transcription and CDK4/6 might have been conflated. They are not:

```
species ... Cd_mRNA, Cd in Cell;                       # transcript and protein are DISTINCT species
CycD1_transcription:    => Cd_mRNA;  (basal+Gli+MYCN) × [EZH2/PRC2 repression]   # EZH2 hits the TRANSCRIPT only
Cd_translation:         Cd_mRNA => Cd_mRNA + Cd;  k_Cd_translation·Cd_mRNA        # separate translation step
Cd_degradation:         Cd => ;  k_Cd_deg·Cd
```

- **CyclinD1 transcript (`Cd_mRNA`) and protein (`Cd`) are separate species**, each with its own reactions.
- **EZH2/PRC2 represses CyclinD1 transcription only** (multiplies `CycD1_transcription`), never the protein.
- **Palbociclib = `kPhRbCd = 0`** — zeroes only the CyclinD-CDK4/6 Rb-phosphorylation rate (code comment:
  *"DISTINCT from CDK4/6i (palbociclib), modeled as kPhRbCd=0"*). It does not touch CyclinD1 transcription or EZH2.
- So Palbo affects EZH2 **only through the natural cascade** CDK4/6 → Rb → E2f → EZH2.

## Rejected fix: gating EZH2 by CDK4/6 activity

A tempting fix is to gate the EZH2 baseline by `kPhRbCd·Cd` (CDK4/6 activity), which is 0 in Palbo. **Rejected
(JP):** that artificially **merges** Palbo's target (CDK4/6) into EZH2 transcription, short-circuiting the
cascade. Palbo must reach EZH2 only via CDK4/6→Rb→E2f, which the model already does.

## The correct resolution: Feature C (commitment dynamics)

The fix is **not in the EZH2 module.** It is that **committed/cycling cells should retain more E2f in G1 than a
sustained CDK4/6-blocked cell** (Spencer 2013 / Cappell 2016 / Fan-Meyer 2021: once past the restriction point,
CDK2/E2f stays elevated and self-sustaining). With that cell-cycle behaviour:

- a **transient within-cycle G1** keeps its E2f → EZH2 stays up → moderate within-cycle gradient (S/G0 ~2);
- a **sustained Palbo arrest** (E2f driven to zero and held) → EZH2 drops the full ~56%.

The two measurements decouple naturally, with **no change to the EZH2 module** and no merging of entities. This
is exactly Feature C — making G0↔G1/S entry CyclinD1/CDK2-controlled (with the population distributions for
stochastic transient-G0 vs immediate re-entry) — so it is deferred to that work.

## Current state / TODO

- Model = clean PRC2 (25/28): H3K27me3 MB/GNP 0.507 (ChIP), CyclinD1 MB/GNP 5.06, EZH2i 2.19, all centered.
- `EZH2 Palbo mRNA drop (~0.44)` left **failing** (0.73) as a visible marker of the writer under-gating.
- Fix via Feature C; re-check this target once the commitment E2f-retention is in.

Related memory: `ezh2-tracks-mitogen-rb-e2f`, `v44-quiescence-bifurcation-direction`, `v44-gli-jmjd3-eraser-reconception`.
