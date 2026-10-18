# CyclinD1 as the Hedgehog memory carrier — the memory is the TRANSCRIPT reservoir (2026-07-21)

*Model corrections from Ho, Tsai & Stearns 2020, Curr Biol 30(14):2829 (PMID 32531277; JP in acknowledgements),
supplied by JP. BAKED into the v44 default; validation 28/29. `simulations/fig_cyclind1_memory.py`.*

## The paper's key result
The Hh proliferation decision sits in the **previous** cell cycle. Single cells inhibited (SANT-1) up to ~11 h before
mitosis M1 still entered the next cycle; only >12 h of inhibition blocked entry — even though 4 h of SANT-1 fully
shuts the pathway, so the decision cycle's G1 had no active signalling. **CyclinD1 is the carrier of that memory**:
Ccnd1 is Hh-dependent, CyclinD1 protein is low in S-arrest but accumulates in G2-arrest, and sister cells make
concordant entry/exit decisions — consistent with inheriting a shared CyclinD1 pool through mitosis. Cilia also
persist further into the cycle than the textbook says: ciliated through G1 **and S**, lost by G2.

## Where the memory lives (JP's mechanism): the mRNA reservoir, not the protein
The **CyclinD1 protein is labile** — degraded for most of the cycle, and forms *rapidly* when there is a lot of
CyclinD1 mRNA. So the protein is a fast **readout**, not the memory. The memory is the **cumulative CyclinD1
transcript reservoir**: a slow mRNA pool that integrates Hh over the cycle, carries through mitosis, and is inherited
by daughters. **Gli raises the pool; EZH2 lowers it.**

Because Gli drives *both* the transcript (up) and EZH2 (up), and EZH2 represses the transcript (down), EZH2 forms an
**incoherent feedforward loop (IFFL)** on Ccnd1. When mitogen falls, the direct Gli drive drops fast, but EZH2 —
being slow (τ ≈ 14 h) — falls only gradually, so its repression **lifts with a lag**, *de-repressing* the transcript
and partly compensating for the lost Gli drive. The more EZH2 was holding the transcript down, the more "room" it has
to de-repress — so the transcript **stays up longer** as mitogen declines.

## What was wrong in the model, and the fix
The model's CyclinD1 mRNA (`Cd_mRNA`) had `k_Cd_mRNA_deg = 0.8` → a ~1-minute half-life → the transcript was **slaved
instantaneously to the mitogen**, so there was no reservoir and the decision was effectively instantaneous.

**Fix (baked):** slow the **mRNA** into a **~11–12 h reservoir** — `k_Cd_mRNA_deg` 0.8→0.001 — with the three
transcription rates (`k_Cd_tx_basal`, `k_Cd_tx_Gli_max`, `k_Cd_tx_MYCN`) scaled ×1/800 so the **steady mRNA and
protein levels are held exactly** (the fold ratios, and the `g_prc2`/H3.3 terms that read `Cd_mRNA`, are unchanged).
The **protein stays fast** (`k_Cd_deg = 1.0`, `k_Cd_translation = 0.801` at their string defaults) — so it remains a
labile readout of the slow reservoir. The protein already survived mitosis (`Cd` is *not* reset in `E_div`); now the
*reservoir* it reads from is the thing that genuinely carries through mitosis. **Validation-neutral (28/29)** —
steady-state levels and division rate (~1.06/day) are unchanged; only the *dynamics* become a memory.

## What the model now reproduces (figure — 6 panels)
- **A — two-tier structure:** slow mRNA **reservoir** (the memory) → fast protein **readout** (labile, forms rapidly
  when transcript is high).
- **B — restriction point in the previous cycle:** after Hh-off the reservoir **coasts ~1 more cycle** (mRNA
  t½ ≈ 13 h ≈ the paper's ~11–12 h); the fast/slaved mRNA arrests immediately.
- **C — EZH2 incoherent feedforward buffers the dose-response:** as steady mitogen falls, EZH2 de-represses (rep
  factor 0.26→0.44), holding the transcript **~1.7× higher** than a frozen-repression control — a flatter, buffered
  mitogen-dependence.
- **D — the reservoir "stays up" after the input drops:** the Gli input resets in <1 h, but the transcript reservoir
  relaxes with t½ ≈ 10 h and the **protein follows the reservoir, not the input** — so the cell coasts.
- **E — widened cilium duty cycle:** `cilium_gate` (default off) — Shh reception ON through G1/S, OFF at G2
  (CyclinB-gated), wider than the textbook G1-only window.
- **F — robustness:** proliferation is flat across the ciliated-fraction scan — intermittent reception is buffered by
  the reservoir (the paper's argument that ~40% ciliated at any moment still lets most cells enter).

## Caveats (do not over-read)
- SMB55 lacks wild-type Ptch1 (constitutive, ligand-independent signalling), so the paper says nothing about the
  Ptch1 negative-feedback/adaptation term — that stays a scanned free parameter (see `realistic_mitogen_input`).
- The ~40% ciliated fraction (SMB55, 2-D) is well below ~94% in vivo (EGL) — treat the duty-cycle fraction as
  system-dependent (a sensitivity scan), not a fixed constant.

## Independent support
The model's **lower CyclinD1 under HU** matches JP's data: HU prolongs S → EZH2 accumulates in the elongated S phase
(the Cdh1-gated S-phase EZH2 boost) → EZH2 represses CyclinD1. Independent support for the EZH2⊣CyclinD1 axis — and
consistent with the IFFL: more EZH2 (from a longer S) pushes the transcript reservoir down.
