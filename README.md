# Ezh2_CyclinD1_Sphase

A real-time ODE model of the **EZH2 / Cyclin D1 / Hedgehog** cell cycle of cerebellar
granule-neuron progenitors (GNPs) and SHH-subgroup medulloblastoma (SHH-MB), for Chahin et al.

The model explains, mechanistically, how **Cyclin D1** integrates a **mitogenic drive**
(Sonic Hedgehog → Gli, plus MYCN in tumour) against an **epigenetic brake** (EZH2 / H3K27me3),
and how that balance sets the drug responses: **vismodegib arrests MB; an EZH2 inhibitor rescues it.**

> 📄 **Start here for the full picture:** `docs/model_writeup.pdf` — a detailed walkthrough of the
> architecture, every equation with a plain-English explanation, and all figures.
> **Current model + outstanding items:** `docs/h3k27me3_repression_model_state.md`.

---

## Current model: v44 (calibrated)

`src/build_model_v44_heldt.py` — a real-time (minutes, no global clock) model built on the
**Heldt et al. 2018 (PNAS)** explicit-replication cell-cycle core, with the EZH2/Hedgehog biology on
top. Built via `build_model_v44(with_ezh2=True, with_hh=True, ...)`. **Validation: 22/27 targets**
(`simulations/validate_v44.py`).

**The two coupled halves:**

1. **Signalling → Cyclin D1 cascade** (the biology): `SHH → Ptch1 ⊣ Smo → Gli → {Gli1, MYCN} →
   CyclinD1`, with EZH2 repressing *Ccnd1*. Inputs: `SHH`, **`HHi`** (Hedgehog inhibitor =
   vismodegib; 0→1; *renamed from `GDC0449`*), `Ptch1_copy_number` (1 GNP / 0.1 MB), `MYCN_amplification`
   (1 GNP / 2.8 MB), `EZH2i` (EZH2 inhibitor = tazemetostat). MYCN keeps a **mitogen-independent floor**
   so vismo is only a *partial* mitogen withdrawal in MB.
2. **Cell-cycle engine** (Heldt 2018): Rb–E2F restriction point (Cyclin-D1-driven, CKI-modulated,
   growth-gated = the Rb-dilution R-point), Skp2–p27 bistable commitment, **explicit DNA replication**,
   CyclinB/CDK1 mitotic switch, APC/C + division reset. EZH2 is a stable, cell-cycle-coupled protein
   (E2F × cyclins × mitogen-dose) that integrates S-phase duration.

### The H3K27me3 epigenetic layer (repression variants; AUM is now the default)

As of June 2026 the **mean-field AUM H3K27me3 module is the default repression** (`with_h3k27_dilution=True`
by default). The two legacy variants remain available behind flags:

| build flag | repression of CyclinD1 | role |
|---|---|---|
| **`with_h3k27_dilution`** (default) | **mean-field AUM mark** + transcription→PRC2 arm + **leaky** Hill `f0 + (1−f0)/(1+(Mk/K)^n)` | the mechanistic model, **default** (**22/27**) |
| `with_h3k27_dilution=False` | `K_rep/(K_rep + EZH2·(1−EZH2i))` — instantaneous | legacy direct-repression model (22/27) |
| `with_h3k27_memory=True` | explicit mark tracks EZH2; EZH2i → mark decays over ~24–48 h | de-repression **kinetics** variant (takes priority over dilution) |

The AUM "feedback OFF" knob (used by the figure scripts, replacing the legacy `K_EZH2_repression→∞`) is
`f0_mk=1.0` (the leaky floor at 1 → `R≡1`, the mark cannot repress).

The dilution module (Purzner spec + literature review, `docs/H3K27me3_CyclinD1_literature_review.md`) is
the **mean-field reduction of the Berry–Howard A/U/M per-nucleosome model**: EZH2-scaled autocatalytic
read-write methylation + de-novo floor (EZH2i-blocked), a **transcription→PRC2 reciprocal eviction arm**
(nascent *Ccnd1* RNA evicts PRC2 → the bivalent set-point *emerges* and vismo arrest is self-reinforcing),
a **replicative halving at early S** (the cycle period sets the dilution frequency; halving lowers the mark
17–25%), lit-review-paced turnover (de-repression t½ ≈ 16 h), and a **leaky** graded repressor —
H3K27me3 dampens *Ccnd1* but never locks it out (Ser5P Pol II stays; residual transcription at full mark).
All EZH2 repression flows through the dilution-sensitive mark. Calibration: `simulations/calibrate_h3k27_aum.py`.

### Key results (figures in `simulations/`, regenerated via `regenerate_figures.sh`)

- **Vismo arrest + EZH2i rescue** (`fig_v44_fig5_rescue`, `fig_v44_vismo_withdrawal`): MB+vismo arrests
  (MYCN-floor CyclinD1 repressed below threshold); +EZH2i de-represses → rescue; CDK4/6i does not.
- **EZH2i fast-tracks CyclinD1 de-repression** (`fig_v44_ezh2i_rescue_kinetics`): rescue is MB/vismo-
  specific (MB sustains EZH2; mitogen-withdrawn GNP loses it).
- **Mitogen withdrawal = vismo** (`fig_v44_mitogen_withdrawal`): same lesion; the EZH2-feedback effect is
  maximal at *partial* withdrawal.
- **Phase-plane analysis** (`fig_v44_ezh2_phaseplane`): the EZH2⊣CyclinD1 loop is a monostable buffer,
  not a switch; it repositions the proliferation bifurcation.
- **Transient-G0 bifurcation** (`fig_v44_g0_bifurcation`): cyclin D1 / birth-p27 heterogeneity.
- **Vismo resistance is MYCN-driven** (`analyze_vismo_resistance.py`): the resistant CyclinD1 floor is
  ~78% MYCN, ~22% basal, ~0.1% Gli1 — even though a real ~3.5%@24h Gli1 residual exists (it sits below
  the steep Gli→CyclinD1 Hill, so it drives ~no CyclinD1).

---

## Older versions (history; not the current model)

| Version | File | What it is | Status |
|---|---|---|---|
| **v42** | `src/build_model_v42_mycn.py` | Published MYCN–EZH2–CyclinD1–HH model on the Gérard–Goldbeter engine. | **Frozen** (`v42-published-baseline`) |
| **v43** | `src/build_model_v43_checkpoint.py` | Intra-S checkpoint via a phase-gated clock; 12/12 validation. | Milestone (`v43-milestone`) |
| **v44** | `src/build_model_v44_heldt.py` | Heldt-based explicit-replication rebuild + the full EZH2/Hedgehog/H3K27me3 biology. | **Current, calibrated 22/27** |

v44 exists because the Gérard–Goldbeter engine's S-phase duration is a structural invariant (only the
global clock ε can change it); Heldt 2018 (`BIOMD0000000700`) makes S-phase emerge from real replication
fork kinetics. v42 is byte-frozen and never modified; all work is on v44.

---

## Repository layout

```
src/
  build_model_v44_heldt.py   # CURRENT model (flags: with_ezh2, with_hh, with_h3k27_memory,
                             #   with_h3k27_dilution, with_two_step_rb, ...). HHi = Hedgehog inhibitor.
  build_model_v42_mycn.py    # FROZEN published v42 (do not modify; keeps its own GDC0449)
  build_model_v43_checkpoint.py
docs/
  model_writeup.pdf / .md                  # full detailed walkthrough (equations + all figures)
  h3k27me3_repression_model_state.md       # current state, the 3 repression variants, outstanding items
  v44_recalibration_and_search_results.md  # the parameter bake + wide-search results
  FIGURES.md                               # which script generates each figure
simulations/
  validate_v44.py                # 22/27 validation suite (H3K27_DILUTION=1 to test the dilution model)
  regenerate_figures.sh          # regenerate the whole figure suite after any model change
  fig_v44_*.py                   # the figure scripts (rescue, feedback, phase-plane, withdrawal, ...)
  sim_*.py / probe_*.py          # ensemble + mechanism experiments (g0 bifurcation, withdrawal, dilution)
  calibrate_h3k27_dilution*.py   # dilution-module calibration + 6-param search
  analyze_vismo_resistance.py    # MYCN-vs-Gli1 decomposition of the resistant floor
```

---

## Getting started

```bash
python -m venv venv && ./venv/bin/pip install -r requirements.txt   # Antimony/Tellurium + libRoadRunner
./venv/bin/python simulations/validate_v44.py                       # 22/27 validation
bash simulations/regenerate_figures.sh                              # regenerate all figures
```

Build the dilution model: `build_model_v44(with_hh=True, with_h3k27_dilution=True)`.

---

## Status / open items

The v44 production model is **calibrated (22/27)** and the H3K27me3 replicative-dilution module is at
**validation parity (22/27)** with biologically-faithful leaky repression. Outstanding (see
`docs/h3k27me3_repression_model_state.md`):

- **5 standing validation fails** (both models): EZH2 transcript S/G0 gradient (the worst), MB HU S/G2
  folds, MYCN GNP+HHi, MB G2+M duration — plus an expanded wide-search to make these a strict win.
- **DONE:** the AUM dilution model is now the **default** repression (validate 22/27 on it; legacy via
  `H3K27_DILUTION=0` or `with_h3k27_dilution=False`). All registry figures regenerated on it — the ~12
  feedback figures were remapped from the legacy `K_EZH2_repression` knob to `f0_mk` (and the phase-plane
  reworked into the (H3K27me3-mark, CyclinD1) plane).
- **Gli→CyclinD1 coupling:** to make the persistent Gli1 residual offload MYCN from the resistant floor,
  lower `K_Gli_act_CycD` / the Hill exponent (a structural change, recalibration needed).
- **Spec §8 (not yet run):** the cycle-length (T_cc) sweep (the replicative-dilution phenotype) and the
  EZH2i × CDK4/6i synergy.

## Citation

Gérard & Goldbeter (2009) *PNAS* (original engine); Heldt, Barr, Cooper, Bakal & Novák (2018) *PNAS*
(`BIOMD0000000700`, the v44 real-time core); Chahin et al. (in preparation) for the EZH2/CyclinD1/MYCN/
Hedgehog integration.
