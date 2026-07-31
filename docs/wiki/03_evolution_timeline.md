# Evolution Timeline

> The chronological spine of the project — from the first published-BioModels mashups (Jan 2026) through the Heldt-core rebuild to the current cell-type-split, two-cyclin, dynamic-CDKI default (Jul 2026). Each milestone records **what changed, why, and the validation state at that moment**. Other wiki pages hang their detail off this backbone: architecture in [Model Architecture](02_model_architecture.md), the surviving mechanisms in [Mechanisms Explored](04_mechanisms_explored.md), the graveyard in [Things Tried & Abandoned](05_things_tried_and_abandoned.md), the scoring history in [Calibration & Validation](06_calibration_and_validation.md), and the sources in [Origins & Foundations](01_origins_and_foundations.md).

There are **two repositories** in the story. The old *foundation* repo (`ezh2_cyclind1_sym`, GitHub `Ezh2_Hh_Ccnd1_model`) held versions **v9–v31** (Jan 2026), built by wiring together published BioModels. The *active* repo (`Ezh2_CyclinD1_Sphase`, GitHub `Ezh2_CyclinD1_Sphase`) holds the S-phase checkpoint rebuild **v42 → v43 → v44 → …** (Apr–Jul 2026, 233 commits). The version numbers do **not** continue across the two repos — v42 is a fresh start on a new core, not "v31 + 11."

A quick validation-count legend (the denominator grows as targets were added): sym era used ad-hoc pass/fail goals; v43 hit **12/12**; v44 ran **22/27 → 28/28**, then **28/29** (target spec formalized), then **30/32** (cell-type split + two-cyclin + CDKI).

---

## Era 0 — Foundations: the published-model mashup (`sym` repo, ~Jan 2026)

The goal from day one: a mechanistic ODE model where **Hedgehog drives CyclinD1**, **EZH2/H3K27me3 gates CyclinD1**, and the readout is **cell-cycle proliferation vs G0 arrest** in GNPs and SHH-MB. The first strategy was to bolt a Hedgehog/Gli front-end and an EZH2/chromatin module onto a *published, already-oscillating* mammalian cell-cycle model. See [Origins & Foundations](01_origins_and_foundations.md) for the full source list (Gérard 2009/2010, Novák–Tyson 2004, Abroudi 2017, Tyson 2001, Yao 2008 Rb-E2F, Schwarz 2018, Goldbeter, Schoeberl 2002).

| Version | Date | What changed | Why | Validation state |
|---|---|---|---|---|
| **v9–v13** | early Jan | Added Gli1/Ptch1 feedback to a cyclin oscillator; fast (6h) Gli kinetics; cell-cycle speed tuning to a 24h period | Get HH to modulate an existing oscillator | Gli stable during proliferation, ~24h period |
| **v14–v16** | Jan | Introduced the **bivalent chromatin model** (H3K27me3 repressive/slow-turnover, H3K4me3 active/fast); tried Gli1 autoregulation at 0.1 (v15) then 0.0 (v16) | Chromatin memory + amplify the Gli signal | v15 runaway Gli; v16 stable but too weak (1 division) |
| **v17** ✓ | 2026-01-15 | **Balanced Gli1 autoregulation = 0.03** (the "critical parameter"); slow H3K27me3 demethylation (k=0.02, ~35h half-life) vs fast H3K4me3 (k=0.4) | 0.03 = ~3% self-amplification, signal boost without a self-sustaining loop | **All 4 goals met**: 0 divisions no-HH; 23.7h period with HH; Gli stable (0.003→0.358); 50× HH fold. `src/build_model_v17_balanced_gli1.py` |
| **v21–v27** | Jan | Reworked Gli architecture repeatedly | Fix biology | Documented as having "fundamental architectural errors" — no proper HH-dependent control; v26 exploded, v27 GliR dominated both conditions |
| **v28** | 2026-01-16 | **Correct Gli biology**: Gli2/Gli3 are *repressors by default* (full-length), HH converts repressor→activator via Smo; Gli1 obligate activator/readout; basal CyclinD1 term added | Prior models had Gli1 driving CyclinD1 directly (wrong); real pathway is a rep→act switch | **BUG**: no Ptch1 differentiation (0.99×) — pathway didn't switch |
| **v29** | 2026-01-16 | **Fixed the SHH-Ptch1 complex degradation** + made SHH a boundary species; 18-param differential-evolution fit | The missing complex-internalization reaction was why the switch didn't fire | Pathway "textbook-perfect": Gli1 5.9×, Smo 3.7×, CycD 7.9× fold. **But S-phase arrest — 0 divisions** (CycB can't accumulate; Cdh1 stays active) |
| **v30** | 2026-01-16 | **CyclinD→p27 control**: Skp2 synthesis Hill on E2F×CycD-CDK4; dual p27 degradation (Skp2 + direct CycD); CycB synthesis n=3 on CycD | Get differential proliferation without a perfect oscillator | Skp2 19× fold, p27 0.47×; **1 division with HH vs 0 without** — proof-of-concept, but oscillations wouldn't sustain |
| **v31** | 2026-01-17 | Tried to swap in the **Goldbeter 3-variable minimal oscillator** | Get sustained cycling | **FAILED**: CVODE `CV_TOO_MUCH_WORK` at t=0.74h. Michaelis-Menten stiffness (K≈0.01) × HH-pathway timescale mismatch. Recommendation: switch to Tyson 2-variable |

**Era-0 verdict.** The `sym` repo proved the *biology* (HH rep→act switch → CyclinD1 → p27 → proliferation, plus bivalent-chromatin memory) but never delivered robust, sustained, concentration-sensitive cycling from a bolted-on published oscillator. The recurring wall — a stiff or bistable oscillator that either stuck in S-phase or wouldn't reset — is what motivated abandoning the mashup and rebuilding on a purpose-built S-phase core. GitHub note: the old clone is **stale**; the local `sym` tree is the authority.

---

## Era 1 — v42: the fresh S-phase rebuild begins (Apr 2026)

The active repo opens. The scientific pivot: the project's *core claim* is a **concentration-dependent intra-S checkpoint** (EZH2 integrates S-phase duration; HU slows forks → longer S → more EZH2), and the Gérard–Goldbeter-style oscillators used in Era 0 **structurally cannot** produce concentration-dependent S-phase. A new core was needed.

| Milestone | Commit date | What changed / why | Validation |
|---|---|---|---|
| **Initial v42 release** | 2026-04-22 | "MYCN–EZH2–CyclinD1–Hedgehog cell cycle model" first commit; hard-coded paths replaced; figure reproduction documented | baseline working model |

v42 carries the HH/MYCN→CyclinD1 module that is later ported forward into v44 (`with_hh=True`). See [Repos, Directories & Reproduction](09_repos_directories_reproduction.md).

---

## Era 2 — v43: the intra-S checkpoint (late May 2026)

| Milestone | Commit date | What changed / why | Validation |
|---|---|---|---|
| Intra-S checkpoint extension plan | 2026-05-29 | Scope the `eps`-in-S mechanism | — |
| **v43 milestone** | 2026-06-02 | **`eps`-in-S intra-S checkpoint** — a phenomenological S-phase progress variable; analysis docs | **12/12 validation** — the retreat/rollback point, tagged `v43-milestone` |

v43 is the last version with the abstract `eps` S-phase variable. It is the safe fallback tag the v44 rebuild retreats to if the explicit-replication branch breaks.

---

## Era 3 — v44: explicit-replication rebuild on the Heldt 2018 core + three structural redesigns (June 2026)

The biggest single re-architecture. Branch `v44-explicit-replication`. The model is rebuilt in **real time (minutes, no `eps`)** on **Heldt, Barr, Cooper, Bakal & Novák 2018** (PNAS, BIOMD0000000700), which has an explicit bistable Rb-E2F restriction point. Builder: `src/build_model_v44_heldt.py`, `build_model_v44(hu=, with_ezh2=, with_hh=)`. Full anatomy in [Model Architecture](02_model_architecture.md).

**Core added on top of Heldt (2026-06-02):**
- **Explicit DNA replication** (`Dna` synthesized by active forks `aRc` at fork speed `kSyDna`).
- **CycB/CDK1 (MPF) mitotic switch** — Cdc25/Wee1 hysteresis.
- **CHK1 intra-S checkpoint** gated on active forks → holds mitosis until `Dna→1`.
- **Mitotic APC/Cdc20 + division-reset event** → sustained cycling.
- **HU → fork speed** (`vfork(HU)` scales `kSyDna`).
- **EZH2 layer** — `Cd` (CyclinD) made dynamic and EZH2-repressible; EZH2 transcription gated on `E2f×(CycE+CycA)`; EZH2 a stable protein reset by dilution at division so it **integrates S-phase duration**.

**Primary goal achieved:** concentration-dependent S-phase — S 5.5h→16.4h across HU 0→2 with **G1 flat**, mitosis always waits for `Dna=1`, no arrest; **EZH2-in-S boost 1.18× (HU=1) / 1.59× (HU=2)** vs experimental 1.31× at 10µM.

Then three **structural redesigns**, each fixing a wall the previous exposed:

| Structural redesign | Commit date | What / why | Validation state |
|---|---|---|---|
| **#1 — real G2 phase** | 2026-06-02 | Reworked mitotic switch so CyclinB-CDK1 is synthesized (inactive, Tyr15-P) only as replication completes (`g2gate` on `Dna`), building *during* G2 instead of pre-stocked and flipping instantly | G2 2%→~20%; period 9→~11–13h; HU redistribution correct (S↑, G2↓) |
| **#2 — growth/size-gated restriction point** | 2026-06-02 | Added a cell `mass` growing exponentially, halving at division, with a **size gate on S-entry** (origins fire only once `mass≥M_size`). Commitment still needs mitogen (CyclinD), but committed cells then *wait for size* → **G1 length = growth time, DECOUPLED from CyclinD level** | Period ~22h; G1 ~48% and nearly identical GNP vs MB despite different Cd; HH-dependence preserved. This is grounded as the Zatulovskiy **Rb-dilution** R-point |
| **#3 — Skp2-p27 feedforward switch** | 2026-06-03 | Skp2 a dynamic E2F target degraded by APC/Cdh1 (Heldt's `C1`); builds the canonical Skp2→p27→Rb→E2F feedforward; p27 reset high at division; `kPhRbCd` raised 0.2→0.5 | Genuine bistable R-point; proliferation-quiescence correct (GNP+SHH cycles, GNP−SHH/+HHi quiesce, MB cycles) |

Why #2 mattered so much: the HH-optimizer had found params giving CyclinD1 MB/GNP ~4× (toward the measured 5.07×) but they made MB cycle too fast (short G1) and crashed CVODE — the model **coupled CyclinD1 level to G1 length**, while the data demand MB have *both* high CyclinD1 *and* a substantial G1. The size gate decoupled them. This is the recurring lesson: level-vs-duration coupling is the model's structural pressure point.

---

## Era 4 — Between-condition calibration: HH saturation, the p16 brake, EZH2 recalibration (June 2026)

With structure complete, a single whole-model calibration pass replaced piecemeal tuning. See [Calibration & Validation](06_calibration_and_validation.md) and [Data & Evidence](07_data_and_evidence.md).

| Milestone | Commit date | What / why | Validation |
|---|---|---|---|
| HH/MYCN optimizer (ratio fit) | 2026-06-02 | 400-sample + Nelder-Mead over 7 HH/MYCN→CyclinD1 params; **CyclinD1 MB/GNP 1.98→4.16** | Revealed the Cd-level/G1-length coupling → drove redesign #2 |
| CKI recalibration | 2026-06-03 | CKIs abundance-proportional, p27→CDK4/6; S-phase shortened 10→3h; phase-proportion validation switched to **count-fractions** | 23/28 |
| GLI1-autoregulation epigenetic memory | 2026-06-10 | `Gli1_epi` slow memory → Gli1 residual after vismodegib (2.3% @24h) without inflating proliferation | 22/28 |
| **p16 competitive CDK4/6 brake** | 2026-06-05 | Strong vismodegib arrest + **EZH2i rescue** matching pRb data; rescue reframed as a **population CyclinD1-threshold effect** | search-tuned to pRb 100/22/69 |
| Gli→Ptch1 negative feedback; MB = broken loop | 2026-06-05 | Canonical Ptch1 feedback; data pins f≈0.1 as the broken-feedback level | high MB Gli1 explained |
| p21 (MB CDK2 brake) + p18/Cdkn2c (INK4, GNP baseline tone) | 2026-06-06 | Widen the CDK-inhibitor panel to INK4 + CIP/KIP family tones | figures regenerated on p16+p18+p21 model |
| **EZH2 synthesis gains a mitogen-dose term** | 2026-06-15 | `EZH2_tx ×= Cd/(Kez_cd+Cd)` — EZH2 tracked only binary commitment before (plateaued); Fig 4H needs a wide dose-response | EZH2 climbs 1.07→2.30 over SHH 0.5→8; **EZH2 MB/GNP 2.05** (was 1.1); MB CyclinD1 fold reaches ~4.65 |
| **Bake wide-search best** (raised threshold) | 2026-06-16 | `kPhRbCd=0.35` (raised commitment threshold, wider GNP sub-threshold range, rescue preserved), 14 baked params; MB+HHi/MB 0.20→0.14 | **22/27**, rescue intact. Docs: `v44_recalibration_and_search_results.md`. Traded MB G2+M duration; HU folds / EZH2-S-gradient / MYCN-GNP+HHi remained unscored misses |
| README overhaul + GDC0449→HHi rename | 2026-06-24 | Cleanup to current state (calibrated v44; HHi; H3K27me3 layer) | — |

Meanwhile, **v45 (stochastic-commitment)** was prototyped in parallel (2026-06-13 → 2026-06-14, then archived to `archive/v45/`): a bistable CDK2-p27 toggle + mitogen-gated bootstrap + mass-cap, **17/17** then **20/20**, with EZH2 re-added as a dynamic species and a MYCN floor protecting MB from Hh-withdrawal arrest. Decision: v44 stays the working model; v45's toggle is "cosmetic" relative to v44's existing size-gate + Skp2-p27. See [Things Tried & Abandoned](05_things_tried_and_abandoned.md).

---

## Era 5 — The H3K27me3-memory explorations (July 2026)

A long, mostly-exploratory arc asking: *what kind of epigenetic memory can the mark carry, and is a strong chromatin memory even compatible with the data?* Most of this did **not** bake; it defined the boundaries.

| Milestone | Commit date | Finding | Bake? |
|---|---|---|---|
| Promote mean-field AUM H3K27me3 module to default repression | 2026-07-02 | Explicit-mark repression becomes the default topology | default |
| H3K27me3→CyclinD1 response characteristics | 2026-07-05 | Register the exploratory set | — |
| "Battery"/charge-persist-discharge; memory phase maps; 2-D `kDeEZ×k_w_mk` sweeps | 2026-07-06 → 07-09 | The mark's ramp effects require a **slow-turnover EZH2** (a threshold); memory is chromatin-intrinsic | exploration |
| **Strong chromatin memory EXCLUDED** by measured folds | 2026-07-09 | Co-optimization: strong memory conflicts with the measured MB/GNP fold; ~450 candidates found **zero feasible**; two data walls (vismo 0.144 + MB/GNP 5.07) | negative result (bound) |
| **Gli→Jmjd3/Kdm6b eraser reconception** | 2026-07-09 | JP ChIP: MB has ~½ GNP H3K27me3 on CyclinD1. EZH2-writer (cycle-coupled) vs Gli-Jmjd3-eraser (mitogen) *race*; `k_jmjd3_gli`. Flips ChIP direction to **MB=0.52× GNP** while keeping 5.07 fold | baked |
| Formal PRC2-complex repression | 2026-07-10 | Replaces `w_ezdir` with a PRC2-Hill; more biological; EZH2i acts via `K_prc2/a_rw_prc2` | **25/27**, all EZH2 anchors nailed |
| **Read-write amplifier** (coop n_rw=2) | 2026-07-10 | Make the mark a functional EZH2 amplifier; loosened noisy targets (JP: don't exact-match noisy RNA-seq/flow/ChIP — band penalties) | **28/28** clean |
| Feature C: commitment carryover | 2026-07-10 | Spencer carryover across division | 27/28 (superseded) |
| **Serial me1/me2/me3 methylation chain PROMOTED** | 2026-07-11 | Rebake around explicit me0→me1→me2→me3; does *not* make steady-state timing matter (read-write buffers) | **28/28**, `with_h3k27_chain` default True |
| Possibility-space: bistable latching mark (v44.1) | 2026-07-11 | Extended bounds unlock a latch, but latch + calibrated are **mutually exclusive** (latch breaks GNP Hh-dependence) | exploration only |

The four surviving mark mechanisms (nerfed Gli-JMJD3 eraser `k_jmjd3_gli=0.005`, coop read-write `n_rw=2`, continuous S-dilution `k_dil_S=0.693`, H3.3 turnover `k_h33=0.003`) become the default; effective decay `δ_eff = del_mk + k_jmjd3_gli·Gli1 + ln2/Tc`.

---

## Era 6 — Two-step Rb baked; the mother-G2 negative result (2026-07-12 → 07-14)

| Milestone | Commit date | What / why | Validation |
|---|---|---|---|
| Set up two-step Rb (mono/hyper) | 2026-07-12 | Literature (Narasimha/Sanidas): CyclinD only *mono*-phosphorylates Rb, CyclinE-CDK2 *hyper*-phosphorylates. Single-step let CyclinD directly release E2F (wrong R-point) and saturated pRb | candidate, default unchanged |
| **PROMOTE two-step Rb to default** | 2026-07-13 | Daughter inherits mono/hyper Rb via `f_commit_carry`; **pRb restored as a valid G0 marker** (norm pRb 0.48 G0 vs 0.82 committed, was saturated 0.88/0.98) | **27/28** default; single-step opt-in stays 28/28. Lone miss: MB HU-S fold (two-step checkpoint interaction) |
| Mother-G2 p27 integrator | 2026-07-13 | Daughter birth-p27 set by mother's G2 mitogen history (Spencer/Min) | **NEGATIVE**: quiescent fraction saturates ≈1.0; birth-p27 *floor* + frequency confound. Fix identified (low-pass tracker with low floor), not built |
| **Joint re-optimize two-step → 28/28** | 2026-07-14 | Genuine **5.6× MB/GNP fold**, de-crutched MB G0 to 12%; strong EZH2i (≥2.2) data-excluded | **28/28** — the "de-crutched" milestone |
| EZH2-concentration convention promoted | 2026-07-13 | Drop the EZH2 division-halving; re-fit `kTlEZ ×0.78` | 26/28 (later reversed in Era 7) |
| Population G0 layer calibrated to flow gate | 2026-07-14 | GNP ~10% / MB ~20% G0 via count-fraction snapshot | matched |

---

## Era 7 — Cdh1-EZH2 S-phase boost; EZH2i catalytic mechanism; transient-G0 / CDK4/6i reframe (2026-07-15 → 07-19)

| Milestone | Commit date | What / why | Validation |
|---|---|---|---|
| **EZH2i catalytic mechanism** (PR #5) | 2026-07-15 | Catalytic inhibitor blocks *writing* only, not the PRC2 complex; repression = occupancy (`PRC2_rep`); de-repression gradual (~24h) + partial (~1.71×) | merged |
| Transient-G0 reframe + CDK4/6i dose-response | 2026-07-15 | CDK4/6-activity-gated reversible G0, **pRb marker**, birth-p27 crutch retired, p27 split; CDK4/6i reframed as a **dose-response** (K≈1.5µM, two measured palbo doses 1µM sub-sat / 5µM sat) — **resistant-clone framing retired** | 27/28 (branch) |
| Structure analysis: v44 as a 3-tier machine | 2026-07-14 | Read the model as fast slaved mitogen input → 35-species oscillator core → slow epigenetic/growth memory; 6 load-bearing dials, `K_EZH2_repression` dead | doc |
| **Cdh1-EZH2 S-phase boost** (default bake) | 2026-07-16→19 | EZH2 = APC/Cdh1 substrate (`kDeEZ_C1`) + dilution convention → elongated/HU-S cells accumulate EZH2; **HU-EZH2-in-S 1.009→1.24**; **reverses** the Era-6 EZH2-concentration promotion | **27/28** |
| PRC2/transcription paradox resolved | 2026-07-17 | Global EZH2 (E2F) vs local eviction; a stable mark via bistable alleles | doc |

---

## Era 8 — mRNA reservoir, noisy Hh input, and the abandoned 16h re-calibration (2026-07-20 → 07-26)

The memory mechanism is reconceived (JP's hypothesis) and the input is made evidence-driven; a hard attempt to abandon the 22h cycle is documented and *not* baked.

| Milestone | Commit date | What / why | Validation |
|---|---|---|---|
| **EZH2 governor re-measure** — `k_Cd_tx_basal ×0.5` baked | 2026-07 | Full de-repression now **Shh-DEPENDENT** (matches cKO: transcript up, no division without Shh, differentiates); thresholds ON 0.18 / EZH2i 0.06 / OFF 0.01; EZH2 raises threshold ~15×; Shh-dependence structural (low basal). Supersedes old 2.6× | **28/28** |
| **CyclinD1 memory = the mRNA transcript RESERVOIR** (default bake) | 2026-07-24 | JP: the Hh memory is the cumulative CyclinD1 *transcript* reservoir, not the protein. `k_Cd_mRNA_deg 0.8→0.001` (mRNA τ≈16h reservoir) + tx ×1/800 to hold levels; protein stays fast; EZH2 = incoherent feedforward → transcript coasts ~1 cycle after Hh-off | **28/29** (target spec now formalized to 29). `fig_cyclind1_memory.py` |
| **Realistic/noisy mitogen input** (all OFF by default) | 2026-07-24 | Tonic + slow OU drift + fast OU noise, cilium duty cycle, weak ~20h Ptch1_slow adaptation. **Key correction**: Gli is quasi-static (carries the noise); the low-pass is the **CyclinD1 mRNA reservoir** (τ=15.9h) + EZH2 IFFL (fast noise attenuated 12–36×) | **28/29** unchanged. `fig_realistic_mitogen.py` |
| **16h-cycle re-calibration attempt** | 2026-07-26 | JP: abandon the fixed 22h cycle, let it fall to Contestabile ~16h. 320-candidate parallel search re-fit to `mu=0.00069` | **25/29** (vs 28/29 at 22h) — **NOT baked**. Casualties: both HU-S+HU-G2 folds break (structural wall, worse at faster S), cKO(OFF) Shh-dependence lost, governor thresholds halve. Fundamental tension: fast cycle → low commitment threshold. Open Q: is 22h this-system's data or a default? |
| **Correct GNP period provenance** | 2026-07-26 | Nakashima GNP period is **~15.9h**, not 22–23h — no conflict with the model | doc/target fix |
| Formal validation-target spec | 2026-07-25 | **29 targets**, cell-type-separated (11 P7-GNP / 13 MB / 5 MB↔GNP ratios) + provenance-tracked. `data/validation_targets.json` + `docs/validation_targets.md` | canonical spec |

---

## Era 9 — The current default: cell-type split → two-cyclin → dynamic-CDKI bake (2026-07-28 → 07-30)

The three most recent structural bakes, each landing as the new default. This is the model as it stands. See [Model Architecture](02_model_architecture.md).

| Milestone | Commit date | What / why | Validation |
|---|---|---|---|
| **Bake the P7-GNP / SHH-MB cell-type split** | 2026-07-28 | MB's longer cycle (24–28h avg) = **same ~16h core + transient-G0 excursions** (Spencer/Cappell), not a slower engine; G0 enter/exit gated on CyclinD1 vs CDKi. GNP low-CDKi → G0~0; MB high-CDKi → G0 → 24–28h. Decouple commitment threshold from growth so governor + cKO survive | **default** |
| **Add CyclinD2 as a separate D-cyclin (two-cyclin model)** | 2026-07-29 | JP + data: CCND2 is **dominant** (83% of GNP D-cyclin pool) and **Hh-buffered** (vismo −35% vs D1 −85%). Add `Cd2` (Gli-independent floor) feeding the same CDK4/6→Rb via `w_Cd2`; total D-cyclin does not collapse under Hh-block → GNP exit ≠ D-cyclin starvation (reinforces EZH2/D1 = governor). Shh-gating survives iff `w_Cd2≤0.2` | **30/32** |
| Fix figure-regen breakages | 2026-07-29 | Two-cyclin manifest rerun surfaced stale figure scripts | — |
| Add dynamic-CDKI module (opt-in, default-neutral) | 2026-07-30 | Individual expressed CDKI species (INK4 p18/p19; CIP-KIP p21/p57) + CDK2 sequestration; reviewer structural fixes (KPC split, `CdP21`, INK4 diff); resolve `p21a` negativity + engage p21/p57 at data levels | — |
| **Bake dynamic-CDKI as the DEFAULT** | 2026-07-30 | Promote the individual-CDKI species to the working model | **30/32** — the current default |

Also in this era, MYCN was reconceived (2026-07-27): SHH-MB is **not** MYCN-amplified (Group-3 is); MYCN is **elevated as a Hh/Gli target** — an additive `MYCN_expr` term (0.54 MB) replaced the amplification multiplier, algebraically neutral. And the **EZH2 two-arms scoping** (2026-07-28): the model holds only the CyclinD1/D2-*governor* arm; the dominant *developmental* arm (EZH2 represses differentiation; cKO → early exit → pool depletion) is deliberately **not** modeled — the arms oppose at tumor initiation. These framings inform interpretation but the governor arm is what's baked. See [Reviewer & Open Questions](10_reviewer_and_open_questions.md).

---

## The through-line (what the timeline shows)

1. **Two rebuilds, not a linear march.** Era 0 (`sym`, bolt-on published oscillators) was abandoned wholesale for the Heldt-core rebuild (v42→v44) because bolted-on oscillators couldn't produce concentration-dependent S-phase and kept sticking/crashing (v31 Goldbeter, CVODE failures).
2. **Structure before parameters.** The v44 breakthroughs were *structural redesigns* (real G2, size-gated R-point, Skp2-p27 feedforward, two-step Rb, mRNA reservoir), each fixing a wall the optimizer had hit — repeatedly the **CyclinD1-level ↔ G1-length coupling**.
3. **A lot of the memory work is boundary-setting, not baking.** Strong chromatin memory, the bistable latch, the mother-G2 integrator, and the 16h re-calibration are documented **negative results** that constrain what the model *can* be — see [Things Tried & Abandoned](05_things_tried_and_abandoned.md).
4. **Validation count is a moving denominator.** 12/12 (v43) → 22/27 → 28/28 → 28/29 (formal spec) → **30/32** (cell-type + two-cyclin + CDKI). Read the count alongside the target spec of its era, detailed in [Calibration & Validation](06_calibration_and_validation.md).
5. **Reframes matter as much as fits.** CDK4/6i resistant-clone → dose-response; EZH2 protein-memory → mRNA-reservoir memory; MYCN amplification → elevated expression; Gli low-pass → CyclinD1-reservoir low-pass; EZH2 lock-in → governor/threshold. Each reframe changed the *story* without necessarily moving the score.
