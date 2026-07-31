# Mechanisms Explored

How each biological mechanism in the EZH2 – CyclinD1 – Hedgehog network is *conceived*, *modeled*, and where it currently stands. This is the mechanism-by-mechanism companion to [Model Architecture](02_model_architecture.md) and [Evolution Timeline](03_evolution_timeline.md); mechanisms that were tried and dropped live in [Things Tried & Abandoned](05_things_tried_and_abandoned.md), and the numbers each is judged against are in [Calibration & Validation](06_calibration_and_validation.md).

**Status vocabulary used throughout**
- **Baked** — engaged in the default `build_model_v44_heldt.py` model and counted in validation (currently **30/32**, two-cyclin + dynamic-CDKI default as of 2026-07-30).
- **Latent** — the machinery exists in the builder but is a no-op at string defaults (a flag or a param at 0), turned on only for a figure or an experiment.
- **Explored** — studied in a standalone reduced model / probe script; never wired into the default.

A recurring design discipline (JP): keep every new mechanism **default-neutral** (verified validation-identical *before* engaging), then engage it in `_ts_bake` — so the calibrated model is never silently perturbed. See [Repos & Reproduction](09_repos_directories_reproduction.md) for the flag/env conventions.

---

## 1. EZH2 / H3K27me3 as a CyclinD1 governor

This is the scientific heart of the project. The through-line: **EZH2's cell-cycle action runs entirely through the CyclinD1 (D-cyclin) transcript** — EZH2 does not gate the kinases directly, and the vismodegib-rescue is the clean readout of exactly this (see §1h).

### 1a. Mark topologies — repressor (setpoint) vs latch (memory)
*Doc:* `docs/prc2_buffering_memory_2026-07-19.md`; *memory:* [prc2-buffering-memory-topologies]. *Files:* `simulations/sim_prc2_buffering_memory.py`, `sim_prc2_writer_dilution_states.py`.

Two jobs a PRC2 mark on the CyclinD1 axis can do, requiring **opposite circuit topologies**:
- **Topology A — mark ⊣ CyclinD1 (negative feedback; what the model has).** Buffers the steady-state *dose-response* (full-v44 CyclinD1(mitogen) gain 0.85→0.47) but **mathematically cannot** create withdrawal memory: a pure repressor only *scales* the mitogen drive, so at mitogen=0 CyclinD1→basal regardless of the mark. Adversarial verification (independent agent, 14 tests across Hill exponents, mark timescales, basal levels): best post-withdrawal banked-division gain of any `g>0` over its own `g=0` was **−0.0001** (never positive). Deep reason: topology A closes a *negative* loop (setpoint regulation), never a latch.
- **Topology C — mark SUSTAINS proliferation (double-negative: PRC2 silences differentiation genes).** This *is* the persistence memory — a bistable latch with hysteresis; division continues through full withdrawal. But it needs the mark to have a mitogen-independent pro-proliferative effect.

**Discriminator (the clean experiment):** EZH2i *during* mitogen withdrawal should **collapse** persistence under topology C but *raise* CyclinD1 (de-repression) under topology A — opposite signs. **Status:** topology A is baked; topology C is explored-only and shown to be **mutually exclusive with calibration** (a full-v44 latch breaks GNP Hh-dependence; see [v44-latch-regime-consequences] and [v44-chromatin-memory-validation-bound]).

### 1b. The mean-field AUM repression module
*Doc:* `docs/h3k27me3_repression_model_state.md`. *Calibration:* `simulations/calibrate_h3k27_aum.py`.

`Mk ∈ [0,1]` = H3K27me3 occupancy over the ~7 kb / ~35-nucleosome bivalent *Ccnd1* domain — a mean-field reduction of the Berry–Dean–Howard A/U/M per-nucleosome model (well-mixed valid at N≈35 because PRC2 read-write is non-local, Dodd 2007). Methylation term:
`Mk_methylation = EZH2·(1−EZH2i)·(k_w·Mk + k0)·(1−Mk)·(1 − g·Cd_mRNA^p/(K_tx^p+Cd_mRNA^p))`
- **read-write** `k_w·Mk` (EED reads me3) + **de-novo nucleation** `k0` (reseeds a halved locus);
- **EZH2-scaled** whole term (so an arrested cell with falling EZH2 stops over-silencing);
- **blocked by EZH2i** (`1−EZH2i`) — all EZH2 repression flows through the mark;
- **transcription→PRC2 reciprocal arm** `(1 − g·tx-Hill)`: nascent transcription evicts PRC2 (a double-negative → positive feedback that makes the bivalent state *emerge* and made the 5.07× fold reachable, previously capped at ~4×).

Readout is **leaky**: `R(Mk) = f0 + (1−f0)/(1+(Mk/K)^n)` — a residual floor (loaded Ser5P Pol II) that dampens rather than locks out. Calibrated: `k_w=0.00160, k0=0.000258, δ=0.00070, K=0.305, n=4.15, f0=0.233, g=0.30`. **Status: baked** (`with_h3k27_dilution=True`).

### 1c. Serial me0→me1→me2→me3 methylation chain
*Memory:* [v44-serial-methylation-chain-promoted]. *Calibration:* `simulations/optimize_chain.py`.

Explicit fractions `m1_me`, `m2_me`, and `Mk (=me3`, the repressive/heritable state). PRC2 catalyses all three forward steps; FAST me0→me1 (`kme1 9.71`), me1→me2 (`kme2 8.27`), SLOW rate-limiting me2→me3 (`kme3 4.56`). Gli→Jmjd3 eraser demethylates stepwise; replication halves all modified fractions.
- **Correction (important, an earlier overclaim retracted):** the chain does **NOT** deliver steady-state division-timing sensitivity — sweeping the period 46h→12h moves settled me3 only ~4%. The ~14h me3 lag is a *de-novo establishment* transient only (Groth-lab consistent: Alabert 2015 / Reveron-Gomez 2018); at cycling steady state me3 halves and refills fast from the me1/me2 precursor pool + strong read-write (rw_frac 0.83). Its real value is biological explicitness + the correct establishment lag, not a timing switch.
- **Data constraint:** speeding me2→me3 saturates me3 in *both* GNP and MB and erases the ChIP MB<GNP contrast (0.70→0.96 as kme3 2.5→15).

**Status: baked** (`with_h3k27_chain=True`, commit 6fa1278; 28/28 at the time). Lumped single-`Mk` fallback via `H3K27_CHAIN=0`.

### 1d. The four mark mechanisms (eraser, cooperative read-write, S-dilution, H3.3 turnover)
*Memory:* [v44-four-mark-mechanisms]. *Probe:* `simulations/probe_four_mechanisms.py`.

Four mechanistically-correct additions, all validation-preserving (28/29), engaged in `_ts_bake`:
1. **Nerfed Gli→JMJD3/KDM6B eraser** `k_jmjd3_gli=0.005` (~21× nerf of the 0.108 chain default) — mitogen-gated active me3 loss, the mechanistically correct source of ChIP MB<GNP. See §1e.
2. **Cooperative read-write** `n_rw=2` — `a_rw_prc2·Mk → a_rw_prc2·Mk^n_rw`; the bistable/latch-capable EED-spreading nonlinearity (near-inert at the saturated mark ~0.98).
3. **Continuous fork-gated S-dilution** `k_dil_S=0.693` replacing the discrete at-S-onset halving, so an elongated/HU S dilutes concurrently (fixes the EZH2-banks-in-S inconsistency — see §1g).
4. **H3.3/HIRA transcription-coupled turnover** `k_h33=0.003` — replication-independent, higher where *Ccnd1* is transcribed (off in quiescent cells → mark protected).

**Net result:** ChIP direction flipped from 1.010 (flat/wrong) to **0.959 (MB<GNP)** while the CyclinD1 fold barely moved. Effective mark decay is now three terms: **δ_eff = del_mk + k_jmjd3_gli·Gli1 + ln2/Tc**. Reaching the ~0.5 ChIP magnitude needs a less-nerfed eraser or the de-saturation re-opt, which is walled by [v44-chromatin-memory-validation-bound]. **Status: baked.**

### 1e. Gli–Jmjd3 eraser reconception (the mark's *direction*)
*Memory:* [v44-gli-jmjd3-eraser-reconception], [v44-h3k27me3-persistence-tension]. JP's ChIP: MB has ~half the GNP H3K27me3 on CyclinD1. Reconceived as a **writer(EZH2, cycle-coupled) vs eraser(Gli-Jmjd3, mitogen-coupled) race** — high-Gli MB strips the mark. `k_jmjd3_gli` sweep 0→0.1 drives ChIP Mk MB/GNP 0.97→0.81. **Tension:** the eraser is the *only* MB<GNP ChIP source, but strengthening it shifts the CyclinD1 fold and shortens mark persistence; no continuous regime simultaneously gives long-mark + ChIP direction + fold. **Status:** baked but deliberately weak (`w_ezdir` default 0).

### 1f. Replicative dilution (the endogenous-period phenotype)
Each S-phase halving raises the readout `R`; a faster cycle dilutes the mark more, so cycle period feeds back on repression — the basis of the T_cc phenotype and the EZH2i × CDK4/6i synergy prediction (CDK4/6i lengthens T_cc → fewer dilutions → mark restored → differentiation). **Status: baked** (via §1b/§1d).

### 1g. Cdh1-gated EZH2 accumulation in S (the concentration convention)
*Doc:* `docs/cdh1_ezh2_sphase_boost_2026-07-19.md`; *memory:* [v44-cdh1-ezh2-sphase-boost] (reverses [v44-ezh2-concentration-promoted]).

JP data: elongating S *should* raise EZH2 per cell; the model showed no effect (HU EZH2-in-S boost stuck at 1.009 vs target 1.31). Root cause: the concentration convention + a single turnover rate. **Fix (baked):** EZH2 is a bona fide APC/C-Cdh1 substrate — (i) **dilution convention** `with_ezh2_conc=False` (EZH2/2 at mitosis, so arrested cells bank it); (ii) low baseline `kDeEZ` 0.00117→0.000281; (iii) **Cdh1(C1)-gated degradation** `EZH2_deg: ... (kDeEZ + kDeEZ_C1·C1)·EZH2` — clears EZH2 in G0 (C1 active) yet accumulates in S/arrest (C1 off). Result: HU boost 1.009→**1.24**, EZH2 G2/G0 1.04→**1.86**, G0/cycling 0.44→**0.62**; 27/28. This is the independent support for the EZH2⊣CyclinD1 axis (longer S → more EZH2 → lower CyclinD1, matching the lower-CyclinD1-under-HU data).

### 1h. Kinetics vs local repression; the poised/hair-trigger reframe
*Memory:* [h3k27me3-kinetics-vs-local-repression] (JP 2026-07-30, foundational, conceptual only).

The mark **spreads over a large domain but its transcriptional-repression impact is very local (proximal promoter)**; the broad domain mostly enhances the *kinetics* of transcriptional control, not the steady-state magnitude. This resolves "marked but expressed" — a gene can be heavily marked yet transcribed when strong activators + elongation machinery are present. Functionally the mark **raises the activation threshold** (responsive only to high activator) and makes the gene **turn off fast** once activators drop — classic bivalency/poising. Maps directly onto the governor finding (EZH2 raises the mitogen threshold for CyclinD1 without decoupling it). **Modeling implication (future, NOT built):** decompose `Mk` into a small local component (steady-state level) + a broad component (rate/time-constant). Status: conceptual capture; JP said "forget it for now, revert."

### 1i. Speed-limiter / rate-sensor framework
*Doc:* `docs/speed_limiter_framework_2026-07-21.md`; *memory:* [ezh2-speed-limiter-framework]. *Files:* `simulations/sim_speed_limiter.py`, `fig_speed_limiter.py`.

JP's synthesis naming what the pieces collectively do. H3K27me3 is a **graded writing-vs-dilution rate-sensor**, not a bistable latch; the bistable element is the **division decision** (R-point), and CyclinD1-being-a-target makes the sensor a **negative-feedback governor**. The writer (EZH2) is proliferation-coupled and **saturating** → there is a *max sustainable division rate*: dividing too fast lets `ln2/Tc` outrun the writer, the mark collapses below θ_diff, and differentiation genes de-repress. GNP sits at the sustainable ceiling (pushed hardest); MB overrides the governor via constitutive Gli/MYCN into the "too-fast / low-mark" regime — consistent with ChIP MB<GNP (MB escapes differentiation by oncogenic means, not the governor). **Status: explored** (conceptual; the one explicit missing piece is an H3K27me3-marked *differentiation-sentinel* locus — only CyclinD1 is marked today).

### 1j. The two EZH2 arms and cKO scoping
*Memory:* [ezh2-two-arms-cko-scoping] (JP 2026-07-28).

EZH2 has **two separable arms**:
- **Arm 1 — developmental, DOMINANT, NOT modeled.** EZH2 represses the *differentiation* program → cKO causes premature differentiation → cell-cycle exit → pool depletion → **Ptch⁺/⁻;EZH2-cKO does NOT increase MB**. EZH2 loss does *not* drive tumors (this corrects an earlier EZH2-as-tumor-suppressor-via-CyclinD1 over-statement).
- **Arm 2 — the CyclinD1(/2) governor, the model's scope.** EZH2 ⊣ CyclinD1 raises the mitogen threshold without decoupling division from Shh; EZH2i lowers the Shh bar (the vismodegib-rescue readout).

The arms **oppose at tumor initiation** — so cKO genetics and acute EZH2i pharmacology can look like opposite functions and both be correct. **cKO mechanism = reduced H3K27me3 *mark*, not reduced PRC2 *occupancy*** (distinct from the drug EZH2i, which blocks writing while the complex stays bound — [v44-ezh2i-catalytic-mechanism]). Model D-cyclin should be read as **CyclinD1/2**.

### 1k. Governor re-measurement (the load-bearing numbers)
*Memory:* [ezh2-governor-reservoir-remeasure]. *File:* `simulations/fig_governor_noise_filter.py`.

Re-measured in the reservoir+IFFL model, then re-baked with `k_Cd_tx_basal ×0.5` (0.0005615→0.0002808, JP-approved; validation-neutral 28/29). Post-bake the Shh-dependence is **structural** (from the *low basal*, not from EZH2). Mitogen thresholds: **ON 0.18 / EZH2i 0.06 / OFF 0.01** → EZH2 repression raises the Shh threshold **~15×** and buffers the dose-response, but all three arrest at Shh=0 (matches the cKO: transcript up, still Shh-dependent, differentiates). **Status: baked.** (This supersedes the old ~2.6× figure in [ezh2-repression-mitogen-requirement].)

---

## 2. CyclinD1 transcript RESERVOIR = the Hh memory + EZH2 incoherent-feedforward

*Doc:* `docs/cyclind1_hh_memory_2026-07-21.md`; *memory:* [v44-cyclind1-hh-memory-carrier]. *File:* `simulations/fig_cyclind1_memory.py`. *Source:* Ho, Tsai & Stearns 2020 (Curr Biol; JP acknowledged).

**The idea:** the Hh proliferation decision sits in the *previous* cell cycle — cells inhibited up to ~11h before mitosis still divide. **CyclinD1 carries that memory, but the memory is the cumulative *transcript* reservoir, not the protein.** The protein is labile (fast readout); the slow mRNA pool integrates Hh over the cycle, carries through mitosis, and is inherited. Gli raises the pool; EZH2 lowers it.

**EZH2 incoherent feedforward (IFFL):** Gli drives *both* the transcript (up) and EZH2 (up), and EZH2 represses the transcript (down). When mitogen falls, direct Gli drive drops fast but slow EZH2 (τ≈14h) lifts its repression with a lag → the transcript **stays up longer**, buffering the dose-response (~1.7× higher than a frozen-repression control).

**How modeled / the fix (baked):** the mRNA had `k_Cd_mRNA_deg = 0.8` (~1-min half-life → slaved to mitogen). Changed to `0.001` (τ≈15.9h reservoir) with the three transcription rates (`k_Cd_tx_basal`, `k_Cd_tx_Gli_max`, `k_Cd_tx_MYCN`) scaled ×1/800 to hold steady levels exactly. Protein stays fast (`k_Cd_deg=1.0`). **Validation-neutral (28/29)** — only the *dynamics* become a memory; after Hh-off the reservoir coasts ~1 cycle. This reservoir is *also the low-pass filter* for §6 (measured first-order low-pass, τ≈16h). **Status: baked.** *Caveat later raised:* [ccnd1-halflife-data-vs-reservoir] — measured CCND1 t½ is ~1.5–2.5h (RNAdecayCafe), ~5–7× faster than the baked 11h reservoir; a fast transcript still coasts ~1 division, so the ~11h memory is arguably committed cell-cycle state (Rb-E2F), not the transcript. Fast-transcript is validation-neutral (28/29); decision pending.

---

## 3. Two-cyclin model: CyclinD1 + CyclinD2

*Memory:* [cyclind1-d2-two-cyclin-model] (JP 2026-07-28). Data (Chahin RNA-seq): **CCND2 is dominant** (83% of the GNP D-cyclin pool) and **Hh-buffered** (drops only ~35% under vismodegib vs CCND1's ~85%). Combined-pool folds: MB/GNP 3.25, GNP+vismo 0.58, MB+vismo 0.42.

**Key implication:** total D-cyclin does *not* collapse under Hh-block (the D2 floor persists) → GNP exit is **not** D-cyclin starvation — reinforcing EZH2/CyclinD1 as a threshold governor, not the on/off switch.

**How modeled:** `species Cd2` with `CycD2_synthesis: => Cd2; (k_Cd2_bas + Cd2_expr + k_Cd2_Gli·Gli/(K_Cd2_Gli+Gli))`, feeding the **same** CDK4/6→Rb drive via `w_Cd2·Cd2`; **no EZH2 repression** on D2. Shh-gating survives the D2 floor iff D2's Gli-independent drive stays below the R-point → constraint `w_Cd2 ≤ 0.2`. **Baked (2026-07-28), the two-cyclin model is now the default, 30/32**: `k_Cd2_bas=2.0, k_Cd2_Gli=2.49, K_Cd2_Gli=0.3, w_Cd2=0.079`; D2 folds MB/GNP 2.39, GNP+HHi 0.66, MB+HHi 0.77. The recal **fixed the persistent CyclinD1 GNP+HHi miss** (0.055→0.100) — with D2 carrying the buffered floor, `k_Cd_tx_basal` could be raised so D1 stays sharply Hh-responsive.

**Is D2 H3K27me3-marked?** JP: it should be. Tested "equally marked as D1" (opt-in `d2_marked`, default OFF): **data-incompatible (27/32) — the vismo direction flips** (D2 goes *up* under vismodegib because its large Gli-independent floor de-represses). A 200-candidate recal confirmed equal-marking is **over-constrained** (best 29/32, only 2/3 D2 folds; three-way trade-off between MB/GNP, the vismo drops, and de-repression). Reconciliation (§1h): D2 can be heavily *marked* while its transcription stays largely mark-*independent*. **Status:** two-cyclin baked; D2-marking latent (flag off), recommended future refinement = a weak/partial D2 mark.

---

## 4. MYCN: elevated expression, not amplification

*Memory:* [mycn-elevated-expression-not-amplified] (JP 2026-07-27). SHH-MB is **not** MYCN-amplified (Group-3 is); MYCN is high because it is a Hh/Gli target. Data verifies the targets: MB/GNP 2.86, GNP+vismo 0.80 (−20%), MB+vismo 0.89 (−11%); escalation wt-GNP 3077 → Ptch⁺/⁻ 4251 → MB 8808.

**Key:** the acute vismodegib response *shrinks* as baseline Hh rises → high-Hh contexts are developmentally **locked / vismo-resistant** (MB MYCN ~86–89% acutely Gli-independent). **How modeled (baked, algebraically neutral 26/29):** `MYCN_synthesis: => MYCN; (k_MYCN_synth_basal + MYCN_expr) + k_MYCN_synth_Gli·Gli/(K_Gli_MYCN+Gli)` — an **additive `MYCN_expr`** term (MB = `0.3·(2.8−1)=0.54`, Gli-set) replaces the old `MYCN_amplification` multiplier. **Abandoned:** the acute *autoregulation latch* (JP wanted MB-only) — the buffer and the GNP/MB discrimination are the same knob `K_auto` with no window between them, plus numerical fragility; the lock is developmental, not an acute switch (`mycn_autoreg` flag deprecated/inert). **Status: baked** as the 5th cell-type knob.

---

## 5. Dynamic-CDKI module (INK4, CIP/KIP, CdP21 buffer, KPC)

*Doc:* `docs/cdki_module_design.md`; *memory:* [dynamic-cdki-module], [cdki-baked-g0-g1-behavior]. Two expert-review rounds applied. Replaces the lumped/mislabeled static `p16`/`p18` params + a single `P21` species standing in for all CIP/KIP.

**RNA-seq scope (P7 GNP → SHH-MB):** modeled = **p18 (Cdkn2c, dominant INK4)**, p19 (Cdkn2d), **p27 (Cdkn1b, dominant CIP/KIP ~85%)**, p21 (Cdkn1a); dropped = p16/p15 (silent in GNP); de-emphasized = p57 (Cdkn1c, a Sox2⁺ quiescent CKI, carried OFF `kSyp57a=0`).

**Structure (10 species added, behind `with_cdki_species`):**
- **INK4 → CDK4/6:** p18 + p19 protein species (p18 stable `kDe_p18=0.003`, p19 short `kDe_p19=0.03`); CDK4/6 competitive brake `1 + p18_prot + p19_prot + w_p27·(p21a+p57a)`.
- **CIP/KIP → CDK2:** p27 + p21 + p57, each with own CyclinE/A-CDK2 sequestration complexes; **p21 carries the PCNA/Rc arm + CRL4^Cdt2** (moved off p27 per reviewer — p27 has no PIP box; later reverted to the `aRc` proxy because the shared iPcna/iRc species couldn't host a near-zero p21 pool → negativity).
- **p27 degradation:** basal + **KPC (constant, mitogen-insensitive, FREE-pool only, `kDeKPC≈0.006`)** + SCF^Skp2. The old CDK4/6-activity-gated clearance is **removed**.
- **CdP21 (the new commitment lever):** p27 buffered on CyclinD1-CDK4/6, `Cd + P21 ⇌ CdP21` (Cd consumed → *saturating*; buffered Cd returned to the Rb-drive numerator = option A; INK4 promotes release). **Not reset at division.** This is the **Shh-gated commitment**: high Shh → high Cd → buffers p27 → commit; Shh=0 → low Cd → free p27 rises → arrest.

**Stage-3 calibration succeeded — 30/32** (best: `kSeqCd=2.132, kRelCd=0.0165, kDeKPC=0.0202, w_ink4=3.62, k_mu_cki=0.362, kSyP21=0.00116`), with GNP-SHH arrest (0 div), cKO@0 = 0.00, and the **GNP 17.6h / MB 26.1h dichotomy**. Fails only the structural HU pair ([v44-hu-s-g2-structural-limit]).

**Two genuine findings from calibration:** (1) **p21/p57 cannot be engaged at data-transcript levels** — at their abundances (~0.49/0.24) they arrest the cell even off the CDK4/6 brake and with p27 near-zero (≈560 sims → 15/32); this is the reviewer's predicted degeneracy *plus* the biological point that **mRNA abundance overstates the active nuclear CDK-inhibitor pool** for these post-translationally-regulated proteins. Resolution: p21/p57 stay expressed-but-inert; functional CDK regulation is **p27 + INK4 + CdP21**. (2) **The baked commitment is SHARP** — no discrete p27-high/pRb-low transient G0 (G0frac=0%); MB's longer cycle is a **longer G1** (21.3h vs GNP 10.8h) set by high INK4, *not* a discrete pause (a change from the prior birth-p27 + slow-clearance model). See §8.

**Status: baked as default (commit ec1a141, JP-approved).** `p57` de-emphasized/inert. Legacy lumped CDKI via `CDKI_SPECIES=0`. Note: CDKI-on validate is slow/stiff (>200s).

---

## 6. Mitogen input (noisy Shh, cilium, Ptch1 adaptation, reservoir low-pass / IFFL)

*Docs:* `docs/realistic_mitogen_input_2026-07-21.md`, `docs/hh-signalling-timing.md`; *memory:* [v44-realistic-mitogen-input]. *Files:* `simulations/sim_realistic_mitogen.py`, `fig_realistic_mitogen.py`. **All new terms OFF by default → validation unchanged.**

Input model `Shh_effective(t) = L(t) · G_cilium(cell-cycle) · A(t)`:
1. **Ligand L(t)** — paracrine, Purkinje-derived, **tonic** (mean 0.45) + a **slow OU drift** (τ≈40h, passband) + **fast OU noise** (τ≈0.5h, stopband); CV≈0.22. *Not* a fast pulse train (that is ERK). Amplitudes are a stated assumption (GNP ligand fluctuation timescale unmeasured); the point is the timescale separation.
2. **Cilium duty cycle** — cilia resorb pre-mitotically, so reception is phase-gated ON in G1/S, OFF through G2/M: `cilium_gate := 1 − w_cilium·Cb^n/(K_cil^n+Cb^n)` (CyclinB-gated; ON-fraction ≈0.86). `w_cilium=0` default.
3. **Adaptation A(t)** — the ~20h temporal adaptation (Ptch1 feedback + Gli2-down): `Ptch1_slow` term (slow, Gli²-driven), tuned **weak** (`k_ptch1_slow_on=0.008` → Gli drops ~35% over ~20h) because GNPs proliferate for days under sustained Shh. **A scanned free parameter** (the central GNP unknown).

**Key correction (2026-07-24):** measuring the transfer function shows **Gli is quasi-static** — the ligand→Gli gain is flat (~0.15) at every period 0.5–160h; Gli *tracks* the ligand and *carries* the noise (Gli protein τ≈40 min–2h ≪ cycle). Gli does **not** low-pass. The actual low-pass filter for the proliferation decision is **downstream: the CyclinD1 mRNA reservoir** (τ≈16h, §2) + the EZH2 IFFL — a clean first-order low-pass that attenuates fast noise (period <3h) **~12×** in amplitude, up to **~36×** in Bode gain, while slow drift passes. This ties the noisy-input story directly to §2. **Status: explored/latent** (figure-only; calibrated model untouched).

---

## 7. Skp2-p27 commitment / two-step Rb R-point / size (Rb-dilution) gate

The bistable commitment switch — "divide or not" (the R-point of §1i).

- **Skp2-p27 & mitogen-dose Skp2.** *Memory:* [v44-skp2-mitogen-dose]. Skp2 is higher in MB (2.3×, JP) but the model had it backwards; fixed via a mitogen(CyclinD1)-dose term `kSySkp2_Cd` → Skp2 MB/GNP 2.31. **Baked (28/29).** Skp2 degrades complexed p27, releasing active CDK2 — the switch that liberates CDK2.
- **Two-step Rb / R-point.** *Memory:* [v44-transient-g0-reframe-and-cdk46i-escape]. A CDK4/6-activity-gated pRb promotion with pRb as the G0 marker (branch, 27/28). Superseded for commitment by the **CdP21 buffer** (§5) — the current Shh-gated R-point.
- **Size gate = Rb-dilution.** *Memory:* [v44-sizegate-is-rb-dilution]. The v44 mass gate is equivalent to the Zatulovskiy Rb-dilution R-point; v44 already has the bistable Skp2-p27 toggle, so porting v45's explicit CDK2-p27 toggle is cosmetic. **Status: baked** (mass/growth gate).
- **CDK4/6i as a dose-response.** *Doc:* `docs/transient_g0_synthesis.md` §7; `simulations/sim_cdk46i_dose_response.py`. Modeling palbo as graded residual CDK4/6 activity `kPhRbCd = KP·resid` over a heterogeneous population reproduces JP's two doses (16–18% pRb⁺ at 1µM → ~2% at 5µM, empirical K≈1.5µM) with **no CDK2-bypass and no senescence compartment**. Falsifiable prediction: the transition is intrinsically **steep** (the bistable R-point makes each cell near-binary), so CyclinD1/CKI spread cannot soften it — intermediate 2–3µM doses should already be ~saturated. Retires the earlier resistant-clone framing.

The related **quiescence bifurcation** (CyclinD1/p27 ratio, Spencer/Cappell; GNP/MB p27-dominant) is in [v44-quiescence-bifurcation-direction]; the **mother-G2 p27 integrator negative result** (does not reproduce the Overton bifurcation) is in [v44-mother-g2-negative-result].

---

## 8. Transient G0, graded G1-lengthening, and mitogen withdrawal

*Docs:* `docs/transient_g0_synthesis.md`, `docs/g1_lengthening_buffering_2026-07-20.md`; *memory:* [transient-g0-literature-constraints], [mb-celltype-transient-g0-parameterization], [gnp-g1-lengthening-buffering], [withdrawal-graded-g1-lengthening], [matched-withdrawal-rate-dependent], [cdki-baked-g0-g1-behavior].

**The biology (firmest first):**
- **GNP outer EGL transient G0 ≈ nil** — obligate cycling with a graded p27 ramp toward *one-way* exit (the Spencer CDK2-low framework is p21-driven cell-line biology; the GNP Cip/Kip is p27).
- **MB G0 is real but p27 does not measure it** — p27 IHC reports total protein; the reversible-G0 remainder (~20% *within the proliferation-competent pool*) is an informed guess, not a measurement. The permanent differentiated exit is deliberately **out of scope** (a separate compartment).
- **CDK4/6i inducibly fills a reversible G0** — cells accumulate p27, go pRb-hypophospho, and **re-enter** (reversible quiescence, not senescence). The functional cycling/arrest switch is **pRb phosphorylation state** (Ser807/811, clone D20B12), not p27 level.

**Transient-G0 in the model, over time:**
- The **prior** model produced a discrete transient G0 via birth-p27 + slow CDK4/6-gated p27 clearance (GNP 32% / MB 19% by p27 marker — [v44-transient-g0-finding]; a deterministic multi-periodic G0 near threshold — [v44-deterministic-multiperiodic-g0]; arrest is **monostable/reversible**, no true bistability — [v44-g0-arrest-monostable]).
- The **baked dynamic-CDKI** model (§5) **removed the discrete G0**: the CdP21 buffer + fast constant KPC keep free p27 low, so commitment is sharp (G0frac=0%). **MB's longer cycle is now a longer G1 set by high INK4**, not a pause. Recovering a discrete transient-G0 *subpopulation* would need the population layer (birth-p27 heterogeneity) + weaker CdP21 buffering — **not** a clearance change.

**Graded G1-lengthening (JP's expectation for GNP — lengthen the cycle without a transient G0):**
- The lever is the latent `mu_eff` growth-slowdown (G1 is growth-limited). Tuning `mu_min_frac 0.4→0.1`, `K_g1len 0.6→0.8` makes the steady period lengthen smoothly 24→84h as Shh falls, then arrest — **validation-compatible 28/29** ([withdrawal-graded-g1-lengthening]). In the baked model it is latent because `K_g1len=0.6` sits below the operating Cd; **raising `K_g1len`→1.8–2.5 recovers clean graded G1** (GNP 26.7h→14.0h as SHH 0.6→1.0, [cdki-baked-g0-g1-behavior]).
- **Important negative result:** JP's *EZH2-reduction* driver is NOT supported — the model banks EZH2 (slow τ~14h, tracks the cycle not the mitogen), so it stays flat during withdrawal; forcing a fast `kDeEZ` breaks validation (23/29, pinned by HU-S/G0). So the phenotype is a **growth-slowdown, not EZH2 reduction**.

**Withdrawal is rate-dependent** ([matched-withdrawal-rate-dependent], [v44-mitogen-withdrawal-depth-gated]): matched to the same Ccnd1, a **fast** decline (<τ~16h) is memory-dominated (all arms arrest identically), while a **slow** decline (120h) is threshold-dominated (the de-repressed OFF cell divides far longer). The EZH2-feedback effect on post-withdrawal divisions is maximal at *partial* withdrawal, ~0 at full; sudden > gradual.

**Complementary piece:** graded G1-lengthening supplies the *buffering* (slow EZH2); the ChIP-seq *direction* (mark maximal at max division, declining with differentiation) still needs the active differentiation erasure from §1a topology C — they are complementary, not alternatives.

**Status:** the sharp-commitment behavior is **baked**; discrete transient-G0 and graded G1-lengthening are **latent** (recoverable via known levers — population birth-p27 / weaker CdP21 for G0; `K_g1len` for graded G1), neither engaged by default.

---

### Cross-references
Foundations of these mechanisms in the published BioModels they were built from: [Origins & Foundations](01_origins_and_foundations.md). Where each sits in the 3-tier structure (fast slaved input → oscillator core → slow epigenetic/growth memory, [v44-structural-map-3tier]): [Model Architecture](02_model_architecture.md). The order they were added and re-baked: [Evolution Timeline](03_evolution_timeline.md). The data each is calibrated against: [Data & Evidence](07_data_and_evidence.md) and [Calibration & Validation](06_calibration_and_validation.md). The figures that visualize each: [Figures Catalog](08_figures_catalog.md).
