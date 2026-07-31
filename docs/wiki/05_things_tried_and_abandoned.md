# Things Tried & Abandoned (dead-ends and negative results)

> The crucial "what did **not** work and **why**" catalog. Over ~4 months JP tried a great many mechanisms; most of the interesting ones are dead-ends, and the *reasons* they died are load-bearing scientific results in their own right — they define the boundaries of what the data will allow. This page is the map of those boundaries. See also [Mechanisms Explored](04_mechanisms_explored.md), [Calibration & Validation](06_calibration_and_validation.md), and [Reviewer & Open Questions](10_reviewer_and_open_questions.md).

**How to read this page.** Each entry gives: (1) what was tried, (2) the result / score, (3) *why* it was abandoned. Almost nothing here was baked into the default model — that's the point. The recurring theme is a small number of **hard data walls** (the measured CyclinD1 MB/GNP ≈ 5.07 fold, the soft vismodegib de-repression ≈ 0.144, the HU S-vs-G2 either-or, the GNP Shh-dependence / cKO constraint) that repeatedly kill otherwise-attractive mechanisms.

The negative results cluster into six families:
- **Bistable latches / strong epigenetic memory** — excluded by the fold + vismo walls.
- **Cell-cycle recalibration** — the 16h attempt and its knock-on losses.
- **Transient-G0 mechanisms** — repeatedly structurally absent in the deterministic single-lineage model.
- **CDKI biology at data abundances** — p21/p57 arrest the cell; mRNA ≠ active inhibitor.
- **Withdrawal-memory dynamics** — the governor is overdamped; no oscillation.
- **Conceptual corrections** — framings JP or the model held that turned out to be backwards.

---

## 1. Bistable latches and strong H3K27me3 memory

### 1.1 Strong self-reinforcing chromatin memory (the bistable latch) — EXCLUDED by two data walls
**Tried.** Re-bake v44 into a memory-competent chromatin regime: strong read-write autocatalysis (`k_w_mk` ≳ 6–16× the baked 0.0045) and/or a durable mark (`del_mk` t½ ≳ 12h, ≤0.001), so the H3K27me3 mark self-latches and sustains a withdrawal memory ([`v44-chromatin-memory-validation-bound`](07_data_and_evidence.md); `sim_mark_memory_chromatin_sweep.py`, `sim_mark_memory_cooptimize.py`, ~450 candidates).

**Result.** Every memory-competent regime **drops CyclinD1 MB+HHi/MB below the pass floor** (e.g. `k_w_mk`=0.025 + `del_mk`=0.001 → 0.077; 23–24/27). More decisively, a constrained search co-optimizing against the **measured MB/GNP = 5.07 fold** found **zero feasible points with memory > 0**: every point with MB/GNP ∈ [4.75, 5.45] and the other targets passing had the mark fully OFF in withdrawal. No Pareto trade-off exists — memory > 0 ⟹ MB/GNP infeasible.

**Why abandoned.** Two independent data walls exclude strong chromatin memory:
1. **Soft:** vismodegib de-repression MB+HHi/MB = 0.144 (~14% residual CyclinD1; a strong memory predicts ~8%).
2. **Hard:** the mitogen-dose fold MB/GNP = 5.07. A memory-competent mark **saturates in both GNP and MB cycling** → loses mitogen-dose sensitivity → cannot repress MB more than GNP → the ratio snaps back to the cascade's raw ~7. The dose sensitivity that pulls the ratio down to 5.07 is *exactly* the graded/weak regime that cannot hold memory. The read-write term `k_w_mk·Mk` locks the arrested MB+HHi mark (residual EZH2 still writes), and this is not fixable by curve shape (`K_mk`, `f0_mk`) or by transcription-eviction `g_prc2` (which pushes the ratio the wrong way). Same physics as the weak-feedback pin from [`v45-feedback-paramspace`] and the monostable G0 arrest below. **Decision: baseline (memory OFF) is the working model; trial re-bake reverted.** A *weak/transient* memory (charged during proliferation, discharges in prolonged arrest) is compatible — but mild by necessity, not tunable up.

### 1.2 G0 arrest as a bistable "EZH2-maintained lock-in" — RETIRED (arrest is monostable)
**Tried.** The framing that sustained high EZH2 holds the H3K27me3 mark → a stable, self-maintained epigenetic lock in arrested G0 cells (`fig_v44_ezh2i_rescue_kinetics`, early `fig_v44_g0_hysteresis`).

**Result / why abandoned.** A definitive two-initial-condition test (`fig_v44_g0_bistability_confirm`): from high-mark or low-mark starts the cell converges to the **same** steady state at every drive → **monostable, reversible**. The earlier `fig_v44_g0_hysteresis` sweep false-positived bistability (a slow-transient artifact) and was superseded. In arrest the writer EZH2 is E2f-gated so it collapses to basal; the mark settles ~0.31 (*below* the cycling mean ~0.39). EZH2 is cell-cycle-gated (Skp2/E2f), so it is **not** MYCN-sustained in a non-cycling cell. The honest paper claim is *reversible/metastable ("sticky")* G0, not an irreversible switch. (The EZH2i-rescue result itself is unchanged — it re-tips a borderline CyclinD1, not an epigenetic lock.)

### 1.3 Long-persistence mark via killing the Gli eraser — blocked by the ChIP direction
**Tried.** JP flagged that the model's mark decays too fast (~0.2h in MB; the Gli→Jmjd3/Kdm6b eraser dominates the me3→me2 rate). Reduce `k_jmjd3_gli` 0.221→0 and `del_mk`→0.0006 to lengthen the mark toward a ~1-cycle memory ([`v44-h3k27me3-persistence-tension`]; `h3k27me3_persistence_exploration.md`).

**Result.** JP's lever *works in isolation* — it lengthens the mark to **18.4h** (~90×).

**Why abandoned (core tension).** The **same eraser term is the ONLY source of the measured MB < GNP ChIP** (MB has ~half the GNP H3K27me3 on CyclinD1). Zeroing it drives MB/GNP mark 0.56 → **1.02** (mark saturates equally in both). Transcription-eviction (`g_prc2`) + replication dilution alone do not separate MB from GNP. **No continuous-model regime gives all three** (long mark + MB<GNP≈0.5 + 5.07 fold). Additionally, "eraser-off gives long functional persistence" was itself over-stated: measured apples-to-apples (EZH2i + cycling), eraser-off persistence is only **~2.8h** (MB), not ~18h — because the mark is drained by **cycling-driven replication dilution**, not the eraser. The ~18h was writing-off *without* cycling. So removing the eraser cannot push functional persistence past ~3h; a true >10h memory needs the bistable latch (data-excluded, §1.1).

### 1.4 The PRC2/transcription paradox → an allele-level stochastic switch — designed, NOT built
**Tried.** JP's paradox: transcription antagonizes PRC2, yet high PRC2 + high expression coexist and peak at max division, and a 0.2h mark "makes no biological sense." A single deterministic locus with the transcription↔PRC2 double-negative is a **bistable switch with no graded intermediate** (weak eviction → MB/GNP≈0.92; strong → all-or-nothing, CyclinD1 jumps ~9× losing dose-grading) — which is precisely *why* the phenomenological eraser existed (`prc2_transcription_paradox.md`, `fig_prc2_transcription_paradox.py`).

**Proposed resolution (validates JP's punctate-transcription idea).** Move the switch to the **single-allele** level (stable latched mark) and drive OFF-flipping by **punctate/stochastic transcription** (rate ∝ mitogen); the *population average* then grades smoothly → MB/GNP ≈ 0.40 + dose-graded CyclinD1, with every allele's mark stable and no fast eraser.

**Why not adopted.** This is a **stochastic module** (Gillespie/tau-leap on allele state) — a structural restructuring, not a parameter re-fit. Default left unchanged; documented as the recommended future direction, not implemented. Later reframed conceptually by [`h3k27me3-kinetics-vs-local-repression`] (broad mark = kinetics, local = repression).

### 1.5 The serial me1/me2/me3 chain as a division-timing switch — OVERCLAIM corrected (but chain kept)
**Tried / result.** The explicit me1→me2→me3 methylation chain was promoted to default expecting it to deliver steady-state division-timing sensitivity. Independent verification (wf_e5ca73e6) showed **it does not**: sweeping period 46h→12h moves settled me3 only ~4%; at cycling steady state me3 merely halves and refills fast from the high me1/me2 pool + read-write autocatalysis (rw_frac 0.83). The ~14h me3 lag is a **de-novo establishment transient only**. The chain was **kept** as default — but for the right reasons (biological explicitness, correct establishment lag, faster me2→me3 JP wanted), *not* as a timing switch. Speeding me2→me3 further erases the ChIP MB<GNP contrast (0.70→0.96), so it's data-bounded.

---

## 2. Cell-cycle recalibration to 16h — reached 28/29 but lost the behaviors

**Tried.** JP: "abandon the 22h fixed cycle, let it fall to the Contestabile ~16.25h value, retest." The 22h period is a **growth-timed constant** (Tc ≈ ln2/μ, μ=0.0005; scales exactly as 1/μ, ~invariant to mitogen). Set μ=0.00069 → Tc 16h and re-optimize (parallel 320-candidate search; `recal_16h.py`) ([`recalibration-16h-cycle-attempt`]).

**Result (global-μ version).** Best set preserving ON-Shh-dependence = **25/29** (vs 28/29 at 22h). Casualties:
- **HU-S (0.73) + HU-G2 (0.37) folds BOTH break** (the structural wall, §4.1; worse at faster S-phase).
- MB G2+M too short (1.0h vs 2.5h); CyclinD1 GNP+HHi (0.047) basal-cut trade-off.

**Then (the fundamental tension).** The faster cycle **lowers the effective commitment CyclinD1 threshold** (cell reaches M_size sooner, commits with less CyclinD1) → **cKO (OFF) Shh-dependence LOST** (OFF divides at SHH=0, *not* basal-fixable — residual Gli drives un-repressed CyclinD1 over the lower threshold) → **governor thresholds HALVE** (ON 0.18→0.09). This is the load-bearing negative result: **fast cycle ↔ low commitment threshold ↔ (cKO Shh-independence + broken HU folds)**.

**Partial rescues, and why the 16h path still wasn't baked.**
- **Cell-type split (`k_mu_cki`, CKI-coupled cycle length):** GNP 16h / MB 22–23h via INK4-slowed MB growth → HU folds recover (they're MB), reaching **26–28/29**. But at the validation-optimal 28/29 the *behaviors* re-weaken (cKO lost, governor ON 0.18→0.07, graded-lengthening weaker) — a **trade-off surface**: validation-optimal (no cKO, weak governor) vs behavior-preserving (26–27/29).
- **`decouple_commit` (remove the size-gate from p27-clearance):** did **not** fix cKO at 16h — the growth→threshold coupling is *also* in the mass-scaled CyclinE/A synthesis, not only the size-gate. Needed an additional threshold re-tune (kPhRbCd 0.35→0.14).
- **Net:** the cell-type split *was* eventually baked (μ=0.000769, k_mu_cki, decouple_commit, GNP core ~17h / MB ~23h; 26/29) — but the pure global-16h recalibration is a documented dead-end, and JP's original "22h is just a default, drop it" premise was answered: 22h is a *growth-timed constant*, and forcing 16h costs the Shh-gating behaviors.

---

## 3. Transient G0 mechanisms — repeatedly structurally absent

### 3.1 MB transient-G0 lengthening at a 16h core — IMPOSSIBLE without a population layer
**Tried.** JP's architecture: MB's 24–28h average is the *same* 16h core cycle plus frequent **transient G0** excursions (Spencer/Cappell), gated on CyclinD1-vs-CDKi — *not* a slower engine. Realize it deterministically at 16h ([`mb-celltype-transient-g0-parameterization`], `characterize_decouple.py`, `pop_celltype_16h.py`).

**Result / why abandoned.** **MB transient G0 ≈ 0%** in every deterministic single-lineage route. Root cause (fundamental): **MB CyclinD1 is ~5× GNP (a measured target) and saturates the Rb Michaelis drive** — at Cd ≫ K_CdRb·brake the CDKi brake is irrelevant, so MB actually commits *more* easily than GNP (cdk_d 0.89 vs 0.79), the opposite of what "more G0" needs. Raising MB birth-p27 to 15× does nothing (huge CyclinD1 clears it every cycle). So the 24–28h-via-transient-G0 picture is intrinsically a **population/stochastic** phenomenon (a low-CyclinD1/high-CDKi subpopulation dwells in G0), not achievable in one parameter set. JP's chosen resolution reverted to **MB core genuinely longer** (k_mu_cki, driven by MB's measured INK4) — contradicting the "transient-G0 is most of it" hypothesis. An arithmetic tension is flagged: flow G0 ~20–25% + 16h core → avg ~20–21h, *not* 24–28h (which would need ~35–40% G0).

### 3.2 Discrete transient G0 via p27-clearance tuning — structurally absent under CdP21
**Tried.** After the dynamic-CDKI module baked, recover a discrete p27-high transient-G0 state in MB by tuning clearance ([`cdki-baked-g0-g1-behavior`], `ws1_transient_g0`).

**Result / why abandoned.** The reviewer's fast-KPC + the CdP21 bimolecular commitment lever produce a **sharp R-point** that **removed the prior model's discrete G0**: MB G0frac = 0%, pRb never < 1.66, free-p27 ~0, birth-p27 has no effect, and it is **NOT recoverable via clearance tuning**. MB's 26h cycle is a longer *G1* (21.3h vs GNP 10.8h) set by high INK4 — **not** a p27-high pause. A discrete transient-G0 now genuinely requires a population-heterogeneity layer.

### 3.3 The `kDeP21Cd` saturated-switch p27 clearance — do NOT tune it
**Tried / result.** The two-step-Rb branch added a mitogen-driven p27 clearance `kDeP21Cd*Cd*commit_gate/(...)`. In MB (Cd~8) it **saturates** and strips committed p27 from 0.6 → ~0.0017, collapsing MB p27-high/G0 from 18–42% to **0%** ([`mb-g0-regression-and-p27-reconception`], git-verified at commit 747d989).

**Why abandoned as a knob.** `kDeP21Cd` is a **saturated switch** — 0.08–0.175 all give 0% G0; only *exactly* 0 flips to 100% arrest. It cannot be used to dial G0. Restored (temporarily) via the MB **birth-p27 crutch** (P21_DIV_MB=1.8), later reduced and then retired (§3.4). Underlying **hard tension:** MB G0-rich and the 5.07 fold are coupled through EZH2 (longer G0 → less EZH2 → de-repressed CyclinD1 → fold breaks toward ~7–10); fold-safe ceiling is ~16–20% G0.

### 3.4 Birth-p27 crutch — RETIRED
**Tried / result.** MB-specific born-high p27 (P21_DIV_MB=1.8, later 1.39) to latch the Skp2-p27 toggle into a growth-timed G0 window.

**Why abandoned.** A non-mechanistic crutch (the model needs ~3× birth-p27 vs the ~1.33× abundance-weighted transcript pool). Retired in the transient-G0 reframe (P21_DIV_MB→0.6, GNP baseline G0 ~0 as a firm anchor); superseded by CDK4/6-activity-gated reversible G0.

### 3.5 The Phase-1 mitogen tracker (birth-p27 = f(mitogen deficit)) — SHELVED (backwards)
**Tried / result.** `with_mitogen_tracker`: set birth p27 from the mitogen deficit. Made MB birth-p27 **~30× BELOW** GNP — **backwards** vs the data (MB p27 is *higher*). Kept opt-in/off as a documented dead-end.

### 3.6 The mother-G2 p27 integrator (Overton bifurcation) — NEGATIVE RESULT
**Tried / result.** A mother-cell-G2 p27 integrator to reproduce the Spencer/Overton birth-p21→CDK2^low-vs-inc fate bifurcation. It **saturates to ≈1** and does not reproduce the bifurcation ([`v44-mother-g2-negative-result`]). Fix identified (low-pass tracker with a low floor) but the integrator itself is a documented dead-end.

---

## 4. Structural walls in the cell-cycle machinery

### 4.1 HU-S vs HU-G2 folds — a STRUCTURAL either-or (28/28 is unreachable honestly)
**Tried.** Pass both hydroxyurea checks simultaneously: MB HU-S fold (cells accumulate in S, target 1.36) and MB HU-G2 fold (HU depletes G2, target 0.23). 50-point Pareto sweep over KmHU_fire × vmin_fork × kSyDna + an independent adversarial refuter agent ([`v44-hu-s-g2-structural-limit`], `_frontier_hu.py`).

**Result.** **0/50 pass both**; closest joint = HU-S 0.740 / HU-G2 0.324 (both miss by a hair). The default is 27/28 with HU-G2 = 0.356 as the lone near-miss.

**Why abandoned (the load-bearing fact).** In-model **BrdU signal and DNA-completion flux are the same quantity** (`BrdU = vfork·aRc`; completion = kSyDna·BrdU). So "counted as S (BrdU+)" ⟺ synthesis flux > 0 ⟺ the cell *is* completing replication and *will* reach G2. "G2 depleted" ⟺ flux ≈ 0 ⟺ BrdU-negative ⟺ not counted as S. **Mutually exclusive by construction.** Every classifier-level fix relocates the failure (SYN_THR↓ balloons DMSO-S to 30.8%; kSyDna↓ fails DMSO-S; pulse-BrdU overshoots HU-S to 7–9×). Decisive biology: real HU-stalled cells incorporate **little** BrdU → which *favors* the flux definition already in the code. So forcing 28/28 via a pulse-BrdU redefinition is **less** biologically faithful. **Recommendation: accept HU-G2 ≈ 0.356 as a noise-consistent flow-cytometry near-miss; do not hack the classifier.** This wall is the perennial casualty of every recalibration (§2) and reappears as the only fail (with HU-S) in the two-cyclin (§6.1) and dynamic-CDKI (§5) baked models.

---

## 5. Individual CDKIs at data abundances — p21/p57 arrest the cell

**Tried.** Model each expressed CDK-inhibitor as its own dynamic species (INK4: p18, p19; CIP/KIP: p27, p21, p57), replacing the lumped static p16/p18 + single "P21," then engage p21 and p57 at their measured transcript abundances (~0.49 / ~0.24) ([`dynamic-cdki-module`], ~560 sims, two expert-review rounds).

**Result.** Engaging p21/p57 at data levels breaks validation to **15/32** — they **over-inhibit CDK2** (via the shared p27 binding constants) and drive p27 to near-zero, arresting the cell *even after* removing them from the CDK4/6 brake and lowering their CDK2 potency.

**Why abandoned.** Two reasons converge: (1) the reviewer's predicted **degeneracy** — p21/p57 are near-unidentifiable against aggregate CDK2 activity, and p21 "lacks the key tyrosine, poor CDK4/6 activator"; (2) the biological point that **mRNA abundance overstates the active nuclear CDK-inhibitor pool** for p21/p57 (heavy post-translational regulation). Resolution: p21/p57 are carried **expressed-but-inert**; the functional CDK regulation is p27 + INK4 + the CdP21 commitment buffer. (p57 further de-emphasized — largely a Sox2+ quiescent-cell CKI, carried OFF.) A related numerical dead-end: moving the PCNA/CRL4-Cdt2 arm onto p21 dragged p21a **negative** (−0.4, corrupting the sim) because iPcna/iRc are replication-coupled *shared* species — reverted; p21 keeps Cdt2 via the aRc proxy. The dynamic-CDKI module *did* bake (30/32) — but with p21/p57 inert, which is itself the negative result.

**Also demoted here — "Option B" CDK2 bypass.** The Guiley/LaBaer/Sherr Cip/Kip-sequestration CDK2-bypass mechanism (`option_b_cdk2_bypass_reading.md`) was read up and prepared as an oscillator-safe module, then **demoted**: it applies only to genuine *acquired* resistance (CCNE1-amp persisting at saturating dose over weeks), not JP's reversible in-vitro palbo data. A naïve constitutive CDK2→Rb term breaks the oscillator (`EZH2_palbo_cycle_gating_tension.md`).

---

## 6. Two-cyclin and MYCN extensions

### 6.1 CyclinD2 marked *equally* to CyclinD1 — DATA-INCOMPATIBLE (vismo direction flips)
**Tried.** JP: "CyclinD2 is also H3K27me3-marked, it should be." Add the *same* PRC2_rep Hill repression to Cd2 synthesis (`d2_marked` flag) ([`cyclind1-d2-two-cyclin-model`]).

**Result.** **27/32 — all 3 D2 folds break and the vismodegib direction FLIPS**: CyclinD2 GNP+HHi 0.66→**1.07** and MB+HHi 0.59→**2.03** (D2 goes *up* under vismodegib, vs the observed ~34% drop).

**Why abandoned.** Under Hh-block the cell arrests → EZH2 falls → PRC2 repression falls → D2's large Gli-*independent* basal floor **de-represses upward**; D1 instead falls because it is Gli-dominated. A 200-candidate re-parameterization confirmed **infeasibility**: best = 29/32, only 2/3 D2 folds — a clean **three-way trade-off** (hit the vismo drops ⟹ D2 Gli-dominated ⟹ MB/GNP collapses to ~1.0; hit MB/GNP ⟹ big Gli-independent elevation ⟹ de-represses up under vismo). "Marked + buffered + drops-under-vismo" is over-constrained for a *strong dynamic* mark. Reconciliation ([`h3k27me3-kinetics-vs-local-repression`]): CCND2 *is* biologically a PRC2 target, but its mark must be **static/developmental** (or a kinetic buffer) that does not dynamically gate transcription. Default: D2 **unmarked** (30/32); `d2_marked` stays OFF.

### 6.2 MYCN acute Gli1-driven autoregulation latch — UNREALIZABLE
**Tried.** JP wanted an MB-only MYCN self-activation latch to explain why the acute vismodegib response *shrinks* as baseline Hh rises (GNP −20%, MB −11%) — a developmental threshold-lock. ~50-candidate search + fine threshold map + low-MYCN-init + robust/slow integrators ([`mycn-elevated-expression-not-amplified`]).

**Result / why abandoned.** **Unrealizable.** The buffer strength (self-sustaining high state) and the GNP/MB discrimination are the **same knob** (K_auto): low K_auto → both latch (MB/GNP ~1.2); high K_auto → GNP stays low but MB no longer holds under HHi (MB+HHi drops to 0.26). **No window between** — GNP's own measured 20% drop puts it too close to any threshold, and GNP/MB MYCN are only ~35% apart at the decision point, too close for a bistable switch to both separate them and sustain a buffer. Also numerically fragile (steep Hill → CVODE failures); slow-integrator + gentle-autoreg fail identically (negative basal analytically). **The lock is developmental, not an acute switch.** Resolution (JP-approved "relabel and move on"): an additive **MYCN_expr** elevated-expression term (MB=0.54, Gli-set) replaces the old `MYCN_amplification` multiplier — algebraically neutral (26/29 unchanged); `mycn_autoreg` flag deprecated/inert. Note the framing correction: SHH-MB is **not** MYCN-amplified (that's Group-3); MYCN is elevated as a Hh/Gli target.

---

## 7. Withdrawal-memory dynamics

### 7.1 A sustained withdrawal oscillation — OVERDAMPED GOVERNOR (0 oscillations in ~200 sims)
**Tried.** JP's proposed delayed negative-feedback loop: CyclinD1↓ → G1 lengthens → H3K27me3 accumulates → transcript drops harder → EZH2 falls → mark decays → transcript rebounds — which might **oscillate**. Probed in the two-cyclin model incl. a 4-agent adversarial-verification workflow (staircase + rate-rule ramps, K_prc2/20, n_prc2→10, del_mk×20, slow EZH2, horizons to 1000h) ([`withdrawal-loop-overdamped-governor`], `fig_withdrawal_reservoir_loop.png`).

**Result.** **Zero sustained oscillations** across ~200 sims. Cycle length is metronomic (CV ≈ 0.00, 17.2h) and **bit-identical** at tight solver tolerance (physical, not numerical). Only ever a single transient overshoot, never ringing.

**Why it can't ring.** transcript→protein is *faithful* (both fall to 38%; the reservoir extends coasting *duration*, not swing), but the overdamp is at **protein→CDK4/6-Rb drive**: the Michaelis term is deep in saturation (Cd 1.4–3.8 ≫ K_CdRb 0.319), so a 62% protein drop compresses to a **16% drive drop**; and the cycle is **growth-pinned** (period = ln2/μ_eff), with G1-lengthening *latent* (K_g1len=0.6 < operating Cd). The residual swing never reaches the period. To make oscillation/graded lengthening observable would need raising K_g1len toward the operating Cd, desaturating the drive, or adding a stochastic/population layer or an explicit fast-feedback delay species — **none of the current continuous machinery can ring.**

### 7.2 Graded G1-lengthening *driven by EZH2 reduction* — NOT supported (EZH2 banks)
**Tried.** JP's overnight hypothesis: the large CyclinD1/D2 pool + **EZH2 reduction** (de-repressing CyclinD1) keeps cells dividing with progressively lengthening G1 as mitogen falls ([`withdrawal-graded-g1-lengthening`], `fig_g1_lengthening_withdrawal.py`).

**Result.** The *phenotype* ("lengthen G1 but not exit") **is** achievable and validation-compatible (28/29) — but via **growth-slowdown** (the latent μ_eff lever: mu_min_frac 0.4→0.1, K_g1len 0.6→0.8 lengthens the period smoothly 24h→84h then arrests, staying Shh-dependent). JP's specific **EZH2-reduction driver is NOT supported**: the model's EZH2 does *not* reduce during withdrawal (it **banks** — slow τ~14h can't track the falling mitogen). Forcing it to reduce (faster turnover) either (i) was an artifact of an uncompensated lower EZH2 *level*, or (ii) with the level held, gives *fewer* divisions, and (iii) **breaks validation to 23/29** (kDeEZ is pinned by the HU-S/G0 data). So the memory is growth-slowdown, not EZH2-reduction dynamics — not baked.

### 7.3 Topology A (mark represses CyclinD1) as a withdrawal memory — PROVEN it cannot
**Tried / result.** Test whether the v44 topology (mark represses CyclinD1, negative feedback) can produce a buffering memory that sustains division through withdrawal ([`prc2-buffering-memory-topologies`], adversarially verified, 14 tests).

**Why abandoned.** Airtight bound: for **any** net-repressive coupling with rep(Mk=0)=1, feedback-off holds the drive at max at every instant → D_on(t) ≤ D_off(t) pointwise → banked divisions and half-life are always ≤ off. Topology A buffers only the *steady-state dose-response*, and actually gives *fewer* post-withdrawal divisions (it sharpens arrest). A true withdrawal memory needs **topology C** (mark *sustains* proliferation by silencing differentiation genes = a positive-feedback latch) — which is the bistable latch that is data-excluded (§1.1). All sub-cycle-dip persistence is the pure Rb-E2F R-point (cell-cycle) memory, mark-independent.

### 7.4 The slow CyclinD1 mRNA reservoir (t½≈11h) — CHALLENGED by measured half-life
**Tried / result.** The baked "memory = slow CyclinD1 transcript reservoir" (`k_Cd_mRNA_deg`=0.001, t½≈11h) was challenged by the RNAdecayCafe atlas: **measured CCND1 t½ = 1.5–2.5h** (28th percentile, moderately labile) — the reservoir is ~5–7× too slow ([`ccnd1-halflife-data-vs-reservoir`]).

**Why (partially) a dead-end.** A *fast* transcript (t½=2h, tx rescaled) is **validation-neutral (28/29)** and *still coasts ~1 division* after Hh-off — because the ~11–12h previous-cycle memory (Ho 2020) is the **committed Rb-E2F state**, not the transcript half-life. So the slow reservoir is a stand-in for "big labile pool + committed cycle." Reducing it to match the data would fix the half-life and the withdrawal-drop timing while keeping the 1-cycle coast — but it **relocates the memory from the transcript to the commitment** and weakens the "reservoir low-passes Shh noise" story (filter τ drops 16h→~2h). Decision pending; not yet re-baked.

---

## 8. Conceptual corrections (framings that turned out backwards)

| Framing (held) | Correction | Where |
|---|---|---|
| **EZH2 is tumor-suppressive via CyclinD1** (EZH2 loss → more proliferation/tumors) | **Wrong.** EZH2 has two arms; the dominant developmental arm (not modeled) represses *differentiation* → cKO causes premature differentiation/exit → progenitor pool depleted → *Ptch⁺/⁻;EZH2-cKO does NOT increase MB*. The model holds only the CyclinD1-*governor* arm (raises the mitogen threshold). The arms **oppose** at tumor initiation. | [`ezh2-two-arms-cko-scoping`] |
| **MB has it backwards: MB mark > GNP mark** (model gave MB 0.65 > GNP 0.43) | ChIP shows **MB < GNP** (~half). Model had no eraser and no Gli→rate coupling; reconceived as a writer(cycle) vs Gli-Jmjd3-eraser(mitogen) race. | [`v44-gli-jmjd3-eraser-reconception`] |
| **SHH-MB is MYCN-amplified** | It is **not** (Group-3 is); MYCN is *elevated expression* as a Hh/Gli target. | §6.2 |
| **CDK4/6i resistance = a resistant clone** (15% low-Rb mixture) | **Retired → dose-response.** Palbo dividing fraction (16–18% @ 1µM, ~2% @ 5µM, both reversible/pRb⁺) is a graded residual-CDK4/6-activity dose-response over a CyclinD1/CKI population — no CDK2-bypass "option B," no senescence, no clone. | [`v44-transient-g0-reframe-and-cdk46i-escape`], `cdk46i_resistance_research.md` |
| **CDK4/6i escape fractions are emergent** | **Honesty caveat:** the 20%/15% subpopulation sizes are JP's culture *inputs* (a soft target); the model verifies the mechanisms + reversibility, **not** the fractions (CyclinD1 heterogeneity moves zero cells). | §7 of the reframe doc; adversarial 14-agent review, commit 499bf7d |
| **Stopping division preserves the mark via less dilution** | **No** — EZH2 collapse dominates, mark lost anyway; division-coupled dilution is a *cost* to the memory, not its source. | [`prc2-buffering-memory-topologies`] |
| **Alternative pocket proteins (p107/p130) as a CDK4/6i resistance route** | **Opposite** — they build the DREAM repressor that *deepens* arrest; "do not build a p107/p130 term." | `cdk46i_resistance_research.md` |

---

## 9. The parallel successor: v45 stochastic-commitment — COMPLETE but SHELVED
**Tried / result.** A parallel prototype (`v45_stochastic_commitment.py`) replacing v44's non-mechanistic size gate with an emergent **bistable CDK2-p27 commitment toggle** + Rb-dilution G1 timer + noisy-birth ensemble. It **works and validates (17/17)** — reproduces the 2N count-fraction (67% vs 68.2% target) from birth noise alone, with the MYCN floor emerging ([`v45-stochastic-commitment`]).

**Why shelved.** Not a failure — a **deliberate parallel track**. v44 remains the committed working model; v45 is the eventual home for exactly the population/stochastic phenomena that v44 cannot express (§3.1, §3.2, §1.4). Notably, [`v44-sizegate-is-rb-dilution`] found v44's mass gate is already ≡ the Zatulovskiy Rb-dilution R-point, so porting v45's toggle into v44 is "cosmetic" — reinforcing that the real v45 payoff is the *stochastic ensemble*, deferred.

---

## Cross-references
- The data walls that do the killing: [Data & Evidence](07_data_and_evidence.md), [Calibration & Validation](06_calibration_and_validation.md).
- The mechanisms that *did* survive: [Mechanisms Explored](04_mechanisms_explored.md), [Model Architecture](02_model_architecture.md).
- How these dead-ends fall in time: [Evolution Timeline](03_evolution_timeline.md).
- Open items several dead-ends point to (population layer, allele-switch, CCND1 half-life re-bake): [Reviewer & Open Questions](10_reviewer_and_open_questions.md).
