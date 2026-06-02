# v44 verification: can the S-checkpoint be made concentration-dependent? — findings

*Run while you slept. Per the plan, I tested step 2 (a concentration-dependent CHK1→Cdc25
checkpoint) BEFORE touching the builder. The verification returned a clear, consistent
negative — and explains why. Bottom line up top:*

> **In the Gérard–Goldbeter model topology, S-phase duration is NOT a concentration-controllable
> quantity. It is a structurally robust ~4.7 h set by the relaxation oscillator. No molecular
> concentration — including the biological checkpoint target (Cdc25/Cdk1-Tyr15) — can lengthen S.
> Only the global clock (eps), or eps gated to the S-window (our v43), can. A genuinely
> concentration-dependent S-checkpoint requires a different model with explicit DNA replication,
> not the GG core.**

This means the v43 eps-in-S checkpoint is **not a hack to be replaced — it is the correct (and
essentially only) abstraction for stretching S in this model class.**

---

## What I tested (3 experiments, all on the v42/v43 core, eps phase-gating disabled so HU drove only the test)

### 1. CHK1 → Cdc25B inhibition (the deep-research "best path") — `diag_v44_checkpoint_test.py`
Added a replication-progress variable `Rep` (integrates during S at fork speed `v_fork(HU)`;
HU slows forks), a `CHK1 := sphase·K^m/(K^m+Rep^m)` (active = replicating AND incomplete), and
gated Cdc25B activation by `K_inh/(K_inh+CHK1)`.

| arm | period | S_dwell | Cb_max |
|---|---|---|---|
| baseline | 19.6 h | 0.241 | 0.19 |
| constant Cdc25B ×0.5 / ×0.2 / ×0.1 | 20.0 / 20.1 / 20.1 h | 0.241 (flat) | 0.75 / 1.53 / 1.81 (inflates) |
| **dynamic CHK1→Cdc25B, HU 0→2** | 19.9–20.0 h (flat) | **0.241 (flat at every HU)** | — |

→ Even maximal, dynamic, replication-gated Cdc25B inhibition does **nothing** to S. Confirms the
relaxation-robustness result: inhibiting a switch *activator* can't hold the transition.

### 2. Hard replication-LICENSING gate on Cdk1 / APC — `diag_v44_licensing_test.py`
Gated Cdk1 (Mb) activation, and separately Cdc20 (APC) activation, **directly** by a steep
license `Rep^8/(0.7^8+Rep^8)` (mitosis literally cannot fire until replication completes).

| arm | period | S_dwell | Cb_max |
|---|---|---|---|
| constant Mb ×0.1 / ×0.02 | 20.1 h | 0.239 / 0.237 (flat) | 1.64 / 1.99 (inflates) |
| dynamic Mb-license, HU 0→2 | 19.7 h (flat) | 0.241 (flat) | 0.21 |
| dynamic Cdc20-license, HU 0→2 | 19.7 h (flat) | 0.241 (flat) | 0.20 |

→ Throttling Cdk1 activation to **2%** does not delay mitosis — it just makes active Cdk1 tiny
while free CycB balloons, and the cycle completes on schedule anyway. A hard licensing gate on the
mitotic module has **zero** leverage on S.

### 3. Where IS S-duration set? — sensitivity scan `diag_v44_sensitivity.py`
Multiplied each rate constant by 0.5× / 2× and measured period + S_dwell:

| module | params | effect on period | effect on S |
|---|---|---|---|
| **Rb–E2F–CyclinE–CyclinA (pacemaker)** | V1e2f, vse2f, kca, Vm1a, kce, V1 | large (or breaks oscillator) | changes S_dwell, but… |
| **CyclinB–Cdk1–Cdc25B–Cdc20 (mitotic)** | vcb, Vm1b | **none (period & S flat)** | **none** |

The mitotic module — *where the biological intra-S checkpoint acts* — has **no leverage**. The
pacemaker is the Rb–E2F–CyclinA loop: the CyclinA/S window is terminated by CyclinA inactivating
E2F (`V1e2f`), a negative feedback, **not** by mitosis.

**And the punchline:** even the pacemaker knobs don't cleanly lengthen S. They lengthen the
*period* by adding time to **G1**, while **S-phase absolute duration stays ~4.7 h**:

| knob | period | S (hours) |
|---|---|---|
| baseline | 19.6 h | 4.7 h |
| kca ×2 (more CyclinA synth) | 24.8 h | 4.7 h (extra time → G1) |
| vse2f ×2 | 27.5 h | 4.7 h (→ G1) |
| Vm1a ×2 | 28.6 h | 4.5 h (S *shorter*) |
| V1 ×2 | 19.0 h | 5.0 h |

S-phase duration is a **structural invariant** of this oscillator. Concentration changes
redistribute time into G1; they cannot stretch S.

---

## Why this happens (the dynamical reason)

S/G2→M is a **fast relaxation discharge slaved to the Rb–E2F–CyclinA pacemaker**, not a
slow-charge leg waiting on cyclin accumulation. There is no slow, rate-limiting variable inside S
for a concentration to throttle. The CyclinA pulse is **self-terminating** via E2F autoinhibition
on a fixed timescale, independent of the mitotic apparatus — so inhibiting Cdc25/Cdk1/APC (which
only affects the slaved mitotic output) leaves S untouched. This is the same robustness that makes
the period insensitive to cyclin amplitude — now shown to extend to *dynamic* and *hard-gated*
concentration levers, not just constant ones.

By contrast, **eps works** precisely because it slows the *pacemaker's own clock* in the S-window —
the one thing that actually sets the traverse time.

---

## Revised options (this changes the 3-step plan)

The 3-step plan's **step 2 (concentration-dependent checkpoint on Cdc25/Cdk1) is not viable in the
GG core.** Step 1 (drop eps → real-time h⁻¹) is still fine but cosmetic — it won't make S
concentration-tunable. So the real fork is:

**Option A — Keep the eps-in-S checkpoint (v43), reframed as principled, not phenomenological.**
The finding *justifies* it: S-duration in a relaxation oscillator is a clock property, and "an
intra-S checkpoint slows the S-phase clock" is the faithful coarse-grained statement. Cheapest,
already validated 12/12, hits the 1.31× anchor. We'd just document *why* it's the right abstraction
(this analysis) and optionally drop the bare `eps` for native h⁻¹ rates so it reads as real time.

**Option B — Rebuild S-phase as an explicit, rate-limited process (the only route to true
concentration-dependence).** Replace the instantaneous G1/S→S→G2 portion with explicit DNA
replication: origin licensing/firing, fork number × fork speed = replication rate, "% genome
replicated" as a slow state variable, and a **hard requirement that mitosis cannot initiate until
replication = 100%**. Then HU lowers fork speed → replication takes longer → S lengthens
concentration-dependently, and CHK1→Cdc25 enforces the gate. *But:* the research found **no
ready-made SBML** for this coupled to a mammalian CDK oscillator — it's a genuine model-building
project (closest scaffolds: Heldt 2018 for an eps-free CDK core, but it has no replication-fork
layer either; the replication layer would be new). This is a v44/v45-scale effort.

**Option C — Hybrid: keep the GG pacemaker, bolt on an explicit replication timer that *hard-gates*
the E2F→S commitment, not the mitotic exit.** Since the pacemaker is Rb–E2F–CyclinA, the
leverage point is the *entry/commitment*, not Cdc25/Cdk1. A replication timer that must complete
before the E2F/CyclinA program can re-arm could, in principle, insert concentration-dependent delay
at a node the pacemaker actually feels. Unproven — would need its own verification test (and risks
just relocating the robustness problem).

---

## My recommendation

**Go with Option A, but with the new framing**, and treat Option B as a separate, scoped future
project if you want true mechanistic replication. The verification did its job: it saved us from
building a concentration-dependent Cdc25 checkpoint that *cannot work* in this model, and it turned
the eps-in-S approach from "a hack we apologized for" into "the demonstrably correct abstraction for
this oscillator class." If you want concentration-dependence badly enough to justify it, Option B is
the real path — but it's a new model, not a checkpoint swap.

Scripts: `diag_v44_checkpoint_test.py`, `diag_v44_licensing_test.py`, `diag_v44_sensitivity.py`
(all untracked, reversible). Nothing committed.
