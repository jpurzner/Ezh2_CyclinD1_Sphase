# How the model runs a cell cycle — phase durations, what rises, what falls

**Model:** v43 = published v42 (MYCN–EZH2–CyclinD1–Hedgehog, Gérard–Goldbeter cell-cycle
core) + the intra-S checkpoint (φ(HU) phase-gated clock) + the EZH2→(Me+Ma) transcription
rewire. All numbers below are from the calibrated DMSO and HU=1.0 simulations
(`calibrate_ezh2_rewire.run`, `validate_v43.py`).

This document answers three questions:
1. **What sets the length of each phase?**
2. **Walk me through one cycle — what goes up, what goes down, in order.**
3. **Where does our checkpoint plug in, and why does it lengthen S but not G1?**

---

## 1. The three layers of the model

The model is best understood as three stacked layers that talk to each other through a
small number of variables:

```
  LAYER A — GROWTH INPUT (Hedgehog + MYCN)        … runs in "real time" (no clock factor)
        SHH → Ptch1 → Smo → Gli → ┐
                          MYCN  →  ┴→  CycD1 mRNA  ───────────┐
                                         ▲                     │ sets how hard
                                         │ repressed by EZH2   │ G1 is pushed
  LAYER B — CELL-CYCLE OSCILLATOR (Gérard–Goldbeter)          ▼
        CycD1 → [Rb–E2F switch] → CycE → CycA → CycB → mitosis → reset
        (every reaction carries the global clock factor `eps`)
                                         │ E2F + (CycE,CycA) activity
  LAYER C — EPIGENETIC READOUT + FEEDBACK                     ▼
        E2F × (Me+Ma) → EZH2 transcription → EZH2 protein ──┘ (represses CycD1)
```

- **Layer A is quasi-static.** None of the Hedgehog/MYCN reactions carry the `eps` factor,
  so they equilibrate fast relative to the ~20 h cycle. Their job is to set a near-constant
  **CycD1 mRNA level** — the "throttle" on G1. With SHH on, CycD1 mRNA ≈ 2.6 (arbitrary units).
- **Layer B is the engine.** Every cell-cycle reaction is multiplied by `eps = 150`. This is a
  pure **time-scaling**: it makes the network's intrinsic period (~0.13 dimensionless units)
  map onto **~19.6 h**. `eps` is the single knob that sets the *speed of the whole clock*.
- **Layer C reads out** the engine (EZH2 transcription is driven by E2F and CDK2 activity) and
  **feeds back** (EZH2 represses CycD1 → throttles Layer A → slows G1).

The intra-S checkpoint we added lives **between B's speed knob and its S-phase state**: it makes
`eps` *drop* specifically while the cell is replicating. More on that in §3.

---

## 2. The cast (what each variable means)

| Symbol | Biology | Role / when it's high |
|---|---|---|
| `Cd_mRNA`, `Cd` | Cyclin D1 mRNA / protein | Growth input; accumulates through G1 |
| `Md` | **CycD1–CDK4/6** active | G1 — starts phosphorylating Rb |
| `pRB`→`pRBp`→`pRBpp` | Retinoblastoma, mono-→hyper-phosphorylated | Rb is the brake; phosphorylation releases it |
| `E2F` | E2F transcription factor (free) | **Spikes at the G1/S restriction point** |
| `Me` | **CycE–CDK2** active | G1/S boundary — trips the Rb switch fully |
| `Ma` | **CycA–CDK2** active | **S/G2 marker** — our "replication window" |
| `Cb` / `Mb` | free Cyclin B / **CycB–CDK1** active | **M-phase marker** / mitotic trigger |
| `Pe`,`Pa`,`Pb` | Cdc25 E/A/B phosphatases | Activate CycE/A/B–CDK (positive feedback) |
| `Cdc20a` | APC/C–Cdc20 | **Destroys CycA & CycB → mitotic exit** |
| `Cdh1a` | APC/C–Cdh1 | Keeps CycB low in G1; switched off in S/G2 |
| `p27` | CDK inhibitor (CKI) | High in quiescence/early G1; degraded to permit S |
| `Skp2` | SCF substrate-recognition | Rises with E2F; degrades p27 |
| `EZH2` | H3K27 methyltransferase | Readout; **represses CycD1** |
| `φ(HU)`, `sphase`, `eps` | our checkpoint | φ slows `eps` while `sphase`≈1 (CycA-high) |

---

## 3. What sets the length of each phase (the core idea)

The Gérard–Goldbeter core is a **relaxation oscillator**: a fast bistable *switch* (the
Rb–E2F–CycE/A positive-feedback loop) that flips ON, followed by a slow *negative feedback*
(CDK → APC/C → cyclin destruction) that resets it. Two consequences are essential to
understanding phase durations — and they are *not* obvious:

**(a) The period is set by the global clock `eps` and the network topology — NOT by the
amplitude of any one cyclin.** We verified this directly: scaling any single activation,
degradation, or synthesis flux (CycA, CycB, Cdc20, …) by 0.6× moves the period by ≤1%
(`diag_phi_levers.py`). The bistable switches fire at the same *time* regardless of how fast
an individual reaction runs — they just reach a different *height*. **Only changing `eps`
moves the period** (halving `eps` ≈ doubles the period). This is why our checkpoint had to act
on `eps`, not on a cyclin flux (the earlier "scale CycA activation by φ" approach left the
period flat — it only shrank the CycA peak).

**(b) Phase lengths are set by how long each transition takes to *arm*, not by the switch
itself.** Concretely:

| Phase | ~Duration (DMSO) | What determines it |
|---|---|---|
| **G1** | **14.8 h (75%)** | How fast CycD1–CDK4/6 + CycE hyperphosphorylate Rb and trip the E2F switch. Faster CycD1 (more SHH/MYCN, less EZH2) → shorter G1. p27 must be degraded first. |
| **S + G2** (CycA-high) | **4.8 h (25%)** | From the E2F spike to mitosis: CycA–CDK2 build-up and the CycA→CycB hand-off. |
| **M** (CycB-high) | **4.3 h (22%)** | CycB–CDK1 active until APC/C–Cdc20 destroys the cyclins. |

(The S and G2 windows are *telescoped* in this simplified core — CycB starts rising only
~0.5 h after CycA, so the strictly "CycA-high / CycB-low" slice is brief. We therefore treat
the whole **CycA-high (`Ma`) window as the operational replication/S–G2 phase**; that is the
window our checkpoint stretches and the window in which EZH2 is read out.)

The **restriction point** (the irreversible commitment) is the E2F spike, which in the DMSO
cycle occurs at **t ≈ 16.1 h — 82% of the way through the cycle**. Everything before it is G1;
the last ~3.5 h is the rapid S→G2→M cascade.

---

## 4. Walk through one cycle (DMSO, period = 19.6 h)

Times are measured from a mitosis (t = 0 = CycB/CycA collapse, two daughter states reset).
"↑/↓" = rising/falling.

### t = 0 – 2.5 h — Mitotic exit → early G1
- **CycB↓, CycA↓ (just destroyed)** by APC/C: `Cdc20a` spiked at mitosis and degraded `Ca`
  (CycA) and `Cb` (CycB). CDK activity collapses → the cell is now 2N, in G1.
- **APC/C–Cdh1 (`Cdh1a`) ON**, holding CycB down. **`p27`↑ to its peak (~65)** — the CKI floods
  in, keeping any residual CDK off. This is the "low-CDK" ground state.
- **Rb (`pRB`) is hypophosphorylated and binds E2F** (`pRBc1`/`pRBc2` complexes) → **free E2F is
  low (~0.2)**. The cell-cycle transcription program is OFF.

### t = 2.5 – 12 h — G1 build-up (the long, rate-limiting stretch)
- **CycD1 protein (`Cd`) and `Md` (CycD1–CDK4/6) accumulate**, fed by the Hedgehog/MYCN input
  through `Cd_mRNA`. `Cd_mRNA`/`Cd` rise to their peak around **t ≈ 12.7 h**.
- **`Md` begins phosphorylating Rb: `pRB → pRBp`.** This is the *first, reversible* dent in the
  brake. Partly-phosphorylated Rb still holds some E2F, so E2F is still low.
- **`Skp2`↑ (tracking the first E2F leak) starts degrading `p27`** → p27 falls from its peak.
  As p27 drops, CDK2 complexes are freed from inhibition. **This is the slow arming of the
  switch** — and the step EZH2 modulates (via CycD1).

### t ≈ 12 – 16 h — The restriction point (G1/S switch flips)
- Once Rb is phosphorylated enough, a little **free E2F appears → drives CycE (`Me`) synthesis →
  CycE–CDK2 phosphorylates `pRBp → pRBpp` (hyperphosphorylated) → releases *more* E2F → drives
  *more* CycE.** This is the **positive-feedback loop** firing. It is fast and switch-like.
- **`E2F` spikes from ~0.2 to ~50 at t ≈ 16.1 h.** `Me` (CycE–CDK2) peaks ~16.5 h, `pRBpp`
  peaks ~16.5 h. **This is commitment — the restriction point.** Cyclin E's Cdc25 (`Pe`) is
  co-activated, sharpening the switch.
- p27 is now fully degraded; the brakes are off.

### t ≈ 16 – 19 h — S / G2 (CycA window)
- **E2F now drives CycA (`Ma`) synthesis.** `Ma` (CycA–CDK2) climbs to its peak (~1.15) at
  **t ≈ 19.6 h**. CycA–CDK2 is the **S-phase engine** (replication) and the marker we read.
- Two negative-feedback arms switch on:
  - **CycA phosphorylates E2F → `E2Fp` (inactive)** → the E2F transcriptional pulse is
    self-terminating (E2F falls back after its spike).
  - **CycA + CycB inactivate APC/C–Cdh1 (`Cdh1a`↓)** → this *permits* CycB to start
    accumulating. (`Cdh1a` actually peaks ~15.6 h then falls as Ma/Mb rise.)
- **EZH2 transcription is ON here** (it is gated on E2F × (CycE+CycA); see §5), so **EZH2 protein
  rises through this window, peaking ~t ≈ 18 h (~0.82)**. This is the high-synthesis window.

### t ≈ 19 – 19.6 h — Mitosis (M)
- **CycB (`Cb`) accumulates and CycB–CDK1 (`Mb`) activates** (its Cdc25, `Pb`, provides positive
  feedback). `Cb` peaks (~0.19) at mitosis. The cell is now 4N, in M.
- **CycB–CDK1 activates APC/C–Cdc20 (`Cdc20a`)**, which **destroys CycA and CycB** →
  CDK activity collapses → **mitotic exit**, and we are back at t = 0 of the next cycle.
- The whole S→G2→M cascade (E2F spike to division) takes only **~3.5 h** — the cycle is
  **G1-dominated**.

**One-line mental model:** *G1 is a slow capacitor charging (CycD1 vs. p27/Rb); the restriction
point is the spark; S/G2/M is a fast discharge (CycA→CycB→APC) that resets the capacitor.*

---

## 5. The EZH2 readout and the feedback loop

**Readout (rewired in v43).** EZH2 transcription is driven by
```
k_basal + k_E2F · E2F/(K+E2F) · [ w_Me · Me/(K_Me+Me)  +  (1−w_Me) · Ma/(K_Ma+Ma) ]
```
i.e. it turns on at the **G1/S boundary** (CycE `Me` arm — restarts transcription as S begins)
and stays on through **S/G2** (CycA `Ma` arm). Because EZH2 protein turns over with τ≈5 h, it
**ramps up within each cycle**, producing the measured within-cycle gradient
**G1 < S < G2** (validated 1.13 / 1.31 / 2.37 × the quiescent G0 level). In quiescence/arrest
both `Me` and `Ma` collapse → the gate closes → EZH2 falls (this reproduces the GDC0449 / HHi
46% reduction).

**Feedback.** EZH2 multiplies the CycD1 transcription rate by `K/(K + EZH2)`. So **more EZH2 →
less CycD1 → slower G1.** This is two things at once:
- a **tonic brake** (removing it shortens the period ~19.6→17.9 h and ~doubles CycD1), and
- the **transducer of the central hypothesis**: *longer S → more EZH2 → stronger CycD1
  repression → longer subsequent G1.*

---

## 6. Where the intra-S checkpoint plugs in, and why it lengthens S but not G1

**Mechanism.** HU depletes dNTPs → slows replication forks. We model fork competence as an
instantaneous function of dose:
```
φ(HU) = φ_min + (1−φ_min)/(1 + (HU/K)^h)        φ ∈ [φ_min, 1],  φ=1 at HU=0
```
and we let φ **slow the cell-cycle clock specifically during replication**:
```
sphase := Ma^n_S /(K_S_phi^n_S + Ma^n_S)        ≈ 1 when CycA is high (S/G2), ≈ 0 otherwise
eps    := eps0 · (1 − sphase·(1−φ(HU)))          eps → eps0·φ in S, eps0 elsewhere
```
- At **HU = 0**, φ = 1 → `eps = eps0` in *every* phase → the model is numerically identical to
  published v42. (This is why all baseline invariants are preserved automatically.)
- At **HU = 1.0** (the 10 µM anchor), φ ≈ 0.5 *inside the CycA window only* → the clock there
  runs at half speed → the cell **dwells ~2.5× longer in S/G2**, while G1 (where sphase≈0, eps
  unchanged) is untouched.

**The data prove the targeting works:**

| Quantity | DMSO | HU = 1.0 | change |
|---|---|---|---|
| Period | 19.6 h | 27.0 h | +38% |
| **G1 (CycA-low, CycB-low)** | **14.8 h** | **14.4 h** | **≈ unchanged** ✅ |
| **CycA/S–G2 window** | **4.8 h** | **12.6 h** | **+160%** ✅ |
| Divisions / 250 h | 13 | 10 | graded, **no arrest** |
| S-phase EZH2 (vs DMSO S) | 1.00× | **1.31×** | matches experiment |
| CycA peak height | 1.15 | 1.15 | **unchanged** (no amplitude artifact) |

So **the entire 7.4 h of added cycle time is added to S/G2**, exactly as an intra-S checkpoint
should — and because EZH2 is synthesized throughout that now-longer window, its S-phase level
rises **1.31×**, emergently (there is *no* HU term in the EZH2 equation). The φ_min floor (0.25)
guarantees the clock only ever *slows*, never stops, so cells keep dividing at every dose.

**Falsifiable prediction.** In this purely-emergent picture, the measured 1.31× EZH2 elevation
is *tied to* a ~2.5–3× lengthening of S-phase at 10 µM HU. If experimental S-phase/doubling-time
data at that dose disagree, φ_min is the knob to re-fit (lower φ_min → stronger slowing & larger
boost; higher → milder).

---

## 7. Quick reference — who drives whom

```
SHH/MYCN ──(+)──▶ CycD1 mRNA ──▶ Cd ──▶ Md(CycD1-CDK4/6) ──(+)──▶ pRB→pRBp
                     ▲                                                  │
            EZH2 ──(−)┘                                                 ▼
                                                  Md+Me ──▶ pRBp→pRBpp ──▶ E2F released
p27 ──(−)──▶ CDK2 complexes        Skp2 ──(−)──▶ p27          │ (RESTRICTION POINT, +fb)
                                      ▲                        ▼
                                    E2F ──(+)──▶ CycE(Me) ──(+)──▶ Cdc25E(Pe) ──(+)──▶ Me
                                      │
                                      ├──(+)──▶ CycA(Ma) ──(+)──▶ Cdc25A(Pa) ──(+)──▶ Ma
                                      │            │
                                      │            ├──(−)──▶ E2F (phosphorylates → off)
                                      │            ├──(−)──▶ Cdh1a (lets CycB rise)
                                      │            └──▶ EZH2 transcription (with Me) ─┐
                                      ▼                                               │
                                   CycB(Cb/Mb) ──(+)──▶ Cdc25B(Pb) ──(+)──▶ Mb        │
                                          │                                           │
                                          └──(+)──▶ Cdc20a ──(−)──▶ CycA, CycB ──▶ EXIT
                                                                                      │
                                   EZH2 protein ◀───────────────────────────────────┘
                                          └──(−)──▶ CycD1 mRNA   (closes the feedback)

  OUR CHECKPOINT:  HU ──▶ φ(HU)  ──▶ slows `eps` while CycA(Ma) is high  ──▶ longer S/G2
```

---

### Files
- Core model: `src/build_model_v42_mycn.py` (unmodified, published v42)
- v43 extension (rewire + φ checkpoint): `src/build_model_v43_checkpoint.py`
- Validation (12/12): `simulations/validate_v43.py`
- Evidence for §3a/§6: `simulations/diag_phi_levers.py`, `diag_phi_clock.py`, `diag_phi_calib.py`
