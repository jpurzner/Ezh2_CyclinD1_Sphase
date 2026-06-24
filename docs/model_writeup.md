---
title: "The EZH2–CyclinD1 Cell-Cycle Model (v44)"
subtitle: "A detailed walkthrough — architecture, equations, and figures"
date: "Chahin et al. — GNP / SHH-medulloblastoma model"
geometry: margin=1in
fontsize: 11pt
colorlinks: true
---

# 1. What the model is and why

This is a real-time (minutes) ordinary-differential-equation (ODE) model of the cell cycle of
**cerebellar granule-neuron progenitors (GNPs)** and **SHH-subtype medulloblastoma (MB)**. Its purpose
is to explain, mechanistically, how a single molecule — **Cyclin D1** — integrates two opposing
inputs: a **mitogenic signal** (Sonic Hedgehog → Gli, plus MYCN in tumour) that drives it up, and an
**epigenetic repressor** (EZH2 / H3K27me3) that holds it down. The balance decides whether a cell keeps
cycling or exits into quiescence/differentiation, and it predicts the drug responses (vismodegib
arrests MB; an EZH2 inhibitor rescues that arrest).

The model has **two coupled halves**:

1. **A signalling → Cyclin D1 cascade** (the part we built): Hedgehog and MYCN drive *Ccnd1*
   transcription; EZH2 represses it. This is where the biology of interest lives.
2. **A cell-cycle engine** (the Heldt 2018 model, adopted): once Cyclin D1 is set, a well-characterised
   network of cyclins, CDKs, the Rb–E2F switch, and the mitotic trigger turns it into actual cell-cycle
   progression in real time, with explicit DNA replication and a size/growth checkpoint.

Cyclin D1 is the hand-off point: the cascade *sets* it, the engine *reads* it.

![Model architecture / wiring diagram. The Hedgehog cascade (left) and MYCN converge on Cyclin D1
transcription; EZH2 (induced by the cell cycle) represses it; Cyclin D1 feeds the Rb–E2F restriction
point of the cell-cycle engine (right).](../simulations/fig_v44_wiring.png)

\newpage

# 2. The signalling → Cyclin D1 cascade

## 2.1 Hedgehog pathway: SHH ⊣ Ptch1 ⊣ Smo → Gli

Hedgehog signalling is a **double-inhibition** relay. The receptor **Patched-1 (Ptch1)** tonically
*inhibits* the transducer **Smoothened (Smo)**. The ligand **SHH** binds Ptch1 and relieves that
inhibition, so Smo becomes active and converts the Gli transcription factor from its **repressor** form
(`Gli_rep`) to its **activator** form (`Gli_act`).

Smo activation rate:
$$ v_{\text{Smo}} = \frac{k_{\text{Smo}}}{1 + \dfrac{P_{\text{free}}\,c_{\text{Ptch1}}}{K_{\text{Ptch,Smo}}}}\;(1 - \text{HHi}) $$

- The denominator is the **Ptch1 brake**: free Ptch1 (`P_free`) suppresses Smo. `c_Ptch1` is the
  *functional Ptch1 copy number* — in GNPs it is 1; in SHH-MB, Ptch1 is mutated/lost, so `c_Ptch1` ≈ 0.1,
  which releases the brake and gives the **tonically high** Hedgehog activity that defines the tumour.
- `(1 − GDC)` is **vismodegib (GDC-0449)**, a Smo antagonist: HHi = 1 fully blocks Smo. (Biologically it
  never blocks 100%, which matters numerically — see §5.)

Gli switching (repressor ↔ activator), driven by Smo with a sharp (Hill-2) switch:
$$ \text{Gli}_{\text{rep}} \xrightarrow{\;k\,\frac{\text{Smo}^2}{K^2+\text{Smo}^2}\;} \text{Gli}_{\text{act}} $$

So more Smo ⇒ more `Gli_act`, less `Gli_rep`. The **read-out of pathway activity** is the *Gli1
transcript level* (`Gli_act + Gli1`), which is how Hedgehog activity is measured experimentally — a GNP
at maximal SHH still has far lower Gli1 than an MB.

## 2.2 MYCN

MYCN is itself a Gli target and is amplified in many MBs:
$$ v_{\text{MYCN}} = k_{\text{bas}}\,A_{\text{MYCN}} + k_{\text{Gli}}\,\frac{\text{Gli}_{\text{act}}+\text{Gli}_1}{K+\text{Gli}_{\text{act}}+\text{Gli}_1} $$

`A_MYCN` is the amplification factor (1 in GNP, ≈2.8 in MB). The first term is a **mitogen-independent
floor** — even with Hedgehog fully off, an MB still makes basal MYCN. This floor is central to the drug
biology: vismodegib removes the Gli-driven part of Cyclin D1 but **cannot remove the MYCN floor**, so
vismo is only a *partial* mitogen withdrawal in MB.

## 2.3 Cyclin D1 transcription — the convergence node

This is the heart of the model. The *Ccnd1* transcription rate is the **(basal + Gli + MYCN) drive**
multiplied by an **EZH2 repression factor**:
$$ v_{Ccnd1} = \Big(\underbrace{k_0}_{\text{basal}} + \underbrace{k_{\text{Gli}}\frac{G^{n}}{K_a^{n}+G^{n}}\cdot\frac{K_r^{n}}{K_r^{n}+\text{Gli}_{\text{rep}}^{n}}}_{\text{Gli activation, repressor-gated}} + \underbrace{k_{\text{MYCN}}\frac{\text{MYCN}^{m}}{K_M^{m}+\text{MYCN}^{m}}}_{\text{MYCN drive}}\Big)\cdot \underbrace{R(\text{EZH2})}_{\text{repression}} $$

where $G=\text{Gli}_{\text{act}}+\text{Gli}_1$. In words: Cyclin D1 mRNA is made in proportion to how
much activating Gli and MYCN are present (each a saturating Hill function), reduced by how much Gli
*repressor* is around, and then **scaled down by the epigenetic repressor**. The repression factor in
the production model is
$$ R(\text{EZH2}) = \frac{K_{\text{rep}}}{K_{\text{rep}} + \text{EZH2}\,(1-\text{EZH2i})} $$

- At EZH2 = 0, $R=1$ (no repression). As EZH2 rises, $R$ falls toward 0 — EZH2 dials Cyclin D1 down.
- `EZH2i` is the EZH2 inhibitor (tazemetostat): EZH2i = 1 sets the repression to 1 (de-repressed).
- (In the newer epigenetic module, §5, this single factor is replaced by an explicit, *leaky*
  H3K27me3-mark term — but the structure "drive × repression" is identical.)

Cyclin D1 protein then follows its mRNA essentially instantaneously (it is short-lived, ~1 min):
$$ \frac{d[\text{Cd}]}{dt} = k_{\text{tl}}[\text{Cd\_mRNA}] - k_{\text{deg}}[\text{Cd}] $$
so $[\text{Cd}] \approx (k_{\text{tl}}/k_{\text{deg}})[\text{Cd\_mRNA}]$ — transcript and protein are
interchangeable read-outs.

## 2.4 EZH2 — the cell-cycle-coupled repressor

EZH2 is **not constant**; it is itself induced by proliferation (it is an Rb–E2F target) and by mitogen
dose. Its transcription:
$$ v_{\text{EZH2}} = k_{\text{bas}} + k_{E2f}\,\frac{E2f}{K+E2f}\Big(w\frac{Ce}{K_e+Ce}+(1-w)\frac{Ca}{K_a+Ca}\Big)\frac{Cd}{K_d+Cd} $$

Reading this term by term: EZH2 is made when (i) **E2F** is active (the commitment marker — saturating,
so EZH2 reports that the cell has committed), gated by (ii) the **S/G2 cyclins** Cyclin E (`Ce`) and
Cyclin A (`Ca`) (so EZH2 peaks in S/G2), and scaled by (iii) a **mitogen-dose** term `Cd/(K_d+Cd)`
(so EZH2 climbs with Cyclin D1 / mitogen across a wide range — reproducing the dose-dependent EZH2
seen over recombinant-SHH titrations, and the ≈2× higher EZH2 in MB vs GNP). EZH2 protein is **long-
lived** (degradation `kDeEZ` gives a half-life of days), which becomes important for the de-repression
kinetics (§5).

The result is a **negative feedback loop**: Cyclin D1 → cell cycle → E2F → EZH2 → ⊣ Cyclin D1. EZH2
buffers Cyclin D1 and ties it to the *current* proliferative/mitogenic state.

\newpage

# 3. The cell-cycle engine (Heldt 2018 core)

Once Cyclin D1 is set, the model uses the **Heldt et al. (2018) real-time mammalian cell-cycle model**
to turn it into progression. We adopt this engine rather than re-deriving it; here is what each part
does and the one equation that couples our biology into it.

- **Restriction point (Rb–E2F).** Retinoblastoma protein (Rb) sequesters E2F. Cyclin D1–CDK4/6
  phosphorylates Rb, releasing E2F, which then drives Cyclin E and the G1/S transition. This is a
  **bistable switch** (E2F activates its own activators): below a Cyclin D1 threshold the cell stays
  in G1 (E2F off); above it, the switch flips and the cell commits. The Cyclin-D1-driven Rb
  phosphorylation, with CDK-inhibitor (CKI) modulation, is:
  $$ v_{\text{Rb-P}} = k_{\text{Ph}}\,\frac{[\text{Cd}]}{K_{CdRb}\,(1+p16+p18+w\,p27) + [\text{Cd}]} $$
  The CKIs **raise the threshold**: INK4 proteins (p16, p18) and CIP/KIP (p27/p21) make it take *more*
  Cyclin D1 to phosphorylate Rb. MB has elevated CKI tone, so MB needs more Cyclin D1 to cycle — the
  mechanism by which the same Cyclin D1 that lets a GNP cycle leaves an MB arrested.

- **Growth / size gate.** Heldt's Rb is conserved in total, so its *concentration* dilutes as the cell
  grows; the model's mass-dependent gate $\text{mass}^n/(M^n+\text{mass}^n)$ is algebraically the
  **Rb-dilution restriction point** (Zatulovskiy 2020): a cell must grow to a critical size (dilute Rb
  enough) before it can pass G1/S. This sets the ~22 h period.

- **Skp2–p27 commitment.** A feed-forward loop: once E2F leaks open, Cyclin E/A build Skp2, which
  degrades p27, which removes the brake on Cyclin E/A — an irreversible bistable commitment (the
  molecular "point of no return").

- **S phase.** Explicit DNA replication: licensed origins (`Rc`) fire into active replication
  complexes (`aRc`) that synthesise DNA (`Dna`: 0 → 1) at fork speed; an intra-S CHK1 checkpoint holds
  mitosis until replication completes. (This explicit S phase is what lets the model know *when* the
  *Ccnd1* locus replicates — essential for the dilution module in §5.)

- **Mitotic switch.** Cyclin B–CDK1 (MPF) with the Cdc25 (activating) / Wee1 (inhibitory) double-
  feedback forms a second bistable switch that triggers mitosis; APC/C–Cdc20 then degrades cyclins and
  resets the cell. Division halves the mass and resets the daughter with **high p27** (born braked),
  which produces a transient post-mitotic G0/G1 window.

The single coupling is: **our Cyclin D1 enters Heldt's $v_{\text{Rb-P}}$**, and Heldt's E2F/cyclins
feed back into our EZH2 synthesis. Everything else in the engine is standard.

\newpage

# 4. Core results — figures

## 4.1 Validation against the data

![Validation grid (Supp. Fig. 7): Cyclin B traces, and model-vs-data for Cyclin D1, MYCN, Gli1, EZH2,
division counts and time-courses across conditions. The model passes 22 of 27 quantitative targets.](../simulations/fig_v44_validation.png)

The model reproduces (among others): the Cyclin D1 MB/GNP fold (~5×), EZH2 MB/GNP (~2×), the
dose-dependent EZH2 rise over SHH, Gli1 reductions under vismodegib, ~22 h periods, and flow-cytometry
phase fractions.

## 4.2 The drug logic: vismodegib arrests, EZH2i rescues

![Fig. 5 — Ptch+/− MB cycles; +vismodegib (HHi) arrests; +vismodegib+EZH2i rescues. Cyclin B traces
show division (peaks) vs arrest (flat).](../simulations/fig_v44_fig5_rescue.png)

This is the central tumour-biology result. Vismodegib removes the Gli drive; with the EZH2 brake
intact, the residual (MYCN-floor) Cyclin D1 falls below the proliferation threshold → **arrest**. Add an
EZH2 inhibitor → the brake is released → Cyclin D1 recovers above threshold → **rescue**. CDK4/6
inhibition does *not* rescue (it acts below the lesion), a key specificity control.

![Population view (Fig. 5C): in a heterogeneous tumour, vismodegib drops the cycling fraction and EZH2i
restores a fraction of it; CDK4/6i does not.](../simulations/fig_v44_fig5_population.png)

## 4.3 The EZH2 ⊣ Cyclin D1 feedback

![EZH2→CyclinD1 negative feedback (Supp. Fig. 8): with vs without the feedback arm — Cyclin B, Cyclin D1
mRNA, EZH2 protein, and their phase relationship.](../simulations/fig_v44_feedback.png)

![The feedback installs a mitogen-gated quiescence threshold: with the EZH2⊣CyclinD1 arm, cells exit the
cycle as mitogen falls; without it, they keep cycling.](../simulations/fig_v44_feedback_role.png)

## 4.4 Hedgehog feedback and Gli1 read-out

![Gli→Ptch1 negative feedback and its disruption in MB: Gli1 step overshoot (GNP vs MB), Ptch1 mRNA
induction, and the broken Smo/Gli loop in tumour.](../simulations/fig_v44_ptch1_feedback.png)

\newpage

# 4.5 Mitogen sensitivity and the operating point

![Mitogen sensitivity expressed in Hh-activity (Gli1) units, with the correct knobs: GNP tuned by SHH
(up), MB tuned by vismodegib (down toward the MYCN floor), ± the EZH2 feedback.](../simulations/fig_v44_mitogen_sensitivity.png)

![Conditions mapped onto the Cyclin D1 × p16 parameter plane with the cycle/arrest boundary: MB+HHi
arrests at the *same* Cyclin D1 a GNP cycles at, because the raised CKI threshold; EZH2i crosses it
back.](../simulations/fig_v44_condition_paramspace.png)

# 4.6 Quiescence (transient G0) and its heterogeneity

![One-cycle anatomy: the transient G0 is the p27-high window. (Phospho-Rb saturates early, so we score
G0 by the p27 marker.)](../simulations/fig_v44_cycle_anatomy.png)

![CyclinD1 sets transient-G0 duration (p27 marker); EZH2i collapses it; below a threshold the cell
arrests.](../simulations/fig_v44_g0_mechanism.png)

![Proliferation–quiescence bifurcation from cyclin D1 / **birth-p27** heterogeneity (N=140 cells): low
birth p27 → immediate re-entry, high → prolonged transient G0 (the Fan–Meyer threshold); low Cyclin D1 →
arrest. SHH and EZH2i tune the split.](../simulations/fig_v44_g0_bifurcation.png)

The cell-to-cell decision is set by the **Cyclin D1 / p27 ratio** at birth crossing one sharp threshold
(Spencer/Cappell/Fan–Meyer). Varying inherited p27 and Cyclin D1 setpoint across a population reproduces
the split between immediate re-entry, transient quiescence, and arrest.

\newpage

# 4.7 The EZH2 feedback as a dynamical system (phase-plane analysis)

![Ferrell-style phase-plane / bifurcation analysis of the EZH2 ⊣ CyclinD1 feedback. (A) The two
nullclines meet at a single stable node — the negative feedback is monostable and *homeostatic*,
buffering Cyclin D1 ~2.2×. (B) Mitogen drive and feedback-removal move the operating point. (C/D) The
feedback shifts the proliferation bifurcation to higher drive / higher mitogen.](../simulations/fig_v44_ezh2_phaseplane.png)

This analysis (in the spirit of James Ferrell's phase-plane teaching) proves the EZH2 loop is **not a
bistable switch** — a single nullcline intersection, monostable and buffering. The bistability that
matters (cycle vs arrest) lives in the Cyclin D1→commitment axis; EZH2's role is to *reposition* that
bifurcation point, consistent with the measured (weak/buffering) feedback strength.

# 4.8 Mitogen withdrawal (GNP) and vismodegib (MB) as the same lesion

![Mitogen-withdrawal experiment (heterogeneous GNP ensemble): post-withdrawal divisions vs final mitogen
level, sudden vs gradual, ± EZH2 feedback. The feedback's effect is **maximal at partial withdrawal**
(the feedback-gated window) and ~0 at full withdrawal; sudden > gradual.](../simulations/fig_v44_mitogen_withdrawal.png)

![Vismodegib on the heterogeneous MB ensemble — the MB analogue. Because the MYCN floor remains, even
full vismo is a *partial* withdrawal, so the EZH2 feedback is decisive: feedback ON → arrest;
no-feedback (≡ EZH2i) → the **bimodal rescue** (a fraction of cells coast).](../simulations/fig_v44_vismo_withdrawal.png)

Mitogen withdrawal (loss of SHH) and vismodegib (Smo block) are the **same molecular lesion** — partial
loss of the Gli-driven Cyclin D1 input — read out through the same EZH2 node. The EZH2-inhibitor rescue
is **MB/vismo-specific**: MB sustains high EZH2, so its arrest is EZH2-maintained and reversible; a
mitogen-withdrawn GNP loses EZH2 (it is cell-cycle-coupled), so its arrest is *not* EZH2-maintained.

\newpage

# 5. The H3K27me3 epigenetic layer (current work)

The production model above lumps EZH2's repression into the single factor $R(\text{EZH2})$. The newer
work makes the **mark itself** explicit, because that is where the interesting kinetics and the
replicative-dilution biology live. Two gated variants exist (both default-off; the production model is
unchanged and still validates 22/27).

## 5.1 Why an explicit mark: de-repression has a real timescale

An EZH2 *inhibitor* (tazemetostat) blocks the methyltransferase but does not remove EZH2 protein. So
de-repression of Cyclin D1 is set by how fast the **H3K27me3 mark** is lost (demethylation + dilution),
not by EZH2 turnover. With an explicit mark `H3K27_Cd` that tracks EZH2 at steady state but decays after
EZH2i, the model reproduces a graded de-repression over ~24–48 h, and the rescue can be triggered at any
time.

![EZH2i fast-tracks Cyclin D1 de-repression → rescues the vismo-arrested MB cell. (A) Cyclin D1 ± EZH2i;
(B) the H3K27me3 mark decaying once deposition is blocked; (C) timed dosing at 0/24/48/72 h — de-represses
whenever given; (D) the rescue is available on-demand.](../simulations/fig_v44_ezh2i_rescue_kinetics.png)

## 5.2 The replicative-dilution module (the spec)

This is the most mechanistic version: the H3K27me3 occupancy `Mk ∈ [0,1]` over the ~7 kb *Ccnd1* domain
obeys its own dynamics, and — crucially — the **cell cycle itself dilutes the mark** at each round of
replication, so the cycle period sets the dilution frequency.

**Continuous dynamics between divisions:**
$$ \frac{d\,Mk}{dt} = \big(\,\underbrace{k_w\,\text{EZH2}\,Mk}_{\text{read-write (EED reads me3)}} + \underbrace{k_0}_{\text{de-novo}}\big)(1-EZH2i)\,(1-Mk) \;-\; \underbrace{\delta\,Mk}_{\text{demethylation/turnover}} $$

- **Read-write** $k_w\,\text{EZH2}\,Mk$ is *autocatalytic*: PRC2 (EED subunit) reads existing me3 and
  writes more, so the mark reinforces itself. It is proportional to EZH2 (the writer) and blocked by
  EZH2i.
- **De-novo** $k_0$ is a small allostery-independent floor so that a *halved* locus can re-seed.
- $(1-Mk)$ limits methylation to unmethylated substrate (saturates at full occupancy).
- $\delta\,Mk$ is active demethylation (KDM6) + histone turnover — the reason the mark is sub-saturating.

**Discrete replicative dilution** — at the locus's replication time in early S (once per cycle):
$$ Mk \;\longrightarrow\; \tfrac{1}{2}\,Mk $$
Replication distributes the parental methylated histones over two daughter strands, halving occupancy.
The cell then has the rest of the cycle to restore it. If the cycle is **fast**, restoration cannot keep
up and the mark ratchets down (de-repressed, proliferative); if **slow**, restoration wins (mark stays,
repressed). This is the loop the design closes: *the cell cycle the mark controls also sets how often
the mark is diluted.*

**Leaky repression (a key biological correction).** H3K27me3 impedes transcriptional *elongation* but
RNA Pol II is still present and substantial transcript is made — so the mark **dampens, it does not lock
out**. The repression factor therefore has a residual floor $f_0$:
$$ R(Mk) = f_0 + \frac{1-f_0}{1+\left(Mk/K\right)^{n}} $$
At $Mk=0$, $R=1$ (full transcription); at high mark, $R$ saturates at the residual $f_0$ (here ≈0.23 —
~23 % of transcription survives even at full occupancy). $R$ multiplies the same Gli+MYCN drive, so the
residual scales with the activator input.

This leaky, graded form turned out to be not only more biological but **better-fitting**: a sharp
lockout Hill capped the MB/GNP Cyclin D1 fold at ~4; the graded leaky repressor reaches the target, and
the module passes the full validation at **22/27 — parity with the production model**.

![Sawtooth dynamics of the replicative-dilution mark (MB): `Mk` is halved at each S-phase and restored
over the cycle, crossing the repression threshold K; Cyclin D1 pulses after each dilution. GNP and MB
both cycle with this module.](../simulations/probe_h3k27_dilution.png)

\newpage

# 6. Additional figures (dose-response, knockdown, ramps)

![Graded EZH2i dose-response of Cyclin D1 and proliferation across GNP / MB / MB+HHi / MB+CDK4/6i — the
rescue threshold and the no-rescue control.](../simulations/fig_v44_ezh2i_doseresponse.png)

![Population EZH2i dose-response: fraction of a heterogeneous tumour pushed back into cycle vs dose
("low dose rescues a fraction").](../simulations/fig_v44_ezh2i_population_dose.png)

![Divisions-to-arrest and the cumulative effect of the EZH2 brake.](../simulations/fig_v44_divisions_to_arrest.png)

![Mitogen ramp-down: gradual vs sudden mitogen decline and the cell-cycle response.](../simulations/fig_v44_mitogen_rampdown.png)

\newpage

# 7. Validation status and what remains

**Validation: 22 / 27 quantitative targets** for both the production model and the new dilution module
(same fails — the module adds no regression). The five outstanding fails are pre-existing:

1. **EZH2 transcript S/G0 gradient** (~5.7; target 1.8–2.5) — the within-cycle EZH2 gradient is too
   steep. The longest-standing open item.
2. **MB hydroxyurea S- and G2-folds** — the HU fork-coupling needs re-tuning for the short S phase.
3. **MYCN GNP+HHi** ratio (0.91 vs 0.78) — a cascade target.
4. **MB G2+M duration** (~3.6 vs 2.5 h) — an accepted trade made during the last recalibration.

**Open decisions / next steps:** (i) whether the replicative-dilution model is promoted to the default
repression mechanism (it is at validation parity); (ii) running the design-spec phenotype simulations —
the **cycle-length (T_cc) sweep** that locates the critical period where dilution outruns restoration
(the replicative-dilution phenotype), and the **EZH2i × CDK4/6i synergy** (CDK4/6i lengthens the cycle →
fewer dilutions → mark restored → differentiation).

# 8. One-paragraph summary

Cyclin D1 sits at the convergence of a mitogenic drive (Hedgehog→Gli, plus MYCN) and an epigenetic brake
(EZH2/H3K27me3). The drive and the brake are each graded; their product, read through a bistable
restriction point with a growth-timed threshold, decides cycle vs exit. Mitogen withdrawal (GNP) and
vismodegib (MB) are the same lesion — partial loss of the Gli drive onto a residual floor — and EZH2
inhibition rescues specifically the MB arrest because only there is EZH2 sustained. Making the H3K27me3
mark explicit, with replicative dilution and a *leaky* (non-lockout) repression, reproduces all of this
and gives the de-repression its real, on-demand timescale — at full validation parity with the original
model.
