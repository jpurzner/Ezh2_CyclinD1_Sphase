# Hedgehog Signaling as a Mitogenic Input: A Modeling-Ready Quantitative Review for a GNP / SHH-Medulloblastoma Cell-Cycle ODE Model

## TL;DR
- **Model the Shh input to cycling GNPs as a slowly-integrating, low-pass-filtered, cell-cycle-gated (duty-cycled) mitogen — NOT as a frequency-encoded pulse train.** Unlike ERK (which activates in discrete asynchronous pulses toggling on/off ~10–20 times per cell cycle, mean burst ~27 min; Albeck et al. 2013) or p53 (fixed-amplitude pulses of period ~5.5 h; Lahav et al. 2004), there is no evidence of intrinsic sub-hour oscillatory Gli dynamics; the pathway's characteristic timescale is the ~20 h "temporal adaptation" of Gli activity seen in the neural tube, gated cyclically by cilium assembly/resorption.
- **The single most important structural feature to include is cilium-dependent gating**: ~94% of outer-EGL GNP progenitors are ciliated (ARL13B+; Ott et al. 2024) but cilia resorb before/at mitosis, so Shh reception is physically interrupted once per cell cycle — a duty cycle, not a continuous input. Adaptation itself (Ptch1 upregulation + Gli2 downregulation + differential Gli-isoform stability) has been characterized in neural tube and NOT in NIH3T3 fibroblasts, and **has never been quantitatively characterized in GNPs specifically — this is the central data gap.**
- **Best-constrained parameters**: Gli1 protein t½ ≈ 40 min; Smo protein/ciliary t½ ≈ 2 h; Gpr161 ciliary exit t½ ≈ 30 min; Gli3FL t½ <30 min (−SuFu) vs >8 h (+SuFu); adaptation timescale ~20 h; SAG EC50 ≈ 3 nM (Chen et al. 2002). **Poorly constrained / must be fitted**: Ptch1 protein t½, all Gli/Ptch1 mRNA half-lives, SuFu turnover, and every GNP-specific dynamic parameter.

## Key Findings

1. **Temporal dynamics, not oscillation.** Under sustained Shh, Gli activity rises to a peak then adapts downward over ~20 h in the neural tube. This is a single slow adaptive pulse, not a repetitive oscillation. The mechanism is threefold (Ptch1 feedback, Gli2 transcriptional downregulation, differential GliA/GliR stability) and is tissue-specific — absent in NIH3T3 fibroblasts.
2. **Hh is a low-pass integrator.** Its effective bandwidth is set by slow protein/ciliary turnover (hours) and the ~20 h adaptation, so it filters out high-frequency ligand fluctuations. This contrasts sharply with the frequency-modulated ERK system.
3. **Cilium gating imposes a cell-cycle duty cycle** on Shh reception, directly coupling the mitogen input to the cell-cycle machinery the user is modeling.
4. **Single-cell Hh response is noisy**; sharp fate/expression domains emerge by downstream GRN filtering, temporal averaging, and population-level cell sorting rather than by precise single-cell readout.
5. **GNP-specific quantitative dynamic data are sparse.** Dose-response and target-gene identity are known, but adaptation kinetics, per-cell noise, and ligand-source dynamics in the EGL are essentially uncharacterized.

## Details

### 1. Temporal adaptation kinetics of Hh/Gli signaling

**Core phenomenon (Dessaud et al. 2007, Nature 450:717; chick/mouse neural [ex vivo] explants and intermediate neural plate cells).** Cells exposed to constant Shh generate an initial burst of Gli activity that then declines; the rate of decline is inversely proportional to Shh concentration, so that ligand concentration is converted into a *duration* of Gli activity. Short exposures to high vs low Shh give similar initial Gli activity and gene expression; divergence appears only over time as cells desensitize. This is the "temporal adaptation" model.

**Quantitative time course (Cohen et al. 2015, Nat Commun 6:6709; mouse neural tube in vivo, GBS-GFP reporter).** Gli activity (GBS-GFP) and Ptch1 protein "initially increase, reaching peak levels early in development (~E9 in mouse), and subsequently decrease to lower levels over time." The adaptation occurs "over a period of ~20 h." Meanwhile the Shh ligand gradient *amplitude increased more than 10-fold* over E8.5–E10.5, while its exponential decay length stayed roughly constant at **19.6 ± 4.2 µm**. So the downstream signal adapts *despite* rising ligand — the hallmark of an adaptive/incoherent system.

**Ptch1 feedback is necessary but not sufficient.** Cohen et al. showed adaptation persists when the pathway is driven downstream of Ptch1 with purmorphamine (10 µM), proving Ptch1 upregulation alone cannot explain it. Their ABC/ODE computational screen identified three mechanisms — Ptch1 transcriptional upregulation, Gli2 transcriptional downregulation, and differential stability of active vs inactive Gli isoforms (active GliA modeled as significantly less stable) — each individually sufficient to reproduce adaptation.

**Gli2 downregulation kinetics (Cohen 2015; mouse neural tube).** Gli2 protein begins decreasing at ~E9 (~15-somite stage), first ventral to the Olig2 domain, then across the entire neural tube by ~25-somite stage. Protein downregulation in dorsal regions precedes mRNA disappearance, implicating both transcriptional and post-translational control.

**Tissue specificity (Cohen 2015; NIH3T3 fibroblasts).** Exposed to 16 nM Shh or 1 µM purmorphamine for 6–48 h, NIH3T3 cells show monotonically *increasing* Gli1 and Ptch1 mRNA — no adaptation — because Shh *upregulates* Gli2 (and Gli3) transcription in these cells, opposite to the neural tube. The mechanistic explanation: the sign of Gli2 transcriptional regulation determines whether the system adapts. This is critical for the user: **adaptation is not a fixed pathway property; it depends on cell-type-specific Gli2 transcriptional wiring.**

**Related dynamics.** Balaskas et al. 2012 (Cell 148:273; mouse neural tube + ODE GRN model) showed the Pax6/Olig2/Nkx2.2 circuit is a bistable, hysteretic network that converts dynamic Gli input into sharp switch-like domain boundaries. Ribes et al. 2010 (Genes Dev 24:1186) and a Shh/Gli three-node timer (Sci Adv 2020, aba8196; chick neural progenitors) documented progressive decay of Gli1/Gli2A/Gli3A activator proteins between ~2.5 and 5.5 "days-different-in-culture," with Gli1 peaking earliest (~3.5DDC).

**GNP gap (explicit):** No study has measured Gli-activity rise/decay kinetics, time-to-peak, or percent adaptation in cerebellar GNPs. GNPs proliferate under sustained Purkinje-derived Shh for days, which is *prima facie* inconsistent with strong neural-tube-style adaptation, suggesting GNPs may resemble the non-adapting (Gli2-upregulating) regime — but this has not been tested. **This is the key missing measurement.**

### 2. Protein and mRNA half-lives / turnover (core parameter table inputs)

| Parameter | Value | System | Source |
|---|---|---|---|
| Gli1 protein t½ | **~40 min** (proteasomal) | *Xenopus* egg extract + mouse NIH3T3 | Huntzicker et al. 2006, Genes Dev 20:276 |
| Gli1 N-degron (aa1–208)-GFP t½ | 180 min | NIH3T3 | Huntzicker 2006 |
| Gli2 (βTrCP-resistant S662A) t½ | **>9 h** (WT shorter, not quantified) | HeLa | Bhatia et al. 2006, J Biol Chem |
| Gli3FL t½ (−SuFu) | **<30 min** | SuFu⁻/⁻ MEF | Humke et al. 2010, Genes Dev 24:670 |
| Gli3FL t½ (+SuFu) | **>8 h** | MEF | Humke 2010 |
| Gli3FL (+Hh) | degrades to ~30% plateau by 6 h | S12/NIH3T3 | Wen et al. 2010, MCB 30:1910 |
| Gli3R | stable ≥6 h | S12 | Wen 2010 |
| SuFu protein t½ | **not measured (GAP)** | — | — |
| Ptch1 protein t½ | **not directly measured**; model-assumed ~7.7 min | — | Lai 2004 (analogy, French & Lauffenburger 1996) |
| Smo protein (total & ciliary) t½ | **~2 h** | smo⁻/⁻ MEF + YFP-Smo | Milenkovic et al. 2009, JCB 187:365 |
| Smo t½ (pulse-chase) | 2.2 h | cultured cells | RanBPM stabilization work |
| Gpr161 ciliary exit t½ | **~30 min** on Shh | mouse cilia | Mukhopadhyay et al. 2013, Cell 152:210 |
| Cyclin D1 (Ccnd1) mRNA t½ | **~30 min** (full-length) vs >3 h (3′UTR-truncated) | Z-138 mantle-cell lymphoma | Wiestner et al. 2007, Blood 109:4599 |
| Gli1/Ptch1/Gli2 mRNA t½ | **not measured in Hh literature (GAP)** — extract from genome-wide sets | — | Sharova 2009 (mouse ES); Schwanhäusser 2011 (NIH3T3) |

**Note on Cyclin D1 mRNA:** the best-documented primary figure is ~30 min for full-length CCND1 (Wiestner et al. 2007: "the half life of the truncated cyclin D1 mRNA (>3 hours) is considerably longer than that of the full-length mRNA (∼30 minutes)"). A frequently-quoted ~15 min value is not traceable to a primary peer-reviewed cyclin-D1 turnover study; use ~30 min.

**Gli processing & degradation machinery.** GliFL→GliR processing is PKA/GSK3/CK1-primed and proteasome-dependent (partial degradation to the repressor). Hh-activated GliA is destabilized by SPOP/Cul3, which competes with SuFu; βTrCP mediates Gli2/Gli3FL degradation. SuFu is the key stability switch: it protects Gli3FL from rapid degradation (t½ <30 min → >8 h; Humke 2010). These give the modeler a mechanistic basis for making GliA less stable than GliFL/GliR (as Cohen 2015 inferred).

**Smo ciliary kinetics.** After Shh/SAG, Smo accumulates in cilia over hours: ~1 h for initial appearance (with reciprocal Ptch1 exit; Rohatgi et al. 2007), reaching near-100% Smo+ cilia by ~6 h (Wilson et al. 2009, PLoS ONE 4:e5182, mouse MEFs). Single-molecule tracking (Milenkovic/Scott/Moerner 2015, PNAS 112:8320; mouse fibroblast cilia, 30 nm/10 ms precision) showed Smo moves primarily by diffusion (D on the order of ~0.1 µm²/s), interrupted by transient binding at sites near the ciliary base; pathway activation (Shh or Ptch1 loss) *reduces* the affinity and frequency of this basal binding, releasing Smo to accumulate. Signal-dependent accumulation is achieved largely by *suppressing constitutive BBSome-dependent exit*.

### 3. Pulsatility, oscillation, and stochastic dynamics

**Hh: no evidence of intrinsic pulsing.** No live-cell Gli biosensor study has reported sub-hour oscillatory or bursty Gli dynamics. Available live-cell work (Xiong et al. 2024, eLife RP96980, zebrafish ptch2:Kaede reporter, ~6 min imaging intervals over ~10 h) shows monotonic, heterogeneous rise trajectories — variable in rate and amplitude between cells — not oscillations. The dominant temporal feature is the slow (~20 h) adaptive pulse. Hh is therefore best characterized as a **slow integrator / low-pass filter**, with an effective filtering timescale set by hours-scale protein and ciliary turnover and the ~20 h adaptation envelope.

**Contrast — ERK/MAPK (frequency-encoded).** Albeck et al. 2013 (Mol Cell 49:249; MCF-10A, EKAR-EV FRET reporter): "ERK is activated in discrete, asynchronous pulses with frequency and duration determined by extracellular concentrations of EGF." At low EGF, isolated bursts last **20–30 min (mean 27 min)** separated by dormant periods of **1–4 h**. Increasing EGF increases pulse *frequency* and duration; individual bursts "becom[e] indistinguishable at higher concentrations of EGF (>200 pg/ml)," with ERK toggling on/off "at least 10–20 times during a single cell cycle." EGFR inhibition reduces pulse frequency; MEK inhibition reduces amplitude. Pulse frequency quantitatively tracks proliferation. This is genuine frequency modulation over minutes-to-hours — the opposite of Hh.

**Contrast — p53 and NF-κB (dynamically encoded).** p53 shows "pulses of fixed amplitude and duration with a period of ~5.5 h" in response to DNA double-strand breaks (Lahav et al. 2004; MCF7), with the *number* of pulses (not amplitude) increasing with damage dose and large cell-to-cell variability in timing. NF-κB shows oscillatory nuclear translocation under TNFα vs sustained under LPS (Covert/Tay). These are the canonical "information in dynamics" systems; Hh conspicuously lacks such fast dynamic encoding.

**Conceptual conclusion for the model:** Hh should be a low-pass/slow-integrator term. It does not carry information in pulse frequency; it integrates ligand exposure over hours and encodes *cumulative duration × level* into Gli activity. Any high-frequency stochasticity in ligand should be strongly attenuated by the pathway.

### 4. Single-cell noise and heterogeneity

**Xiong et al. 2013 (Cell 153:550; zebrafish neural tube, in toto imaging, >200 tracked cells).** Early Shh response is heterogeneous and "fuzzy": progenitors of different fates are initially spatially intermixed with widely overlapping gene-expression domains. Sharp domains form *afterward* by **cell sorting**, which corrects the imprecision of noisy inductive signaling. Ectopically induced motor-neuron progenitors robustly sort to correct positions.

**Xiong et al. 2024 (eLife RP96980).** Single-cell Shh-response trajectories show "significant heterogeneity" and only partial correlation with final fate; motor-neuron progenitors (pMN) are especially variable, and in the posterior neural tube their response dynamics overlap those of more-ventral lateral-floor-plate cells. The authors frame this as a **precision limit** of morphogen-interpretation GRNs in small, dynamic cell fields, arguing precision emerges at the population (not single-cell) level.

**Noise sources:** intrinsic transcriptional noise (bursting), plus extrinsic — positional noise from cell movement, gradient fluctuation, variable ciliation state, and cell-cycle phase. Xiong's work highlights positional noise (cell movement) as a major extrinsic contributor requiring sorting-based correction (adhesion code, Tsai et al. 2020, Science 370:113).

**Morphogen precision / positional information (Zagorski et al. 2017, Science 356:1379; mouse neural tube).** Combined readout of antiparallel Shh (ventral) and BMP (dorsal) gradients achieves domain-boundary positional precision on the order of **~1 cell diameter**. Zagorski et al. concluded single gradients are too noisy beyond ~30 h and that opposing gradients are needed. This was contested: **Vetter & Iber 2022 (Nat Commun 13:1145) and the 2024 exchange (Nat Commun 15:929 / 45148 & 45149)** argue the positional error was overestimated (partly due to 8-bit imaging limits and exponential-fit artifacts) and that a *single* gradient can in principle achieve the observed ~1-cell-diameter boundary precision. **Report the discrepancy, don't pick a side:** positional precision is ~1 cell diameter; whether one or two gradients are required to explain it is unresolved. Information-theoretic (Gregor/Tkačik "positional information in bits") framing has been applied primarily to Drosophila Bicoid (Dubuis et al. 2013, PNAS 110:16301); a rigorous bits-of-information value specific to the Hh gradient is not established.

### 5. Cilium-mediated gating and cell-cycle coupling (high priority)

**Cilium/cell-cycle coupling (general).** Primary cilia assemble in G1/G0 and resorb on cell-cycle re-entry; resorption occurs in two waves — one at G1/S after serum stimulation and a robust wave at G2/M (Tucker; Pugacheva et al. 2007, Cell — Aurora A/HDAC6 axis). The cilium is incompatible with the mitotic spindle (shared centriole/basal body), enforcing resorption before mitosis. Cilia disassembly is triggered by Aurora A → HDAC6 → axonemal microtubule destabilization; Nek2–Kif24 and Tctex-1 (phospho-T94) accelerate resorption and shorten G1 (Li et al. 2011). Forced ciliary resorption shortens G1 and hastens S-phase entry (Tctex-1 phosphomimetic converts hours→minutes). Nde1 (Kim et al. 2011) and CPAP/CENPJ (Gabriel et al. 2016) negatively regulate ciliary length/cell-cycle coupling. Resorption rates in Chlamydomonas ~0.3–0.7 µm/min; mammalian disassembly can also occur by whole-cilium shedding/decapitation (Phua et al. 2017, Cell — decapitation precedes resorption and can itself be mitogenic).

**GNP-specific ciliation (Ott et al. 2024, PNAS 121:e2408083121 / JCB 223:e202404038; mouse cerebellum, volume EM + immunostaining).** In the outer EGL, **~94% of GNP progenitors are ciliated (ARL13B+)**; this "decreased to 75% in the P27KIP1 positive differentiating cells in the inner EGL," then to **12% in the molecular layer and 10% in the IGL**. GNP progenitor cilia average **~1 µm** long. GNP progenitors resorb cilia pre-mitotically (after S-phase, in vivo and in vitro); differentiating GCs additionally undergo a distinct "cilia deconstruction" (internalized/concealed cilia, centriole capping) that permanently removes Hh responsiveness. Aurka is expressed highly in G2/M progenitors but not differentiating cells.

**GNP ciliogenesis control.** Atoh1 controls primary-cilium formation in GNPs to permit Shh-triggered proliferation (Chang et al. 2019, Dev Cell 48:184; via Cep131/centriolar satellites) — Atoh1 loss impairs ciliogenesis and Hh response. TTBK2 (Tau-tubulin kinase 2) is highly expressed in proliferating GNPs under Shh and *stabilizes* cilia (inhibits disassembly), promoting proliferation; the E3 ligase HUWE1 degrades TTBK2 to trigger disassembly and differentiation (PMC11445238 / Cell Death Differ 2024). TTBK2 depletion inhibits SHH-MB proliferation — a therapeutic node.

**Implication (duty cycle).** Because ~94% of progenitors are ciliated but cilia resorb before mitosis, **Shh reception is phase-gated: available through G1/S, extinguished around G2/M**, restored after division. The duty cycle is high (most of the cycle is ciliated) but non-zero-fraction interrupted. For the model this argues for a **cell-cycle-gated Shh input** — a multiplicative gating function on the mitogen term that drops during mitosis and depends on ciliary state, rather than a constant.

### 6. Hh as a mitogen in cerebellar GNPs

**Ligand source and system.** Purkinje cells secrete Shh; GNPs in the EGL proliferate in response (Wechsler-Reya & Scott 1999). Major murine GNP expansion is P0–P14, peaking ~P4–P8, completed ~P15–P21. In vitro, GNPs require insulin/IGF for survival; Shh (or SAG) permits continued cycling. Shh acts as a **paracrine** mitogen (Purkinje→GNP); autocrine contribution is minor.

**Dose-response.** SAG EC50 ≈ **3 nM** (Chen et al. 2002, PNAS 99:14071, Shh-LIGHT2 reporter: "an EC50 of ≈3 nM"), with apparent KD ≈ 59 nM ("significantly higher than the SAG EC50 for pathway activation"); response is biphasic (agonist at low nM, inhibitory >1 µM). In GNP cultures, 5 nM SAG or Shh-conditioned medium sustains proliferation and BrdU incorporation for up to 96 h (Nakashima et al. 2015, Dev Dynamics). A precise EC50 for Shh on *GNP proliferation* and cell-cycle-length-vs-dose curves are not well quantified.

**Downstream targets.** Shh→Gli induces N-Myc (Mycn), Cyclin D1 (Ccnd1) and Cyclin D2 in GNPs (Kenney & Rowitch 2000, MCB 20:9055). N-Myc is necessary and sufficient for GNP proliferation and drives D-cyclins. Cyclin D1 both accelerates the cycle and stabilizes ATOH1 (via CDK4-mediated S309 phosphorylation blocking degradative S328 phosphorylation), coupling proliferation to the immature state (EMBO J 2021, 40:e105712). Precise time-to-induction and fold-changes are only semi-quantified.

**Signal persistence / memory.** Withdrawal of Shh causes **irreversible** growth arrest and differentiation of GNPs (cell-cycle exit). Even *with* continued Shh, GNPs exit the cycle on a cell-intrinsic ~1-week (in vitro) schedule — an intrinsic timer overrides sustained mitogen. This is a crucial modeling constraint: the mitogen is necessary but a cell-autonomous differentiation program imposes a finite proliferative window.

### 7. Noise buffering and robustness mechanisms

**Bistability/hysteresis (Lai, Robertson & Schaffer 2004, Biophys J 86:2748).** The Shh network = positive Gli autoregulation nested in negative Ptch feedback → bistable switch with hysteresis; flips fate at threshold Shh. Stochastic (Gillespie) analysis: noise-induced flipping near the threshold, first-passage time ~30 h (high Shh) to >300 h (low Shh). Most sensitive parameter: Gli promoter strength (v_max,G) — 10% change → >1000% change in switching threshold. Parameter set (Table 1): K1=K2=8.3×10⁻¹⁰ M; k_deg(all Gli)=0.009 min⁻¹ (t½≈77 min); K_Shh=0.58–2.0 nM; K_Ptc=8.3×10⁻¹¹ M; k_deg,P=0.09 min⁻¹; v_max,G=0.24 nM/min; k_g3rc(Gli3→Gli3R)=0.012 min⁻¹; v_max,P=7.5×10⁻¹⁰ M/min. Not in BioModels.

**GRN filtering.** The downstream GRN (Pax6/Olig2/Nkx2.2 in neural tube; Rb/E2F in cycling cells) is itself bistable/hysteretic (Balaskas 2012) and thresholds/filters noisy upstream Gli input — the analog of the user's Rb/E2F switch.

**Temporal averaging.** Slow mRNA/protein turnover integrates and averages transcriptional-burst noise: noise attenuation scales with the ratio of averaging (half-life) time to the noise correlation time — longer half-lives → lower relative variance (standard telegraph-model/Paulsson result). This is directly relevant to the user's bursting (k_on) module: the hours-scale Hh readout half-lives buffer fast burst noise.

**Population-level error correction.** Cell sorting via an adhesion code (Tsai et al. 2020, Science 370:113; zebrafish), apoptosis of mis-specified cells, and opposing BMP gradients all correct single-cell imprecision at the tissue level (Xiong 2013).

### 8. Synthesis and modeling recommendations

**Recommended functional form for the Shh input term in a GNP cell-cycle ODE model:**

The evidence supports a **composite input**, not a single simple form:

`Shh_effective(t) = L(t) · G_cilium(cell-cycle phase) · A(t)`

where:
- **L(t) = tonic/slowly-varying ligand** (paracrine, Purkinje-derived). Treat as approximately **constant / slowly ramping** at the population level; there is no evidence for fast pulsatile ligand delivery in the EGL (though this is unmeasured — flag as assumption).
- **G_cilium = cell-cycle-gated multiplier (duty cycle).** Set ≈1 during G1/S (ciliated), →0 through G2/M (cilium resorbed), restored after division. This is the highest-value structural addition and directly couples the mitogen to the Rb/E2F module.
- **A(t) = adaptation/low-pass filter.** Implement as Ptch1 negative feedback + optional GliA-destabilization, giving a first-order adaptation with **τ ≈ 20 h** IF GNPs adapt like neural tube. BUT: because GNPs proliferate for days under sustained Shh and may resemble the non-adapting NIH3T3 (Gli2-upregulating) regime, **treat the adaptation strength as a free parameter and run sensitivity analysis**; the sign of Gli2 transcriptional regulation is the switch.

**Do NOT** implement Shh as a frequency-encoded pulse train (that is ERK's mode, not Hh's). Do treat the pathway as a low-pass filter of ligand fluctuations (cutoff timescale ~hours).

**Existing models to adapt/couple:**
- **Lai, Robertson & Schaffer 2004** (Biophys J 86:2748): full ODE bistable-switch with published parameter table (above); ideal skeleton for the Gli/Ptch core. Not in BioModels — transcribe from Table 1.
- **Saha & Schaffer 2006** (Development 133:889): spatial extension of Lai; useful if gradient/transport matters. Not in BioModels.
- **Cohen et al. 2015** (Nat Commun 6:6709): ABC-fitted adaptation ODE with three mechanisms; posteriors in supplement; open-access. Best source for adaptation structure.
- **Balaskas et al. 2012** (Cell 148:273): downstream GRN ODE — couple to Rb/E2F analog.
- **BioModels: BIOMD0000000395** (Sivakumar et al. 2011, OMICS 15:729) — curated Hedgehog SBML with parameter set; the main downloadable Hh model.

**Well-constrained parameters:** Gli1 protein t½ (~40 min); Smo protein/ciliary t½ (~2 h); Smo ciliary accumulation (~1 h onset, ~6 h to saturate); Gpr161 exit t½ (~30 min); Gli3FL SuFu-dependent stability (<30 min vs >8 h); adaptation timescale (~20 h, neural tube); GNP ciliation fractions (94/75/12/10%); SAG EC50 (~3 nM); Lai 2004 full parameter set; ERK pulse numbers (for contrast).

**Poorly constrained / must be fitted or sensitivity-tested:**
- Ptch1 protein half-life (never directly measured; Lai's ~7.7 min is a generic-receptor analogy).
- All Gli1/Ptch1/Gli2 mRNA half-lives (extract approximate values from genome-wide sets: Sharova 2009 mouse ES, Schwanhäusser 2011 NIH3T3).
- SuFu protein turnover (unmeasured).
- Whether/how strongly GNPs adapt (the central gap).
- GNP-specific single-cell CV/Fano factors for Gli1/Ptch1 (unmeasured — only zebrafish neural tube heterogeneity exists qualitatively).
- Temporal dynamics of Purkinje Shh secretion in the EGL (assumed constant; unmeasured).
- Cyclin D1/N-Myc induction time constants and fold-changes downstream of Shh in GNPs (only semi-quantified).

## Recommendations

1. **Build the Gli/Ptch core from Lai 2004** (transcribe Table 1 parameters), then add cilium gating and adaptation. Use Gli1 t½≈40 min and Smo t½≈2 h as anchored values; set Ptch1 protein t½ and all mRNA half-lives as free parameters with priors from genome-wide datasets.
2. **Implement cilium-dependent cell-cycle gating** of the Shh input (G_cilium multiplier tied to your Rb/E2F phase variable, dropping through G2/M). This is the most defensible mechanistic coupling to the cell cycle and is GNP-supported (94% ciliated progenitors; pre-mitotic resorption).
3. **Parameterize adaptation as tunable.** Default τ≈20 h (neural tube), but scan adaptation strength from strong (neural-tube-like) to zero (NIH3T3/GNP-proliferative-like). Report model behavior across this axis — it is the largest structural uncertainty.
4. **Treat Shh ligand as tonic/constant** at baseline; add low-frequency noise and verify the pathway low-pass-filters it (a validation that your parameterization reproduces Hh's integrator behavior).
5. **Use half-life-dependent temporal averaging** as your noise-buffering mechanism in the bursting module: longer readout half-lives attenuate k_on burst noise. Couple upstream Gli noise into the Rb/E2F bistable switch and test noise-induced flipping (Lai 2004 predicts this near threshold).
6. **Benchmarks that would change the recommendations:** (a) If a GNP Gli live-reporter shows adaptation with measurable τ → fix the adaptation term. (b) If GNP ligand is shown to fluctuate on <1 h timescales → revisit the low-pass assumption. (c) If GNP single-cell Gli CV is measured → replace assumed noise amplitude. (d) If cilium duty cycle is quantified per phase → replace the gating function with data.

## Caveats
- **The single biggest gap: adaptation has never been quantified in GNPs.** All ~20 h adaptation kinetics come from neural tube; transfer to GNPs is an assumption, and GNP biology (days-long proliferation under sustained Shh) hints they may not adapt strongly.
- Many "half-lives" come from heterologous systems (Xenopus extract, HeLa, MEF, NIH3T3), not GNPs; transferability is uncertain and flagged per row.
- Ptch1 protein t½, SuFu t½, and Hh-target mRNA half-lives are effectively unmeasured — do not present fitted values as data.
- The Zagorski vs Vetter/Iber dispute over morphogen positional precision (~1 cell diameter; one vs two gradients) is unresolved; both are cited.
- No bits-of-positional-information value is rigorously established for the Hh gradient specifically (only for Drosophila Bicoid).
- Smo single-molecule diffusion coefficient (~0.1 µm²/s order) should be confirmed directly from Milenkovic 2015 before use.
- Cohen 2015 model parameters are ABC *posteriors over wide priors* (spanning ~4 orders of magnitude), not tight measurements — use as structural guidance, not precise constants.
- The Gli2 S662A ">9 h" figure was retrieved via a secondary restatement of Bhatia et al. 2006 J Biol Chem — verify against the primary paper before use.