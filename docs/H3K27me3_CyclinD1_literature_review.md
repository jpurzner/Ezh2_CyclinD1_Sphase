# H3K27me3 / Cyclin D1 Model — Literature Review & Evidence Base

*Consolidated summary of the papers reviewed, their main findings in relation to our model, and the extracted timing/kinetics tables. Companion to the model design spec and to `h3k27me3_repression_model_state.md` (the model that implements this evidence base = the mean-field AUM `with_h3k27_dilution` module). Citation details marked † should be verified against PMID/DOI before publication; the underlying finding is well-supported.*

---

## 1. Model context (recap)

We are modeling H3K27me3-mediated repression of *Ccnd1* and its replicative dilution, added as a compact module to our existing cell-cycle model. The defining feature is that the **replication frequency that dilutes the mark is endogenous** — set by the Cyclin D1–driven cycle, which the mark itself controls via repression of *Ccnd1*. The module is one state variable (`M`, fractional H3K27me3 occupancy of the ~7 kb / ~35-nucleosome domain), one early-S halving event, and one repressive factor on the existing *Ccnd1* synthesis term.

Key locus facts anchoring the literature: *Ccnd1* is **bivalent** (broad H3K4me3 over the H3K27me3 domain), **non-CpG-methylated** (early-replicating, CpG-island-like), **large** (read-write sufficient, no structural-memory term needed), and carries **initiation-competent Ser5P Pol II** (poised, elongation-limited).

---

## 2. Key conclusions from the literature (executive summary)

1. **The complex, not the mark, is the proximal repressor.** PRC1 (and complex occupancy generally) throttles transcription at initiation / burst frequency; H3K27me3 is the slower targeting/memory layer. → Our readout should modulate initiation/burst frequency with a leaky floor, not an elongation block.
2. **Repression operates on three separable timescales.** Mark removal (days) ≫ PRC2-complex removal (~12 h) > PRC1 removal (~hours). In the physiological dilution scenario nothing is actively removed, so **mark loss is rate-limiting** and complex/transcription are fast followers — which justifies the instantaneous Hill readout.
3. **The coupling is bidirectional.** Nascent RNA from ongoing transcription inhibits/evicts PRC2 in real time. → A reciprocal transcription→PRC2 arm should gate the methylation rate; the bivalent sub-saturating set-point M_ss emerges from this antagonism rather than being imposed.
4. **PRC2 is dynamic, not parked.** ~20% chromatin-bound, seconds-scale residence; the mark is an allosteric *effector* of catalysis, not a binding stabilizer. → Justifies treating read-write as a catalytic-rate term scaled by local mark, not an occupancy term.
5. **Restoration is slow and density-dependent.** H3K27me3 level recovery takes ~1 day (faster at large dense domains), setting τ_restore ≈ 10–18 h and the dilution-vs-restoration race.
6. **Absolute output is anchorable.** Genome-average ~2 mRNA/h; robust genes ~10–50/h; ceiling ~1 Pol II/s. → v_init anchor for *Ccnd1*.

---

## 3. Findings by theme (papers + model relevance)

### 3.1 Mechanism of repression — mark vs complex
**Finding.** Polycomb represses via poised-Pol II holding, PRC1 compaction, and elongation interference; H3K27me3 marks chromatin for reader/effector recruitment rather than physically blocking Pol II. Acute PRC1 (RING1B) removal derepresses targets *while H3K27me3 persists*, and PRC1 controls transcription by limiting PIC formation / burst frequency.
**Model relevance.** The repressive term acts at **initiation/burst frequency** with a leaky floor η₀ (loaded Ser5P Pol II → residual firing). `M` is a valid proxy because complex occupancy tracks the mark at steady state.
**Papers.** Blackledge & Klose 2021 (Nat Rev Mol Cell Biol); Dobrinić et al. 2021 (NSMB)†; Klose-lab "deep OFF / PIC" 2024†; Finogenova et al. 2020 (eLife).

### 3.2 Two-timescale derepression — PRC1 ≠ PRC2 ≠ mark
**Finding.** Derepression kinetics depend on what is removed: PRC1 degron ~hours, PRC2 (SUZ12) degron ~12 h for a subset, EZH2 catalytic inhibition days (replication-coupled mark dilution).
**Model relevance.** Confirms the **two-timescale architecture**: slow mark (rate-limiting, our dilution event) + fast complex/transcription follower (justifies quasi-steady-state readout). Predicts a falsifiable signature — degron-fast vs EZH2i-slow derepression.
**Papers.** Dobrinić et al. 2021 (NSMB)†; Petracovici & Bonasio 2021 (Mol Cell)†; McCabe et al. 2012 (Nature); Knutson et al. 2012 (Nat Chem Biol)†; Jadhav et al. 2020 (Mol Cell).

### 3.3 Transcription ↔ PRC2 bidirectional coupling — the reciprocal arm
**Finding.** PRC2 continuously samples promoters; at active genes it binds nascent RNA (5′ ezRNAs), which inhibits its HMTase activity and/or evicts it. PRC2–RNA and PRC2–chromatin binding are mutually antagonistic (higher affinity for G-tract RNA than chromatin). RNA–PRC2 interactions also tune Pol II pausing/elongation.
**Model relevance.** Adds a **transcription→PRC2 feedback arm** (double-negative ⇒ positive feedback) on top of read-write. Implemented as a transcription-gated methylation rate `k_w·Z(t)·(1 − g(v_Ccnd1))`. This lets the **bivalent M_ss emerge** as the antagonism balance and sharpens the switch / adds hysteresis. Magnitude is debated → include with tunable strength.
**Papers.** Kaneko et al. 2013, 2014; Davidovich et al. 2013 (NSMB); Cifuentes-Rojas et al. 2014; Beltran et al. 2016 (Genome Res)†; Wang et al. 2017 (Mol Cell); Rosenberg et al. 2021 (NSMB)†; Song et al. 2023 (Science)†. Co-occupancy with poised Ser5P Pol II: Stock et al. 2007; Brookes et al. 2012; Ferrai et al. 2017 (Mol Syst Biol).

### 3.4 H3K27me3 dilution & restoration
**Finding.** Positional inheritance is faithful (histone recycling) but *level* restoration by de novo methylation is slow, locus-specific, spread across the cycle, and faster at dense/repressed domains (H1-facilitated). ~40% residual still silences; poised promoters need more erasure; derepression takes multiple divisions, proliferation-rate dependent.
**Model relevance.** Sets the **dilution event** (½ at early S), **τ_restore ≈ 10–18 h** (fast end, given the large dense *Ccnd1* domain), and the **promoter-specific threshold K** (high, because bivalent/poised). Replication-rate dependence is the core dilution phenotype.
**Papers.** Reverón-Gómez et al. 2018 (Mol Cell); Escobar et al. 2023 (Nat Commun)†; Alabert et al. 2015 (Genes Dev)†; Jadhav et al. 2020 (Mol Cell); Howard-lab hybrid model 2021 (eLife)†.

### 3.5 Absolute transcription rates
**Finding.** Genome-average ~2 mRNA/h; highly transcribed genes pack ~1 Pol II/270 bp; peak initiation ~1 Pol II/s; burst size set by core promoter/Pol II availability, frequency by enhancers/initiation probability; re-initiation in seconds–minutes.
**Model relevance.** Anchors **v_init ≈ 10–50 mRNA/h** for fully-active *Ccnd1*; PRC1 throttles the **frequency** arm (matches our initiation/burst-frequency readout). Short *Ccnd1* mRNA half-life makes output track repression quickly.
**Papers.** Schwanhäusser et al. 2011 (Nature); *Fast transcription rates…* EMBO Rep 2011†; bursting reviews (Mol Cell 2025†, PNAS 2025†).

### 3.6 PRC2 recruitment, read-write, allostery
**Finding.** EED reads H3K27me3 and allosterically activates EZH2 (read-write); accessory subunits (AEBP2, PCL/MTF2) drive chromatin binding; the mark acts as an allosteric effector of catalysis, not a dwell-time stabilizer; PRC2 is ~20% bound / seconds residence.
**Model relevance.** Justifies the **read-write term** `k_w·Z·M·(1−M)` as a catalytic-rate term scaled by existing mark, and the **de novo floor k₀** (accessory-driven nucleation) needed to reseed after halving.
**Papers.** Margueron et al. 2009 (Nature); Youmans et al. 2018 (Genes Dev); Finogenova et al. 2020 (eLife); single-molecule review NAR 2021†; DIPG K27M SMT 2018 (Nat Commun)†.

### 3.7 E2F–EZH2 cell-cycle coupling
**Finding.** EZH2 and EED are E2F targets downstream of pRB-E2F (four E2F sites in the EZH2 promoter); the cell-cycle-coupled EZH2 isoform is consistent with maintaining H3K27me3 during replication.
**Model relevance.** Supplies the **EZH2 capacity waveform Z(t)** that gates methylation, phase-locked to G1/S — the restoration arm of the network. (Already represented in our existing model.)
**Papers.** Bracken et al. 2003 (EMBO J); Chen et al. 2018 (Epigenetics & Chromatin).

### 3.8 Pol II pausing & bivalency (readout context)
**Finding.** PRC2 targets carry poised Ser5P-only Pol II (no Ser2P/Ser7P), initiation-competent but elongation-blocked, persisting through differentiation; H3K4me3 promotes pause-release via INTS11; in cortical NPC→neuron, pause-release is decoupled from activation and H3K27me3 positively correlates with paused Pol II.
**Model relevance.** Decides whether the M→Cyclin D1 readout needs a **transcriptional delay term**. Because the developmental context shows K27me3-coincident paused Pol II, a delay is possible — flagged as checkable via a *Ccnd1* Pol II pausing index in our data. Default: instantaneous Hill (delay ≪ τ_restore).
**Papers.** Brookes et al. 2012; Ferrai et al. 2017 (Mol Syst Biol); Wang et al. 2023 (Nature); Liu et al. 2017 (Cell Reports).

### 3.9 Modeling lineage
**Finding.** Per-nucleosome stochastic models (Dodd) require cooperativity + non-local recruitment for replication-resistant bistability; the Howard A–U–M PRC2 model fits mass-spec H3K27me3 kinetics and shows slow chromatin dynamics filter TF fluctuations; ODE reductions exist.
**Model relevance.** We **lift the A–U–M kernel** (calibrated to real restoration kinetics) and add the endogenous replication clock; Dodd sets the cooperativity/non-locality requirements; mean-field is justified at N≈35.
**Papers.** Dodd et al. 2007 (Cell); Angel et al. 2011 (Nature); Berry, Dean & Howard 2017 (Cell Systems); Howard-lab hybrid model 2021 (eLife)†; Klingel et al. 2021 (FEBS J)†.

---

## 4. Timing & kinetics tables

*Organized by process category; within the model the axis runs fast (seconds) → slow (days).*

### A. Mechanical transit times (seconds–minutes) — justify the instantaneous dilution event
| Process | Timing | Reference | Evidence (method) | Notes |
|---|---|---|---|---|
| Replication fork transit over 10 kb | ~5–7 min single-fork (~2.5–3.5 bidirectional) | Standard fork rate ~1.5–2 kb/min | DNA fiber / combing (textbook) | Halving ≈ instantaneous (≪ T_cc) |
| Pol II elongation over 10 kb | ~2.5–5 min | *Fast transcription rates…* EMBO Rep 2011† | Live-cell MS2 imaging | 1.3–4.3 kb/min |
| Pol II elongation, IEGs | 3.5–4 kb/min | IEG kinetics paper† | Native nascent-transcript time course | Rapidly inducible gene |

### B. PRC2 complex chromatin dynamics (seconds) — the writer is not parked
| Process | Timing | Reference | Evidence (method) | Notes |
|---|---|---|---|---|
| PRC2 chromatin-bound fraction | ~20% bound / ~80% diffusing | Youmans et al. 2018 (Genes Dev) | Endogenous HaloTag-EZH2/SUZ12 + SPT | Continuous sampling |
| PRC2 residence time | Seconds (FRAP recovery tens of sec) | Youmans et al. 2018 | FRAP + SPT at marked locus | Greater immobile fraction at marked loci |
| Mark effect on dwell time | No immediate change when EED–H3K27me3 blocked | Youmans et al. 2018 | A-395 inhibitor + SPT | Mark = allosteric effector, not stabilizer |
| Residence time ↔ function | Correlates with silencing | DIPG K27M SMT 2018 (Nat Commun)† | Live-cell SMT | K27M prolongs EZH2 residence |

### C. Transcription ↔ PRC2 antagonism (real-time, per-transcript) — reciprocal arm
| Process | Timing | Reference | Evidence (method) | Notes |
|---|---|---|---|---|
| Nascent RNA inhibits/evicts PRC2 | Real-time, per-transcript | Kaneko 2013/2014; Davidovich 2013 | EZH2 PAR-CLIP/CLIP; in vitro HMTase | "Sensing" of transcription |
| PRC2–RNA vs PRC2–chromatin | Mutually antagonistic | Beltran et al. 2016 (Genome Res)† | In vitro binding + RIP | Higher affinity for G-tract RNA |
| G-quadruplex inactivation | mechanism (not a rate) | Song et al. 2023 (Science)† | Cryo-EM | Structural basis |
| RNA–PRC2 → Pol II | Blocks productive elongation | Rosenberg et al. 2021 (NSMB)† | CLIP + Pol II-S5P ChIP-seq | Links axis to pausing |

> Caveat: in vivo significance/specificity of PRC2–RNA binding is debated (Brockdorff 2013; RIP artifacts; "junk-mail" model not fully validated).

### D. Transcription output rates — the v_init anchor
| Process | Rate | Reference | Evidence (method) | Notes |
|---|---|---|---|---|
| Average mammalian transcription | ~2 mRNA/h | Schwanhäusser et al. 2011 (Nature) | 4sU labeling + RNA-seq + ribo-prof + MS | Genome average |
| Peak Pol II initiation | ~1 Pol II/s (~3,600/h) | GAL10, via review† | Live MS2 imaging (yeast) | Burst ceiling |
| Highly transcribed density | ~1 Pol II/270 bp; ~20 nascent/5.4 kb | EMBO Rep 2011† | Live-cell MS2 | Continuous firing |
| Re-initiation interval | Seconds–minutes | Bursting review 2025† | smFISH / live imaging | Burst structure |
| Median mRNA half-life | ~9 h | Schwanhäusser et al. 2011 | Labeling decay | **CCND1 shorter (~30 min–2 h); use measured** |

### E. Complex-removal → derepression (hours) — fast-follower ceiling
| Removed | Timing | Reference | Evidence (method) | Notes |
|---|---|---|---|---|
| PRC1 (RING1B) | gone ~2 h; derepress ~4 h | Dobrinić 2021†; Klose 2024† | RING1B-AID + nascent RNA-seq; Zic2 smFISH | Derepresses **while H3K27me3 persists** |
| PRC1, real-time onset | captured live | Szczurek 2024† | Live MS2 post-degron | Single-transcript sensitivity |
| PRC2 (SUZ12) | subset ~12 h | Petracovici & Bonasio 2021† | SUZ12-AID + RNA-seq | Intermediate timescale |

### F. Mark removal after EZH2 inhibition (days) — replication-coupled
| Process | Timing | Reference | Evidence (method) | Notes |
|---|---|---|---|---|
| H3K27me3 t½ (EPZ-6438) | ~1 day; 90% by 3–4 days | Knutson et al. 2012† | Western/MS time course | Passive, replication-coupled |
| H3K27me3 (GSK126) | <24 h onset, max ~2 days | McCabe et al. 2012 (Nature) | Western time course | — |
| H3K27me3 (EPZ6438, bivalent) | ↓ by 12 h, persists ≥72 h | osteoblast paper† | Western/ChIP time course | Bivalent context |

### G. H3K27me3 restoration after replication (hours–~1 day) — sets τ_restore
| Process | Timing | Reference | Evidence (method) | Notes |
|---|---|---|---|---|
| Positional inheritance | Immediate (recycling) | Reverón-Gómez et al. 2018 | ChOR-seq | Pattern fast; levels slow |
| Level restoration | Across cycle, locus-specific | Reverón-Gómez et al. 2018 | ChOR-seq | Heterogeneity source |
| Restoration vs density | Faster at dense domains; H1 accelerates | Escobar et al. 2023† | ChOR-seq + H1 perturbation | Large dense Ccnd1 → faster |
| Mammalian timescale | ~1 day (vs ~2 h plant) | Howard-lab hybrid 2021 (eLife)† | Modeling + cross-species | τ_restore anchor |
| Two propagation modes | distinct old/new kinetics | Alabert et al. 2015† | SILAC MS | A–U–M calibration |

### H. Dilution-driven derepression, physiological (days/divisions) — THE MODEL TIMESCALE
| Process | Timing | Reference | Evidence (method) | Notes |
|---|---|---|---|---|
| Dilution-driven derepression | Multiple divisions; proliferation-dependent | Jadhav et al. 2020 (Mol Cell) | EZH2-null + ChIP + expression | **Rate-limiting in the model** |
| Replication-rate dependence | 0.3 vs up to 4 div/day | Jadhav et al. 2020 | Lineage-resolved | T_cc sets dilution phenotype |
| Derepression threshold | ~40% residual silences; poised need more erasure | Jadhav et al. 2020 | Quantitative ChIP vs expression | Sets K (high for bivalent Ccnd1) |

### Timescale hierarchy (figure axis)
| Band | Processes | Governs |
|---|---|---|
| Seconds | PRC2 sampling; Pol II re-initiation; nascent-RNA sensing | Writer accessibility; transcription↔PRC2 gate |
| Minutes | Fork transit; Pol II elongation | Instantaneous dilution event |
| Hours | PRC1 removal (~4 h); PRC2 removal (~12 h); transcript accumulation | Fast-follower readout (degron) |
| ~1 day | Mark restoration (τ_restore); EZH2i mark t½ | Recovery vs dilution race |
| Days/divisions | Physiological dilution-driven derepression | **Rate-limiting derepression timescale** |

---

## 5. Findings → model parameters (synthesis)

| Finding | Constrains | Value / form |
|---|---|---|
| Complex is proximal; controls initiation/burst frequency | Readout form | `v = v_init·(η₀ + (1−η₀)/(1+(M/K)ⁿ))`, frequency-acting |
| Mark rate-limiting; complex fast-follower | Readout dynamics | Instantaneous Hill in `M` (quasi-steady-state OK) |
| Nascent RNA inhibits PRC2 (real-time) | Methylation gating | `k_w·Z·(1 − g(v))`; M_ss emerges from antagonism |
| ~40% residual silences; poised = higher threshold | K, n | K ≈ 0.5–0.6·M_ss; n ≈ 4–8 |
| Restoration ~1 day, faster at dense domains | τ_restore | 10–18 h (fast end) |
| EZH2i t½ ~1 day, multi-division derepression | dilution/maintenance | ½ at early S; days to cross K |
| Read-write allosteric; de novo floor needed | k_w, k₀ | catalytic-rate term + small k₀ > 0 |
| ~2/h average, robust 10–50/h, ceiling ~1/s | v_init | ~10–50 mRNA/h for active Ccnd1 |
| Poised Ser5P Pol II; K27me3-coincident pausing | readout delay | optional small delay; check pausing index |
| EZH2 is E2F target | Z(t) | E2F-gated capacity waveform (existing model) |

---

## 6. Open gaps

- **No direct measurement of H3K27me3 restoration kinetics at a large bivalent locus like *Ccnd1*** — τ_restore is currently inferred from genome-wide / other-locus data. A locus-resolved ChOR-seq or nascent-chromatin ChIP at *Ccnd1* would pin it.
- **PRC2-complex (degron) derepression timing at a large bivalent locus** is uncharacterized — the SUZ12 ~12 h figure is a genome-wide "subset" number in ESCs.
- **Strength of the transcription→PRC2 gate in vivo** is debated; constrain it from our own bivalent M_ss and EZH2i/degron kinetics rather than assuming.
- ***Ccnd1* Pol II pausing index** in our cells decides whether a transcriptional-delay term is needed.
- ***Ccnd1* mRNA half-life** (gene-specific, ARE-regulated) should be measured, not taken from the genome median.

---

## 7. Consolidated reference list

*Grouped by theme. † = verify bibliographic details before publication.*

**Mechanism / mark vs complex**
- Blackledge NP & Klose RJ (2021). *The molecular principles of gene regulation by Polycomb repressive complexes.* Nat Rev Mol Cell Biol 22:815–833.
- Dobrinić P, Szczurek AT, Klose RJ (2021). *PRC1 drives Polycomb-mediated gene repression by controlling transcription initiation and burst frequency.* Nat Struct Mol Biol 28:811–824.†
- Klose lab (2024). *The Polycomb system sustains promoters in a deep OFF state by limiting pre-initiation complex formation.*†
- Szczurek AT et al. (2024). Live-cell MS2 imaging after RING1B degron.†
- Finogenova K et al. (2020). *Structural basis for PRC2 decoding of active histone marks H3K36me2/3.* eLife.

**Two-timescale derepression**
- Petracovici A & Bonasio R (2021). *Distinct PRC2 subunits regulate maintenance and establishment of Polycomb repression during differentiation.* Mol Cell.†
- McCabe MT et al. (2012). *EZH2 inhibition as a therapeutic strategy (GSK126).* Nature 492:108–112.
- Knutson SK et al. (2012). *A selective inhibitor of EZH2 (EPZ-6438).* Nat Chem Biol.†

**Transcription ↔ PRC2 coupling / poised Pol II**
- Stock JK et al. (2007). *Ring1-mediated ubiquitination of H2A restrains poised RNAPII at bivalent genes.* Nat Cell Biol 9:1428–1435.
- Brookes E et al. (2012). *Polycomb associates genome-wide with a specific RNAPII variant in ESCs.* Cell Stem Cell.†
- Ferrai C et al. (2017). *RNA polymerase II primes Polycomb-repressed developmental genes throughout terminal neuronal differentiation.* Mol Syst Biol.
- Davidovich C et al. (2013). *Promiscuous RNA binding by PRC2.* Nat Struct Mol Biol 20:1250–1257.
- Kaneko S et al. (2013; 2014). PRC2 binds nascent RNA at active promoters / keeps PRC2 in check. Genes Dev.†
- Cifuentes-Rojas C et al. (2014). RNA regulation of PRC2.
- Beltran M et al. (2016). *The interaction of PRC2 with RNA or chromatin is mutually antagonistic.* Genome Res 26:896.†
- Wang X et al. (2017). *Targeting of PRC2 to RNA by short repeats of consecutive guanines.* Mol Cell 65:1056–1067.
- Rosenberg M et al. (2021). *Motif-driven interactions between RNA and PRC2 are rheostats that regulate transcription elongation.* Nat Struct Mol Biol 28:103–117.†
- Song J et al. (2023). *Structural basis for inactivation of PRC2 by G-quadruplex RNA.* Science.†

**PRC2 dynamics / recruitment / read-write**
- Margueron R et al. (2009). *Role of the Polycomb protein EED in propagation of repressive histone marks.* Nature 461:762–767.
- Youmans DT et al. (2018). *Live-cell imaging reveals the dynamics of PRC2 and recruitment to chromatin by SUZ12-associated subunits.* Genes Dev 32:794.
- DIPG H3.3K27M single-molecule dynamics (2018). Nat Commun.†

**Dilution / restoration**
- Reverón-Gómez N et al. (2018). *Accurate recovery of the chromatin landscape after DNA replication.* Mol Cell 72:239–249.
- Escobar TM et al. (2023). H1 / ChOR-seq restoration.†
- Alabert C et al. (2015). *Two distinct modes for propagation of histone PTMs across the cell cycle.* Genes Dev.†
- Jadhav U et al. (2020). *Replicational dilution of H3K27me3 in mammalian cells and the role of poised promoters.* Mol Cell 78:141–151.

**Transcription rates**
- Schwanhäusser B et al. (2011). *Global quantification of mammalian gene expression control.* Nature 473:337–342.
- *Fast transcription rates of RNA polymerase II in human cells* (2011). EMBO Rep.†

**E2F–EZH2 coupling**
- Bracken AP et al. (2003). *EZH2 is downstream of the pRB-E2F pathway, essential for proliferation and amplified in cancer.* EMBO J 22:5323–5335.
- Chen J et al. (2018). EZH2 splice variant in the cell cycle. Epigenetics & Chromatin.

**Pausing / bivalency**
- Wang H et al. (2023). *H3K4me3 regulates RNAPII promoter-proximal pause-release.* Nature 615:339–348.
- Liu J et al. (2017). *Dynamics of RNAPII pausing and bivalent H3 methylation during neuronal differentiation.* Cell Reports 20:1307–1318.

**Modeling lineage**
- Dodd IB et al. (2007). *Theoretical analysis of epigenetic cell memory by nucleosome modification.* Cell 129:813–822.
- Angel A et al. (2011). *A Polycomb-based switch underlying quantitative epigenetic memory.* Nature 476:105–108.
- Berry S, Dean C, Howard M (2017). *Slow chromatin dynamics allow Polycomb target genes to filter fluctuations in TF activity.* Cell Systems 4:445–457.
- Howard-lab hybrid protein-assembly/histone-modification model (2021). eLife.†
- Klingel V et al. (2021). *Model-based robustness and bistability analysis for methylation-based epigenetic memory.* FEBS J.†

---

## 8. How this maps to the implemented model

The mean-field AUM `with_h3k27_dilution` module (see `h3k27me3_repression_model_state.md` §2) implements
this evidence base as a single-state reduction:

| Lit-review feature (§) | Implementation |
|---|---|
| Read-write + de-novo floor, allosteric (§3.6) | `EZH2·(k_w·Mk + k0)·(1−Mk)` |
| EZH2 = E2F-coupled writer capacity Z(t) (§3.7) | whole methylation term ∝ EZH2 (falls when the cell exits the cycle) |
| Transcription→PRC2 reciprocal eviction (§3.3) | `·(1 − g·Cd_mRNA^p/(K_tx^p+Cd_mRNA^p))` — M_ss emerges, arrest self-reinforces |
| τ_restore 10–18 h; EZH2i t½ ~1 day (§3.4, F/G) | `δ=0.0007` → de-repression t½ ≈ 16 h |
| Replicative dilution, ½ at S; T_cc phenotype (§3.4, H) | halving event at `Dna>0.05`; dilution lowers Mk 17–25 % |
| Leaky, frequency-acting readout w/ floor (§3.1, §3.8) | `R = f0 + (1−f0)/(1+(Mk/K)^n)` |
