# PRC2 in vitro binding literature — census of affinity, kinetics and substrate preference

Compiled 2026-10-23.

## Reading standard applied

Every number below was read in the **Results, Methods, Table or figure-legend text of the primary
paper**. Numbers available only from an abstract or a review are marked `UNVERIFIED — abstract only`.
Numbers that exist only as a curve inside a figure *image* (not in any machine-readable text) are
marked `value in figure image only`, and the fold-change or qualitative statement from the text is
given instead of a fabricated value.

Confidence grades:

| Grade | Meaning |
|---|---|
| **A** | Direct measurement, read in the primary paper's own text/table |
| **B** | Model-derived (simulation, fit-extrapolation, ODE model) |
| **C** | Secondary source, review, or abstract only |
| **X** | Not found / not measured anywhere I could reach |

Reading was split across three readers (this agent plus two delegated readers) working to the same
standard. Rows marked † were read and re-verified directly by the compiling agent against the
primary full text; the rest were read once by a delegated reader from the primary full text.

**The single most important caveat in this whole file:** these K<sub>d</sub> values span 10–200 mM
monovalent salt, 0–35 bp vs 60–100 bp linker DNA, 4–30 °C, and presence/absence of yeast tRNA
competitor. Values measured under different conditions are **not interchangeable**. Reading the salt
column is not optional.

---

## 1. PRC2 binding to free DNA

| # | Quantity | Value | Technique | Substrate | Enzyme | Salt | Temp | n | Paper | Where | Conf |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1.1 † | K<sub>d</sub>, PRC2 core | 64 ± 4 nM | EMSA (1.2% agarose), 45 nM Flc probe | **30 bp dsDNA** "PRE 11L" (GC-rich) | human PRC2 4-mer EZH2/EED/SUZ12/RBBP4 | **50 mM NaCl**, 25 mM HEPES pH 7.9, 5 mM MgCl₂, 0.05% Tween, 4% glycerol | 5 min **on ice**; gel at 4 °C | 3 | Choi 2017 NSMB 24:1039 (PMID 29058710) | Results "PRC2 binds dsDNA…"; Fig 2b | A |
| 1.2 † | K<sub>d</sub>, PHF1–PRC2 | **30 ± 3 nM** | EMSA | 30 bp PRE 11L | human PHF1–PRC2 5-mer | 50 mM NaCl | on ice | 3 | Choi 2017 | Fig 2b | A |
| 1.3 † | K<sub>d</sub>, PHF1<sup>C</sup>–PRC2 | 49 ± 4 nM | EMSA | 30 bp PRE 11L | PRC2 core + PHF1(515–567) only | 50 mM NaCl | on ice | 3 | Choi 2017 | Fig 2b | A |
| | | *Authors' own caveat: probe was 45 nM, i.e. comparable to K<sub>d</sub>, so these may **underestimate** the true affinities (stated in Fig 2b legend).* | | | | | | | | | |
| 1.4 † | K<sub>d,app</sub>, PRC2 5-mer on long DNA, **AFM** | **150 ± 12 nM**, Hill n = 1.1 ± 0.1 (all-points fit 140 ± 8 nM) | AFM in liquid, single-molecule counting | **2.5 kb** DNA (12×601), 20 nM | human 5-mer +AEBP2 | **25 mM KCl**, 38 mM HEPES pH 7.5, 0.75 mM ZnCl₂, 1 mM TCEP; deposition 25 mM KCl + 10 mM MgCl₂ | bind 30 min at **30 °C**; image at **19 °C** | >70 molecules per conc.; 1152 molecules total | Heenan 2020 NAR 48:2969 (PMC7102998) | Abstract + Results "Quantifying PRC2 binding by AFM"; Fig 3A | A |
| 1.5 † | K<sub>d,app</sub>, same complex, **EMSA control** | **75 ± 5 nM**, Hill n = 1.7 ± 0.2 (global fit to three 2.5-kb sequences) | EMSA, ³²P | 2.5 kb 12×601, GC-poor, GC-poor+76% GC island | human 5-mer +AEBP2 | **100 mM KCl**, 50 mM Tris pH 7.5, 0.1 mM ZnCl₂, 0.05% NP-40, 0.5 mM EDTA, 5% glycerol, **no tRNA** | 1 h at 30 °C | 3 sequences | Heenan 2020 | Methods EMSA; Fig 3F | A |
| 1.6 † | K<sub>d,app</sub> for **DNA looping** (as distinct from binding) | **900 ± 400 nM**, Hill n = 1.4 ± 0.3 | AFM, fraction of looped molecules | 2.5 kb 12×601 | human 5-mer +AEBP2 | 25 mM KCl | 19 °C | as 1.4 | Heenan 2020 | Fig 3A (brown curve) | A |
| 1.7 | K<sub>d,app</sub>, 60 bp dsDNA **vs salt** | **5.1 ± 0.60 nM** (10 mM KCl); **82 ± 13 nM** (25 mM); **170 ± 15 nM** (100 mM); **not measurable** (200 mM). log K<sub>d</sub> vs log[KCl] slope **m ≈ 1.4 ± 0.68** | fluorescence polarization | ds-d(N)₆₀ | human PRC2 5-mer | 10 / 25 / 100 / 200 mM KCl, no tRNA | 25 °C | ≥3 | Hemphill 2023 PNAS (PMC10266057) | Table 1; Fig 3B | A |
| 1.8 | K<sub>d,app</sub>, 50 bp dsDNA | 5.0 ± 0.46 nM (10 mM KCl); **390 ± 31 nM** (25 mM KCl) | FP | ds-d(N)₅₀ | human PRC2 5-mer | as above | 25 °C | ≥3 | Hemphill 2023 | Table 1 | A |
| 1.9 † | K<sub>d,app</sub>, (mCpG)₂₄ vs (CpG)₂₄ | **11.3 ± 6.2 nM** methylated; **~50-fold weaker** unmethylated; requires AEBP2 (PRC2 ΔAEBP2 4-mer does not discriminate) | EMSA | 48 bp dsDNA | human PRC2 5-mer +AEBP2 | **100 mM KCl** | 30 °C, 30 min, no tRNA | 3 | Wang 2017 NSMB 24:1028 (PMC5771497) | Results; Fig 5b,c | A |
| 1.10 † | (CG)₃₀ vs (TA)₃₀ | CG-rich bound with "substantially lower K<sub>d,app</sub>" | fluorescence polarization | 60 bp dsDNA | PRC2 5-mer | **25 mM KCl** | 30 °C | 3 | Wang 2017 NSMB | Fig 5a | A (qualitative); **value in figure image only** |
| 1.11 | K<sub>d</sub>, PRC2–MTF2 holo on CpG-island DNA | **6.3 nM** (no competitor); **15.0 nM** (+200 ng/µl yeast tRNA) | EMSA, 200 pM probe | 100 bp CGI_Lhx6 genomic DNA | human PRC2.1 5-mer +MTF2 | 100 mM KCl | on ice, 1 h | ≥3 | Chen 2020 Mol Cell 77:1265 (PMC7571800) | Results; STAR Methods | A |
| 1.12 | K<sub>d</sub>, PRC2–AEBP2<sup>short</sup> vs core vs AEBP2<sup>long</sup> | **45.5 ± 3.6 nM** (AEBP2-S, Hill 2.3 ± 0.4); core PRC2, PRC2–AEBP2-L(iso1), (iso2) all **> 4000 nM** | EMSA, 5 nM probe | 24 bp CpG hairpin (CDKN2B CGI) | human core PRC2 ± AEBP2 isoform | **100 mM KCl** | 30 °C, 30 min | 3 | Mucha 2025 EMBO J (PMC12669776) | Fig 1C,D | A |
| 1.13 | K<sub>d</sub>, core PRC2 vs PRC2–PALI1 | core **> 4 µM**; PRC2–PALI1<sup>PIR</sup> **155 ± 26 nM**; GC-poor control DNA **73.7 ± 10.0 nM** → not CpG-specific | fluorescence anisotropy, 5 nM probe | 46 bp dsDNA hairpin | human core PRC2 ± PALI1<sup>PIR</sup> | **100 mM KCl** | 30 °C | 3 | Zhang 2021 Nat Commun 12:4592 (PMC8319299) | Fig 6a | A |
| 1.14 | PRC2.1<sup>MTF2</sup> ± EPOP | 3.4 nM (−EPOP) → **9.2 nM** (+EPOP), i.e. EPOP **weakens** DNA binding 2–3× | EMSA, 0.2 nM probe, +0.2 mg/ml tRNA | 100 bp CGI_Lhx6 DNA | human PRC2.1 ± EPOP | 100 mM KCl | on ice, 30 min | 3 | Gong 2026 Nat Commun (PMC12901234) | Results; Methods | A |
| 1.15 | Composition dependence (AFM) | Full 5-mer binds robustly; **catalytic moiety** (EZH2 + SUZ12<sup>VEFS</sup> + EEDΔ81) and **regulatory moiety** (SUZ12ΔVEFS + AEBP2 + RBBP4) both at background | AFM, 30 and 90 nM | 2.5 kb 12×601 | three PRC2 variants | 25 mM KCl | 19 °C | ≥70–200 molecules per point | Heenan 2020 | Fig 5 | A (qualitative) |
| 1.16 † | **Sequence specificity on long DNA — negative result** | No spatial preference for a 200 bp 76%-GC CpG island embedded in 2.5 kb low-GC DNA, by either position histogram or EMSA. Positive control BspMI localised to 320 ± 2 nm vs expected 323 nm (N = 91) | AFM localisation + EMSA | 2.5 kb constructs | human 5-mer +AEBP2 | 25 mM KCl (AFM) / 100 mM KCl (EMSA) | 19 / 30 °C | N = 2081 molecules | Heenan 2020 | Figs 3F, 4B,C | A |
| 1.17 † | Isolated Polycomblike winged-helix domain | Pcl<sup>PHD2–WH</sup> **9.7 ± 0.2 µM**; bxd-PRE 20-mer panel **3–12 µM** with no sequence preference, CpG-containing and CpG-free probes equal; human PHF1<sup>PHD2–WH</sup> "midmicromolar"; wing-1 and helix-α3 mutants lose binding | fluorescence polarization, 45 nM Flc probe | 30 bp PRE 11L; 20 bp bxd PRE probes | *Drosophila* Pcl(491–694), human PHF1(165–363) — **isolated domains, not in PRC2** | 50 mM NaCl | **20 °C**, 5 min | 3 | Choi 2017 | Fig 4c,d; Suppl Fig 4 | A |
| 1.18 | Isolated PHF1/MTF2 reader cassettes, CpG DNA (ITC) | PHF1(26–360) **1.2 ± 0.3 µM**; PHF1(165–360) **0.5 ± 0.1 µM**; MTF2(180–378) **2.1 ± 0.3 µM**. Sequence series NCpGN: PHF1 0.8 → 31 µM, MTF2 6.4 → 33 µM; ACpGA/ACpGT/TCpGA no binding | ITC | 12–14 bp dsDNA | isolated domains, not in PRC2 | **150 mM NaCl** | **10 °C** | ≥2 | Li 2017 Nature 549:287 (PMC5937281) | Ext. Data Table 2 | A |

> **Conflict worth noting.** Choi 2017 found the isolated Pcl/PHF1 WH domain to be *sequence-non-specific*
> (rows 1.17) and explicitly contradicts Li 2017's CpG-specific binding model in its Discussion. Li 2017
> (row 1.18) reports a ~40-fold spread across NCpGN sequences by ITC at 150 mM NaCl. Different technique,
> different salt, different construct boundaries. Both are grade A; they disagree.

---

## 2. PRC2 binding to nucleosomes — the linker-DNA dependence

This is the cleanest quantitative story in the field and it is **all from one paper's internal series**,
so the comparisons are self-consistent.

| # | Quantity | Value | Technique | Substrate | Enzyme | Salt | Temp | n | Paper | Where | Conf |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2.1 † | K<sub>d,app</sub>, **nucleosome core particle, NO linker** | **40 ± 16 µM** | EMSA, H2B-Cy5 and ³²P | NCP, 147 bp | human PRC2 5-mer +AEBP2 | **100 mM KCl**, 50 mM Tris pH 7.5, 2.5 mM MgCl₂, 0.1 mM ZnCl₂, 0.05% NP-40, 0.1 mg/ml BSA, 5% glycerol, **no tRNA** | **30 °C**, 30 min | 3 | Wang 2017 NSMB | Results; Fig 3a,b | A |
| 2.2 † | K<sub>d,app</sub>, **mononucleosome with linker** | **280 ± 19 nM** | EMSA, ³²P | 207 bp template (147 bp 601 + 60 bp linker) | same | 100 mM KCl | 30 °C | 3 | Wang 2017 NSMB | Fig 1b,c | A |
| 2.3 † | K<sub>d,app</sub>, free 147 bp DNA | strongest of the series (tighter than any nucleosome) | EMSA | protein-free 147 bp | same | 100 mM KCl | 30 °C | 3 | Wang 2017 NSMB | Fig 3b | A (ordering); **value in figure image only** |
| 2.4 † | K<sub>d,app</sub>, **trinucleosome** | **23 ± 9 nM** | EMSA, ³²P | 3×601, **60 bp linkers** | same | 100 mM KCl | 30 °C | 3 | Wang 2017 NSMB | Fig 1b,c | A |
| 2.5 † | K<sub>d,app</sub>, **dodecanucleosome** | **9 ± 2 nM** | EMSA, ³²P | 12×601, 60 bp linkers | same | 100 mM KCl | 30 °C | 3 | Wang 2017 NSMB | Fig 1b,c | A |
| | | *Not normalised per nucleosome. Affinity per nucleosome rises mono→tri then plateaus. Hill coefficient rises substantially from mono to arrays.* | | | | | | | | | |
| 2.6 † | K<sub>d,app</sub>, **dinucleosome, 10 bp linker** | **276 ± 16 nM** (substantially weaker than 50 or 100 bp linker) | EMSA, fluorescent dinucleosome | CpG-island promoter dinucleosome mimetic, 10 / 50 / 100 bp linker | PRC2 5-mer | 100 mM KCl | 30 °C | 3 | Wang 2017 NSMB | Results "PRC2 prefers to bind and methylate long-linker dinucleosomes"; Fig 4a,b | A (10 bp value); 50/100 bp **values in figure image only** |
| 2.7 | K<sub>d</sub>, core PRC2, **NCP147 vs NCP182** | NCP147 (no linker) **2560 ± 1250 nM**; NCP182 (35 bp linker) **330 ± 30 nM** → ~8× | EMSA, 5 nM Cy5 probe + salmon-sperm DNA | NCP147 vs NCP182 | human core PRC2 4-mer | **100 mM KCl** | **4 °C**, 30 min | 3 | Zhang 2021 Nat Commun (PMC8319299) | Figs 6b, 7a,b | A |
| 2.8 † | K<sub>d</sub>, PRC2 core on mononucleosome **with 34 bp linkers both sides** | **163 ± 13 nM** | EMSA, 45 nM Flc-mononucleosome | 215 bp (147 bp 601 + 34 bp each side) | human PRC2 4-mer (no AEBP2) | **50 mM NaCl** | on ice, 5 min | 3 | Choi 2017 | Fig 2a | A |
| 2.9 † | Same, PHF1–PRC2 / PHF1<sup>C</sup>–PRC2 | **66 ± 4 nM** / **211 ± 9 nM** | EMSA | 215 bp mononucleosome | human 5-mer | 50 mM NaCl | on ice | 3 | Choi 2017 | Fig 2a | A |
| 2.10 | K<sub>d</sub>, PRC2–PALI1<sup>PIR</sup> on mononucleosome | NCP182 **19.0 ± 0.6 nM** (>15× tighter than core); NCP147 **48.5 ± 1.2 nM** → PALI1 gives >50× boost **even with no linker** | EMSA | NCP147, NCP182 | human PRC2 + PALI1<sup>PIR</sup> | 100 mM KCl | 4 °C | 3 | Zhang 2021 Nat Commun | Figs 6b, 7a,b | A |
| 2.11 | K<sub>d</sub>, PRC2.2 (AEBP2-L + JARID2) on mononucleosome | **24 nM** | EMSA, 5 nM Cy5 | NCP182 (35 bp linker) | PRC2–AEBP2<sup>L</sup>–JARID2 | **100 mM NaCl**, 2.5 mM MgCl₂, 0.1 mM ZnCl₂ | gel 4 °C | 3 | Mucha 2025 EMBO J | Fig EV8A,B | A |
| 2.12 | PRC2.1<sup>MTF2</sup> ± EPOP on mononucleosome | 122.9 nM (−EPOP) → **295.4 nM** (+EPOP) | EMSA, 0.5 nM probe | 601 + 30 bp extension both ends (~207 bp) | human PRC2.1 ± EPOP ± Elongin BC | **100 mM NaCl** | on ice, 30 min | 3 | Gong 2026 Nat Commun | Results; Methods | A |
| 2.13 | K<sub>d</sub>, PRC2 catalytic core (PRC2c) on mononucleosome | **25.9 ± 10.7 nM** (unmodified, 147 bp) and **29.9 ± 3.9 nM** (H3K27M, 147 bp) | fluorescence-intensity titration (MDCC-H3), 10 nM nucleosome | 147 bp core particle | human PRC2c (EZH2–EED–SUZ12<sup>VEFS</sup>) | **40 mM KCl**, 400 µM SAM, no tRNA | **20 °C** | 3 | Beltran 2019 NSMB (PMC6778522) | Suppl Fig 4c; Fig 3c,d | A |
| 2.14 | K<sub>d</sub>, PRC2c on mononucleosome (EMSA) | **0.4 ± 0.1 µM** WT | EMSA, 20 nM nucleosome | recombinant mononucleosome | human PRC2c | 100 mM KCl | not stated | 3 | Houslay 2026 **bioRxiv preprint** (doi 10.64898/2026.10.15.755990) | Fig 4b | A (preprint, not peer-reviewed) |
| 2.15 | EPOP vs MTF2 on **dinucleosomes** and linker sequence | Both strongly increase PRC2 dinucleosome affinity over core. PRC2–EPOP > PRC2–MTF2 on **(GA)ₙ** linker; PRC2–MTF2 > PRC2–EPOP on **(GCN)ₙ** linker. On matched naked DNA both help but **without** sequence specificity → preference is chromatin-dependent | EMSA, Hill fit | dinucleosome, two 601 cores + **40 bp linker** of defined sequence | human PRC2-core, PRC2–EPOP, PRC2–MTF2 | **100 mM NaCl**, 2.5 mM MgCl₂, 0.1 mM ZnCl₂ | **30 °C**, 30 min | 3 biological replicates | Kwon/Reinberg 2026 PNAS (PMC12890814) | Fig 4 | A (comparisons); **absolute K<sub>d</sub> in figure image only** |

> **Note on the two sub-µM "no-linker" numbers.** Wang 2017 gets 40 µM for a 147 bp NCP with the
> 5-mer at 100 mM KCl; Zhang 2021 gets 2.56 µM for NCP147 with the 4-mer at 100 mM KCl and 4 °C;
> Beltran 2019 gets 26 nM for a 147 bp NCP with the minimal 3-subunit PRC2c at 40 mM KCl with SAM.
> That is a ~1500-fold spread for nominally the same substrate class. Salt (40 vs 100 mM), temperature
> (20 vs 4 vs 30 °C), enzyme composition and readout (gel vs solution fluorescence) all differ. Do not
> average these. If you need one number for a model, pick the one whose conditions match your use case
> and carry the condition with it.

---

## 3. Kinetics: k<sub>on</sub>, k<sub>off</sub>, residence times

### 3.1 Single-molecule TIRF — Choi et al. 2017, full Table 1 †

Conditions for all rows: human PRC2 complexes, EZH2 N-terminally labelled with DY-547 (>50% labelled),
injected at **1–2 nM**. Imaging buffer **50 mM HEPES pH 7.5, 50 mM NaCl, 10% glycerol, 0.005% Tween-20,
2 mM Trolox, 3.2% glucose, 5 mg/ml BSA, glucose oxidase/catalase**. 50 ms frames, 5000 frames,
532 nm at 40 W/cm². Dwell-time histograms corrected for DY-547 photobleaching (single-exponential,
τ = 18 s) then fit to y(t) = A₁e^(−t/τ_off,1) + A₂e^(−t/τ_off,2). Temperature not stated (room temperature
implied). Substrates: 12×601 arrays with **30 bp linker**; 177 bp naked 601 DNA; 601 mononucleosomes.

| Substrate | Protein | τ<sub>off,1</sub> (s) | % A₁ | **τ<sub>off,2</sub> (s)** | % A₂ | λ<sub>on</sub> (s) | k<sub>on</sub> (×10⁶ M⁻¹s⁻¹) | n |
|---|---|---|---|---|---|---|---|---|
| 12-mer array | PRC2 | 0.29 ± 0.03 | 62 ± 6 | **2.28 ± 0.44** | 38 ± 6 | 15.53 ± 6.48 | 32.05 ± 13.44 | 3 |
| 12-mer array | PHF1–PRC2 | 0.12 ± 0.08 | 46 ± 4 | **6.48 ± 0.52** | 54 ± 4 | 34.21 ± 29.66 | 14.62 ± 12.67 | 3 |
| 12-mer array | PHF1<sup>C</sup>–PRC2 | 0.22 ± 0.05 | 62 ± 4 | 2.86 ± 0.07 | 38 ± 4 | 13.82 ± 3.55 | 36.19 ± 9.29 | 3 |
| 12-mer array | PHF1<sup>WH>E</sup>–PRC2 | 0.17 ± 0.11 | 35 ± 11 | 3.11 ± 0.91 | 65 ± 11 | 22.59 ± 8.22 | 44.27 ± 16.10 | 3 |
| Mononucleosome | PRC2 | 0.30 ± 0.10 | 75 ± 10 | 3.53 ± 0.90 | 25 ± 10 | 17.68 ± 7.42 | 28.28 ± 11.87 | 4 |
| Mononucleosome | PHF1–PRC2 | 0.16 ± 0.11 | 46 ± 15 | **9.25 ± 4.87** | 54 ± 15 | 37.03 ± 10.77 | 13.5 ± 3.93 | 5 |
| Mononucleosome | PHF1<sup>C</sup>–PRC2 | 0.22 ± 0.31 | 72 ± 7 | 4.17 ± 1.46 | 28 ± 7 | 25.28 ± 21.08 | 19.78 ± 16.49 | 3 |
| Mononucleosome | PHF1<sup>WH>E</sup>–PRC2 | 0.15 ± 0.11 | 50 ± 19 | 5.34 ± 0.13 | 50 ± 19 | 20.48 ± 5.30 | 24.41 ± 6.31 | 4 |
| **601-DNA (177 bp)** | PRC2 | 0.38 ± 0.06 | 68 ± 5 | **4.45 ± 1.90** | 32 ± 5 | 11.72 ± 3.63 | 42.67 ± 13.24 | 3 |
| **601-DNA** | PHF1–PRC2 | 0.25 ± 0.04 | 44 ± 8 | **10.38 ± 2.46** | 56 ± 8 | 18.39 ± 4.55 | 27.19 ± 6.72 | 5 |
| **601-DNA** | PHF1<sup>C</sup>–PRC2 | 0.27 ± 0.08 | 70 ± 5 | 6.83 ± 1.63 | 30 ± 5 | 12.21 ± 0.51 | 40.95 ± 1.71 | 4 |
| **601-DNA** | PHF1<sup>WH>E</sup>–PRC2 | 0.23 ± 0.10 | 49 ± 11 | 5.71 ± 1.15 | 51 ± 11 | 19.02 ± 4.01 | 26.30 ± 5.55 | 4 |

Grade **A** throughout. Statistics: two-sided heteroscedastic t-test. Chromatin τ_off,2 PRC2 vs PHF1–PRC2
P = 0.00049; DNA P = 0.011; **mononucleosome P = 0.057, i.e. the PHF1 effect on mononucleosomes is not
significant at α = 0.05** — worth knowing, since the 3.53 → 9.25 s number looks like the biggest effect
in the table but has n = 4/5 and a wide SD.

Readings that follow from this table:

- τ_off,2 is ~2× longer on **naked DNA** than on 12-mer chromatin for every complex → DNA, not histone,
  carries the affinity.
- k<sub>on</sub> is **not significantly different** between complexes (all 1.4–4.4 × 10⁷ M⁻¹s⁻¹). PHF1 works
  entirely through k<sub>off</sub>, not k<sub>on</sub>.
- 46–75% of all binding events are sub-second (τ_off,1 ≈ 0.1–0.4 s) — "the" residence time is bimodal.
- The authors' own equilibrium cross-check: K<sub>d</sub> estimated from k<sub>on</sub>·τ_off would be
  low-nM, but they explicitly decline to do this ("uncertainties in the on-rate constants preclude an
  exact determination of equilibrium affinities") and instead measure K<sub>d</sub> by EMSA/MST.

### 3.2 Bulk dissociation kinetics

| # | Quantity | Value | Technique | Conditions | Paper | Conf |
|---|---|---|---|---|---|---|
| 3.2.1 † | Half-life of pre-formed PRC2–dsDNA complex on RNA challenge | **t½ = 91 ± 5 s** (main text); **τ = 94.6 ± 4.4 s** (Fig 5e legend) | FP competitive dissociation, (GGAA)₁₀ competitor | 5 nM A488-dsDNA, **25 mM KCl**, 30 °C, no tRNA, n = 3 | Wang 2017 NSMB, Fig 5e | A |
| 3.2.2 | Intrinsic k₋₁ from 60 bp dsDNA | **2.4 ± 0.25 × 10⁻³ s⁻¹** (t½ ≈ 4.8 min) | FP competitive dissociation, full competitor series | **10 mM KCl**, 25 °C, no tRNA, ≥3 | Hemphill 2023 PNAS, Table 1 | A |
| 3.2.3 | Intrinsic k₋₁ from G4 RNA | **1.4 ± 1.1 × 10⁻³ s⁻¹** (t½ ≈ 8 min) | as above | 10 mM KCl, 25 °C | Hemphill 2023, Table 1 | A |
| 3.2.4 | **Direct-transfer** rate constant k<sub>θ</sub>, dsDNA → G4 RNA | **170 ± 41 M⁻¹s⁻¹** (60 bp ligand); **170 ± 8.2** (50 bp ligand) | FP competitive dissociation | 10 mM KCl, 25 °C | Hemphill 2023, Table 1 | A |
| 3.2.5 | k<sub>θ</sub>, G4 RNA → dsDNA | **91 ± 31** (ds-d(N)₅₀); **100 ± 19 M⁻¹s⁻¹** (ds-d(N)₆₀) | as above | 10 mM KCl, 25 °C | Hemphill 2023, Table 1 | A |
| 3.2.6 | DNA→DNA self-exchange | k₋₁ **9.1 × 10⁻⁵ s⁻¹**, k<sub>θ</sub> **260 M⁻¹s⁻¹** (n = 1, flagged in paper); 5.3 ± 2.3 × 10⁻⁴ s⁻¹ / 150 ± 21 M⁻¹s⁻¹ (alternative protocol) | FP | 10 mM KCl, 25 °C | Hemphill 2023, Table 1 | A |
| 3.2.7 | Crossover for direct transfer to dominate | k<sub>θ</sub>/k₋₁ ≈ **10⁵ M⁻¹** → direct-transfer flux exceeds simple dissociation above ~10 µM competitor | derived from Table 1 | — | Hemphill 2023, Discussion | A (ratio) / B (biological extrapolation) |

> **Resolved conflict — trapping vs direct transfer.** Wang 2017 NSMB concluded that RNA **traps** PRC2
> after spontaneous dissociation rather than actively invading, because a 10-fold drop in RNA competitor
> did not change the 91 s half-life (Suppl Fig 8c). Hemphill 2023 PNAS repeated this over a full competitor
> concentration series and found k<sub>off</sub><sup>obs</sup> rises **linearly** with competitor with no
> plateau, with the direct-transfer model beating classic competition by BIC in every pairing. Hemphill
> explicitly reconciles: the earlier papers sampled only one or two competitor concentrations (200-fold and
> 2000-fold excess) and are consistent with the direct-transfer model at those points. **Use Hemphill's
> rate constants; the 91 s figure is a single-concentration observable, not an intrinsic k<sub>off</sub>.**

### 3.3 The residence-time trap, stated explicitly

Three different quantities are called "residence time" in this literature and they are **not comparable**:

1. **τ_off from smTIRF** (Choi 2017, §3.1): the lifetime of a fluorescent spot at zero applied force,
   biexponential, 0.1–10 s. This is a true dwell time under equilibrium.
2. **t½ / k₋₁ from competitive dissociation** (§3.2): the lifetime of the complex under a competitor
   challenge. 91 s (Wang) or 5–12 min (Hemphill k₋₁). These are 1–2 orders of magnitude longer than the
   smTIRF τ_off,2 — different substrates, different salt, and a competitor-dependent observable.
3. **Rupture force from optical tweezers** (§4): the force at which a mechanically loaded complex fails
   at a 0.1 µm/s pull rate. This is **not a lifetime at all** — it is a loading-rate-dependent quantity,
   and converting it to an off-rate requires a Bell–Evans analysis that none of these papers performs.

Heenan 2020 itself cross-references "the ~100 s residency time of PRC2 on dsDNA" from Wang 2017 to argue
its 30 min incubation is equilibrated — a legitimate use of quantity (2), but note it is not the smTIRF number.

---

## 4. Bridging, looping and mechanical stability

### 4.1 Optical tweezers — Leicher et al. 2020 PNAS †

Conditions: LUMICKS C-Trap dual-trap, **room temperature**, imaging buffer **10 mM Tris pH 8.0, 0.1 mM EDTA,
200 mM KCl, 0.5 mM MgCl₂, 0.1% Tween, 1 mM DTT**. 12×601 arrays with **30 bp linker** (or 50 bp where
stated), 4 kbp DNA handles. PRC2 at **500 nM** unless noted, free in solution throughout; 5–10 s incubation
at zero force; pull at **0.1 µm/s**. ≥15 single tethers per condition; first pull only. Force resolution
0.2 pN, stability <0.3 pN.

| # | Quantity | Value | Conf |
|---|---|---|---|
| 4.1.1 † | Nucleosome unwrapping force, **no PRC2** | **16.0 ± 0.3 pN**; 12 ± 1 transitions per array; ΔL₀ = 76 ± 1 bp (SD 11 bp) | A |
| 4.1.2 † | C1 cluster — PRC2 sequestering free/linker DNA | ΔL₀ < 51 bp; withstands ~17 pN; seen on bare DNA alone | A |
| 4.1.3 † | C2 — mononucleosome unwrapping, ± PRC2 | 16.2 ± 0.3 pN (no PRC2) vs **16.6 ± 0.3 pN** (+PRC2); ΔL₀ 77 ± 1 bp; **not significant** | A |
| 4.1.4 † | C3 — adjacent-nucleosome bridging (Nuc 1–2) | 105–179 bp; **16.3 ± 0.4 pN** | A |
| 4.1.5 † | **C4 — non-adjacent bridging** | ΔL₀ > 180 bp (86/958 transitions with PRC2; **0/234 without**); **21.4 ± 1.0 pN** — significantly higher than C2 and C3 | A |
| 4.1.6 † | **Preferred geometry** | 78% of C4 transitions are 180–300 bp → **Nuc 1–3** (one spacer nucleosome) is the dominant bridging mode | A |
| 4.1.7 † | Linker-length test | 50 bp linker array: C4 mean ΔL₀ is **~40 bp larger** than 30 bp linker array (= 2 × linker increment, confirming Nuc 1–3). C3 shifts only ~7 bp | A |
| 4.1.8 † | Concentration dependence of bridging | C4 population rises from **5 nM to 100 nM** PRC2 and **saturates above ~100 nM**. C1 (free-DNA sequestration) keeps rising past that — consistent with the ~900 nM looping K<sub>d</sub> from AFM (row 1.6) | A |
| 4.1.9 † | C4 fraction vs linker | 9% of transitions (30 bp linker) vs 15% (50 bp linker) | A |
| 4.1.10 † | **Dinucleosome control** | C4 essentially absent on a 2-mer substrate (P = 5.9 × 10⁻²⁷, KS test) → C4 needs ≥3 nucleosomes, rules out free-DNA-only or one-nucleosome-plus-DNA bridging | A |
| 4.1.11 † | **dCas9 roadblock control** | Two dCas9 roadblocks on bare DNA + PRC2 produced **no** large transitions → the "PRC2 slides until it hits steric blocks" model is **disfavoured** | A |
| 4.1.12 † | SAM effect | SAM raises bare-nucleosome force 16.0 → **17.7 ± 0.2 pN** independent of PRC2; C4 force **unaffected** by SAM | A |
| 4.1.13 † | **H3K27me3 effect on mechanical stability** | **50% H3K27me3** arrays: C3 **14.6 ± 0.2** vs 16.3 ± 0.4 pN unmodified; C4 **16.9 ± 0.8** vs 21.4 ± 1.0 pN → destabilised. **100% H3K27me3** returns to unmodified levels. Non-monotonic | A |
| 4.1.14 † | **H3K27M effect** (with SAM) | Dose-dependent **stabilisation**: overall 17.8 ± 0.3 (WT) → 18.3 ± 0.3 (50%) → **19.6 ± 0.5 pN** (100%). C3 17.8 ± 0.5 → 19.6 ± 0.7; **C4 19.8 ± 1.2 → 26.1 ± 2.6 pN**. H3K27R control: no effect | A |
| 4.1.15 † | AEBP2 | Abolishes C1 (free-DNA sequestration; P = 0.00071 KS); **lowers mean transition force 17.0 ± 0.3 → 14.5 ± 0.3 pN** → AEBP2 destabilises PRC2–chromatin engagement | A |
| 4.1.16 † | JARID2 | **No significant change** in cluster distribution or transition force | A |
| 4.1.17 † | Compaction (negative-stain EM) | PRC2 condenses 12-mer and 4-mer arrays in a concentration-dependent manner (1:1 and 1:4 array:PRC2, 150 nM array / 600 nM PRC2); 95–371 masks per condition | A |
| 4.1.18 | In silico bridging populations | Reweighted: **Nuc 1–3 = 18%**, Nuc 1–2 = 4%, Nuc 1–4 = 7% (30 bp linker tetranucleosome, zero tension) | **B (model)** |

> **Correction to a prior note.** The statement that "the majority of nucleosomes are bound at 500 nM" is
> **not a measurement in Leicher 2020** — it is the authors' justification for their working concentration,
> cited to Wang 2017 NSMB (ref 7). The measured concentration dependence in Leicher 2020 is row 4.1.8
> (bridging saturates ~100 nM). The 21.4 ± 1.0 vs 16.6 ± 0.3 pN comparison is correct, but note 16.6 pN is
> the **C2 (mononucleosome) cluster with PRC2 present**, not a PRC2-free baseline (which is 16.0–16.2 pN).

### 4.2 AFM — Heenan et al. 2020 NAR †

| # | Quantity | Value | Conf |
|---|---|---|---|
| 4.2.1 † | DNA bend angle at bound PRC2 | **84 ± 3°** (mean ± SEM, N<sub>protein</sub> = 331) vs **29.3 ± 0.4°** for naked DNA (N<sub>segments</sub> = 2844; WLC prediction 28.9°) — a **3-fold** increase, tangent vectors 20 nm apart | A |
| 4.2.2 † | DNA wrapping | **No** evidence of wrapping — no correlation between higher bend angle and shortened contour length | A |
| 4.2.3 † | Naked DNA persistence length (internal control) | **49.1 ± 0.4 nm** (N = 640), matching the canonical 50 nm → the surface preparation is properly 2D-equilibrated | A |
| 4.2.4 † | Loop anchoring | Loops "typically anchored by **two or more** PRC2 molecules"; volumetric analysis resolves two peaks at **400 ± 200 nm³** and **1300 ± 600 nm³** (monomer and dimer/higher, N = 63/64), plus dimers seen without DNA | A |
| 4.2.5 † | Compaction vs occupancy | Mass-weighted area decreases monotonically with PRC2–DNA complex volume (N = 300) → compaction scales with number of bound PRC2 | A |

### 4.3 Cryo-EM geometry — Poepsel, Kasinath & Nogales 2018 NSMB (PMC5805599)

| # | Quantity | Value | Conf |
|---|---|---|---|
| 4.3.1 † | Bridging geometry | A **single** PRC2 engages **two** nucleosomes: the H3K27me3-bearing nucleosome's H3 tail into the EED aromatic cage, the substrate nucleosome's H3 tail into the EZH2 SET active site. Both contacts are made **via nucleosomal DNA near the DNA exit point**, coinciding with where H3 emerges | A |
| 4.3.2 † | Linker DNA in the 35 bp construct | A 35 bp DNA model fits the linker density with a **bend angle of 52°**; nucleosome cores are not displaced from their positioning sequences | A |
| 4.3.3 † | **Linker-length tolerance** | Structures solved for DiNcl**30**, DiNcl**35**, DiNcl**40**. PRC2 accommodates ±half a helical turn by **flipping the modified nucleosome ~180°** (DiNcl30) so the other H3 copy presents its tail, rather than by changing the PRC2 contacts | A |
| 4.3.4 † | Contacts | EZH2 SBD (basic patches 16–20 RKRVK, 27–30 RQLKR) contacts the **minor groove of the upper DNA gyre** of the modified nucleosome; EED N-terminal stretch + WD40 lateral surface contact both gyres, with orientation varying across 3D classes | A |

> This structure is the molecular basis of Leicher's C3 (Nuc 1–2) cluster. Leicher's C4 (Nuc 1–3) mode has
> **no** experimental structure — only the coarse-grained docking model (row 4.1.18, grade B).

### 4.4 Simulation of DNA looping — Leicher/Lin 2021 NAR (PMC8216278)

Mostly computational. Experimental part: 10 kbp bare DNA, same C-Trap and same buffer as §4.1
(200 mM KCl, room temperature, 500 nM PRC2, 0.1 µm/s).

| # | Quantity | Value | Conf |
|---|---|---|---|
| 4.4.1 † | PRC2 4-mer core on bare DNA | **73 rupture events** detected across the dataset (force and ΔL₀ distributions in Fig 6B) | A (experiment); **values in figure image only** |
| 4.4.2 † | PRC2 5-mer (+AEBP2) on bare DNA | **No rupture events detectable at all** → AEBP2 abolishes the stable DNA-sequestering mode. Consistent with 4.1.15 | A |
| 4.4.3 | Computed binding free energy, PRC2 on 147 bp 601 | PRC2b (C2-bound, AEBP2-like) **19.1 ± 0.9 k<sub>B</sub>T** vs experimental EMSA 17.2–17.8 k<sub>B</sub>T; PRC2f (C2-free, core-like) **24.7 ± 1.1 k<sub>B</sub>T** | **B (coarse-grained MD)** |
| 4.4.4 | Cooperative looping | ΔF(open→looped, 294 bp DNA): bare +4.7 k<sub>B</sub>T; one PRC2f +0.6; **two PRC2f −10.8 k<sub>B</sub>T**. One PRC2b +2.7, two PRC2b −2.2. Non-additive → cooperativity via DNA-mediated allostery, not protein–protein contact. Forcing a tight dimer (PRC2d) **abolishes** the cooperativity | **B (model)** |

---

## 5. Effect of H3K27me3 and other marks on BINDING (not catalysis)

This is the weakest-evidenced of the seven questions and the answer is largely negative.

| # | Mark | Effect on **binding affinity** | Technique / conditions | Paper | Conf |
|---|---|---|---|---|---|
| 5.1 † | **H3K27me3** on mononucleosomes and trinucleosomes | "PRC2 bound all of these nucleosome variants with **similar affinity**… H3K27me3-nucleosomes bound only **slightly more tightly**" than the 280 nM unmodified. Authors' own conclusion: "histone marks contribute little to the targeting of PRC2 to nucleosomes" | EMSA, 100 mM KCl, 30 °C, n = 3 | Wang 2017 NSMB, Results + Fig 2b,c | A (qualitative); **absolute K<sub>d</sub> in figure image only** |
| 5.2 † | **H3K4me3** | affinity "only marginally weakened"; catalysis separately reduced 4-fold | same | Wang 2017 NSMB, Fig 2c | A (qualitative) |
| 5.3 † | **H3K27M** on nucleosomes | bound "only slightly more tightly" than unmodified — i.e. **no** high-affinity sequestration at the level of equilibrium binding | same | Wang 2017 NSMB, Fig 2b,c | A (qualitative) |
| 5.4 † | H3K27M, **mechanical** | By contrast, H3K27M **does** stabilise the complex mechanically: C4 rupture 19.8 → 26.1 pN at 100% K27M (row 4.1.14). Grade A, but this is a force, not an affinity |  | Leicher 2020 | A |
| 5.5 † | H3K27me3, **mechanical** | 50% H3K27me3 **destabilises** (C4 21.4 → 16.9 pN); 100% restores baseline (row 4.1.13). Non-monotonic | | Leicher 2020 | A |
| 5.6 | H3K4me3 / H3K36me3 / tailless H3Δ38 | "PRC2 bound… with **similar affinity**", with or without JARID2-K116me3; "the histone H3 tail contributes minimally to the affinity of PRC2 for nucleosomes" | EMSA, 50 mM KCl, room temp, 187 bp nucleosome, n = 3 | Cookis 2025 NSMB (PMC11832421) | A |
| 5.7 | H3Kc36me3 | nucleosome affinity **indistinguishable** from unmodified, while H3K27me1/me3 catalysis is **>10-fold inhibited** | EMSA, 45 nM 6-FAM mononucleosome, n = 3 (buffer cited to Choi 2017) | Finogenova 2020 eLife (PMC7725500), Fig 3A,B | A |
| 5.8 | PRC2 interface mutants on nucleosomes | PRC2<sup>CXC>A</sup> or <sup>EED>A</sup> 2–3× weaker; double mutant ~5× weaker, still sub-µM. WT stated only as "mid-nanomolar" | EMSA, 215 bp mononucleosome, n = 3 | Finogenova 2020, Fig 2A,B | A (folds); **X for absolute WT K<sub>d</sub> — never printed** |
| 5.9 | **H2AK119ub1** | JARID2(1–530) **alone** (not in PRC2): **2.8 ± 0.5 µM** on H2AK119ub1 nucleosome vs **5.0 ± 0.7 µM** unmodified vs 5.1 ± 1 µM ub-I44A → ~2×, driven by k<sub>on</sub> | **BLI** (Octet, streptavidin), 150 mM NaCl, **25 °C**, n = 3 | Cooper 2016 Nat Commun 7:13661 (PMC5133711), Fig 4b | A |
| 5.10 | H2AK119ub1 in assembled PRC2.2 | Cryo-EM + HMTase only; JARID2-K116me3→EED is the primary allosteric driver with a "modest" contribution from JARID2/AEBP2 contacting ubiquitin. **No affinity measurement of any kind** | — | Kasinath 2021 Science 371:eabc3393 | **X** |
| 5.11 | **EED–peptide affinities** (a different quantity — free peptide to isolated subunit) | H3K27me3 **41.7 ± 2.8 µM**; JARID2-K116me3 8.07 ± 0.49 µM; PALI1-K1241me3 7.49 ± 0.59, me2 19.0 ± 2.3, me1 241 ± 38 µM | fluorescence-anisotropy displacement, 100 mM KCl, 30 °C, n = 3 | Zhang 2021 Nat Commun, Fig 3a/5b | A |
| 5.12 | EED–trimethyl-lysine panel | H3K9me3, H4K20me3, H1K26me3, H3K27me3 all **10–45 µM**; ~4× weaker per methyl removed; H3K4me3/K36me3/K79me3 not appreciably bound | competition fluorescence + ITC, 150 mM NaCl, 20 °C | Margueron 2009 Nature 461:762 | A for the range; per-peptide values in Suppl Table S1, **not retrieved** |

> **Bottom line for Q3:** there is, as far as I can find, **no direct quantitative measurement of a
> K<sub>d</sub> shift caused by H3K27me3 on a nucleosome**. Three independent papers (Wang 2017, Cookis 2025,
> Finogenova 2020) report that histone marks have little or no effect on PRC2 **binding** while having large
> effects on **catalysis**. Every H3K27me3 number in the literature is either (a) a free-peptide K<sub>d</sub>
> to isolated EED (~40 µM) or (b) an allosteric stimulation of k<sub>cat</sub>. Leicher 2020's rupture forces
> are the only mark-dependent mechanical measurement, and they are non-monotonic in methylation level.

---

## 6. Accessory subunits — summary of the affinity effects

| Subunit | Effect on **free DNA** | Effect on **nucleosomes** | Effect on **residence time** | Best source |
|---|---|---|---|---|
| **PHF1 (PCL1)** | 64 → **30 nM** (2.1×) | 163 → **66 nM** (2.5×) | τ_off,2 chromatin 2.28 → **6.48 s**; DNA 4.45 → **10.38 s**; via k<sub>off</sub> only. WH-domain mutant (R304E/S307E/N310E/K323E/K324E) reverts to core PRC2 | Choi 2017 † |
| **MTF2 (PCL2)** | holo PRC2.1<sup>MTF2</sup> **6.3 nM** on 100 bp CGI DNA; isolated cassette 2.1 µM by ITC | strongly increases dinucleosome affinity; prefers (GCN)ₙ linkers | not measured | Chen 2020; Li 2017; Kwon 2026 |
| **PHF19 (PCL3)** | — | — | not measured; tightest G4-RNA binder of the PCLs (34.0 nM) | Zhang 2019 |
| **AEBP2 (short)** | **confers** CpG-methyl specificity (11.3 nM mCpG, ~50× over CpG); core PRC2 > 4 µM → **45.5 nM** with AEBP2-S on CpG24 | present in the 280 nM 5-mer measurement | **Destabilises mechanically**: abolishes free-DNA sequestration (C1), lowers transition force 17.0 → 14.5 pN; abolishes detectable ruptures on bare DNA entirely | Wang 2017 †; Mucha 2025; Leicher 2020 †; Leicher 2021 † |
| **AEBP2 (long isoform)** | **auto-inhibitory**: K<sub>d</sub> > 4 µM on CpG24; removing/neutralising its two acidic tracts raises PRC2.2 affinity for 46 bp CpG DNA **>15-fold** but changes NCP182 K<sub>d</sub> ≤2-fold | 24 nM for PRC2.2 on NCP182 | — | Mucha 2025 EMBO J |
| **JARID2** | **no effect**: "JARID2 and EZH1 do not alter PRC2–nucleosome binding or RNA inhibition" | no effect | **No significant change** in optical-tweezers cluster distribution or transition force | Wang 2017 †; Leicher 2020 † |
| **EPOP** | **weakens** DNA binding 2–3× (3.4 → 9.2 nM) | **weakens** mononucleosome binding (122.9 → 295.4 nM); but **increases** dinucleosome affinity over core, preferring (GA)ₙ linkers | — | Gong 2026; Kwon 2026 |
| **PALI1** | core > 4 µM → **155 ± 26 nM** (>20×), **not** CpG-specific | NCP182 330 → **19.0 nM**; NCP147 2560 → **48.5 nM** | — | Zhang 2021 |

Two negative controls worth keeping: Son 2013 (Genes Dev) and Lee 2018 (Mol Cell) both report
JARID2/AEBP2 effects on nucleosome binding **qualitatively only** — no K<sub>d</sub> anywhere in either
paper. Do not let a review turn "markedly enhanced" into a number.

---

## 7. RNA

### 7.1 PRC2 affinity for RNA

| # | Ligand | K<sub>d</sub> | Technique | Enzyme | Salt | tRNA? | Temp | n | Paper | Conf |
|---|---|---|---|---|---|---|---|---|---|---|
| 7.1.1 | **G4 24 = (UUAGGG)₄** | **129 ± 6 nM** (Hill 0.92) | FA, 5 nM probe | core PRC2 4-mer | **200 mM KCl** | **yes**, 0.1 mg/ml | 30 °C | 3 | Zhang 2019 NSMB (PMC6736635) Table 1 | A |
| 7.1.2 | G4 24 | 54.6 ± 3.8 / **34.0 ± 1.8** / 82.4 ± 9.8 / 40.9 ± 3.9 nM | FA | +AEBP2 / **+PHF19** / +MTF2 / +MTF2+EPOP | 200 mM KCl | yes | 30 °C | 3 | Zhang 2019 Table 1 | A |
| 7.1.3 | **(UGAGUG)₄** size-matched non-G-tract control | **no detectable binding**, all complexes | FA | all above | 200 mM KCl | yes | 30 °C | 3 | Zhang 2019 Table 1 | A |
| 7.1.4 | G4 24, **HMTase-mimicking buffer** | **19.2 ± 2.7 nM**; +JARID2-K116me3 **46.7 ± 6.9 nM**; +H3K27me3 22.5 ± 3.1 nM | FA | PRC2-AEBP2 | **100 mM KCl**, no Mg/Zn | **no** | 30 °C | 11 / 6 / 4 | Zhang 2019 Table 3 | A |
| 7.1.5 | **(G4A4)₄** (40 nt) | **16.7 ± 1.2 nM** | FA, 20 nM probe | PRC2c (EZH2–EED–SUZ12<sup>VEFS</sup>) | **40 mM KCl** + 400 µM SAM | no | 20 °C | 3 | Beltran 2019 NSMB | A |
| 7.1.6 | PIM1 mRNA G4 (24 nt) | **22.5 ± 1.8 nM** | FA | PRC2c | 40 mM KCl | no | 20 °C | 3 | Beltran 2019 | A |
| 7.1.7 | K⁺ vs Li⁺ | binds in 100 mM KCl, **no significant binding in 100 mM LiCl** → the G-quadruplex fold itself is required | FA | PRC2c | 100 mM KCl / LiCl | no | 20 °C | 3 | Beltran 2019 Fig 3a,b | A (qualitative) |
| 7.1.8 | **r(G₃A₂)₄ G4 vs salt** | **2.3 ± 0.35** (10 mM KCl), **4.4 ± 0.34** (25), **1.4 ± 0.15** (100), **6.0 ± 0.93 nM** (200). log K<sub>d</sub> vs log[KCl] slope **m ≈ 0 ± 0.34** — essentially **salt-independent** | FP, 5 nM probe | PRC2 5-mer / 5me | 10–200 mM KCl | no | 25 °C | ≥3 | Hemphill 2023 PNAS Table 1 | A |
| 7.1.9 | r(GGAA)₁₀ | **78 ± 12 nM** (25 mM KCl) | FP | PRC2 5-mer | 25 mM KCl | no | 25 °C | ≥3 | Hemphill 2023 Table 1 | A |
| 7.1.10 | TERRA (50 nt, 4× UUAGGG) | **18 ± 9 nM**; TERRA<sup>mut</sup> (G-tracts removed) **66 ± 15 nM**; Poly(U) 28 ± 4; Poly(C) 28 ± 8 nM | EMSA | human 6-mer (+AEBP2 210–503, +JARID2 119–450) | 100 mM KCl | yes, 0.05 mg/ml | 30 °C | — | Song 2026 Nat Commun (PMC13316113) | A |
| 7.1.11 | HOTAIR 1–400 | 136 ± 22 nM (5-mer); 255 ± 3 nM (4-mer, Hill 1.7) | EMSA | human PRC2 | **buffer not in the accessible text** (Online Methods absent from PMC) | unknown | unknown | 2 | Davidovich 2013 NSMB (PMC3823624) | A (value) / **X (conditions)** |
| 7.1.12 | MBP mRNA 1–300 (non-cognate) | **110 ± 10 nM** (4-mer) — the promiscuity result | EMSA | human PRC2 | as above, unknown | unknown | unknown | 2 | Davidovich 2013 | A / X |
| 7.1.13 | RepA I–V (434 nt) vs MBP 1–300 vs P4–P6 | **14.6 ± 5.3** / **245 ± 77** / **> 5000 nM** | EMSA | human PRC2 5-mer | as 2013 | — | — | 3 | Davidovich 2015 Mol Cell (PMC4320675) Fig 1A | A |
| 7.1.14 | Same three RNAs, **alternative buffer** | all bound with **similar** K<sub>d</sub> regardless of length → specificity collapses | EMSA on TBE gel | PRC2 5-mer | 100 mM KCl, 2.5 mM MgCl₂, 0.1 mM ZnCl₂, **+0.1 mg/ml tRNA** | yes | RT fold | 3 | Davidovich 2015 Fig 3B | A |
| 7.1.15 | RepA I–IV / hHOTAIR 1–300 / Tsix I–IV / MBP / P4–P6 | **81** / **93** / **321** / **≫1000** / **≫1000 nM** | double-filter binding, 2 nM probe | mouse PRC2 4-mer | **100 mM NaCl**, 5 mM MgCl₂ | **yes**, 50 ng/µl | **30 °C** | — | Cifuentes-Rojas 2014 Mol Cell (PMC4107928) Fig 2B,C | A |
| 7.1.16 | Isolated subunits | EZH2 alone **28–32 nM**; SUZ12 alone 86–235 nM; EED **≫1000 nM**; JARID2 > 250 nM | double-filter binding | individual FLAG-subunits | 100 mM NaCl | yes | 30 °C | — | Cifuentes-Rojas 2014 | A |
| 7.1.17 | Endogenous CLIP-derived G-rich motifs | Aebp2 93.3 ± 4.8, HspA5 11.3 ± 0.3, Fbxo5 28.5 ± 1.4, Cbx7 35.3 ± 4.7, Mcm2 77.3 ± 2.8, Erlec1 12.5 ± 0.6 nM. **G→C mutation: Fbxo5 → 2058 ± 1274; Mcm2 → 3258 ± 2133 nM** | EMSA | mouse PRC2 | conditions by cross-reference to Cifuentes-Rojas 2014 | — | 30 °C | — | Rosenberg 2021 NSMB (PMC8050941) Fig 5 | A (K<sub>d</sub>) / B (conditions) |
| 7.1.18 | **poly(A)₄₀** | The consistent negative across every lab: inert as competitor, no inhibition, lowest affinity of the 40-mer homopolymers | multiple | multiple | multiple | — | — | — | Wang 2017 NSMB †, Long 2017, Song 2026 | A |

> **The tRNA variable.** Same lab, same G4-class ligand: ~1.4 nM at 10 mM KCl without tRNA (Hemphill)
> vs ~129 nM at 200 mM KCl with 0.1 mg/ml tRNA (Zhang). Treat yeast tRNA as a first-order experimental
> variable, not a buffer detail. Wang 2017 NSMB states outright that without tRNA competitor the RNA
> K<sub>d,app</sub> values are "in the **picomolar** range" (no number given).

### 7.2 RNA vs DNA/nucleosome competition

| # | Observation | Value | Paper | Conf |
|---|---|---|---|---|
| 7.2.1 † | RNA strips PRC2 from DNA; DNA cannot strip PRC2 from RNA | (GGAA)₁₀ competes PRC2 off A488-dsDNA with t½ = 91 ± 5 s; **10 µM dsDNA fails** to compete PRC2 off A488-RNA; poly(A)₄₀ inert | Wang 2017 NSMB, Fig 5e,f | A |
| 7.2.2 † | Sub-stoichiometric eviction from chromatin | Pre-formed PRC2–trinucleosome complexes **fully disrupted at PRC2:RNA ≈ 2:1** by (GGAA)₁₀. Same for H3K27me3-trinucleosome and for PRC2–JARID2 6-mer. Poly(A)₄₀ fails | Wang 2017 NSMB, Figs 1f, 2d, 6e | A |
| 7.2.3 | G4 RNA blocks nucleosome binding | 500 nM (G4A4)₄ or PIM1 G4 **blocks** PRC2c binding to 147 bp nucleosomes (K<sub>d</sub> 25.9 nM) and to H3K27M nucleosomes (29.9 nM); non-G4 PIM1 control has no effect | Beltran 2019, Fig 3c,d | A |
| 7.2.4 | Sequence-selective eviction | TERRA, TERRA<sup>mut</sup>, pUG-fold all strip ³²P-trinucleosome from PRC2; Poly(U), Poly(C), P10, Poly(A) do **not** | Song 2026, Fig 5a | A |
| 7.2.5 | Shared binding surface | EZH2 GR-loop and K-loop basic patches bind **both** G4 RNA and nucleosome: mutating them costs 5–12× in RNA affinity and 2–31× in nucleosome affinity (PRC2c mononucleosome 0.4 → 13.8 µM for the GR+K double mutant) | Houslay 2026 **bioRxiv preprint** | A (preprint) |
| 7.2.6 | RNA depletion increases nucleosome binding | RNase-A treatment of ESC nuclear extract increases PRC2 binding to both 147 bp and 183/185 bp nucleosomes, independent of PCL2/AEBP2/JARID2 | Beltran 2019, Fig 3g | A |
| 7.2.7 | **Rate constants** | See §3.2 rows 3.2.2–3.2.7 — the direct-transfer constants are the only measured rate constants for RNA↔DNA exchange on PRC2 | Hemphill 2023 | A |
| 7.2.8 | RNA effect on single-molecule residence time | **Not measured experimentally anywhere.** Hemphill 2023 simulates it | — | **X** (experiment) / B (simulation) |

### 7.3 Catalysis, kept separate (these are NOT affinities)

| Quantity | Value | Paper |
|---|---|---|
| IC₅₀ of RNA for PRC2 HMTase | **~500 nM** for TERRA, TERRA<sup>mut</sup> and pUG-fold alike (at 600 nM PRC2); **no** inhibition by Poly(U)/Poly(C)/P10/Poly(A) | Song 2026 |
| Explicit K<sub>d</sub>≠IC₅₀ caveat | K<sub>d</sub>(TERRA) = 18 nM but IC₅₀ ≈ 500 nM, "because the inhibition of PRC2 is determined not only by the K<sub>d</sub>… but by the RNA-PRC2 stoichiometry, given that 600 nM PRC2 was used" | Song 2026, Results |
| RNA inhibits methylation of **DNA-free H3** | G4 256 RNA inhibits PRC2 toward free H3.1, not only toward nucleosomes → part of the catalytic inhibition is independent of DNA competition | Zhang 2019, Fig 5a |
| RNA does **not** inhibit EZH2 automethylation | automethylation persists while H3K27 methylation is abolished → RNA is not an active-site inhibitor | Wang 2017 NSMB †, Fig 1d,e |
| Allosteric activation weakens G4 binding | PRC2c–G4A4 K<sub>d</sub>: 8.0 ± 1.1 nM alone → 31.9 ± 12.9 nM (+SAM) → 59.2 ± 4.9 nM (+SAM+JARID2-K116me3+H3K27M). SAM-binding mutant N688A unchanged | Houslay 2026 preprint |
| PRC2 catalysis, linker-length "Goldilocks" | HMTase on dinucleosomes is **non-monotonic** in linker length — 50 bp stimulates more than 10 or 100 bp, even though **binding** increases monotonically to 100 bp. Binding optimum ≠ catalytic optimum | Wang 2017 NSMB †, Fig 4c,d |
| H3Kc36me3 | binding unchanged, H3K27me1/me3 formation **>10-fold** inhibited | Finogenova 2020 |

---

## 8. Does PRC2 slide or diffuse on DNA?

**Answer: no measurement exists. Grade X.**

- No published in vitro single-molecule measurement of a PRC2 1D diffusion coefficient on DNA or chromatin
  (no DNA-curtain, no sliding-trajectory, no MSD analysis). Searches across Europe PMC for
  `facilitated diffusion`, `one-dimensional diffusion`, `diffusion coefficient`, `sliding`, `target search`
  and `DNA curtain` combined with PRC2 return nothing primary.
- The **one direct experimental test** is a negative result: Leicher 2020 placed two well-separated dCas9
  roadblocks on bare DNA and reasoned that if PRC2 engaged two DNA segments and **slid** until it hit steric
  blocks, large force-induced transitions should appear. None did, despite dCas9 being stably bound. The
  authors conclude this "disfavors the nonspecific roadblock model" (row 4.1.11). That argues against
  long-range sliding as the mechanism behind the bridged states, but it was not designed to measure
  diffusion and does not exclude short-range sliding. Grade A as a negative result.
- Heenan 2020 found PRC2 binds **promiscuously along** 2.5 kb DNA with no positional preference for an
  embedded 76% GC island, with a localisation assay verified to 2 nm by a BspMI positive control. Consistent
  with either fast sliding or sequence-blind direct binding; it does not distinguish them.
- Choi 2017's λ_on (12–37 s at 1–2 nM) and k_on (1.4–4.4 × 10⁷ M⁻¹s⁻¹) describe the 3D search rate, not a
  1D diffusion constant. Note k_on near 10⁷–10⁸ M⁻¹s⁻¹ is at or above the typical 3D diffusion-limited
  range, which is sometimes taken as a hint of facilitated diffusion — but Choi makes no such claim and the
  k_on errors are 40–90% of the value.

---

## 9. Verification of the three priors I was asked to check

| Prior | Verdict |
|---|---|
| **Choi 2017** — "mid-nanomolar affinity and multisecond residence on free DNA, PHF1 prolonging 2–3×" | **Confirmed and now fully quantified** (full text obtained via nature.com; EPFL's green-OA record holds no file). DNA K<sub>d</sub>: PRC2 64 ± 4 nM, PHF1–PRC2 **30 ± 3 nM** (2.1×). Mononucleosome: 163 ± 13 → **66 ± 4 nM** (2.5×). Residence τ_off,2 on 601-DNA: 4.45 ± 1.90 → **10.38 ± 2.46 s** (2.3×); on 12-mer chromatin 2.28 ± 0.44 → **6.48 ± 0.52 s** (2.8×). Full Table 1 reproduced in §3.1. Three refinements to the prior: (i) the effect is entirely on k<sub>off</sub>, k<sub>on</sub> is unchanged; (ii) the kinetics are **biexponential** — 46–75% of events are sub-second; (iii) the salt is **50 mM NaCl**, which is low, and the EMSA probe (45 nM) is comparable to the K<sub>d</sub>, which the authors say may make these an **underestimate** of affinity. |
| **Wang 2017 NSMB** — "DNA dominates; linker length increases affinity; H3K27me3 has little effect" | **Confirmed, with numbers.** NCP (no linker) 40 ± 16 µM → mononucleosome with 60 bp linker 280 ± 19 nM → trinucleosome 23 ± 9 nM → dodecanucleosome 9 ± 2 nM; free 147 bp DNA tighter still. Dinucleosome with 10 bp linker 276 ± 16 nM, substantially weaker than 50/100 bp. H3K27me3, H3K27M and H3K4me3 nucleosomes all bound with "similar affinity" to unmodified. Conditions: 100 mM KCl, 30 °C, no tRNA, n = 3. One thing to correct if it is carried anywhere: **276 ± 16 nM is a dinucleosome with a 10 bp linker, not a mononucleosome**. |
| **Leicher 2020 PNAS** — "bridged states saturate 5–100 nM, majority of nucleosomes bound at 500 nM, rupture 21.4 ± 1.0 vs 16.6 ± 0.3 pN" | **Mostly confirmed, one correction.** Bridging (C4) population rises 5 → 100 nM and saturates: correct. 21.4 ± 1.0 pN (C4) vs 16.6 ± 0.3 pN (C2): correct, but 16.6 pN is the mononucleosome cluster **in the presence of PRC2**; the PRC2-free baseline is 16.0–16.2 ± 0.3 pN. **"Majority of nucleosomes bound at 500 nM" is not a measurement in this paper** — it is the authors' justification for their working concentration, cited to Wang 2017 (ref 7). Also add: the experiment is at **200 mM KCl**, the highest salt of any study in this file, which matters when comparing its conclusions to the 25–100 mM KCl affinity work. |

---

## 10. Conflicts and traps to carry forward

1. **Salt dependence is the dominant variable for DNA, and is absent for RNA.** Hemphill 2023 measured both
   on the same instrument: log K<sub>d</sub> vs log[KCl] slope is **~1.4 for 60 bp dsDNA** (5 nM at 10 mM KCl
   → unmeasurable at 200 mM) but **~0 for G4 RNA**. Any model that treats PRC2–DNA affinity as a constant
   across ionic conditions is wrong by two orders of magnitude over the physiological range.
2. **Linker DNA is the second dominant variable.** 40 µM (no linker) vs 280 nM (60 bp linker) in the same
   paper — a ~150-fold effect, larger than anything any accessory subunit does.
3. **Nucleosome arrays are not N× mononucleosomes.** Wang's array K<sub>d</sub> values were deliberately
   **not** normalised for binding-site number, and the Hill coefficient rises substantially from mono to
   arrays. Do not divide by 12.
4. **Binding optimum ≠ catalytic optimum.** Binding rises monotonically with linker length to 100 bp;
   HMTase activity peaks at 50 bp. Several papers slide between the two.
5. **Three incompatible "residence times."** See §3.3.
6. **"Trapping" vs "direct transfer"** for RNA-mediated eviction: superseded in favour of direct transfer
   (§3.2 note). Use k₋₁ and k<sub>θ</sub>, not the 91 s half-life.
7. **IC₅₀ is not K<sub>d</sub>.** Song 2026 states the gap explicitly: 18 nM vs ~500 nM for the same RNA,
   because inhibition is stoichiometry-limited at 600 nM enzyme.
8. **AEBP2 pulls in opposite directions depending on readout.** It confers mCpG specificity and tightens
   binding to short CpG DNA (Wang 2017, Mucha 2025) but **destabilises** PRC2–chromatin engagement
   mechanically and abolishes the DNA-sequestering mode (Leicher 2020/2021). These are not contradictory —
   affinity and binding geometry are different things — but a single "AEBP2 effect" number is meaningless.
9. **Isolated domains vs holo-complex.** PCL winged-helix domains alone bind DNA at 0.5–10 µM
   (Choi 2017, Li 2017). The same protein inside PRC2 contributes a 2–3× change against a 30–60 nM
   background. Do not transfer the isolated-domain K<sub>d</sub> into a holo-complex model.
10. **PCL sequence specificity is contested.** Choi 2017 (FP, 50 mM NaCl, 20 °C) finds none; Li 2017
    (ITC, 150 mM NaCl, 10 °C) finds a ~40-fold NCpGN preference. Choi explicitly rebuts Li in its Discussion.
11. **Preprint flag.** Houslay 2026 (G4 sequesters the EZH2 GR/K loops, monomeric PRC2c) directly contradicts
    Song 2026's RNA-induced-dimerisation model. Houslay is a bioRxiv preprint, not peer-reviewed.

---

## 11. What is still missing

**Papers I could not obtain:**

- **Wang X. et al. 2017 Mol Cell 65:1056, "Targeting of PRC2 to RNA by Short Repeats of Consecutive
  Guanines"** (PMID 28306504). Not in PMC. ScienceDirect and cell.com both return bot-detection blocks.
  This is the largest single gap: it holds the primary (G4A4)₄ / (GGAA)₁₀ / G-tract-length K<sub>d</sub>
  series that later papers cite. **No numbers from it appear in this file.** It is likely obtainable
  manually through a Queen's library session.
- **Schmitges 2011 Mol Cell 42:330** — no PMCID, no route. So there is **no verified H3K4me3/H3K36me3
  peptide K<sub>d</sub>** here, only Finogenova 2020's secondary summary (grade C).
- **Nekrasov 2005 EMBO Rep 6:348** (PMC1299286) — PDF-only, 403 to automated access. The main *Drosophila*
  PRC2 nucleosome-binding gap. The abstract quotes no K<sub>d</sub> either.
- **Jani 2019 PNAS 116:8295** — PDF-only; the "5-fold higher k<sub>cat</sub> for unmodified H3K36" is
  abstract-only (grade C).
- **Yuan 2012 Science 337:971** (nucleosome density activates PRC2) — abstract only, no PMCID, grade C.

**Numbers that exist only inside figure images** (would need the figure read visually or the source-data
files downloaded):

- Wang 2017 NSMB: Fig 2c (H3K27me3 / H3K27M / H3K4me3 mono- and trinucleosome K<sub>d</sub> values and their
  naked-DNA controls), Fig 3b (the full linker-length series), Fig 4b (50 and 100 bp dinucleosome K<sub>d</sub>),
  Fig 5a ((CG)₃₀ vs (TA)₃₀), Fig 7a (the summary affinity ladder). **These are the single most valuable
  remaining numbers for question 1 and question 3.** Source data are available with the paper online.
- Finogenova 2020 Fig 2B: the absolute WT PRC2–mononucleosome K<sub>d</sub> ("mid-nanomolar", never printed).
- Kwon/Reinberg 2026 PNAS Fig 4B/D/F: all nine dinucleosome K<sub>d</sub> values for core/EPOP/MTF2 ×
  R/GA/GCN linkers.
- Long 2017 eLife Fig 1D: the 40-mer homopolymer K<sub>d</sub> table including poly(A).
- Mucha 2025 Figs 1D, 6I, EV8; Leicher 2021 Fig 6B (rupture force and ΔL₀ distributions on bare DNA);
  Beltran 2019 Fig 3a,b (FA K<sub>d</sub> at 100 mM KCl); Cifuentes-Rojas 2014 Figs 4E/6F (IC₅₀);
  Long 2020 Ext Data Fig 3 (absolute k<sub>off</sub>).
- Margueron 2009 Suppl Table S1: per-peptide EED K<sub>d</sub>. (Zhang 2021's 41.7 ± 2.8 µM for
  H3K27me3–EED at 100 mM KCl / 30 °C is the best verified substitute.)

**Quantities nobody appears to have measured:**

- **A K<sub>d</sub> shift caused by H3K27me3 on a nucleosome.** Three papers say "no effect" qualitatively;
  none gives a paired number. This is the clean experiment the field is missing.
- **A K<sub>d</sub> for H3K27M nucleosomes** that would support the sequestration model. Wang 2017 says
  "slightly more tightly"; Beltran 2019 gets 29.9 nM for H3K27M vs 25.9 nM unmodified with PRC2c — i.e.
  no difference. Only the mechanical measurement (Leicher 2020) supports stabilisation.
- **Any 1D diffusion coefficient for PRC2** (§8).
- **Any SPR, BLI or ITC measurement of holo-PRC2 binding to a nucleosome.** The entire nucleosome-affinity
  literature is EMSA plus a handful of fluorescence-anisotropy experiments. The only BLI number in this
  file is for free JARID2(1–530), not PRC2.
- **Single-molecule residence time as a function of RNA** (7.2.8).
- **Drosophila PRC2 quantitative binding**, of any kind.
- **Any structure of the Nuc 1–3 bridged state.** It is the dominant bridging mode by force spectroscopy
  (78% of C4) and the mechanically strongest state (21.4 pN), but exists structurally only as a
  coarse-grained docking model.
