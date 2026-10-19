# Paper table + selection paragraph — half-lives

Drop-in table and paragraph for the manuscript. Full provenance in
[stability_data_provenance.md](stability_data_provenance.md) and
[protein_stability_provenance.md](protein_stability_provenance.md).

---

## Table — mRNA and protein half-lives used to parameterise the model

| Gene (symbol) | mRNA t½ (h) | protein t½ (h) | source (mRNA / protein) | in model |
|---|:--:|:--:|---|---|
| Cyclin D1 (*Ccnd1*) | 3.5 | 0.4 | Serdar 2025 / Li 2021 | mRNA |
| p18 (*Cdkn2c*) | 3.7 | 23.9 | Serdar 2025 / Zecha 2018 | mRNA; protein (sensitivity) |
| p19 (*Cdkn2d*) | 3.6 | n.d. | Serdar 2025 / — | mRNA |
| p21 (*Cdkn1a*) | 5.2 | 1.0 | Serdar 2025 / Li 2021 | mRNA; protein (sensitivity) |
| p27 (*Cdkn1b*) | 3.6 | 44.5 | Serdar 2025 / Zecha 2018 | mRNA |
| Skp2 (*Skp2*) | 3.4 | 6.8 | Serdar 2025 / Li 2021 | mRNA |
| EZH2 (*Ezh2*) | 4.2 | 57 | Serdar 2025 / Zecha 2018 | mRNA |

**Caption.** Messenger-RNA half-lives (degradation-only) parameterise the two-step
transcription→mRNA→protein modules and are used throughout. Protein half-lives are the measured values for
reference; they are imposed only as a sensitivity analysis for p18 and p21 (the short-lived/unregulated CDK
inhibitors), whereas p27 and Skp2 are cleared by the model's cell-cycle-regulated proteolysis (SCF–Skp2, KPC,
APC/Cdh1) and EZH2 protein turnover is fit to EZH2-coupled data. mRNA sources: Serdar et al. (2025), mouse
embryonic cortex E14.5. Protein sources: Li et al. (2021), RPE1/HCT116; Zecha et al. (2018), HeLa. n.d., no
protein measurement available.

---

## Short paragraph — how the values were selected

> Transcript and protein half-lives were drawn from a curated compilation of published mammalian turnover data
> (24 datasets, 120 experimental arms; degradation-only estimates). Rather than averaging, a single value per
> gene was chosen by a fixed preference ranking that favoured, first, the most physiologically relevant
> dataset — for mRNA, mouse embryonic cortex (Serdar et al., 2025), the closest available neural-progenitor
> arm, which supplied all six transcripts; for protein, neural-progenitor and well-characterised cell-line
> turnover proteomes ahead of post-mitotic tissue — and, second, degradation-only measurements, because the
> model carries its own division/dilution term; censored bounds and non-comparable model fits were excluded.
> The selected mRNA values were representative of the relevant datasets (within ~2-fold of an independent
> iPSC-derived neural-progenitor measurement; Mufteev et al., 2026), and the wide across-dataset ranges
> reflect actinomycin-D transcription-arrest and whole-tissue measurements that the ranking deprioritises. As
> the transcription→mRNA→protein cascades were parameterised to preserve steady-state abundance, the model's
> calibration is insensitive to these half-lives (identical validation with normalized values), which set only
> the kinetics of the transcriptional-burst response.
