# Protein stability (half-life) data — provenance & methods

Companion to [stability_data_provenance.md](stability_data_provenance.md) (which leads with the **mRNA**
half-lives that are *baked* into the model). This file is the full **protein** half-life record. Same source:
the lab's **stability-catalog** (AMOS external-parameter arm), a curated aggregation of 24 published datasets /
120 arms; no value is a lab measurement; all are *degradation-only* estimates. Half-lives in hours; verified
against the model code 2026-10-19.

> **Headline:** the baseline model uses **no** catalog protein half-life directly. CDK-inhibitor and Skp2
> protein turnover is **cell-cycle-regulated** in the model (SCF–Skp2, KPC, APC/Cdh1), and EZH2 protein
> turnover is **fit** to EZH2-coupled targets. Measured protein half-lives enter only as a **sensitivity
> analysis** on the two short-lived/unregulated inhibitors p18 and p21 (`cdki_measured_protein`, default OFF),
> which still validated 30/32.

---

## 0. How one protein value was selected per gene (the ranking rule)

Identical machinery to the mRNA selection (`stability.py:_half_life`, `prefer="degradation_only"`): a **fixed,
deterministic preference ranking**, first qualifying arm wins, value = median of its replicates.

1. **Curated dataset order** (tissue-relevance first) for protein:
   Nakanoh 2026 ESC-derived **neural progenitor** (mouse, human) → Li 2021 cell-line turnover
   (RPE1 → U2OS → HCT116 → HEK293T) → Zecha 2018 HeLa → Fornasiero 2018 / Li 2025 brain tissue →
   Mathieson 2018 immune/neuron → Schwanhäusser 2011 / McShane 2016 NIH3T3 → Swovick 2021 fibroblast.
2. **Degradation-only** preferred: an arm reporting a *removal rate that includes growth dilution* is demoted
   (+500). This matters here — e.g. Zecha's p18 is **23.9 h** (pure degradation) vs **14.1 h** (turnover
   *incl.* cell doubling) from the *same* study; we take 23.9 h because the model adds division separately.
3. **Exclusions**: censored bounds (+1000); McShane **1-state** fits (+100; they measure a different quantity
   than the 2-state steady-state fit — visible in EZH2 below, where the 1-state 4.98 h sits far below the rest).
4. Lowest rank wins.

The protein winners differ by gene (unlike mRNA, where Serdar E14.5 won everything) because the top-ranked
neural-progenitor arm (Nakanoh) only measures a few of these genes — so the ranking falls through to the first
arm that has the gene: **p18, p27 → Zecha 2018 HeLa; p21, Skp2, Cyclin D1 → Li 2021.**

---

## 1. Chosen protein half-lives + how each is used

| Gene / symbol | chosen protein t½ | source arm | class | **use in model** |
|---|---|---|---|---|
| p18 / *Cdkn2c* | **23.87 h** | zecha2018:HeLa (n=19) | — | sensitivity variant only (`cdki_measured_protein`) |
| p21 / *Cdkn1a* | **0.98 h** | li2021:RPE1 (n=5) | short-lived | sensitivity variant only (`cdki_measured_protein`) |
| p27 / *Cdkn1b* | **44.47 h** | zecha2018:HeLa (n=15) | — | examined, **not used** (regulated clearance) |
| Skp2 / *Skp2* | **6.77 h** | li2021:HCT116 (n=8) | short-lived | examined, **not used** (regulated clearance) |
| CyclinD1 / *Ccnd1* | **0.42 h** | li2021:RPE1 (n=8) | short-lived | consistent w/ model short Cd; level not catalog-set |
| EZH2 / *Ezh2* | **57.47 h** | zecha2018:HeLa (n=8) | — | **not used** — EZH2 protein t½ is model-fit (~25 h effective) |
| p19 / *Cdkn2d* | **no data** | — | — | no protein arm exists; model value assumed (~0.4 h, reviewer) |

## 2. Complete per-gene arm list (what the other datasets said)

**p18 / CDKN2C** — chosen **23.9 h** (Zecha HeLa). Only cell-line turnover value; all other arms are
post-mitotic **brain/organ tissue** where p18 is long-lived: Li 2025 amygdala 28 h, thalamus 66 h, frontal
cortex 96 h, hippocampus 217 h, **cerebellum 228 h**, kidney 998 h; Schwanhäusser NIH3T3 189 h. (The dividing
HeLa value is the appropriate turnover for a cycling progenitor; the tissue values reflect no-dilution neurons.)

**p19 / CDKN2D** — **no protein measurement in the catalog.** Model value is a reviewer estimate (~0.4 h, p19^INK4d is a short-lived S-peaking INK4).

**p21 / CDKN1A** — chosen **0.98 h** (Li RPE1). Tight, all 5 arms short-lived: RPE1 0.98, U2OS 1.36, HCT116 1.52, Zecha HeLa 2.15 h. Robust short-lived CDKI.

**p27 / CDKN1B** — chosen **44.5 h** (Zecha HeLa). **Strongly context-dependent:** Nakanoh ESC-NP (mouse) **10.2 h**, Mathieson mouse-neuron 15.6 h, Zecha HeLa 44.5 h, brain tissue 55–143 h (Fornasiero cerebellum 84 h, Li 2025 cerebellum 96 h). The 4× spread between a neural-progenitor value (10 h) and HeLa (44 h) is exactly why we keep the model's **regulated** p27 clearance rather than one scalar.

**Skp2 / SKP2** — chosen **6.77 h** (Li HCT116). Range 3.3–13.7 h: McShane NIH3T3 3.3, Schwanhäusser 3.4, Zecha HeLa 5.6, Li HCT116 6.8, Li HEK293T 9.5, Li U2OS 13.7 h. Short-lived; kept as regulated (E2F-synthesis + APC/Cdh1) in the model.

**CyclinD1 / CCND1** — chosen **0.42 h** (Li RPE1). Cell lines tightly short: RPE1 0.42, HeLa 0.45, HCT116 0.47, U2OS 0.55, Schwanhäusser 1.31 h; tissue/fibroblast outliers 12.8 h (hamster), 33.9 h (kidney). Consistent short-lived Cyclin D1.

**EZH2 / EZH2** — chosen **57.5 h** (Zecha HeLa). Wide (12–57 h): McShane 2-state 11.6, Zecha-incl-doubling 19.0, **Nakanoh ESC-NP mouse 21.2 / human 39.2**, Schwanhäusser 22.0, Li RPE1 26.1, Zecha HeLa 57.5 h. **Not used** — the model's EZH2 turnover is fit; its ~25 h effective value matches the neural-progenitor arms (21–39 h) better than the HeLa chosen scalar.

## 3. Why measured protein half-lives are not imposed on the baseline model

The CDK inhibitors (p21, p27) and Skp2 are cleared by **cell-cycle-regulated proteolysis** in the model
(SCF–Skp2 for p21/p27; APC/Cdh1 for Skp2; KPC for p27), giving effective cycling-cell half-lives of ~1–2 h that
drive the Skp2–p27 commitment switch. The catalog's values are **degradation-only, from non-cycling or
steady-state cell lines**, and (for p27) are themselves cell-type-dependent (10 h in neural progenitors → 44 h
in HeLa). A single measured scalar would therefore mis-specify the regulated clearance. The two exceptions —
**p18** (an unregulated INK4 brake level) and **p21** (short-lived, where the measured 0.98 h ≈ the model's
~0.6 h) — are the only ones imposed, and only as a sensitivity check (still 30/32). EZH2 protein turnover is
pinned by EZH2-coupled dynamics (HU-arrest in-S boost, G0-withdrawal), not by the catalog.

## 4. Source citations (protein arms)

- **Zecha J, Meng C, Zolg DP, Samaras P, Wilhelm M, Kuster B (2018)** — HeLa pulsed-SILAC proteome turnover (p18, p27, EZH2). *Mol. Cell. Proteomics* 17:974–992. doi:10.1074/mcp.RA118.000583.
- **Li J, Cai Z, Vaites LP, … Gygi SP (2021)** — RPE1/U2OS/HCT116/HEK293T turnover, non-linear model (p21, Skp2, Cyclin D1). *Molecular Cell* 81:4722–4735.e5. doi:10.1016/j.molcel.2021.09.015.
- **Nakanoh S, et al. (2026)** — ESC-derived neural-progenitor turnover (p27, EZH2 alternatives). *Developmental Cell* 61:1867–1882.e8. doi:10.1016/j.devcel.2026.07.014.
- Alternative/comparison arms: Li W et al. (2025) *Cell* 188:2267 (brain/organ tissue); Fornasiero et al. (2018) *Nat. Commun.* 9:4230 (cerebellum); Mathieson et al. (2018) *Nat. Commun.* 9:689 (neuron/immune); Schwanhäusser et al. (2011) *Nature* 473:337 (NIH3T3); McShane et al. (2016) *Cell* 167:803 (NIH3T3/RPE1); Dörrbaum et al. (2018) *eLife* 7:e34202; Swovick et al. (2021) *MCP* 20:100041.

## 5. Draft methods paragraph (protein)

> **Protein half-lives.** Protein turnover values were compiled from the same curated catalog using its
> degradation-only preference ranking (neural-progenitor arms first, then cell-line turnover proteomes, then
> brain tissue), adopting the degradation constant rather than the growth-inclusive removal rate because the
> model supplies its own dilution term. The catalog yielded p18 23.9 h and p27 44.5 h (Zecha et al., 2018,
> HeLa), p21 0.98 h, Skp2 6.8 h and Cyclin D1 0.42 h (Li et al., 2021, RPE1/HCT116), and EZH2 57 h (Zecha et
> al., 2018); no protein half-life was available for p19 (*Cdkn2d*). For the baseline model these were not
> imposed directly: p21, p27 and Skp2 are cleared by the model's cell-cycle-regulated proteolysis (SCF–Skp2,
> KPC, APC/Cdh1), whose effective cycling-cell half-lives (~1–2 h) differ from the degradation-only
> measurements made in non-cycling cell lines — p27 itself ranged from 10 h in ESC-derived neural progenitors
> (Nakanoh et al., 2026) to 44 h in HeLa — and EZH2 turnover was instead constrained by EZH2-coupled cell-cycle
> observations. Measured protein half-lives were imposed only as a sensitivity analysis on the two
> short-lived/unregulated inhibitors p18 and p21, which left the validation suite unchanged (30/32 targets).
