# All measured PROTEIN half-lives (catalog) and what the model used

Every arm in the stability-catalog for each gene (hours), sorted low→high. **✓ = catalog's chosen arm**
(degradation-only preference ranking). The last column of each header says what the MODEL does with it.

## p18 / CDKN2C  —  chosen **23.87 h**  (zecha2018:HeLa)
*19 arms, 14.1–997.9 h (71× spread). Model: sensitivity variant (cdki_measured_protein) — 23.9 h*

| t½ (h) | dataset : cell-type / variant | |
|---:|---|:--:|
| 14.12 | zecha2018:HeLa/turnover_incl_cell_doubling |  |
| 23.87 | zecha2018:HeLa | ✓ |
| 27.66 | li2025:amygdala |  |
| 39.99 | li2025:liver |  |
| 50.56 | li2025:lung |  |
| 65.55 | li2025:spleen |  |
| 66.48 | li2025:thalamus |  |
| 78.92 | li2025:olfactory_bulb |  |
| 79.17 | li2025:heart |  |
| 95.71 | li2025:frontal_cortex |  |
| 97.05 | dorrbaum2018:glia_enriched |  |
| 103.82 | li2025:substantia_nigra |  |
| 125.41 | li2025:striatum |  |
| 134.39 | li2025:entorhinal_cortex |  |
| 189.49 | schwanhausser2011:NIH3T3 |  |
| 216.84 | li2025:hippocampus |  |
| 227.91 | li2025:cerebellum |  |
| 554.97 | li2025:gut |  |
| 997.85 | li2025:kidney |  |

## p19 / CDKN2D

**No protein measurement in the catalog.** Model: NO DATA — model assumes ~0.4 h (reviewer)

## p21 / CDKN1A  —  chosen **0.98 h**  (li2021:RPE1)
*5 arms, 1.0–2.1 h (2× spread). Model: sensitivity variant (cdki_measured_protein) — 0.98 h*

| t½ (h) | dataset : cell-type / variant | |
|---:|---|:--:|
| 0.98 | li2021:RPE1/non-linear model | ✓ |
| 1.36 | li2021:U2OS/non-linear model |  |
| 1.52 | li2021:HCT116/non-linear model |  |
| 2.01 | zecha2018:HeLa/turnover_incl_cell_doubling |  |
| 2.15 | zecha2018:HeLa |  |

## p27 / CDKN1B  —  chosen **44.47 h**  (zecha2018:HeLa)
*15 arms, 10.2–143.5 h (14× spread). Model: NOT used — regulated clearance (SCF-Skp2/KPC), ~1-2 h effective*

| t½ (h) | dataset : cell-type / variant | |
|---:|---|:--:|
| 10.22 | nakanoh2026:ESC_derived_NP_mouse |  |
| 15.57 | mathieson2018:mouse_neuron |  |
| 23.10 | zecha2018:HeLa/turnover_incl_cell_doubling |  |
| 41.71 | li2025:spleen |  |
| 44.47 | zecha2018:HeLa | ✓ |
| 47.74 | mathieson2018:B_cell |  |
| 54.67 | li2025:olfactory_bulb |  |
| 67.55 | li2025:frontal_cortex |  |
| 84.42 | fornasiero2018:cerebellum_homogenate |  |
| 85.24 | li2025:entorhinal_cortex |  |
| 95.77 | li2025:cerebellum |  |
| 96.19 | li2025:heart |  |
| 98.35 | li2025:substantia_nigra |  |
| 123.11 | li2025:kidney |  |
| 143.49 | li2025:liver |  |

## Skp2 / SKP2  —  chosen **6.77 h**  (li2021:HCT116)
*8 arms, 3.3–6.8 h (2× spread). Model: NOT used — regulated clearance (E2F synth + APC/Cdh1)*

| t½ (h) | dataset : cell-type / variant | |
|---:|---|:--:|
| 3.30 | mcshane2016:NIH3T3/2state |  |
| 3.38 | schwanhausser2011:NIH3T3 |  |
| 3.57 | mcshane2016:NIH3T3/1state |  |
| 4.62 | zecha2018:HeLa/turnover_incl_cell_doubling |  |
| 5.59 | zecha2018:HeLa |  |
| 6.77 | li2021:HCT116/non-linear model | ✓ |
| 9.47 | li2021:HEK293T/linear model |  |
| 13.73 | li2021:U2OS/linear model |  |

## CyclinD1 / CCND1  —  chosen **0.42 h**  (li2021:RPE1)
*8 arms, 0.4–33.9 h (81× spread). Model: NOT used as a scalar — model Cd is short-lived, consistent*

| t½ (h) | dataset : cell-type / variant | |
|---:|---|:--:|
| 0.42 | li2021:RPE1/non-linear model | ✓ |
| 0.44 | zecha2018:HeLa/turnover_incl_cell_doubling |  |
| 0.45 | zecha2018:HeLa |  |
| 0.47 | li2021:HCT116/non-linear model |  |
| 0.55 | li2021:U2OS/manual inspection |  |
| 1.31 | schwanhausser2011:NIH3T3 |  |
| 12.77 | swovick2021:primary_fibroblast_hamster |  |
| 33.87 | li2025:kidney |  |

## EZH2 / EZH2  —  chosen **57.47 h**  (zecha2018:HeLa)
*8 arms, 11.6–57.5 h (5× spread). Model: NOT used — EZH2 protein t½ model-fit (~25 h effective)*

| t½ (h) | dataset : cell-type / variant | |
|---:|---|:--:|
| 4.98 | mcshane2016:NIH3T3/1state |  |
| 11.60 | mcshane2016:NIH3T3/2state |  |
| 18.99 | zecha2018:HeLa/turnover_incl_cell_doubling |  |
| 21.19 | nakanoh2026:ESC_derived_NP_mouse |  |
| 22.02 | schwanhausser2011:NIH3T3 |  |
| 26.11 | li2021:RPE1/linear model |  |
| 39.18 | nakanoh2026:ESC_derived_NP_human |  |
| 57.47 | zecha2018:HeLa | ✓ |

