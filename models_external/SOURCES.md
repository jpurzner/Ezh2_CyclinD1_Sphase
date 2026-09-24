# External model provenance

Downloaded 2026-09-24 via `.claude/skills/biomodels/scripts/bmfetch.py get <ID> -o models_external/biomodels`.
Each sha256 was checked against the manifest the BioModels API returns for that entry
(`bmfetch files <ID>`), so these are byte-verified, not merely "the download succeeded".

Load check: `tellurium.loadSBMLModel(f)` followed by `simulate(0, 100, 201)`, asserting all
outputs finite and non-negative.

| file | BioModels ID | what it is | sha256 (first 16) | loads + integrates |
|---|---|---|---|---|
| `BIOMD0000000154_url.xml` | BIOMD0000000154 | Geva-Zatorsky2006 p53 Model 3 | `582f640acfb6a4da…` | **no** — uses `delay()`; libRoadRunner cannot integrate DDEs |
| `BIOMD0000000155_url.xml` | BIOMD0000000155 | Geva-Zatorsky2006 p53 Model 6 | `e376475269e30be4…` | **no** — uses `delay()`; libRoadRunner cannot integrate DDEs |
| `BIOMD0000000156_url.xml` | BIOMD0000000156 | Geva-Zatorsky2006 p53 Model 5 | `70ca5dd300e73162…` | yes |
| `BIOMD0000000157_url.xml` | BIOMD0000000157 | Geva-Zatorsky2006 p53 Model 4 | `3346c709a413e31f…` | yes |
| `BIOMD0000000158_url.xml` | BIOMD0000000158 | Geva-Zatorsky2006 p53 Model 2 | `0c6c6337632df87e…` | yes |
| `BIOMD0000000159_url.xml` | BIOMD0000000159 | Geva-Zatorsky2006 p53 Model 1 | `18917194cccf661e…` | yes |
| `BIOMD0000000646_url.xml` | BIOMD0000000646 | Barr2016 - all-or-nothing G1/S transition | `37b8ce92f47b495f…` | yes |
| `BIOMD0000000660_url.xml` | BIOMD0000000660 | Barr2017 - p21 dynamics in hTert-RPE1 | `401b3adfc7d04058…` | yes |
| `BIOMD0000000700_url.xml` | BIOMD0000000700 | Heldt2018 - proliferation-quiescence decision (the v44 core) | `56ada2bc9084381c…` | yes |
| `Pandey2018_wildtype.xml` | BIOMD0000000954 | Pandey2018 - reversible quiescence<->proliferation | `a2556db19d998a0e…` | yes |
| `Yan2012.xml` | BIOMD0000000720 | Yan2012 - Rb-E2F dynamics with miR449 | `10f32e364304522c…` | yes |

## Notes

* **`MODEL1703030000` and `BIOMD0000000700` are the same model.** The brief listed the
  `MODEL…` submission ID; the API 301-redirects it to the curated `BIOMD…` accession. The
  downloaded file is byte-identical (md5 `bd4e342d43664c8e0b13b12d2789e90e`) to the
  `models_external/heldt2018_BIOMD700.xml` already in this repo, so v44's core is confirmed
  to be the published Heldt model, unmodified.
* **Geva-Zatorsky Models 3 and 6 (BIOMD…154, …155) cannot be run here.** They are delay
  differential equation models and libRoadRunner raises
  `Unable to support delay differential equations`. Models 1, 2, 4 and 5 integrate fine.
  If a DDE variant is ever needed, it has to be re-implemented (e.g. with a linear chain of
  intermediate species approximating the delay), not loaded.
* **Pandey2018 (BIOMD0000000954) and Yan2012 (BIOMD0000000720) were not in the brief.** They
  were found by the UniProt-annotation probe (see `docs/protein_noise_plan.md`) and added
  because both model the quiescence<->proliferation transition that the noise work is meant
  to drive. Pandey2018's file is named after the model, not the accession — that is how
  BioModels serves it.
* **EARM** (Albeck 2008 / Lopez 2013) is cloned at `models_external/repos/earm`
  (github.com/sorgerlab/earm, shallow clone). It is the reference for sampling measured
  protein distributions into an ODE ensemble; it is PySB, not SBML, and is read, not run.
* **HARISSA** — resolved: `harissa` 3.1.2 on PyPI, github.com/ulysseherbach/harissa,
  "Simulation and inference of gene regulatory networks based on transcriptional bursting".
  Confirms the brief's guess at the author. Not installed; it is a design reference for the
  PDMP burst module.
