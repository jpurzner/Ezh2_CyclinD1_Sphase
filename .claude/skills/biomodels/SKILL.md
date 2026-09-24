---
name: biomodels
description: Find, download, verify and load models from the EBI BioModels repository (SBML), and check which published models cover a given protein. Use when asked to fetch a BIOMD/MODEL accession, look for a published ODE model of a pathway or protein, resolve a paper to its deposited model, or check that an external SBML actually integrates in Tellurium/libRoadRunner.
---

# BioModels

The EBI BioModels REST API, as it actually behaves. Every endpoint and quirk below was probed
live against `https://www.ebi.ac.uk/biomodels` — where the documented behaviour and the real
behaviour differ, the real behaviour is what is written here.

Use `scripts/bmfetch.py` rather than hand-rolling requests; it handles the redirects, the
retries and the checksum verification that this API needs.

```bash
python3 .claude/skills/biomodels/scripts/bmfetch.py search "cell cycle" -n 10 --facets
python3 .claude/skills/biomodels/scripts/bmfetch.py info   BIOMD0000000700
python3 .claude/skills/biomodels/scripts/bmfetch.py files  BIOMD0000000700
python3 .claude/skills/biomodels/scripts/bmfetch.py get    BIOMD0000000700 -o models_external/biomodels
python3 .claude/skills/biomodels/scripts/bmfetch.py probe  P24385,Q15910,P04198
```

## The five things that will bite you

**1. Everything 301-redirects.** `curl` without `-L` returns an empty body and exit code 0, so a
naive download looks like a silent failure. Use `curl -sL`, or `urllib`, which follows redirects
by default.

**2. `MODEL…` and `BIOMD…` accessions can be the same model.** A submission gets a `MODEL…` ID; if
it is later manually curated it also gets a `BIOMD…` accession, and the API redirects the former
to the latter. `MODEL1703030000` and `BIOMD0000000700` are one model (Heldt 2018). Always resolve
an accession with `info` before concluding you have two different models.

**3. `format=json` is not optional and the search endpoint is slow.** Plain keyword search can
take 30–60 s and intermittently times out. `bmfetch` retries with backoff; anything you write
yourself should too.

**4. Query syntax: `field:value` works, quoted multi-word filters return HTTP 400.**

| query | result |
|---|---|
| `Heldt` | works |
| `publication:Heldt` | works |
| `UNIPROT:P24385` | works — this is the useful one |
| `GO:0000278` | works |
| `cell cycle AND curationstatus:"Manually curated"` | **HTTP 400** |
| `cell cycle` + separate `filter=` param | **HTTP 400** |

To filter by curation status or modelling approach, **read the facets off an unfiltered response**
instead (`--facets`). Every search returns `facets` for `curationstatus`, `modelformat`,
`modellingapproach`, `TAXONOMY`, `GO`, `UNIPROT` and `CHEBI`, with counts.

**5. Downloading without `filename=` gives you an OMEX zip, not the SBML.**
`/model/download/<ID>` alone returns a COMBINE archive (`.omex`). To get the SBML, pass the main
file's name from the manifest — which is what `bmfetch get` does automatically.

## Endpoints

| purpose | endpoint |
|---|---|
| search | `/search?query=<q>&numResults=<n>&format=json` |
| metadata | `/<ID>?format=json` |
| file manifest | `/model/files/<ID>?format=json` |
| download one file | `/model/download/<ID>?filename=<name>` |
| download OMEX archive | `/model/download/<ID>` |

The manifest splits files into `main` (the SBML you want) and `additional` (auto-generated
BioPAX/MATLAB/Octave/XPP, the curation figure, SED-ML, COPASI `.cps`, RDF metadata). Each entry
carries `fileSize`, `md5sum`, `sha1sum` and `sha256sum` — **verify the download against the
manifest sha256**; `bmfetch get` does this and prints `sha256 OK` or `!! SHA256 MISMATCH`.

## Finding a model for a protein

`UNIPROT:<accession>` search is the high-value query: it hits the curated annotations rather than
free text, so it finds models that never mention the gene symbol. Use `bmfetch probe` for a list.

Two caveats:

* **Annotations are predominantly human accessions**, even for models built on mouse data. Probe
  with the human accession and note the model's actual organism from `info`.
* **A zero result is real information.** It means no curated BioModels entry annotates that
  protein — so there is no published model to borrow structure or parameters from, and any
  parameter for it has to come from primary literature.

Get accessions from the UniProt REST API, and **read the gene symbol off the returned entry**:

```bash
curl -s 'https://rest.uniprot.org/uniprotkb/search?query=gene:SMO+AND+organism_id:9606+AND+reviewed:true&fields=accession,gene_primary,protein_name&format=json&size=5'
```

`gene_exact:` **is not exact** — `gene_exact:SMO` ranks `SMOX` (spermine oxidase) above `SMO`
(smoothened). Filter the results on `genes[0].geneName.value` yourself rather than taking the
first hit. Three accessions in this repo's own species map were initially wrong for exactly this
reason and were caught only by reading the returned symbols back.

## Loading and checking an SBML

Downloading is not validating. A deposited model can fail to integrate.

```python
import tellurium as te, numpy as np
r = te.loadSBMLModel(path)
res = r.simulate(0, 100, 201)
assert np.all(np.isfinite(res))
```

Known failure mode: **delay differential equations**. Models using SBML's `delay(x, tau)` raise
`RuntimeError: Unable to support delay differential equations` from libRoadRunner. Two of the six
Geva-Zatorsky p53 models (BIOMD0000000154, BIOMD0000000155) are DDEs and cannot be run here — they
would need re-implementing with a linear chain of intermediates approximating the delay.

Also check, before trusting a converted model:

* `hasOnlySubstanceUnits` on each species, and the compartment volumes — amounts and
  concentrations are easy to mix up when converting a model to molecule counts for stochastic
  simulation.
* Whether a quantity that should be a **parameter with a rate rule** has been encoded as a
  species in a compartment. Making it a compartment species silently rescales it by the volume.

## Recording provenance

Anything downloaded goes under `models_external/` with a `SOURCES.md` recording, per file: the
accession, the URL, the date, the sha256, and whether it loads and integrates. See
`models_external/SOURCES.md` in this repo for the format.
