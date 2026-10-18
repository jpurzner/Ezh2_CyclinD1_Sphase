# Results — model section (restructured draft, 2026-08-09)

> Restructured from the original single paragraph into **Skeleton / Parameterization / Population modeling**,
> per JP. Results-level altitude only — a full specification goes in Methods. Placeholders in [brackets].

---

Having established that Ezh2 directly represses Cyclin D1 through H3K27me3 (Fig. 4L), we sought to understand
the quantitative consequences of this repression for cell-cycle control and drug response. We therefore
constructed an ordinary differential equation (ODE) model that couples Hedgehog signaling and its target
genes (*Mycn*, *Ezh2*) to the mammalian cell-cycle engine (Fig. S8A). The model was assembled and used in
three stages — a cell-cycle **skeleton** built from established models, **parameterization** first to
generate cycling cells and then against our own data, and a **population** layer — which we describe in turn.
A complete specification of the equations, parameters, and fitting procedure is provided in the Methods; the
account below is intended to convey how the model was constructed.

## Model skeleton

The cell-cycle core is the real-time mammalian G1/S model of Heldt et al. (2018), in which DNA replication is
represented explicitly so that S-phase has a mechanistic, variable duration — a feature that proves essential
for modeling replicative dilution of the histone mark (below). To this core we added four modules, each
adapted from an established model of a specific cell-cycle behavior (Table 1): a Cyclin B/CDK1 mitotic switch
that makes entry into and exit from mitosis abrupt and all-or-none through Cdc25–Wee1 feedback (Pomerening et
al., 2003); an intra-S checkpoint by which ongoing DNA replication restrains mitotic entry, enforcing the
dependence of mitosis on completed replication (Hartwell and Weinert, 1989); a growth-gated restriction point
that requires the cell to reach a critical size before licensing commitment and S-phase entry (Ginzberg et
al., 2015; Cadart et al., 2018); and a Skp2–p27 feedforward switch, in which E2f-induced Skp2 degrades p27,
that renders the G0→G1 commitment a sharp, all-or-none decision (Carrano et al., 1999; Yao et al., 2008).
Because these components are drawn from the literature, the skeleton reproduces canonical cell-cycle behavior
independently of the tumor-specific regulation described below.

Onto this skeleton we grafted the Hedgehog–Ezh2–Cyclin D1 axis. Cyclin D1 mRNA transcription serves as the
central integration node, receiving three inputs — basal transcription, Gli-driven transcription (via
Shh/Ptch1/Smo/Gli), and Mycn-driven transcription — all of which are gated by an explicit H3K27me3 mark at the
Cyclin D1 locus. We modeled the mark as a fractional promoter occupancy that Ezh2 deposits and that is removed
both by demethylation and by replicative dilution — that is, halved at each round of DNA replication. Because
the mark throttles rather than silences transcription, a residual level of Cyclin D1 expression persists even
at maximal occupancy, so the mark acts as a graded rheostat rather than an on/off switch. Its steady level —
and hence the degree of Cyclin D1 repression — reflects a balance between deposition, which tracks Ezh2
(itself coupled to both mitogen dose and cell-cycle activity), and dilution, which tracks how rapidly the
cell divides. Capturing both arms of this balance is precisely why the core must represent DNA replication in
real time: Ezh2 is a stable protein that accumulates in proportion to S-phase duration, whereas the mark is
diluted in proportion to the frequency of S-phase.

**Table 1. Cell-cycle skeleton: modules and their literature sources.**

| Module | Behavior it provides | Source |
|---|---|---|
| Real-time G1/S core | explicit DNA replication → mechanistic, variable S-phase duration | Heldt et al., 2018 |
| Cyclin B/CDK1 mitotic switch | abrupt, all-or-none mitotic entry and exit (Cdc25–Wee1 bistability) | Pomerening et al., 2003 |
| Intra-S / replication checkpoint | ongoing replication restrains mitotic entry | Hartwell and Weinert, 1989 |
| Growth-gated restriction point | critical-size requirement licenses commitment and S-phase entry | Ginzberg et al., 2015; Cadart et al., 2018 |
| Skp2–p27 feedforward switch | sharp, all-or-none G0→G1 commitment (E2f → Skp2 → p27 degradation) | Carrano et al., 1999; Yao et al., 2008 |

*The Hedgehog→Gli/Mycn inputs and the H3K27me3 rheostat on Cyclin D1 are specific to this study (Fig. 4);
the replicative-dilution behavior of H3K27me3 is supported by [H3K27me3 dilution/restoration ref].*

## Parameterization

We parameterized the model in two stages. First, we tuned the cell-cycle skeleton to generate normally
cycling cells: parameters were set so that the model produces sustained oscillations of the core regulators,
a physiological cell-cycle length, and the expected succession of phases, confirming that the engine runs
before any tumor-specific tuning. Second, we constrained the Hedgehog–Ezh2–Cyclin D1 layer against our own
measurements, which we organized into a structured panel of quantitative targets (Table 2): the relative
transcript and protein levels of Cyclin D1, Ezh2, Gli1 and Mycn in granule-neuron precursors (GNPs) versus
medulloblastoma (MB); the de-repression of Cyclin D1 upon Ezh2 inhibition; the proliferative response to
Hedgehog blockade (vismodegib) and to CDK4/6 inhibition; the accumulation of cells under replication stress
(hydroxyurea); and the cell-cycle period and phase fractions. We required a single parameter set to reproduce
this entire panel simultaneously ([N of M targets]; Fig. S8[B–M]), so that the model accounts for all of the
data jointly rather than any one observation in isolation. Parameter values, their sources and bounds, and
the fitting procedure are detailed in the Methods.

**Table 2. Structured data used to parameterize the regulatory layer.**

| Data class | Measurement | Comparison / assay |
|---|---|---|
| Transcript level | Cyclin D1, Gli1 | MB vs GNP (RNA-seq) |
| Protein level | Ezh2, CDK6 | MB vs GNP (TMT proteomics) |
| Epigenetic de-repression | Cyclin D1 after Ezh2 inhibition | fold change |
| Mitogen dependence | Cyclin D1, proliferation after Hedgehog blockade | vismodegib dose |
| Downstream dependence | proliferation after CDK4/6 inhibition | palbociclib |
| Replication stress | phase accumulation | hydroxyurea |
| Cell-cycle timing | cycle period | GNP |
| Phase distribution | G0/G1/S/G2 fractions | flow cytometry (GNP, MB) |

## Population modeling

The equations above describe a single, average cell, but several of the phenomena we set out to explain are
properties of a population: only a fraction of MB cells reside in a transient G0 at any given time, and
Hedgehog or Ezh2 inhibition arrests some cells while others continue to divide. To represent this, we
simulated a population of cells whose key regulatory setpoints — for example, the levels of the CDK inhibitors
and of Cyclin D1 — were drawn from distributions rather than held fixed. Simulating this ensemble reproduces
the graded, single-cell-resolved responses measured by flow cytometry — the fraction of cycling versus
arrested cells, and how that fraction shifts with each perturbation — which a single deterministic cell
cannot capture. We used this population layer for the drug-response analyses in Fig. [X].

---

### Notes / placeholders for JP
- **[N of M targets]** — insert the final validation count for the version cited in the paper (current baked
  model: 31/32; earlier calibrations reported 28/28 and 28/29 — pick the one that matches the figures).
- **[Fig. X], [Fig. S8B–M]** — figure numbers.
- **H3K27me3 replicative-dilution reference** — a citation for the two-fold dilution of histone marks at
  replication (I did not invent one; e.g. the Groth-lab histone-recycling / mark-restoration literature).
- Table 2 lists data *classes*; if you prefer, I can expand each row to the specific numeric target and
  tolerance (that content is essentially the validation-target table and may fit better in Methods/Supp).
