#!/usr/bin/env python3
"""T3 -- build the turnover / burst prior table for the v44 model species.

Joins three published datasets onto the genes encoding v44's species and writes
`data/priors/turnover_burst_priors.csv` (long/tidy: one row per gene x quantity x source).

Sources (raw files under data/priors/raw/, provenance in data/priors/SOURCES.md):
  Schwanhausser 2011 Nature 473:337, Suppl. Table 3  -- mouse NIH3T3: protein & mRNA copy
      number, protein & mRNA half-life, transcription rate (vsr), translation rate (ksp).
      NB: the hosted file is the 2013 CORRIGENDUM version (nature11848).
  Larsson 2019 Nature 565:251, Suppl. Tables 1 & 3   -- mouse primary fibroblasts / ES cells:
      allele-resolved transcriptional burst frequency and burst size (telegraph model).
  Mathieson 2018 Nat Commun 9:689, Suppl. Data 2     -- primary cells incl. MOUSE NEURONS:
      protein half-life. Neurons are POST-MITOTIC, so these are dilution-free (pure
      degradation) -- the class-(b) reference for the T4 dilution audit.
  McShane 2016 Cell 167:803, mmc1 (NIH3T3) and mmc4 (RPE) -- pulse-chase time course fitted
      with BOTH a 1-state exponential and a 2-state model (fast pool A -> stable pool B).
      Gives the degradation PROFILE (ED / NED / UN) and, for the 2-state fit, kA, kB and
      ktransfer. RPE is the line Heldt 2018 was fitted in, so mmc4 is the direct reference
      for v44's Heldt-inherited constants.
  Eden 2011 Science 331:764, SOM Table S4              -- bleach-chase in LIVING human H1299,
      100 proteins. NB bleach-chase measures the TOTAL removal rate alpha = degradation +
      dilution, so these are class (a) EFFECTIVE values, unlike Schwanhausser/Mathieson.

UNITS / CAVEATS that the CSV carries per row:
  * Larsson burst frequency is reported PER mRNA-DEGRADATION-RATE (delta), not per hour.
    To get per-hour: f[1/h] = f[/delta] * ln2 / t_half_mRNA[h]. Burst size is dimensionless
    (mRNAs per burst). Values are per ALLELE; a diploid cell has 2x the frequency.
  * Schwanhausser half-lives are the well-known long-tail estimates; Li, Bickel & Biggin 2014
    (PeerJ 2:e270) argue copy numbers and translation rates are several-fold higher. Treat
    single-gene half-lives as order-of-magnitude, and sweep the ones the model is sensitive to.
  * growth_corrected: 'yes' for Mathieson (post-mitotic primary cells -> no dilution term),
    'unknown' for Schwanhausser (cycling NIH3T3; the Methods describe a correction but this
    has NOT been re-verified here -- see SOURCES.md).

Usage:  python3 simulations/build_turnover_priors.py
"""
import ast
import csv
import os
import re
import statistics
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "priors", "raw")
OUT = os.path.join(ROOT, "data", "priors", "turnover_burst_priors.csv")

# ---------------------------------------------------------------- species map
# v44 species -> the gene that encodes it. Accessions verified against the UniProt REST API
# (rest.uniprot.org, reviewed entries only). `role` flags what the species does in the model,
# which is what decides whether it needs explicit bursts (option 4) or drift (option 3).
#
# NB `gene_exact:` in the UniProt query language is NOT exact -- gene_exact:SMO ranks SMOX
# first. Every accession below was read off the returned entry, not the query.
SPECIES = [
    # v44 species,          gene,      mouse,    human,    role
    ("Rb/Rbm/pRb/RbE2f",    "Rb1",     "P13405", "P06400", "R-point; size-diluted (Zatulovskiy)"),
    ("E2f/RbE2f/RbmE2f",    "E2f1",    "Q61501", "Q01094", "R-point TF; bistable switch"),
    ("E1/E1C1",             "Fbxo5",   "Q7TSG3", "Q9UKT4", "Emi1, APC/C inhibitor"),
    ("C1/pC1/E1C1",         "Fzr1",    "Q9R1K5", "Q9UM11", "Cdh1, APC/C co-activator"),
    ("P21/CeP21/CaP21",     "Cdkn1a",  "P39689", "P38936", "CKI; BURST CANDIDATE (Barr 2017)"),
    ("Ce/CeP21",            "Ccne1",   "Q61457", "P24864", "G1/S cyclin"),
    ("Ca/CaP21",            "Ccna2",   "P51943", "P20248", "S/G2 cyclin"),
    ("aPcna/iPcna",         "Pcna",    "P17918", "P12004", "replication clamp"),
    ("P53",                 "Trp53",   "P02340", "P04637", "damage response"),
    ("Cd/Cd_mRNA",          "Ccnd1",   "P25322", "P24385", "MITOGEN INTEGRATOR; BURST CANDIDATE"),
    ("Skp2",                "Skp2",    "Q9Z0Z3", "Q13309", "bistable commitment (Skp2-p27)"),
    ("MPF/preMPF",          "Ccnb1",   "P24860", "P14635", "mitotic cyclin"),
    ("MPF/preMPF",          "Cdk1",    "P11440", "P06493", "mitotic kinase"),
    ("Cdc20",               "Cdc20",   "Q9JJ66", "Q12834", "APC/C co-activator, mitotic exit"),
    ("EZH2/EZH2m",          "Ezh2",    "Q61188", "Q15910", "PRC2 catalytic; THE EPIGENETIC BRAKE"),
    ("(PRC2, implicit)",    "Suz12",   "Q80U70", "Q15022", "PRC2 core (not an explicit species)"),
    ("(PRC2, implicit)",    "Eed",     "Q921E6", "O75530", "PRC2 core (not an explicit species)"),
    ("Gli1/Gli1_mRNA",      "Gli1",    "P47806", "P08151", "HH output; MB marker"),
    ("Gli_act/Gli_rep",     "Gli2",    "Q0VGT2", "P10070", "lumped into Gli_act"),
    ("Gli_act/Gli_rep",     "Gli3",    "Q61602", "P10071", "lumped into Gli_rep"),
    ("MYCN",                "Mycn",    "P03966", "P04198", "MB amplification; mitogen floor"),
    ("Ptch1_free/_mRNA",    "Ptch1",   "Q61115", "Q13635", "HH receptor; GNP=1 / MB=0.1 dosage"),
    ("Smo_active",          "Smo",     "P56726", "Q99835", "vismodegib target"),
    ("(Cdt2 param)",        "Dtl",     "Q3TLR7", "Q9NZJ0", "CRL4-Cdt2, p21 degron (param, not species)"),
    ("(Wee1a assignment)",  "Wee1",    "P47810", "P30291", "mitotic switch (algebraic, not a species)"),
    ("(Cdc25a assignment)", "Cdc25a",  "P48964", "P30304", "mitotic switch (algebraic, not a species)"),
    ("(kPhRbCd param)",     "Cdk4",    "P30285", "P11802", "CyclinD partner (lumped into Cd)"),
    ("(kPhRbCd param)",     "Cdk6",    "Q64261", "Q00534", "CyclinD partner (lumped into Cd)"),
    ("(Ce/Ca partner)",     "Cdk2",    "P97377", "P24941", "lumped into Ce/Ca"),
    ("(P21 pool)",          "Cdkn1b",  "P46414", "P46527", "p27; v44 lumps p21/p27 into P21"),
    ("(p18 param)",         "Cdkn2c",  "Q60772", "P42773", "INK4c; constitutive GNP brake"),
    ("(p19 param)",         "Cdkn2d",  "Q60773", "P55273", "INK4d"),
    ("(Chk1 assignment)",   "Chek1",   "O35280", "O14757", "intra-S checkpoint (algebraic)"),
]


def _f(v):
    try:
        x = float(v)
        return None if x != x else x
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------- McShane
def load_mcshane(fname, sheet):
    """gene -> (1-state half-life h, profile, kA, kB, ktransfer). Profile is ED/NED/UN."""
    import openpyxl
    wb = openpyxl.load_workbook(os.path.join(RAW, fname), read_only=True)
    ws = wb[sheet]
    it = ws.iter_rows(values_only=True)
    hdr = list(next(it))
    out = {}
    for row in it:
        if not row or not row[1]:
            continue
        d = dict(zip(hdr, row))
        rec = (_f(d.get("Half-life (exponential) 1-state-model [h]")),
               str(d.get("Degradation profile") or "").strip(),
               _f(d.get("2-state-model_kA [h-1]")),
               _f(d.get("2-state-model_kB [h-1]")),
               _f(d.get("2-state-model_ktransfer [h-1]")))
        for g in str(row[1]).split(";"):
            g = g.strip()
            if g and g not in out:
                out[g] = rec
    return out


# ------------------------------------------------------------------- Eden
def load_eden():
    """gene -> bleach-chase half-life (h). Parsed from the SOM PDF, Table S4."""
    import csv as _csv
    path = os.path.join(ROOT, "data", "priors", "eden2011_bleachchase_halflives.csv")
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return {r["gene"].strip().upper(): float(r["half_life_h"]) for r in _csv.DictReader(fh)}


# ------------------------------------------------------------ Schwanhausser
def load_schwanhausser():
    import xlrd
    path = os.path.join(RAW, "schwanhausser2011_TableS3.xls")
    sheet = xlrd.open_workbook(path).sheet_by_index(0)
    # col -> (quantity, units)
    COLS = {13: ("protein_copy_number", "molecules/cell"),
            16: ("mrna_copy_number", "molecules/cell"),
            19: ("protein_half_life", "h"),
            22: ("mrna_half_life", "h"),
            25: ("transcription_rate_vsr", "molecules/(cell*h)"),
            28: ("translation_rate_ksp", "molecules/(mRNA*h)")}
    out = {}
    for r in range(1, sheet.nrows):
        row = sheet.row_values(r)
        for g in str(row[2]).split(";"):
            g = g.strip()
            if g and g not in out:
                out[g] = {q: (_f(row[i]), u) for i, (q, u) in COLS.items()}
    return out


# ------------------------------------------------------------------- Larsson
_NUM = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def _first(cell):
    """Larsson stores numpy-repr arrays as text, e.g. '[0.965 0.745 1.249]'."""
    if cell is None:
        return None
    m = _NUM.findall(str(cell))
    return float(m[0]) if m else None


def load_larsson():
    import openpyxl
    out = {}
    for fname, tissue in (("larsson2019_TableS1_fibroblast_burst.xlsx", "mouse primary fibroblast"),
                          ("larsson2019_TableS3_ESC_burst.xlsx", "mouse embryonic stem cell")):
        wb = openpyxl.load_workbook(os.path.join(RAW, fname), read_only=True)
        per_gene = {}
        for allele in ("C57", "CAST"):
            if allele not in wb.sheetnames:
                continue
            for row in wb[allele].iter_rows(min_row=2, values_only=True):
                if not row or not row[0]:
                    continue
                f, b, mu = _first(row[2]), _first(row[3]), _f(row[4])
                if f is None or b is None:
                    continue
                per_gene.setdefault(str(row[0]), []).append((f, b, mu))
        # average the two alleles
        out[tissue] = {g: (statistics.mean(v[0] for v in vs),
                           statistics.mean(v[1] for v in vs),
                           statistics.mean(v[2] for v in vs if v[2] is not None) if any(v[2] is not None for v in vs) else None,
                           len(vs))
                       for g, vs in per_gene.items()}
    return out


# ----------------------------------------------------------------- Mathieson
def load_mathieson():
    import openpyxl
    wb = openpyxl.load_workbook(os.path.join(RAW, "mathieson2018_SupplData2.xlsx"), read_only=True)
    ws = wb["protein half lives high qual"]
    hdr = [str(h) if h else "" for h in next(ws.iter_rows(max_row=1, values_only=True))]
    # keep the two mouse-neuron replicate half-life columns (post-mitotic -> dilution-free)
    cols = [i for i, h in enumerate(hdr) if h.startswith("Mouse Neurons") and h.endswith("half_life")]
    out = {}
    for row in ws.iter_rows(min_row=2, values_only=True):
        if not row or not row[0]:
            continue
        vals = [_f(row[i]) for i in cols]
        vals = [v for v in vals if v is not None]
        if vals:
            out[str(row[0]).strip().upper()] = (statistics.mean(vals), len(vals))
    return out


# ---------------------------------------------------------------------- main
def main():
    try:
        sch, lar, mat = load_schwanhausser(), load_larsson(), load_mathieson()
        mcs_3t3 = load_mcshane("mcshane2016_mmc1_3T3.xlsx", "FinalTable_PostNorm_3T3_")
        mcs_rpe = load_mcshane("mcshane2016_mmc4_RPE.xlsx", "FinalTable_RPE_")
        eden = load_eden()
    except ImportError as e:
        sys.exit(f"missing reader dependency ({e}); need xlrd (for .xls) and openpyxl (.xlsx)")

    rows = []

    def add(sp, gene, quantity, value, units, source, table, celltype, growth_corrected, notes=""):
        rows.append(dict(v44_species=sp, gene=gene, quantity=quantity,
                         value=("NA" if value is None else f"{value:.6g}"), units=units,
                         source=source, table=table, cell_type=celltype,
                         growth_corrected=growth_corrected, notes=notes))

    for sp, gene, mouse, human, role in SPECIES:
        add(sp, gene, "uniprot_mouse", None, "accession", "UniProt", "REST api", "-", "n/a",
            notes=mouse)
        add(sp, gene, "uniprot_human", None, "accession", "UniProt", "REST api", "-", "n/a",
            notes=human)
        add(sp, gene, "model_role", None, "-", "this repo", "src/build_model_v44_heldt.py",
            "-", "n/a", notes=role)

        s = sch.get(gene)
        for q, units in (("protein_copy_number", "molecules/cell"),
                         ("mrna_copy_number", "molecules/cell"),
                         ("protein_half_life", "h"), ("mrna_half_life", "h"),
                         ("transcription_rate_vsr", "molecules/(cell*h)"),
                         ("translation_rate_ksp", "molecules/(mRNA*h)")):
            val = s[q][0] if s else None
            add(sp, gene, q, val, units, "Schwanhausser 2011 (2013 corrigendum)",
                "Suppl. Table 3", "mouse NIH3T3", "unknown",
                notes="" if s else "gene absent from the 5029-protein table")

        for tissue, tbl in lar.items():
            hit = tbl.get(gene)
            note = ("mean of C57+CAST alleles" if hit and hit[3] == 2 else
                    "single allele only" if hit else "gene absent from table")
            add(sp, gene, "burst_frequency", hit[0] if hit else None,
                "per mRNA-degradation-rate (delta), per allele",
                "Larsson 2019", "Suppl. Table 1 (fibro) / 3 (ESC)", tissue, "n/a",
                notes=note + "; f[1/h] = value*ln2/t_half_mRNA[h]")
            add(sp, gene, "burst_size_B", hit[1] if hit else None, "mRNAs per burst",
                "Larsson 2019", "Suppl. Table 1 (fibro) / 3 (ESC)", tissue, "n/a", notes=note)
            add(sp, gene, "mrna_mean_per_allele", hit[2] if hit else None, "molecules/allele",
                "Larsson 2019", "Suppl. Table 1 (fibro) / 3 (ESC)", tissue, "n/a", notes=note)

        m = mat.get(gene.upper())
        add(sp, gene, "protein_half_life", m[0] if m else None, "h",
            "Mathieson 2018", "Suppl. Data 2 (high qual)", "mouse primary neuron", "yes",
            notes=("mean of %d replicate(s); POST-MITOTIC -> pure degradation, no dilution" % m[1])
                  if m else "gene absent from the high-quality table")

        # McShane: 1-state half-life + the 2-state (fast/stable pool) decomposition
        for tbl, key, ct in ((mcs_3t3, gene, "mouse NIH3T3"), (mcs_rpe, gene.upper(), "human RPE")):
            hit = tbl.get(key) or tbl.get(gene)
            prof = hit[1] if hit else ""
            note = (f"degradation profile {prof}" if hit else "gene absent from table")
            if ct == "human RPE":
                note += "; RPE is the line Heldt 2018 was fitted in"
            add(sp, gene, "protein_half_life_1state", hit[0] if hit else None, "h",
                "McShane 2016", "mmc1 (3T3) / mmc4 (RPE)", ct, "unknown", notes=note)
            add(sp, gene, "degradation_profile", None, "ED|NED|UN",
                "McShane 2016", "mmc1 (3T3) / mmc4 (RPE)", ct, "n/a",
                notes=(prof or "NA") + " (NED = non-exponential: fast unassembled pool then stable)")
            add(sp, gene, "deg_rate_kA_fast_pool", hit[2] if hit else None, "1/h",
                "McShane 2016", "mmc1 (3T3) / mmc4 (RPE)", ct, "unknown",
                notes="2-state model: newly made / unassembled pool")
            add(sp, gene, "deg_rate_kB_stable_pool", hit[3] if hit else None, "1/h",
                "McShane 2016", "mmc1 (3T3) / mmc4 (RPE)", ct, "unknown",
                notes="2-state model: assembled / complexed pool -- match model species to THIS "
                      "when the species represents the functional complexed form")
            add(sp, gene, "deg_rate_ktransfer", hit[4] if hit else None, "1/h",
                "McShane 2016", "mmc1 (3T3) / mmc4 (RPE)", ct, "unknown",
                notes="2-state model: rate of A -> B (maturation / assembly)")

        e = eden.get(gene.upper())
        add(sp, gene, "protein_half_life_effective", e, "h",
            "Eden 2011", "SOM Table S4", "human H1299 (living)", "no",
            notes=("bleach-chase; INCLUDES dilution (alpha = degradation + dilution)"
                   if e else "gene absent from the 100-protein set"))

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    n_na = sum(1 for r in rows if r["value"] == "NA" and r["quantity"] not in
               ("uniprot_mouse", "uniprot_human", "model_role"))
    n_num = sum(1 for r in rows if r["quantity"] not in
                ("uniprot_mouse", "uniprot_human", "model_role"))
    print(f"wrote {OUT}")
    print(f"  {len(rows)} rows, {len(SPECIES)} gene entries")
    print(f"  numeric slots: {n_num - n_na}/{n_num} filled, {n_na} explicit NA (never imputed)")


if __name__ == "__main__":
    main()
