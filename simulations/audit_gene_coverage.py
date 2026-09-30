#!/usr/bin/env python3
"""Per-gene data coverage audit for the v44 model species.

Answers "what do we actually have for each gene in the model?" across every dataset held,
rather than per-dataset. Reads data/priors/turnover_burst_priors.csv plus Battich's
cell-cycle table directly (Battich is relative-rate only, so it is not in the priors CSV).

Usage:  python3 simulations/audit_gene_coverage.py [--csv out.csv]
"""
import argparse
import csv
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRIORS = os.path.join(ROOT, "data", "priors", "turnover_burst_priors.csv")
BATTICH = os.path.join(ROOT, "data", "priors", "raw", "battich2020", "aax3072_table-s1.csv")

# What each column needs, as (label, quantity, source-prefix, cell_type-substring)
COLS = [
    ("prot#",    "protein_copy_number",     "Schwanhausser", "NIH3T3"),
    ("mRNA#",    "mrna_copy_number",        "Schwanhausser", "NIH3T3"),
    ("t_p Sch",  "protein_half_life",       "Schwanhausser", "NIH3T3"),
    ("t_p McS",  "protein_half_life_1state","McShane",       "NIH3T3"),
    ("t_p neur", "protein_half_life",       "Mathieson",     "neuron"),
    ("t_m Sch",  "mrna_half_life",          "Schwanhausser", "NIH3T3"),
    ("t_m Her",  "mrna_half_life",          "Herzog",        "ESC"),
    ("B fib",    "burst_size_B",            "Larsson",       "fibroblast"),
    ("B esc",    "burst_size_B",            "Larsson",       "embryonic stem"),
    ("ksp",      "translation_rate_ksp",    "Schwanhausser", "NIH3T3"),
]


def load_priors():
    d = {}
    for r in csv.DictReader(open(PRIORS)):
        d.setdefault(r["gene"], []).append(r)
    return d


def pick(rows, quantity, src, ct):
    for r in rows:
        if r["quantity"] == quantity and r["source"].startswith(src) and ct in r["cell_type"]:
            return r["value"]
    return "NA"


def load_battich():
    """gene -> strategy group. Battich gives cycle-RELATIVE rates only (log2 vs cycle median)."""
    try:
        rows = list(csv.reader(open(BATTICH)))
    except FileNotFoundError:
        return {}
    hdr = rows[1]
    gi, si = hdr.index("gene_symbol"), hdr.index("strategy_group")
    return {r[gi].strip().upper(): r[si] for r in rows[2:] if len(r) > max(gi, si) and r[gi]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", help="also write the matrix here")
    a = ap.parse_args()

    P, B = load_priors(), load_battich()
    genes = list(P.keys())

    hdr = f"{'gene':9s} " + " ".join(f"{c[0]:>8s}" for c in COLS) + f" {'Batt':>5s}  role"
    print("### Per-gene coverage across all held datasets")
    print("### '.' = present, '-' = absent from that dataset\n")
    print(hdr)
    print("-" * len(hdr))

    out = []
    for g in genes:
        rows = P[g]
        vals = [pick(rows, q, s, c) for _, q, s, c in COLS]
        marks = "".join(f"{('.' if v != 'NA' else '-'):>9s}" for v in vals)
        bat = B.get(g.upper(), "")
        role = next((r["notes"] for r in rows if r["quantity"] == "model_role"), "")
        n = sum(1 for v in vals if v != "NA")
        print(f"{g:9s}{marks} {(bat or '-'):>5s}  {role[:42]}")
        out.append(dict(gene=g, n_filled=n, battich_strategy=bat, role=role,
                        **{c[0]: v for c, v in zip(COLS, vals)}))

    # summary
    print()
    tot = len(genes)
    for i, (lab, q, s, c) in enumerate(COLS):
        n = sum(1 for g in genes if pick(P[g], q, s, c) != "NA")
        print(f"  {lab:10s} {n:2d}/{tot}  ({100*n/tot:3.0f}%)")
    nb = sum(1 for g in genes if g.upper() in B)
    print(f"  {'Battich':10s} {nb:2d}/{tot}  ({100*nb/tot:3.0f}%)")

    complete = [d["gene"] for d in out if d["n_filled"] >= 8]
    empty = [d["gene"] for d in out if d["n_filled"] <= 2]
    print(f"\n  best covered (>=8/10): {', '.join(complete) or 'none'}")
    print(f"  worst covered (<=2/10): {', '.join(empty) or 'none'}")

    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
            w.writeheader(); w.writerows(out)
        print(f"\nwrote {a.csv}")


if __name__ == "__main__":
    main()
