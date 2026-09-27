#!/usr/bin/env python3
"""Independent test of Ochiai 2020's claim that PRC2 lowers burst frequency and raises burst size.

Ochiai state the result in one sentence pointing at a correlation heatmap (Fig. 2E), with no
coefficients, no effect sizes and no statistics in the text, and their deposited tables are the
raw count matrices rather than the fitted kinetics. So the claim cannot be checked from their
supplement. This script rebuilds it from independent public data:

    ChIP-Atlas gene-level Suz12 / Ezh2 binding in ES_cells   (333 / 133 experiments, mm10)
      x Larsson 2019 Table S3 burst kinetics in mouse ESC     (same cell type)

THE CONFOUNDER this is built to expose: PRC2 marks lowly expressed genes, and BOTH burst
frequency and burst size correlate positively with mean expression (r ~ +0.5 / +0.3). So a naive
marginal correlation of any repressive mark against f and B will be negative for both, and a
"frequency down, size up" pattern can only be assessed after conditioning on expression.

Marginally f and B are uncorrelated (r ~ 0.00), but at fixed mean they trade off
(partial r ~ -0.20 to -0.30), so conditioning is not optional.

Requires: /tmp/ca_Suz12.tsv, /tmp/ca_Ezh2.tsv  (see SOURCES for the ChIP-Atlas URLs)
Usage:    python3 simulations/test_prc2_bursting.py
"""
import os
import re
import statistics
import sys

import numpy as np
import openpyxl

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CHIP = {"Suz12": "/tmp/ca_Suz12.tsv", "Ezh2": "/tmp/ca_Ezh2.tsv"}
BURST = os.path.join(ROOT, "data", "priors", "raw", "larsson2019_TableS3_ESC_burst.xlsx")
_NUM = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def load_chip(path, cell_type):
    """gene -> (fraction of experiments with a called peak, mean MACS2 score)."""
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        cols = [i for i, h in enumerate(hdr) if "|" in h and h.split("|", 1)[1] == cell_type]
        if not cols:
            raise SystemExit(f"no '{cell_type}' columns in {path}")
        out = {}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(hdr):
                continue
            v = []
            for i in cols:
                try:
                    v.append(float(p[i]))
                except ValueError:
                    v.append(0.0)
            a = np.array(v)
            out[p[0].strip()] = (float((a > 0).mean()), float(a.mean()))
    return out, len(cols)


def load_burst(path):
    """gene -> (burst frequency, burst size, mean expression), averaged over both alleles."""
    def first(c):
        m = _NUM.findall(str(c)) if c is not None else []
        return float(m[0]) if m else None

    wb = openpyxl.load_workbook(path, read_only=True)
    acc = {}
    for allele in ("C57", "CAST"):
        if allele not in wb.sheetnames:
            continue
        for row in wb[allele].iter_rows(min_row=2, values_only=True):
            if not row or not row[0]:
                continue
            f, b = first(row[2]), first(row[3])
            try:
                mu = float(row[4])
            except (TypeError, ValueError):
                mu = None
            if None in (f, b, mu) or min(f, b, mu) <= 0:
                continue
            acc.setdefault(str(row[0]), []).append((f, b, mu))
    return {g: tuple(statistics.mean(v[i] for v in vs) for i in range(3))
            for g, vs in acc.items()}


def resid(y, x):
    A = np.vstack([x, np.ones_like(x)]).T
    return y - A @ np.linalg.lstsq(A, y, rcond=None)[0]


def main():
    missing = [p for p in CHIP.values() if not os.path.exists(p)]
    if missing:
        sys.exit("missing ChIP-Atlas files: " + ", ".join(missing) +
                 "\n  curl -sL -o /tmp/ca_Suz12.tsv "
                 "https://chip-atlas.dbcls.jp/data/mm10/target/Suz12.10.tsv")
    L = load_burst(BURST)

    for antigen, path in CHIP.items():
        C, n_exp = load_chip(path, "ES_cells")
        sh = sorted(set(L) & set(C))
        fb = np.array([C[g][0] for g in sh])
        f = np.log([L[g][0] for g in sh])
        b = np.log([L[g][1] for g in sh])
        m = np.log([L[g][2] for g in sh])

        print(f"\n{'='*78}")
        print(f"### {antigen}  (ChIP-Atlas ES_cells, {n_exp} experiments)  x  Larsson ESC bursting")
        print(f"    n = {len(sh)} genes")
        print(f"\n  {'':28s} {'marginal r':>11s} {'partial r | mean expr':>23s}")
        for nm, y in (("burst FREQUENCY  log f", f),
                      ("burst SIZE       log B", b),
                      ("mean expression  log m", m)):
            mr = np.corrcoef(fb, y)[0, 1]
            if "mean expression" in nm:
                print(f"  {nm:28s} {mr:+11.3f} {'—':>23s}")
            else:
                pr = np.corrcoef(resid(fb, m), resid(y, m))[0, 1]
                print(f"  {nm:28s} {mr:+11.3f} {pr:+23.3f}")

        # effect size: bound vs unbound, matched on expression
        BOUND, UNB = fb > 0.5, fb < 0.05
        fr = np.exp(f); br = np.exp(b)
        print(f"\n  bound (>50% of experiments) n={BOUND.sum()}, unbound (<5%) n={UNB.sum()}")
        print(f"  RAW medians   f {np.median(fr[BOUND])/np.median(fr[UNB]):.2f}x   "
              f"B {np.median(br[BOUND])/np.median(br[UNB]):.2f}x   "
              f"expr {np.median(np.exp(m)[BOUND])/np.median(np.exp(m)[UNB]):.2f}x"
              "   <- dominated by PRC2 targets being lowly expressed")

        bins = np.quantile(m, np.linspace(0, 1, 9))
        rows = []
        print(f"\n  EXPRESSION-MATCHED (within expression octiles, >=15 genes each side):")
        print(f"    {'expr bin':>20s} {'n_bound':>8s} {'f ratio':>9s} {'B ratio':>9s}")
        for i in range(len(bins) - 1):
            sel = (m >= bins[i]) & (m < bins[i + 1])
            if (sel & BOUND).sum() < 15 or (sel & UNB).sum() < 15:
                continue
            rf = np.median(fr[sel & BOUND]) / np.median(fr[sel & UNB])
            rb = np.median(br[sel & BOUND]) / np.median(br[sel & UNB])
            rows.append((rf, rb, (sel & BOUND).sum()))
            print(f"    {np.exp(bins[i]):8.2f}-{np.exp(bins[i+1]):<11.2f} "
                  f"{(sel & BOUND).sum():8d} {rf:9.2f} {rb:9.2f}")
        if rows:
            w = np.array([r[2] for r in rows], float)
            w /= w.sum()
            print(f"\n    weighted: f ratio = {np.sum(w*[r[0] for r in rows]):.2f}, "
                  f"B ratio = {np.sum(w*[r[1] for r in rows]):.2f}   (1.00 = no effect)")
            print("    Ochiai predicts f ratio < 1 and B ratio > 1.")
        print("\n    NB only the lowest expression bins have enough bound genes to match on —"
              "\n    PRC2 targets are so skewed to low expression that this is a"
              "\n    low-expression-only result.")


if __name__ == "__main__":
    main()
