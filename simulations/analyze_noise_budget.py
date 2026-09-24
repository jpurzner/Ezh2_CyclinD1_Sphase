#!/usr/bin/env python3
"""Noise budget: what CV and what correlation time should each v44 species carry?

Reads data/priors/turnover_burst_priors.csv (built by build_turnover_priors.py) and, for
every gene with enough data, evaluates the two-stage bursty-transcription result

    CV_p^2  ~  1/<p>  +  (1 + B)/<m> * g_p/(g_p + g_m)                      (intrinsic)
    tau_corr ~ max(1/g_m, 1/g_p_eff)

with g_p_eff = g_p + mu the removal rate INCLUDING dilution (mu = ln2/T_cc). It then puts
those next to the DECAY RATE THE MODEL ACTUALLY USES for the same species, which is what
decides whether injected noise survives to the cell-cycle engine or is averaged away.

The comparison is the point: a species whose model lifetime is much shorter than its real
lifetime cannot carry realistic noise, because the downstream (hour-scale) engine integrates
over many independent fluctuations.

Usage:  python3 simulations/analyze_noise_budget.py [--tcc 16]
"""
import argparse
import csv
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV = os.path.join(ROOT, "data", "priors", "turnover_burst_priors.csv")
LN2 = math.log(2.0)

# ------------------------------------------------------------------------------------------
# What v44 actually uses, in MODEL UNITS (1/min -- the Heldt core and everything bolted onto
# it run in real-time minutes; see validate_v44.run() "Real-time minutes").
#
# `kind`: 'protein' or 'mrna'. `note` records where the constant came from, which is exactly
# the class (a)/(b)/(c) question the dilution audit (T4) has to answer.
MODEL_RATES = {
    # gene      species        kind        k [1/min]   note
    "E2f1":   [("E2f",        "protein",   0.05,     "Heldt kDeE2f (fitted in cycling cells -> class (a))")],
    "Cdkn1a": [("P21",        "protein",   0.0025,   "Heldt kDeP21 basal; plus Skp2/Cdt2 terms (class (a))")],
    "Ccne1":  [("Ce",         "protein",   0.004,    "Heldt kDeCe basal (+kDeCeCa*Ca) (class (a))")],
    "Ccna2":  [("Ca",         "protein",   0.01,     "Heldt kDeCa basal (+kDeCaC1*C1) (class (a))")],
    "Fbxo5":  [("E1",         "protein",   0.0005,   "Heldt kDeE1 (class (a))")],
    "Ccnd1":  [("Cd",         "protein",   1.0,      "k_Cd_deg -- INHERITED FROM v42, NOT rescaled to minutes"),
               ("Cd_mRNA",    "mrna",      0.8,      "k_Cd_mRNA_deg -- INHERITED FROM v42, NOT rescaled")],
    "Ezh2":   [("EZH2",       "protein",   0.0011685, "kDeEZ, fitted to HU-arrest + G0 withdrawal (class (a)/(c))"),
               ("EZH2m",      "mrna",      0.0098,   "kDeEZm (class (c))")],
    "Gli1":   [("Gli1",       "protein",   0.8,      "k_Gli1_deg -- v42 units"),
               ("Gli1_mRNA",  "mrna",      0.8,      "k_Gli1_mRNA_deg -- v42 units")],
    "Mycn":   [("MYCN",       "protein",   1.0,      "k_MYCN_deg -- v42 units")],
    "Ptch1":  [("Ptch1_free", "protein",   0.5,      "k_Ptch1_deg -- v42 units"),
               ("Ptch1_mRNA", "mrna",      0.8,      "k_Ptch1_mRNA_deg -- v42 units")],
    "Smo":    [("Smo_active", "protein",   1.2,      "k_Smo_inact -- v42 units (inactivation, not degradation)")],
}


def load_priors():
    """gene -> {quantity: {cell_type: value}}"""
    d = {}
    with open(CSV) as fh:
        for r in csv.DictReader(fh):
            if r["value"] == "NA":
                continue
            try:
                v = float(r["value"])
            except ValueError:
                continue
            d.setdefault(r["gene"], {}).setdefault(r["quantity"], {})[r["cell_type"]] = v
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tcc", type=float, default=16.0, help="cell-cycle length in hours")
    a = ap.parse_args()
    mu = LN2 / a.tcc  # 1/h dilution rate

    pri = load_priors()
    print(f"# Noise budget   (T_cc = {a.tcc} h, mu = {mu:.4f}/h)")
    print(f"# CV_p^2 = 1/<p> + (1+B)/<m> * g_p_eff/(g_p_eff+g_m);  g_p_eff includes dilution\n")
    hdr = (f"{'gene':8s} {'<p>':>9s} {'<m>':>6s} {'B':>5s} {'t_m':>6s} {'t_p':>7s} "
           f"{'filt':>5s} {'CV_p':>6s} {'tau_c/h':>8s} | {'model species':13s} {'model tau':>10s} {'ratio':>7s}")
    print(hdr); print("-" * len(hdr))

    for gene in sorted(pri):
        q = pri[gene]
        def g1(name, prefer=None):
            m = q.get(name, {})
            if not m: return None
            if prefer:
                for k in m:
                    if prefer in k: return m[k]
            return next(iter(m.values()))

        p   = g1("protein_copy_number")
        tp  = g1("protein_half_life", prefer="NIH3T3")      # cycling-cell value
        tpn = q.get("protein_half_life", {}).get("mouse primary neuron")
        tm  = g1("mrna_half_life")
        B   = g1("burst_size_B", prefer="fibroblast")
        m_allele = g1("mrna_mean_per_allele", prefer="fibroblast")
        m   = g1("mrna_copy_number") or (2 * m_allele if m_allele else None)

        rows = MODEL_RATES.get(gene, [(None, None, None, None)])
        if not (tm and B and m):
            for sp, kind, k, _ in rows:
                if sp:
                    tau_model = 1.0 / k / 60.0
                    print(f"{gene:8s} {'insufficient prior data':>48s} | {sp:13s} {tau_model:9.2f}h {'--':>7s}")
            continue

        g_m = LN2 / tm
        g_p = (LN2 / tp) if tp else None
        g_p_eff = (g_p + mu) if g_p else mu      # dilution-only if no measured half-life
        filt = g_p_eff / (g_p_eff + g_m)
        cv2 = ((1.0 / p) if p else 0.0) + (1.0 + B) / m * filt
        tau_c = max(1.0 / g_m, 1.0 / g_p_eff)

        first = True
        for sp, kind, k, _ in rows:
            lead = (f"{gene:8s} {p if p else float('nan'):9.0f} {m:6.1f} {B:5.2f} "
                    f"{tm:6.1f} {(tp if tp else float('nan')):7.1f} {filt:5.2f} "
                    f"{math.sqrt(cv2):6.2f} {tau_c:8.1f}") if first else " " * 63
            first = False
            if sp is None:
                print(lead + " | (no v44 species)")
                continue
            tau_model = 1.0 / k / 60.0                      # hours
            tau_bio = (1.0 / g_m) if kind == "mrna" else (1.0 / g_p_eff)
            ratio = tau_bio / tau_model
            flag = "  <-- model >20x too fast" if ratio > 20 else ("  <-- model too fast" if ratio > 3 else "")
            print(f"{lead} | {sp:13s} {tau_model:9.2f}h {ratio:6.1f}x{flag}")

    print("\n# Legend: t_m/t_p = mRNA/protein half-life (h). filt = g_p_eff/(g_p_eff+g_m), the fraction")
    print("# of mRNA noise reaching the protein. tau_c = protein-noise correlation time.")
    print("# ratio = (biological correlation time) / (model relaxation time) for that species.")

    # ---------------------------------------------------------------- sensitivity band
    # The mRNA half-life is the load-bearing and least trustworthy number: it sets BOTH the
    # filter factor and the correlation time, and Schwanhausser's estimates are systematically
    # long (Ccnd1 = 14.0 h here, vs the ~0.5-1 h usually quoted from actinomycin-D chase).
    # Rather than impute a value, sweep it and show what the choice buys.
    print(f"\n\n# mRNA-lifetime sensitivity for the burst candidates (T_cc = {a.tcc} h)")
    print("# NOT an imputation -- this is the range the conclusion has to survive.\n")
    band = [("Ccnd1", 0.5, 14.0), ("Cdkn1a", 0.5, 9.0), ("Mycn", 0.5, 9.0), ("Ezh2", 5.0, 9.85)]
    print(f"{'gene':8s} {'t_p/h':>6s} {'<m>':>6s} {'B':>6s} | " +
          "  ".join(f"t_m={t:<5g}" for t in (0.5, 1, 2, 4, 9, 14)))
    print(f"{'':8s} {'':6s} {'':6s} {'':6s} | " +
          "  ".join(f"{'CV/tau':<10s}" for _ in range(6)))
    for gene, _lo, _hi in band:
        q = pri.get(gene, {})
        def g1(name, prefer=None):
            mm = q.get(name, {})
            if not mm: return None
            if prefer:
                for k in mm:
                    if prefer in k: return mm[k]
            return next(iter(mm.values()))
        p_ = g1("protein_copy_number")
        tp = g1("protein_half_life", prefer="NIH3T3") or g1("protein_half_life")
        B  = g1("burst_size_B", prefer="fibroblast") or g1("burst_size_B")
        ma = g1("mrna_mean_per_allele", prefer="fibroblast") or g1("mrna_mean_per_allele")
        m_ = g1("mrna_copy_number") or (2 * ma if ma else None)
        if not (B and m_ and tp):
            print(f"{gene:8s} {'-- needs B, <m> and a protein half-life --':>50s}")
            continue
        g_p_eff = LN2 / tp + mu
        cells = []
        for tm_ in (0.5, 1, 2, 4, 9, 14):
            g_m_ = LN2 / tm_
            f_ = g_p_eff / (g_p_eff + g_m_)
            cv = math.sqrt(((1.0 / p_) if p_ else 0.0) + (1.0 + B) / m_ * f_)
            tau = max(1.0 / g_m_, 1.0 / g_p_eff)
            cells.append(f"{cv:.2f}/{tau:<5.1f}")
        print(f"{gene:8s} {tp:6.1f} {m_:6.1f} {B:6.2f} | " + "  ".join(f"{c:<10s}" for c in cells))
    print("\n# cell = CV_p / tau_c[h].  Where CV barely moves across the row, the mRNA half-life")
    print("# does not matter and drift (option 3) is safe. Where tau_c moves from ~1 h to ~20 h,")
    print("# the SAME snapshot CV means either fast fluctuation or near-frozen heterogeneity --")
    print("# opposite predictions for basin escape, and the choice must be justified, not defaulted.")


if __name__ == "__main__":
    main()
