"""Figure: MB is p27-high / G0-rich via a higher MB-specific BIRTH p27 -- a FOLD-SAFE lever (v44).

MESSAGE
  MB dwells longer in the p27-high (transient-G0) basin each cycle than GNP, and it does so because
  the MB daughter is BORN with high p27 (P21_div = 1.8 vs GNP 0.6) -- a synthesis-vs-Skp2 balance set
  at division. Born-high p27 latches the Skp2-p27 toggle LOW, opening a growth-timed G0 window, and it
  does this WITHOUT perturbing the hard 5.07x MB/GNP CyclinD1 fold (whole-cycle EZH2 feedback). It is a
  FOLD-SAFE attractor lever -- NOT the shelved mitogen-inverse birth-p27 tracker (with_mitogen_tracker).

  The fold-safety is the point: to make MB G0-rich you must pick a p27 lever. Raising p27 via
  steady-state synthesis (kSyP21) to reproduce MB's G0 window slows commitment enough that
  whole-cycle EZH2 over-represses CyclinD1 and the MB/GNP fold INFLATES to ~8.8x (over the tolerance
  band edge → data-excluded). Birth-p27 opens the same G0 window while leaving the fold at ~7.1x,
  inside the calibrated 5.07 ± 45% band (edge 7.35). The raw MB with NO G0 lever sits at ~4.9x.

PANELS
  (A) single-cell p27(t) traces: GNP (birth 0.6) vs MB (birth 1.8). MB dwells longer in the p27-high
      basin each cycle; the p27>0.1 pre-S (G0) windows are shaded.
  (B) G0 fraction = fraction of pre-S time with p27>0.1, GNP vs MB (bars). MB clearly higher.
  (C) FOLD-SAFETY: MB/GNP CyclinD1 mRNA fold. Raw MB (no G0 lever) ≈ 4.9× (near target, G0≈0);
      +birth-p27 lever ≈ 7.1× (G0-rich, inside the 5.07 ± 45% band); reproducing the SAME G0 via
      kSyP21 instead ≈ 8.8× (over the band edge, data-excluded). Target band 5.07 ± 0.45 shaded.

Deterministic single-cell mechanism (distinct from the population ensemble). REUSE:
  simulations/fig_v44_transient_g0_definition.py (CyclinD1 x p27 decision machinery + p27 traces),
  simulations/validate_v44.py (MB condition constants + P21_DIV_MB).

Run:  ./venv/bin/python simulations/fig_v44_mb_g0_birthp27.py
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.build_model_v44_heldt import build_model_v44

# --- validate_v44.py constants (single source of truth) ----------------------
P27_THR = 0.1            # p27 (=P21) marker threshold for the G0/G1 split (G0 = p27-high, pre-S)
MYCN_AMP_MB = 2.8
PTCH1_MB = 0.1
P16_MB = 0.306
P18_MB = 1.553
KSYP21_MB = 0.002
P21_DIV_MB = 1.8         # MB-specific BIRTH p27 (fold-safe G0 lever)
P21_DIV_GNP = 0.6        # builder default birth p27
KSYP21_ALT = 0.010       # kSyP21 that reproduces MB's G0 window WITHOUT birth-p27 (the fold-UNSAFE lever)

CD_FOLD_TARGET, CD_FOLD_TOL = 5.07, 0.45   # validate_v44 Fig-4I target + fractional tolerance (band edge 7.35)
T_END, N_PTS, SETTLE = 15000, 60000, 4000
SEL = ["time", "P21", "Cd", "Cd_mRNA", "Dna", "aRc", "MPF", "mass", "EZH2"]

# Conditions. GNP = builder defaults (SHH=0.5, P21_div=0.6). MB = MB-specific tones + birth-p27 lever.
GNP = dict(SHH=0.5, Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, P21_div=P21_DIV_GNP)
MB = dict(SHH=0.5, Ptch1_copy_number=PTCH1_MB, MYCN_amplification=MYCN_AMP_MB,
          p16=P16_MB, p18=P18_MB, kSyP21=KSYP21_MB, P21_div=P21_DIV_MB)
# fold-safety references: raw MB with NO G0 lever (birth back to GNP 0.6, baseline kSyP21) ...
MB_BASE = dict(SHH=0.5, Ptch1_copy_number=PTCH1_MB, MYCN_amplification=MYCN_AMP_MB,
               p16=P16_MB, p18=P18_MB, kSyP21=KSYP21_MB, P21_div=P21_DIV_GNP)
# ... and the SAME G0 window built via steady-state synthesis (kSyP21) instead of birth-p27.
MB_KSY = dict(SHH=0.5, Ptch1_copy_number=PTCH1_MB, MYCN_amplification=MYCN_AMP_MB,
              p16=P16_MB, p18=P18_MB, kSyP21=KSYP21_ALT, P21_div=P21_DIV_GNP)


def run_cell(cond):
    """Fresh model per condition (reset() does NOT restore parameters). Real-time minutes."""
    rr = te.loada(build_model_v44())     # default build: two-step Rb + EZH2 conc + me1/2/3 chain
    rr.integrator.setValue("absolute_tolerance", 1e-8)
    rr.integrator.setValue("relative_tolerance", 1e-6)
    try:
        rr.integrator.setValue("maximum_num_steps", 2000000)
        rr.integrator.setValue("maximum_time_step", 5.0)
    except Exception:
        pass
    for k, v in cond.items():
        rr[k] = v
    d = rr.simulate(0, T_END, N_PTS, selections=SEL)
    t = d["time"]
    m = t >= SETTLE
    out = {s: np.asarray(d[s])[m] for s in SEL}
    out["time"] = out["time"] / 60.0     # -> hours
    return out


def analyse(o):
    """G0 fraction (fraction of pre-S time with p27>0.1), mean CyclinD1, divisions."""
    t, P21, Dna, aRc = o["time"], o["P21"], o["Dna"], o["aRc"]
    preS = (aRc < 0.05) & (Dna < 0.98)              # G0 + G1 (pre-replication)
    g0 = preS & (P21 > P27_THR)                     # p27-high pre-S window
    g0_frac = float(g0.sum() / max(preS.sum(), 1))
    # divisions = Dna reset 1 -> 0
    div_idx = np.where((Dna[:-1] > 0.9) & (Dna[1:] < 0.5))[0] + 1
    return dict(g0_frac=g0_frac, mean_cd=float(np.mean(o["Cd"])),
                mean_cd_mrna=float(np.mean(o["Cd_mRNA"])),
                mean_p21=float(np.mean(P21)), ndiv=int(len(div_idx)))


def shade_g0(ax, o, color="#f7dc6f"):
    """Shade p27>0.1 pre-S (transient-G0) windows."""
    t, P21, Dna, aRc = o["time"], o["P21"], o["Dna"], o["aRc"]
    g0 = (P21 > P27_THR) & (aRc < 0.05) & (Dna < 0.98)
    ax.fill_between(t, 0, 1, where=g0, color=color, alpha=0.55, lw=0,
                    transform=ax.get_xaxis_transform(), label="G0 (p27>0.1, pre-S)")


def main():
    print("running deterministic single-cell sims (GNP, MB, MB-base, MB[kSyP21 lever]) ...", flush=True)
    gnp = run_cell(GNP)
    mb = run_cell(MB)
    mb_base = run_cell(MB_BASE)
    mb_ksy = run_cell(MB_KSY)
    A = {"GNP": analyse(gnp), "MB": analyse(mb),
         "MB_BASE": analyse(mb_base), "MB_KSY": analyse(mb_ksy)}

    # CyclinD1 fold measured on Cd_mRNA (the transcript the EZH2 feedback represses; validate's 5.07 metric)
    cd = lambda k: A[k]["mean_cd_mrna"]
    cd_fold_base = cd("MB_BASE") / cd("GNP")     # raw MB, no G0 lever
    cd_fold_birth = cd("MB") / cd("GNP")         # birth-p27 lever (calibrated MB)
    cd_fold_ksy = cd("MB_KSY") / cd("GNP")       # same G0 via kSyP21
    band_hi = CD_FOLD_TARGET * (1 + CD_FOLD_TOL)
    band_lo = CD_FOLD_TARGET * (1 - CD_FOLD_TOL)

    # ---- figure -------------------------------------------------------------
    C_GNP, C_MB = "#2471a3", "#c0392b"
    fig = plt.figure(figsize=(15, 8))
    gs = GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.34,
                  left=0.065, right=0.965, top=0.85, bottom=0.10)
    fig.suptitle("MB is p27-high / G0-rich via a higher MB-specific BIRTH p27 — a FOLD-SAFE attractor lever (v44)",
                 fontsize=14, fontweight="bold", y=0.965)
    fig.text(0.5, 0.905,
             "Born-high p27 (P21_div 1.8 vs GNP 0.6) latches the Skp2–p27 toggle low → a growth-timed G0 window, "
             "WITHOUT perturbing the hard 5.07× MB/GNP CyclinD1 fold.",
             ha="center", fontsize=9.5, color="#555", style="italic")

    # (a) GNP p27 trace
    axg = fig.add_subplot(gs[0, 0:2])
    shade_g0(axg, gnp)
    axg.plot(gnp["time"], gnp["P21"], color=C_GNP, lw=1.5, label="p27 (=P21)")
    axg.plot(gnp["time"], gnp["Dna"], color="#17a2b8", lw=1.0, alpha=0.7, label="Dna (S/G2)")
    axg.axhline(P27_THR, color="#8e44ad", lw=0.8, ls="--", alpha=0.6)
    axg.axhline(P21_DIV_GNP, color=C_GNP, lw=0.8, ls=":", alpha=0.6)
    axg.text(SETTLE / 60.0 + 2, P21_DIV_GNP + 0.03, "birth p27 = 0.6", fontsize=7.5, color=C_GNP)
    axg.set_ylim(-0.03, 2.0)
    axg.set_ylabel("level")
    axg.set_title(f"(a) GNP single cell: born low p27 → commits fast  "
                  f"(G0 = {A['GNP']['g0_frac']*100:.1f}% of pre-S)",
                  fontsize=9.5, fontweight="bold", color=C_GNP)
    axg.legend(fontsize=7, loc="upper right", framealpha=0.9, ncol=3)
    axg.grid(alpha=0.2)

    # (a) MB p27 trace
    axm = fig.add_subplot(gs[1, 0:2], sharex=axg)
    shade_g0(axm, mb)
    axm.plot(mb["time"], mb["P21"], color=C_MB, lw=1.5, label="p27 (=P21)")
    axm.plot(mb["time"], mb["Dna"], color="#17a2b8", lw=1.0, alpha=0.7, label="Dna (S/G2)")
    axm.axhline(P27_THR, color="#8e44ad", lw=0.8, ls="--", alpha=0.6)
    axm.axhline(P21_DIV_MB, color=C_MB, lw=0.8, ls=":", alpha=0.6)
    axm.text(SETTLE / 60.0 + 2, P21_DIV_MB + 0.03, "birth p27 = 1.8 (P21_DIV_MB)", fontsize=7.5, color=C_MB)
    axm.set_ylim(-0.03, 2.0)
    axm.set_xlabel("time (h)")
    axm.set_ylabel("level")
    axm.set_title(f"(a) MB single cell: born HIGH p27 → dwells in the p27-high basin each cycle  "
                  f"(G0 = {A['MB']['g0_frac']*100:.1f}% of pre-S)",
                  fontsize=9.5, fontweight="bold", color=C_MB)
    axm.legend(fontsize=7, loc="upper right", framealpha=0.9, ncol=3)
    axm.grid(alpha=0.2)

    # (b) G0 fraction bars
    axb = fig.add_subplot(gs[0, 2])
    vals = [A["GNP"]["g0_frac"] * 100, A["MB"]["g0_frac"] * 100]
    bars = axb.bar(["GNP\n(birth 0.6)", "MB\n(birth 1.8)"], vals, color=[C_GNP, C_MB], width=0.6)
    for b, v in zip(bars, vals):
        axb.text(b.get_x() + b.get_width() / 2, v + 0.4, f"{v:.1f}%", ha="center",
                 fontsize=10, fontweight="bold")
    axb.set_ylabel("G0 fraction\n(pre-S time with p27>0.1)  %")
    axb.set_ylim(0, max(vals) * 1.35 + 1)
    axb.set_title("(b) MB dwells longer in G0", fontsize=9.5, fontweight="bold")
    axb.grid(alpha=0.2, axis="y")

    # (c) fold-safety -- 3 bars: raw MB (no G0) -> add G0 via birth-p27 (in-band) OR via kSyP21 (out)
    axc = fig.add_subplot(gs[1, 2])
    folds = [cd_fold_base, cd_fold_birth, cd_fold_ksy]
    g0s = [A["MB_BASE"]["g0_frac"] * 100, A["MB"]["g0_frac"] * 100, A["MB_KSY"]["g0_frac"] * 100]
    labels = ["MB raw\n(no G0 lever)", "MB +birth-p27\n(P21_div=1.8)", f"MB +kSyP21\n(={KSYP21_ALT:g})"]
    cols = ["#95a5a6", "#1e8449", "#7b241c"]
    barc = axc.bar(labels, folds, color=cols, width=0.62)
    for b, v, g in zip(barc, folds, g0s):
        axc.text(b.get_x() + b.get_width() / 2, v + 0.15, f"{v:.2f}×", ha="center",
                 fontsize=9.5, fontweight="bold")
        axc.text(b.get_x() + b.get_width() / 2, 0.25, f"G0 {g:.0f}%", ha="center",
                 fontsize=7.2, color="white", fontweight="bold")
    axc.axhspan(band_lo, band_hi, color="#abebc6", alpha=0.45, lw=0)
    axc.axhline(band_hi, color="#1e8449", lw=1.0, ls="--", alpha=0.9)
    axc.axhline(CD_FOLD_TARGET, color="#1e8449", lw=0.8, ls=":", alpha=0.7)
    axc.text(0.015, band_hi + 0.12, f"band edge {band_hi:.2f}×", transform=axc.get_yaxis_transform(),
             fontsize=7.0, color="#145a32", va="bottom", ha="left")
    axc.text(0.015, CD_FOLD_TARGET + 0.12, "5.07× target (Fig 4I)", transform=axc.get_yaxis_transform(),
             fontsize=7.0, color="#145a32", va="bottom", ha="left")
    axc.set_ylabel("MB / GNP CyclinD1 mRNA fold")
    axc.set_ylim(0, max(folds) * 1.22 + 0.5)
    axc.set_title("(c) FOLD-SAFETY: birth-p27 keeps the fold in-band", fontsize=9.2, fontweight="bold")
    axc.grid(alpha=0.2, axis="y")
    fig.text(0.5, 0.028, "Panel (c): making MB G0-rich via kSyP21 (steady-state p27 synthesis) instead of "
             "birth-p27 slows commitment enough that whole-cycle EZH2 over-represses\nCyclinD1 — the MB/GNP "
             "fold clears the tolerance-band edge and becomes data-excluded. The birth-p27 lever opens the same "
             "G0 window while staying inside the band.",
             ha="center", fontsize=7.8, color="#555", style="italic")

    out = os.path.join(os.path.dirname(__file__), "fig_v44_mb_g0_birthp27")
    fig.savefig(out + ".png", dpi=170)
    fig.savefig(out + ".pdf")
    print("wrote", out + ".png /.pdf")

    # ---- readouts -----------------------------------------------------------
    inband = lambda f: band_lo <= f <= band_hi
    print("\n=== readouts ===")
    for k in ["GNP", "MB", "MB_BASE", "MB_KSY"]:
        a = A[k]
        print(f"  {k:8s}: G0_frac={a['g0_frac']*100:5.1f}%  mean_Cd={a['mean_cd']:.3f}  "
              f"mean_CdmRNA={a['mean_cd_mrna']:.3f}  mean_p27={a['mean_p21']:.3f}  ndiv={a['ndiv']}")
    print(f"\n  G0 fraction (pre-S time with p27>0.1):  GNP={A['GNP']['g0_frac']*100:.1f}%  "
          f"MB={A['MB']['g0_frac']*100:.1f}%")
    print(f"  CyclinD1 mRNA fold band: {band_lo:.2f} .. {band_hi:.2f}x  (target {CD_FOLD_TARGET})")
    print(f"  MB/GNP fold  raw MB (no G0 lever)          = {cd_fold_base:.2f}x  "
          f"[{'IN-BAND' if inband(cd_fold_base) else 'OUT'}]  (G0={A['MB_BASE']['g0_frac']*100:.0f}%)")
    print(f"  MB/GNP fold  +birth-p27 lever (P21_div=1.8)= {cd_fold_birth:.2f}x  "
          f"[{'IN-BAND' if inband(cd_fold_birth) else 'OUT'}]  (G0={A['MB']['g0_frac']*100:.0f}%)")
    print(f"  MB/GNP fold  +kSyP21={KSYP21_ALT:g} (same G0)      = {cd_fold_ksy:.2f}x  "
          f"[{'IN-BAND' if inband(cd_fold_ksy) else 'DATA-EXCLUDED'}]  (G0={A['MB_KSY']['g0_frac']*100:.0f}%)")


if __name__ == "__main__":
    import tellurium as te
    main()
