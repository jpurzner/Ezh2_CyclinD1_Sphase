"""Supp figure (transcript version): G0/G1 commitment in the plane of Cyclin D1 TRANSCRIPT x EZH2 REPRESSION.

Same map as fig_v44_cd_ezh2_commitment.py, but the x-axis is the Cyclin D1 mRNA TRANSCRIPT (the RNA-seq
observable) rather than an abstract drive.  Because Ezh2 represses Cyclin D1 *transcription*, the actual
transcript present = (transcriptional drive) x (1 - repression); commitment fires when the resulting
Cyclin D1 protein (proportional to transcript) clears the threshold, which the MB CDK-inhibitor tone
raises above GNP.  So the cycle/arrest boundary in (transcript drive, repression r) is  drive = mRNA*/(1-r).

We measure from the model (mean Cd_mRNA):
  mRNA*_MB , mRNA*_GNP : transcript cycling threshold in MB vs GNP CDK-inhibitor context
  mRNA_MB              : MB operating Cyclin D1 transcript (full dynamic model)
  F                    : EZH2i de-repression fold (transcript)  ->  r_phys = 1 - 1/F
  DROP                 : vismodegib transcript drop (MB+HHi / MB)

The transcription drive is swept by scaling the three Cyclin D1 transcription terms
(k_Cd_tx_basal / k_Cd_tx_Gli_max / k_Cd_tx_MYCN) by a common factor, with EZH2i=1 (repression OFF).

Run:  ./venv/bin/python simulations/fig_v44_cdmrna_ezh2_commitment.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

MB_TONE  = dict(p16=0.15, p18=1.5, ksyp21=0.004)
GNP_TONE = dict(p16=0.0,  p18=0.4, ksyp21=0.002)
TX0 = dict(basal=0.3086, gli=46.31, mycn=35.22)          # default Cyclin D1 transcription terms


def _runner():
    rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
    rr.integrator.setValue("relative_tolerance", 1e-6)
    return rr


def simulate(rr, tx=1.0, SHH=0.5, f=0.1, mycn=2.86, ezh2i=0.0, hhi=0.0,
             p16=0.0, p18=0.4, ksyp21=0.002):
    """mean settled Cd_mRNA (transcript) + division count.  tx scales the transcription terms;
    EZH2i=1 turns Ezh2 repression OFF, giving a clean transcript-drive scan."""
    for te_end, te_pts in ((12000, 48000), (8000, 32000), (6000, 24000)):
        for atol in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
            rr.reset()
            rr['SHH'] = SHH; rr['Ptch1_copy_number'] = f; rr['MYCN_amplification'] = mycn
            rr['HHi'] = hhi; rr['EZH2i'] = ezh2i
            rr['p16'] = p16; rr['p18'] = p18; rr['kSyP21'] = ksyp21
            rr['k_Cd_tx_basal'] = tx * TX0['basal']
            rr['k_Cd_tx_Gli_max'] = tx * TX0['gli']
            rr['k_Cd_tx_MYCN'] = tx * TX0['mycn']
            rr.integrator.setValue("absolute_tolerance", atol)
            try: rr.integrator.setValue("maximum_num_steps", 300000)
            except Exception: pass
            try:
                r = rr.simulate(0, te_end, te_pts, selections=["time", "MPF", "Cd_mRNA"]); break
            except Exception:
                r = None
        if r is not None:
            break
    if r is None:
        return np.nan, -1
    t = r['time']; m = t >= 4000; tt = t[m]; dt = tt[1] - tt[0]
    pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=int(200 / dt))
    return float(np.mean(r['Cd_mRNA'][m])), len(pk)


def threshold_scan(tone, f, mycn, label):
    """Scan the transcription drive (EZH2i=1, repression OFF) and return (transcript, cyc) points and
    the transcript threshold set by this CDK-inhibitor tone."""
    rr = _runner()
    pts = []
    for tx in [0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.18, 0.27, 0.40, 0.60, 0.85, 1.2]:
        m, div = simulate(rr, tx=tx, f=f, mycn=mycn, ezh2i=1.0, **tone)
        if np.isnan(m):
            continue
        pts.append((m, div >= 2))
        print(f"   {label:4s}  tx={tx:5.2f}  Cd_mRNA={m:7.3f}  div={div:2d}  {'CYCLE' if div>=2 else 'arrest'}")
    pts.sort()
    arr = [m for m, c in pts if not c]; cyc = [m for m, c in pts if c]
    if arr and cyc:
        hi_arr = max(arr); lo_cyc = min([c for c in cyc if c > hi_arr] or cyc)
        thr = 0.5 * (hi_arr + lo_cyc)
    else:
        thr = np.nan
    return pts, thr


print("=" * 70)
print("Measuring commitment thresholds in TRANSCRIPT units (Ezh2 repression OFF):")
mb_pts,  mRNA_MB_thr  = threshold_scan(MB_TONE,  f=0.1, mycn=2.86, label="MB")
gnp_pts, mRNA_GNP_thr = threshold_scan(GNP_TONE, f=1.0, mycn=1.0,  label="GNP")
print(f"\n  mRNA* (MB tone) = {mRNA_MB_thr:.3f}   mRNA* (GNP tone) = {mRNA_GNP_thr:.3f}"
      f"   -> MB brake raises transcript threshold {mRNA_MB_thr/mRNA_GNP_thr:.2f}x")

print("\nMeasuring MB operating transcript + EZH2i de-repression fold (full dynamic model):")
rr_full = _runner()
mRNA_MB, _     = simulate(rr_full, f=0.1, mycn=2.86, **MB_TONE)                 # repression ON
mRNA_MB_nr, _  = simulate(rr_full, f=0.1, mycn=2.86, ezh2i=1.0, **MB_TONE)      # repression OFF, same drive
mRNA_HHi, _    = simulate(rr_full, f=0.1, mycn=2.86, hhi=1.0, **MB_TONE)        # vismodegib
F = (mRNA_MB_nr / mRNA_MB) if (mRNA_MB and not np.isnan(mRNA_MB_nr) and mRNA_MB > 1e-6) else 2.5
F = float(np.clip(F, 1.3, 6.0))
r_phys = 1.0 - 1.0 / F
DROP = float(mRNA_HHi / mRNA_MB) if (mRNA_MB and not np.isnan(mRNA_HHi)) else 0.14
print(f"  mRNA_MB={mRNA_MB:.3f}  mRNA_MB(no repression)={mRNA_MB_nr:.3f}  mRNA_MB+HHi={mRNA_HHi:.3f}")
print(f"  de-repression fold F={F:.2f} -> r_phys={r_phys:.2f}   vismo transcript drop={DROP:.3f} ({100*(1-DROP):.0f}%)")

# ---- operating points in (transcript drive, repression).  transcript_present = drive*(1-r). ----
drive_MB = mRNA_MB / (1 - r_phys)
drive_HHi = DROP * drive_MB
OPS = {
    "MB":             dict(drive=drive_MB,  r=r_phys, cyc=True),
    "MB +HHi":        dict(drive=drive_HHi, r=r_phys, cyc=False),
    "MB +HHi +EZH2i": dict(drive=drive_HHi, r=0.0,    cyc=(drive_HHi > mRNA_MB_thr)),
}
print("\nOperating points (drive, r, transcript=drive*(1-r)):")
for n, o in OPS.items():
    tr = o['drive'] * (1 - o['r'])
    print(f"  {n:16s} drive={o['drive']:.3f}  r={o['r']:.2f}  transcript={tr:.3f}  "
          f"{'CYCLE' if tr>mRNA_MB_thr else 'arrest'}  (mRNA*_MB={mRNA_MB_thr:.3f})")

# ======================================================================
# Figure
# ======================================================================
fig, ax = plt.subplots(figsize=(9.4, 7.2))
DMAX = drive_MB * 1.4
XMIN = max(0.3, mRNA_GNP_thr * 0.6)
rr_grid = np.linspace(0.0, 0.85, 260)
dd_grid = np.geomspace(XMIN, DMAX, 260)
DD, RR = np.meshgrid(dd_grid, rr_grid)
NET = DD * (1 - RR)
ax.contourf(DD, RR, (NET > mRNA_MB_thr).astype(float), levels=[-0.5, 0.5, 1.5],
            colors=["#fdecea", "#e8f6ee"], zorder=0)

rline = np.linspace(0.0, 0.85, 200)
ax.plot(mRNA_MB_thr / (1 - rline), rline, color="#2c3e50", lw=2.4, zorder=3,
        label=f"MB commitment boundary (mRNA*={mRNA_MB_thr:.2f})")
ax.plot(mRNA_GNP_thr / (1 - rline), rline, color="#7f8c8d", lw=1.8, ls="--", zorder=3,
        label=f"GNP boundary (no brake, mRNA*={mRNA_GNP_thr:.2f}; {mRNA_MB_thr/mRNA_GNP_thr:.1f}× lower)")

ax.text(mRNA_MB_thr * 0.62, 0.74, "ARREST  (G0)", color="#c0392b", fontsize=12, fontweight="bold", ha="center")
ax.text(DMAX * 0.45, 0.06, "CYCLING", color="#1e8449", fontsize=12, fontweight="bold", ha="center")

PCOL = {"MB": "#4A1486", "MB +HHi": "#C2185B", "MB +HHi +EZH2i": "#1e8449"}
LOFF = {"MB": (8, 8), "MB +HHi": (-6, 12), "MB +HHi +EZH2i": (8, -14)}
P = {}
for n, o in OPS.items():
    x, y = o['drive'], o['r']; P[n] = (x, y)
    ax.scatter([x], [y], s=170, c=PCOL[n], edgecolor="white", lw=2, zorder=6)
    dx, dy = LOFF[n]
    ax.annotate(n, (x, y), textcoords="offset points", xytext=(dx, dy),
                fontsize=10, fontweight="bold", color=PCOL[n], zorder=7)

# (1) vismodegib = mitogen withdrawal: transcript drive collapses MB -> MB+HHi; Ezh2 dips in step.
t = np.linspace(0, 1, 60)
mw_x = drive_MB * (drive_HHi / drive_MB) ** t
mw_y = r_phys - 0.30 * np.sin(np.pi * t)
ax.plot(mw_x, mw_y, color="#E08214", lw=2.6, zorder=4)
ax.add_patch(FancyArrowPatch((mw_x[-3], mw_y[-3]), (mw_x[-1], mw_y[-1]), arrowstyle="-|>",
             mutation_scale=18, color="#E08214", lw=2.6, zorder=5))
ax.annotate("vismodegib = mitogen withdrawal\n(transcript ↓; Ezh2 dips in step →\na few extra divisions before G0)",
            (mw_x[30], mw_y[30]), textcoords="offset points", xytext=(-6, -42),
            fontsize=8.6, color="#9C5400", fontweight="bold", ha="center", zorder=7)

# (2) EZH2i rescue: from MB+HHi straight DOWN (repression -> 0) across the boundary
ax.add_patch(FancyArrowPatch(P["MB +HHi"], P["MB +HHi +EZH2i"], arrowstyle="-|>",
             mutation_scale=18, color="#1e8449", lw=2.8, zorder=5))
ax.annotate("EZH2i rescue\n(repression ↓ → transcript\nde-repressed → re-enter cycle)",
            P["MB +HHi"], textcoords="offset points", xytext=(14, -34),
            fontsize=8.6, color="#1e8449", fontweight="bold", zorder=7)

# CDK4/6i: off-plane
ax.text(DMAX * 0.98, 0.83,
        "MB +CDK4/6i:\nVmax kinase block downstream\nof Cyclin D1 — off this plane,\nnot re-crossable → no rescue",
        fontsize=8.0, color="#7b241c", ha="right", va="top",
        bbox=dict(boxstyle="round,pad=0.35", fc="#fdedec", ec="#c0392b", lw=1.2), zorder=7)

ax.set_xscale("log")
ax.set_xlim(XMIN, DMAX); ax.set_ylim(0, 0.85)
ax.set_xlabel("Cyclin D1 transcript  (Cd mRNA, mitogen/Gli/MYCN-set, a.u.)", fontsize=11)
ax.set_ylabel("Ezh2 repression strength  (0 = none → strong)", fontsize=11)
ax.set_title("Commitment in the Cyclin D1-transcript × Ezh2-repression plane\n"
             "transcript present = drive × (1 − repression); MB brake raises the threshold",
             fontsize=12, fontweight="bold")
ax.legend(loc="lower right", fontsize=9, framealpha=0.92)
ax.grid(True, alpha=0.15)
fig.text(0.5, -0.02,
         "x = Cyclin D1 transcriptional drive (mRNA before Ezh2 repression); the transcript actually "
         "present is drive × (1 − repression). Two perturbations cross the SAME boundary: vismodegib "
         "lowers the drive (out of cycle), EZH2i lifts repression (back in).",
         ha="center", fontsize=8.0, color="#555", style="italic")

plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), "fig_v44_cdmrna_ezh2_commitment")
fig.savefig(out + ".png", dpi=165, bbox_inches="tight")
fig.savefig(out + ".pdf", bbox_inches="tight")
print("\nwrote", out + ".png /.pdf")
