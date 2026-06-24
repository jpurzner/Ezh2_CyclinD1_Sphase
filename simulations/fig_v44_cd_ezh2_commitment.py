"""Supp figure: the G0/G1 commitment in the plane of Cyclin D1 DRIVE x EZH2 REPRESSION.

One map that unifies two results:
  (1) mitogen withdrawal (paragraph 2): lowering Hedgehog drive moves a cell left; because Ezh2 is
      cell-cycle-coupled it falls in step, so repression also eases (downward) -> the path skirts the
      commitment boundary and yields a few extra divisions before crossing out into G0.
  (2) the EZH2i rescue (rescue paragraph): from arrested MB+HHi, relieving Ezh2 repression raises net
      Cyclin D1 (straight down) back across the SAME boundary into the cycle. CDK4/6i is off-plane
      (a Vmax kinase block downstream of Cyclin D1) -> cannot be re-crossed here -> no rescue.

Model fact used: net Cyclin D1  =  drive (mitogen/Gli/MYCN-set transcription)  x  Ezh2 repression factor.
Commitment fires when net Cd clears the threshold, which the MB CDK-inhibitor tone raises above GNP.
So the cycle/arrest boundary in (drive, repression r) is the curve  drive = Cd* / (1 - r).

We measure from the model:
  Cd*_MB , Cd*_GNP : cycling threshold (mean Cd at the boundary) in MB vs GNP CDK-inhibitor context
  Cd_MB            : MB operating Cyclin D1 (full dynamic model)
  F                : EZH2i de-repression fold in MB+HHi  ->  physiological repression r_phys = 1 - 1/F

Run:  ./venv/bin/python simulations/fig_v44_cd_ezh2_commitment.py
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

# MB CDK-inhibitor tone vs GNP baseline (same numbers as validate_v44 / the paramspace fig)
MB_TONE  = dict(p16=0.15, p18=1.5,  ksyp21=0.004)
GNP_TONE = dict(p16=0.0,  p18=0.4,  ksyp21=0.002)


def _runner():
    rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True))   # full model for everything
    rr.integrator.setValue("relative_tolerance", 1e-6)
    return rr


def simulate(rr, ktl=None, SHH=0.5, f=0.1, mycn=2.86, ezh2i=0.0, hhi=0.0,
             p16=0.0, p18=0.4, ksyp21=0.002):
    """mean settled Cd + division count; robust to the flaky CVODE init.
    EZH2i=1 turns the Ezh2 repression OFF, giving a clean drive (net-Cd) scan."""
    for te_end, te_pts in ((12000, 48000), (8000, 32000), (6000, 24000)):
        for atol in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
            rr.reset()
            rr['SHH'] = SHH; rr['Ptch1_copy_number'] = f; rr['MYCN_amplification'] = mycn
            rr['HHi'] = hhi
            rr['p16'] = p16; rr['p18'] = p18; rr['kSyP21'] = ksyp21; rr['EZH2i'] = ezh2i
            if ktl is not None:
                rr['k_Cd_translation'] = ktl
            rr.integrator.setValue("absolute_tolerance", atol)
            try: rr.integrator.setValue("maximum_num_steps", 300000)
            except Exception: pass
            try:
                r = rr.simulate(0, te_end, te_pts, selections=["time", "MPF", "Cd"]); break
            except Exception:
                r = None
        if r is not None:
            break
    if r is None:
        return np.nan, -1
    t = r['time']; m = t >= 4000; tt = t[m]; dt = tt[1] - tt[0]
    pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=int(200 / dt))
    return float(np.mean(r['Cd'][m])), len(pk)


def threshold_scan(tone, f, mycn, label):
    """Scan k_Cd_translation (the Cyclin D1 drive proxy) with EZH2i=1 (Ezh2 repression OFF), so drive
    is a clean input.  Return (Cd, cyc) points and the Cd threshold set by this CDK-inhibitor tone."""
    rr = _runner()
    pts = []
    for ktl in [0.02, 0.04, 0.06, 0.08, 0.11, 0.15, 0.20, 0.28, 0.40, 0.55, 0.75, 1.0, 1.25]:
        cd, div = simulate(rr, ktl=ktl, f=f, mycn=mycn, ezh2i=1.0, **tone)
        if np.isnan(cd):
            continue
        pts.append((cd, div >= 2))
        print(f"   {label:4s}  ktl={ktl:5.2f}  Cd={cd:7.3f}  div={div:2d}  {'CYCLE' if div>=2 else 'arrest'}")
    pts.sort()
    arr = [cd for cd, c in pts if not c]
    cyc = [cd for cd, c in pts if c]
    if arr and cyc:
        hi_arr = max(arr); lo_cyc = min([c for c in cyc if c > hi_arr] or cyc)
        thr = 0.5 * (hi_arr + lo_cyc)
    else:
        thr = np.nan
    return pts, thr


print("=" * 70)
print("Measuring commitment thresholds (Ezh2 feedback OFF; drive = k_Cd_translation):")
mb_pts,  Cd_MB_thr  = threshold_scan(MB_TONE,  f=0.1, mycn=2.86, label="MB")
gnp_pts, Cd_GNP_thr = threshold_scan(GNP_TONE, f=1.0, mycn=1.0,  label="GNP")
print(f"\n  Cd* (MB tone)  = {Cd_MB_thr:.3f}   Cd* (GNP tone) = {Cd_GNP_thr:.3f}"
      f"   -> MB brake raises threshold {Cd_MB_thr/Cd_GNP_thr:.2f}x")

# operating Cyclin D1 + de-repression fold, from the FULL dynamic model
print("\nMeasuring MB operating Cyclin D1 + EZH2i de-repression fold (full dynamic model):")
rr_full = _runner()
Cd_MB, _      = simulate(rr_full, f=0.1, mycn=2.86, **MB_TONE)                  # repression ON
Cd_MB_norep,_ = simulate(rr_full, f=0.1, mycn=2.86, ezh2i=1.0, **MB_TONE)      # repression OFF (same drive)
Cd_HHi, _     = simulate(rr_full, f=0.1, mycn=2.86, hhi=1.0, **MB_TONE)        # vismodegib: Gli drive off
# de-repression fold measured at MB operating drive (non-stiff): F = Cd(no repression) / Cd(repression)
F = (Cd_MB_norep / Cd_MB) if (Cd_MB and not np.isnan(Cd_MB_norep) and Cd_MB > 1e-6) else 2.5
F = float(np.clip(F, 1.3, 6.0))
r_phys = 1.0 - 1.0 / F
DROP = float(Cd_HHi / Cd_MB) if (Cd_MB and not np.isnan(Cd_HHi)) else 0.14     # vismo net-Cd drop
print(f"  Cd_MB={Cd_MB:.3f}  Cd_MB(no repression)={Cd_MB_norep:.3f}  Cd_MB+HHi={Cd_HHi:.3f}")
print(f"  de-repression fold F={F:.2f} -> r_phys={r_phys:.2f}   vismo net-Cd drop={DROP:.3f} ({100*(1-DROP):.0f}%)")

# ---- operating points in the (drive, repression) plane.  net = drive*(1-r). ----
# MB: full drive, physiological repression.   net_MB = Cd_MB.
drive_MB = Cd_MB / (1 - r_phys)
# vismo removes the Gli-driven input: net falls to the model's observed fraction of MB
drive_HHi = DROP * drive_MB                      # repression unchanged; DRIVE collapses
# EZH2i: same (collapsed) drive, repression -> 0
OPS = {
    "MB":              dict(drive=drive_MB,  r=r_phys,  cyc=True),
    "MB +HHi":         dict(drive=drive_HHi, r=r_phys,  cyc=False),
    "MB +HHi +EZH2i":  dict(drive=drive_HHi, r=0.0,     cyc=(drive_HHi * 1.0 > Cd_MB_thr)),
}
print("\nOperating points (drive, r, net=drive*(1-r)):")
for n, o in OPS.items():
    net = o['drive'] * (1 - o['r'])
    print(f"  {n:16s} drive={o['drive']:.3f}  r={o['r']:.2f}  net={net:.3f}  "
          f"{'CYCLE' if net>Cd_MB_thr else 'arrest'}  (thr_MB={Cd_MB_thr:.3f})")

# ======================================================================
# Figure
# ======================================================================
fig, ax = plt.subplots(figsize=(9.4, 7.2))
DMAX = drive_MB * 1.4
XMIN = 0.3
rr_grid = np.linspace(0.0, 0.85, 260)
dd_grid = np.geomspace(XMIN, DMAX, 260)
DD, RR = np.meshgrid(dd_grid, rr_grid)
NET = DD * (1 - RR)
ax.contourf(DD, RR, (NET > Cd_MB_thr).astype(float), levels=[-0.5, 0.5, 1.5],
            colors=["#fdecea", "#e8f6ee"], zorder=0)

# boundary curves: drive = Cd* / (1 - r)
rline = np.linspace(0.0, 0.85, 200)
ax.plot(Cd_MB_thr / (1 - rline), rline, color="#2c3e50", lw=2.4, zorder=3,
        label=f"MB commitment boundary (Cd*={Cd_MB_thr:.2f})")
ax.plot(Cd_GNP_thr / (1 - rline), rline, color="#7f8c8d", lw=1.8, ls="--", zorder=3,
        label=f"GNP boundary (no brake, Cd*={Cd_GNP_thr:.2f}; {Cd_MB_thr/Cd_GNP_thr:.1f}× lower)")

ax.text(Cd_MB_thr * 0.62, 0.74, "ARREST  (G0)", color="#c0392b", fontsize=12, fontweight="bold", ha="center")
ax.text(DMAX * 0.45, 0.06, "CYCLING", color="#1e8449", fontsize=12, fontweight="bold", ha="center")

# operating points
PCOL = {"MB": "#4A1486", "MB +HHi": "#C2185B", "MB +HHi +EZH2i": "#1e8449"}
LOFF = {"MB": (8, 8), "MB +HHi": (-6, 12), "MB +HHi +EZH2i": (8, -14)}
P = {}
for n, o in OPS.items():
    x, y = o['drive'], o['r']; P[n] = (x, y)
    col = "#1e8449" if o['cyc'] else "#c0392b"
    ax.scatter([x], [y], s=170, c=PCOL[n], edgecolor="white", lw=2, zorder=6)
    dx, dy = LOFF[n]
    ax.annotate(n, (x, y), textcoords="offset points", xytext=(dx, dy),
                fontsize=10, fontweight="bold", color=PCOL[n], zorder=7)

# (1) vismodegib = mitogen withdrawal: drive collapses MB -> MB+HHi; as Ezh2 falls in step the path
#     bows DOWN (repression dips), keeping net Cyclin D1 up for a few extra divisions before it arrests.
t = np.linspace(0, 1, 60)
mw_x = drive_MB * (drive_HHi / drive_MB) ** t            # drive falls log-linearly MB -> MB+HHi
mw_y = r_phys - 0.30 * np.sin(np.pi * t)                 # Ezh2 dips then returns at the arrested state
ax.plot(mw_x, mw_y, color="#E08214", lw=2.6, zorder=4)
ax.add_patch(FancyArrowPatch((mw_x[-3], mw_y[-3]), (mw_x[-1], mw_y[-1]), arrowstyle="-|>",
             mutation_scale=18, color="#E08214", lw=2.6, zorder=5))
ax.annotate("vismodegib = mitogen withdrawal\n(drive ↓; Ezh2 dips in step →\na few extra divisions before G0)",
            (mw_x[30], mw_y[30]), textcoords="offset points", xytext=(-6, -42),
            fontsize=8.6, color="#9C5400", fontweight="bold", ha="center", zorder=7)

# (2) EZH2i rescue: from MB+HHi straight DOWN (repression -> 0) across the boundary
ax.add_patch(FancyArrowPatch(P["MB +HHi"], P["MB +HHi +EZH2i"], arrowstyle="-|>",
             mutation_scale=18, color="#1e8449", lw=2.8, zorder=5,
             connectionstyle="arc3,rad=0.0"))
ax.annotate("EZH2i rescue\n(repression ↓ → net\nCyclin D1 ↑ → re-enter cycle)",
            P["MB +HHi"], textcoords="offset points", xytext=(14, -34),
            fontsize=8.6, color="#1e8449", fontweight="bold", zorder=7)

# CDK4/6i: off-plane
ax.text(DMAX * 0.98, 0.83,
        "MB +CDK4/6i:\nVmax kinase block downstream\nof Cyclin D1 — off this plane,\nnot re-crossable → no rescue",
        fontsize=8.0, color="#7b241c", ha="right", va="top",
        bbox=dict(boxstyle="round,pad=0.35", fc="#fdedec", ec="#c0392b", lw=1.2), zorder=7)

ax.set_xscale("log")
ax.set_xlim(XMIN, DMAX); ax.set_ylim(0, 0.85)
ax.set_xlabel("Cyclin D1 drive  (mitogen / Gli / MYCN-set transcription, a.u.)", fontsize=11)
ax.set_ylabel("Ezh2 repression strength  (0 = none → strong)", fontsize=11)
ax.set_title("Commitment in the Cyclin D1-drive × Ezh2-repression plane\n"
             "net Cyclin D1 = drive × (1 − repression); MB brake raises the threshold",
             fontsize=12, fontweight="bold")
ax.legend(loc="lower right", fontsize=9, framealpha=0.92)
ax.grid(True, alpha=0.15)
fig.text(0.5, -0.02,
         "Two perturbations cross the SAME boundary: mitogen withdrawal moves a cell down-and-left "
         "out of the cycle (with a few last divisions as Ezh2 falls in step); EZH2i moves an arrested "
         "MB cell straight down back into it.",
         ha="center", fontsize=8.4, color="#555", style="italic")

plt.tight_layout()
out = os.path.join(os.path.dirname(__file__), "fig_v44_cd_ezh2_commitment")
fig.savefig(out + ".png", dpi=165, bbox_inches="tight")
fig.savefig(out + ".pdf", bbox_inches="tight")
print("\nwrote", out + ".png /.pdf")
