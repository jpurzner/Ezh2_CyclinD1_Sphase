"""Figure: experimental conditions as INPUT PARAMETERS, and where they sit in PARAMETER SPACE.

Panel A — input-parameter matrix: the context + drug knobs we set for each condition
  (SHH, functional Ptch1 f, MYCN amplification, p16, and the drug switches HHi/EZH2i/HU/CDK4/6i).
Panel B — parameter space: the commitment decision is governed by CyclinD1 (the drive) vs p16 (the
  competitive CDK4/6 brake). Each condition is placed at its model steady-state CyclinD1 and its p16;
  the cycle/arrest boundary is computed by scanning CyclinD1 (via k_Cd_translation) x p16. Shows how
  vismodegib drops MB below the p16-raised threshold and EZH2i (CyclinD1 up) rescues it.

Run:  ./venv/bin/python simulations/fig_v44_condition_paramspace.py
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

rr = te.loada(build_model_v44())
rr.integrator.setValue("relative_tolerance", 1e-6)
DEF_kPhRbCd = rr['kPhRbCd']; DEF_ktl = rr['k_Cd_translation']

# condition -> input knobs.  f = functional Ptch1 (Ptch1_copy_number); cdk46i sets kPhRbCd=0;
# starve sets k_Cd_translation=0.  p16: 0 in GNP (EZH2-silenced), 3.3 in MB; ~1 for the Ptch+/- intermediate.
CONDS = [
    ("GNP +SHH",          dict(SHH=0.5, f=1.0, mycn=1.0,  p16=0.0, hhi=0, ezh2i=0, hu=0, cdk46i=0, starve=0)),
    ("GNP +HHi",          dict(SHH=0.5, f=1.0, mycn=1.0,  p16=0.0, hhi=1, ezh2i=0, hu=0, cdk46i=0, starve=0)),
    ("GNP +EZH2i",        dict(SHH=0.5, f=1.0, mycn=1.0,  p16=0.0, hhi=0, ezh2i=1, hu=0, cdk46i=0, starve=0)),
    ("GNP serum-starve",  dict(SHH=0.5, f=1.0, mycn=1.0,  p16=0.0, hhi=0, ezh2i=0, hu=0, cdk46i=0, starve=1)),
    ("Ptch+/- (P7)",      dict(SHH=0.5, f=0.5, mycn=1.38, p16=1.0, hhi=0, ezh2i=0, hu=0, cdk46i=0, starve=0)),
    ("MB",                dict(SHH=0.5, f=0.1, mycn=2.86, p16=3.3, hhi=0, ezh2i=0, hu=0, cdk46i=0, starve=0)),
    ("MB +HHi",           dict(SHH=0.5, f=0.1, mycn=2.86, p16=3.3, hhi=1, ezh2i=0, hu=0, cdk46i=0, starve=0)),
    ("MB +HHi +EZH2i",    dict(SHH=0.5, f=0.1, mycn=2.86, p16=3.3, hhi=1, ezh2i=1, hu=0, cdk46i=0, starve=0)),
    ("MB +CDK4/6i",       dict(SHH=0.5, f=0.1, mycn=2.86, p16=3.3, hhi=0, ezh2i=0, hu=0, cdk46i=1, starve=0)),
    ("MB +HU",            dict(SHH=0.5, f=0.1, mycn=2.86, p16=3.3, hhi=0, ezh2i=0, hu=1, cdk46i=0, starve=0)),
]


def simulate(SHH, f, mycn, p16, hhi, ezh2i, hu, cdk46i, starve, ktl_override=None):
    for atol in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
        rr.reset()
        rr['kPhRbCd'] = DEF_kPhRbCd; rr['k_Cd_translation'] = DEF_ktl
        rr['SHH'] = SHH; rr['Ptch1_copy_number'] = f; rr['MYCN_amplification'] = mycn; rr['p16'] = p16
        rr['GDC0449'] = hhi; rr['EZH2i'] = ezh2i; rr['HU'] = hu
        if cdk46i: rr['kPhRbCd'] = 0.0
        if starve: rr['k_Cd_translation'] = 0.0
        if ktl_override is not None: rr['k_Cd_translation'] = ktl_override
        rr.integrator.setValue("absolute_tolerance", atol)
        try:
            rr.integrator.setValue("maximum_num_steps", 300000)
        except Exception:
            pass
        try:
            r = rr.simulate(0, 12000, 48000, selections=["time", "MPF", "Cd"]); break
        except Exception:
            r = None
    if r is None:
        return np.nan, -1
    t = r['time']; m = t >= 4000; tt = t[m]; dt = tt[1] - tt[0]
    pk, _ = find_peaks(r['MPF'][m], prominence=0.15, distance=int(200 / dt))
    return float(np.mean(r['Cd'][m])), len(pk)


def cd_short(SHH, f, mycn, p16, hhi, ezh2i, hu, cdk46i, starve, **_):
    """Robust steady-state CyclinD1 from a short settle (less likely to crash than the full run)."""
    for atol in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
        rr.reset()
        rr['kPhRbCd'] = DEF_kPhRbCd; rr['k_Cd_translation'] = DEF_ktl
        rr['SHH'] = SHH; rr['Ptch1_copy_number'] = f; rr['MYCN_amplification'] = mycn; rr['p16'] = p16
        rr['GDC0449'] = hhi; rr['EZH2i'] = ezh2i; rr['HU'] = hu
        if cdk46i: rr['kPhRbCd'] = 0.0
        if starve: rr['k_Cd_translation'] = 0.0
        rr.integrator.setValue("absolute_tolerance", atol)
        try:
            r = rr.simulate(0, 3000, 9000, selections=["time", "Cd"]); break
        except Exception:
            r = None
    if r is None:
        return np.nan
    t = r['time']; return float(np.mean(r['Cd'][t >= 1500]))


def boundary_cd(p16):
    """Critical CyclinD1 at this p16 (interp of the cycle/arrest grid)."""
    if not boundary:
        return np.nan
    bx = [b[0] for b in boundary]; bp = [b[1] for b in boundary]
    return float(np.interp(p16, bp, bx))


# --- cycle/arrest boundary in (CyclinD1, p16): scan ktl x p16 in the MB context (computed first) ---
print("computing cycle/arrest boundary ...")
P16_GRID = [0.0, 0.75, 1.5, 2.25, 3.0, 4.0, 5.0]
KTL_GRID = [0.03, 0.05, 0.07, 0.09, 0.12, 0.16, 0.22, 0.32, 0.5, 0.75, 1.05]
bx, by, bc = [], [], []   # Cd, p16, cycling(bool) of grid points
boundary = []             # (Cd_crit, p16)
for p in P16_GRID:
    pts = []
    for ktl in KTL_GRID:
        cd, div = simulate(0.5, 0.1, 2.86, p, 0, 0, 0, 0, 0, ktl_override=ktl)
        if np.isnan(cd):
            continue
        bx.append(cd); by.append(p); bc.append(div >= 2)
        pts.append((cd, div >= 2))
    pts.sort()
    arr = [cd for cd, c in pts if not c]
    cyc = [cd for cd, c in pts if c]
    if arr and cyc:
        boundary.append((0.5 * (max(arr) + min([c for c in cyc if c > max(arr)] or cyc)), p))
boundary.sort(key=lambda z: z[1])

print("computing per-condition CyclinD1 + cycling ...")
rows = []
for name, k in CONDS:
    cd, div = simulate(**k)
    if np.isnan(cd):                       # full cycling run was too stiff -> robust short Cd + boundary call
        cd = cd_short(**k)
        cyc = (not k['cdk46i']) and (not k['starve']) and (cd > boundary_cd(k['p16']))
        note = "  (Cd via short settle; cycle inferred from boundary)"
    else:
        cyc = div >= 2
        note = ""
    rows.append(dict(name=name, **k, cd=cd, cyc=cyc))
    print(f"  {name:18s} Cd={cd:6.2f}  {'CYCLE' if cyc else 'arrest'}{note}")

# ======================================================================
# Figure
# ======================================================================
fig = plt.figure(figsize=(16, 7.2))
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.22)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])

# ---------------- Panel A: input-parameter matrix ----------------
cols = [("SHH", "SHH"), ("Ptch1 f", "f"), ("MYCN", "mycn"), ("p16", "p16"),
        ("HHi", "hhi"), ("EZH2i", "ezh2i"), ("HU", "hu"), ("CDK4/6i", "cdk46i")]
hue = {"SHH": "#f1c40f", "f": "#27ae60", "mycn": "#e67e22", "p16": "#8e44ad",
       "hhi": "#e74c3c", "ezh2i": "#e74c3c", "hu": "#e74c3c", "cdk46i": "#e74c3c"}
# normalizer per column: intensity 0..1 (perturbation strength)
def intensity(key, v):
    if key == "SHH":   return 0.0 if v >= 0.5 else 0.85          # SHH present = baseline (light)
    if key == "f":     return (1.0 - v) / 0.9                    # loss of functional Ptch1
    if key == "mycn":  return (v - 1.0) / (2.86 - 1.0)
    if key == "p16":   return v / 3.3
    return float(v)                                              # binary drugs

nC, nR = len(cols), len(rows)
for j, (clabel, key) in enumerate(cols):
    for i, row in enumerate(rows):
        v = row[key]
        inten = np.clip(intensity(key, v), 0, 1)
        from matplotlib.colors import to_rgb
        base = np.array(to_rgb(hue[key]))
        cell = (1 - inten) * np.ones(3) + inten * base
        axA.add_patch(plt.Rectangle((j, nR - 1 - i), 1, 1, fc=cell, ec="white", lw=2))
        # text
        if key in ("hhi", "ezh2i", "hu", "cdk46i"):
            txt = "✓" if v else "·"
        elif key == "SHH":
            txt = f"{v:g}"
        elif key == "f":
            txt = f"{v:g}"
        else:
            txt = f"{v:g}"
        axA.text(j + 0.5, nR - 1 - i + 0.5, txt, ha="center", va="center",
                 fontsize=9, color="#222" if inten < 0.55 else "white",
                 fontweight="bold" if inten > 0.05 else "normal")
# row labels (conditions) + cycle/arrest tag
for i, row in enumerate(rows):
    y = nR - 1 - i + 0.5
    axA.text(-0.15, y, row['name'], ha="right", va="center", fontsize=9.5,
             fontweight="bold" if row['name'] in ("GNP +SHH", "MB") else "normal")
    tag = "●" if row['cyc'] else "○"
    axA.text(nC + 0.2, y, tag, ha="center", va="center", fontsize=11,
             color="#1e8449" if row['cyc'] else "#c0392b")
# column headers
for j, (clabel, key) in enumerate(cols):
    axA.text(j + 0.5, nR + 0.15, clabel, ha="center", va="bottom", fontsize=9, fontweight="bold", rotation=0)
axA.text(nC + 0.2, nR + 0.15, "cyc?", ha="center", va="bottom", fontsize=8.5, style="italic", color="#555")
# group separator: context | drugs (after column 4)
axA.plot([4, 4], [0, nR], color="#34495e", lw=1.8)
axA.text(2.0, -0.55, "context inputs", ha="center", fontsize=9, style="italic", color="#34495e")
axA.text(6.0, -0.55, "drug switches", ha="center", fontsize=9, style="italic", color="#34495e")
axA.set_xlim(-3.0, nC + 0.6); axA.set_ylim(-0.9, nR + 0.7); axA.axis("off")
axA.set_title("A   Input parameters per condition", fontsize=12.5, fontweight="bold", loc="left")
axA.text(-3.0, -0.55, "serum-starve = k_Cd_transl→0;  CDK4/6i = kPhRbCd→0", fontsize=7.3, color="#777", ha="left")

# ---------------- Panel B: parameter space (CyclinD1 x p16) ----------------
bx = np.array(bx); by = np.array(by); bc = np.array(bc)
axB.scatter(bx[bc], by[bc],  s=26, c="#abebc4", edgecolor="none", zorder=1)
axB.scatter(bx[~bc], by[~bc], s=26, c="#f5b7b1", edgecolor="none", zorder=1)
if len(boundary) >= 2:
    bX = [b[0] for b in boundary]; bY = [b[1] for b in boundary]
    axB.plot(bX, bY, color="#2c3e50", lw=2.0, ls="--", zorder=2, label="cycle / arrest boundary")
axB.text(0.17, 4.6, "ARREST", color="#c0392b", fontsize=11, fontweight="bold", ha="center")
axB.text(5.2, 0.35, "CYCLING", color="#1e8449", fontsize=11, fontweight="bold", ha="center")

# place conditions (skip CDK4/6i -- off-plane; skip starve Cd~0)
lbl_off = {"GNP +SHH": (-2, 11), "GNP +HHi": (-2, 11), "GNP +EZH2i": (6, 8),
           "Ptch+/- (P7)": (9, -3), "MB": (-2, 12), "MB +HHi": (-8, -15),
           "MB +HHi +EZH2i": (2, 11), "MB +HU": (6, -13)}
pts = {}
for row in rows:
    if row['cdk46i'] or row['starve'] or np.isnan(row['cd']):
        continue
    x, y = max(row['cd'], 0.05), row['p16']
    pts[row['name']] = (x, y)
    col = "#1e8449" if row['cyc'] else "#c0392b"
    axB.scatter([x], [y], s=130, c=col, edgecolor="white", lw=1.8, zorder=5)
    dx, dy = lbl_off.get(row['name'], (8, 8))
    axB.annotate(row['name'], (x, y), textcoords="offset points", xytext=(dx, dy),
                 fontsize=8.5, fontweight="bold", color=col, zorder=6)
# rescue arrow MB+HHi -> MB+HHi+EZH2i
if "MB +HHi" in pts and "MB +HHi +EZH2i" in pts:
    p1, p2 = pts["MB +HHi"], pts["MB +HHi +EZH2i"]
    axB.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=16,
                                  color="#4A1486", lw=2.2, zorder=4,
                                  connectionstyle="arc3,rad=-0.2"))
    axB.text(1.6, 3.65, "EZH2i rescue\n(CyclinD1 ↑ ~3×)", color="#4A1486", fontsize=8.2,
             fontweight="bold", ha="center")
# CDK4/6i annotation (off-plane)
axB.text(6.6, 4.3, "MB +CDK4/6i:\nkPhRbCd→0 (Vmax block)\noff-plane → arrest,\nno EZH2i rescue",
         fontsize=7.6, color="#7b241c", ha="left", va="top",
         bbox=dict(boxstyle="round,pad=0.3", fc="#fdedec", ec="#c0392b", lw=1))

axB.set_xscale("log")
axB.set_xlim(0.1, 11); axB.set_ylim(-0.3, 5.3)
axB.set_xlabel("CyclinD1 level  (the drive)  — log scale", fontsize=10.5)
axB.set_ylabel("p16  (competitive CDK4/6 brake)", fontsize=10.5)
axB.set_title("B   Parameter space: CyclinD1 vs p16 commitment threshold", fontsize=12.5, fontweight="bold", loc="left")
axB.legend(loc="lower left", fontsize=8.5, framealpha=0.9)
axB.grid(True, which="both", alpha=0.18)
axB.text(0.115, -0.05, "points = deterministic mean cell · tumour CyclinD1/p27 heterogeneity straddles the "
         "boundary → fractional response (MB+HHi 22% → +EZH2i 69% cycling)",
         fontsize=7.2, color="#555", ha="left", va="top", style="italic")

fig.suptitle("v44 model — experimental conditions: input parameters & where they sit in the proliferation parameter space",
             fontsize=13, fontweight="bold", y=1.0)
plt.tight_layout(rect=[0, 0, 1, 0.97])
out = os.path.join(os.path.dirname(__file__), "fig_v44_condition_paramspace")
fig.savefig(out + ".png", dpi=160, bbox_inches="tight")
fig.savefig(out + ".pdf", bbox_inches="tight")
print("wrote", out + ".png /.pdf")
