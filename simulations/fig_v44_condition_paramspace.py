"""Figure: experimental conditions as INPUT PARAMETERS, and where they sit in PARAMETER SPACE.

Panel A — input-parameter matrix: the context + drug knobs we set for each condition
  (SHH, functional Ptch1 f, MYCN amplification, and the CDK-inhibitor tone p16 + p21, plus the
  drug switches HHi/EZH2i/HU/CDK4/6i).
Panel B — parameter space: the commitment decision is governed by CyclinD1 (the drive) vs the
  CDK-inhibitor tone (the brake: p16 raises the CyclinD1->Rb half-max via CDK4/6, p21 inhibits CDK2;
  both are co-elevated GNP->MB). Each condition is placed at its model steady-state CyclinD1 and its
  brake tone; the cycle/arrest boundary is computed by scanning CyclinD1 (k_Cd_translation) x tone.
  Shows how vismodegib drops MB below the brake-raised threshold and EZH2i (CyclinD1 up) rescues it.

Run:  ./venv/bin/python simulations/fig_v44_condition_paramspace.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.colors import to_rgb
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

rr = te.loada(build_model_v44())
rr.integrator.setValue("relative_tolerance", 1e-6)
DEF_kPhRbCd = rr['kPhRbCd']; DEF_ktl = rr['k_Cd_translation']
P16_MB, KSY_BASE, KSY_MB = 1.2, 0.002, 0.004     # MB brake = p16 1.2 + p21 (kSyP21) 2x baseline

# brake "tone" tau: GNP=0 (p16=0, p21 baseline), MB=1 (p16=1.2, p21 2x).  p16(tau)=1.2*tau, kSyP21(tau)=0.002*(1+tau)
def p16_of(t):  return P16_MB * t
def ksy_of(t):  return KSY_BASE * (1.0 + t)

# condition -> input knobs.  tone = combined CDK-inhibitor tone; f = functional Ptch1.
CONDS = [
    ("GNP +SHH",         dict(SHH=0.5, f=1.0, mycn=1.0,  tone=0.0, hhi=0, ezh2i=0, hu=0, cdk46i=0, starve=0)),
    ("GNP +HHi",         dict(SHH=0.5, f=1.0, mycn=1.0,  tone=0.0, hhi=1, ezh2i=0, hu=0, cdk46i=0, starve=0)),
    ("GNP +EZH2i",       dict(SHH=0.5, f=1.0, mycn=1.0,  tone=0.0, hhi=0, ezh2i=1, hu=0, cdk46i=0, starve=0)),
    ("GNP serum-starve", dict(SHH=0.5, f=1.0, mycn=1.0,  tone=0.0, hhi=0, ezh2i=0, hu=0, cdk46i=0, starve=1)),
    ("Ptch+/- (P7)",     dict(SHH=0.5, f=0.5, mycn=1.38, tone=0.4, hhi=0, ezh2i=0, hu=0, cdk46i=0, starve=0)),
    ("MB",               dict(SHH=0.5, f=0.1, mycn=2.86, tone=1.0, hhi=0, ezh2i=0, hu=0, cdk46i=0, starve=0)),
    ("MB +HHi",          dict(SHH=0.5, f=0.1, mycn=2.86, tone=1.0, hhi=1, ezh2i=0, hu=0, cdk46i=0, starve=0)),
    ("MB +HHi +EZH2i",   dict(SHH=0.5, f=0.1, mycn=2.86, tone=1.0, hhi=1, ezh2i=1, hu=0, cdk46i=0, starve=0)),
    ("MB +CDK4/6i",      dict(SHH=0.5, f=0.1, mycn=2.86, tone=1.0, hhi=0, ezh2i=0, hu=0, cdk46i=1, starve=0)),
    ("MB +HU",           dict(SHH=0.5, f=0.1, mycn=2.86, tone=1.0, hhi=0, ezh2i=0, hu=1, cdk46i=0, starve=0)),
]


def _set(SHH, f, mycn, tone, hhi, ezh2i, hu, cdk46i, starve):
    rr['kPhRbCd'] = DEF_kPhRbCd; rr['k_Cd_translation'] = DEF_ktl
    rr['SHH'] = SHH; rr['Ptch1_copy_number'] = f; rr['MYCN_amplification'] = mycn
    rr['p16'] = p16_of(tone); rr['kSyP21'] = ksy_of(tone)
    rr['GDC0449'] = hhi; rr['EZH2i'] = ezh2i; rr['HU'] = hu
    if cdk46i: rr['kPhRbCd'] = 0.0
    if starve: rr['k_Cd_translation'] = 0.0


def simulate(ktl_override=None, **k):
    for te_end, te_pts in ((12000, 48000), (8000, 32000)):
        for atol in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
            rr.reset(); _set(**k)
            if ktl_override is not None: rr['k_Cd_translation'] = ktl_override
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


def cd_short(**k):
    for atol in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
        rr.reset(); _set(**k)
        rr.integrator.setValue("absolute_tolerance", atol)
        try:
            r = rr.simulate(0, 3000, 9000, selections=["time", "Cd"]); break
        except Exception:
            r = None
    if r is None:
        return np.nan
    t = r['time']; return float(np.mean(r['Cd'][t >= 1500]))


# ---- cycle/arrest boundary in (CyclinD1, tone): scan ktl x tone in the MB-genetics context ----
print("computing cycle/arrest boundary ...")
TONE_GRID = [0.0, 0.25, 0.5, 0.75, 1.0, 1.3]
KTL_GRID = [0.03, 0.05, 0.07, 0.09, 0.12, 0.16, 0.22, 0.32, 0.5, 0.75, 1.05]
bx, byy, bc = [], [], []
boundary = []
for tn in TONE_GRID:
    pts = []
    for ktl in KTL_GRID:
        cd, div = simulate(ktl_override=ktl, SHH=0.5, f=0.1, mycn=2.86, tone=tn,
                           hhi=0, ezh2i=0, hu=0, cdk46i=0, starve=0)
        if np.isnan(cd):
            continue
        bx.append(cd); byy.append(tn); bc.append(div >= 2); pts.append((cd, div >= 2))
    pts.sort()
    arr = [cd for cd, c in pts if not c]; cyc = [cd for cd, c in pts if c]
    if arr and cyc:
        boundary.append((0.5 * (max(arr) + min([c for c in cyc if c > max(arr)] or cyc)), tn))
boundary.sort(key=lambda z: z[1])


def boundary_cd(tone):
    if not boundary:
        return np.nan
    return float(np.interp(tone, [b[1] for b in boundary], [b[0] for b in boundary]))


print("computing per-condition CyclinD1 + cycling ...")
rows = []
for name, k in CONDS:
    cd, div = simulate(**k)
    if np.isnan(cd):
        cd = cd_short(**k)
        cyc = (not k['cdk46i']) and (not k['starve']) and (cd > boundary_cd(k['tone']))
        note = "  (Cd via short settle; cycle from boundary)"
    else:
        cyc = div >= 2; note = ""
    rows.append(dict(name=name, **k, p16=p16_of(k['tone']), ksy=ksy_of(k['tone']), cd=cd, cyc=cyc))
    print(f"  {name:18s} Cd={cd:6.2f}  p16={p16_of(k['tone']):.2f} p21x={ksy_of(k['tone'])/KSY_BASE:.1f}  "
          f"{'CYCLE' if cyc else 'arrest'}{note}")

# ======================================================================
# Figure
# ======================================================================
fig = plt.figure(figsize=(16.5, 7.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1.0], wspace=0.24)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])

# ---------------- Panel A: input-parameter matrix ----------------
cols = [("SHH", "SHH"), ("Ptch1 f", "f"), ("MYCN", "mycn"), ("p16", "p16"), ("p21", "p21x"),
        ("HHi", "hhi"), ("EZH2i", "ezh2i"), ("HU", "hu"), ("CDK4/6i", "cdk46i")]
hue = {"SHH": "#f1c40f", "f": "#27ae60", "mycn": "#e67e22", "p16": "#8e44ad", "p21x": "#7d3c98",
       "hhi": "#e74c3c", "ezh2i": "#e74c3c", "hu": "#e74c3c", "cdk46i": "#e74c3c"}
def intensity(key, v):
    if key == "SHH":   return 0.0 if v >= 0.5 else 0.85
    if key == "f":     return (1.0 - v) / 0.9
    if key == "mycn":  return (v - 1.0) / (2.86 - 1.0)
    if key == "p16":   return v / P16_MB
    if key == "p21x":  return (v - 1.0) / 1.0
    return float(v)
for r in rows:
    r["p21x"] = r["ksy"] / KSY_BASE
nC, nR = len(cols), len(rows)
for j, (clabel, key) in enumerate(cols):
    for i, row in enumerate(rows):
        v = row[key]; inten = np.clip(intensity(key, v), 0, 1)
        cell = (1 - inten) * np.ones(3) + inten * np.array(to_rgb(hue[key]))
        axA.add_patch(plt.Rectangle((j, nR - 1 - i), 1, 1, fc=cell, ec="white", lw=2))
        if key in ("hhi", "ezh2i", "hu", "cdk46i"):
            txt = "✓" if v else "·"
        elif key == "p21x":
            txt = f"{v:g}×"
        else:
            txt = f"{v:g}"
        axA.text(j + 0.5, nR - 1 - i + 0.5, txt, ha="center", va="center", fontsize=8.8,
                 color="#222" if inten < 0.55 else "white", fontweight="bold" if inten > 0.05 else "normal")
for i, row in enumerate(rows):
    y = nR - 1 - i + 0.5
    axA.text(-0.15, y, row['name'], ha="right", va="center", fontsize=9.5,
             fontweight="bold" if row['name'] in ("GNP +SHH", "MB") else "normal")
    axA.text(nC + 0.2, y, "●" if row['cyc'] else "○", ha="center", va="center",
             fontsize=11, color="#1e8449" if row['cyc'] else "#c0392b")
for j, (clabel, key) in enumerate(cols):
    axA.text(j + 0.5, nR + 0.12, clabel, ha="center", va="bottom", fontsize=8.8, fontweight="bold")
axA.text(nC + 0.2, nR + 0.12, "cyc?", ha="center", va="bottom", fontsize=8.3, style="italic", color="#555")
axA.plot([5, 5], [0, nR], color="#34495e", lw=1.8)
axA.text(2.5, -0.55, "context inputs  (p16, p21 = CDK-inhibitor tone)", ha="center", fontsize=8.8, style="italic", color="#34495e")
axA.text(7.0, -0.55, "drug switches", ha="center", fontsize=8.8, style="italic", color="#34495e")
axA.set_xlim(-3.0, nC + 0.6); axA.set_ylim(-0.95, nR + 0.7); axA.axis("off")
axA.set_title("A   Input parameters per condition", fontsize=12.5, fontweight="bold", loc="left")
axA.text(-3.0, -0.9, "p16: competitive CDK4/6 brake (raises K_CdRb) · p21: CDK2 inhibitor (kSyP21, ×baseline) · "
         "serum-starve k_Cd_transl→0 · CDK4/6i kPhRbCd→0", fontsize=6.8, color="#777", ha="left")

# ---------------- Panel B: parameter space (CyclinD1 x CDK-inhibitor tone) ----------------
bx = np.array(bx); byy = np.array(byy); bc = np.array(bc)
axB.scatter(bx[bc], byy[bc], s=26, c="#abebc4", edgecolor="none", zorder=1)
axB.scatter(bx[~bc], byy[~bc], s=26, c="#f5b7b1", edgecolor="none", zorder=1)
if len(boundary) >= 2:
    axB.plot([b[0] for b in boundary], [b[1] for b in boundary], color="#2c3e50", lw=2.0, ls="--",
             zorder=2, label="cycle / arrest boundary")
axB.text(0.17, 1.18, "ARREST", color="#c0392b", fontsize=11, fontweight="bold", ha="center")
axB.text(5.2, 0.08, "CYCLING", color="#1e8449", fontsize=11, fontweight="bold", ha="center")

lbl_off = {"GNP +SHH": (-2, 11), "GNP +HHi": (-2, 11), "GNP +EZH2i": (6, 8),
           "Ptch+/- (P7)": (9, -3), "MB": (-2, 12), "MB +HHi": (-8, -15),
           "MB +HHi +EZH2i": (2, 11), "MB +HU": (6, -13)}
pts = {}
for row in rows:
    if row['cdk46i'] or row['starve'] or np.isnan(row['cd']):
        continue
    x, y = max(row['cd'], 0.05), row['tone']
    pts[row['name']] = (x, y)
    col = "#1e8449" if row['cyc'] else "#c0392b"
    axB.scatter([x], [y], s=130, c=col, edgecolor="white", lw=1.8, zorder=5)
    dx, dy = lbl_off.get(row['name'], (8, 8))
    axB.annotate(row['name'], (x, y), textcoords="offset points", xytext=(dx, dy),
                 fontsize=8.5, fontweight="bold", color=col, zorder=6)
if "MB +HHi" in pts and "MB +HHi +EZH2i" in pts:
    p1, p2 = pts["MB +HHi"], pts["MB +HHi +EZH2i"]
    axB.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=16, color="#4A1486",
                                  lw=2.2, zorder=4, connectionstyle="arc3,rad=-0.18"))
    axB.text(1.5, 0.83, "EZH2i rescue\n(CyclinD1 ↑ ~3×)", color="#4A1486", fontsize=8.2,
             fontweight="bold", ha="center")
axB.text(6.4, 1.18, "MB +CDK4/6i:\nkPhRbCd→0 (Vmax block)\noff-plane → arrest,\nno EZH2i rescue",
         fontsize=7.6, color="#7b241c", ha="left", va="top",
         bbox=dict(boxstyle="round,pad=0.3", fc="#fdedec", ec="#c0392b", lw=1))
# right axis: what the tone means
for tn, lab in [(0.0, "GNP\np16 0 / p21 1×"), (0.4, "Ptch+/-\np16 0.5 / p21 1.4×"), (1.0, "MB\np16 1.2 / p21 2×")]:
    axB.text(12.0, tn, lab, fontsize=7.0, color="#555", va="center", ha="left")

axB.set_xscale("log")
axB.set_xlim(0.1, 11); axB.set_ylim(-0.12, 1.45)
axB.set_xlabel("CyclinD1 level  (the drive)  — log scale", fontsize=10.5)
axB.set_ylabel("CDK-inhibitor tone  (p16 + p21, co-elevated GNP→MB)", fontsize=10.5)
axB.set_title("B   Parameter space: CyclinD1 vs the CDK-inhibitor brake", fontsize=12.5, fontweight="bold", loc="left")
axB.legend(loc="lower left", fontsize=8.5, framealpha=0.9)
axB.grid(True, which="both", alpha=0.18)
axB.text(0.115, -0.20, "points = deterministic mean cell · tumour heterogeneity straddles the boundary "
         "→ fractional response (MB+HHi 23% → +EZH2i 73% cycling)",
         fontsize=7.2, color="#555", ha="left", va="top", style="italic")

fig.suptitle("v44 model — experimental conditions: input parameters & where they sit in the proliferation parameter space",
             fontsize=13, fontweight="bold", y=1.0)
plt.tight_layout(rect=[0, 0, 1, 0.97])
out = os.path.join(os.path.dirname(__file__), "fig_v44_condition_paramspace")
fig.savefig(out + ".png", dpi=160, bbox_inches="tight")
fig.savefig(out + ".pdf", bbox_inches="tight")
print("wrote", out + ".png /.pdf")
