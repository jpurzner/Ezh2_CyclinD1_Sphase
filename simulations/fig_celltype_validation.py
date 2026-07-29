"""Split GNP vs MB validation figure with the cell-type parameter differences shown.

Runs validate_v44 at the cell-type-split operating point (GNP 16h / MB ~23.5h core; JP 2026-07-27:
'MB core genuinely longer via k_mu_cki' + decouple_commit). Left panel = the parameter differences
(only the biological identities differ; the cycle-length split EMERGES from the MB CDKi via k_mu_cki).
Middle/right = per-cell-type validation forest (model/target ratio, tolerance band). Bottom = MB<->GNP
cross-ratios. Writes fig_celltype_validation.pdf + .png.
"""
import os, sys, json, subprocess, tempfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(ROOT, 'venv/bin/python'); VAL = os.path.join(ROOT, 'simulations/validate_v44.py')

# ---- the cell-type-split operating point (two-cyclin recal, cKO-preserving) ----
# These MIRROR the baked _ts_bake defaults (two-cyclin model, 30/32); passing them explicitly
# documents the operating point. If a future recal changes _ts_bake, sync these too.
SPLIT = dict(mu=0.000769, k_mu_cki=0.307, kPhRbCd=0.1403, K_CdRb=0.4032,
             k_Cd_tx_Gli_max=0.152118, k_Cd_tx_basal=0.000529, k_Cd_tx_MYCN=0.031157,
             w_Cd2=0.079448, k_Cd2_bas=2.0, k_Cd2_Gli=2.49, K_Cd2_Gli=0.3)
GNP_CORE_H, MB_CORE_H = 17.2, 23.5   # measured core periods at this point

# palette (status + cell-type identity)
C_GNP, C_MB, C_CROSS = '#2c7fb8', '#c0392b', '#7d5ba6'
C_PASS, C_FAIL = '#2a9d5c', '#d1495b'
INK, MUT, FAINT, GRID = '#2d2d2d', '#666666', '#999999', '#d9d9d9'
BAND = '#e8ede9'
plt.rcParams.update({'font.size': 9, 'font.family': 'DejaVu Sans', 'axes.linewidth': 0.8,
                     'svg.fonttype': 'none', 'pdf.fonttype': 42})

# name -> (celltype, tol, kind, target_override, short_label)
M = {
 'CyclinD1 GNP+HHi/GNP':        ('GNP', 0.40, 'fold', None, 'CyclinD1 +HHi/ctrl'),
 'CyclinD2 GNP+HHi/GNP':        ('GNP', 0.30, 'fold', None, 'CyclinD2 +HHi/ctrl'),
 'MYCN GNP+HHi/GNP':            ('GNP', 0.30, 'fold', None, 'MYCN +HHi/ctrl'),
 'Gli1 GNP+HHi reduction':      ('GNP', 0.05, 'fold', None, 'Gli1 +HHi reduction'),
 'EZH2i CycD1 fold (GNP)':      ('GNP', 0.42, 'fold', None, 'EZH2i CyclinD1 fold'),
 'EZH2 G0/cycling (~0.6)':      ('GNP', 0.40, 'fold', None, 'EZH2 G0/cycling'),
 'EZH2 Palbo mRNA drop (~0.44)':('GNP', 0.42, 'fold', None, 'EZH2 Palbo drop'),
 'Period GNP (~22h)':           ('GNP', 0.30, 'value', 16.0, 'cell-cycle period (h)'),
 'GNP+SHH cycles (>0 div)':     ('GNP', None, 'binary', None, 'proliferates +Shh'),
 'GNP-SHH arrest (0 div)':      ('GNP', None, 'binary', None, 'arrests −Shh'),
 'GNP+HHi arrest (0 div)':      ('GNP', None, 'binary', None, 'arrests +HHi'),
 'GNP serum-starve arrest':     ('GNP', None, 'binary', None, 'arrests serum-starve'),
 'CyclinD1 MB+HHi/MB':          ('MB', 0.45, 'fold', None, 'CyclinD1 +HHi/ctrl'),
 'CyclinD2 MB+HHi/MB':          ('MB', 0.45, 'fold', None, 'CyclinD2 +HHi/ctrl'),
 'MYCN MB+HHi/MB':              ('MB', 0.25, 'fold', None, 'MYCN +HHi/ctrl'),
 'MB 2N (G0+G1) count% (flow)': ('MB', 0.30, 'value', None, '2N (G0+G1) %'),
 'MB S count% (flow; BrdU Ts~3h)':('MB', 0.42, 'value', None, 'S %'),
 'MB G2+M duration ~2.5h (direct)':('MB', 0.55, 'value', None, 'G2+M duration (h)'),
 'MB HU S fold':                ('MB', 0.40, 'fold', None, 'HU S fold'),
 'MB HU G2 fold':               ('MB', 0.40, 'fold', None, 'HU G2 fold'),
 'EZH2 transcript S/G0 (1.8-2.5)':('MB', 0.45, 'fold', None, 'EZH2 mRNA S/G0'),
 'EZH2 protein G2/G0 (1.48)':   ('MB', 0.45, 'fold', None, 'EZH2 protein G2/G0'),
 'HU EZH2-in-S boost (1.31)':   ('MB', 0.35, 'fold', None, 'EZH2 in HU-S'),
 'MB cycles (>0 div)':          ('MB', None, 'binary', None, 'proliferates'),
 'MB+CDK4/6i arrest (0)':       ('MB', None, 'binary', None, 'arrests +CDK4/6i'),
 'MB+CDK4/6i+EZH2i no rescue':  ('MB', None, 'binary', None, '+EZH2i no rescue'),
 'CyclinD1 MB/GNP':             ('cross', 0.20, 'fold', None, 'CyclinD1 MB/GNP'),
 'CyclinD2 MB/GNP':             ('cross', 0.25, 'fold', None, 'CyclinD2 MB/GNP'),
 'MYCN MB/GNP':                 ('cross', 0.30, 'fold', None, 'MYCN MB/GNP'),
 'Gli1 MB/GNP':                 ('cross', 0.40, 'fold', None, 'Gli1 MB/GNP'),
 'EZH2 MB/GNP':                 ('cross', 0.35, 'fold', None, 'EZH2 MB/GNP'),
 'Skp2 MB/GNP':                 ('cross', 0.30, 'fold', None, 'Skp2 MB/GNP'),
}

def run_validate():
    tf = tempfile.NamedTemporaryFile('w', suffix='.json', delete=False); tf.close()
    env = dict(os.environ, VALIDATE_PARAMS=json.dumps(SPLIT), VALIDATE_RESULTS=tf.name,
               TWO_STEP_RB='1', H3K27_CHAIN='1', EZH2_CONC='0', DECOUPLE_COMMIT='1')
    subprocess.run([PY, VAL], cwd=ROOT, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=400)
    R = json.load(open(tf.name)); os.unlink(tf.name)
    return {c['name']: c for c in R['checks']}

def compute_extra():
    """Two targets that apply to BOTH cell types but the harness scores only once (JP 2026-07-27):
    EZH2 mRNA S/G0 for GNP (scRNA pools GNP Fig4A + MB Fig4B), Gli1 +HHi reduction for MB (98% down)."""
    os.environ.update(VALIDATE_PARAMS=json.dumps(SPLIT), DECOUPLE_COMMIT='1',
                      TWO_STEP_RB='1', H3K27_CHAIN='1', EZH2_CONC='0')
    sys.path.insert(0, os.path.join(ROOT, 'simulations'))
    import validate_v44 as V
    gnp = V.run(**V.CONDITIONS['GNP + SHH']); mb = V.run(**V.CONDITIONS['MB']); mbh = V.run(**V.CONDITIONS['MB + HHi'])
    thr = V.PRB_G0_THR
    return dict(ez_gnp=float(V.classify_grad(gnp, thr, 'EZH2m')['S']),
                gli_mb=float(1 - V.mean_settled(mbh, 'Gli1') / V.mean_settled(mb, 'Gli1')))

def mk_row(lab, act, tgt, tol):
    ratio = act / tgt
    return dict(lab=lab, ratio=ratio, tol=tol, passed=abs(ratio - 1) <= tol, kind='fold', act=act, tgt=tgt)

def insert_after(rows, after_lab, new):
    for i, r in enumerate(rows):
        if r['lab'] == after_lab:
            rows.insert(i + 1, new); return
    rows.append(new)

def prep(chk):
    """Return per-celltype ordered lists of rows: (label, ratio, tol, pass, kind, actual, target)."""
    groups = {'GNP': [], 'MB': [], 'cross': []}
    for name, (ct, tol, kind, tover, lab) in M.items():
        if name not in chk: continue
        c = chk[name]; act = c['actual']; tgt = tover if tover is not None else c['target']
        if kind == 'binary':
            passed = bool(c['pass']); ratio = None
        else:
            ratio = act / tgt if tgt else np.nan
            passed = abs(ratio - 1.0) <= tol
        groups[ct].append(dict(lab=lab, ratio=ratio, tol=tol, passed=passed, kind=kind, act=act, tgt=tgt))
    return groups

def draw_forest(ax, rows, title, accent):
    """Horizontal model/target ratio forest with tolerance band; binaries as check/cross glyphs."""
    folds = [r for r in rows if r['kind'] != 'binary']
    bins  = [r for r in rows if r['kind'] == 'binary']
    n = len(folds) + (1 if bins else 0) + len(bins)
    y = n - 1
    ax.axvline(1.0, color=FAINT, lw=1.0, zorder=1)
    yticks, ylabs = [], []
    for r in folds:
        ax.add_patch(plt.Rectangle((1 - r['tol'], y - 0.34), 2 * r['tol'], 0.68,
                     facecolor=BAND, edgecolor='none', zorder=0))
        col = C_PASS if r['passed'] else C_FAIL
        xr = min(max(r['ratio'], 0.02), 2.15)
        ax.plot([1, xr], [y, y], color=col, lw=1.4, alpha=0.55, zorder=2)
        ax.plot(xr, y, 'o', ms=8, color=col, mec='white', mew=1.1, zorder=3,
                clip_on=False if r['ratio'] > 2.15 else True)
        # value label
        av = r['act']
        txt = f"{av:.2f}" if abs(av) < 100 else f"{av:.0f}"
        ax.text(2.28, y, txt, fontsize=7.6, va='center', ha='left', color=INK if r['passed'] else C_FAIL,
                fontweight='normal' if r['passed'] else 'bold')
        if r['ratio'] > 2.15:
            ax.text(2.16, y, '›', fontsize=11, va='center', ha='left', color=col)
        yticks.append(y); ylabs.append(r['lab']); y -= 1
    if bins:
        y -= 0.2
        for r in bins:
            col = C_PASS if r['passed'] else C_FAIL
            ax.text(1.0, y, '✓' if r['passed'] else '✗', fontsize=12, va='center', ha='center',
                    color=col, fontweight='bold')
            yticks.append(y); ylabs.append(r['lab']); y -= 1
    ax.set_yticks(yticks); ax.set_yticklabels(ylabs, fontsize=8.2)
    ax.set_ylim(-0.6, n - 0.4); ax.set_xlim(0, 2.55)
    ax.set_xticks([0, 0.5, 1, 1.5, 2]); ax.set_xticklabels(['0', '', '1', '', '2'], fontsize=7.5)
    ax.tick_params(length=0)
    for s in ('top', 'right', 'left'): ax.spines[s].set_visible(False)
    ax.spines['bottom'].set_color(GRID)
    npass = sum(r['passed'] for r in rows)
    ax.set_title(f"{title}   ({npass}/{len(rows)})", fontsize=11, fontweight='bold', color=accent, loc='left', pad=8)

def main():
    chk = run_validate()
    passed = sum(1 for n in M if n in chk and (
        chk[n]['pass'] if M[n][3] is None or M[n][2] == 'binary'
        else abs(chk[n]['actual'] / M[n][3] - 1) <= M[n][1]))
    G = prep(chk)
    # add the two targets that apply to BOTH cell types (shown, not re-scored)
    ex = compute_extra()
    insert_after(G['GNP'], 'EZH2 Palbo drop', mk_row('EZH2 mRNA S/G0', ex['ez_gnp'], 2.0, 0.45))
    insert_after(G['MB'], 'MYCN +HHi/ctrl', mk_row('Gli1 +HHi reduction', ex['gli_mb'], 0.99, 0.05))

    fig = plt.figure(figsize=(15.5, 9.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[0.92, 1.15, 1.15], height_ratios=[1, 0.34],
                          hspace=0.42, wspace=0.42, left=0.055, right=0.965, top=0.86, bottom=0.075)
    axP = fig.add_subplot(gs[0, 0])
    axG = fig.add_subplot(gs[0, 1]); axM = fig.add_subplot(gs[0, 2])
    axC = fig.add_subplot(gs[1, 1:])
    axL = fig.add_subplot(gs[1, 0]); axL.axis('off')

    fig.suptitle('v44 cell-type parameterization — P7 GNP vs SHH-medulloblastoma', x=0.055, ha='left',
                 fontsize=15.5, fontweight='bold', color=INK, y=0.965)
    fig.text(0.055, 0.905, f"one engine, two cell types · only the biological identities differ · "
             f"the cycle-length split emerges from the MB CDKi · {passed}/{len(M)} harness-scored targets",
             fontsize=10, color=MUT, ha='left')

    # ---- Panel: parameter differences ----
    axP.axis('off')
    axP.set_title('cell-type parameters', fontsize=11, fontweight='bold', color=INK, loc='left', pad=8)
    rows = [('Hedgehog', None, None, None),
            ('Ptch1 copy #', '1.0', '0.1', 'Ptch loss'),
            ('MYCN expression', '1.0×', '2.86×', 'Gli-set'),
            ('CyclinD2 expr', '1.0×', '3.6×', 'Hh-buffered'),
            ('CDK inhibitors', None, None, None),
            ('p16  (INK4)', '0', '0.31', ''),
            ('p18  (INK4)', '0.46', '1.55', 'dominant'),
            ('birth p27', '0', '0.6', 'CIP/KIP')]
    y = 1.0; dy = 0.083
    axP.text(0.46, y, 'GNP', color=C_GNP, fontsize=10, fontweight='bold', ha='center', transform=axP.transAxes)
    axP.text(0.78, y, 'MB', color=C_MB, fontsize=10, fontweight='bold', ha='center', transform=axP.transAxes)
    y -= dy * 0.7
    for lab, g, m, note in rows:
        if g is None:  # section header
            y -= dy * 0.35
            axP.text(0.0, y, lab, fontsize=8.6, fontweight='bold', color=MUT, transform=axP.transAxes)
            y -= dy; continue
        axP.text(0.0, y, lab, fontsize=9, color=INK, transform=axP.transAxes)
        axP.text(0.46, y, g, fontsize=9.2, color=C_GNP, ha='center', transform=axP.transAxes)
        axP.text(0.78, y, m, fontsize=9.2, color=C_MB, ha='center', fontweight='bold', transform=axP.transAxes)
        if note: axP.text(0.99, y, note, fontsize=7.4, color=FAINT, ha='right', style='italic', transform=axP.transAxes)
        y -= dy
    # emergent cycle-length bar (the outcome)
    y -= dy * 0.5
    axP.text(0.0, y, 'core cell cycle', fontsize=8.6, fontweight='bold', color=MUT, transform=axP.transAxes)
    axP.text(0.99, y, 'emergent, via k_mu_cki', fontsize=7.4, color=FAINT, ha='right', style='italic', transform=axP.transAxes)
    y -= dy * 1.05
    bx0, bw = 0.0, 0.72
    for lab, val, col, maxv in [('GNP', GNP_CORE_H, C_GNP, MB_CORE_H), ('MB', MB_CORE_H, C_MB, MB_CORE_H)]:
        w = bw * val / (MB_CORE_H * 1.08)
        axP.add_patch(FancyBboxPatch((bx0, y - 0.028), w, 0.05, boxstyle='round,pad=0,rounding_size=0.008',
                      facecolor=col, edgecolor='none', transform=axP.transAxes, clip_on=False))
        axP.text(bx0 + w + 0.02, y, f"{lab}  {val:.1f} h", fontsize=9, color=col, va='center',
                 fontweight='bold', transform=axP.transAxes)
        y -= dy * 1.02
    y -= dy * 0.6
    axP.text(0.0, y, 'shared engine: mu 7.7e-4, k_mu_cki 0.31, decouple_commit,\nkPhRbCd 0.14, K_CdRb 0.40, w_Cd2 0.08 (two-cyclin)',
             fontsize=7.2, color=FAINT, transform=axP.transAxes, va='top', family='DejaVu Sans')

    # ---- validation forests ----
    draw_forest(axG, G['GNP'], 'P7 GNP', C_GNP)
    draw_forest(axM, G['MB'], 'SHH-medulloblastoma', C_MB)
    draw_forest(axC, G['cross'], 'MB ↔ GNP ratios', C_CROSS)
    axG.set_xlabel('model / target', fontsize=8.5, color=MUT)
    axM.set_xlabel('model / target', fontsize=8.5, color=MUT)
    axC.set_xlabel('model / target', fontsize=8.5, color=MUT)

    # legend
    axL.text(0.0, 0.9, 'validation', fontsize=8.6, fontweight='bold', color=MUT, transform=axL.transAxes)
    axL.plot(0.06, 0.62, 'o', ms=8, color=C_PASS, mec='white', mew=1.1, transform=axL.transAxes)
    axL.text(0.14, 0.62, 'within tolerance', fontsize=8.4, color=INK, va='center', transform=axL.transAxes)
    axL.plot(0.06, 0.40, 'o', ms=8, color=C_FAIL, mec='white', mew=1.1, transform=axL.transAxes)
    axL.text(0.14, 0.40, 'outside tolerance', fontsize=8.4, color=INK, va='center', transform=axL.transAxes)
    axL.add_patch(plt.Rectangle((0.02, 0.12), 0.09, 0.12, facecolor=BAND, edgecolor='none', transform=axL.transAxes))
    axL.text(0.14, 0.18, 'tolerance band', fontsize=8.4, color=INK, va='center', transform=axL.transAxes)

    out = os.path.join(ROOT, 'simulations', 'fig_celltype_validation')
    fig.savefig(out + '.pdf', bbox_inches='tight')
    fig.savefig(out + '.png', dpi=150, bbox_inches='tight')
    print('wrote', out + '.pdf / .png  |  passed', passed, '/', len(M))

if __name__ == '__main__':
    main()
