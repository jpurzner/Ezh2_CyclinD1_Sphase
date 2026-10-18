"""CCND1 transcript half-life vs the model's reservoir assumption (JP 2026-07-24). Uses the RNAdecayCafe atlas
(Vock; 16 nucleotide-recoding RNA-seq datasets, uniform EZbakR pipeline; data/AvgKdegs_genes_v1.1.csv). Question:
is the model's SLOW CyclinD1 mRNA reservoir (t½~11h) consistent with the measured CCND1 half-life? 3 panels:
 A  half-life ladder: CCND1 vs labile (MYC/FOS/MYCN) and stable (ACTB/GAPDH) anchors + D-cyclins + Hh genes;
    model reservoir (11h) marked — it sits in the STABLE tier, far from real CCND1.
 B  model Ccnd1 drop after Hh withdrawal at the baked 11h vs data-realistic 1.5-2.5h — the model coasts ~12h,
    the data implies a ~2h drop.
 C  CCND1 per-cell-line spread (donorm vs raw), MDA-MB amplification outlier, 3'end-only vs total flag.
Run: ./venv/bin/python simulations/fig_ccnd1_halflife.py"""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1'); os.environ.setdefault('EZH2_CONC', '0')
import tellurium as te
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from src.build_model_v44_heldt import build_model_v44
HERE = os.path.dirname(os.path.abspath(__file__)); UPH = 62.8; LN2 = np.log(2)
df = pd.read_csv(os.path.join(os.path.dirname(HERE), 'data', 'AvgKdegs_genes_v1.1.csv'))
GNP = dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464)
MODEL_THALF = LN2 / (te.loada(build_model_v44())['k_Cd_mRNA_deg'] * UPH)   # baked k_Cd_mRNA_deg -> t½ in h (now the MEASURED 3.4h, swapped 2026-10-16)

def med(g, col='avg_donorm_halflife'):
    s = df[df.feature_ID == g][col].replace(0, np.nan).dropna(); return (np.median(s), s.values)
GENES = [('MYC', 'labile'), ('FOS', 'labile'), ('MYCN', 'labile'), ('CCND1', 'ccnd1'), ('GLI1', 'hh'),
         ('CCND3', 'dcyc'), ('EZH2', 'hh'), ('CCND2', 'dcyc'), ('ACTB', 'stable'), ('GAPDH', 'stable')]
COL = {'labile': '#c0392b', 'ccnd1': '#0e6251', 'dcyc': '#16a085', 'hh': '#8e44ad', 'stable': '#2c3e50'}

# model Ccnd1 decline after Hh withdrawal at several half-lives
def wd_traj(thalf_h, Twarm=130, Tpost=44):
    kmu = LN2 / (thalf_h * UPH); scale = kmu / 0.001
    ov = dict(k_Cd_mRNA_deg=kmu, k_Cd_tx_basal=0.00028075924684120905 * scale,
              k_Cd_tx_Gli_max=0.0940819151141191 * scale, k_Cd_tx_MYCN=0.020257706030264946 * scale)
    for atol in (1e-8, 1e-7, 1e-6, 1e-5):                 # stiffness retry: some half-lives need a looser tol
        try:
            rr = te.loada(build_model_v44(params=ov)); rr.reset()
            for k, v in GNP.items(): rr[k] = v
            rr['SHH'] = 0.5; rr.integrator.setValue('absolute_tolerance', atol); rr.integrator.setValue('maximum_num_steps', 12000000)
            rr.simulate(0, Twarm * UPH, 300); rr['HHi'] = 1
            r = rr.simulate(Twarm * UPH, (Twarm + Tpost) * UPH, 2500, selections=['time', 'Cd_mRNA'])
            return r['time'] / UPH - Twarm, r['Cd_mRNA'] / r['Cd_mRNA'][0]
        except Exception:
            continue
    th = np.linspace(0, Tpost, 2500)                       # all tolerances failed -> flat NaN (skip this half-life)
    return th, np.full_like(th, np.nan)

fig, ax = plt.subplots(1, 3, figsize=(21.0, 6.6))
fig.subplots_adjust(left=0.045, right=0.99, top=0.82, bottom=0.13, wspace=0.28)

# A: half-life ladder (log-x)
a = ax[0]
gm_all = df[df.avg_donorm_halflife > 0].groupby('feature_ID').avg_donorm_halflife.median()
a.axvspan(gm_all.quantile(.1), gm_all.quantile(.9), color='#ecf0f1', alpha=0.9, zorder=0)
a.axvline(gm_all.median(), color='#bbb', ls='-', lw=1, zorder=1); a.text(gm_all.median(), 10.4, 'genome median', fontsize=7, color='#999', ha='center')
for i, (g, cat) in enumerate(GENES):
    dn, dnv = med(g, 'avg_donorm_halflife'); rw, _ = med(g, 'avg_halflife')
    y = len(GENES) - i
    a.plot(dnv, [y] * len(dnv), 'o', color=COL[cat], ms=4, alpha=0.35, zorder=2)
    a.plot(dn, y, 'o', color=COL[cat], ms=11, mec='k', mew=0.6, zorder=3)
    a.plot(rw, y, 'D', color=COL[cat], ms=6, mfc='none', mew=1.4, zorder=3)
    a.text(0.42, y, g, fontsize=9.5, ha='right', va='center', fontweight=('bold' if g == 'CCND1' else 'normal'), color=COL[cat])
a.axvline(MODEL_THALF, color='#e67e22', lw=2.4, ls='--', zorder=4)
a.text(MODEL_THALF, len(GENES) + 0.6, 'model mRNA\n(baked %.1f h, measured)' % MODEL_THALF, fontsize=8.5, color='#e67e22', ha='center', fontweight='bold')
a.set_xscale('log'); a.set_xlim(0.4, 40); a.set_ylim(0.3, len(GENES) + 1.4)
a.set_xlabel('mRNA half-life (h, log)  ·  ● donorm  ◇ raw  (small = per cell line)'); a.set_yticks([])
a.set_title('A · CCND1 is a LABILE transcript (t½≈1.5 h donorm),\nfar below the model\'s 11 h reservoir', fontsize=11, fontweight='bold')

# B: model Ccnd1 drop after Hh withdrawal
b = ax[1]
for th, col, lab in ((MODEL_THALF, '#e67e22', 'model as baked (11 h)'), (2.5, '#16a085', 'data raw (2.5 h)'), (1.5, '#0e6251', 'data donorm (1.5 h)')):
    t, mr = wd_traj(th); w = t >= -2
    b.plot(t[w], mr[w], color=col, lw=2.6, label=lab)
    i50 = np.argmax(mr < 0.5); t50 = t[i50] if mr[i50] < 0.5 else np.nan
    if t50 == t50: b.plot(t50, 0.5, 'v', color=col, ms=9, mec='k', mew=0.5)
b.axhline(0.5, color='#bbb', ls=':', lw=1); b.axvline(0, color='#333', ls=':', lw=1.2)
b.text(1, 0.94, 'Hh OFF', fontsize=8.5, color='#333')
b.text(13, 0.55, 'model coasts\n~12 h to 50%', fontsize=8.3, color='#e67e22', fontweight='bold')
b.text(3.2, 0.30, 'data → ~2 h to 50%', fontsize=8.3, color='#0e6251', fontweight='bold')
b.set_xlabel('time after Hh withdrawal (h)'); b.set_ylabel('Ccnd1 transcript (norm. to start)'); b.set_xlim(-2, 40); b.set_ylim(0, 1.05)
b.set_title('B · "How long until Ccnd1 drops?" — the baked reservoir\npredicts ~12 h; the measured half-life predicts ~2 h', fontsize=11, fontweight='bold')
b.legend(fontsize=8.6, loc='upper right'); b.grid(alpha=0.2)

# C: CCND1 per cell line
c = ax[2]
cc = df[df.feature_ID == 'CCND1'].copy().sort_values('avg_donorm_halflife')
y = np.arange(len(cc))
for k, (_, r) in enumerate(cc.iterrows()):
    src = "3'end" if r.avg_RPKM_total == 0 else ("total" if r.avg_RPKM_3pend == 0 else "both")
    c.plot([r.avg_donorm_halflife, r.avg_halflife], [k, k], color='#bbb', lw=1, zorder=1)
    c.plot(r.avg_donorm_halflife, k, 'o', color='#0e6251', ms=8, zorder=3)
    c.plot(r.avg_halflife, k, 'D', color='#0e6251', ms=6, mfc='none', mew=1.4, zorder=3)
    c.text(0.42, k, '%s (%s)' % (r.cell_line, src), fontsize=8, ha='right', va='center', color='#444')
c.axvline(np.median(cc.avg_donorm_halflife), color='#0e6251', ls='--', lw=1.4)
c.text(np.median(cc.avg_donorm_halflife), len(cc) - 0.3, 'median %.1f h' % np.median(cc.avg_donorm_halflife), fontsize=8, color='#0e6251', ha='center')
c.set_xscale('log'); c.set_xlim(0.4, 25); c.set_yticks([])
c.set_xlabel('CCND1 half-life (h, log)  ·  ● donorm  ◇ raw'); c.set_ylabel('')
c.set_title('C · CCND1 per cell line (donorm ● vs raw ◇)\nMDA-MB raw 18 h = its CCND1 amplification (→3.9 h donorm)', fontsize=11, fontweight='bold')

fig.suptitle('CCND1 transcript half-life (RNAdecayCafe atlas) vs the model\'s reservoir:  measured CCND1 is LABILE (~1.5–2.5 h, 28th pct — near MYC/FOS, far below ACTB/GAPDH), '
             'so the baked 11 h reservoir is ~5–7× too slow — the ~11 h memory is likely the committed cell-cycle state, not the transcript',
             fontsize=11.4, fontweight='bold', y=0.965)
out = os.path.join(HERE, 'fig_ccnd1_halflife.png')
fig.savefig(out, dpi=140); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
print('CCND1 median donorm t½ = %.2f h (raw %.2f h); model baked = %.1f h' % (med('CCND1')[0], med('CCND1', 'avg_halflife')[0], MODEL_THALF))
