"""v44 -- mitogen sensitivity in Hh-activity (Gli1) units, with the CORRECT knobs:
  * GNP: the mitogen is SHH (ligand). Sweep SHH up from 0; GNP proliferation is SHH-driven.
    (To reach BELOW the SHH=0 basal-Smo-leak floor -- needed to see the no-feedback threshold --
     we extend the axis downward with a little HHi; SHH remains the physiological mitogen.)
  * MB: Hh is tonically high (Ptch1 loss, ligand-independent). The knob is vismodegib (HHi),
    which titrates the pathway DOWN. Sweep GDC up from 0 at the MB baseline.
Both knobs are tunable; the common x-axis is Gli1 transcript (the Hh-activity readout). Each panel
shows CyclinD1 and proliferation vs Gli1, with the EZH2->CyclinD1 feedback ON vs OFF.

Run:  ./venv/bin/python simulations/fig_v44_mitogen_sensitivity.py   (~12 min cold; caches)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
from simulations.validate_v44 import _new_rr, count_divisions, P16_MB, P18_MB, KSYP21_MB

XSEL = ["time", "MPF", "Cd", "Cd_mRNA", "Gli1_mRNA", "Gli1", "EZH2"]
T_END, N_PTS = 9000, 12000
WIN_H = (T_END - 4000) / 60.0
_rr = _new_rr()


def run(shh, hhi, ezh2i, ptch, mycn, p16=0.0, p18=None, ksyp21=None):
    for te_end, te_pts in ((T_END, N_PTS), (6000, 9000), (5000, 7000)):   # horizon retry for stiff (MB+vismo) conditions
        for atol in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
            _rr.reset()
            try: _rr.integrator.setValue("maximum_num_steps", 400000)
            except Exception: pass
            _rr['SHH'] = shh; _rr['Ptch1_copy_number'] = ptch; _rr['HHi'] = hhi
            _rr['EZH2i'] = ezh2i; _rr['MYCN_amplification'] = mycn
            _rr['p16'] = p16                          # MB CKI tones (INK4/CIP-KIP) raise the commitment threshold
            if p18 is not None: _rr['p18'] = p18
            if ksyp21 is not None: _rr['kSyP21'] = ksyp21
            _rr.integrator.setValue("absolute_tolerance", atol)
            try:
                r = _rr.simulate(0, te_end, te_pts, selections=XSEL)
                m = r['time'] >= 4000; nd, _, _ = count_divisions(r)
                return {s: float(np.mean(r[s][m])) for s in XSEL[1:]} | {'div': nd / ((te_end - 4000) / 60.0) * 168.0}
            except Exception:
                continue
    return {s: np.nan for s in XSEL[1:]} | {'div': np.nan}


# GNP: SHH is the mitogen (GDC=0); a few GDC>0 points extend the axis below the SHH=0 floor
GNP_PTS = [(0.5, 1.0), (0.5, 0.55), (0.5, 0.25)] + [(s, 0.0) for s in (0.0, 0.12, 0.25, 0.45, 0.7, 1.0, 1.4)]
# MB: vismodegib (GDC) titrates the tonic-high Hh down
MB_GDC = [0.0, 0.25, 0.45, 0.6, 0.72, 0.82, 0.9, 0.96, 1.0]


_MBCKI = dict(p16=P16_MB, p18=P18_MB, ksyp21=KSYP21_MB)


def sweep(points, ezh2i, ptch, mycn, is_mb):
    if is_mb:                                    # points = GDC (vismo) values at MB baseline, SHH=0.5; MB CKI tones on
        rows = [run(0.5, hhi, ezh2i, ptch, mycn, **_MBCKI) for hhi in points]
    else:                                        # points = (SHH, GDC) tuples in GNP context (GNP CKI defaults)
        rows = [run(shh, hhi, ezh2i, ptch, mycn) for (shh, hhi) in points]
    d = {k: np.array([r[k] for r in rows]) for k in rows[0]}
    o = np.argsort(d['Gli1_mRNA'])               # order by Hh activity (Gli1)
    return {k: v[o] for k, v in d.items()}


_CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_v44_mitogen_sensitivity_cache.npz')
if os.path.exists(_CACHE) and '--fresh' not in sys.argv:
    z = np.load(_CACHE); D = {k: z[k] for k in z.files}
    print('loaded cache; --fresh to recompute')
else:
    print('GNP +fb (SHH sweep) ...');  gfb = sweep(GNP_PTS, 0.0, 1.0, 1.0, False)
    print('GNP -fb ...');              gnf = sweep(GNP_PTS, 1.0, 1.0, 1.0, False)
    print('MB +fb (vismo sweep) ...'); mfb = sweep(MB_GDC, 0.0, 0.1, 2.8, True)
    print('MB -fb ...');               mnf = sweep(MB_GDC, 1.0, 0.1, 2.8, True)
    D = {}
    for tag, d in [('gfb', gfb), ('gnf', gnf), ('mfb', mfb), ('mnf', mnf)]:
        for k, v in d.items(): D[f'{tag}_{k}'] = v
    np.savez(_CACHE, **D)


def g(tag, k): return D[f'{tag}_{k}']
CFB, CNF = '#762A83', '#1B7837'
fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.6))


def cross(tag):                                  # Gli1 where it starts cycling (div>1)
    G, d = g(tag, 'Gli1_mRNA'), g(tag, 'div')
    a = np.where(d > 1)[0]
    return (G[0] if a[0] == 0 else 0.5 * (G[a[0] - 1] + G[a[0]])) if len(a) else np.nan


def panel(ax, tag_fb, tag_nf, key, ylabel, title, shade=True):
    ax.plot(g(tag_nf, 'Gli1_mRNA'), g(tag_nf, key), '--s', color=CNF, lw=2, ms=4, label='no feedback')
    ax.plot(g(tag_fb, 'Gli1_mRNA'), g(tag_fb, key), '-o', color=CFB, lw=2.4, ms=4, label='EZH2 feedback')
    tFB = cross(tag_fb)
    if shade and not np.isnan(tFB):
        ax.axvline(tFB, color=CFB, ls=':', lw=1.3, alpha=0.8)
        ax.axvspan(ax.get_xlim()[0], tFB, color='#F2EEF6', alpha=0.45)
    ax.set_xlabel('Hh activity = Gli1 transcript (a.u.)'); ax.set_ylabel(ylabel)
    ax.set_title(title, loc='left', fontweight='bold', fontsize=10.5)
    ax.legend(fontsize=8, loc='upper left')


panel(axes[0, 0], 'gfb', 'gnf', 'Cd', 'CyclinD1 protein (a.u.)', 'A  GNP CyclinD1 (mitogen = SHH $\\to$)')
panel(axes[0, 1], 'mfb', 'mnf', 'Cd', 'CyclinD1 protein (a.u.)', 'B  MB CyclinD1 (vismodegib $\\leftarrow$)')
panel(axes[1, 0], 'gfb', 'gnf', 'EZH2', 'EZH2 (a.u.)', 'C  GNP EZH2 dose-response (Fig 4H)')
panel(axes[1, 1], 'mfb', 'mnf', 'EZH2', 'EZH2 (a.u.)', 'D  MB EZH2 vs Hh activity')
axes[1, 0].text(0.97, 0.05, 'EZH2 now tracks the mitogen\nDOSE (saturating) -> rises across\nthe wide rShh range, not just at\nthe commitment threshold',
                transform=axes[1, 0].transAxes, fontsize=6.8, color='#B26500', va='bottom', ha='right')
axes[0, 0].text(0.97, 0.05, 'feedback represses & flattens\nCyclinD1 (lower mitogen sensitivity);\narrest region shaded',
                transform=axes[0, 0].transAxes, fontsize=6.8, color='#444', va='bottom', ha='right')

fig.suptitle('v44 mitogen sensitivity in Hh activity (Gli1): CyclinD1 (top) and the EZH2 dose-response (bottom), GNP via SHH / MB via vismodegib, +/- the EZH2$\\dashv$CyclinD1 feedback',
             fontsize=10, y=1.0)
fig.tight_layout()
fig.savefig('simulations/fig_v44_mitogen_sensitivity.png', dpi=160, bbox_inches='tight')
fig.savefig('simulations/fig_v44_mitogen_sensitivity.pdf', bbox_inches='tight')
print('wrote simulations/fig_v44_mitogen_sensitivity.{png,pdf}')
print(f"  GNP Gli1 range {g('gfb','Gli1_mRNA').min():.3f}-{g('gfb','Gli1_mRNA').max():.3f} (SHH); MB {g('mfb','Gli1_mRNA').min():.3f}-{g('mfb','Gli1_mRNA').max():.3f} (vismo)")
