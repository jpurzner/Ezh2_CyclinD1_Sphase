"""v44 -- mitogen sensitivity in Hh-activity (Gli1) units, with the CORRECT knobs:
  * GNP: the mitogen is SHH (ligand). Sweep SHH up from 0; GNP proliferation is SHH-driven.
    (To reach BELOW the SHH=0 basal-Smo-leak floor -- needed to see the no-feedback threshold --
     we extend the axis downward with a little GDC0449; SHH remains the physiological mitogen.)
  * MB: Hh is tonically high (Ptch1 loss, ligand-independent). The knob is vismodegib (GDC0449),
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


def run(shh, gdc, ezh2i, ptch, mycn, p16=0.0, p18=None, ksyp21=None):
    for te_end, te_pts in ((T_END, N_PTS), (6000, 9000), (5000, 7000)):   # horizon retry for stiff (MB+vismo) conditions
        for atol in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
            _rr.reset()
            try: _rr.integrator.setValue("maximum_num_steps", 400000)
            except Exception: pass
            _rr['SHH'] = shh; _rr['Ptch1_copy_number'] = ptch; _rr['GDC0449'] = gdc
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
        rows = [run(0.5, gdc, ezh2i, ptch, mycn, **_MBCKI) for gdc in points]
    else:                                        # points = (SHH, GDC) tuples in GNP context (GNP CKI defaults)
        rows = [run(shh, gdc, ezh2i, ptch, mycn) for (shh, gdc) in points]
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
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0))


def cross(tag):                                  # Gli1 where it starts cycling (div>1)
    G, d = g(tag, 'Gli1_mRNA'), g(tag, 'div')
    a = np.where(d > 1)[0]
    return (G[0] if a[0] == 0 else 0.5 * (G[a[0] - 1] + G[a[0]])) if len(a) else np.nan


def panel(ax, tag_fb, tag_nf, title, drive_note):
    # CyclinD1 protein vs Gli1 = the mitogen dose-response (primary)
    ax.plot(g(tag_nf, 'Gli1_mRNA'), g(tag_nf, 'Cd'), '--s', color=CNF, lw=2, ms=4, label='CyclinD1, no feedback')
    ax.plot(g(tag_fb, 'Gli1_mRNA'), g(tag_fb, 'Cd'), '-o', color=CFB, lw=2.4, ms=4, label='CyclinD1, EZH2 feedback')
    tFB, tNF = cross(tag_fb), cross(tag_nf)
    # shade the ARRESTED region (Gli1 below the cycling threshold) for each condition
    for t, c in [(tNF, CNF), (tFB, CFB)]:
        if not np.isnan(t):
            ax.axvline(t, color=c, ls=':', lw=1.4, alpha=0.8)
    lo = ax.get_xlim()[0]
    if not np.isnan(tFB):
        ax.axvspan(lo, tFB, color='#F2EEF6', alpha=0.45)
    ax.set_xlabel('Hh activity = Gli1 transcript (a.u.)'); ax.set_ylabel('CyclinD1 protein (a.u.)')
    ax.set_title(title, loc='left', fontweight='bold', fontsize=11)
    ax.legend(fontsize=8, loc='upper left')
    ax.text(0.97, 0.05, drive_note, transform=ax.transAxes, fontsize=7.0, color='#444', va='bottom', ha='right')
    ax.text(0.97, 0.93, f'cycling threshold (div>0):\nno-fb Gli1 {tNF:.3f}  |  +fb {tFB:.3f}\n(dotted; arrest = shaded)',
            transform=ax.transAxes, fontsize=6.6, color='#666', va='top', ha='right')


panel(axes[0], 'gfb', 'gnf', 'A  GNP — mitogen = SHH (drives Hh UP $\\to$)',
      'SHH $\\to$ raises Gli1 (the GNP mitogen);\nvismo only extends the axis below the\nSHH=0 basal floor')
panel(axes[1], 'mfb', 'mnf', 'B  MB — knob = vismodegib ($\\leftarrow$ Hh tonically high)',
      'MB Hh tonically HIGH (right edge);\nvismo (GDC0449) titrates it DOWN\ntoward the MYCN floor')

fig.suptitle('v44 mitogen sensitivity in Hh activity (Gli1): GNP driven UP by SHH (ligand), MB tuned DOWN by vismodegib -- with vs without the EZH2$\\dashv$CyclinD1 feedback (arrest region shaded)',
             fontsize=9.6, y=1.0)

fig.suptitle('v44 mitogen sensitivity in Hh activity (Gli1): GNP driven by SHH (ligand), MB tuned by vismodegib (suppresses tonic Hh) -- both with vs without the EZH2$\\dashv$CyclinD1 feedback',
             fontsize=9.8, y=1.0)
fig.tight_layout()
fig.savefig('simulations/fig_v44_mitogen_sensitivity.png', dpi=160, bbox_inches='tight')
fig.savefig('simulations/fig_v44_mitogen_sensitivity.pdf', bbox_inches='tight')
print('wrote simulations/fig_v44_mitogen_sensitivity.{png,pdf}')
print(f"  GNP Gli1 range {g('gfb','Gli1_mRNA').min():.3f}-{g('gfb','Gli1_mRNA').max():.3f} (SHH); MB {g('mfb','Gli1_mRNA').min():.3f}-{g('mfb','Gli1_mRNA').max():.3f} (vismo)")
