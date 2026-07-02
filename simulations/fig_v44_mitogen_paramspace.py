"""v44 -- 2D mitogen-sensitivity parameter planes on the committed Heldt-engine model.

How does the mitogen (SHH/Hh) threshold to proliferate move as a second parameter is tuned? Two planes,
each coloured by proliferation (divisions / 168 h) with the cycle/arrest boundary drawn, GNP context:
  (A) SHH x EZH2/H3K27me3-feedback strength (AUM f0_mk leaky floor; LOWER f0 = STRONGER repression).
      Stronger feedback pushes the Hh threshold to the right (needs more mitogen) -- the mitogen-sensitivity knob.
  (B) SHH x p16 (the MB-specific INK4 / CDK4/6 brake). More CKI also raises the Hh threshold.

The default (f0_mk=0.233) is marked; the boundary is the cycle/arrest line.

Caches the grids (simulations/*_cache.npz, gitignored; --fresh recomputes). ~20-30 min cold.
Run:  ./venv/bin/python simulations/fig_v44_mitogen_paramspace.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from simulations.validate_v44 import _new_rr, SEL, count_divisions

T_END, N_PTS = 8000, 12000           # 133 h window (enough to resolve cycle vs arrest)
WIN_H = (T_END - 3000) / 60.0
_rr = _new_rr()


def divrun(shh, f0=None, p16=0.0, ptch1_cn=1.0, mycn_amp=1.0):
    """divisions/168h for one (mitogen, parameter) point; mirrors validate_v44.run + sets f0_mk (AUM floor)."""
    for te_end, te_pts in ((T_END, N_PTS), (6000, 9000)):
        for atol in (1e-9, 1e-8, 1e-7, 1e-6):
            _rr.reset()
            try: _rr.integrator.setValue("maximum_num_steps", 300000)
            except Exception: pass
            _rr['SHH'] = shh; _rr['Ptch1_copy_number'] = ptch1_cn; _rr['HHi'] = 0.0
            _rr['EZH2i'] = 0.0; _rr['MYCN_amplification'] = mycn_amp; _rr['p16'] = p16
            if f0 is not None:
                _rr['f0_mk'] = f0
            _rr.integrator.setValue("absolute_tolerance", atol)
            try:
                r = _rr.simulate(0, te_end, te_pts, selections=SEL)
                n, _, _ = count_divisions(r)
                return n / ((te_end - 3000) / 60.0) * 168.0
            except Exception:
                continue
    return np.nan


def grid(shh_vals, p_vals, kind):
    Z = np.full((len(p_vals), len(shh_vals)), np.nan)
    for i, p in enumerate(p_vals):
        for j, s in enumerate(shh_vals):
            Z[i, j] = divrun(s, f0=p) if kind == 'K' else divrun(s, p16=p)
        print(f'  [{kind}] {("f0_mk" if kind=="K" else "p16")}={p:<7.3g} row done: '
              + ' '.join(f'{v:.0f}' for v in Z[i]))
    return Z


SHH = np.array([0.0, 0.1, 0.18, 0.28, 0.4, 0.55, 0.72, 0.9])
F0S = np.array([0.05, 0.15, 0.233, 0.4, 0.65, 0.9])      # AUM leaky floor; lower = stronger repression
P16 = np.array([0.0, 0.15, 0.4, 0.8, 1.5])
_CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fig_v44_mitogen_paramspace_cache.npz')
if os.path.exists(_CACHE) and '--fresh' not in sys.argv:
    z = np.load(_CACHE); ZK, ZP = z['ZK'], z['ZP']; SHH = z['SHH']; F0S = z['F0S']; P16 = z['P16']
    print('loaded cached grids; pass --fresh to recompute')
else:
    print('computing SHH x f0_mk grid ...'); ZK = grid(SHH, F0S, 'K')
    print('computing SHH x p16 grid ...'); ZP = grid(SHH, P16, 'p16')
    np.savez(_CACHE, ZK=ZK, ZP=ZP, SHH=SHH, F0S=F0S, P16=P16)

fig, axes = plt.subplots(1, 2, figsize=(13, 5.0))
vmax = np.nanmax([np.nanmax(ZK), np.nanmax(ZP)])


def plane(ax, X, Y, Z, ylabel, ydefault, title, ylog=False):
    pm = ax.pcolormesh(X, Y, Z, shading='nearest', cmap='viridis', vmin=0, vmax=vmax)
    cs = ax.contour(X, Y, np.nan_to_num(Z), levels=[1.0], colors='white', linewidths=2.2)
    ax.clabel(cs, fmt='cycle/arrest', fontsize=7.5)
    if ydefault is not None:
        ax.axhline(ydefault, color='#E41A1C', ls='--', lw=1.4)
        ax.text(X[-1], ydefault, ' default', color='#E41A1C', fontsize=7.5, va='bottom', ha='right')
    if ylog: ax.set_yscale('log')
    ax.set_xlabel('mitogen (SHH / Hh dose)'); ax.set_ylabel(ylabel)
    ax.set_title(title, loc='left', fontweight='bold', fontsize=11)
    fig.colorbar(pm, ax=ax, label='proliferation (divisions / 168 h)')


plane(axes[0], SHH, F0S, ZK, 'H3K27me3 leaky floor $f_0$  (lower = STRONGER feedback)', 0.233,
      'A  SHH x EZH2/H3K27me3-feedback strength', ylog=False)
axes[0].invert_yaxis()  # stronger feedback (lower f0) at top
axes[0].annotate('stronger feedback\n$\\to$ higher Hh threshold', xy=(0.28, 0.5), xytext=(0.45, 0.38),
                 fontsize=7.5, color='white', arrowprops=dict(arrowstyle='->', color='white', lw=1.2))
plane(axes[1], SHH, P16, ZP, 'p16 (INK4 / CDK4/6 brake)', 0.15, 'B  SHH x p16 (CKI brake)')
axes[1].annotate('more CKI\n$\\to$ higher Hh threshold', xy=(0.4, 0.8), xytext=(0.5, 1.15),
                 fontsize=7.5, color='white', arrowprops=dict(arrowstyle='->', color='white', lw=1.2))

fig.suptitle('v44 (Heldt-engine): 2D mitogen-sensitivity planes -- the Hh threshold to proliferate vs EZH2-feedback strength and the p16 CKI brake',
             fontsize=11, y=1.0)
fig.tight_layout()
fig.savefig('simulations/fig_v44_mitogen_paramspace.png', dpi=160, bbox_inches='tight')
fig.savefig('simulations/fig_v44_mitogen_paramspace.pdf', bbox_inches='tight')
print('wrote simulations/fig_v44_mitogen_paramspace.{png,pdf}')
