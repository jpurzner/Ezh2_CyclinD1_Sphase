"""Memory phase diagrams: how H3K27me3 BUFFERING MEMORY changes across pairs of the mark-equation
parameters (JP 2026-07-21). Self-contained REDUCED equation (no full model needed -- this is the
summary equation itself):

    dM/dt = lam * W0 * (a0_prc2 + a_rw_prc2 * M^n_rw) * (1 - M)   -   delta_eff * M
    delta_eff = del_mk + k_jmjd3_gli*Gli1 + ln2/Tc + k_h33*tx       (the three-term effective decay)
    R(M) = f0 + (1-f0)/(1 + (M/K_prc2)^n_prc2)  ;  buffering B = 1 - R

Two FACES of memory:
  * PASSIVE persistence  tau = ln2 / delta_eff   (hours) -- how long a transient/withdrawal is remembered.
    Set purely by the DECAY params (del_mk, eraser k_jmjd3_gli*Gli1, dilution ln2/Tc, H3.3 k_h33).
  * LATCH memory = BISTABILITY: cooperative read-write (n_rw>=2) + low nucleation (a0) makes the mark
    self-sustaining -> HYSTERESIS in the mitogen(writing lam) response = the mark REMEMBERS its ON/OFF state.
    Metric = width in lam over which two stable states coexist (0 = monostable = no latch).

W0 is calibrated so the baked config sits at M* ~ 0.96 (the saturated mark). The baked operating point is
marked (*) on every panel. Run: ./venv/bin/python simulations/fig_memory_parameter_pairs.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); UPH = 62.8

# ---- baked v44 baseline (2026-07-21 four-mechanism config) ----
BASE = dict(a0=0.0028749753507882796, a_rw=0.0022657109194688194, n_rw=2.0,
            del_mk=0.00044000677598027846, k_jmjd3_gli=0.005, Gli=0.5,
            Tc=22.0, k_h33=0.003, tx=0.14,
            K_prc2=0.006856946735209077, n_prc2=4.015744557919901, f0=0.17482554666158287)
MSTAR = 0.96  # target baked ON-state mark used to calibrate the writing scale W0

def delta_eff(p):
    return p['del_mk'] + p['k_jmjd3_gli'] * p['Gli'] + np.log(2) / (p['Tc'] * UPH) + p['k_h33'] * p['tx']

# calibrate W0 so baseline lam=1 gives a fixed point at MSTAR
_d0 = delta_eff(BASE)
W0 = _d0 * MSTAR / ((BASE['a0'] + BASE['a_rw'] * MSTAR ** BASE['n_rw']) * (1 - MSTAR))
BASE['W0'] = W0

MG = np.concatenate([np.logspace(-6, -3, 120), np.linspace(1e-3, 1 - 1e-4, 320)])

def n_stable_fp(p, lam):
    d = delta_eff(p)
    f = lam * p['W0'] * (p['a0'] + p['a_rw'] * MG ** p['n_rw']) * (1 - MG) - d * MG
    idx = np.where(np.diff(np.sign(f)) != 0)[0]
    n = int(np.sum((f[idx] > 0) & (f[idx + 1] < 0)))     # interior + -> - crossings = stable
    if f[0] < 0:      # flow at M~0 points DOWN -> M=0 (fully-erased) is a stable OFF attractor
        n += 1
    return n

LAM = np.linspace(0.04, 6.0, 70)

def hysteresis_width(p):
    bist = [lam for lam in LAM if n_stable_fp(p, lam) >= 2]
    return (max(bist) - min(bist)) if bist else 0.0

def tau_hours(p):
    return np.log(2) / delta_eff(p) / UPH

def with_over(**ov):
    p = dict(BASE); p.update(ov); return p

# ---- panel specs: (metric, xkey, xlabel, xrange, ykey, ylabel, yrange) ----
LATCH = 'latch'; TAU = 'tau'
LOGK = {'a0'}   # parameters plotted on a log axis
PANELS = [
    (LATCH, 'a_rw', 'a_rw_prc2 (read-write)', (0.0, 0.02), 'a0', 'a0_prc2 (nucleation)', (1e-6, 1e-2), {}),
    (LATCH, 'n_rw', 'n_rw (cooperativity)', (1.0, 5.0), 'a0', 'a0_prc2 (nucleation)', (1e-6, 1e-2), {}),
    (LATCH, 'n_rw', 'n_rw (cooperativity)', (1.0, 5.0), 'a_rw', 'a_rw_prc2 (read-write)', (0.0, 0.02), {'a0': 1e-5}),
    (LATCH, 'Tc', 'Tc cell-cycle (h)', (8.0, 60.0), 'a0', 'a0_prc2 (nucleation)', (1e-6, 1e-2), {}),
    (TAU, 'k_jmjd3_gli', 'k_jmjd3_gli (eraser)', (0.0, 0.03), 'del_mk', 'del_mk (turnover)', (0.0, 0.004), {}),
    (TAU, 'k_jmjd3_gli', 'k_jmjd3_gli (eraser)', (0.0, 0.03), 'Tc', 'Tc cell-cycle (h)', (8.0, 60.0), {}),
    (TAU, 'k_h33', 'k_h33 (H3.3 turnover)', (0.0, 0.01), 'Tc', 'Tc cell-cycle (h)', (8.0, 60.0), {}),
    (TAU, 'k_h33', 'k_h33 (H3.3 turnover)', (0.0, 0.01), 'del_mk', 'del_mk (turnover)', (0.0, 0.004), {}),
]

def axis_vals(key, rng, n):
    return np.logspace(np.log10(rng[0]), np.log10(rng[1]), n) if key in LOGK else np.linspace(rng[0], rng[1], n)

NX = NY = 46
fig, axes = plt.subplots(2, 4, figsize=(21, 10.2))
fig.subplots_adjust(left=0.045, right=0.98, top=0.87, bottom=0.075, hspace=0.30, wspace=0.32)

for k, (metric, xk, xl, xr, yk, yl, yr, fixed) in enumerate(PANELS):
    ax = axes[k // 4, k % 4]
    xs = axis_vals(xk, xr, NX); ys = axis_vals(yk, yr, NY)
    Z = np.zeros((NY, NX))
    for iy, yv in enumerate(ys):
        for ix, xv in enumerate(xs):
            p = with_over(**{xk: xv, yk: yv, **fixed})
            Z[iy, ix] = hysteresis_width(p) if metric == LATCH else tau_hours(p)
    if metric == LATCH:
        cmap = 'viridis'; lab = 'bistable-latch memory\n(hysteresis width in mitogen λ)'
        vmax = max(Z.max(), 1e-6)
        im = ax.pcolormesh(xs, ys, Z, cmap=cmap, shading='auto', vmin=0, vmax=vmax)
        # boundary of the bistable (memory) region
        if Z.max() > 0:
            ax.contour(xs, ys, (Z > 1e-6).astype(float), levels=[0.5], colors='white', linewidths=1.8, linestyles='--')
        tag = 'LATCH'
    else:
        cmap = 'magma'; lab = 'passive memory  τ = ln2/δ_eff  (h)'
        im = ax.pcolormesh(xs, ys, Z, cmap=cmap, shading='auto')
        cs = ax.contour(xs, ys, Z, levels=[3, 6, 10, 15, 20], colors='white', linewidths=0.8, alpha=0.7)
        ax.clabel(cs, inline=True, fontsize=7, fmt='%g h')
        tag = 'τ'
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02); cb.ax.tick_params(labelsize=7)
    cb.set_label(lab, fontsize=7.5)
    # baked operating point
    bx, by = BASE[xk], BASE[yk]
    if xr[0] <= bx <= xr[1] and yr[0] <= by <= yr[1]:
        ax.plot(bx, by, marker='*', ms=20, mfc='#ff2d55', mec='white', mew=1.4, zorder=5)
    if xk in LOGK: ax.set_xscale('log')
    if yk in LOGK: ax.set_yscale('log')
    ax.set_xlabel(xl, fontsize=9); ax.set_ylabel(yl, fontsize=9)
    note = '  (at ' + ', '.join(f'{kk}={vv:g}' for kk, vv in fixed.items()) + ')' if fixed else ''
    if metric == TAU and 'k_jmjd3_gli' in (xk, yk):
        note += '  (Gli1=0.5, high-mitogen/MB)'
    ax.set_title(f'{chr(65+k)}  [{tag}]  {yl.split(" (")[0]}  vs  {xl.split(" (")[0]}{note}', fontsize=9.0, fontweight='bold')

fig.suptitle(
    'H3K27me3 buffering MEMORY across mark-equation parameter pairs   —   dM/dt = λ·W0·(a0 + a_rw·M^n_rw)·(1−M) − δ_eff·M,   '
    'δ_eff = del_mk + k_jmjd3_gli·Gli1 + ln2/Tc + k_h33·tx\n'
    'Row 1 = LATCH memory (bistable hysteresis; white dashed = onset of self-sustaining ON state — needs n_rw≥2 + low a0).   '
    'Row 2 = PASSIVE memory τ (eraser & fast cycling ⇒ shorter memory).   ★ = baked v44 operating point.',
    fontsize=11.5, fontweight='bold', y=0.975)
out = os.path.join(HERE, 'fig_memory_parameter_pairs.png')
fig.savefig(out, dpi=140); fig.savefig(out.replace('.png', '.pdf'))
print('->', out)
print(f'  W0 calibrated = {W0:.3f}   baseline δ_eff = {_d0:.5f}   baseline τ = {tau_hours(BASE):.2f} h')
print(f'  baked bistable? hysteresis width = {hysteresis_width(BASE):.3f}  (0 => monostable, no latch)')
# a0->low probe to show the latch opening
plow = with_over(a0=1e-5, a_rw=0.012, n_rw=3)
print(f'  low-a0 / high-a_rw / n_rw=3 probe: hysteresis width = {hysteresis_width(plow):.3f}  (latch opens)')
