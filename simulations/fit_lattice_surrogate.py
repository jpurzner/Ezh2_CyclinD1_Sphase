"""Fit an ANTIMONY-EXPRESSIBLE two-box surrogate to the nucleosome lattice.

The lattice is ground truth but cannot live in Antimony (no arrays, no
per-nucleosome replication, and v44 cannot be restarted mid-arrest). This fits
the cheapest thing that CAN: two scalars, core and flank, with cross-coupled
read-write. Exactly the shape of the existing with_proximal_distal module, so a
good fit means the module survives the redesign with new parameters; a bad fit
means the geometry has to be carried some other way.

    dm_c/dt = k_eff*P_c*(1-m_c) - (del + k_kdm6b*Gli1*psi_c + exch_c)*m_c
    dm_f/dt = k_eff*P_f*(1-m_f) - (del + k_kdm6b*Gli1*psi_f + exch_f)*m_f
    P_c = a0 + a_rw*(w_cc*m_c + w_cf*m_f)^n     (w_cc = 1 - w_cf)
    P_f = a0 + a_rw*(w_ff*m_f + w_fc*m_c)^n     (w_ff = 1 - w_fc)
    at S phase:  m -> (1 - dil)*m

Targets: the lattice's (core, flank-plateau) across the Gli1 dose-response and
the cell-cycle-length sweep.
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src', 'lattice'))
import numpy as np
from dataclasses import replace
from scipy.optimize import least_squares
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from promoter_lattice import Params, Inputs, ensemble

OUT = os.path.dirname(os.path.abspath(__file__))
PAR = Params(a_rw_prc2=3.5)
INP = Inputs(Gli1=0.3, activity=0.4, tx=0.4)
EZI_GRID = [0.0, 0.2, 0.4, 0.6, 0.8, 0.95]
GLI_GRID = [0.0, 0.3, 0.6, 1.0, 1.5, 2.2, 3.0]
CYC_GRID = [None, 12.0, 18.0, 24.0, 36.0, 48.0, 96.0]


# ----------------------------------------------------------------- surrogate
def surrogate_ss(theta, Gli1, cycle_h, hours=600.0, dt=0.02, ezh2i=0.0):
    k_eff, a0, a_rw, n, w_cf, w_fc, psi_c, psi_f, ex_c, ex_f = theta
    dl, kk, dil = PAR.del_mk, PAR.k_kdm6b, PAR.dil_frac
    m_c, m_f = 0.8, 0.8
    nS = int(hours / dt)
    nxt = cycle_h if cycle_h else np.inf
    t = 0.0
    acc_c, acc_f, nacc = 0.0, 0.0, 0
    for i in range(nS):
        Pc = a0 + a_rw * max(0.0, (1 - w_cf) * m_c + w_cf * m_f) ** n
        Pf = a0 + a_rw * max(0.0, (1 - w_fc) * m_f + w_fc * m_c) ** n
        dc = (1 - ezh2i) * k_eff * Pc * (1 - m_c) - (dl + kk * Gli1 * psi_c + ex_c) * m_c
        df = (1 - ezh2i) * k_eff * Pf * (1 - m_f) - (dl + kk * Gli1 * psi_f + ex_f) * m_f
        m_c = min(1.0, max(0.0, m_c + dc * dt))
        m_f = min(1.0, max(0.0, m_f + df * dt))
        t += dt
        if t >= nxt:
            m_c *= (1 - dil); m_f *= (1 - dil); nxt += cycle_h
        if i > nS // 2:
            acc_c += m_c; acc_f += m_f; nacc += 1
    return acc_c / nacc, acc_f / nacc


# ----------------------------------------------------------------- lattice targets
def lattice_targets(n_cells=20, cache=os.path.join(OUT, '.lattice_targets.json')):
    if os.path.exists(cache):
        return json.load(open(cache))
    tg = dict(gli=[], cyc=[])
    for g in GLI_GRID:
        e = ensemble(n_cells=n_cells, hours=500, dt=0.05, par=PAR,
                     inp=replace(INP, Gli1=g), init_me3=0.8)
        x, pf, ic = e['x'], e['profile'], e['in_core']
        tg['gli'].append(dict(Gli1=g, core=float(pf[ic].mean()),
                              plateau=float(np.percentile(pf[x > 2000], 75))))
        print(f'  lattice Gli1={g:4.1f} core={tg["gli"][-1]["core"]:.3f} '
              f'plateau={tg["gli"][-1]["plateau"]:.3f}')
    for ch in CYC_GRID:
        e = ensemble(n_cells=n_cells, hours=600, dt=0.05, par=PAR, inp=INP,
                     init_me3=0.8, cycle_h=ch)
        x, pf, ic = e['x'], e['profile'], e['in_core']
        tg['cyc'].append(dict(cycle_h=ch, core=float(pf[ic].mean()),
                              plateau=float(np.percentile(pf[x > 2000], 75))))
        print(f'  lattice cycle={str(ch):>5} core={tg["cyc"][-1]["core"]:.3f} '
              f'plateau={tg["cyc"][-1]["plateau"]:.3f}')
    tg['ezi'] = []
    for z in EZI_GRID:
        e = ensemble(n_cells=n_cells, hours=500, dt=0.05, par=PAR,
                     inp=replace(INP, EZH2i=z), init_me3=0.8)
        x, pf, ic = e['x'], e['profile'], e['in_core']
        tg['ezi'].append(dict(EZH2i=z, core=float(pf[ic].mean()),
                              plateau=float(np.percentile(pf[x > 2000], 75))))
        print(f'  lattice EZH2i={z:4.2f} core={tg["ezi"][-1]["core"]:.3f} '
              f'plateau={tg["ezi"][-1]["plateau"]:.3f}')
    json.dump(tg, open(cache, 'w'), indent=2)
    return tg


# ----------------------------------------------------------------- fit
def residuals(theta, tg):
    r = []
    for row in tg['gli']:
        c, f = surrogate_ss(theta, row['Gli1'], None)
        r += [c - row['core'], f - row['plateau']]
    for row in tg['cyc']:
        c, f = surrogate_ss(theta, INP.Gli1, row['cycle_h'])
        r += [c - row['core'], f - row['plateau']]
    return np.array(r)


def residuals_axis(theta, tg, axis):
    r = []
    if axis in ('gli', 'both'):
        for row in tg['gli']:
            c, f = surrogate_ss(theta, row['Gli1'], None)
            r += [c - row['core'], f - row['plateau']]
    if axis in ('cyc', 'both'):
        for row in tg['cyc']:
            c, f = surrogate_ss(theta, INP.Gli1, row['cycle_h'])
            r += [c - row['core'], f - row['plateau']]
    return np.array(r)


def surrogate_fixed_points(theta, Gli1, grid=1201):
    """Is the fitted surrogate bistable? The lattice is NOT (see
    lattice_switching). If the surrogate is, it has invented a behaviour."""
    k_eff, a0, a_rw, n, w_cf, w_fc, psi_c, psi_f, ex_c, ex_f = theta
    dl, kk = PAR.del_mk, PAR.k_kdm6b
    ms = np.linspace(0, 1, grid)
    # flank self-consistency with m_c slaved (w_fc small): treat flank alone
    g = np.array([k_eff * (a0 + a_rw * ((1 - w_fc) * m) ** n) * (1 - m)
                  - (dl + kk * Gli1 * psi_f + ex_f) * m for m in ms])
    roots = []
    for i in range(len(ms) - 1):
        if g[i] * g[i + 1] < 0:
            roots.append(float(0.5 * (ms[i] + ms[i + 1])))
    return roots


def cross_validate(tg, x0, lo, hi):
    print('\n=== OUT-OF-SAMPLE: fit one axis, predict the other ===')
    out = {}
    for fit_axis, pred_axis in [('gli', 'cyc'), ('cyc', 'gli')]:
        res = least_squares(residuals_axis, x0, bounds=(lo, hi), args=(tg, fit_axis),
                            xtol=1e-9, ftol=1e-9)
        rp = residuals_axis(res.x, tg, pred_axis)
        obs = np.array([v for row in tg[pred_axis]
                        for v in ((row['core'], row['plateau']))])
        rmse_in = float(np.sqrt((res.fun ** 2).mean()))
        rmse_out = float(np.sqrt((rp ** 2).mean()))
        r2_out = float(1 - (rp ** 2).sum() / ((obs - obs.mean()) ** 2).sum())
        print(f'  fit {fit_axis:>3} -> predict {pred_axis:>3}:  '
              f'in-sample RMSE {rmse_in:.4f} | OUT-OF-SAMPLE RMSE {rmse_out:.4f}  R2 {r2_out:+.3f}')
        out[f'{fit_axis}->{pred_axis}'] = dict(rmse_in=rmse_in, rmse_out=rmse_out, r2_out=r2_out,
                                               theta=res.x.tolist())
    return out


def main():
    free_n = '--free-n' in sys.argv
    print('=== lattice targets ===')
    tg = lattice_targets()
    #        k_eff  a0    a_rw   n    w_cf  w_fc  psi_c psi_f ex_c  ex_f
    x0 = np.array([0.25, 0.10, 3.50, 2.0, 0.40, 0.15, 1.00, 0.25, 0.06, 0.06])
    lo = np.array([0.02, 0.00, 0.20, 0.5, 0.00, 0.00, 0.10, 0.00, 0.00, 0.00])
    hi = np.array([3.00, 1.00, 20.0, 4.0, 1.00, 1.00, 3.00, 2.00, 0.60, 0.60])
    if not free_n:
        lo[3] = hi[3] = x0[3] = 2.0          # pin n to the lattice's n_rw
        lo[3], hi[3] = 1.999, 2.001
    print('\n=== fitting ===')
    t0 = time.time()
    res = least_squares(residuals, x0, bounds=(lo, hi), args=(tg,),
                        xtol=1e-10, ftol=1e-10, verbose=0)
    names = ['k_eff', 'a0', 'a_rw', 'n_rw', 'w_cf', 'w_fc', 'psi_c', 'psi_f', 'ex_c', 'ex_f']
    print(f'  converged in {time.time()-t0:.0f}s, cost={res.cost:.5f}')
    for nm, v in zip(names, res.x):
        print(f'    {nm:7s} = {v:.4f}')
    r = res.fun
    obs = np.array([v for row in tg['gli'] for v in (row['core'], row['plateau'])]
                   + [v for row in tg['cyc'] for v in (row['core'], row['plateau'])])
    rmse = float(np.sqrt((r ** 2).mean()))
    r2 = float(1 - (r ** 2).sum() / ((obs - obs.mean()) ** 2).sum())
    print(f'  RMSE = {rmse:.4f}   R2 = {r2:.4f}   (observable range {obs.min():.3f}-{obs.max():.3f})')

    # per-axis breakdown
    print('\n=== fidelity ===')
    print(f'  {"Gli1":>5} {"core lat":>9} {"core sur":>9} {"plat lat":>9} {"plat sur":>9}')
    for row in tg['gli']:
        c, f = surrogate_ss(res.x, row['Gli1'], None)
        print(f'  {row["Gli1"]:5.1f} {row["core"]:9.3f} {c:9.3f} {row["plateau"]:9.3f} {f:9.3f}')
    print(f'  {"cycle":>5} {"core lat":>9} {"core sur":>9} {"plat lat":>9} {"plat sur":>9}')
    for row in tg['cyc']:
        c, f = surrogate_ss(res.x, INP.Gli1, row['cycle_h'])
        lab = 'none' if row['cycle_h'] is None else f'{row["cycle_h"]:g}'
        print(f'  {lab:>5} {row["core"]:9.3f} {c:9.3f} {row["plateau"]:9.3f} {f:9.3f}')

    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    g = [r_['Gli1'] for r_ in tg['gli']]
    sc = [surrogate_ss(res.x, gg, None) for gg in g]
    ax[0].plot(g, [r_['core'] for r_ in tg['gli']], 'o-', color='#b0402f', label='lattice core')
    ax[0].plot(g, [s[0] for s in sc], 'o--', color='#b0402f', mfc='none', label='surrogate core')
    ax[0].plot(g, [r_['plateau'] for r_ in tg['gli']], 's-', color='#2f6fb0', label='lattice flank')
    ax[0].plot(g, [s[1] for s in sc], 's--', color='#2f6fb0', mfc='none', label='surrogate flank')
    ax[0].set_xlabel('Gli1 (KDM6B drive)'); ax[0].set_ylabel('me3'); ax[0].legend(fontsize=7.5)
    ax[0].set_title('Gli1 axis', fontweight='bold', fontsize=10)
    cyc = [200 if r_['cycle_h'] is None else r_['cycle_h'] for r_ in tg['cyc']]
    scy = [surrogate_ss(res.x, INP.Gli1, r_['cycle_h']) for r_ in tg['cyc']]
    ax[1].plot(cyc, [r_['core'] for r_ in tg['cyc']], 'o-', color='#b0402f', label='lattice core')
    ax[1].plot(cyc, [s[0] for s in scy], 'o--', color='#b0402f', mfc='none', label='surrogate core')
    ax[1].plot(cyc, [r_['plateau'] for r_ in tg['cyc']], 's-', color='#2f6fb0', label='lattice flank')
    ax[1].plot(cyc, [s[1] for s in scy], 's--', color='#2f6fb0', mfc='none', label='surrogate flank')
    ax[1].set_xlabel('cell-cycle length (h; 200 = none)'); ax[1].set_ylabel('me3'); ax[1].legend(fontsize=7.5)
    ax[1].set_title('cell-cycle axis', fontweight='bold', fontsize=10)
    fig.suptitle(f'Two-box Antimony surrogate vs nucleosome lattice  '
                 f'(R$^2$={r2:.3f}, RMSE={rmse:.3f}{", n free" if free_n else ", n pinned to 2"})',
                 fontweight='bold')
    fig.tight_layout()
    tag = '_freen' if free_n else ''
    fig.savefig(os.path.join(OUT, f'lattice_surrogate_fit{tag}.png'), dpi=140, bbox_inches='tight')
    json.dump(dict(theta=dict(zip(names, res.x.tolist())), rmse=rmse, r2=r2, free_n=free_n),
              open(os.path.join(OUT, f'lattice_surrogate_fit{tag}.json'), 'w'), indent=2)
    print(f'\n  wrote lattice_surrogate_fit{tag}.png/.json')

    cv = cross_validate(tg, x0, lo, hi)

    print('\n=== HELD-OUT PERTURBATION: fit on Gli1+cycle, predict EZH2i ===')
    rr_ = []
    for row in tg.get('ezi', []):
        c, f = surrogate_ss(res.x, INP.Gli1, None, ezh2i=row['EZH2i'])
        rr_ += [c - row['core'], f - row['plateau']]
        print(f'  EZH2i={row["EZH2i"]:4.2f}  core lat={row["core"]:.3f} sur={c:.3f} '
              f'| plateau lat={row["plateau"]:.3f} sur={f:.3f}')
    if rr_:
        rr_ = np.array(rr_)
        obs_z = np.array([v for row in tg['ezi'] for v in (row['core'], row['plateau'])])
        rmse_z = float(np.sqrt((rr_ ** 2).mean()))
        r2_z = float(1 - (rr_ ** 2).sum() / ((obs_z - obs_z.mean()) ** 2).sum())
        print(f'  --> HELD-OUT RMSE {rmse_z:.4f}   R2 {r2_z:+.3f}')
        cv['heldout_ezh2i'] = dict(rmse=rmse_z, r2=r2_z)

    print('\n=== does the surrogate invent bistability the lattice does not have? ===')
    for g_ in [0.0, 0.3, 1.0]:
        rts = surrogate_fixed_points(res.x, g_)
        print(f'  Gli1={g_:4.1f}  flank fixed points: {[round(r_,3) for r_ in rts]}'
              f'   -> {"BISTABLE" if len(rts) >= 3 else "monostable"}')

    json.dump(dict(theta=dict(zip(names, res.x.tolist())), rmse=rmse, r2=r2,
                   free_n=free_n, cross_validation=cv),
              open(os.path.join(OUT, f'lattice_surrogate_fit{tag}.json'), 'w'), indent=2)


if __name__ == '__main__':
    main()
