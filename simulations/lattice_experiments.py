"""Experiment driver for the nucleosome-lattice promoter model.

  python simulations/lattice_experiments.py <exp> [...]

  bifurcation   mean-field fixed points vs read-write strength
  branches      lattice vs mean-field in a homogeneous flank (finite-size switching)
  profile       full promoter: depleted core, bounded domain, boundary position
  versions      model variants: local vs 3D-contact kernel, chain vs collapsed, no read-write
  gli           Gli1 (KDM6B) dose-response = the GNP -> MB axis
  cycle         replication dilution: steady-state me3 vs cell-cycle length
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src', 'lattice'))
import numpy as np
from dataclasses import replace
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from promoter_lattice import (PromoterLattice, Geometry, Params, Inputs, ensemble,
                              U, AC, ME1, ME2, ME3)
import mean_field as mf

OUT = os.path.dirname(os.path.abspath(__file__))
# working point: inside the bistable window on the high branch
PAR = Params(a_rw_prc2=3.5)
INP = Inputs(Gli1=0.3, activity=0.4, tx=0.4)


def _save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, name + '.png'), dpi=140, bbox_inches='tight')
    plt.close(fig)
    print('  wrote', name + '.png')


def boundary_metrics(x, prof, in_core):
    """Core depletion, plateau height, and the downstream boundary position/width."""
    core = float(prof[in_core].mean())
    dn = x > 0
    plateau = float(np.percentile(prof[dn & (x > 2000)], 75)) if (dn & (x > 2000)).any() else np.nan
    res = dict(core_me3=core, plateau=plateau, depletion=plateau - core)
    if plateau > 0.08:
        half = 0.5 * (plateau + core)
        xs, ps = x[dn], prof[dn]
        above = np.where(ps > half)[0]
        res['edge_bp'] = float(xs[above[0]]) if len(above) else np.nan
        # 20-80% rise width
        lo_t, hi_t = core + 0.2 * (plateau - core), core + 0.8 * (plateau - core)
        i_lo = np.where(ps > lo_t)[0]; i_hi = np.where(ps > hi_t)[0]
        res['width_bp'] = (float(xs[i_hi[0]] - xs[i_lo[0]])
                           if len(i_lo) and len(i_hi) else np.nan)
    else:
        res['edge_bp'] = np.nan; res['width_bp'] = np.nan
    return res


# ----------------------------------------------------------------- experiments
def exp_bifurcation():
    grid = np.arange(0.5, 8.01, 0.1)
    b = mf.bifurcation(PAR, INP, grid)
    bw = grid[b['bistable']]
    print(f'  bistable window a_rw = {bw.min():.2f} .. {bw.max():.2f}'
          if len(bw) else '  no bistable window')
    print(f'  kappa_param at saddle-node = {b["kappa"][b["bistable"]][0]:.2f}')
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.plot(grid, b['high'], color='#b0402f', lw=2, label='stable (high / marked)')
    ax.plot(grid, b['low'], color='#2f6fb0', lw=2, label='stable (low / unmarked)')
    ax.plot(grid, b['unstable'], color='#888', lw=1.6, ls='--', label='unstable (separatrix)')
    if len(bw):
        ax.axvspan(bw.min(), bw.max(), color='#f2e0c9', alpha=.5, zorder=0, label='bistable')
    ax.axvline(PAR.a_rw_prc2, color='k', ls=':', lw=.9)
    ax.set_xlabel('$a_{rw}$ (EED read-write strength)'); ax.set_ylabel('flank me3 fixed point')
    ax.set_title('Mean-field: the flank is bistable', fontweight='bold')
    ax.legend(fontsize=8)
    _save(fig, 'lattice_bifurcation')
    json.dump({k: np.asarray(v).tolist() for k, v in b.items()},
              open(os.path.join(OUT, 'lattice_bifurcation.json'), 'w'))


def exp_branches():
    """Homogeneous flank: does the lattice sit on the mean-field branches?"""
    # a genuinely homogeneous flank: no nucleation, no p300, no KDM6B, uniform
    # exchange -- so the lattice and the mean-field see IDENTICAL rates.
    geom = Geometry(up_bp=2000, dn_bp=2000, core_up_bp=-1e9, core_dn_bp=-1e9,
                    cgi_sigma_bp=1e-6, body_start_bp=-1e6)
    print(f'  {"a_rw":>5} {"MF low":>7} {"MF high":>8} {"lat<-0":>8} {"lat<-1":>8}')
    rows = []
    for arw in [2.0, 2.75, 3.5, 4.5, 6.0]:
        par = replace(PAR, a_rw_prc2=arw, k_kdm6b=0.0, k_ac=0.0,
                      nucleation_gain=0.0)
        fps = mf.fixed_points(par, INP, grid=401, body=1.0,
                              psi_kdm=0.0, psi_p300=0.0, nucl=0.0)
        st = sorted([m for m, s in fps if s])
        lo_mf, hi_mf = st[0], st[-1]
        got = {}
        for init, lab in [(0.0, 'from 0'), (1.0, 'from 1')]:
            e = ensemble(n_cells=6, hours=500, dt=0.05, geom=geom, par=par,
                         inp=INP, init_me3=init, burn_frac=0.6)
            got[lab] = e['profile'].mean()
        print(f'  {arw:5.2f} {lo_mf:7.3f} {hi_mf:8.3f} {got["from 0"]:8.3f} {got["from 1"]:8.3f}')
        rows.append(dict(a_rw=arw, mf_low=lo_mf, mf_high=hi_mf,
                         lat_from0=got['from 0'], lat_from1=got['from 1']))
    json.dump(rows, open(os.path.join(OUT, 'lattice_branches.json'), 'w'), indent=2)


def exp_profile():
    e = ensemble(n_cells=30, hours=500, dt=0.05, par=PAR, inp=INP, init_me3=0.8)
    m = boundary_metrics(e['x'], e['profile'], e['in_core'])
    print('  ' + '  '.join(f'{k}={v:.3g}' for k, v in m.items()))
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    ax.fill_between(e['x'], e['profile'] - e['profile_sd'], e['profile'] + e['profile_sd'],
                    color='#b0402f', alpha=.18, lw=0)
    ax.plot(e['x'], e['profile'], color='#b0402f', lw=2, label='me3 (ensemble mean $\\pm$ sd)')
    ax.axvspan(e['x'][e['in_core']].min(), e['x'][e['in_core']].max(),
               color='#cfe3f5', alpha=.6, zorder=0, label='depleted core')
    if not np.isnan(m['edge_bp']):
        ax.axvline(m['edge_bp'], color='k', ls=':', lw=1)
        ax.text(m['edge_bp'] + 120, 0.05, f"edge {m['edge_bp']:.0f} bp", fontsize=8)
    ax.set_xlabel('position relative to TSS (bp)'); ax.set_ylabel('H3K27me3 occupancy')
    ax.set_title('Depleted core and bounded flanking domain emerge from local antagonism',
                 fontweight='bold', fontsize=10.5)
    ax.legend(fontsize=8)
    _save(fig, 'lattice_profile')
    json.dump(m, open(os.path.join(OUT, 'lattice_profile.json'), 'w'), indent=2)


def exp_versions():
    variants = [
        ('v1 local kernel',      replace(PAR, kernel='local')),
        ('v2 3D contact 1/d',    replace(PAR, kernel='contact', contact_exp=1.0)),
        ('v2b 3D contact 1/d^2', replace(PAR, kernel='contact', contact_exp=2.0)),
        ('v3 collapsed me chain', replace(PAR, k_me1=12.0, k_me2=12.0)),
        ('v4 no read-write',     replace(PAR, a_rw_prc2=0.0, nucleation_gain=3.0)),
        ('v5 no KDM6B',          replace(PAR, k_kdm6b=0.0)),
        ('v6 no p300',           replace(PAR, k_ac=0.0)),
    ]
    rows = []
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    for lab, par in variants:
        e = ensemble(n_cells=16, hours=500, dt=0.05, par=par, inp=INP, init_me3=0.8)
        m = boundary_metrics(e['x'], e['profile'], e['in_core'])
        rows.append(dict(variant=lab, **{k: (None if np.isnan(v) else float(v))
                                         for k, v in m.items()}))
        print(f'  {lab:22s} core={m["core_me3"]:.3f} plateau={m["plateau"]:.3f} '
              f'depl={m["depletion"]:+.3f} edge={m["edge_bp"]:.0f} width={m["width_bp"]:.0f}'
              .replace('nan', '  - '))
        ax.plot(e['x'], e['profile'], lw=1.7, label=lab)
    ax.set_xlabel('position relative to TSS (bp)'); ax.set_ylabel('H3K27me3 occupancy')
    ax.set_title('Model variants', fontweight='bold'); ax.legend(fontsize=7.5)
    _save(fig, 'lattice_versions')
    json.dump(rows, open(os.path.join(OUT, 'lattice_versions.json'), 'w'), indent=2)


def exp_gli():
    rows = []
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
    for g in [0.0, 0.3, 0.6, 1.0, 1.5, 2.2, 3.0]:
        inp = replace(INP, Gli1=g)
        e = ensemble(n_cells=16, hours=500, dt=0.05, par=PAR, inp=inp, init_me3=0.8)
        m = boundary_metrics(e['x'], e['profile'], e['in_core'])
        rows.append(dict(Gli1=g, **{k: (None if np.isnan(v) else float(v)) for k, v in m.items()}))
        print(f'  Gli1={g:4.1f} core={m["core_me3"]:.3f} plateau={m["plateau"]:.3f} '
              f'edge={m["edge_bp"]:.0f}'.replace('nan', ' - '))
        ax[0].plot(e['x'], e['profile'], lw=1.7, label=f'Gli1={g:g}')
    ax[0].set_xlabel('position vs TSS (bp)'); ax[0].set_ylabel('me3'); ax[0].legend(fontsize=7.5)
    ax[0].set_title('KDM6B dose (Gli1) widens the depleted core', fontweight='bold', fontsize=10)
    gg = [r['Gli1'] for r in rows]
    ax[1].plot(gg, [r['core_me3'] for r in rows], 'o-', color='#b0402f', label='core')
    ax[1].plot(gg, [r['plateau'] for r in rows], 's-', color='#2f6fb0', label='flank plateau')
    ax[1].set_xlabel('Gli1'); ax[1].set_ylabel('me3'); ax[1].legend(fontsize=8)
    ax[1].set_title('core vs flank response', fontweight='bold', fontsize=10)
    _save(fig, 'lattice_gli')
    json.dump(rows, open(os.path.join(OUT, 'lattice_gli.json'), 'w'), indent=2)


def exp_cycle():
    rows = []
    for ch in [None, 12.0, 18.0, 24.0, 36.0, 48.0, 96.0]:
        e = ensemble(n_cells=16, hours=600, dt=0.05, par=PAR, inp=INP,
                     init_me3=0.8, cycle_h=ch)
        m = boundary_metrics(e['x'], e['profile'], e['in_core'])
        lab = 'no division' if ch is None else f'{ch:g} h'
        rows.append(dict(cycle_h=ch, **{k: (None if np.isnan(v) else float(v)) for k, v in m.items()}))
        print(f'  cycle={lab:>11s} core={m["core_me3"]:.3f} plateau={m["plateau"]:.3f}')
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    xs = [r['cycle_h'] if r['cycle_h'] else 200 for r in rows]
    ax.plot(xs, [r['plateau'] for r in rows], 'o-', color='#b0402f', label='flank plateau')
    ax.plot(xs, [r['core_me3'] for r in rows], 's-', color='#2f6fb0', label='core')
    ax.set_xlabel('cell-cycle length (h; 200 = no division)'); ax.set_ylabel('me3')
    ax.set_title('Replicative dilution vs cycle length', fontweight='bold', fontsize=10.5)
    ax.legend(fontsize=8)
    _save(fig, 'lattice_cycle')
    json.dump(rows, open(os.path.join(OUT, 'lattice_cycle.json'), 'w'), indent=2)


def exp_switching():
    """Is the lattice unimodal in the mean-field's bistable window, or is it
    hopping between the two branches? Mean-field cannot tell these apart, and
    they mean completely different things biologically."""
    geom = Geometry(up_bp=2000, dn_bp=2000, core_up_bp=-1e9, core_dn_bp=-1e9,
                    cgi_sigma_bp=1e-6, body_start_bp=-1e6)
    fig, axes = plt.subplots(2, 3, figsize=(13, 6.4))
    rows = []
    for ax, arw in zip(axes.ravel(), [2.75, 3.5, 4.5]):
        par = replace(PAR, a_rw_prc2=arw, k_kdm6b=0.0, k_ac=0.0, nucleation_gain=0.0)
        fps = mf.fixed_points(par, INP, grid=401, body=1.0, psi_kdm=0.0,
                              psi_p300=0.0, nucl=0.0)
        st = sorted([m for m, s_ in fps if s_])
        traces = []
        for c in range(12):
            lat = PromoterLattice(geom, par, seed=100 + c, init_me3=0.5)
            out = lat.run(900, 0.05, INP, record_every=2.0)
            traces.append(out['profile'].mean(1))
        T = np.array(traces); tt = out['t']
        for tr in T:
            ax.plot(tt, tr, lw=.7, alpha=.75)
        for m in st:
            ax.axhline(m, color='k', ls='--', lw=.9)
        ax.set_title(f'$a_{{rw}}$={arw}  (MF stable: {", ".join(f"{m:.2f}" for m in st)})',
                     fontsize=9.5)
        ax.set_xlabel('h'); ax.set_ylabel('flank me3'); ax.set_ylim(0, 1)
        # distribution over the second half
        half = T[:, len(tt)//2:].ravel()
        axd = axes.ravel()[3 + list([2.75, 3.5, 4.5]).index(arw)]
        axd.hist(half, bins=np.linspace(0, 1, 36), color='#b0402f', alpha=.8)
        for m in st:
            axd.axvline(m, color='k', ls='--', lw=.9)
        axd.set_xlabel('flank me3'); axd.set_ylabel('count'); axd.set_title('distribution', fontsize=9)
        # crude bimodality: dip between the two MF branches
        hist, edges = np.histogram(half, bins=np.linspace(0, 1, 36), density=True)
        rows.append(dict(a_rw=arw, mf_stable=st, mean=float(half.mean()),
                         sd=float(half.std()),
                         frac_above_mid=float((half > 0.45).mean())))
        print(f'  a_rw={arw:4.2f}  MF stable={[round(m,3) for m in st]}  '
              f'lattice mean={half.mean():.3f} sd={half.std():.3f}  '
              f'frac>0.45={rows[-1]["frac_above_mid"]:.2f}')
    fig.suptitle('Does the lattice occupy the mean-field branches, or smear through them?',
                 fontweight='bold')
    _save(fig, 'lattice_switching')
    json.dump(rows, open(os.path.join(OUT, 'lattice_switching.json'), 'w'), indent=2)


EXPS = dict(bifurcation=exp_bifurcation, switching=exp_switching, branches=exp_branches, profile=exp_profile,
            versions=exp_versions, gli=exp_gli, cycle=exp_cycle)

if __name__ == '__main__':
    todo = sys.argv[1:] or list(EXPS)
    for name in todo:
        print(f'\n=== {name} ===')
        t0 = time.time()
        EXPS[name]()
        print(f'  ({time.time()-t0:.0f}s)')
