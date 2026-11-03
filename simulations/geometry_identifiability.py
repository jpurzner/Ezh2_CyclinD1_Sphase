"""Which lattice parameters can be recovered from an average H3K27me3 PROFILE?

Answers the question before the data arrives: a metagene profile gives 7 shape
features; the lattice has many more parameters. This computes the sensitivity
Jacobian dF/dln(p), then its SVD, to show which parameter COMBINATIONS the
profile constrains and which are degenerate.

  python simulations/geometry_identifiability.py [--cells N]
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src', 'lattice'))
import numpy as np
from dataclasses import replace
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from promoter_lattice import Geometry, Params, Inputs, ensemble
from geometry_fit import profile_features, features_vector, FEATURES

OUT = os.path.dirname(os.path.abspath(__file__))
# symmetric, generous extents so neither side is truncated
GEOM = Geometry(up_bp=6000, dn_bp=9000, core_up_bp=400, core_dn_bp=700)
PAR = Params(a_rw_prc2=3.5)
INP = Inputs(Gli1=0.3, activity=0.4, tx=0.4)
PARAMS = ['a_rw_prc2', 'k_kdm6b', 'kdm_core_width', 'p_local',
          'del_mk', 'k_h33', 'nucleation_gain', 'a0_prc2']


def run(par, cells, seed0=0):
    e = ensemble(n_cells=cells, hours=500, dt=0.05, geom=GEOM, par=par,
                 inp=INP, init_me3=0.8, seed0=seed0)
    return profile_features(e['x'], e['profile'])


def main():
    cells = 28
    if '--cells' in sys.argv:
        cells = int(sys.argv[sys.argv.index('--cells') + 1])

    print(f'baseline ({cells} cells/point)')
    base = run(PAR, cells)
    for k in FEATURES:
        print(f'  {k:11s} {base[k]:9.3f}')

    # noise floor: same parameters, different seeds
    reps = [features_vector(run(PAR, cells, seed0=1000 * i)) for i in range(4)]
    noise = np.nanstd(np.array(reps), axis=0)
    print('\nnoise floor (sd over 4 independent ensembles):')
    for k, v in zip(FEATURES, noise):
        print(f'  {k:11s} {v:9.3f}')

    print('\nsensitivity  dF/dln(p)  (+/-25% perturbation)')
    J, used = [], []
    hdr = '  ' + f'{"parameter":16s}' + ''.join(f'{k[:9]:>10s}' for k in FEATURES) + '   max|S/noise|'
    print(hdr)
    for pname in PARAMS:
        v0 = getattr(PAR, pname)
        hi = run(replace(PAR, **{pname: v0 * 1.25}), cells)
        lo = run(replace(PAR, **{pname: v0 * 0.75}), cells)
        d = (features_vector(hi) - features_vector(lo)) / np.log(1.25 / 0.75)
        snr = np.nanmax(np.abs(d) / np.where(noise > 0, noise, np.nan))
        J.append(d); used.append(pname)
        print('  ' + f'{pname:16s}' + ''.join(
            ('      nan' if not np.isfinite(x) else f'{x:10.1f}') for x in d) + f'   {snr:8.1f}')
    J = np.array(J)

    # scale rows by the noise floor so the SVD is in units of "detectability"
    Jn = J / np.where(noise > 0, noise, np.nan)
    Jn = np.nan_to_num(Jn)
    U, S, Vt = np.linalg.svd(Jn, full_matrices=False)
    print('\nSVD of the noise-scaled Jacobian (parameter space):')
    print(f'  singular values: {np.array2string(S, precision=1)}')
    print(f'  condition number: {S[0]/max(S[-1],1e-12):.3g}')
    print('\n  identifiable combinations (top 3):')
    for i in range(min(3, len(S))):
        terms = sorted(zip(used, U[:, i]), key=lambda t: -abs(t[1]))[:3]
        print(f'    #{i+1} (sv {S[i]:7.1f}): ' +
              '  '.join(f'{n}{"+" if w>0 else "-"}{abs(w):.2f}' for n, w in terms))
    print('\n  WORST-determined combinations (bottom 2):')
    for i in range(len(S) - 2, len(S)):
        terms = sorted(zip(used, U[:, i]), key=lambda t: -abs(t[1]))[:3]
        print(f'    (sv {S[i]:7.2f}): ' +
              '  '.join(f'{n}{"+" if w>0 else "-"}{abs(w):.2f}' for n, w in terms))

    fig, ax = plt.subplots(figsize=(9, 4.4))
    im = ax.imshow(np.clip(Jn, -30, 30), cmap='RdBu_r', aspect='auto', vmin=-30, vmax=30)
    ax.set_xticks(range(len(FEATURES))); ax.set_xticklabels(FEATURES, rotation=35, ha='right', fontsize=8)
    ax.set_yticks(range(len(used))); ax.set_yticklabels(used, fontsize=8)
    ax.set_title('Profile sensitivity, scaled by the noise floor\n'
                 '(how many sd a 25% parameter change moves each feature)',
                 fontweight='bold', fontsize=10)
    fig.colorbar(im, label='sd')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'geometry_identifiability.png'), dpi=140, bbox_inches='tight')
    json.dump(dict(baseline={k: float(base[k]) for k in FEATURES},
                   noise=dict(zip(FEATURES, noise.tolist())),
                   params=used, jacobian=J.tolist(),
                   singular_values=S.tolist()),
              open(os.path.join(OUT, 'geometry_identifiability.json'), 'w'), indent=2)
    print('\n  wrote geometry_identifiability.png/.json')


if __name__ == '__main__':
    main()
