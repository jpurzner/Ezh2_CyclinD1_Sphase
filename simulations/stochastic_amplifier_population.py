"""Does a STOCHASTIC chromatin layer turn the two-compartment amplifier's all-or-none
cliff into a graded, tunable locked fraction?  (JP question, 2026-10-22)

BACKGROUND.  `population_amplifier.py` concluded that the amplifier is either inert
(validation-safe regime) or a permanent lock (data-excluded regime), with "no robust
moderately-prolonged transient dwell in between".  But in that script the only
variability is a lognormal over EXCURSION MAGNITUDE -- the driving input.  The
chromatin layer itself (Mk, P_prox) is a deterministic ODE, and a deterministic
bistable switch can only ever give OFF-or-locked.

The single-cell chromatin literature says the in-between regime exists but is realised
as HETEROGENEITY ACROSS CELLS: silencing is all-or-none per cell and the regulator
tunes the FRACTION of cells silenced, not the level (Bintu 2016; Anink-Groenen 2014's
boundary diffuses).  If that is right, an intermediate a_P that reads as "inert"
deterministically should give a tunable locked FRACTION once P_prox fluctuates.

WHAT THIS SCRIPT ADDS.  A chemical-Langevin noise term on P_prox, representing a
promoter with a FINITE number of proximal nucleosome sites (Nsites):

    dP_prox/dt = lam_on - lam_off + sqrt((lam_on + lam_off)/Nsites) * eta(t)

the standard CLE for a birth-death process on Nsites sites (stationary variance
P(1-P)/Nsites = the binomial result).  Nsites -> inf recovers the current model exactly.

IMPLEMENTATION NOTE -- why the noise is smooth, not a jump process.
The obvious implementations (chunked integration + SSA, or timed SBML events) both
force a CVODE restart, and v44 CANNOT BE RESTARTED during deep arrest: chunked
integration dies at t~6250 min with CV_CONV_FAILURE.  This is PRE-EXISTING -- it
reproduces with with_proximal_distal=False, and no tolerance or integrator setting
rescues it (see the docstring note on "deep-arrest integrator stiffness").  So eta(t)
is synthesised spectrally as a band-limited Gaussian process (random-phase sum of
cosines), which is smooth and needs no restart.  Band-limited rather than white is
also the more physical choice: chromatin fluctuations are correlated, not delta-
correlated.  --selftest calibrates the amplitude against the exact binomial variance.

Run:
  PYTHONPATH=. ./venv/bin/python simulations/stochastic_amplifier_population.py --selftest
  PYTHONPATH=. ./venv/bin/python simulations/stochastic_amplifier_population.py
"""
import sys, os, time, argparse, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

BP = dict(kPhRbCd=0.1748299602292031, K_CdRb=0.4025007594938155, w_ink4=4.646376316224687,
          kSyP21=0.0004043793202695936, K_cdk6_sink=3.879643001142466, w_p18=1.003950915185932,
          w_p19=0.6655007434943913, k_Cd_tx_Gli_max=0.24197652013556434, k_Cd_tx_MYCN=0.11031706915457935)

T_EXC, T_POST = 22000, 90000          # excursion / post-excursion window (min), as in population_amplifier.py
M_MODES = 48                           # spectral modes in the noise path
TAU_LO, TAU_HI = 40.0, 1500.0          # noise band: periods 40..1500 min (P_prox relaxation is ~200-500 min)


# ---------------------------------------------------------------- noise path
def mode_freqs():
    """Angular frequencies, log-spaced across the band."""
    return 2 * np.pi / np.exp(np.linspace(np.log(TAU_LO), np.log(TAU_HI), M_MODES))


def eta_numpy(t, phases, amps):
    """The same band-limited path the Antimony model evaluates (for the selftest)."""
    w = mode_freqs()
    return (amps[None, :] * np.cos(w[None, :] * np.atleast_1d(t)[:, None] + phases[None, :])).sum(1)


CAL_CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), '.noise_calibration.json')


def calibrate_amplitude(seed=0, n_real=400):
    """Scale the random-phase cosine sum so it acts as unit-intensity white noise.

    For a birth-death process the CLE stationary variance is P(1-P)/N.  We integrate
    the scalar OU-like system dP = (lam_on - lam_off)dt + sig*eta dt around a fixed
    operating point and solve for the scale that reproduces it.
    """
    if os.path.exists(CAL_CACHE):
        with open(CAL_CACHE) as fh:
            c = json.load(fh)
        if c.get('modes') == M_MODES and c.get('band') == [TAU_LO, TAU_HI]:
            return c['scale'], c['got'], c['target']
    rng = np.random.default_rng(seed)
    w = mode_freqs()
    # unit-amplitude cosines, scale solved below
    base = np.ones(M_MODES)
    # operating point: symmetric birth-death, P* = 0.5, relaxation rate k
    k, N, Pstar = 0.01, 12.0, 0.5
    sig = np.sqrt((k * Pstar + k * (1 - Pstar)) / N)        # CLE amplitude at P*
    target = Pstar * (1 - Pstar) / N                         # binomial variance
    dt, T = 1.0, 400000
    tt = np.arange(0, T, dt)
    var = []
    for _ in range(n_real):
        ph = rng.uniform(0, 2 * np.pi, M_MODES)
        e = (base[None, :] * np.cos(w[None, :] * tt[:, None] + ph[None, :])).sum(1)
        P = np.empty_like(tt); P[0] = Pstar
        for i in range(1, len(tt)):                          # Euler, linearised
            P[i] = P[i-1] + (-k * (P[i-1] - Pstar)) * dt + sig * e[i-1] * dt
        var.append(P[len(tt)//4:].var())
    scale = float(np.sqrt(target / np.mean(var)))
    with open(CAL_CACHE, 'w') as fh:
        json.dump(dict(scale=scale, got=float(np.mean(var)), target=float(target),
                       modes=M_MODES, band=[TAU_LO, TAU_HI]), fh)
    return scale, float(np.mean(var)), float(target)


# ---------------------------------------------------------------- model
def noise_block(nscale):
    w = mode_freqs()
    terms = " + ".join(f"cos({wj:.10g}*time + ph_{j})" for j, wj in enumerate(w))
    phs = "".join(f"\n  ph_{j} = 0;" for j in range(M_MODES))
    return (
        "\n  // ---- stochastic chromatin layer: CLE noise on P_prox (Nsites proximal sites) ----"
        "\n  Nsites = 12; nscale = %.10g; eps_P = 0.01;" % nscale
        + phs +
        "\n  lam_on_P := k_onP*a_P*Mk^n_anch*(1 - elong_gate)*(1 - P_prox);"
        "\n  lam_off_P := (k_offP*elong_gate + basal_offP)*P_prox;"
        "\n  gclamp := (P_prox/(P_prox + eps_P))*((1 - P_prox)/(1 - P_prox + eps_P));"
        "\n  xi_t := " + terms + ";"
        "\n  noise_amp := nscale*sqrt((lam_on_P + lam_off_P)/Nsites)*gclamp;"
        "\n  P_noise: => P_prox; Cell*noise_amp*xi_t;"
    )


def build(nscale, a_P, f_P=0.80):
    import tellurium as te
    from src.build_model_v44_heldt import build_model_v44
    P = dict(BP); P.update(dict(a_P=a_P, f_P=f_P))
    src = build_model_v44(with_p27_optionB=True, with_cdk6_gli=True,
                          with_proximal_distal=True, params=P)
    src = src.replace("\nend", noise_block(nscale) + "\nend", 1)
    rr = te.loada(src)
    rr.integrator.setValue('relative_tolerance', 1e-7)
    rr.integrator.setValue('absolute_tolerance', 1e-9)
    rr.integrator.setValue('maximum_num_steps', 15000000)
    return rr, P


def setp(rr, P, mult):
    for k, v in P.items():
        try: rr[k] = v
        except Exception: pass
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 0.3
    try: rr['Cd2_expr'] = rr['CD2_EXPR_MB']
    except Exception: pass
    rr['p18'] = 1.73 * mult; rr['p19'] = 0.58 * mult; rr['kSyP21'] = BP['kSyP21'] * mult


def set_cell(rr, rng, nsites, noise_on):
    for j in range(M_MODES):
        rr[f'ph_{j}'] = float(rng.uniform(0, 2 * np.pi))
    rr['Nsites'] = float(nsites)
    if not noise_on:
        rr['nscale'] = 0.0


def dwell(rr, P, mult, rng, nsites, noise_on, nscale):
    """Excursion to `mult` x CDKI then back to 1x. Returns dwell (h) or None if locked."""
    setp(rr, P, mult); rr.reset()
    rr['nscale'] = nscale
    set_cell(rr, rng, nsites, noise_on)
    setp(rr, P, mult)
    try:
        rr.simulate(0, T_EXC, T_EXC // 10, selections=['time', 'Dna'])
        setp(rr, P, 1.0)
        r2 = rr.simulate(T_EXC, T_EXC + T_POST, T_POST // 10,
                         selections=['time', 'Dna', 'P_prox'])
    except Exception:
        return 'ERR'
    dna = r2['Dna']; dd = np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0]
    return float((r2['time'][dd[0]] - T_EXC) / 60) if len(dd) else None


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--cells', type=int, default=24)
    ap.add_argument('--nsites', type=float, default=12.0)
    ap.add_argument('--mult', type=float, default=5.0, help='CDKI excursion magnitude')
    ap.add_argument('--a_P', type=float, default=0.10)
    ap.add_argument('--grid', type=float, nargs='+',
                    default=[0.72, 0.74, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.85])
    ap.add_argument('--out', default=None)
    a = ap.parse_args()

    t0 = time.time()
    scale, got, target = calibrate_amplitude()
    print(f'[calibration] random-phase scale = {scale:.6g}  '
          f'(raw var {got:.4g} -> target binomial var {target:.4g})')
    if a.selftest:
        # verify: scaled path reproduces the binomial variance at a second operating point
        rng = np.random.default_rng(7); w = mode_freqs()
        k, N, Pstar = 0.004, 20.0, 0.3
        sig = np.sqrt((k * Pstar + k * (1 - Pstar)) / N); tgt = Pstar * (1 - Pstar) / N
        dt = 1.0; tt = np.arange(0, 400000, dt); vs = []
        for _ in range(200):
            ph = rng.uniform(0, 2 * np.pi, M_MODES)
            e = scale * np.cos(w[None, :] * tt[:, None] + ph[None, :]).sum(1)
            P = np.empty_like(tt); P[0] = Pstar
            for i in range(1, len(tt)):
                P[i] = P[i-1] + (-k * (P[i-1] - Pstar)) * dt + sig * e[i-1] * dt
            vs.append(P[len(tt)//4:].var())
        print(f'[selftest] second operating point: var {np.mean(vs):.4g} vs binomial {tgt:.4g} '
              f'(ratio {np.mean(vs)/tgt:.3f})')
        print(f'[selftest] done in {time.time()-t0:.1f}s')
        return

    GRID = a.grid          # f_P values; f_P is the controlling axis (a_P alone does not reach the cliff)
    rows = []
    for f_P in GRID:
        rr, P = build(scale, a.a_P, f_P)
        for noise_on, lab in [(False, 'deterministic'), (True, f'stochastic N={a.nsites:g}')]:
            rng = np.random.default_rng(1234)
            ds = [dwell(rr, P, a.mult, rng, a.nsites, noise_on, scale if noise_on else 0.0)
                  for _ in range(1 if not noise_on else a.cells)]
            errs = sum(1 for d in ds if d == 'ERR')
            ds = [d for d in ds if d != 'ERR']
            locked = sum(1 for d in ds if d is None)
            fin = [d for d in ds if d is not None]
            frac = 100.0 * locked / len(ds) if ds else float('nan')
            rows.append(dict(f_P=f_P, a_P=a.a_P, mode=lab, n=len(ds), errs=errs, pct_locked=frac,
                             median_dwell=(float(np.median(fin)) if fin else None)))
            md = rows[-1]['median_dwell']
            print(f'f_P={f_P:<6g} {lab:<18s} n={len(ds):<3d} locked={frac:5.1f}%  '
                  f'median dwell={"lock" if md is None else f"{md:.0f}h"}' +
                  (f'  [{errs} integration errors]' if errs else ''))

    out = a.out or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'fig_v44_stochastic_amplifier')
    with open(out + '.json', 'w') as fh:
        json.dump(dict(rows=rows, nsites=a.nsites, cells=a.cells, mult=a.mult, a_P=a.a_P,
                       scale=scale, band=[TAU_LO, TAU_HI], modes=M_MODES), fh, indent=2)

    det = {r['f_P']: r for r in rows if r['mode'] == 'deterministic'}
    sto = {r['f_P']: r for r in rows if r['mode'] != 'deterministic'}
    xs = GRID
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.step(xs, [det[x]['pct_locked'] for x in xs], where='mid', color='#888888',
            lw=2, label='deterministic (current model)')
    ax.plot(xs, [sto[x]['pct_locked'] for x in xs], 'o-', color='#b0402f', lw=2, ms=6,
            label=f'stochastic chromatin (Nsites={a.nsites:g}, {a.cells} cells)')
    ax.axvline(0.80, color='k', ls=':', lw=0.9)
    ax.text(0.802, 55, 'BAKED default\n$f_P$=0.80', fontsize=7.5)
    ax.set_xlabel('$f_P$  (repression floor; LOWER = stronger proximal repression)')
    ax.invert_xaxis()
    ax.set_ylabel('% of cells permanently locked')
    ax.set_title('Does chromatin noise smear the lock cliff into a graded fraction?', fontweight='bold')
    ax.legend(fontsize=8.5); ax.set_ylim(-3, 103)
    fig.tight_layout()
    fig.savefig(out + '.png', dpi=140, bbox_inches='tight'); fig.savefig(out + '.pdf', bbox_inches='tight')
    print('wrote', out + '.png/.pdf/.json', f'({time.time()-t0:.0f}s)')


if __name__ == '__main__':
    main()
