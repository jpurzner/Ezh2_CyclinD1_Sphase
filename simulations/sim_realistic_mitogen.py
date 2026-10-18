"""Biologically realistic, EVIDENCE-DRIVEN Shh mitogen input + a DECOMPOSITION of the Gli signal, and the
demonstration that the pathway LOW-PASS-filters it (JP 2026-07-21; grounded in docs/hh-signalling-timing.md).

Evidence-driven design (see the doc):
  * LIGAND L(t): tonic / slowly-varying paracrine (Purkinje-derived); NO fast pulse train (that is ERK, not Hh).
    Modelled as tonic + a SLOW OU drift (tau_lo ~ 40 h, the "real" slow variation, PASSBAND) + a FAST OU noise
    (tau_hi ~ 0.5 h, sub-hour, STOPBAND) that the pathway should filter out. (GNP ligand fluctuation timescale is
    UNMEASURED -> the amplitudes are a stated assumption; the point is the timescale SEPARATION.)
  * GLI is quasi-static in v44: it tracks the ligand with a FLAT gain (~0.15) at all periods (measured, panel D) --
    consistent with the doc (Gli responds on the ~40 min-2 h protein timescale = "instant" vs the ~22 h cycle; no
    intrinsic Gli oscillation). So the Gli signal itself is NOISY -- it carries the fast ligand fluctuations.
  * THE LOW-PASS FILTER IS DOWNSTREAM: the slow CyclinD1 mRNA RESERVOIR (tau ~ 11-12 h, the baked Hh memory) plus the
    EZH2 incoherent feedforward. The measured transfer function ligand->Cd_mRNA is a clean low-pass: fast fluctuations
    (period < 2 h) attenuated ~20-40x, slow drift (period > 40 h) passed. This is the doc's "temporal averaging by
    slow mRNA turnover" (sec 7) and it is what buffers the PROLIFERATION decision against Gli noise.
  * CILIUM duty cycle G(t): ~94% of EGL GNP progenitors ciliated, cilia RESORB pre-mitotically -> Shh reception ON in
    G1/S, OFF at G2/M (cilium_gate on CyclinB). A cell-cycle-coupled multiplier; ON-fraction tuned to ~0.9.
  * ADAPTATION A(t): the ~20 h temporal adaptation (Ptch1 feedback + Gli2-down + GliA-destabilization); STRENGTH is
    the central GNP unknown (GNPs proliferate for days under sustained Shh) -> a scanned free parameter (Ptch1_slow).

Outputs (-> _cache.pkl; figure = fig_realistic_mitogen.py):
  (1) the realistic ligand + its slow/fast components;
  (2) the Gli-signal DECOMPOSITION (tonic -> +noise -> +cilium -> +adaptation);
  (3) the reservoir low-passing it (Gli noisy vs Cd_mRNA smooth, same input);
  (4) the BODE transfer function ligand->{Gli, Cd_mRNA} (the quantitative low-pass);
  (5) cilium duty-cycle waveform + ON fraction;
  (6) cell-to-cell heterogeneity.
Run: ./venv/bin/python simulations/sim_realistic_mitogen.py"""
import os, sys, pickle
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1'); os.environ.setdefault('EZH2_CONC', '0')
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44
HERE = os.path.dirname(os.path.abspath(__file__)); UPH = 62.8
GNP = dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464)
ADAPT = dict(k_ptch1_slow_on=0.008, w_ptch1_slow=1.5, k_ptch1_slow_off=0.0005)   # ~20 h WEAK adaptation (free strength; GNPs proliferate for days -> weak)
CILIUM = dict(w_cilium=0.9, K_cil=0.70)                                          # duty cycle; K_cil -> ON-fraction ~0.9
KTL0, KTLEZ0 = 0.801, 0.0150988626303029


def ligand(T_h, dt_h=0.1, tonic=0.45, sd_lo=0.075, tau_lo=40.0, sd_hi=0.060, tau_hi=0.5,
           lo=0.16, hi=0.85, seed=7):
    """Tonic Shh + SLOW OU drift (paracrine, passband) + FAST OU noise (sub-hour, stopband -> to be filtered).
    OU stationary sd = sigma_step*sqrt(tau/2)  ->  sigma_step = sd*sqrt(2/tau). Components re-centred so the
    realized ligand mean == tonic (a controlled proliferating setpoint, not a per-realization offset)."""
    rng = np.random.default_rng(seed)
    t = np.arange(0, T_h, dt_h); n = len(t)
    def ou(tau, sd):
        s = sd * np.sqrt(2.0 / tau); x = np.zeros(n)
        for i in range(1, n):
            x[i] = x[i-1] - x[i-1] * dt_h / tau + s * np.sqrt(dt_h) * rng.normal()
        return x - x.mean()                                # re-centre so tonic sets the mean
    slow = ou(tau_lo, sd_lo); fast = ou(tau_hi, sd_hi)
    L = np.clip(tonic + slow + fast, lo, hi)
    L_slow = np.clip(tonic + slow, lo, hi)                  # slow-only (what SHOULD pass the filter)
    return t, L, L_slow, slow, fast


def run(L_t, L_arr, over, T_h, sel, chunk=6.0):
    """Piecewise-integrate the model under a time-varying SHH = L(t) (chunk in model units)."""
    rr = te.loada(build_model_v44(params=over)); rr.reset()
    for k, v in GNP.items(): rr[k] = v
    rr.integrator.setValue('relative_tolerance', 1e-6); rr.integrator.setValue('maximum_num_steps', 2000000)
    def Lf(th): return float(np.interp(th, L_t, L_arr))
    T, out = [], {k: [] for k in sel if k != 'time'}
    t = 0.0
    while t < T_h * UPH - 1e-6:
        rr['SHH'] = Lf(t / UPH)
        r = None
        for atol in (1e-8, 1e-7, 1e-6, 1e-5):
            try: rr.integrator.setValue('absolute_tolerance', atol); r = rr.simulate(t, t + chunk, 3, selections=sel); break
            except Exception: rr.integrator.setValue('maximum_time_step', 2.0)
        if r is None: t += chunk; continue
        T += list(np.asarray(r['time'])[1:])
        for k in out: out[k] += list(np.asarray(r[k])[1:])
        t += chunk
    return np.array(T) / UPH, {k: np.array(v) for k, v in out.items()}


def lockin(period_h, tonic=0.45, amp=0.08, over=None, n_settle=4, n_meas=10, ppc=40):
    """Transfer function: drive SHH = tonic + amp*sin(2pi t/period); extract the Fourier component of each readout
    at the drive frequency (lock-in) -> gain = |readout|/|ligand|. Separates the driven response from cycle noise."""
    rr = te.loada(build_model_v44(params=(over or dict(w_cilium=0)))); rr.reset()
    for k, v in GNP.items(): rr[k] = v
    rr['SHH'] = tonic; rr.integrator.setValue('relative_tolerance', 1e-7); rr.integrator.setValue('absolute_tolerance', 1e-9)
    rr.integrator.setValue('maximum_num_steps', 3000000)
    rr.simulate(0, 180 * UPH, 50)
    per = period_h * UPH; dt = per / ppc; t = 180 * UPH
    ts = []; sig = {'Gli_act': [], 'Cd_mRNA': [], 'Cd': []}; lig = []
    for _ in range(int((n_settle + n_meas) * ppc)):
        shh = tonic + amp * np.sin(2 * np.pi * (t - 180 * UPH) / per); rr['SHH'] = max(shh, 0.001)
        r = rr.simulate(t, t + dt, 2, selections=['time', 'Gli_act', 'Cd_mRNA', 'Cd'])
        ts.append(r['time'][-1])
        for k in sig: sig[k].append(r[k][-1])
        lig.append(shh); t += dt
    ts = np.array(ts); lig = np.array(lig); m = ts > (180 * UPH + n_settle * per)
    ph = 2 * np.pi * (ts - 180 * UPH) / per
    def famp(x):
        x = np.array(x); c = 2 * np.mean((x[m] - x[m].mean()) * np.cos(ph[m])); s = 2 * np.mean((x[m] - x[m].mean()) * np.sin(ph[m]))
        return np.hypot(c, s)
    la = famp(lig)
    return {k: famp(v) / la for k, v in sig.items()}


if __name__ == '__main__':
    T_h = 160.0
    t, L_full, L_slow, slow, fast = ligand(T_h)
    res = dict(t_lig=t, L_full=L_full, L_slow=L_slow, slow=slow, fast=fast,
               lig_stats=dict(mean=float(L_full.mean()), cv=float(L_full.std() / L_full.mean()),
                              cv_slow=float(slow.std() / L_full.mean()), cv_fast=float(fast.std() / L_full.mean())))
    print(f"ligand: mean={L_full.mean():.3f} cv={L_full.std()/L_full.mean():.3f} "
          f"(slow {slow.std()/L_full.mean():.3f} + fast {fast.std()/L_full.mean():.3f})", flush=True)
    sel = ['time', 'Gli_act', 'Gli1', 'Cd_mRNA', 'Cd', 'cilium_gate', 'Cb', 'Dna', 'Ptch1_slow']

    # (2) Gli-signal decomposition: add one component at a time
    print("decomposition (tonic -> +noise -> +cilium -> +adaptation)...", flush=True)
    tonic = float(np.mean(L_full))
    res['d1_tonic'] = run(t, np.full_like(t, tonic), dict(w_cilium=0), T_h, sel)
    res['d2_noise'] = run(t, L_full, dict(w_cilium=0), T_h, sel)
    res['d3_cilium'] = run(t, L_full, CILIUM, T_h, sel)
    res['d4_full'] = run(t, L_full, dict(**CILIUM, **ADAPT), T_h, sel)
    for k in ('d1_tonic', 'd2_noise', 'd3_cilium', 'd4_full'):
        g = res[k][1]['Gli_act']; cm = res[k][1]['Cd_mRNA']
        print(f"  {k}: Gli mean={g.mean():.3f} cv={g.std()/g.mean():.3f} | Cd_mRNA cv={cm.std()/cm.mean():.3f}", flush=True)

    # (3) low-pass check on the realistic input: Gli CV vs Cd_mRNA CV (transcript smooths the Gli noise)
    gd = res['d3_cilium'][1]
    res['lowpass_realistic'] = dict(gli_cv=float(gd['Gli_act'].std() / gd['Gli_act'].mean()),
                                    cd_mRNA_cv=float(gd['Cd_mRNA'].std() / gd['Cd_mRNA'].mean()),
                                    cd_cv=float(gd['Cd'].std() / gd['Cd'].mean()))
    print(f"  realistic-input low-pass: Gli CV={res['lowpass_realistic']['gli_cv']:.3f} -> "
          f"Cd_mRNA CV={res['lowpass_realistic']['cd_mRNA_cv']:.3f} "
          f"({res['lowpass_realistic']['gli_cv']/max(res['lowpass_realistic']['cd_mRNA_cv'],1e-9):.1f}x smoother)", flush=True)

    # (4) BODE transfer function ligand -> {Gli, Cd_mRNA, Cd}
    rr0 = te.loada(build_model_v44()); res['reservoir_tau_h'] = float(1.0 / rr0['k_Cd_mRNA_deg'] / UPH)   # first-order low-pass tau
    print(f"Bode transfer function (lock-in, cilium off); reservoir tau={res['reservoir_tau_h']:.1f}h...", flush=True)
    periods = [0.5, 1, 2, 4, 8, 12, 16, 24, 32, 48, 72, 110, 160]
    bode = {'period': [], 'Gli_act': [], 'Cd_mRNA': [], 'Cd': []}
    for p in periods:
        try:
            g = lockin(p); bode['period'].append(p)
            for k in ('Gli_act', 'Cd_mRNA', 'Cd'): bode[k].append(g[k])
            print(f"  period={p:6.1f}h  Gli={g['Gli_act']:.3f}  Cd_mRNA={g['Cd_mRNA']:.3f}  Cd={g['Cd']:.3f}", flush=True)
        except Exception as e:
            print(f"  period={p:6.1f}h ERROR {e}", flush=True)
    res['bode'] = {k: np.array(v) for k, v in bode.items()}

    # (5) cilium duty-cycle ON fraction (tune reference)
    cg = res['d3_cilium'][1]['cilium_gate']; res['cilium_on_fraction'] = float((cg > 0.5).mean())
    print(f"  cilium duty cycle ON-fraction = {res['cilium_on_fraction']:.2f} (target ~0.9; 94% ciliated - G2/M off)", flush=True)

    # (6) cell-to-cell heterogeneity (different ligand realizations + per-cell CyclinD1/EZH2 tone)
    print("cell-to-cell heterogeneity...", flush=True)
    rng = np.random.default_rng(2028); cells = []
    for c in range(6):
        _, Lc, Lc_slow, _, _ = ligand(T_h, seed=100 + c)
        cd = float(np.clip(np.exp(rng.normal(0, 0.5)), 0.4, 2.5)); ez = float(np.clip(np.exp(rng.normal(0, 0.4)), 0.5, 2.0))
        ov = dict(k_Cd_translation=KTL0 * cd, kTlEZ=KTLEZ0 * ez, **CILIUM, **ADAPT)
        tc, oc = run(t, Lc, ov, T_h, ['time', 'Gli_act', 'Cd_mRNA', 'Dna'])
        cells.append((tc, oc['Gli_act'], oc['Cd_mRNA'], Lc))
    res['cells'] = cells
    pickle.dump(res, open(os.path.join(HERE, 'sim_realistic_mitogen_cache.pkl'), 'wb'))
    print("-> cached; figure: fig_realistic_mitogen.py", flush=True)
