"""v45 PROTOTYPE — stochastic commitment kernel (parallel to v44, NOT a replacement yet).

Replaces the v44 size-gate (M_commit) with a BISTABLE CDK2-p27 toggle + NOISY birth partitioning, so
the G0-vs-fast-G1 split EMERGES from a bifurcation (Spencer/Cappell/Yao) instead of being a phase every
cell traverses. Cyclin D1 = mitogen sensor that sets the basin probability (p27 synthesis), NOT a timer.
Rb is made size-independently so [Rb]=Rb/mass DILUTES as the cell grows = the committed fast-G1 timer
(Zatulovskiy 2020). Division partitions p27 (inherited brake, noisy) and CDK2 (Spencer carryover, noisy)
across daughters -> they straddle the commitment separatrix -> a born-committed fraction + a G0 tail.

Design goal: a committed-cycler period ~16-18h (which supplies ~68/16/16 count-fractions) PLUS a G0 dwell
on a fraction of daughters lifting the population MEAN toward ~22h, with the MB-vs-GNP difference emerging
from mitogen (more CyclinD1 -> less p27 synth -> fewer G0-born). Falsifiable: right-skewed period
distribution; born-committed vs G0 split = CDK2inc/CDK2low (CDK2-reporter); EZH2i -> shift to fast mode.

STATUS (prototype; infrastructure DONE, bifurcation calibration is open). The kernel + ensemble harness +
age-weighted count-fraction/Euler-Lotka readout + a fixed-point (nullcline) analyzer (`fixed_points`) are
built and run. Calibration findings (from the analyzer + ensemble sweeps), i.e. what a working fit needs:
  1. The CDK2-p27 toggle IS bistable (3 fixed points) for RbC ~ 5-11 with current params -- the structure
     exists. Cd OUT of the E2F drive (basin-probability setter via p27) is correct; Cd IN it collapses
     the MB period (timer bug).
  2. Birth RbC must sit INSIDE the bistable window (lowered Rb_ss to 9 so birth RbC~9). Good.
  3. BUT the separatrix is HIGH (CDK2act ~0.63 at RbC=9) vs the birth CDK2act (~0.1, because birth p27
     buffers CDK2 heavily), so ~no cell is born committed -> no fast mode. FIX: lower/sharpen the
     separatrix and increase its Cd-sensitivity (Kp down / ksp0 up) so (a) a real born-committed fraction
     exists and (b) it differs GNP vs MB. The separatrix barely moves with Cd now -> weak mitogen effect.
  4. The saddle-node (where G0 cells are FORCED to commit) is at RbC~3 = over-growth (mass~3) -> the ~25h
     deterministic G0. FIX: move the SN to RbC~6-7 so G0 cells commit within ~1 doubling.
  5. *** KEY STRUCTURAL FINDING ***: birth noise alone gives a BIMODAL outcome (born-committed-at-birth vs
     G0-until-SN), NOT the continuous right-skewed G0 TAIL the proposal wants. A continuous tail requires
     INTRINSIC noise (stochastic p27/CDK2 escape across the separatrix) = a Langevin/Gillespie SDE, not a
     deterministic ODE with only noisy births. That is a real additional build, not a parameter.
So: the mechanism is implementable and now mapped, but the fit is a multi-stage calibration (separatrix
placement + Cd-leverage + SN position + Rb-dilution G1) AND an SDE for the tail. v44 remains committed.

Run:  ./venv/bin/python simulations/v45_stochastic_commitment.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from scipy.optimize import brentq

# ---------------------------------------------------------------------------
# The bistable commitment kernel (minutes). Reactions for species; assignments for the algebra.
# ---------------------------------------------------------------------------
KERNEL = """
model commitment_kernel
  species mass, Rb, CDK2, p27, Cdh1, Cb, Dna;
  mass = 1; Rb = 1.0; CDK2 = 0.04; p27 = 0.6; Cdh1 = 1; Cb = 0; Dna = 0;

  // ---- tunable parameters ----
  mu = 0.0009;                       // growth rate (sets the committed-cycler timescale via Rb dilution)
  ksRb = 0.0135; kdRb = 0.0015;      // Rb made size-INDEPENDENTLY -> [Rb]=Rb/mass dilutes (Rb_ss=9 -> birth
                                     // RbC~9 sits INSIDE the bistable window 5-11, so birth noise can scatter cells)
  KfireRb = 5.0; hRb = 8;            // origins fire only once Rb CONCENTRATION dilutes below this = the G1 timer
  Cd = 1.0;                          // mitogen (CyclinD1): 1 = GNP, ~7 = MB (sets p27 synthesis)
  wE = 2.6; d0 = 0.06; K_CdRb = 0.45; Km = 0.42; nE = 4;   // E2F release driven by CDK2 toggle (+ small basal d0)
                                                           // -- Cd is NOT in the drive: it's a basin-probability setter (via p27), not a timer
  ksE = 0.040; ksE0 = 0.0006; kdE = 0.013; Ki = 0.22;   // CDK2 (E2F-driven) ; p27 buffers CDK2
  ksp0 = 0.013; Kp = 1.2; kdp = 1.6; kdp0 = 0.0022;     // p27 synth (mitogen-suppressed) + CDK2(Skp2) clearance
  kon = 0.015; koff = 0.8;           // APC/C-Cdh1: CDK2act inactivates it (point of no return)
  ksCb = 0.024; kdCb = 0.010; KG2 = 0.85; nG2 = 6;      // CyclinB builds when Cdh1 off AND replication done
  kFire = 0.014; Kfire = 0.55; nFire = 6;               // origins fire when CDK2act crosses (G1->S); -> S ~3.6h
  // (G2 is short in this prototype -- a slow CyclinB build / Erlang chain is the refinement)

  // ---- algebraic ----
  RbC     := Rb/mass;                              // Rb CONCENTRATION (dilutes with growth)
  CDK2act := CDK2/(1 + p27/Ki);                    // p27 buffers/inhibits CDK2
  drive   := d0 + wE*CDK2act;                      // CDK2 (the toggle) drives Rb-P; Cd enters only via p27 (ksp)
  RbP     := drive/(K_CdRb*RbC + drive);           // fraction Rb inactivated (0..1); low [Rb] lowers the bar
  E2F     := RbP^nE/(Km^nE + RbP^nE);              // ultrasensitive E2F release
  ksp     := ksp0/(1 + Cd/Kp);                     // MITOGEN sets the brake: more Cd -> less p27 synth
  fire    := (CDK2act^nFire/(Kfire^nFire + CDK2act^nFire)) * (KfireRb^hRb/(KfireRb^hRb + RbC^hRb));
             // origins fire only when COMMITTED (CDK2act high) AND Rb diluted (RbC low) -> Rb-dilution G1 timer
  g2gate  := Dna^nG2/(KG2^nG2 + Dna^nG2);          // CyclinB only after replication ~complete

  // ---- dynamics ----
  Growth:  => mass; mu*mass;
  RbSyn:   => Rb; ksRb;
  RbDeg:   Rb => ; kdRb*Rb;
  CDK2syn: => CDK2; ksE*E2F + ksE0;
  CDK2deg: CDK2 => ; kdE*CDK2;
  p27syn:  => p27; ksp;
  p27deg:  p27 => ; (kdp*CDK2act + kdp0)*p27;        // CDK2 (Skp2) clears p27 = the toggle
  Cdh1on:  => Cdh1; kon*(1 - Cdh1);
  Cdh1off: Cdh1 => ; koff*CDK2act*Cdh1;
  DnaRep:  => Dna; kFire*fire*(1 - Dna);             // replicate once committed
  CbSyn:   => Cb; ksCb*(1 - Cdh1)*g2gate;
  CbDeg:   Cb => ; kdCb*Cb;
end
"""

Cb_DIV = 0.5                                   # CyclinB-CDK1 (MPF) threshold = mitosis/division
COMMIT_CDK2 = 0.30                             # CDK2act above this = committed (for G0/G1 classification)
SEL = ["time", "mass", "Rb", "CDK2", "p27", "Cdh1", "Cb", "Dna"]
_rr = te.loada(KERNEL)


def _cdk2act(CDK2, p27, Ki=0.22):
    return CDK2 / (1 + p27 / Ki)


def run_cell(p27_0, CDK2_0, Cd, mass_0=1.0, t_max=4000):
    """Simulate one cell from a (noisy) birth state to its first division. Returns trajectory dict
    + premitotic CDK2 (for daughter carryover), or None if it fails to divide in t_max."""
    for atol in (1e-9, 1e-8, 1e-7, 1e-6):
        _rr.reset()
        _rr['Cd'] = Cd; _rr['mass'] = mass_0
        _rr['p27'] = p27_0; _rr['CDK2'] = CDK2_0
        _rr['Cdh1'] = 1.0; _rr['Cb'] = 0.0; _rr['Dna'] = 0.0
        _rr['Rb'] = _rr['ksRb'] / _rr['kdRb'] * mass_0    # Rb at its size-scaled steady level
        _rr.integrator.setValue("absolute_tolerance", atol)
        _rr.integrator.setValue("relative_tolerance", 1e-7)
        try: _rr.integrator.setValue("maximum_num_steps", 200000)
        except Exception: pass
        try:
            r = _rr.simulate(0, t_max, t_max, selections=SEL)
        except Exception:
            continue
        Cb = r['Cb']
        idx = np.argmax(Cb > Cb_DIV)
        if Cb[idx] <= Cb_DIV:                 # never divided
            return None, np.nan
        sub = {k: r[k][:idx + 1] for k in SEL}
        return sub, float(r['CDK2'][idx])
    return None, np.nan


def phase_at(traj):
    """Per-time-point phase + key ages. 2N=Dna<0.05 (G0 if CDK2act<COMMIT else G1); S=0.05..0.95; G2M>=0.95."""
    Dna = traj['Dna']; cdk2a = _cdk2act(traj['CDK2'], traj['p27'])
    is_S = (Dna >= 0.05) & (Dna < 0.95)
    is_G2M = Dna >= 0.95
    is_2N = ~is_S & ~is_G2M
    is_G0 = is_2N & (cdk2a < COMMIT_CDK2)
    is_G1 = is_2N & (cdk2a >= COMMIT_CDK2)
    return dict(G0=is_G0, G1=is_G1, S=is_S, G2M=is_G2M)


def ensemble(Cd_base, N=200, sig_cd=0.0, P21_div=0.7, sig_p=0.5, phi=0.12, sig_c=0.4, seed=0):
    """Draw N noisy births, simulate each to division. Return periods, per-cell phase ages, G0 durations."""
    rng = np.random.default_rng(seed)
    cells = []
    for _ in range(N):
        Cd = Cd_base * (rng.lognormal(0, sig_cd) if sig_cd else 1.0)
        p27_0 = P21_div * rng.lognormal(0, sig_p)
        CDK2_0 = phi * rng.lognormal(0, sig_c)
        traj, _ = run_cell(p27_0, CDK2_0, Cd)
        if traj is None:
            continue
        t = traj['time']; per = t[-1] / 60.0          # hours
        ph = phase_at(traj)
        durs = {k: float(np.sum(v) * (t[1] - t[0]) / 60.0) for k, v in ph.items()}  # phase durations (h)
        g0_dur = durs['G0']
        cells.append(dict(per=per, durs=durs, g0=g0_dur, traj=traj, ph=ph))
    return cells


def population_growth_rate(periods):
    """Euler-Lotka for symmetric division: find lam s.t. <2 exp(-lam Tc)> = 1."""
    Tc = np.array(periods) * 60.0
    f = lambda lam: np.mean(2 * np.exp(-lam * Tc)) - 1.0
    try:
        return brentq(f, 1e-5, 1.0)
    except Exception:
        return np.log(2) / np.mean(Tc)


def count_fractions(cells):
    """Age-weighted COUNT-fractions across the population (the flow quantity)."""
    periods = [c['per'] for c in cells]
    lam = population_growth_rate(periods)                  # per minute
    acc = dict(G0=0.0, G1=0.0, S=0.0, G2M=0.0); W = 0.0
    for c in cells:
        t = c['traj']['time']; dt = t[1] - t[0]
        w = 2 * lam * np.exp(-lam * t)                     # age density within this cell's cycle
        for k, sel in c['ph'].items():
            acc[k] += float(np.sum(w * sel) * dt)
        W += float(np.sum(w) * dt)
    return {k: 100 * v / W for k, v in acc.items()}, lam


def fixed_points(RbC, Cd, p=None):
    """Nullcline analysis of the CDK2-p27 toggle at fixed (RbC, Cd): roots of g(x), x=CDK2act.
    3 roots = BISTABLE (low/separatrix/high). The calibration tool for placing the bifurcation."""
    if p is None:
        rr = te.loada(KERNEL)
        p = {k: rr[k] for k in ['wE', 'd0', 'K_CdRb', 'Km', 'nE', 'ksE', 'ksE0', 'kdE', 'Ki',
                                'ksp0', 'Kp', 'kdp', 'kdp0']}
    ksp = p['ksp0'] / (1 + Cd / p['Kp'])
    xs = np.linspace(1e-4, 3.0, 4000)
    p27 = ksp / (p['kdp'] * xs + p['kdp0'])
    cdk2_from_p27 = xs * (1 + p27 / p['Ki'])
    drive = p['d0'] + p['wE'] * xs
    RbP = drive / (p['K_CdRb'] * RbC + drive)
    E2F = RbP ** p['nE'] / (p['Km'] ** p['nE'] + RbP ** p['nE'])
    cdk2_from_E2F = (p['ksE'] * E2F + p['ksE0']) / p['kdE']
    g = cdk2_from_p27 - cdk2_from_E2F
    sc = np.where(np.diff(np.sign(g)) != 0)[0]
    return [0.5 * (xs[i] + xs[i + 1]) for i in sc]


if __name__ == "__main__":
    print("BISTABILITY (fixed points of the CDK2-p27 toggle vs RbC):")
    for RbC in [11, 9, 7, 5, 3]:
        print(f"  RbC={RbC:4.1f}  GNP {[round(x,2) for x in fixed_points(RbC,1.0)]}  "
              f"MB {[round(x,2) for x in fixed_points(RbC,7.0)]}")
    print()
    print("=" * 70)
    print("SINGLE-CELL sanity (deterministic, low-p27 birth = born committed):")
    for nm, Cd in [("GNP", 1.0), ("MB", 7.0)]:
        traj, _ = run_cell(p27_0=0.2, CDK2_0=0.15, Cd=Cd)
        if traj is None:
            print(f"  {nm}: did not divide"); continue
        ph = phase_at(traj); t = traj['time']; dt = (t[1] - t[0]) / 60
        d = {k: float(np.sum(v) * dt) for k, v in ph.items()}
        print(f"  {nm} (Cd={Cd}): period {t[-1]/60:.1f}h  G0 {d['G0']:.1f} G1 {d['G1']:.1f} "
              f"S {d['S']:.1f} G2M {d['G2M']:.1f}")
    print("\nENSEMBLE (noisy births):")
    for nm, Cd in [("GNP", 1.0), ("MB", 7.0)]:
        cells = ensemble(Cd, N=150)
        if not cells:
            print(f"  {nm}: no cells divided"); continue
        pers = np.array([c['per'] for c in cells]); g0s = np.array([c['g0'] for c in cells])
        committed = np.mean(g0s < 2.0)                     # born-committed = short G0
        cf, lam = count_fractions(cells)
        print(f"  {nm} (Cd={Cd}, n={len(cells)}/{150}): period mean {pers.mean():.1f}h median {np.median(pers):.1f}h "
              f"mode~{pers.min():.1f}h  born-committed {100*committed:.0f}%")
        print(f"       count-fractions G0/G1/S/G2M = {cf['G0']:.0f}/{cf['G1']:.0f}/{cf['S']:.0f}/{cf['G2M']:.0f}  "
              f"(2N={cf['G0']+cf['G1']:.0f}%; target ~68/16/16)")
