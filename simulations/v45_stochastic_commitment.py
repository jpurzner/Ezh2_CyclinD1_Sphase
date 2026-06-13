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

STATUS (prototype; MECHANISM DEMONSTRATED, quantitative fit ~half done). The kernel + ensemble + readout
+ fixed-point analyzer run, and the core claim is now PROVEN with the current params:
  * The mixture RESOLVES the period-vs-count-fraction over-constraint that broke v44. GNP: 51% born
    committed -> 2N count-fraction = 67% (target 68.2%) -- because fast committed cyclers supply ~45% 2N
    and the G0 tail lifts it to ~68% WITHOUT overshoot. This is the whole point, and it works.
  * The period distribution is RIGHT-SKEWED (fast mode ~11.6h + a G0 tail), from BIRTH NOISE ALONE --
    so the earlier "needs an SDE" worry was wrong: noisy p27/CDK2 partitioning across the separatrix
    gives variable G0 dwell. (An SDE would only smooth the tail.)
  * Mitogen effect emerges: MB 62% born committed / mean 13h vs GNP 51% / 15h. Cd acts via the p27
    CLEARANCE dynamics (MB clears p27 faster -> escapes G0 sooner), NOT via the separatrix position
    (which is on the E2F nullcline, Cd-independent).
Working regime: toggle bistable RbC~5-11; Rb_ss=9 (birth in window); KfireRb=3.5 (Rb-dilution G1 ~17h);
ksp0=0.04, Kp=0.4 (strong Cd leverage), Ki=0.12, kdp=2.5; birth phi=0.55, sig_p=0.55, P21_div=0.42.
G2/S REFINEMENT DONE (was REMAINING #1): G2p accumulator timer gives a real ~2.3h G2M (was ~0); a Dna>0.15
LATCH decouples S DURATION from the slow RbC-firing ramp -> clean ~3h S (was a throttled ~5h). Result:
GNP 2N/S/G2M = 74/17/9, MB = 70/20/10 (target 2N=68, S=16, G2M~2.5h DURATION not the unreliable 4N 16%).
v44-BIOLOGY WIRING DONE (was REMAINING #2): `condition_cd()` runs the v44 HH/MYCN/EZH2->CyclinD1 cascade
to steady state per condition and feeds the GNP-normalized Cd into v45's mitogen (two-timescale: the
cascade is quasi-static, and v44's steady EZH2 is ~condition-independent so the EZH2->Cd feedback doesn't
dynamically lengthen the period). `rescue_panel()` shows the paper's story EMERGES: HH withdrawal raises
the quiescent (G0) fraction (MB 29->47%, GNP 44->67%); EZH2i RESCUES (MB+HHi G0 47->35%); the effect is
larger in GNP than MB (MYCN floors MB's CyclinD1); period lengthens under HHi (the slow extra divisions).
ARREST MECHANISM DONE (was REMAINING #2): the CyclinD bootstrap d0:=dmax*Cd/(Kd0+Cd) makes Rb mono-
phospho MITOGEN-gated, and a mass cap (Mmax,nM) stops growth in the low basin so quiescent cells can't
dilute their way to commitment. Together they give a SHARP commitment threshold ~Cd 0.2-0.3: Cd>0.3 ->
all cycle; Cd<0.2 -> mostly PERMANENT arrest (only born-committed daughters divide). Emergent result:
GNP+HHi (Cd 0.16) -> 68% arrest, but MB+HHi (Cd 0.78, MYCN-floored ABOVE threshold) -> 0% permanent
arrest (only transient G0 48%) -- i.e. the MYCN floor PROTECTS MB from Hh-withdrawal arrest, the paper's
claim. EZH2i rescues (MB+HHi quiescent 48->31%). MB baseline 21% quiescent matches DMSO data (21.7%).
REMAINING (mechanical, not structural): (1) S a few % high (MB 22 vs 16) + MEAN period short of ~22h
(GNP 18, MB 16) -- lower mu lengthens the (growth-coupled) G1 while S/G2 (rate-fixed) shrink in fraction.
(2) Embed EZH2 as a v45 cycle species for the (small) within-cycle feedback + the Gli1-residual / EZH2
MB/GNP validation checks. v44 stays the committed working model.

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
  species mass, Rb, CDK2, p27, Cdh1, G2p, Dna;
  mass = 1; Rb = 1.0; CDK2 = 0.04; p27 = 0.6; Cdh1 = 1; G2p = 0; Dna = 0;

  // ---- tunable parameters ----
  mu = 0.0009;                       // growth rate (sets the committed-cycler timescale via Rb dilution)
  Mmax = 4.2; nM = 16;               // mass cap: cycling cells divide at mass ~3.5 (below cap), but a cell stuck in the
                                     // low basin grows to Mmax, floors RbC at 9/Mmax~2.1, stays bistable-low -> PERMANENT
                                     // arrest (quiescent cells don't grow). Without it, unbounded growth dilutes RbC->0
                                     // and forces commitment at ANY Cd>0 (low-Cd cells only slow, never arrest).
  ksRb = 0.0135; kdRb = 0.0015;      // Rb made size-INDEPENDENTLY -> [Rb]=Rb/mass dilutes (Rb_ss=9 -> birth
                                     // RbC~9 sits INSIDE the bistable window 5-11, so birth noise can scatter cells)
  KfireRb = 3.5; hRb = 8;            // origins fire only once Rb CONCENTRATION dilutes below this = the G1 timer
                                     // (3.5: cells must dilute Rb further -> committed G1 ~17h -> period ~22h, lowers S count-fraction)
  Cd = 1.0;                          // mitogen (CyclinD1): 1 = GNP, ~7 = MB (sets p27 synthesis)
  wE = 2.6; dmax = 0.12; Kd0 = 1.0; K_CdRb = 0.45; Km = 0.42; nE = 4;  // E2F release: CDK2 toggle + a SATURATING CyclinD bootstrap
                                                           // Cd enters ONLY the BASAL bootstrap d0 (CyclinD-CDK4/6 mono-phospho of Rb),
                                                           // NOT the wE*CDK2act toggle gain (that gain-coupling collapses the MB period)
  ksE = 0.040; ksE0 = 0.0015; kdE = 0.013; Ki = 0.12;   // CDK2 (E2F-driven, faster bootstrap) ; p27 buffers CDK2
  ksp0 = 0.040; Kp = 0.4; kdp = 2.5; kdp0 = 0.0022;     // p27 synth (STRONG mitogen suppression Kp=0.4) + faster CDK2(Skp2) clearance
  // NB: Cd does NOT move the separatrix (it's on the E2F nullcline); Cd acts via the p27-CLEARANCE DYNAMICS
  // -- MB (low ksp) clears p27 faster -> CDK2act rises across the fixed separatrix sooner -> shorter G0.
  kon = 0.015; koff = 0.8;           // APC/C-Cdh1: CDK2act inactivates it (point of no return)
  kG2 = 0.0067; KG2 = 0.85; nG2 = 6;                    // G2 timer: G2p accumulates while replicated+committed
                                                        // -> G2 duration = 1/kG2 ~ 2.5h (mitotic CyclinB-build delay, abstracted)
  kFire = 0.018; Kfire = 0.55; nFire = 6;               // origins fire when CDK2act crosses (G1->S); -> S ~3h

  // ---- algebraic ----
  RbC     := Rb/mass;                              // Rb CONCENTRATION (dilutes with growth)
  CDK2act := CDK2/(1 + p27/Ki);                    // p27 buffers/inhibits CDK2
  d0      := dmax*Cd/(Kd0 + Cd);                    // CyclinD-CDK4/6 mono-phospho = MITOGEN-GATED bootstrap, SATURATING:
                                                   // 0 at Cd=0 (no basal Rb-P -> dilution alone can't start E2F -> ARREST);
                                                   // ~0.06 at GNP Cd=1; saturates ~0.105 at MB Cd=7 (keeps MB's bistable birth -> retains G0)
  drive   := d0 + wE*CDK2act;                      // CyclinD bootstrap (d0) + CyclinE/CDK2 toggle (wE*CDK2act, the commitment switch)
  RbP     := drive/(K_CdRb*RbC + drive);           // fraction Rb inactivated (0..1); low [Rb] lowers the bar
  E2F     := RbP^nE/(Km^nE + RbP^nE);              // ultrasensitive E2F release
  ksp     := ksp0/(1 + Cd/Kp);                     // MITOGEN sets the brake: more Cd -> less p27 synth
  fire    := (CDK2act^nFire/(Kfire^nFire + CDK2act^nFire)) * (KfireRb^hRb/(KfireRb^hRb + RbC^hRb));
             // origins fire only when COMMITTED (CDK2act high) AND Rb diluted (RbC low) -> Rb-dilution G1 timer
  latch   := Dna^8/(0.15^8 + Dna^8);               // once replication has truly STARTED (Dna>0.15, above the pre-commit
                                                   // fire-leak ~0.03) it runs to completion -> decouples S DURATION from
                                                   // the slow RbC-firing ramp (clean ~3h S) WITHOUT bypassing the G1 timer
  g2gate  := Dna^nG2/(KG2^nG2 + Dna^nG2);          // CyclinB only after replication ~complete

  // ---- dynamics ----
  Growth:  => mass; mu*mass*(1 - (mass/Mmax)^nM);   // exponential growth with a soft cap at Mmax (quiescent cells stall)
  RbSyn:   => Rb; ksRb;
  RbDeg:   Rb => ; kdRb*Rb;
  CDK2syn: => CDK2; ksE*E2F + ksE0;
  CDK2deg: CDK2 => ; kdE*CDK2;
  p27syn:  => p27; ksp;
  p27deg:  p27 => ; (kdp*CDK2act + kdp0)*p27;        // CDK2 (Skp2) clears p27 = the toggle
  Cdh1on:  => Cdh1; kon*(1 - Cdh1);
  Cdh1off: Cdh1 => ; koff*CDK2act*Cdh1;
  DnaRep:  => Dna; kFire*(fire + latch - fire*latch)*(1 - Dna);   // fire INITIATES (G1 timer); latch COMPLETES (clean S)
  G2acc:   => G2p; kG2*(1 - Cdh1)*g2gate;               // G2 clock: runs only when committed (Cdh1 off) AND replication done
end
"""

G2_DIV = 1.0                                   # G2 clock threshold = mitosis/division
COMMIT_CDK2 = 0.30                             # CDK2act above this = committed (for G0/G1 classification)
SEL = ["time", "mass", "Rb", "CDK2", "p27", "Cdh1", "G2p", "Dna"]
_rr = te.loada(KERNEL)


def _cdk2act(CDK2, p27, Ki=0.22):
    return CDK2 / (1 + p27 / Ki)


def run_cell(p27_0, CDK2_0, Cd, mass_0=1.0, t_max=4000):
    """Simulate one cell from a (noisy) birth state to its first division. Returns
    (trajectory dict, premitotic CDK2, status) with status in {'divided','arrested','failed'}.
    'arrested' = integrated cleanly but never divided in t_max (mass-capped, stuck in the low basin)."""
    for atol in (1e-9, 1e-8, 1e-7, 1e-6):
        _rr.reset()
        _rr['Cd'] = Cd; _rr['mass'] = mass_0
        _rr['p27'] = p27_0; _rr['CDK2'] = CDK2_0
        _rr['Cdh1'] = 1.0; _rr['G2p'] = 0.0; _rr['Dna'] = 0.0
        _rr['Rb'] = _rr['ksRb'] / _rr['kdRb'] * mass_0    # Rb at its size-scaled steady level
        _rr.integrator.setValue("absolute_tolerance", atol)
        _rr.integrator.setValue("relative_tolerance", 1e-7)
        try: _rr.integrator.setValue("maximum_num_steps", 200000)
        except Exception: pass
        try:
            r = _rr.simulate(0, t_max, t_max, selections=SEL)
        except Exception:
            continue
        G2p = r['G2p']
        idx = np.argmax(G2p > G2_DIV)
        if G2p[idx] <= G2_DIV:                # integrated but never divided -> permanent arrest (mass cap + low basin)
            return None, np.nan, 'arrested'
        sub = {k: r[k][:idx + 1] for k in SEL}
        return sub, float(r['CDK2'][idx]), 'divided'
    return None, np.nan, 'failed'             # integrator gave up (rare); not counted as a biological arrest


def phase_at(traj):
    """Per-time-point phase + key ages. 2N=Dna<0.05 (G0 if CDK2act<COMMIT else G1); S=0.05..0.95; G2M>=0.95."""
    Dna = traj['Dna']; cdk2a = _cdk2act(traj['CDK2'], traj['p27'])
    is_S = (Dna >= 0.05) & (Dna < 0.95)
    is_G2M = Dna >= 0.95
    is_2N = ~is_S & ~is_G2M
    is_G0 = is_2N & (cdk2a < COMMIT_CDK2)
    is_G1 = is_2N & (cdk2a >= COMMIT_CDK2)
    return dict(G0=is_G0, G1=is_G1, S=is_S, G2M=is_G2M)


def ensemble(Cd_base, N=200, sig_cd=0.0, P21_div=0.42, sig_p=0.55, phi=0.55, sig_c=0.5, seed=0,
             return_stats=False):
    """Draw N noisy births, simulate each to division. Returns the list of DIVIDED cells
    (or, if return_stats, also a dict with arrested/failed counts -> the population arrest fraction)."""
    rng = np.random.default_rng(seed)
    cells = []; n_arrested = 0; n_failed = 0
    for _ in range(N):
        Cd = Cd_base * (rng.lognormal(0, sig_cd) if sig_cd else 1.0)
        p27_0 = P21_div * rng.lognormal(0, sig_p)
        CDK2_0 = phi * rng.lognormal(0, sig_c)
        traj, _, status = run_cell(p27_0, CDK2_0, Cd)
        if traj is None:
            n_arrested += status == 'arrested'; n_failed += status == 'failed'
            continue
        t = traj['time']; per = t[-1] / 60.0          # hours
        ph = phase_at(traj)
        durs = {k: float(np.sum(v) * (t[1] - t[0]) / 60.0) for k, v in ph.items()}  # phase durations (h)
        g0_dur = durs['G0']
        cells.append(dict(per=per, durs=durs, g0=g0_dur, traj=traj, ph=ph))
    if return_stats:
        resolved = len(cells) + n_arrested                 # exclude integrator failures from the denominator
        arrest_frac = n_arrested / resolved if resolved else np.nan
        return cells, dict(N=N, divided=len(cells), arrested=n_arrested, failed=n_failed,
                           arrest_frac=arrest_frac)
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
        p = {k: rr[k] for k in ['wE', 'dmax', 'Kd0', 'K_CdRb', 'Km', 'nE', 'ksE', 'ksE0', 'kdE', 'Ki',
                                'ksp0', 'Kp', 'kdp', 'kdp0']}
    ksp = p['ksp0'] / (1 + Cd / p['Kp'])
    xs = np.linspace(1e-4, 3.0, 4000)
    p27 = ksp / (p['kdp'] * xs + p['kdp0'])
    cdk2_from_p27 = xs * (1 + p27 / p['Ki'])
    drive = p['dmax'] * Cd / (p['Kd0'] + Cd) + p['wE'] * xs   # saturating Cd-gated CyclinD bootstrap
    RbP = drive / (p['K_CdRb'] * RbC + drive)
    E2F = RbP ** p['nE'] / (p['Km'] ** p['nE'] + RbP ** p['nE'])
    cdk2_from_E2F = (p['ksE'] * E2F + p['ksE0']) / p['kdE']
    g = cdk2_from_p27 - cdk2_from_E2F
    sc = np.where(np.diff(np.sign(g)) != 0)[0]
    return [0.5 * (xs[i] + xs[i + 1]) for i in sc]


# ===========================================================================
# v44 BIOLOGY WIRING -- replace the bare `Cd` parameter with the v44 HH/MYCN/EZH2 -> CyclinD1
# cascade, so the v45 commitment toggle is driven by the REAL mitogen module and the paper's
# perturbations (HHi/vismo, EZH2i, Ptch1-loss = MB, MYCN amplification) flow through to commitment.
#
# Two-timescale separation (justified): the HH cascade is QUASI-STATIC vs the ~hours cell cycle
# (v44's own framing), and v44's steady EZH2 is nearly condition-independent (GNP 1.17 vs MB 1.21
# -> the EZH2->CyclinD1 feedback does NOT dynamically lengthen the period, per the v44 calibration
# outcome). So we run the cascade to steady state per condition, read its CyclinD1 (Cd), and feed
# the GNP-NORMALIZED value (GNP=1) into v45's mitogen scale (where Kp=0.4 was tuned for Cd~1..7).
# This reproduces the rescue ORDERING from the actual biology; the (small) within-cycle EZH2
# feedback is a later refinement once EZH2 is embedded in the v45 cycle.
# ===========================================================================
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

# Condition knobs map to the v44 HH-cascade boundary inputs (mirrors simulations/validate_v44.py).
PTCH1_MB = 0.1           # MB = Ptch1 loss (constitutive Hedgehog)
MYCN_AMP_MB = 2.8        # MB MYCN amplification
CONDITIONS = {           # SHH, Ptch1_copy_number, GDC0449(=HHi/vismo), EZH2i, MYCN_amplification
    'GNP':              dict(SHH=0.5, Ptch1_copy_number=1.0,      GDC0449=0, EZH2i=0, MYCN_amplification=1.0),
    'GNP+HHi':          dict(SHH=0.5, Ptch1_copy_number=1.0,      GDC0449=1, EZH2i=0, MYCN_amplification=1.0),
    'GNP+EZH2i':        dict(SHH=0.5, Ptch1_copy_number=1.0,      GDC0449=0, EZH2i=1, MYCN_amplification=1.0),
    'MB':               dict(SHH=0.5, Ptch1_copy_number=PTCH1_MB, GDC0449=0, EZH2i=0, MYCN_amplification=MYCN_AMP_MB),
    'MB+HHi':           dict(SHH=0.5, Ptch1_copy_number=PTCH1_MB, GDC0449=1, EZH2i=0, MYCN_amplification=MYCN_AMP_MB),
    'MB+EZH2i':         dict(SHH=0.5, Ptch1_copy_number=PTCH1_MB, GDC0449=0, EZH2i=1, MYCN_amplification=MYCN_AMP_MB),
    'MB+HHi+EZH2i':     dict(SHH=0.5, Ptch1_copy_number=PTCH1_MB, GDC0449=1, EZH2i=1, MYCN_amplification=MYCN_AMP_MB),
}

_V44_RR = None
_CD_CACHE = {}


def _v44_cascade():
    global _V44_RR
    if _V44_RR is None:
        from build_model_v44_heldt import build_model_v44
        _V44_RR = te.loada(build_model_v44(with_growth=True, with_hh=True, with_two_step_rb=True))
    return _V44_RR


def condition_cd(name, normalize=True):
    """Steady-state CyclinD1 (Cd) from the v44 HH/MYCN/EZH2 cascade for a named condition.
    normalize=True returns it on the v45 mitogen scale (GNP=1)."""
    if name in _CD_CACHE:
        return _CD_CACHE[name]
    rr = _v44_cascade(); rr.reset()
    for k, v in CONDITIONS[name].items():
        rr[k] = v
    try:
        rr.simulate(0, 3000, 1500)
    except Exception:
        pass
    cd = float(rr['Cd'])
    _CD_CACHE[name] = cd
    if normalize:
        if 'GNP' not in _CD_CACHE:
            condition_cd('GNP', normalize=False)
        return cd / _CD_CACHE['GNP']
    return cd


def rescue_panel(N=300):
    """Drive the v45 ensemble with the v44-derived Cd for each condition.
    arrest% = newborns that PERMANENTLY exit the cycle (mass-capped, low basin) = the pRb-/Ki67-
    fraction analog. Quiescent% ~ arrest + (1-arrest)*cyclerG0 (snapshot estimate; ignores growth-
    dilution of the arrested pool). EZH2i should LOWER both arrest% and Quiescent%."""
    print("RESCUE PANEL (v45 commitment driven by the v44 HH/MYCN/EZH2 -> CyclinD1 cascade):")
    print(f"  {'condition':16s} {'Cd(norm)':>9s} {'arrest%':>8s} {'cyclerG0%':>10s} {'Quiesc%':>8s} "
          f"{'meanT(h)':>9s} {'2N/S/G2M':>10s}")
    for name in CONDITIONS:
        cd = condition_cd(name)
        cells, st = ensemble(cd, N=N, return_stats=True)
        arr = 100 * st['arrest_frac']
        if not cells:
            print(f"  {name:16s} {cd:9.2f} {arr:7.0f}% {'--':>10s} {arr:7.0f}%  (all arrest)")
            continue
        pers = np.array([c['per'] for c in cells])
        cf, _ = count_fractions(cells)
        quiesc = arr + (1 - st['arrest_frac']) * cf['G0']
        print(f"  {name:16s} {cd:9.2f} {arr:7.0f}% {cf['G0']:9.0f}% {quiesc:7.0f}% "
              f"{pers.mean():8.1f}  {cf['G0']+cf['G1']:.0f}/{cf['S']:.0f}/{cf['G2M']:.0f}")


if __name__ == "__main__":
    print("BISTABILITY (fixed points of the CDK2-p27 toggle vs RbC):")
    for RbC in [11, 9, 7, 5, 3]:
        print(f"  RbC={RbC:4.1f}  GNP {[round(x,2) for x in fixed_points(RbC,1.0)]}  "
              f"MB {[round(x,2) for x in fixed_points(RbC,7.0)]}")
    print()
    print("=" * 70)
    print("SINGLE-CELL sanity (deterministic, low-p27 birth = born committed):")
    for nm, Cd in [("GNP", 1.0), ("MB", 7.0)]:
        traj, _, _ = run_cell(p27_0=0.2, CDK2_0=0.15, Cd=Cd)
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
    print("\n" + "=" * 70)
    rescue_panel(N=200)
