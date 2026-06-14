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

STATUS: COMPLETE & VALIDATED (simulations/validate_v45.py -> 20/20), WITH a dynamic EZH2 species. v45 is
a self-consistent standalone model (parallel to v44, which stays the committed working model). EZH2 is a
DYNAMIC v45 species and a Rb-E2f TARGET, so its synthesis tracks the CyclinD1-CDK4/6-Rb-E2f axis:
MITOGEN-DOSE-dependent (saturating in CyclinD1: kEZsyn*(Cd/(Kcd+Cd))) AND cycle-gated (ezgate, ~0 in G0 /
under CDK4/6i, ~1 in S/G2). STABLE -> integrates over the cycle (longer S -> more EZH2); HALVED at
division. EZH2 REPRESSES CyclinD1: Cd := Cd_drive*K/(K+EZH2*(1-EZH2i)) (Cd_drive = the EZH2-independent
Gli/MYCN drive from the v44 cascade, condition_drive) -> the CyclinD1<->EZH2 NEGATIVE FEEDBACK (Fig 4K)
is INTERNAL to v45. Calibrated to the paper folds: EZH2 MB/GNP ~2.1x (Fig4J 2.05x), CyclinD1 MB/GNP ~4.4x
(Fig4I 5.07x -- the feedback pulls MB CyclinD1 DOWN from the cascade's 7.6x toward the measured fold),
EZH2i de-repression ~2.6x (Fig3C 2.2x), G2/G0 ~1.7x (Fig4A/B). The EZH2i rescue is EMERGENT (remove v45's
EZH2 -> Cd_eff jumps to Cd_drive). The kernel + ensemble + readout + fixed-point analyzer run, and the
core claim is PROVEN with the current params:
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
SOFT FIT TUNED: mu 0.0009->0.0008 + kFire 0.018->0.022 -> GNP single-cell period 22.8h, ensemble MEDIAN
22.6h (on the ~22-23h literature value; the mean 19h is lower due to the fast committed mode + G0 tail).
MB phase-fractions 70/19/11 (2N target 68, S 16, G2M~2.5h DURATION) -- closer to target than v44's own
fit (v44 over-counts 2N~80, under-counts S~10). MB baseline 24% quiescent ~ DMSO data 21.7%.
REMAINING: Embed EZH2 as a v45 cycle species for the (small) within-cycle feedback + the Gli1-residual /
EZH2 MB/GNP validation checks; write validate_v45.py. v44 stays the committed working model.

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
  species mass, Rb, CDK2, p27, Cdh1, G2p, Dna, EZH2;
  mass = 1; Rb = 1.0; CDK2 = 0.04; p27 = 0.6; Cdh1 = 1; G2p = 0; Dna = 0; EZH2 = 1.5;

  // ---- tunable parameters ----
  mu = 0.0008;                       // growth rate (sets the committed-cycler timescale via Rb dilution)
                                     // 0.0008 (was 0.0009): lengthens the growth-coupled G1 -> GNP mean period ~20h (toward
                                     // the ~22h literature value) while the rate-fixed S/G2 shrink as count-fractions
  Mmax = 4.2; nM = 16;               // mass cap: cycling cells divide at mass ~3.5 (below cap), but a cell stuck in the
                                     // low basin grows to Mmax, floors RbC at 9/Mmax~2.1, stays bistable-low -> PERMANENT
                                     // arrest (quiescent cells don't grow). Without it, unbounded growth dilutes RbC->0
                                     // and forces commitment at ANY Cd>0 (low-Cd cells only slow, never arrest).
  ksRb = 0.0135; kdRb = 0.0015;      // Rb made size-INDEPENDENTLY -> [Rb]=Rb/mass dilutes (Rb_ss=9 -> birth
                                     // RbC~9 sits INSIDE the bistable window 5-11, so birth noise can scatter cells)
  KfireRb = 3.5; hRb = 8;            // origins fire only once Rb CONCENTRATION dilutes below this = the G1 timer
                                     // (3.5: cells must dilute Rb further -> committed G1 ~17h -> period ~22h, lowers S count-fraction)
  // ---- EZH2 -> CyclinD1 feedback (the paper's mechanism) ----
  Cd_drive = 2.56;                   // EZH2-INDEPENDENT CyclinD1 transcription drive (Gli/MYCN), set per condition
                                     // from the v44 cascade's UN-repressed Cd (GNP drive = GNP+EZH2i/GNP = 2.56).
  EZH2i = 0;                         // EZH2->CyclinD1 repression toggle (1 = EZH2 inhibitor -> removes the brake -> de-represses Cd)
  K_EZH2_Cd = 1.0;                   // EZH2 repression half-max; calibrated (with kEZsyn/Kcd) to EZH2 MB/GNP ~2.1x (Fig4J 2.05x),
                                     // CyclinD1 MB/GNP ~4.4x (Fig4I 5.07x), EZH2i de-repression ~2.6x (Fig3C 2.2x), GNP Cd~1.0 (period ~21h)
  kEZbas = 0.00020; kEZsyn = 0.0045; kDeEZ = 0.00010; Kez = 0.35; Kcd = 2.0;  // EZH2 is a Rb-E2f target -> synthesis tracks
                                     // the CyclinD1-CDK4/6-Rb-E2f axis: MITOGEN-DOSE-dependent (proportional to CyclinD1 = Cd) AND
                                     // cycle-gated (ezgate, high in S/G2, ~0 in G0 / under CDK4/6i). Fig 4A/B/H-J. STABLE (small
                                     // kDeEZ) -> INTEGRATES over the cycle (longer S -> more EZH2, the HU result); HALVED at division.
  wE = 2.6; dmax = 0.10; Kd0 = 0.82; nd0 = 2; K_CdRb = 0.45; Km = 0.42; nE = 4;  // E2F release: CDK2 toggle + a SATURATING CyclinD bootstrap
                                                           // Cd enters ONLY the BASAL bootstrap d0 (CyclinD-CDK4/6 mono-phospho of Rb),
                                                           // NOT the wE*CDK2act toggle gain (that gain-coupling collapses the MB period)
  ksE = 0.040; ksE0 = 0.0015; kdE = 0.013; Ki = 0.12;   // CDK2 (E2F-driven, faster bootstrap) ; p27 buffers CDK2
  ksp0 = 0.040; Kp = 0.4; kdp = 2.5; kdp0 = 0.0022;     // p27 synth (STRONG mitogen suppression Kp=0.4) + faster CDK2(Skp2) clearance
  // NB: Cd does NOT move the separatrix (it's on the E2F nullcline); Cd acts via the p27-CLEARANCE DYNAMICS
  // -- MB (low ksp) clears p27 faster -> CDK2act rises across the fixed separatrix sooner -> shorter G0.
  kon = 0.015; koff = 0.8;           // APC/C-Cdh1: CDK2act inactivates it (point of no return)
  kG2 = 0.0067; KG2 = 0.85; nG2 = 6;                    // G2 timer: G2p accumulates while replicated+committed
                                                        // -> G2 duration = 1/kG2 ~ 2.5h (mitotic CyclinB-build delay, abstracted)
  kFire = 0.022; Kfire = 0.55; nFire = 6;               // origins fire when CDK2act crosses (G1->S); -> S ~2.5h
                                                        // (period-scaled BrdU anchor: at Tc~16-18h, Ts/Tc~15% -> Ts~2.5h, not 3h)

  // ---- algebraic ----
  RbC     := Rb/mass;                              // Rb CONCENTRATION (dilutes with growth)
  CDK2act := CDK2/(1 + p27/Ki);                    // p27 buffers/inhibits CDK2
  Cd      := Cd_drive*K_EZH2_Cd/(K_EZH2_Cd + EZH2*(1 - EZH2i));   // EZH2 REPRESSES CyclinD1 (the paper's feedback);
                                                   // EZH2i=1 removes it -> Cd jumps to Cd_drive (the de-repression / rescue)
  d0      := dmax*Cd^nd0/(Kd0^nd0 + Cd^nd0);        // CyclinD-CDK4/6 mono-phospho = MITOGEN-GATED bootstrap, SATURATING (Hill-2):
                                                   // ~0 at Cd<0.3 (no basal Rb-P -> dilution alone can't start E2F -> ARREST, incl. GNP-SHH Cd~0.25);
                                                   // ~0.06 at GNP Cd=1; saturates ~0.10 at MB Cd=7 (keeps MB's bistable birth -> retains G0).
                                                   // Hill-2 (vs Hill-1) sharpens the low-Cd cutoff so Cd~0.25 arrests while Cd~0.78 (MYCN-floored MB+HHi) still cycles
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
  ezgate  := CDK2act^2/(Kez^2 + CDK2act^2);        // E2f/commitment gate for EZH2 synthesis: ~0 in G0 (and under CDK4/6i,
                                                   // which blocks commitment -> abolishes cycle-driven EZH2, Fig 4F), ~1 once committed (S/G2)

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
  EZH2syn: => EZH2; kEZbas + kEZsyn*(Cd/(Kcd + Cd))*ezgate;   // MITOGEN-DOSE (SATURATING in Cd: high-mitogen MB doesn't over-drive
                                                       // EZH2) x CYCLE-GATE (ezgate): Rb-E2f-driven, CyclinD1-dependent EZH2 -> tracks Hh/CyclinD1
                                                       // dose AND peaks S/G2; the CyclinD1<->EZH2 NEGATIVE feedback (Fig 4K)
  EZH2deg: EZH2 => ; kDeEZ*EZH2;                        // STABLE -> EZH2 integrates synthesis over the cycle (halved at division -> dilution)
end
"""

G2_DIV = 1.0                                   # G2 clock threshold = mitosis/division
COMMIT_CDK2 = 0.30                             # CDK2act above this = committed (for G0/G1 classification)
SEL = ["time", "mass", "Rb", "CDK2", "p27", "Cdh1", "G2p", "Dna", "EZH2", "Cd"]
_rr = te.loada(KERNEL)


def _cdk2act(CDK2, p27, Ki=0.22):
    return CDK2 / (1 + p27 / Ki)


def run_cell(p27_0, CDK2_0, Cd_drive, EZH2_0=1.5, EZH2i=0, mass_0=1.0, t_max=4000):
    """Simulate one cell from a (noisy) birth state to its first division. Cd_drive = the EZH2-independent
    CyclinD1 drive (Gli/MYCN); EZH2_0 = inherited EZH2 at birth (the stable, diluted-at-division integrator).
    Returns (trajectory dict, premitotic CDK2, premitotic EZH2, status) with status in
    {'divided','arrested','failed'}. 'arrested' = integrated cleanly but never divided in t_max."""
    for atol in (1e-9, 1e-8, 1e-7, 1e-6):
        _rr.reset()
        _rr['Cd_drive'] = Cd_drive; _rr['EZH2i'] = EZH2i; _rr['mass'] = mass_0
        _rr['p27'] = p27_0; _rr['CDK2'] = CDK2_0; _rr['EZH2'] = EZH2_0
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
            return None, np.nan, np.nan, 'arrested'
        sub = {k: r[k][:idx + 1] for k in SEL}
        return sub, float(r['CDK2'][idx]), float(r['EZH2'][idx]), 'divided'
    return None, np.nan, np.nan, 'failed'     # integrator gave up (rare); not counted as a biological arrest


def phase_at(traj):
    """Per-time-point phase + key ages. 2N=Dna<0.05 (G0 if CDK2act<COMMIT else G1); S=0.05..0.95; G2M>=0.95."""
    Dna = traj['Dna']; cdk2a = _cdk2act(traj['CDK2'], traj['p27'])
    is_S = (Dna >= 0.05) & (Dna < 0.95)
    is_G2M = Dna >= 0.95
    is_2N = ~is_S & ~is_G2M
    is_G0 = is_2N & (cdk2a < COMMIT_CDK2)
    is_G1 = is_2N & (cdk2a >= COMMIT_CDK2)
    return dict(G0=is_G0, G1=is_G1, S=is_S, G2M=is_G2M)


def equilibrate_ezh2(Cd_drive, EZH2i=0, iters=6):
    """EZH2 is inherited (halved) across divisions, so its birth value is a fixed point.
    Iterate a committed cell deterministically: EZH2_birth -> run -> EZH2_div -> birth = EZH2_div/2.
    Returns the self-consistent birth EZH2 for a CYCLING cell in this condition."""
    E = 1.5
    for _ in range(iters):
        traj, _, E_div, status = run_cell(0.20, 0.55, Cd_drive, EZH2_0=E, EZH2i=EZH2i)
        if status != 'divided':
            return E                                       # non-cycling here; birth value moot
        E = E_div / 2.0
    return E


def ensemble(Cd_drive, EZH2i=0, N=200, P21_div=0.42, sig_p=0.55, phi=0.55, sig_c=0.5, sig_ez=0.25,
             seed=0, return_stats=False):
    """Draw N noisy births, simulate each to division. Cd_drive = EZH2-independent CyclinD1 drive;
    EZH2 birth value is equilibrated (self-consistent halving) then drawn with lognormal noise.
    Returns the list of DIVIDED cells (+ a stats dict with arrest fraction if return_stats)."""
    rng = np.random.default_rng(seed)
    E_eq = equilibrate_ezh2(Cd_drive, EZH2i)               # steady birth EZH2 for a cycling cell
    cells = []; n_arrested = 0; n_failed = 0
    for _ in range(N):
        p27_0 = P21_div * rng.lognormal(0, sig_p)
        CDK2_0 = phi * rng.lognormal(0, sig_c)
        EZH2_0 = E_eq * rng.lognormal(0, sig_ez)
        traj, _, _, status = run_cell(p27_0, CDK2_0, Cd_drive, EZH2_0=EZH2_0, EZH2i=EZH2i)
        if traj is None:
            n_arrested += status == 'arrested'; n_failed += status == 'failed'
            continue
        t = traj['time']; per = t[-1] / 60.0          # hours
        ph = phase_at(traj)
        durs = {k: float(np.sum(v) * (t[1] - t[0]) / 60.0) for k, v in ph.items()}  # phase durations (h)
        cells.append(dict(per=per, durs=durs, g0=durs['G0'], traj=traj, ph=ph))
    if return_stats:
        resolved = len(cells) + n_arrested                 # exclude integrator failures from the denominator
        arrest_frac = n_arrested / resolved if resolved else np.nan
        return cells, dict(N=N, divided=len(cells), arrested=n_arrested, failed=n_failed,
                           arrest_frac=arrest_frac, E_eq=E_eq)
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
        p = {k: rr[k] for k in ['wE', 'dmax', 'Kd0', 'nd0', 'K_CdRb', 'Km', 'nE', 'ksE', 'ksE0', 'kdE', 'Ki',
                                'ksp0', 'Kp', 'kdp', 'kdp0']}
    ksp = p['ksp0'] / (1 + Cd / p['Kp'])
    xs = np.linspace(1e-4, 3.0, 4000)
    p27 = ksp / (p['kdp'] * xs + p['kdp0'])
    cdk2_from_p27 = xs * (1 + p27 / p['Ki'])
    d0 = p['dmax'] * Cd ** p['nd0'] / (p['Kd0'] ** p['nd0'] + Cd ** p['nd0'])   # saturating Cd-gated CyclinD bootstrap
    drive = d0 + p['wE'] * xs
    RbP = drive / (p['K_CdRb'] * RbC + drive)
    E2F = RbP ** p['nE'] / (p['Km'] ** p['nE'] + RbP ** p['nE'])
    cdk2_from_E2F = (p['ksE'] * E2F + p['ksE0']) / p['kdE']
    g = cdk2_from_p27 - cdk2_from_E2F
    sc = np.where(np.diff(np.sign(g)) != 0)[0]
    return [0.5 * (xs[i] + xs[i + 1]) for i in sc]


# ===========================================================================
# v44 BIOLOGY WIRING -- drive the v45 commitment toggle from the v44 HH/MYCN -> CyclinD1 cascade,
# so the paper's perturbations (HHi/vismo, Ptch1-loss = MB, MYCN amplification) flow through.
#
# EZH2 is NOT taken from the cascade -- it is a DYNAMIC v45 cycle species (see the kernel) that
# represses CyclinD1 inside the v45 cell cycle (the paper's feedback). From the cascade we take only
# the EZH2-INDEPENDENT drive `Cd_drive` (the Gli/MYCN transcription), obtained by running the cascade
# with EZH2i=1 (repression removed) and GNP-normalizing. v45's own EZH2 then re-applies the repression:
# Cd = Cd_drive * K/(K + EZH2*(1-EZH2i)). So the EZH2i RESCUE is emergent (remove v45's EZH2 -> Cd
# rises to Cd_drive), and EZH2's phase-dependence + MB/GNP elevation are v45 OUTPUTS. The HH cascade is
# quasi-static vs the ~hours cycle, so steady Cd_drive is principled.
# ===========================================================================
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))

# Condition knobs map to the v44 HH-cascade boundary inputs (mirrors simulations/validate_v44.py).
PTCH1_MB = 0.1           # MB = Ptch1 loss (constitutive Hedgehog)
MYCN_AMP_MB = 2.8        # MB MYCN amplification
CONDITIONS = {           # SHH, Ptch1_copy_number, GDC0449(=HHi/vismo), EZH2i, MYCN_amplification
    'GNP':              dict(SHH=0.5, Ptch1_copy_number=1.0,      GDC0449=0, EZH2i=0, MYCN_amplification=1.0),
    'GNP-SHH':          dict(SHH=0.0, Ptch1_copy_number=1.0,      GDC0449=0, EZH2i=0, MYCN_amplification=1.0),
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


def _cascade_cd(cond_params, force_ezh2i=None):
    """Run the v44 cascade to steady state and return its CyclinD1 (raw, un-normalized)."""
    rr = _v44_cascade(); rr.reset()
    for k, v in cond_params.items():
        rr[k] = v
    if force_ezh2i is not None:
        rr['EZH2i'] = force_ezh2i
    try:
        rr.simulate(0, 3000, 1500)
    except Exception:
        pass
    return float(rr['Cd'])


def condition_cd(name, normalize=True):
    """Steady-state REPRESSED CyclinD1 from the v44 cascade (reference value; GNP-normalized).
    This is the effective Cd v45 should land at the cycle level once its own EZH2 represses Cd_drive."""
    key = ('cd', name)
    if key not in _CD_CACHE:
        _CD_CACHE[key] = _cascade_cd(CONDITIONS[name])
    cd = _CD_CACHE[key]
    if normalize:
        if ('cd', 'GNP') not in _CD_CACHE:
            _CD_CACHE[('cd', 'GNP')] = _cascade_cd(CONDITIONS['GNP'])
        return cd / _CD_CACHE[('cd', 'GNP')]
    return cd


def condition_drive(name):
    """(Cd_drive, EZH2i) for a condition. Cd_drive = the EZH2-INDEPENDENT Gli/MYCN CyclinD1 drive
    (cascade with EZH2i=1, repression removed), GNP-normalized to the cascade's REPRESSED GNP Cd so
    the v45 scale matches (GNP drive ~2.56). EZH2i = the condition's flag (whether v45 re-represses)."""
    key = ('drive', name)
    if key not in _CD_CACHE:
        _CD_CACHE[key] = _cascade_cd(CONDITIONS[name], force_ezh2i=1)      # naked Gli/MYCN drive
    if ('cd', 'GNP') not in _CD_CACHE:
        _CD_CACHE[('cd', 'GNP')] = _cascade_cd(CONDITIONS['GNP'])
    drive = _CD_CACHE[key] / _CD_CACHE[('cd', 'GNP')]
    return drive, CONDITIONS[name]['EZH2i']


def ezh2_stats(cells):
    """EZH2 read-outs over cyclers: cycle-mean EZH2, the G2/G0 EZH2 ratio (Fig 4A/B ~2x),
    and the cycle-mean effective CyclinD1 (Cd = Cd_drive * EZH2-repression)."""
    allE, cdm, g0, g2 = [], [], [], []
    for c in cells:
        E = c['traj']['EZH2']; Cd = c['traj']['Cd']; ph = c['ph']
        allE.append(float(E.mean())); cdm.append(float(Cd.mean()))
        if ph['G0'].any():  g0.append(float(E[ph['G0']].mean()))
        if ph['G2M'].any(): g2.append(float(E[ph['G2M']].mean()))
    return dict(mean=float(np.mean(allE)) if allE else np.nan,
                cd_mean=float(np.mean(cdm)) if cdm else np.nan,
                g2_g0=(float(np.mean(g2)) / float(np.mean(g0))) if (g0 and g2) else np.nan)


def rescue_panel(N=300):
    """Drive the v45 ensemble from the v44 Gli/MYCN drive; EZH2 represses CyclinD1 INSIDE v45.
    arrest% = newborns that PERMANENTLY exit (mass-capped, low basin) = the pRb-/Ki67- analog.
    Quiescent% ~ arrest + (1-arrest)*cyclerG0. EZH2(rel) = cycle-mean EZH2 / GNP's."""
    print("RESCUE PANEL (v45; CyclinD1 = v44 Gli/MYCN drive repressed by v45's OWN dynamic EZH2):")
    print(f"  {'condition':16s} {'drive':>6s} {'Cd_eff':>7s} {'EZH2rel':>8s} {'arrest%':>8s} "
          f"{'Quiesc%':>8s} {'meanT':>6s} {'2N/S/G2M':>10s}")
    gnp_E = None
    for name in CONDITIONS:
        drive, ezh2i = condition_drive(name)
        cells, st = ensemble(drive, EZH2i=ezh2i, N=N, return_stats=True)
        arr = 100 * st['arrest_frac']
        if not cells:
            print(f"  {name:16s} {drive:6.2f} {'--':>7s} {'--':>8s} {arr:7.0f}%  (all arrest)")
            continue
        ez = ezh2_stats(cells)
        if name == 'GNP':
            gnp_E = ez['mean']
        erel = ez['mean'] / gnp_E if gnp_E else np.nan
        cf, _ = count_fractions(cells)
        quiesc = arr + (1 - st['arrest_frac']) * cf['G0']
        pers = np.array([c['per'] for c in cells])
        print(f"  {name:16s} {drive:6.2f} {ez['cd_mean']:7.2f} {erel:8.2f} {arr:7.0f}% "
              f"{quiesc:7.0f}% {pers.mean():6.1f}  {cf['G0']+cf['G1']:.0f}/{cf['S']:.0f}/{cf['G2M']:.0f}")


if __name__ == "__main__":
    print("SINGLE-CELL sanity (deterministic, low-p27 birth = born committed):")
    for nm in ("GNP", "MB"):
        drive, ezh2i = condition_drive(nm)
        E0 = equilibrate_ezh2(drive, ezh2i)
        traj, _, _, _ = run_cell(p27_0=0.2, CDK2_0=0.15, Cd_drive=drive, EZH2_0=E0, EZH2i=ezh2i)
        if traj is None:
            print(f"  {nm}: did not divide"); continue
        ph = phase_at(traj); t = traj['time']; dt = (t[1] - t[0]) / 60
        d = {k: float(np.sum(v) * dt) for k, v in ph.items()}
        E = traj['EZH2']; Cd = traj['Cd']
        print(f"  {nm} (drive={drive:.2f}): period {t[-1]/60:.1f}h  G0 {d['G0']:.1f} G1 {d['G1']:.1f} "
              f"S {d['S']:.1f} G2M {d['G2M']:.1f}  EZH2 {E[0]:.2f}->{E[-1]:.2f}  Cd_eff {Cd.mean():.2f}")
    print("\nENSEMBLE (noisy births):")
    for nm in ("GNP", "MB"):
        drive, ezh2i = condition_drive(nm)
        cells = ensemble(drive, EZH2i=ezh2i, N=150)
        if not cells:
            print(f"  {nm}: no cells divided"); continue
        pers = np.array([c['per'] for c in cells]); g0s = np.array([c['g0'] for c in cells])
        cf, lam = count_fractions(cells); ez = ezh2_stats(cells)
        print(f"  {nm} (drive={drive:.2f}, n={len(cells)}/150): period mean {pers.mean():.1f}h "
              f"median {np.median(pers):.1f}h  born-committed {100*np.mean(g0s < 2.0):.0f}%  "
              f"EZH2 mean {ez['mean']:.2f} (G2/G0 {ez['g2_g0']:.2f})  Cd_eff {ez['cd_mean']:.2f}")
        print(f"       count-fractions G0/G1/S/G2M = {cf['G0']:.0f}/{cf['G1']:.0f}/{cf['S']:.0f}/{cf['G2M']:.0f}  "
              f"(2N={cf['G0']+cf['G1']:.0f}%; target ~68/16/16)")
    print("\n" + "=" * 70)
    rescue_panel(N=200)
