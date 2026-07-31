"""Comprehensive validation of the v44 model (Heldt-based explicit-replication rebuild).

Single source of truth for calibration: runs every experimental condition, classifies
cell-cycle phases, counts divisions, measures steady levels, and scores ALL targets from
docs/Ezh2_CcnD1_model_targets.md. Mirrors the role of validate_v42.py for the old model.

Unlike v42 (Gerard-Goldbeter, eps-scaled "hours"), v44 is REAL TIME in MINUTES. Cyclin B is
`Cb := MPF + preMPF`; the cell cycle readouts are pRb / P21(=p27) / aRc / Dna / Skp2 / mass.

Phase classification (per the manuscript IF assay: phospho-Rb Ser807/811 + DAPI/PCNA):
  G0  : pre-replication, Rb HYPO-phosphorylated (pRb < thr)   -> quiescent-like / p27-high
  G1  : pre-replication, Rb phosphorylated (pRb >= thr)        -> committed, pre-S
  S   : origins fired & replicating (aRc high, Dna < 0.98)
  G2  : replication complete (Dna >= 0.98) until division
The pRb threshold is anchored so DMSO MB G0 matches the renormalized experimental value;
the HU redistribution and all other conditions are then predictions.

Edit PARAMS to iterate calibration, then:  ./venv/bin/python simulations/validate_v44.py
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

# ---------------------------------------------------------------------------
# Calibration parameter set (edit here; empty {} = builder defaults).
# ---------------------------------------------------------------------------
PARAMS = json.loads(os.environ.get('VALIDATE_PARAMS', '') or '{}')  # param overrides (wide-search driver); {} = builder defaults
                    # conditions in the full model -> needs the saturating CyclinD1->Rb drive fix.

P27_THR = 0.1        # p27 (P21) marker threshold (legacy; the G0/G1 split now uses pRb, below)
PRB_G0_THR = float(os.environ.get('PRB_G0_THR', '1.5'))   # G0 = pre-S cells with pRb(hyper) BELOW this = hypophospho/arrested
                     # (2026-07-15 transient-G0 reframe: pRb-hypophospho is the functional arrest switch, not p27 level)
MYCN_AMP_MB = 2.8
CD2_EXPR_MB = 3.59       # CyclinD2 MB-specific developmental elevation (added to the shared Cd2 basal); GNP = 0
PTCH1_MB = 0.1            # MB = Ptch1 loss (constitutive Hedgehog); v44 uses 0.1 (0 -> species->0)
P16_MB = 0.306           # MB CDK-inhibitor tones -- wide-search baked (were 0.15/1.5/0.004). p18 still the
P18_MB = float(os.environ.get('P18_MB', '1.73'))    # dominant INK4 (GNP 0.464 -> MB); CDKI data tone (was 1.553 legacy)
P19_MB = float(os.environ.get('P19_MB', '0.58'))    # Cdkn2d(p19) MB tone (GNP 0.36 builder default); used only when CDKI_SPECIES=1
KSYP21_MB = 0.002        #   (search co-tuned with the raised commitment threshold kPhRbCd=0.35 to preserve the rescue.)
P21_DIV_MB = float(os.environ.get('P21_DIV_MB', '0.6'))   # env-overridable. 2026-07-15 transient-G0 reframe: the
                         # MB-specific BIRTH-p27 "crutch" is RETIRED (set to the GNP baseline 0.6). MB's baseline
                         # reversible G0 must now EMERGE from real INK4 biology: the p27-clearance is gated on CDK4/6
                         # ACTIVITY (kPhRbCd*Cd), and MB's high p16/p18 brakes that activity -> slower clearance ->
                         # a modest baseline p27-high/pRb-low G0, GNP ~0. And CDK4/6i (kPhRbCd=0) fills the G0 in both.
                         # (docs/transient_g0_synthesis.md). NB the INK4 brake needs the joint re-fit to reach MB ~20%.
                         #   p27 (CIP/KIP) also brakes CDK4/6 via w_p27=1 (model default).
                         #   CIP/KIP = p21/p27 (kSyP21 2x the GNP baseline 0.002).
                         # GNP: p16=0, p18=0.4 (builder default), p21 baseline. Together they raise the
                         # CyclinD1 commitment threshold so vismo's CyclinD1 drop arrests MB; EZH2i (CyclinD1
                         # up) overcomes the COMPETITIVE INK4 brake -> rescue. CDK4/6i (kPhRbCd=0, Vmax) NOT rescuable.

SEL = ["time", "Cb", "MPF", "Cd", "Cd2", "Cd_mRNA", "MYCN", "Gli1", "EZH2", "EZH2m",
       "E2f", "pRb", "P21", "Skp2", "aRc", "Dna", "mass", "vfork", "Mk"]

# Build once; reset + set runtime inputs per condition (fast).
# H3K27_DILUTION=1 env var swaps in the replicative-dilution (leaky H3K27me3) repression module.
_MODEL = build_model_v44(with_ezh2=True, with_hh=True,
                         with_h3k27_dilution=(os.environ.get('H3K27_DILUTION', '1') == '1'),
                         with_h3k27_chain=(os.environ.get('H3K27_CHAIN', '1') == '1'),
                         with_two_step_rb=(os.environ.get('TWO_STEP_RB', '1') == '1'),
                         with_ezh2_conc=(os.environ.get('EZH2_CONC', '0') == '1'),   # 2026-07-19: DILUTION convention is now default (Cdh1-EZH2 re-bake); set EZH2_CONC=1 for the legacy concentration convention
                         with_mitogen_tracker=(os.environ.get('MITOGEN_TRACKER', '0') == '1'),
                         decouple_commit=(os.environ.get('DECOUPLE_COMMIT', '1') == '1'),  # 2026-07-27: BAKED default (cell-type split); G0/commit gated on CyclinD1/CDKi (not size). DECOUPLE_COMMIT=0 for the legacy size-gated commitment.
                         mycn_autoreg=(os.environ.get('MYCN_AUTOREG', '0') == '1'),  # 2026-07-27: MYCN from Gli1 + bistable self-activation (SHH-MB is NOT MYCN-amplified); HHi conditions use establish-then-withdraw
                         with_cdki_species=(os.environ.get('CDKI_SPECIES', '1') == '1'),  # 2026-07-30 BAKED default: individual dynamic CDKI species (set CDKI_SPECIES=0 for the legacy lumped-CDKI model)
                         params=PARAMS or None)

_MYCN_AUTOREG = os.environ.get('MYCN_AUTOREG', '0') == '1'


def _new_rr():
    rr = te.loada(_MODEL)
    rr.integrator.setValue("absolute_tolerance", 1e-9)
    rr.integrator.setValue("relative_tolerance", 1e-6)
    return rr


def run(shh=0.5, ptch1_cn=1.0, hhi=0.0, ezh2i=0.0, mycn_amp=1.0, p16=0.0, p18=None, p19=None, ksyp21=None, p21div=None,
        hu=0.0, cdk46i=False, serum_starve=False, t_end=12000, n_pts=48000):
    """Run one condition. Real-time minutes."""
    rr = _new_rr()
    last = None
    # tolerance + horizon + max-step retry: high-CyclinD1 MB conditions (esp. MB+HHi+EZH2i) and the
    # faster-EZH2 (kDeEZ) arrests can be stiff; capping the internal step cures CV_CONV/ILL_INPUT.
    for maxstep in (1e9, 20.0, 5.0):
      for te_end, te_pts in ((t_end, n_pts), (8000, 32000), (6000, 24000)):
        for atol in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5):
            rr.reset()
            try:
                rr.integrator.setValue("maximum_num_steps", 1000000)
                rr.integrator.setValue("maximum_time_step", maxstep)
            except Exception:
                pass
            rr['SHH'] = shh
            rr['Ptch1_copy_number'] = ptch1_cn
            rr['HHi'] = hhi
            rr['EZH2i'] = ezh2i
            # MYCN elevated expression (JP 2026-07-27): SHH-MB is NOT MYCN-amplified; MB's higher MYCN is a
            # developmentally-Gli-set elevated EXPRESSION, modeled as an additive MB-specific term (0 in GNP).
            # mycn_amp is kept as the interface (an expression FOLD); MYCN_expr = basal*(fold-1) is algebraically
            # equal to the old basal*amp multiplier, so all folds are unchanged.
            rr['MYCN_expr'] = rr['k_MYCN_synth_basal'] * (mycn_amp - 1.0)
            rr['Cd2_expr'] = CD2_EXPR_MB if ptch1_cn < 0.5 else 0.0   # CyclinD2 MB developmental elevation (MB has Ptch loss)
            rr['p16'] = p16              # INK4 (Cdkn2a): competitive CDK4/6 brake (0 in GNP, elevated in MB)
            if p18 is not None:
                rr['p18'] = p18          # INK4 (Cdkn2c): constitutive (GNP default 0.4), ~3x in MB
            _p19 = p19 if p19 is not None else (P19_MB if ptch1_cn < 0.5 else None)
            if _p19 is not None:
                try: rr['p19'] = _p19    # INK4 (Cdkn2d): CDKI-species mode; MB tone auto-applied (no-op if p19 param absent)
                except Exception: pass
            if ksyp21 is not None:
                rr['kSyP21'] = ksyp21    # p21/p27 (CDK2 inhibitor) synthesis; elevated in MB
            if p21div is not None:
                rr['P21_div'] = p21div   # MB-specific BIRTH p27 (fold-safe G0 lever; see P21_DIV_MB)
            rr['HU'] = hu
            if cdk46i:
                rr['kPhRbCd'] = 0.0          # block CycD-CDK4/6-mediated Rb phosphorylation
            if serum_starve:
                rr['k_Cd_translation'] = 0.0  # remove mitogen (CyclinD1) -> G0
            rr.integrator.setValue("absolute_tolerance", atol)
            try:
                if _MYCN_AUTOREG and hhi > 0:
                    # establish-then-withdraw: the MYCN bistable latch requires the cell to establish with
                    # Hh ON (MYCN latched high in MB) BEFORE vismodegib. A fresh HHi start can't reproduce
                    # cell-type-specific MYCN buffering (both types begin identical with Gli->0). So settle at
                    # HHi=0, then apply HHi and measure the withdrawn phase.
                    rr['HHi'] = 0.0
                    rr.simulate(0, 8000, 8000, selections=["time", "MYCN"])
                    rr['HHi'] = hhi
                    return rr.simulate(0, te_end, te_pts, selections=SEL)
                return rr.simulate(0, te_end, te_pts, selections=SEL)
            except Exception as e:
                last = e
    raise last


def count_divisions(res, settle=3000):
    """Divisions = MPF (mitotic) peaks after settle. Returns (n, peak_times_h, periods_h)."""
    t = res['time']; mpf = res['MPF']; m = t >= settle
    tt = t[m]
    dt = tt[1] - tt[0]
    peaks, _ = find_peaks(mpf[m], prominence=0.15, distance=int(200 / dt))  # >=200 min apart
    pk_t = tt[peaks] / 60.0
    per = np.diff(pk_t) if len(pk_t) > 1 else np.array([])
    return len(peaks), pk_t, per


def mean_last(res, sp, n=6000):
    return float(np.mean(res[sp][-n:]))


def mean_settled(res, sp, settle=4000):
    t = res['time']; return float(np.mean(res[sp][t >= settle]))


def classify(res, pRb_thr, settle=4000):
    """Cell-cycle phase fractions.  Returns (count_frac, ezph, dur_frac).

    *** count_frac is the COUNT-fraction = what flow cytometry actually measures ***, i.e. the
    fraction of cells *sitting* in each phase in an asynchronous, exponentially growing population.
    It is NOT the duration-fraction (phase time / period): in exp. growth there are ~2x as many
    just-divided cells as about-to-divide cells, so late phases hold fewer cells than their duration
    share.  We apply the standard age-density (lambda) correction: the normalized age density is
    n(a) = 2*lam*exp(-lam*a), a = time since the last division, lam = ln2/Tc.  Each settled time point
    is re-weighted by n(age); the weighted phase fraction is the count-fraction (= Tpot=lam*Ts/LI
    kinetics literature).  dur_frac (the old time-fraction) is kept for reference/figures only.

    G0 is now the pRb-HYPOPHOSPHORYLATION marker (2026-07-15 transient-G0 reframe): G0 = pre-S cells whose Rb is
    NOT hyperphosphorylated (pRb < PRB_G0_THR) -- the genuinely arrested/quiescent pre-commitment pool. This is the
    functional cycling/arrest switch (what CDK4/6i and the marker panel actually track), replacing the old
    p27-level marker (p27 IHC reports total protein -> not disjoint from cycling; see docs/transient_g0_synthesis.md).
    NOTE: flow CANNOT resolve G0 from G1 (both 2N) -- this G0/G1 split is a model-internal partition; the permanent
    differentiation exit (the larger non-cycling fraction) is a SEPARATE compartment, not this reversible G0."""
    t = res['time']; m = t >= settle; tt = t[m]; dt = tt[1] - tt[0]
    P21 = res['P21'][m]; aRc = res['aRc'][m]; Dna = res['Dna'][m]; MPF = res['MPF'][m]; pRbh = res['pRb'][m]
    # S-phase = ACTIVE DNA SYNTHESIS (BrdU/EdU analog: the synthesis flux vfork*aRc above a threshold),
    # NOT mere aRc presence. Under HU (slow forks) vfork ~ 0.1 so aRc-high but synthesis-low cells are
    # BrdU-NEGATIVE (the data measures BrdU) and belong in 2N/G1, not S. At HU=0 vfork=1 EXACTLY, so
    # vfork*aRc = aRc and DMSO is identical to the old aRc>0.05 gate (SYN_THR ~ that scale).
    SYN_THR = 0.02
    syn = res['vfork'][m] * aRc
    in_S = (syn > SYN_THR) & (Dna < 0.98)
    in_G2 = (Dna >= 0.98)
    preS = ~in_S & ~in_G2
    G0 = preS & (pRbh < PRB_G0_THR); G1 = preS & (pRbh >= PRB_G0_THR)   # G0 = pRb-hypophospho (arrested/pre-commit)
    sels = {"G0": G0, "G1": G1, "S": in_S, "G2": in_G2}
    # duration-fraction (time-fraction = single-cell phase-duration share)
    dur = {k: float(v.mean()) for k, v in sels.items()}; sd = sum(dur.values()) or 1.0
    dur = {k: 100 * v / sd for k, v in dur.items()}
    # COUNT-fraction (flow) via the age-density re-weighting n(a)=2*lam*exp(-lam*a)
    pk, _ = find_peaks(MPF, prominence=0.15, distance=int(200 / dt))
    if len(pk) >= 2:
        divt = tt[pk]; Tc = float(np.mean(np.diff(divt))); lam = np.log(2) / Tc
        idx = np.searchsorted(divt, tt, side="right") - 1            # last division before each point
        age = np.where(idx >= 0, tt - divt[np.clip(idx, 0, None)], np.nan)
        w = np.where(idx >= 0, 2 * lam * np.exp(-lam * np.clip(age, 0, None)), 0.0)
        W = np.sum(w) or 1.0
        cnt = {k: float(np.sum(w * v) / W) for k, v in sels.items()}; sc = sum(cnt.values()) or 1.0
        cnt = {k: 100 * v / sc for k, v in cnt.items()}
    else:
        cnt = dict(dur)                                              # not cycling -> no correction
    ez = res['EZH2'][m]
    def _phase_ez(ph, sel):
        if sel.any():
            return float(ez[sel].mean())
        if ph == 'G0':                       # commitment carryover: no distinct G0 -> pre-S EZH2 trough
            ref = ez[preS] if preS.any() else ez
            return float(np.mean(np.sort(ref)[:max(1, len(ref) // 4)])) if len(ref) else np.nan
        return np.nan
    ezph = {ph: _phase_ez(ph, sel) for ph, sel in sels.items()}
    return cnt, ezph, dur


# ---------------------------------------------------------------------------
# Run all conditions
# ---------------------------------------------------------------------------
CONDITIONS = {
    'GNP + SHH':         dict(shh=0.5, ptch1_cn=1.0, hhi=0.0, ezh2i=0.0, mycn_amp=1.0),
    'GNP + CDK4/6i':     dict(shh=0.5, ptch1_cn=1.0, hhi=0.0, ezh2i=0.0, mycn_amp=1.0, cdk46i=True),
    'GNP - SHH':         dict(shh=0.0, ptch1_cn=1.0, hhi=0.0, ezh2i=0.0, mycn_amp=1.0),
    'GNP + HHi':         dict(shh=0.5, ptch1_cn=1.0, hhi=1.0, ezh2i=0.0, mycn_amp=1.0),
    'GNP + EZH2i':       dict(shh=0.5, ptch1_cn=1.0, hhi=0.0, ezh2i=1.0, mycn_amp=1.0),
    'GNP + HHi + EZH2i': dict(shh=0.5, ptch1_cn=1.0, hhi=1.0, ezh2i=1.0, mycn_amp=1.0),
    'GNP Serum-starved': dict(shh=0.5, ptch1_cn=1.0, hhi=0.0, ezh2i=0.0, mycn_amp=1.0, serum_starve=True),
    'MB':                dict(shh=0.5, ptch1_cn=PTCH1_MB, hhi=0.0, ezh2i=0.0, mycn_amp=MYCN_AMP_MB, p16=P16_MB, p18=P18_MB, p19=P19_MB, ksyp21=KSYP21_MB, p21div=P21_DIV_MB),
    'MB + HHi':          dict(shh=0.5, ptch1_cn=PTCH1_MB, hhi=1.0, ezh2i=0.0, mycn_amp=MYCN_AMP_MB, p16=P16_MB, p18=P18_MB, ksyp21=KSYP21_MB, p21div=P21_DIV_MB),
    'MB + EZH2i':        dict(shh=0.5, ptch1_cn=PTCH1_MB, hhi=0.0, ezh2i=1.0, mycn_amp=MYCN_AMP_MB, p16=P16_MB, p18=P18_MB, ksyp21=KSYP21_MB, p21div=P21_DIV_MB),
    'MB + HHi + EZH2i':  dict(shh=0.5, ptch1_cn=PTCH1_MB, hhi=1.0, ezh2i=1.0, mycn_amp=MYCN_AMP_MB, p16=P16_MB, p18=P18_MB, ksyp21=KSYP21_MB, p21div=P21_DIV_MB),
    'MB + CDK4/6i':      dict(shh=0.5, ptch1_cn=PTCH1_MB, hhi=0.0, ezh2i=0.0, mycn_amp=MYCN_AMP_MB, p16=P16_MB, p18=P18_MB, ksyp21=KSYP21_MB, p21div=P21_DIV_MB, cdk46i=True),
    'MB + CDK4/6i+EZH2i':dict(shh=0.5, ptch1_cn=PTCH1_MB, hhi=0.0, ezh2i=1.0, mycn_amp=MYCN_AMP_MB, p16=P16_MB, p18=P18_MB, ksyp21=KSYP21_MB, p21div=P21_DIV_MB, cdk46i=True),
    'MB + HU':           dict(shh=0.5, ptch1_cn=PTCH1_MB, hhi=0.0, ezh2i=0.0, mycn_amp=MYCN_AMP_MB, p16=P16_MB, p18=P18_MB, ksyp21=KSYP21_MB, p21div=P21_DIV_MB, hu=1.0),
}


def main():
    print("=" * 92)
    print("v44 MODEL VALIDATION  (Heldt-based explicit-replication; real-time minutes)")
    print("PARAMS:", PARAMS if PARAMS else "(builder defaults)")
    print("=" * 92)

    sims = {name: run(**kw) for name, kw in CONDITIONS.items()}

    # G0 split by the p27 marker (principled; no anchoring needed)
    mb = sims['MB']
    pRb_thr = P27_THR

    # ---- summary table ----
    print(f"\n{'Condition':20s} | {'Div':>4} {'Period':>7} | "
          f"{'CycD1':>6} {'CdmRNA':>7} {'MYCN':>6} {'Gli1':>6} {'EZH2':>6} {'pRb':>6} {'p27':>6}")
    print("-" * 92)
    div = {}
    for name in CONDITIONS:
        r = sims[name]
        n, _, per = count_divisions(r)
        div[name] = n
        p = f"{np.mean(per):.1f}h" if len(per) else "arrest"
        print(f"{name:20s} | {n:4d} {p:>7s} | {mean_settled(r,'Cd'):6.3f} "
              f"{mean_settled(r,'Cd_mRNA'):7.3f} {mean_settled(r,'MYCN'):6.3f} "
              f"{mean_settled(r,'Gli1'):6.3f} {mean_settled(r,'EZH2'):6.3f} "
              f"{mean_settled(r,'pRb'):6.2f} {mean_settled(r,'P21'):6.3f}")

    # ---- target checks ----
    passed = [0]; total = [0]; rows = []

    def check(name, actual, target, tol=0.20, unit='', zero_ok=False):
        total[0] += 1
        if zero_ok:
            ok = abs(actual) <= 0.001
        elif target == 0:
            ok = abs(actual) <= 1
        else:
            ok = abs(actual - target) / abs(target) <= tol
        if ok:
            passed[0] += 1
        rows.append((ok, name, actual, target, unit))

    # transcript (Cd_mRNA matches RNA-seq); levels = mean settled
    cd = lambda c: mean_settled(sims[c], 'Cd_mRNA')
    cd2 = lambda c: mean_settled(sims[c], 'Cd2')   # CyclinD2 (separate, Hh-buffered D-cyclin)
    my = lambda c: mean_settled(sims[c], 'MYCN')
    gl = lambda c: mean_settled(sims[c], 'Gli1')
    ez = lambda c: mean_settled(sims[c], 'EZH2')
    ezm = lambda c: mean_settled(sims[c], 'EZH2m')

    # Section A — between-condition ratios
    check("CyclinD1 GNP+HHi/GNP", cd('GNP + HHi')/cd('GNP + SHH'), 0.157, 0.40)
    check("CyclinD1 MB+HHi/MB",   cd('MB + HHi')/cd('MB'),         0.144, 0.45)  # MB_HHi: 86% drop
    check("CyclinD1 MB/GNP",      cd('MB')/cd('GNP + SHH'),        5.07, 0.20)  # Fig 4I. Model 5.6 (11% high), TIGHT band -- the 2026-07-14 joint re-opt genuinely lands this (was 7.1 passing only on a wide 0.45 band; the old 7.1 was the inflated birth-p27=1.8 crutch, since removed)
    # CyclinD2 (Chahin RNA-seq): DOMINANT D-cyclin, Hh-BUFFERED (drops only ~35% vs CyclinD1 ~85%). Separate species.
    if cd2('GNP + SHH') > 1e-6:   # only when the two-cyclin model is engaged (w_Cd2/k_Cd2_bas > 0)
        check("CyclinD2 MB/GNP",       cd2('MB')/cd2('GNP + SHH'),      2.39, 0.25)
        check("CyclinD2 GNP+HHi/GNP",  cd2('GNP + HHi')/cd2('GNP + SHH'), 0.66, 0.30)
        check("CyclinD2 MB+HHi/MB",    cd2('MB + HHi')/cd2('MB'),       0.59, 0.45)  # MB drop non-sig (trend, padj 0.076) -> wide band
    check("MYCN GNP+HHi/GNP",     my('GNP + HHi')/my('GNP + SHH'), 0.78, 0.30)
    check("MYCN MB+HHi/MB",       my('MB + HHi')/my('MB'),         0.86, 0.25)
    check("MYCN MB/GNP",          my('MB')/my('GNP + SHH'),        2.80, 0.30)
    check("Gli1 GNP+HHi reduction", 1 - gl('GNP + HHi')/gl('GNP + SHH'), 0.99, 0.05)
    check("Gli1 MB/GNP",          gl('MB')/gl('GNP + SHH'),        6.90, 0.40)  # raw RNA-seq
    check("EZH2 MB/GNP",          ez('MB')/ez('GNP + SHH'),        2.05, 0.35)
    # Skp2 (JP scRNA timecourse + MB table): HIGHER in MB (~2.3x GNP-P7, mitogen-induced via the kSySkp2_Cd dose
    # term). The differentiation DECLINE (log2fc 1.73 down, t50 P12) is captured mechanistically by that same term
    # (Skp2 falls as mitogen/CyclinD1 drops); not checked here because serum-starve is a deeper arrest than mild
    # differentiation and the commitment needs Skp2~0 in the p27-high G0.
    sk = lambda c: mean_settled(sims[c], 'Skp2')
    check("Skp2 MB/GNP",          sk('MB')/sk('GNP + SHH'),        2.3, 0.30)
    # EZH2i de-repression of CyclinD1 (~2x)
    check("EZH2i CycD1 fold (GNP)", cd('GNP + EZH2i')/cd('GNP + SHH'), 2.2, 0.42)
    # EZH2 G0/cycling
    ez_cyc = ez('GNP + SHH'); ez_g0 = ez('GNP Serum-starved')
    check("EZH2 G0/cycling (~0.6)", ez_g0/ez_cyc if ez_cyc else 0, 0.6, 0.40)

    # proliferation-quiescence (divisions)
    check("GNP+SHH cycles (>0 div)", 1 if div['GNP + SHH'] > 0 else 0, 1, zero_ok=False, tol=0.01)
    check("GNP-SHH arrest (0 div)",  div['GNP - SHH'], 0, zero_ok=True)
    check("GNP+HHi arrest (0 div)",  div['GNP + HHi'], 0, zero_ok=True)
    check("GNP serum-starve arrest", div['GNP Serum-starved'], 0, zero_ok=True)
    check("MB cycles (>0 div)",      1 if div['MB'] > 0 else 0, 1, tol=0.01)
    # NOTE: with the data-matched HH, MB+HHi retains ~3x a cycling GNP's CyclinD1, so single cells
    # still cycle; vismodegib's effect (and the EZH2i rescue) is a POPULATION/fractional effect
    # (sim_ezh2i_population_dose.py / fig_v44_fig5_population.py), matching Fig 5D/E. The single-cell
    # arrest/rescue checks are therefore replaced by the CDK4/6i (kinase-block) contrast below.
    check("MB+CDK4/6i arrest (0)",   div['MB + CDK4/6i'], 0, zero_ok=True)
    check("MB+CDK4/6i+EZH2i no rescue", div['MB + CDK4/6i+EZH2i'], 0, zero_ok=True)

    # period
    _, _, per_gnp = count_divisions(sims['GNP + SHH'])
    check("Period GNP (~22h)", np.mean(per_gnp) if len(per_gnp) else 0, 22.0, 0.30)

    # ---- phase proportions: COUNT-fractions (what flow measures), via the lambda-correction ----
    # EVIDENCE / PROVENANCE of these targets (and why the comparison quantity changed):
    #  * DMSO G0/G1/S/G2 = 24.8/43.4/15.7/16.1 are FLOW-CYTOMETRY count-fractions (DNA-content gating):
    #    the fraction of cells SITTING in each phase. The model must therefore report the age-weighted
    #    COUNT-fraction (classify -> count_frac), NOT the duration-fraction (phase time / period).
    #    In exp. growth ~2x more just-divided than about-to-divide cells, so duration% over-counts late
    #    phases; n(a)=2*lam*exp(-lam*a), lam=ln2/Tc (the Tpot=lam*Ts/LI lambda-correction). [Comparing
    #    duration% to flow count% was a real bug -- it spuriously "passed" S and G2.]
    #  * Flow CANNOT resolve G0 vs G1 (both 2N). The 24.8% "G0" is really a p27-/Ki67- QUIESCENT
    #    fraction (a marker, not a per-cycle duration); for MB (growth fraction < 1) the 2N pool also
    #    holds truly-quiescent cells this single-cycle classifier misses (a growth-fraction term, not
    #    modeled). So the G0/G1 split is SOFT/model-internal -- we test the resolvable 2N pool (G0+G1).
    #  * S: reliable anchor is cumulative-BrdU (control GCP Tc-Ts ~20h -> Ts ~3h at Tc~23h), NOT the
    #    15.7% gate; S duration tracks the assumed period (true Tc~28h -> Ts ~8h). Checked vs flow (soft).
    #  * G2/M: the 16.1% 4N gate is an UNRELIABLE duration proxy (late-S near-4N content, doublets,
    #    tetraploidy). Direct methods (Fujita G2~2h + M~0.5h; pHH3/BrdU mitotic index peaking <2h) give
    #    G2+M ~2.5h. We anchor G2/M to that DIRECT DURATION, not the 4N count.
    f0, ez0, f0_dur = classify(mb, pRb_thr)
    f1, ez1, f1_dur = classify(sims['MB + HU'], pRb_thr)
    TGT_DMSO = dict(G0=24.8, G1=43.4, S=15.7, G2=16.1)              # flow count-fractions (soft for G0/G1, G2)
    _, _, per_mb = count_divisions(mb); Tc_mb = float(np.mean(per_mb)) if len(per_mb) else 23.0
    check("MB 2N (G0+G1) count% (flow)", f0['G0'] + f0['G1'], TGT_DMSO['G0'] + TGT_DMSO['G1'], 0.30, '%')
    check("MB S count% (flow; BrdU Ts~3h)", f0['S'], TGT_DMSO['S'], 0.42, '%')
    check("MB G2+M duration ~2.5h (direct)", (f0_dur['G2'] + f0_dur.get('M', 0)) / 100.0 * Tc_mb, 2.5, 0.55, 'h')
    TGT_HU_FOLD = dict(G0=1.19, G1=1.16, S=1.36, G2=0.23)
    for ph in ('S', 'G2'):  # the directionally clear ones (count-fraction folds)
        mf = f1[ph]/f0[ph] if f0[ph] else 0
        check(f"MB HU {ph} fold", mf, TGT_HU_FOLD[ph], 0.40)

    # EZH2 within-cycle gradient + HU boost (MB)
    check("EZH2 transcript S/G0 (1.8-2.5)", classify_grad(mb, pRb_thr, 'EZH2m')['S'], 2.0, 0.45)
    check("EZH2 protein G2/G0 (1.48)", ez0['G2']/ez0['G0'] if ez0['G0'] else 0, 1.48, 0.45)
    check("EZH2 Palbo mRNA drop (~0.44)", ezm('GNP + CDK4/6i')/ezm('GNP + SHH'), 0.44, 0.42)   # Fig 4: CDK4/6i drops EZH2 mRNA 56% (E2f-gated); tests the writer is cycle-gated (Point 2)
    check("HU EZH2-in-S boost (1.31)", ez1['S']/ez0['S'] if ez0['S'] else 0, 1.31, 0.35)

    # ---- print results ----
    print("\n" + "=" * 92)
    print("TARGET CHECKS")
    print("=" * 92)
    for ok, name, actual, target, unit in rows:
        sym = "[+]" if ok else "[-]"
        tstr = f"{target}{unit}" if not isinstance(target, float) else f"{target:g}{unit}"
        print(f"  {sym} {name:34s}: {actual:8.3f}  (target {tstr})")

    print(f"\nPhase COUNT-fractions (MB, flow-comparable):  DMSO {fmt(f0)}   HU {fmt(f1)}")
    print(f"  duration-fractions (single-cell, NOT flow):  DMSO {fmt(f0_dur)}")
    print(f"  flow targets G0/G1/S/G2 = 24.8/43.4/15.7/16.1 (G0/G1 unresolvable -> 2N=68.2; G2/M 4N soft, direct ~2.5h)")
    print(f"\n{'=' * 92}\nVALIDATION SUMMARY: {passed[0]}/{total[0]} targets passed\n{'=' * 92}")

    try:
        marks = {c: float(mean_settled(sims[c], 'Mk')) for c in sims}
    except Exception:
        marks = {}
    out = dict(model='v44_heldt', params=PARAMS, passed=passed[0], total=total[0],
               divisions=div, pRb_thr=float(pRb_thr), marks=marks,
               phase_dmso_count=f0, phase_dmso_duration=f0_dur, phase_hu_count=f1,
               checks=[{'name': n, 'actual': float(a), 'target': t, 'pass': bool(o)}
                       for o, n, a, t, u in rows])
    _rp = os.environ.get('VALIDATE_RESULTS') or os.path.join(os.path.dirname(__file__), 'validation_v44_results.json')
    with open(_rp, 'w') as fh:
        json.dump(out, fh, indent=2, default=float)
    return out


def classify_grad(res, pRb_thr, var, settle=4000):
    t = res['time']; m = t >= settle
    P21 = res['P21'][m]; aRc = res['aRc'][m]; Dna = res['Dna'][m]; v = res[var][m]
    in_S = (aRc > 0.05) & (Dna < 0.98); in_G2 = (Dna >= 0.98); preS = ~in_S & ~in_G2
    G0 = preS & (P21 > P27_THR)
    # G0 reference = the low-EZH2 quiescent phase (high-p27 pre-S). With commitment carryover
    # (f_commit_carry>0) committed cells lack a distinct high-p27 G0 -> fall back to the pre-S TROUGH
    # (lowest quartile of pre-S), which is the G0-equivalent low point. f_commit_carry=0 keeps v[G0].mean().
    if G0.sum() >= max(3, int(0.02 * len(v))):
        g0v = v[G0].mean()
    else:
        ref = v[preS] if preS.any() else v
        g0v = float(np.mean(np.sort(ref)[:max(1, len(ref) // 4)])) if len(ref) else np.nan
    e = lambda sel: (float(v[sel].mean())/g0v if sel.any() and g0v else np.nan)
    return dict(G1=e(preS & (P21 <= P27_THR)), S=e(in_S), G2=e(in_G2))


def fmt(f):
    return "/".join(f"{f[p]:.0f}" for p in ('G0', 'G1', 'S', 'G2'))


if __name__ == "__main__":
    main()
