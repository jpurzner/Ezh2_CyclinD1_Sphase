"""
Phase-3 calibration of the rewired EZH2 module (v43, EZH2 transcription -> Ma).

Why this exists
---------------
v42 gates EZH2 transcription on pRBpp, a CycE/CDK2 spike present only ~3% of the
cycle, so EZH2 is effectively flat and cannot show the measured within-cycle
gradient nor integrate S-phase duration. v43 (rewire_ezh2=True) drives EZH2
transcription off Ma = active CycA/CDK2 (the S/G2 marker) and speeds protein
turnover. This script:

  1. Re-checks the BETWEEN-condition EZH2 invariants v42 already satisfies
     (so the rewire does not regress the published validation):
         - mean cycling EZH2          ~0.49   (sets EZH2i->CycD1 ~2x, repression)
         - EZH2 G0/cycling ratio       ~0.60
         - GDC (HHi) reduces EZH2       25-47%
  2. Measures the NEW within-cycle gradient the rewire is meant to unlock:
         DMSO  G1/S/G2 vs G0 = 1.16 / 1.32 / 1.48   (Fig. 4 data)

Operational phase calls (data-free, from the cell-cycle state itself):
    G2/M : Cb (CycB) high                       -> 4N, mitotic
    S    : Ma (CycA/CDK2) high AND Cb low        -> replicating
    G1   : Ma low AND Cb low (post-division)     -> 2N, pre-S
    G0   : separate quiescent run (no SHH)       -> non-cycling reference

Run:  ./venv/bin/python simulations/calibrate_ezh2_rewire.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from src.build_model_v43_checkpoint import (
    build_model_v43, EZH2_REWIRE_DEFAULTS)
import tellurium as te


# EZH2 transcription params that live in the inherited v42 string and can be tuned
# at RUNTIME (no rebuild needed). basal sets the G0/G1 floor; E2F_amp sets the
# Ma-driven amplitude. Lowering basal drops quiescent G0 below the G1 trough so the
# within-cycle gradient (G0<G1<S<G2) emerges.
def run(shh=0.5, ptch1_cn=1.0, hhi=0.0, ezh2i=0.0, hu=0.0,
        t_end=200, n_pts=20000, rewire=True, ezh2_rewire_params=None,
        checkpoint=None, serum_starve=False, basal=None, e2f_amp=None):
    m = build_model_v43(checkpoint_params=checkpoint, rewire_ezh2=rewire,
                        ezh2_rewire_params=ezh2_rewire_params)
    rr = te.loada(m)
    rr['SHH'] = shh
    rr['Ptch1_copy_number'] = ptch1_cn
    rr['HHi'] = hhi
    rr['EZH2i'] = ezh2i
    rr['HU'] = hu
    if basal is not None:
        rr['k_EZH2_mRNA_synth_basal'] = basal
    if e2f_amp is not None:
        rr['k_EZH2_mRNA_synth_E2F'] = e2f_amp
    if shh == 0.0:
        rr['SHH_Ptch'] = 0.0
        rr['Ptch1_free'] = 1.0
        rr['Smo_active'] = 0.01
    if serum_starve:
        rr['k_Cd_translation'] = 0.0
    return rr.simulate(0, t_end, n_pts)


def phase_call(result, t_start=40):
    """Per-timepoint phase from Ma (CycA/CDK2, S/G2 marker) and Cb (CycB, M marker).

    Thresholds are fractions of each species' own peak (NOT distribution midpoints):
    the cyclins are sharply bimodal, so a midpoint sits at the spike top and collapses
    S. Fraction-of-max cleanly separates the S-onset (Ma~0.4 of peak, Cb still low)
    from the long low-Ma G1 and the Cb-high G2/M.
    """
    t = result['time']
    ma = result['[Ma]']
    cb = result['[Cb]']
    mask = t >= t_start
    thr_ma = 0.20 * ma[mask].max()    # S-onset Ma ~0.45 vs G1 Ma <0.13
    thr_cb = 0.25 * cb[mask].max()    # G2/M when CycB rises past ~1/4 peak
    phase = np.empty(len(t), dtype=object)
    for i in range(len(t)):
        if not mask[i]:
            phase[i] = None
        elif cb[i] > thr_cb:
            phase[i] = 'G2M'
        elif ma[i] > thr_ma:
            phase[i] = 'S'
        else:
            phase[i] = 'G1'
    return phase, dict(thr_ma=float(thr_ma), thr_cb=float(thr_cb))


def mean_ezh2_by_phase(result, t_start=40):
    phase, thr = phase_call(result, t_start)
    ez = result['[EZH2]']
    out = {}
    for ph in ('G1', 'S', 'G2M'):
        sel = np.array([p == ph for p in phase])
        out[ph] = float(np.mean(ez[sel])) if sel.sum() else np.nan
    return out, thr


def s_phase_duration(result, t_start=40):
    """Mean contiguous S-phase (Ma-high, Cb-low) duration per cycle."""
    t = result['time']
    phase, _ = phase_call(result, t_start)
    durs, in_s, start = [], False, None
    for i in range(len(t)):
        if phase[i] == 'S' and not in_s:
            in_s, start = True, t[i]
        elif phase[i] != 'S' and in_s:
            in_s = False
            durs.append(t[i] - start)
    return (float(np.mean(durs)) if durs else 0.0), durs


def mean_last(result, sp, n=1500):
    return float(np.mean(result[f'[{sp}]'][-n:]))


def period(result, t_start=40):
    t = result['time']; cb = result['[Cb]']
    mask = t >= t_start
    pk, _ = find_peaks(cb[mask], prominence=0.05, distance=10)
    pt = t[mask][pk]
    return float(np.mean(np.diff(pt))) if len(pt) > 1 else float('nan')


def evaluate(ezh2_rewire_params=None, basal=None, e2f_amp=None, verbose=True):
    """Full report for one rewire parameter set (HU=0)."""
    kw = dict(ezh2_rewire_params=ezh2_rewire_params, basal=basal, e2f_amp=e2f_amp)
    cyc = run(shh=0.5, **kw)          # GNP + SHH (cycling)
    g0 = run(shh=0.0, **kw)           # quiescent reference
    hhi = run(shh=0.5, hhi=1.0, **kw)  # HHi

    ez_cyc = mean_last(cyc, 'EZH2')
    ez_g0 = mean_last(g0, 'EZH2')
    ez_gdc = mean_last(hhi, 'EZH2')
    by, thr = mean_ezh2_by_phase(cyc)
    sdur, _ = s_phase_duration(cyc)
    per = period(cyc)

    # within-cycle ratios vs G0
    rG1, rS, rG2 = by['G1'] / ez_g0, by['S'] / ez_g0, by['G2M'] / ez_g0
    invariants = dict(
        mean_cyc_EZH2=ez_cyc,
        G0_over_cyc=ez_g0 / ez_cyc if ez_cyc else np.nan,
        GDC_reduction=1 - ez_gdc / ez_cyc if ez_cyc else np.nan,
    )
    gradient = dict(G1=rG1, S=rS, G2=rG2)

    if verbose:
        p = dict(EZH2_REWIRE_DEFAULTS); p.update(ezh2_rewire_params or {})
        print(f"  params: K_Ma={p['K_Ma_EZH2']:.3f} k_deg={p['k_EZH2_deg']:.3f} "
              f"k_tx={p['k_EZH2_translation']:.3f}")
        print(f"  period={per:.1f}h  S-dur={sdur:.1f}h  thr_ma={thr['thr_ma']:.3f} "
              f"thr_cb={thr['thr_cb']:.3f}")
        print(f"  [INVARIANTS]  mean_cyc_EZH2={invariants['mean_cyc_EZH2']:.3f} (~0.49)"
              f"   G0/cyc={invariants['G0_over_cyc']:.2f} (~0.60)"
              f"   GDC_red={invariants['GDC_reduction']*100:.0f}% (25-47%)")
        print(f"  [GRADIENT]    G1={rG1:.2f} (1.16)  S={rS:.2f} (1.32)  G2={rG2:.2f} (1.48)")
    return dict(invariants=invariants, gradient=gradient, period=per,
                s_duration=sdur, by_phase=by, ez_g0=ez_g0)


if __name__ == "__main__":
    print("=" * 78)
    print("v43 EZH2-REWIRE CALIBRATION (HU=0)")
    print("=" * 78)
    # Fixed rewire shape (gives correct G2/G1 rise); search the G0/G1 floor via basal.
    REW = dict(k_EZH2_deg=0.20, K_Ma_EZH2=0.35, k_EZH2_translation=0.80)
    print(f"Fixed rewire shape: {REW}")
    print("Searching basal floor + Ma amplitude to set G0/cyc and gradient normalization.\n")
    for basal in (0.030, 0.020, 0.012):
        for e2f_amp in (2.5, 3.5):
            print(f"* basal={basal}  E2F_amp={e2f_amp}")
            evaluate(REW, basal=basal, e2f_amp=e2f_amp)
            print()
