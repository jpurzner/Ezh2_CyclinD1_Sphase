"""
Diagnostics for the v43 intra-S checkpoint extension.

Goals:
  1. Confirm v43 with HU=0 reproduces the v42 baseline (checkpoint inert).
  2. Characterise the cell-cycle phase structure (define S-phase operationally).
  3. Measure the effect of HU (replication stress) on:
       - cell-cycle period
       - S-phase duration
       - mean EZH2 protein, and EZH2 in S-phase vs G0
   against the experimental target: 10 uM HU -> 1.31-fold EZH2 in S-phase (Fig.4G).

Run:  ./venv/bin/python simulations/sim_checkpoint_diagnostics.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from src.build_model_v43_checkpoint import build_model_v43, CHECKPOINT_DEFAULTS
import tellurium as te

MODEL_STR = build_model_v43()


def run_sim(shh=0.5, ptch1_cn=1.0, gdc=0.0, ezh2i=0.0, mycn_amp=1.0,
            hu=0.0, t_end=168, n_pts=None, serum_starve=False, cdk46i=False,
            checkpoint=None):
    """Run a v43 simulation. checkpoint = dict of checkpoint param overrides."""
    if n_pts is None:
        n_pts = int(t_end * 100)
    rr = te.loada(MODEL_STR)
    rr['SHH'] = shh
    rr['Ptch1_copy_number'] = ptch1_cn
    rr['GDC0449'] = gdc
    rr['EZH2i'] = ezh2i
    rr['MYCN_amplification'] = mycn_amp
    rr['HU'] = hu
    if ptch1_cn == 0.0:
        rr['Ptch1_mRNA'] = 0.0
        rr['Ptch1_free'] = 0.0
        rr['SHH_Ptch'] = 0.0
        rr['Smo_active'] = 1.0
    elif shh == 0.0:
        rr['SHH_Ptch'] = 0.0
        rr['Ptch1_free'] = 1.0
        rr['Smo_active'] = 0.01
    if serum_starve:
        rr['k_Cd_translation'] = 0.0
    if cdk46i:
        rr['V1'] = 0.0
    if checkpoint:
        for k, v in checkpoint.items():
            rr[k] = v
    return rr.simulate(0, t_end, n_pts)


def count_divisions(result, t_start=10):
    t = result['time']
    cb = result['[Cb]']
    mask = t >= t_start
    peaks, _ = find_peaks(cb[mask], prominence=0.05, distance=10)
    peak_times = t[mask][peaks]
    periods = np.diff(peak_times) if len(peak_times) > 1 else np.array([])
    return len(peaks), peak_times, periods


def phase_assignment(result):
    """Operational phase call per timepoint from cyclin state.

    Uses Ma = active CycA/CDK2 (the S/G2 marker, high ~16% of cycle) and
    Cb = CycB (the M marker). pRBpp is NOT used: it is a brief CycE/CDK2 spike
    present only ~3% of the cycle and is a poor time-in-S integrator.

    G2/M  : Cb (CycB) high                        -> 4N, mitotic
    S     : Ma (CycA/CDK2) high AND Cb low         -> replicating
    G0G1  : Ma low AND Cb low                      -> pre-S / quiescent
    Returns array of {'G0G1','S','G2M'} and the thresholds used.
    """
    ma = result['[Ma]']
    cb = result['[Cb]']
    # thresholds: low/high midpoints (robust, data-free)
    thr_ma = 0.5 * (np.percentile(ma, 95) + np.percentile(ma, 5))
    thr_cb = 0.5 * (np.percentile(cb, 95) + np.percentile(cb, 5))
    phase = np.empty(len(ma), dtype=object)
    for i in range(len(ma)):
        if cb[i] > thr_cb:
            phase[i] = 'G2M'
        elif ma[i] > thr_ma:
            phase[i] = 'S'
        else:
            phase[i] = 'G0G1'
    return phase, dict(thr_ma=float(thr_ma), thr_cb=float(thr_cb))


def s_phase_duration_per_cycle(result, t_start=20):
    """Mean contiguous S-phase duration (h) per cycle after t_start."""
    t = result['time']
    phase, _ = phase_assignment(result)
    mask = t >= t_start
    t, phase = t[mask], phase[mask]
    durations, in_s, start = [], False, None
    for i in range(len(t)):
        if phase[i] == 'S' and not in_s:
            in_s, start = True, t[i]
        elif phase[i] != 'S' and in_s:
            in_s = False
            durations.append(t[i] - start)
    # drop incomplete trailing interval already handled (only closed intervals appended)
    return float(np.mean(durations)) if durations else 0.0, durations


def mean_in_phase(result, species, target_phase, t_start=20):
    t = result['time']
    vals = result[f'[{species}]']
    phase, _ = phase_assignment(result)
    mask = (t >= t_start)
    sel = mask & (phase == target_phase)
    return float(np.mean(vals[sel])) if sel.sum() > 0 else np.nan


def mean_last(result, species, n=1500):
    return float(np.mean(result[f'[{species}]'][-n:]))


# =====================================================================
if __name__ == "__main__":
    print("=" * 78)
    print("v43 CHECKPOINT DIAGNOSTICS")
    print("=" * 78)
    print("Checkpoint defaults:", CHECKPOINT_DEFAULTS)

    # ---- 1. DMSO (HU=0) vs HU in cycling GNPs ----
    dmso = run_sim(shh=0.5, ptch1_cn=1.0, hu=0.0)
    hu1 = run_sim(shh=0.5, ptch1_cn=1.0, hu=1.0)

    for label, r in [('DMSO (HU=0)', dmso), ('HU=1', hu1)]:
        n, _, per = count_divisions(r)
        period = np.mean(per) if len(per) else float('nan')
        sdur, _ = s_phase_duration_per_cycle(r)
        ezh2_mean = mean_last(r, 'EZH2')
        ezh2_S = mean_in_phase(r, 'EZH2', 'S')
        ezh2_G0 = mean_in_phase(r, 'EZH2', 'G0G1')
        chk1 = mean_last(r, 'CHK1_active')
        dstress = mean_last(r, 'D_stress')
        print(f"\n--- {label} ---")
        print(f"  divisions(168h)={n}  period={period:.1f}h  S-phase={sdur:.1f}h")
        print(f"  EZH2 mean={ezh2_mean:.3f}  EZH2_S={ezh2_S:.3f}  EZH2_G0={ezh2_G0:.3f}"
              f"  S/G0={ezh2_S/ezh2_G0 if ezh2_G0==ezh2_G0 and ezh2_G0>0 else float('nan'):.2f}")
        print(f"  CHK1_active={chk1:.3f}  D_stress={dstress:.3f}")

    # EZH2-in-S boost: HU vs DMSO
    ezh2_S_dmso = mean_in_phase(dmso, 'EZH2', 'S')
    ezh2_S_hu = mean_in_phase(hu1, 'EZH2', 'S')
    boost = ezh2_S_hu / ezh2_S_dmso if ezh2_S_dmso > 0 else float('nan')
    print(f"\n>>> HU/DMSO EZH2-in-S boost = {boost:.3f}  (target 1.31)")

    # ---- 2. Time-course figure ----
    fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharex=True)
    species_rows = [('Cb', 'Cyclin B (division marker)'),
                    ('pRBpp', 'pRBpp (EZH2 driver)'),
                    ('EZH2', 'EZH2 protein')]
    for col, (r, title) in enumerate([(dmso, 'DMSO (HU=0)'), (hu1, 'HU=1')]):
        t = r['time']
        for row, (sp, lbl) in enumerate(species_rows):
            ax = axes[row, col]
            ax.plot(t, r[f'[{sp}]'], color='#1b9e77', lw=1)
            if row == 0:
                ax.set_title(title, fontsize=11, fontweight='bold')
            ax.set_ylabel(lbl, fontsize=8)
            ax.grid(alpha=0.2)
            # overlay checkpoint signal on Cb row
            if row == 0:
                ax2 = ax.twinx()
                ax2.plot(t, r['[CHK1_active]'], color='#d95f02', lw=0.8, alpha=0.7)
                ax2.set_ylabel('CHK1_active', color='#d95f02', fontsize=7)
                ax2.set_ylim(-0.05, 1.05)
        axes[-1, col].set_xlabel('Time (h)', fontsize=9)
    fig.suptitle('v43 checkpoint: DMSO vs HU (cycling GNP + SHH)',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    out = 'simulations/fig_checkpoint_diagnostics.png'
    plt.savefig(out, dpi=180, bbox_inches='tight')
    plt.close()
    print(f"\nSaved: {out}")
