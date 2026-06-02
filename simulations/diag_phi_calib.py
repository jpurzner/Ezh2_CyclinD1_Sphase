"""Joint calibration: can ONE EZH2 turnover satisfy BOTH targets?
  (i)  within-cycle gradient at HU=0:  G2/G0 ~1.48 (>=1.30), monotone G1<S<G2
  (ii) emergent HU boost at HU=1.0:    S-phase EZH2 ~1.31x DMSO

Tension: (i) wants FAST EZH2 turnover (track phase); (ii) wants SLOW turnover
(integrate S-dwell). Scan k_EZH2_deg x phi-strength (K_HU_phi) on the phase-gated
clock model and print both metrics so we can see if a sweet spot exists.

Run:  ./venv/bin/python simulations/diag_phi_calib.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from src.build_model_v43_checkpoint import build_model_v43
import tellurium as te

MA_ON = "Ma_activation: Mai -> Ma; Vm1a * (Mai / (K1a + Mai)) * Pa * phi_HU * eps;"
MA_OFF = "Ma_activation: Mai -> Ma; Vm1a * (Mai / (K1a + Mai)) * Pa * eps;"
K_S_PHI, N_S = 0.20, 6


def build(k_deg, K_HU_phi, phi_min=0.45, h=3.0):
    rew = dict(K_Ma_EZH2=0.35, K_Me_EZH2=0.35, w_Me=0.65,
               k_EZH2_deg=k_deg, k_EZH2_translation=0.80)
    m = build_model_v43(checkpoint_params=dict(phi_min=phi_min, K_HU_phi=K_HU_phi, h_HU=h),
                        rewire_ezh2=True, ezh2_rewire_params=rew).replace(MA_ON, MA_OFF)
    block = (f"eps0 = 150;\nK_S_phi = {K_S_PHI};\nn_S = {N_S};\n"
             f"sphase := Ma^n_S / (K_S_phi^n_S + Ma^n_S);\n"
             f"eps := eps0 * (1 - sphase * (1 - phi_HU));\n")
    return m.replace("eps = 150;", block)


def sim(model, hu, shh=0.5, t_end=300, n_pts=30000):
    rr = te.loada(model)
    rr['SHH'] = shh; rr['Ptch1_copy_number'] = 1.0; rr['HU'] = hu
    rr['k_EZH2_mRNA_synth_basal'] = 0.02; rr['k_EZH2_mRNA_synth_E2F'] = 3.5
    if shh == 0.0:
        rr['SHH_Ptch'] = 0.0; rr['Ptch1_free'] = 1.0; rr['Smo_active'] = 0.01
    return rr.simulate(0, t_end, n_pts)


def period(r, t_start=50):
    t = r['time']; cb = r['[Cb]']; m = t >= t_start
    pk, _ = find_peaks(cb[m], prominence=0.05, distance=10)
    pt = t[m][pk]
    return float(np.mean(np.diff(pt))) if len(pt) > 1 else float('nan')


def phases(r, thr_ma, thr_cb, t_start=50):
    t = r['time']; ma = r['[Ma]']; cb = r['[Cb]']; m = t >= t_start
    g1 = m & (ma <= thr_ma) & (cb <= thr_cb)
    s = m & (ma > thr_ma) & (cb <= thr_cb)
    g2 = m & (cb > thr_cb)
    return g1, s, g2


def metrics(k_deg, K_HU_phi):
    model = build(k_deg, K_HU_phi)
    cyc = sim(model, 0.0)
    g0 = sim(model, 0.0, shh=0.0)
    m0 = cyc['time'] >= 50
    thr_ma = 0.20 * cyc['[Ma]'][m0].max()
    thr_cb = 0.25 * cyc['[Cb]'][m0].max()
    ez = cyc['[EZH2]']
    g1, s, g2 = phases(cyc, thr_ma, thr_cb)
    ez_g0 = float(g0['[EZH2]'][-1500:].mean())
    rG1 = ez[g1].mean() / ez_g0
    rS = ez[s].mean() / ez_g0
    rG2 = ez[g2].mean() / ez_g0
    mono = rG1 < rS < rG2
    # HU=1.0 emergent boost: S/G2 EZH2 (Ma>thr) vs DMSO
    def ez_sg2(r):
        mm = r['time'] >= 50
        sel = mm & (r['[Ma]'] > thr_ma)
        return float(r['[EZH2]'][sel].mean())
    ez0 = ez_sg2(cyc)
    hu1 = sim(model, 1.0)
    boost = ez_sg2(hu1) / ez0
    return dict(per=period(cyc), per_hu=period(hu1), rG1=rG1, rS=rS, rG2=rG2,
                mono=mono, boost=boost, meanEZH2=float(ez[m0].mean()))


if __name__ == "__main__":
    print("=" * 84)
    print("JOINT CALIB: within-cycle gradient (HU=0) vs emergent HU=1.0 boost")
    print(f"phase-gated clock K_S_phi={K_S_PHI} n_S={N_S}; phi_min=0.45 h=3")
    print("TARGETS: G2/G0>=1.30 & monotone  AND  HU=1.0 boost~1.31")
    print("=" * 84)
    print(f"{'k_deg':>7}{'K_HUphi':>8}{'meanEZH2':>9}{'G1/G0':>7}{'S/G0':>7}"
          f"{'G2/G0':>7}{'mono':>6}{'per0':>6}{'perHU':>6}{'boost':>7}")
    for k_deg in (0.20, 0.12, 0.08, 0.05):
        for K_HU_phi in (1.0, 0.7, 0.5):
            r = metrics(k_deg, K_HU_phi)
            print(f"{k_deg:>7.2f}{K_HU_phi:>8.2f}{r['meanEZH2']:>9.3f}{r['rG1']:>7.2f}"
                  f"{r['rS']:>7.2f}{r['rG2']:>7.2f}{str(r['mono']):>6}"
                  f"{r['per']:>6.1f}{r['per_hu']:>6.1f}{r['boost']:>7.2f}")
