"""Find the phi(HU) coupling point that genuinely slows the cycle and prolongs S.

Background
----------
Scaling Ma_activation (CycA/CDK2 buildup) by phi only SHRINKS CycA amplitude; the
Gerard-Goldbeter period is robust to CycA amplitude, so period/S-dwell stay flat
(see diag_phi.py). The canonical intra-S -> G2/M effector is CDK1 (CycB/Mb)
activation: an unfinished S-phase (low phi) inhibits CDK1 activation (via ATR-CHK1
-> CDC25) to DELAY mitosis. So phi should scale Mb_activation, not Ma_activation.

This script reverts the builder's Ma coupling and re-injects phi at alternative
points via pure string surgery on the built model, then compares HU dose-response:
  - period (CycB peaks)            -> should LENGTHEN with HU
  - S_dwell (Ma above FIXED DMSO thr) -> should LENGTHEN with HU (longer S/G2)
  - Cb (free CycB) stats           -> watch for phase-detector inflation artifact
  - emergent EZH2 in the Ma-high window -> should RISE ~1.31x at HU=1.0

Run:  ./venv/bin/python simulations/diag_phi_coupling.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from src.build_model_v43_checkpoint import build_model_v43
import tellurium as te

REW = dict(K_Ma_EZH2=0.35, K_Me_EZH2=0.35, w_Me=0.65,
           k_EZH2_deg=0.20, k_EZH2_translation=0.80)

MA_RXN = "Ma_activation: Mai -> Ma; Vm1a * (Mai / (K1a + Mai)) * Pa * phi_HU * eps;"
MA_OFF = "Ma_activation: Mai -> Ma; Vm1a * (Mai / (K1a + Mai)) * Pa * eps;"
MB_OFF = "Mb_activation: Mbi -> Mb; Vm1b * (Mbi / (K1b + Mbi)) * Pb * eps;"
MB_ON  = "Mb_activation: Mbi -> Mb; Vm1b * (Mbi / (K1b + Mbi)) * Pb * phi_HU * eps;"

# Base model: built with builder (phi on Ma), then REVERT Ma coupling so phi_HU is
# defined but unused; coupling variants are injected from this clean base.
_BASE = build_model_v43(rewire_ezh2=True, ezh2_rewire_params=REW)
assert MA_RXN in _BASE
_BASE_NOCPL = _BASE.replace(MA_RXN, MA_OFF)
assert MB_OFF in _BASE_NOCPL


def make(coupling):
    m = _BASE_NOCPL
    if coupling == "Ma":
        m = m.replace(MA_OFF, MA_RXN)
    elif coupling == "Mb":
        m = m.replace(MB_OFF, MB_ON)
    elif coupling == "MaMb":
        m = m.replace(MA_OFF, MA_RXN).replace(MB_OFF, MB_ON)
    elif coupling == "none":
        pass
    else:
        raise ValueError(coupling)
    return m


def sim(coupling, hu, t_end=300, n_pts=30000):
    rr = te.loada(make(coupling))
    rr['SHH'] = 0.5
    rr['Ptch1_copy_number'] = 1.0
    rr['HU'] = hu
    rr['k_EZH2_mRNA_synth_basal'] = 0.02
    rr['k_EZH2_mRNA_synth_E2F'] = 3.5
    return rr.simulate(0, t_end, n_pts)


def period(r, t_start=50):
    t = r['time']; cb = r['[Cb]']; m = t >= t_start
    pk, _ = find_peaks(cb[m], prominence=0.05, distance=10)
    pt = t[m][pk]
    return float(np.mean(np.diff(pt))) if len(pt) > 1 else float('nan')


def report(coupling):
    base = sim(coupling, 0.0)
    m0 = base['time'] >= 50
    thr_ma = 0.20 * base['[Ma]'][m0].max()
    # emergent EZH2 in Ma-high (S/G2) window, FIXED DMSO threshold
    def ezh2_shigh(r):
        mm = r['time'] >= 50
        sel = mm & (r['[Ma]'] > thr_ma)
        return float(r['[EZH2]'][sel].mean()) if sel.sum() else np.nan
    ez0 = ezh2_shigh(base)
    print(f"\n=== coupling = {coupling}   (fixed thr_ma={thr_ma:.3f}) ===")
    print(f"{'HU':>5}{'phi':>6}{'period':>8}{'S_dwell':>9}{'Ma_max':>8}"
          f"{'Cb_max':>8}{'Mb_max':>8}{'EZH2_SG2':>9}{'boost':>7}")
    for hu in (0.0, 0.5, 1.0, 2.0):
        r = sim(coupling, hu)
        mm = r['time'] >= 50
        phi = 0.45 + 0.55 / (1 + (hu / 1.0) ** 3)
        sdw = (r['[Ma]'][mm] > thr_ma).mean()
        ez = ezh2_shigh(r)
        print(f"{hu:>5.1f}{phi:>6.2f}{period(r):>8.1f}{sdw:>9.3f}"
              f"{r['[Ma]'][mm].max():>8.3f}{r['[Cb]'][mm].max():>8.3f}"
              f"{r['[Mb]'][mm].max():>8.3f}{ez:>9.3f}{ez/ez0:>7.2f}")


if __name__ == "__main__":
    print("=" * 78)
    print("phi(HU) COUPLING-POINT SEARCH: which flux genuinely slows S / lengthens cycle")
    print("=" * 78)
    for c in ("Ma", "Mb", "MaMb"):
        report(c)
