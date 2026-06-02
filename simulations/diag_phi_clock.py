"""Phase-gated clock slowdown: the intra-S checkpoint slows the cell-cycle CLOCK
specifically during S/G2 (replication window) as a function of phi(HU).

Why: no single-flux amplitude scaling moves the GG period (diag_phi_levers.py);
only the global timescale eps does. So we make eps phase-dependent -- reduced toward
eps0*phi ONLY while replication is ongoing (Ma high), normal elsewhere:

    sphase  := Ma^n / (K_S_phi^n + Ma^n)        # ~1 in S/G2, ~0 in G1/M
    eps     := eps0 * (1 - sphase * (1 - phi_HU))

At phi=1 (HU=0): eps = eps0 everywhere -> identical to v42.
At phi<1 in S:   eps -> eps0*phi  -> clock ticks slower -> S/G2 prolonged.
phi >= phi_min > 0 -> graded slowing, never arrest. No single switch is throttled
relative to its antagonist, so no free-CycB inflation.

Biologically: replication stress (low dNTP, low phi) holds/slows the cell-cycle
engine within S until replication can complete (intra-S / ATR-CHK1 checkpoint).

Tests phi_min=0.45,K=1.0,h=3 across HU and a couple K_S_phi/n_S settings.
Run:  ./venv/bin/python simulations/diag_phi_clock.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
from src.build_model_v43_checkpoint import build_model_v43
import tellurium as te

REW = dict(K_Ma_EZH2=0.35, K_Me_EZH2=0.35, w_Me=0.65,
           k_EZH2_deg=0.20, k_EZH2_translation=0.80)

# Build base, strip the builder's phi-on-Ma coupling.
MA_ON = "Ma_activation: Mai -> Ma; Vm1a * (Mai / (K1a + Mai)) * Pa * phi_HU * eps;"
MA_OFF = "Ma_activation: Mai -> Ma; Vm1a * (Mai / (K1a + Mai)) * Pa * eps;"
BASE = build_model_v43(rewire_ezh2=True, ezh2_rewire_params=REW).replace(MA_ON, MA_OFF)
assert "eps = 150;" in BASE


def make(K_S, n_S):
    """Convert eps to a phase-gated assignment rule."""
    block = (f"eps0 = 150;\n"
             f"K_S_phi = {K_S};\n"
             f"n_S = {n_S};\n"
             f"sphase := Ma^n_S / (K_S_phi^n_S + Ma^n_S);\n"
             f"eps := eps0 * (1 - sphase * (1 - phi_HU));\n")
    return BASE.replace("eps = 150;", block)


def period(r, t_start=50):
    t = r['time']; cb = r['[Cb]']; m = t >= t_start
    pk, _ = find_peaks(cb[m], prominence=0.05, distance=10)
    pt = t[m][pk]
    return float(np.mean(np.diff(pt))) if len(pt) > 1 else float('nan')


def n_div(r, t_start=50):
    t = r['time']; cb = r['[Cb]']; m = t >= t_start
    pk, _ = find_peaks(cb[m], prominence=0.05, distance=10)
    return len(pk)


def sim(model, hu, t_end=300, n_pts=30000):
    rr = te.loada(model)
    rr['SHH'] = 0.5; rr['Ptch1_copy_number'] = 1.0; rr['HU'] = hu
    rr['k_EZH2_mRNA_synth_basal'] = 0.02; rr['k_EZH2_mRNA_synth_E2F'] = 3.5
    return rr.simulate(0, t_end, n_pts)


def report(K_S, n_S):
    model = make(K_S, n_S)
    base = sim(model, 0.0)
    m0 = base['time'] >= 50
    thr_ma = 0.20 * base['[Ma]'][m0].max()
    def ez_sg2(r):
        mm = r['time'] >= 50
        sel = mm & (r['[Ma]'] > thr_ma)
        return float(r['[EZH2]'][sel].mean()) if sel.sum() else np.nan
    ez0 = ez_sg2(base)
    print(f"\n=== phase-gated clock: K_S_phi={K_S} n_S={n_S}  (thr_ma={thr_ma:.3f}) ===")
    print(f"{'HU':>5}{'phi':>6}{'period':>8}{'div':>5}{'S_dwell':>9}{'S_hrs':>7}"
          f"{'Ma_max':>8}{'Cb_max':>8}{'EZH2_SG2':>9}{'boost':>7}")
    for hu in (0.0, 0.5, 1.0, 1.5, 2.0):
        r = sim(model, hu)
        mm = r['time'] >= 50
        phi = 0.45 + 0.55 / (1 + (hu / 1.0) ** 3)
        sdw = (r['[Ma]'][mm] > thr_ma).mean()
        per = period(r)
        print(f"{hu:>5.1f}{phi:>6.2f}{per:>8.1f}{n_div(r):>5d}{sdw:>9.3f}"
              f"{sdw*per:>7.1f}{r['[Ma]'][mm].max():>8.3f}{r['[Cb]'][mm].max():>8.3f}"
              f"{ez_sg2(r):>9.3f}{ez_sg2(r)/ez0:>7.2f}")


if __name__ == "__main__":
    print("=" * 80)
    print("INTRA-S CHECKPOINT via phase-gated cell-cycle CLOCK slowdown (phi on eps in S)")
    print("=" * 80)
    for K_S, n_S in ((0.25, 4), (0.35, 4), (0.20, 6)):
        report(K_S, n_S)
