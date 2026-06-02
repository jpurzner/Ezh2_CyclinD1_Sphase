"""v44 CORE HYPOTHESIS TEST (the reason for the whole project):
   longer S-phase (HU) -> more EZH2 -> stronger CyclinD repression -> longer SUBSEQUENT G1.

On the explicit-replication Heldt-based v44 model, sweep HU and measure per-cycle:
  S duration, G1 duration, mean EZH2 in S-phase, mean Cd.
Run with EZH2->Cd feedback ON (EZH2i=0) and OFF (EZH2i=1):
  - EZH2-in-S should rise with HU in BOTH (longer dwell in the E2f x CDK2 window).
  - G1 should lengthen with HU ONLY when feedback is ON (EZH2 -> represses Cd -> slows
    Rb phosphorylation -> longer G1). With feedback OFF, G1 stays flat -> isolates the loop.

Run:  ./venv/bin/python simulations/v44_ezh2_hypothesis.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

M = build_model_v44(with_ezh2=True)


def run(hu, ezh2i, t_end=12000, n_pts=48000):
    rr = te.loada(M)
    rr['HU'] = hu
    rr['EZH2i'] = ezh2i
    return rr.simulate(0, t_end, n_pts,
                       selections=["time", "Dna", "MPF", "aRc", "Ce", "Ca", "Cd", "EZH2"])


def per_cycle(res, settle=4000):
    # S-onset is detected by ORIGIN FIRING (aRc>0.1), which depends on CDK2 not fork speed,
    # so G1 is measured cleanly. (Using Dna>0.05 leaks slow early-S into G1 under HU.)
    t = res["time"]; dna = res["Dna"]; arc = res["aRc"]
    div = np.where((dna[:-1] > 0.9) & (dna[1:] < 0.1))[0] + 1
    div = div[t[div] > settle]
    g1, sdur, ez_s, ez_all, cd = [], [], [], [], []
    for k in range(len(div) - 1):
        d0, d1 = div[k], div[k + 1]
        st = t[d0:d1]
        af = np.where(arc[d0:d1] > 0.1)[0]          # S-onset = origins fired
        a95 = np.where(dna[d0:d1] > 0.95)[0]        # S-end = replication complete
        if not (len(af) and len(a95)):
            continue
        ts0, ts1 = st[af[0]], st[a95[0]]
        g1.append(ts0 - st[0]); sdur.append(ts1 - ts0)
        smask = (t >= ts0) & (t <= ts1)             # S-phase window
        ez_s.append(res["EZH2"][smask].mean())
        cyc = (t >= st[0]) & (t <= st[-1])
        ez_all.append(res["EZH2"][cyc].mean()); cd.append(res["Cd"][cyc].mean())
    f = lambda x: float(np.mean(x)) if x else float('nan')
    return dict(g1=f(g1) / 60, s=f(sdur) / 60, ez_s=f(ez_s), ez_all=f(ez_all), cd=f(cd),
                ncyc=len(g1))


if __name__ == "__main__":
    print("=" * 80)
    print("v44 HYPOTHESIS: HU -> longer S -> more EZH2 -> CyclinD repression -> longer G1")
    print("=" * 80)
    for ezh2i, lab in ((0, "feedback ON (EZH2->Cd)"), (1, "feedback OFF (EZH2i=1)")):
        print(f"\n--- {lab} ---")
        print(f"{'HU':>5}{'G1_h':>8}{'S_h':>8}{'EZH2_inS':>10}{'EZH2_mean':>11}"
              f"{'Cd_mean':>9}{'ncyc':>6}")
        base_ezs = None
        for hu in (0.0, 0.5, 1.0, 2.0):
            r = per_cycle(run(hu, ezh2i))
            if base_ezs is None:
                base_ezs = r['ez_s']
            boost = r['ez_s'] / base_ezs
            print(f"{hu:>5.1f}{r['g1']:>8.2f}{r['s']:>8.2f}{r['ez_s']:>10.3f}"
                  f"{r['ez_all']:>11.3f}{r['cd']:>9.3f}{r['ncyc']:>6d}   "
                  f"EZH2-in-S boost={boost:.2f}")
    print("\nExpect: EZH2-in-S rises with HU (both); G1 rises with HU only when feedback ON.")
