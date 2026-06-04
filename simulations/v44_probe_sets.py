"""Quick probe of a few hand-designed param sets: check stability + arrest/rescue guard +
figure-relevant ratios. Goal: lift CyclinD1 MB/GNP & lower MB+HHi/MB via Gli de-saturation while
keeping MYCN moderate (preserve MB+HHi arrest & numerical stability), and set the EZH2i fold."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

COND = {
    "GNP":     dict(SHH=0.5, Ptch1_copy_number=1.0, MYCN_amplification=1.0, GDC0449=0.0, EZH2i=0.0),
    "GNP+HHi": dict(SHH=0.5, Ptch1_copy_number=1.0, MYCN_amplification=1.0, GDC0449=1.0, EZH2i=0.0),
    "GNP+EZH2i":dict(SHH=0.5,Ptch1_copy_number=1.0, MYCN_amplification=1.0, GDC0449=0.0, EZH2i=1.0),
    "MB":      dict(SHH=0.5, Ptch1_copy_number=0.1, MYCN_amplification=2.8, GDC0449=0.0, EZH2i=0.0),
    "MB+HHi":  dict(SHH=0.5, Ptch1_copy_number=0.1, MYCN_amplification=2.8, GDC0449=1.0, EZH2i=0.0),
    "MB+HHi+EZH2i": dict(SHH=0.5, Ptch1_copy_number=0.1, MYCN_amplification=2.8, GDC0449=1.0, EZH2i=1.0),
}
SPEC = ["time", "Cd_mRNA", "EZH2", "MYCN", "MPF"]

SETS = {
 "default": {},
 # Gli de-saturation: smaller K_Ptch_Smo & K_Gli_act_CycD + higher Gli max -> MB(Ptch1 loss) gets
 # more Gli-driven CyclinD1 than GNP; MYCN kept moderate (n=3) so MB+HHi still drops below threshold.
 "gliDesat_A": dict(K_Ptch_Smo=0.18, K_Gli_act_CycD=0.6, k_Cd_tx_Gli_max=10.0, k_Cd_tx_basal=0.2,
                    k_Cd_tx_MYCN=15.0, K_MYCN_Cd=1.3, n_MYCN_Cd=3, K_EZH2_repression=0.8),
 "gliDesat_B": dict(K_Ptch_Smo=0.15, K_Gli_act_CycD=0.6, k_Cd_tx_Gli_max=13.0, k_Cd_tx_basal=0.2,
                    k_Cd_tx_MYCN=18.0, K_MYCN_Cd=1.3, n_MYCN_Cd=4, K_EZH2_repression=0.8),
 "gliDesat_C": dict(K_Ptch_Smo=0.15, K_Gli_act_CycD=0.55, k_Cd_tx_Gli_max=12.0, k_Cd_tx_basal=0.15,
                    k_Cd_tx_MYCN=16.0, K_MYCN_Cd=1.2, n_MYCN_Cd=4, K_EZH2_repression=0.85),
}


def ndiv(r):
    t = r["time"]; m = t >= 3000; tt = t[m]; dt = tt[1]-tt[0]
    pk, _ = find_peaks(r["MPF"][m], prominence=0.15, distance=int(200/dt))
    return len(pk)


def measure(params):
    try:
        model = build_model_v44(with_ezh2=True, with_hh=True, params=(params or None))
    except Exception as e:
        return None, f"build fail {e}"
    out = {}
    for name, inp in COND.items():
        try:
            rr = te.loada(model)
            rr.integrator.setValue("absolute_tolerance", 1e-9)
            rr.integrator.setValue("relative_tolerance", 1e-6)
            for k, v in inp.items(): rr[k] = v
            r = rr.simulate(0, 9000, 18000, selections=SPEC)
            sel = r["time"] >= 4500
            out[name] = {s: max(float(np.mean(r[s][sel])), 1e-9) for s in SPEC[1:]}
            out[name]["div"] = ndiv(r)
        except Exception as e:
            return None, f"{name} CRASH"
    return out, "ok"


for label, P in SETS.items():
    o, status = measure(P)
    print("=" * 64)
    print(f"[{label}]  {status}")
    if o is None:
        continue
    guard = (o["GNP"]["div"] > 0 and o["MB"]["div"] > 0 and o["GNP+HHi"]["div"] == 0 and
             o["MB+HHi"]["div"] == 0 and o["MB+HHi+EZH2i"]["div"] > 0)
    print("  div: " + "  ".join(f"{c}={o[c]['div']}" for c in COND) + f"   GUARD={'OK' if guard else 'FAIL'}")
    r = lambda sp,a,b: o[a][sp]/o[b][sp]
    print(f"  CyclinD1 GNP+HHi/GNP={r('Cd_mRNA','GNP+HHi','GNP'):.3f} (0.14)  "
          f"MB+HHi/MB={r('Cd_mRNA','MB+HHi','MB'):.3f} (0.40)  MB/GNP={r('Cd_mRNA','MB','GNP'):.2f} (5.07)")
    print(f"  EZH2i fold (GNP)={r('Cd_mRNA','GNP+EZH2i','GNP'):.2f} (2.2)  "
          f"EZH2 MB/GNP={r('EZH2','MB','GNP'):.2f} (2.05)  MYCN MB/GNP={r('MYCN','MB','GNP'):.2f} (2.8)")
