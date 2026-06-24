"""v44 HH/MYCN/CyclinD1 between-condition calibration (Section A targets).

Measures steady/mean levels of CyclinD1 (Cd), MYCN, Gli1, EZH2 in the four conditions and
scores the cross-condition ratios:
  GNP    = SHH 0.5, MYCN_amplification 1.0, HHi 0
  GNP+HHi= SHH 0.5, MYCN_amplification 1.0, HHi 1
  MB     = SHH 0.5, MYCN_amplification 2.8, HHi 0
  MB+HHi = SHH 0.5, MYCN_amplification 2.8, HHi 1

Targets (Section A):
  CyclinD1  GNP+HHi/GNP=0.14   MB+HHi/MB=0.40   MB/GNP=5.07
  MYCN      GNP+HHi/GNP=0.78   MB+HHi/MB=0.86   MB/GNP=2.8
  Gli1      GNP+HHi/GNP<0.01 (>99% reduction)
  EZH2      MB/GNP=2.05        GNP+HHi/GNP=0.53-0.75 (25-47% reduction; cell-cycle coupled)

Edit PARAMS, re-run.  ./venv/bin/python simulations/v44_calibrate_hh.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

# HH/MYCN/CyclinD1 parameters under calibration (others stay at builder defaults)
PARAMS = {
    "k_Cd_tx_basal": 0.3, "k_Cd_tx_Gli_max": 4.0, "k_Cd_tx_MYCN": 15,
    "K_MYCN_Cd": 1.5, "n_MYCN_Cd": 3,
    "k_MYCN_synth_basal": 0.3, "k_MYCN_synth_Gli": 0.102, "k_MYCN_deg": 1.0,
    "k_Cd_translation": 0.26,
    # HH-pathway de-saturation (so MB Ptch1-loss -> higher Gli -> graded CyclinD1)
    "K_Ptch_Smo": 0.4, "K_Gli_act_CycD": 0.4,
}
# MB = Ptch1 loss (constitutive Hedgehog) + MYCN amplification; GNP = WT (Ptch1=1), SHH-driven.
COND = {
    "GNP":     dict(SHH=0.5, Ptch1_copy_number=1.0, MYCN_amplification=1.0, HHi=0.0),
    "GNP+HHi": dict(SHH=0.5, Ptch1_copy_number=1.0, MYCN_amplification=1.0, HHi=1.0),
    "MB":      dict(SHH=0.5, Ptch1_copy_number=0.1, MYCN_amplification=2.8, HHi=0.0),
    "MB+HHi":  dict(SHH=0.5, Ptch1_copy_number=0.1, MYCN_amplification=2.8, HHi=1.0),
}
SPECIES = ["Cd", "MYCN", "Gli1", "EZH2"]


def measure(params):
    m = build_model_v44(with_ezh2=True, with_hh=True, params=params)
    out = {}
    for name, inp in COND.items():
        try:
            rr = te.loada(m)
            for k, v in inp.items():
                rr[k] = v
            rr.integrator.setValue("absolute_tolerance", 1e-9)
            rr.integrator.setValue("relative_tolerance", 1e-6)
            r = rr.simulate(0, 16000, 32000, selections=["time"] + SPECIES)
            sel = r["time"] >= 6000
            out[name] = {sp: float(np.mean(r[sp][sel])) for sp in SPECIES}
        except Exception as e:
            print(f"  [{name} FAILED: {str(e)[:50]}]")
            out[name] = {sp: float('nan') for sp in SPECIES}
    return out


def report(o):
    def ratio(sp, a, b):
        return o[a][sp] / o[b][sp] if o[b][sp] else float('nan')
    print(f"\n{'species':<9}" + "".join(f"{c:>9}" for c in COND))
    for sp in SPECIES:
        print(f"{sp:<9}" + "".join(f"{o[c][sp]:>9.3f}" for c in COND))
    print("\nRATIO                       model     target")
    rows = [
        ("CyclinD1 GNP+HHi/GNP", ratio("Cd", "GNP+HHi", "GNP"), "0.14"),
        ("CyclinD1 MB+HHi/MB",   ratio("Cd", "MB+HHi", "MB"),   "0.40"),
        ("CyclinD1 MB/GNP",      ratio("Cd", "MB", "GNP"),      "5.07"),
        ("MYCN     GNP+HHi/GNP", ratio("MYCN", "GNP+HHi", "GNP"), "0.78"),
        ("MYCN     MB+HHi/MB",   ratio("MYCN", "MB+HHi", "MB"),   "0.86"),
        ("MYCN     MB/GNP",      ratio("MYCN", "MB", "GNP"),      "2.8"),
        ("Gli1     GNP+HHi/GNP", ratio("Gli1", "GNP+HHi", "GNP"), "<0.01"),
        ("EZH2     MB/GNP",      ratio("EZH2", "MB", "GNP"),      "2.05"),
        ("EZH2     GNP+HHi/GNP", ratio("EZH2", "GNP+HHi", "GNP"), "0.53-0.75"),
    ]
    for lab, val, tgt in rows:
        print(f"{lab:<24}{val:>9.3f}     {tgt}")


if __name__ == "__main__":
    print("=" * 60)
    print("v44 HH/MYCN/CyclinD1 BETWEEN-CONDITION CALIBRATION")
    print("PARAMS:", PARAMS)
    report(measure(PARAMS))
