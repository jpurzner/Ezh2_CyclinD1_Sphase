"""v44 build test: Heldt 2018 core + CycB/CDK1 (MPF) mitotic switch + division reset.

Adds to the Heldt one-shot G1/S/replication model:
  - CyclinB/CDK1 hysteresis switch: MPF (active) <-> preMPF (Tyr15-P), with Cdc25
    (activated by MPF, inhibited by Chk1) and Wee1 (inhibited by MPF, activated by Chk1).
  - CHK1 replication checkpoint gated on aRc (active replication forks = unfinished S):
    holds MPF inactive until DNA replication completes (Dna -> 1, forks disassemble).
  - Mitotic APC/Cdc20: activated by MPF, degrades MPF/preMPF and CycA -> mitotic exit.
  - Division event at mitotic exit: reset Dna=0, re-license origins, dephosphorylate Rb,
    halve cyclins -> next G1. Makes the one-shot model CYCLE.

Goal of this test: confirm sustained, clean cycling with sane phase structure.
Run:  ./venv/bin/python simulations/v44_build_test.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.signal import find_peaks
import tellurium as te

HELDT = open(os.path.join(os.path.dirname(__file__), "..",
             "models_external", "heldt2018.ant")).read()

MITOSIS_BLOCK = """
  // ===== v44 ADDITION: CyclinB/CDK1 (MPF) mitotic switch + Cdc20 + Chk1 checkpoint =====
  species MPF in Cell, preMPF in Cell, Cdc20 in Cell;
  MPF = 0; preMPF = 0; Cdc20 = 0;

  Cb := MPF + preMPF;                                  // total CyclinB-CDK1 (readout)
  Chk1 := aRc/(jChk + aRc);                            // checkpoint: active forks = unfinished S
  // basal a25/aWee let the switch ignite once Chk1 clears (avoids cold-start lock)
  Cdc25a := (a25 + (1 - a25)*MPF^nMpf/(KmMpf^nMpf + MPF^nMpf))/(1 + wChk*Chk1);
  Wee1a  := (aWee + (1 - aWee)*KmMpf^nMpf/(KmMpf^nMpf + MPF^nMpf))*(1 + wChkW*Chk1);

  SynCycB: => MPF; Cell*kSyCb*E2f;                     // CycB-CDK1 synthesized (active) via E2f
  Wee1phos: MPF => preMPF; Cell*kWee*Wee1a*MPF;        // Tyr15 phosphorylation (inactivate)
  Cdc25dephos: preMPF => MPF; Cell*k25*Cdc25a*preMPF;  // Tyr15 dephosphorylation (activate)
  DegMPF: MPF => ; Cell*(kDeCbBas + kDeCb*Cdc20)*MPF;
  DegpreMPF: preMPF => ; Cell*(kDeCbBas + kDeCb*Cdc20)*preMPF;
  Cdc20act: => Cdc20; Cell*kaCdc20*(MPF^nCdc20/(KmCdc20^nCdc20 + MPF^nCdc20))*(1 - Cdc20);
  Cdc20inact: Cdc20 => ; Cell*kiCdc20*Cdc20;
  DegCycACdc20: Ca => ; Cell*kDeCaCdc20*Cdc20*Ca;      // mitotic CycA destruction

  kSyCb = 0.01; kDeCb = 0.08; kDeCbBas = 0.004;
  kWee = 0.4; k25 = 0.7; KmMpf = 0.35; nMpf = 4;
  a25 = 0.15; aWee = 0.1;
  wChk = 12; wChkW = 4; jChk = 0.03;
  kaCdc20 = 0.3; kiCdc20 = 0.12; KmCdc20 = 0.5; nCdc20 = 8; kDeCaCdc20 = 1.0;

  // Division at mitotic entry (MPF crosses high): reset to a G1 daughter.
  E_div: at (MPF > 1): Dna = 0, Rc = 1, pRc = 0, aRc = 0, iRc = 0, Rb = Rb + pRb, pRb = 0, Ce = Ce/2, Ca = Ca/2, E1 = E1/2, MPF = 0, preMPF = 0, Cdc20 = 0 ;
"""


def build():
    assert "\nend" in HELDT
    return HELDT.replace("\nend", MITOSIS_BLOCK + "\nend")


def main():
    m = build()
    rr = te.loada(m)
    res = rr.simulate(0, 6000, 12000, selections=[
        "time", "E2f", "Ce", "Ca", "Dna", "MPF", "Cb", "Cdc20", "pRb", "aRc", "Chk1"])
    t = res["time"]
    # detect cycles via MPF peaks
    mpf = res["MPF"]
    pk, props = find_peaks(mpf, prominence=0.2, distance=100)
    print(f"MPF peaks (mitoses): {len(pk)} at t(min) = {np.round(t[pk]).astype(int).tolist()}")
    if len(pk) > 1:
        per = np.diff(t[pk])
        print(f"period(min) = {np.round(per).astype(int).tolist()}  ~{np.mean(per)/60:.2f} h")
    # Dna completes?
    print(f"Dna range: {res['Dna'].min():.2f} .. {res['Dna'].max():.2f}")
    print(f"MPF range: {mpf.min():.3f} .. {mpf.max():.3f}   Cdc20 max: {res['Cdc20'].max():.3f}")
    print(f"Ca range: {res['Ca'].min():.3f} .. {res['Ca'].max():.3f}")
    print(f"E2f range: {res['E2f'].min():.3f} .. {res['E2f'].max():.3f}")
    print(f"pRb range: {res['pRb'].min():.3f} .. {res['pRb'].max():.3f}")
    # downsampled trajectory every 250 min for first 4000 min
    print("\n t(min)  E2f   Ce    Ca    Dna   MPF   Cdc20  pRb   Chk1")
    for i in range(0, min(len(t), 8000), 500):
        print(f"  {t[i]:5.0f} {res['E2f'][i]:5.2f} {res['Ce'][i]:5.2f} {res['Ca'][i]:5.2f} "
              f"{res['Dna'][i]:5.2f} {res['MPF'][i]:5.2f} {res['Cdc20'][i]:5.2f} "
              f"{res['pRb'][i]:5.2f} {res['Chk1'][i]:5.2f}")


if __name__ == "__main__":
    main()
