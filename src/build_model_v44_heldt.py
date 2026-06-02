"""
Model v44 - explicit-replication cell cycle on the Heldt 2018 core.

WHY THIS EXISTS
---------------
v42/v43 used the Gerard-Goldbeter relaxation oscillator, whose S-phase duration is a
STRUCTURAL INVARIANT: no molecular concentration can lengthen S (verified exhaustively
in simulations/diag_v44_*; only the global eps clock could, which is unbiological). To
make S-phase duration mechanistic and concentration-dependent, v44 rebuilds the engine
on Heldt, Barr, Cooper, Bakal & Novak 2018 (PNAS 115:2532, BIOMD0000000700): a real-time
(minutes, NO eps) model with EXPLICIT DNA replication -- DNA is synthesized by active
replication complexes (aRc) at fork speed kSyDna until Dna reaches 1.

Heldt is a one-shot G1/S model (proliferation-quiescence decision); it has CycE/CycA but
no mitosis and does not cycle. v44 adds:

  1. A CyclinB/CDK1 (MPF) mitotic switch with Cdc25/Wee1 hysteresis (Cdc25 activated by
     MPF & inhibited by Chk1; Wee1 inhibited by MPF & activated by Chk1; basal a25/aWee
     let the switch ignite once the checkpoint clears).
  2. A CHK1 intra-S checkpoint gated on active replication forks (aRc): holds MPF inactive
     until DNA replication completes (Dna -> 1, forks disassemble, aRc -> 0).
  3. A mitotic APC/Cdc20 (MPF-activated) that destroys CycA/CycB -> mitotic exit, and a
     division-reset EVENT (reset Dna, re-license origins, dephosphorylate Rb, halve
     cyclins) -> sustained cycling.
  4. HU coupling: hydroxyurea / dNTP depletion slows replication forks. vfork(HU) scales
     kSyDna in the DNA-synthesis flux, so HU lengthens S concentration-dependently while
     G1 (set by Rb-E2F/CycD) is unchanged, and mitosis waits for Dna=1 (no arrest).

Later stages re-attach the Hedgehog/MYCN -> CyclinD (Cd) input and the EZH2 layer.

Verified (simulations/v44_fork_doseresponse.py): fork speed x2..x0.25 -> S 3.7..12.4 h,
G1 ~flat 2.8..3.6 h, mitosis always at Dna~1.0, no arrest.
"""
import os

_HELDT_PATH = os.path.join(os.path.dirname(__file__), "..",
                           "models_external", "heldt2018.ant")

# ---- CyclinB/CDK1 (MPF) mitotic switch + Cdc20 + Chk1 checkpoint + division event ----
# Calibrated inline in simulations/v44_build_test.py (sustained cycling, replication-gated).
MITOSIS_BLOCK = """
  // ===== v44: CyclinB/CDK1 (MPF) mitotic switch + Cdc20 + Chk1 intra-S checkpoint =====
  species MPF in Cell, preMPF in Cell, Cdc20 in Cell;
  MPF = 0; preMPF = 0; Cdc20 = 0;

  Cb := MPF + preMPF;                                   // total CyclinB-CDK1 (readout)
  Chk1 := aRc/(jChk + aRc);                             // checkpoint: active forks = unfinished S
  Cdc25a := (a25 + (1 - a25)*MPF^nMpf/(KmMpf^nMpf + MPF^nMpf))/(1 + wChk*Chk1);
  Wee1a  := (aWee + (1 - aWee)*KmMpf^nMpf/(KmMpf^nMpf + MPF^nMpf))*(1 + wChkW*Chk1);

  SynCycB: => MPF; Cell*kSyCb*E2f;                      // CycB-CDK1 synthesized (active) via E2f
  Wee1phos: MPF => preMPF; Cell*kWee*Wee1a*MPF;         // Tyr15 phosphorylation (inactivate)
  Cdc25dephos: preMPF => MPF; Cell*k25*Cdc25a*preMPF;   // Tyr15 dephosphorylation (activate)
  DegMPF: MPF => ; Cell*(kDeCbBas + kDeCb*Cdc20)*MPF;
  DegpreMPF: preMPF => ; Cell*(kDeCbBas + kDeCb*Cdc20)*preMPF;
  Cdc20act: => Cdc20; Cell*kaCdc20*(MPF^nCdc20/(KmCdc20^nCdc20 + MPF^nCdc20))*(1 - Cdc20);
  Cdc20inact: Cdc20 => ; Cell*kiCdc20*Cdc20;
  DegCycACdc20: Ca => ; Cell*kDeCaCdc20*Cdc20*Ca;       // mitotic CycA destruction

  kSyCb = 0.01; kDeCb = 0.08; kDeCbBas = 0.004;
  kWee = 0.4; k25 = 0.7; KmMpf = 0.35; nMpf = 4;
  a25 = 0.15; aWee = 0.1;
  wChk = 12; wChkW = 4; jChk = 0.03;
  kaCdc20 = 0.3; kiCdc20 = 0.12; KmCdc20 = 0.5; nCdc20 = 8; kDeCaCdc20 = 1.0;

  // Division at mitotic entry (MPF crosses high): reset to a G1 daughter.
  E_div: at (MPF > 1): Dna = 0, Rc = 1, pRc = 0, aRc = 0, iRc = 0, Rb = Rb + pRb, pRb = 0, Ce = Ce/2, Ca = Ca/2, E1 = E1/2, MPF = 0, preMPF = 0, Cdc20 = 0 ;
"""

# ---- HU -> replication fork speed (dNTP depletion slows forks) ----
HU_BLOCK = """
  // ===== v44: HU -> fork speed (hydroxyurea depletes dNTPs -> slower replication) =====
  HU = 0;                 // hydroxyurea dose (0 = none; 1 ~ 10 uM experiment)
  vmin_fork = 0.1;        // residual fork speed at saturating HU (>0 -> S finite, no arrest)
  KmHU_fork = 1.0;        // HU IC50 for fork slowing
  hHU_fork = 3;           // Hill coefficient
  vfork := vmin_fork + (1 - vmin_fork)*KmHU_fork^hHU_fork/(KmHU_fork^hHU_fork + HU^hHU_fork);
"""

_DNA_RXN_OLD = "Synthesis_of_DNA: aRc => aRc + Dna; Cell*kSyDna*aRc;"
_DNA_RXN_NEW = "Synthesis_of_DNA: aRc => aRc + Dna; Cell*kSyDna*vfork*aRc;"

# ---- EZH2 epigenetic layer + CyclinD (Cd) regulation ----
# Heldt's Cd (CyclinD mitogen drive on Rb) is a constant 0.65. v44 makes it dynamic and
# repressible by EZH2 (H3K27me3 at CCND1), so the central hypothesis closes:
#   longer S -> EZH2 accumulates more (gated on E2f x CDK2 activity) -> stronger Cd
#   repression -> slower Rb phosphorylation -> longer subsequent G1.
# EZH2 transcription tracks E2f x (CycE + CycA) -- the Heldt analogs of v43's E2F x (Me+Ma)
# gate -- so EZH2 integrates time spent in the high-CDK2 (S) window. EZH2 protein turnover
# (kDeEZ) is set for tau ~ a few hours so EZH2 ramps within a cycle (within-cycle gradient).
# Calibrated in simulations/v44_ezh2_calib.py. mitogen is a placeholder for the HH/MYCN
# drive (stage 12c). EZH2i toggles the EZH2->Cd feedback (1 = OFF).
EZH2_BLOCK = """
  // ===== v44: EZH2 epigenetic layer; CyclinD (Cd) repressed by EZH2 =====
  species EZH2m in Cell, EZH2 in Cell;
  EZH2m = 0.1; EZH2 = 0.5;

  mitogen = 1.0;            // HH/MYCN proliferative drive (placeholder; stage 12c)
  Cd_max = 1.30;           // CyclinD scale (sets baseline Cd ~0.65 given EZH2 repression)
  K_EZH2_Cd = 0.5;         // EZH2 repression half-constant on CyclinD
  EZH2i = 0;               // EZH2->Cd feedback toggle (1 = OFF)
  Cd := Cd_max*mitogen*(K_EZH2_Cd/(K_EZH2_Cd + EZH2*(1 - EZH2i)));

  kEZbas = 0.0003;         // basal EZH2 transcription
  kEZE2f = 0.010;          // E2f-driven EZH2 transcription
  K_E2f_EZ = 0.3;
  K_Ce_EZ = 0.5; K_Ca_EZ = 0.8; wCe = 0.5;   // CycE(S-onset)/CycA(S-G2) gate weights
  kDeEZm = 0.02;           // EZH2 mRNA turnover
  kTlEZ = 0.004;           // EZH2 translation
  // EZH2 protein is STABLE (slow turnover) and RESET BY DILUTION at division (halved in
  // E_div). So EZH2 INTEGRATES synthesis over the cycle: a longer S (gated synthesis on
  // longer) -> higher EZH2 peak. This is what makes "longer S -> more EZH2" quantitative.
  kDeEZ = 0.0008;          // EZH2 protein turnover (tau ~20h >> cycle -> integrator)

  EZH2_tx: => EZH2m; Cell*(kEZbas + kEZE2f*E2f/(K_E2f_EZ + E2f)*(wCe*Ce/(K_Ce_EZ + Ce) + (1 - wCe)*Ca/(K_Ca_EZ + Ca)));
  EZH2m_deg: EZH2m => ; Cell*kDeEZm*EZH2m;
  EZH2_tl: EZH2m => EZH2m + EZH2; Cell*kTlEZ*EZH2m;
  EZH2_deg: EZH2 => ; Cell*kDeEZ*EZH2;
"""


def _load_heldt():
    with open(_HELDT_PATH) as f:
        return f.read()


def build_model_v44(hu=None, with_ezh2=True):
    """Build v44 = Heldt 2018 core + mitotic switch + HU->fork-speed coupling.

    with_ezh2=True (default) also makes CyclinD (Cd) dynamic and EZH2-repressible and adds
    the EZH2 epigenetic layer. with_ezh2=False keeps Heldt's constant Cd=0.65 (core engine
    only). `hu` optionally overrides the default HU=0 in the string.
    """
    m = _load_heldt()
    if _DNA_RXN_OLD not in m:
        raise RuntimeError("Heldt DNA-synthesis reaction not found in expected form.")
    if "\nend" not in m:
        raise RuntimeError("Heldt model 'end' marker not found.")

    # 1. HU-scale the replication fork flux
    m = m.replace(_DNA_RXN_OLD, _DNA_RXN_NEW)
    # 2. inject mitotic switch + HU blocks
    blocks = MITOSIS_BLOCK + HU_BLOCK

    if with_ezh2:
        # Make Cd dynamic (assignment rule) + add EZH2: remove Cd's const + init so the
        # EZH2 block's `Cd := ...` assignment rule governs it.
        if "const Cell, Cd, Skp2," not in m:
            raise RuntimeError("Heldt const declaration for Cd not in expected form.")
        m = m.replace("const Cell, Cd, Skp2,", "const Cell, Skp2,")
        m = m.replace("\n  Cd = 0.65;", "")   # drop the constant init; Cd := ... governs it
        blocks += EZH2_BLOCK
        # EZH2/EZH2m are diluted 2x at division (stable protein -> integrates S-duration)
        blocks = blocks.replace("MPF = 0, preMPF = 0, Cdc20 = 0 ;",
                                "MPF = 0, preMPF = 0, Cdc20 = 0, EZH2 = EZH2/2, EZH2m = EZH2m/2 ;")

    m = m.replace("\nend", blocks + "\nend")
    if hu is not None:
        m = m.replace("HU = 0;", f"HU = {hu};")
    return m


if __name__ == "__main__":
    m = build_model_v44()
    print("v44 model built:", len(m), "chars")
    for tok in ["MPF", "Chk1 := aRc", "Cdc25a", "Wee1a", "vfork :=", "kSyDna*vfork*aRc",
                "E_div: at (MPF > 1)"]:
        assert tok in m, f"missing {tok}"
    print("mitotic switch + HU fork coupling present.")
    import tellurium as te
    rr = te.loada(m)
    res = rr.simulate(0, 4000, 8000, selections=["time", "Dna", "MPF"])
    from scipy.signal import find_peaks
    import numpy as np
    pk, _ = find_peaks(res["MPF"], prominence=0.2, distance=50)
    print(f"sanity: {len(pk)} mitoses, Dna in [{res['Dna'].min():.2f},{res['Dna'].max():.2f}]")
