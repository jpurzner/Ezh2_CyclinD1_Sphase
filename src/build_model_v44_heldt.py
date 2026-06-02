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

  // Division at mitotic entry (MPF crosses high): reset to a TRUE G1 daughter.
  // CyclinA (Ca) and CyclinE (Ce) are set low (mitotic APC/SCF degradation), so the daughter
  // starts Rb-hypophosphorylated with low CDK2 -> G1 length is the slow CyclinD-driven Rb
  // phosphorylation + cyclin rebuild time (otherwise inherited Ca/2 fires the R-point instantly).
  Ca_div = 0.02; Ce_div = 0.10;
  E_div: at (MPF > 1): Dna = 0, Rc = 1, pRc = 0, aRc = 0, iRc = 0, Rb = Rb + pRb, pRb = 0, Ce = Ce_div, Ca = Ca_div, E1 = E1/2, MPF = 0, preMPF = 0, Cdc20 = 0 ;
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

# ---- EZH2 epigenetic layer (transcription + protein; the Cd repression lives downstream) ----
# EZH2 transcription tracks E2f x (CycE + CycA) -- the Heldt analogs of v43's E2F x (Me+Ma)
# gate. EZH2 protein is STABLE (slow turnover) and reset by DILUTION at division (halved in
# E_div), so it INTEGRATES synthesis over the cycle: a longer S -> higher EZH2. EZH2i toggles
# the EZH2->CyclinD1 feedback (1 = OFF / EZH2 inhibitor).
EZH2_CORE_BLOCK = """
  // ===== v44: EZH2 epigenetic layer =====
  species EZH2m in Cell, EZH2 in Cell;
  EZH2m = 0.1; EZH2 = 0.5;
  EZH2i = 0;               // EZH2->CyclinD1 feedback toggle (1 = OFF)
  kEZbas = 0.0003; kEZE2f = 0.010; K_E2f_EZ = 0.3;
  K_Ce_EZ = 0.5; K_Ca_EZ = 0.8; wCe = 0.5;        // CycE(S-onset)/CycA(S-G2) gate weights
  kDeEZm = 0.02; kTlEZ = 0.004; kDeEZ = 0.0003;   // stable EZH2 -> integrates S-duration
  // kDeEZ=0.0003 calibrated (v44_calibrate_ezh2.py): HU->EZH2-in-S boost 1.28x (target 1.31,
  // alt 1.22); transcript gradient S/G0=2.0, G2/G0=2.1 (Section D 1.8-2.5). Decoupled from gradient.
  EZH2_tx: => EZH2m; Cell*(kEZbas + kEZE2f*E2f/(K_E2f_EZ + E2f)*(wCe*Ce/(K_Ce_EZ + Ce) + (1 - wCe)*Ca/(K_Ca_EZ + Ca)));
  EZH2m_deg: EZH2m => ; Cell*kDeEZm*EZH2m;
  EZH2_tl: EZH2m => EZH2m + EZH2; Cell*kTlEZ*EZH2m;
  EZH2_deg: EZH2 => ; Cell*kDeEZ*EZH2;
"""

# ---- CyclinD placeholder (used when with_hh=False): constant mitogen x EZH2 repression ----
CD_PLACEHOLDER_BLOCK = """
  // CyclinD (Cd) placeholder drive (no HH module): constant mitogen x EZH2 repression
  mitogen = 1.0; Cd_max = 1.30; K_EZH2_Cd = 0.5;
  Cd := Cd_max*mitogen*(K_EZH2_Cd/(K_EZH2_Cd + EZH2*(1 - EZH2i)));
"""

# ---- Hedgehog + MYCN -> CyclinD1 (drives Heldt's Cd) ----
# Ported from v42 (build_model_v42_mycn): SHH->Ptch1->Smo->Gli->CyclinD1; MYCN drives CyclinD1
# GDC-independently; EZH2 represses CyclinD1 transcription. Signaling runs fast (quasi-static
# mitogen input) vs the ~hours cell cycle. Cd is now a dynamic species (Cd_mRNA -> Cd).
# NOTE: parameter MAGNITUDES are inherited from v42 and not yet rescaled to the Heldt
# minute/Cd~0.65 frame -- deferred to the whole-data calibration pass (k_Cd_translation/k_Cd_deg
# set the Cd scale). SHH, GDC0449 are boundary inputs; Ptch1_copy_number, MYCN_amplification params.
HH_MYCN_BLOCK = """
  // ===== v44: Hedgehog + MYCN -> CyclinD1 (drives Cd) =====
  // HH species initialized near the PROLIFERATING steady state so Cd is ~0.65 from t=0
  // (starting from low Gli/Cd lets the cell miss its restriction-point window -> false G0).
  species $SHH = 0.5, $GDC0449 = 0.0;
  species SHH_Ptch = 0.0, Ptch1_free = 0.3, Ptch1_mRNA = 0.6, Smo_active = 0.8;
  species Gli_rep = 0.1, Gli_act = 0.5, Gli1_mRNA = 1.0, Gli1 = 1.4;
  species MYCN = 0.4, Cd_mRNA = 2.6, Cd in Cell;
  Cd = 0.70;
  Ptch1_copy_number = 1.0; MYCN_amplification = 1.0;

  k_Ptch1_tx = 0.5; k_Ptch1_mRNA_deg = 0.8; k_Ptch1_translation = 1.0; k_Ptch1_deg = 0.5;
  k_SHH_Ptch_bind = 5.0; k_SHH_Ptch_release = 0.1; k_SHH_Ptch_deg = 0.8;
  k_Smo_act = 1.5; k_Smo_inact = 1.2; K_Ptch_Smo = 0.4;
  k_Gli_rep_to_act = 3.0; k_Gli_act_to_rep = 1.5; K_Smo_Gli_switch = 0.6;
  Vmax_Gli1_tx = 1.2; K_Gli_act_Gli1 = 0.3; n_Gli_act = 2; K_Gli_rep_Gli1 = 0.4; n_Gli_rep = 2;
  k_Gli1_mRNA_deg = 0.8; k_Gli1_translation = 1.2; k_Gli1_deg = 0.8;
  k_Cd_tx_basal = 0.3; k_Cd_tx_Gli_max = 4.0; K_Gli_act_CycD = 0.4; K_Gli_rep_CycD = 0.3;
  k_Cd_mRNA_deg = 0.8;
  k_MYCN_synth_basal = 0.3; k_MYCN_synth_Gli = 0.102; K_Gli_MYCN = 0.5; k_MYCN_deg = 1.0;
  k_Cd_tx_MYCN = 15.0; K_MYCN_Cd = 1.5; n_MYCN_Cd = 3;
  K_EZH2_repression = 0.5;
  k_Cd_translation = 0.26; k_Cd_deg = 1.0;        // Cd scale (~0.65 baseline; tune later)

  Ptch1_transcription: => Ptch1_mRNA; k_Ptch1_tx*Ptch1_copy_number;
  Ptch1_mRNA_degradation: Ptch1_mRNA => ; k_Ptch1_mRNA_deg*Ptch1_mRNA;
  Ptch1_translation: Ptch1_mRNA => Ptch1_mRNA + Ptch1_free; k_Ptch1_translation*Ptch1_mRNA;
  Ptch1_degradation: Ptch1_free => ; k_Ptch1_deg*Ptch1_free;
  SHH_Ptch_binding: Ptch1_free => SHH_Ptch; k_SHH_Ptch_bind*SHH*Ptch1_free;
  SHH_Ptch_release: SHH_Ptch => Ptch1_free; k_SHH_Ptch_release*SHH_Ptch;
  SHH_Ptch_degradation: SHH_Ptch => ; k_SHH_Ptch_deg*SHH_Ptch;
  Smo_activation: => Smo_active; k_Smo_act/(1 + Ptch1_free/K_Ptch_Smo)*(1 - GDC0449);
  Smo_inactivation: Smo_active => ; k_Smo_inact*Smo_active;
  Gli_rep_to_act: Gli_rep => Gli_act; k_Gli_rep_to_act*Smo_active^2/(K_Smo_Gli_switch^2 + Smo_active^2)*Gli_rep;
  Gli_act_to_rep: Gli_act => Gli_rep; k_Gli_act_to_rep*(1 - Smo_active^2/(K_Smo_Gli_switch^2 + Smo_active^2))*Gli_act;
  Gli1_transcription: => Gli1_mRNA; Vmax_Gli1_tx*Gli_act^n_Gli_act/(K_Gli_act_Gli1^n_Gli_act + Gli_act^n_Gli_act)*(1 - Gli_rep^n_Gli_rep/(K_Gli_rep_Gli1^n_Gli_rep + Gli_rep^n_Gli_rep));
  Gli1_mRNA_degradation: Gli1_mRNA => ; k_Gli1_mRNA_deg*Gli1_mRNA;
  Gli1_translation: Gli1_mRNA => Gli1_mRNA + Gli1; k_Gli1_translation*Gli1_mRNA;
  Gli1_degradation: Gli1 => ; k_Gli1_deg*Gli1;
  MYCN_synthesis: => MYCN; k_MYCN_synth_basal*MYCN_amplification + k_MYCN_synth_Gli*(Gli_act + Gli1)/(K_Gli_MYCN + Gli_act + Gli1);
  MYCN_degradation: MYCN => ; k_MYCN_deg*MYCN;
  CycD1_transcription: => Cd_mRNA; (k_Cd_tx_basal + k_Cd_tx_Gli_max*(Gli_act + Gli1)^n_Gli_act/(K_Gli_act_CycD^n_Gli_act + (Gli_act + Gli1)^n_Gli_act)*(K_Gli_rep_CycD^n_Gli_rep/(K_Gli_rep_CycD^n_Gli_rep + Gli_rep^n_Gli_rep)) + k_Cd_tx_MYCN*MYCN^n_MYCN_Cd/(K_MYCN_Cd^n_MYCN_Cd + MYCN^n_MYCN_Cd))*(K_EZH2_repression/(K_EZH2_repression + EZH2*(1 - EZH2i)));
  CycD1_mRNA_degradation: Cd_mRNA => ; k_Cd_mRNA_deg*Cd_mRNA;
  Cd_translation: Cd_mRNA => Cd_mRNA + Cd; k_Cd_translation*Cd_mRNA;
  Cd_degradation: Cd => ; k_Cd_deg*Cd;
"""


def _load_heldt():
    with open(_HELDT_PATH) as f:
        return f.read()


import re as _re


def _apply_overrides(model, overrides):
    """Substitute scalar parameter initial values (`name = value;`) in the Antimony string.
    Used for calibration sweeps over `const` kinetic parameters that cannot be set at runtime.
    """
    for name, val in overrides.items():
        pat = _re.compile(rf"(\b{_re.escape(name)}\s*=\s*)[-0-9.eE]+(\s*;)")
        new, n = pat.subn(rf"\g<1>{val}\g<2>", model, count=1)
        if n == 0:
            raise RuntimeError(f"override parameter not found: {name}")
        model = new
    return model


def build_model_v44(hu=None, with_ezh2=True, with_hh=True, params=None):
    """Build v44 = Heldt 2018 core + mitotic switch + HU->fork-speed coupling.

    with_ezh2=True (default): add the EZH2 epigenetic layer and make CyclinD (Cd) dynamic.
    with_hh=True (default; requires with_ezh2): drive Cd from the full Hedgehog/MYCN module
      (SHH, GDC0449, MYCN, Ptch1_copy_number inputs) with EZH2 repression of CyclinD1.
    with_hh=False: Cd is a placeholder (constant mitogen x EZH2 repression).
    with_ezh2=False: Heldt's constant Cd=0.65 (core engine only).
    `hu` optionally overrides the default HU=0 in the string.
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
        # Cd becomes dynamic: free it from Heldt's const + constant init.
        if "const Cell, Cd, Skp2," not in m:
            raise RuntimeError("Heldt const declaration for Cd not in expected form.")
        m = m.replace("const Cell, Cd, Skp2,", "const Cell, Skp2,")
        m = m.replace("\n  Cd = 0.65;", "")   # drop Heldt's constant Cd init
        blocks += EZH2_CORE_BLOCK
        blocks += HH_MYCN_BLOCK if with_hh else CD_PLACEHOLDER_BLOCK
        # EZH2/EZH2m are diluted 2x at division (stable protein -> integrates S-duration)
        blocks = blocks.replace("MPF = 0, preMPF = 0, Cdc20 = 0 ;",
                                "MPF = 0, preMPF = 0, Cdc20 = 0, EZH2 = EZH2/2, EZH2m = EZH2m/2 ;")

    m = m.replace("\nend", blocks + "\nend")
    if hu is not None:
        m = m.replace("HU = 0;", f"HU = {hu};")
    if params:
        m = _apply_overrides(m, params)
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
