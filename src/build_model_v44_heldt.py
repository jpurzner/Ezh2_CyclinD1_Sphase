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
  g2gate := Dna^nG2/(KG2^nG2 + Dna^nG2);               // ~1 only when replication ~complete -> G2
  Cdc25a := (a25 + (1 - a25)*MPF^nMpf/(KmMpf^nMpf + MPF^nMpf))/(1 + wChk*Chk1);
  Wee1a  := (aWee + (1 - aWee)*KmMpf^nMpf/(KmMpf^nMpf + MPF^nMpf))*(1 + wChkW*Chk1);

  // CyclinB-CDK1 is synthesized (inactive, Tyr15-P = preMPF) only as replication completes
  // (g2gate on Dna), so it BUILDS in G2 rather than being pre-stocked during S -> a real
  // multi-hour G2. Cdc25/Wee1 hysteresis then gives sharp mitotic entry once it accumulates.
  SynCycB: => preMPF; Cell*kSyCb*g2gate;               // CycB-CDK1 made (inactive) in G2
  Wee1phos: MPF => preMPF; Cell*kWee*Wee1a*MPF;         // Tyr15 phosphorylation (inactivate)
  Cdc25dephos: preMPF => MPF; Cell*k25*Cdc25a*preMPF;   // Tyr15 dephosphorylation (activate)
  DegMPF: MPF => ; Cell*(kDeCbBas + kDeCb*Cdc20)*MPF;
  DegpreMPF: preMPF => ; Cell*(kDeCbBas + kDeCb*Cdc20)*preMPF;
  Cdc20act: => Cdc20; Cell*kaCdc20*(MPF^nCdc20/(KmCdc20^nCdc20 + MPF^nCdc20))*(1 - Cdc20);
  Cdc20inact: Cdc20 => ; Cell*kiCdc20*Cdc20;
  DegCycACdc20: Ca => ; Cell*kDeCaCdc20*Cdc20*Ca;       // mitotic CycA destruction

  kSyCb = 0.020; kDeCb = 0.08; kDeCbBas = 0.004;       // kSyCb=0.020 -> G2+M ~3.3h (within tol of the
                                                        // ~2.5h direct anchor; shortening it just shifts time into G0+G1 via the growth gate -> worse 2N count)
  kWee = 0.4; k25 = 0.7; KmMpf = 0.35; nMpf = 4;
  a25 = 0.02; aWee = 0.1;
  wChk = 12; wChkW = 4; jChk = 0.03;
  KG2 = 0.85; nG2 = 6;                                  // G2 gate: CycB synthesis ramps as Dna->1
  kaCdc20 = 0.3; kiCdc20 = 0.12; KmCdc20 = 0.5; nCdc20 = 8; kDeCaCdc20 = 1.0;

  // Division at mitotic entry (MPF crosses high): reset to a TRUE G1 daughter.
  // CyclinA (Ca) and CyclinE (Ce) are set low (mitotic APC/SCF degradation), so the daughter
  // starts Rb-hypophosphorylated with low CDK2 -> G1 length is the slow CyclinD-driven Rb
  // phosphorylation + cyclin rebuild time (otherwise inherited Ca/2 fires the R-point instantly).
  // Daughter is born with HIGH p27 (P21_div) + low Skp2 -> transient G0 (p27-positive / phospho-Rb-
  // negative) until the Skp2-p27-E2F feedforward commits it. (G0 = phospho-Rb(Ser807/811)- OR p27+.)
  Ca_div = 0.30; Ce_div = 0.10; MPF_div = 0.5; P21_div = 0.6;   // Ca_div=0.30 keeps GNP & MB cycling
  E_div: at (MPF > MPF_div): Dna = 0, Rc = 1, pRc = 0, aRc = 0, iRc = 0, Rb = Rb + pRb, pRb = 0, P21 = P21_div, CeP21 = 0, CaP21 = 0, Skp2 = 0.05, Ce = Ce_div, Ca = Ca_div, E1 = E1/2, MPF = 0, preMPF = 0, Cdc20 = 0 ;
"""

# ---- HU -> replication fork speed (dNTP depletion slows forks) ----
HU_BLOCK = """
  // ===== v44: HU -> fork speed (hydroxyurea depletes dNTPs -> slower replication) =====
  HU = 0;                 // hydroxyurea dose (0 = none; 1 ~ 10 uM experiment)
  vmin_fork = 0.1;        // residual fork speed at saturating HU (>0 -> S finite, no arrest)
  KmHU_fork = 0.2;        // HU IC50 (wide-search baked, was 0.35)
  hHU_fork = 3;           // Hill coefficient
  vfork := vmin_fork + (1 - vmin_fork)*KmHU_fork^hHU_fork/(KmHU_fork^hHU_fork + HU^hHU_fork);
"""

_DNA_RXN_OLD = "Synthesis_of_DNA: aRc => aRc + Dna; Cell*kSyDna*aRc;"
_DNA_RXN_NEW = "Synthesis_of_DNA: aRc => aRc + Dna; Cell*kSyDna*vfork*aRc;"

# ---- Structural redesign #2: cell-growth-gated restriction point ----
# Cell mass grows exponentially and halves at division (size homeostasis). A SIZE GATE is
# placed on S-ENTRY (origin firing, Phosphorylation_priming_of_replication_complexes): origins
# fire only once mass >= a critical size. Commitment (E2f release) STILL requires mitogen
# (CyclinD), so low-CyclinD cells (GNP without SHH, or +HHi) never commit and stay quiescent --
# HH-dependence preserved. Committed cells then WAIT for size before replicating, so G1 length =
# the GROWTH TIME, DECOUPLED from CyclinD level: high CyclinD (MB) no longer collapses G1 to ~0,
# letting the model have MB with high CyclinD1 AND a substantial G1 (reconciling the
# phase-proportion and HH between-condition calibrations).
#
# MECHANISTIC GROUNDING -- the size gate IS the Rb-DILUTION / titration restriction point
# (Zatulovskiy et al. 2020 Science 369:466; cf. the v45 successor's explicit Rb-dilution timer):
# Rb is made at a ~constant, size-INDEPENDENT amount per cell, so its CONCENTRATION [Rb] = tRb/mass
# is titrated DOWN by growth; commitment / S-entry are licensed once [Rb] dilutes below a threshold
# that CyclinD-CDK4/6 can overcome. Heldt's total Rb is CONSERVED (tRb := Rb+pRb+RbE2f, no synth/deg;
# the division event converts pRb->Rb without halving the total), so [Rb] is exactly proportional to
# 1/mass and oscillates tRb -> tRb/2 over the cycle. Hence the mass gate mass^n/(M^n+mass^n) is
# ALGEBRAICALLY the Rb-dilution gate Krb^n/(Krb^n + [Rb]^n) up to a reparametrization (M = tRb/Krb).
# So the "size checkpoint" here is not a postulated absolute-size sensor -- it is the molecular Rb
# titration. (To make Rb-dilution behave DIFFERENTLY from a mass gate one must break Rb conservation
# with size-independent, CONDITION-REGULATED Rb synthesis -- currently unconstrained by data.)
GROWTH_BLOCK = """
  // ===== v44 structural #2: cell-growth-gated restriction point (size gate on S-entry) =====
  species mass in Cell;
  mass = 1.0;
  mu = 0.0005;            // specific growth rate (1/min); ~ ln2/period for size homeostasis (~22h)
  M_size = 2.5;           // critical cell size for S-entry (origin firing)
  M_commit = 1.05;        // critical cell size for COMMITMENT (G0->G1; Skp2-p27 feedforward fires) ->
                          // transient G0 ~20% (MB); cell grows in G0 (p27 high) until mass>=M_commit
  n_size = 6;             // steepness of the size gates
  size_gate := mass^n_size/(M_size^n_size + mass^n_size);
  commit_gate := mass^n_size/(M_commit^n_size + mass^n_size);
  Growth: => mass; Cell*mu*mass;
"""
_FIRE_OLD = ("Phosphorylation_priming_of_replication_complexes: Rc => pRc; "
             "Cell*((kPhRc*(Ce + Ca)^n/(jCy^n + (Ce + Ca)^n))*Rc);")
_FIRE_NEW = ("Phosphorylation_priming_of_replication_complexes: Rc => pRc; "
             "Cell*((kPhRc*(Ce + Ca)^n/(jCy^n + (Ce + Ca)^n))*size_gate*Rc);")

# ---- Structural redesign #3: Skp2-p27 feedforward restriction-point switch ----
# Heldt's Skp2 is constant; here it becomes a dynamic E2F target degraded by APC/C-Cdh1 (C1).
# G0: E2F off -> low Skp2 synthesis AND active Cdh1 (C1) degrades Skp2 -> Skp2 low -> p27 (P21) NOT
# degraded -> p27 stays HIGH (constitutive), CDK2 inhibited, Rb hypophosphorylated, E2F off -- a
# self-reinforcing, p27-POSITIVE transient G0. CyclinD partially phosphorylates Rb -> some E2F ->
# Skp2 rises -> p27 degraded -> CyclinE/CDK2 active -> phosphorylates Cdh1 OFF (C1->pC1) -> Skp2 no
# longer destroyed -> more Skp2 -> p27 gone -> full Rb-P -> more E2F. This Skp2-p27-Rb-E2F feedforward
# IS the bistable R-point and yields the transient G0 (high p27 / phospho-Rb-negative) seen in all
# GNP and MB cells. p27 is reset HIGH at division (P21_div).
SKP2_BLOCK = """
  // ===== v44 structural #3: Skp2-p27 feedforward R-point (Skp2 = dynamic E2F target, Cdh1-degraded) =====
  kSySkp2 = 0.02;          // E2F-driven Skp2 transcription (Skp2 is a direct E2F target)
  kSySkp2bas = 0.002;      // basal Skp2 synthesis
  kDeSkp2C1 = 1.0;         // APC/C-Cdh1 (C1)-mediated Skp2 degradation (keeps Skp2 low in G0/G1)
  kDeSkp2bas = 0.02;       // basal Skp2 turnover
  Skp2_synthesis: => Skp2; Cell*(kSySkp2bas + kSySkp2*E2f);
  Skp2_degradation_Cdh1: Skp2 => ; Cell*kDeSkp2C1*C1*Skp2;
  Skp2_decay: Skp2 => ; Cell*kDeSkp2bas*Skp2;
"""

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
  kEZbas = 0.00027; kEZE2f = 0.022; K_E2f_EZ = 0.3; Kez_cd = 2.94;   // wide-search baked (EZH2 mitogen-dose term)
  K_Ce_EZ = 0.5; K_Ca_EZ = 0.8; wCe = 0.582;       // CycE(S-onset)/CycA(S-G2) gate weights (wCe baked)
  kDeEZm = 0.02; kTlEZ = 0.004; kDeEZ = 0.00015;   // stable EZH2 -> integrates S-duration (kDeEZ baked)
  // EZH2 is a Rb-E2f target driven by CyclinD1-CDK4/6: synthesis is cycle-gated (E2f x CycE/CycA, peaks S/G2)
  // AND MITOGEN-DOSE dependent (the *Cd/(Kez_cd+Cd) factor). The Cd term reproduces the dose-dependent EZH2
  // increase over a WIDE rShh range (Fig 4H) and the MB/GNP=2.05x (Fig 4J), while SATURATING (Kez_cd=2) so the
  // ~7x CyclinD1 fold compresses to ~2x EZH2. Without it EZH2 only tracks the binary commitment (MB/GNP~1.1).
  EZH2_tx: => EZH2m; Cell*(kEZbas + kEZE2f*E2f/(K_E2f_EZ + E2f)*(wCe*Ce/(K_Ce_EZ + Ce) + (1 - wCe)*Ca/(K_Ca_EZ + Ca))*Cd/(Kez_cd + Cd));
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
  species Gli_rep = 0.1, Gli_act = 0.5, Gli1_mRNA = 1.0, Gli1 = 1.4, Gli1_epi = 0.0001;
  species MYCN = 0.4, Cd_mRNA = 2.6, Cd in Cell;
  Cd = 0.70;
  // Ptch1_copy_number is the FUNCTIONAL Ptch1 fraction (gene dosage x functional competence):
  //   GNP = 1.0 (two wt alleles), Ptch1+/- = 0.5, MB-with-LOH = 0.1 (residual / non-functional).
  //   It gates Ptch1's REPRESSION of Smo (function), NOT Ptch1 transcription (production).
  Ptch1_copy_number = 1.0; MYCN_amplification = 1.0;

  // Ptch1 is a Gli TARGET: transcription = basal floor + Gli-induced (closes the canonical HH
  // negative feedback Gli->Ptch1-|Smo-|Gli). In MB the induced Ptch1 is non-functional (low
  // Ptch1_copy_number) so the loop is BROKEN -> constitutive Gli AND high Ptch1 mRNA (SHH-MB marker).
  k_Ptch1_basal = 0.1493; k_Ptch1_Gli = 1.873; K_Gli_Ptch = 0.3978;
  k_Ptch1_mRNA_deg = 0.8; k_Ptch1_translation = 1.0; k_Ptch1_deg = 0.5;
  k_SHH_Ptch_bind = 5.0; k_SHH_Ptch_release = 0.1; k_SHH_Ptch_deg = 0.8;
  // HH/MYCN->CyclinD1 DE-SATURATED + refit WITH the Gli->Ptch1 negative feedback (v44_recalibrate_gli.py):
  // Gli1 MB/GNP 7.2x, CyclinD1 6.8x, vismo crashes MB CyclinD1 to 0.142 of MB (~1.0x a cycling GNP).
  // Gli supplies ~86% of MB CyclinD1; Mycn/basal the HHi-resistant residual. Used WITH the saturating
  // CyclinD1->Rb drive (with_cd_sat). The Gli->Ptch1 loop (k_Ptch1_basal/Gli, K_Gli_Ptch above) damps
  // GNP Gli (transient overshoot, adaptive) and is BROKEN in MB (low functional Ptch1 -> constitutive Gli).
  k_Smo_act = 1.573; k_Smo_inact = 1.2; K_Ptch_Smo = 0.1025;
  k_Gli_rep_to_act = 3.0; k_Gli_act_to_rep = 1.5; K_Smo_Gli_switch = 0.8879;
  Vmax_Gli1_tx = 1.27; K_Gli_act_Gli1 = 1.182; n_Gli_act = 2; K_Gli_rep_Gli1 = 0.4; n_Gli_rep = 2;
  k_Gli1_mRNA_deg = 0.8; k_Gli1_translation = 1.2; k_Gli1_deg = 0.8;
  // GLI1 autoregulation via a SLOW epigenetic memory (Gli1_epi, a 0..1 chromatin state at the Gli
  // locus). Gli1_epi is CHARGED by Smo activity (the drug target; Hill K_Gli1_auto, n_Gli1_auto),
  // gated to BROKEN feedback (1 - Ptch1_copy_number) so it is GNP-silent and MB-active, and DISCHARGES
  // slowly (k_epi_off). It feeds GLI1 TRANSCRIPTION Smo-independently (k_Gli1_auto*Gli1_epi below).
  // Charging from Smo (not Gli1) makes it a non-bistable capacitor: vismo blocks Smo -> charging stops
  // -> the memory discharges, so vismo on cycling MB leaves a Gli1 RESIDUAL that decays over ~a day.
  // WHY it feeds Gli1 (not Gli_act): MYCN floors CyclinD1 but has NO path to Gli1, yet the data show a
  // large Gli1 residual after vismo (MB_GDC0449 Gli1 413 / MB 17881 = 2.3% at 24h, 41x GNP+vismo's 10)
  // -- a Gli-INTRINSIC residual. Feeding Gli1 reproduces that residual WITHOUT inflating CyclinD1/
  // proliferation (the residual Gli1 is small, so CyclinD1 after vismo stays MYCN-floored ~ pRb 1/4).
  // The no-memory model gives Gli1~0 by 24h, so ONLY the memory reproduces the residual.
  // k_epi_off=0.0008 set by the 24h timepoint: model Gli1(24h)/baseline ~0.035 (data 0.023; same order,
  // decay beyond 24h unmeasured). k_Gli1_auto=0 = legacy (off).
  k_Gli1_auto = 0.012; K_Gli1_auto = 0.4; n_Gli1_auto = 4; k_epi_on = 0.02; k_epi_off = 0.0008;
  // Re-attribution: with the memory ON, reduce the Smo-driven Gli1 transcription in MB (broken feedback)
  // so the memory SUPPLIES part of MB's constitutive Gli1 rather than inflating baseline. Gates by
  // (Ptch1_copy_number + (1 - Ptch1_copy_number)*g_smo_Gli1_broken): GNP (copy=1) unchanged.
  // 1.0 = legacy (no re-attribution); 0.873 = calibrated (keeps MB Gli1 6.9x, Cd 7.58x with memory on).
  g_smo_Gli1_broken = 0.873;
  k_Cd_tx_basal = 0.483; k_Cd_tx_Gli_max = 59.2; K_Gli_act_CycD = 0.4568; K_Gli_rep_CycD = 0.3;   // wide-search baked (basal, Gli_max)
  k_Cd_mRNA_deg = 0.8;
  k_MYCN_synth_basal = 0.3; k_MYCN_synth_Gli = 0.102; K_Gli_MYCN = 0.5; k_MYCN_deg = 1.0;
  k_Cd_tx_MYCN = 21.66; K_MYCN_Cd = 1.655; n_MYCN_Cd = 3.658;   // wide-search baked (was 35.22)
  K_EZH2_repression = 0.539;   // wide-search baked (was 0.75)
  k_Cd_translation = 0.801; k_Cd_deg = 1.0;   // wide-search baked (was 0.75). Cd protein scale: GNP (and MB+HHi == cycling-GNP level by
  // the data) must CLEANLY clear the cycling threshold. The desaturated Gli->Cd recalibration to
  // MB_GDC0449 dropped GNP Cd toward the bistable knife-edge (0.4 hysteretic, 0.5/0.65 left MB+HHi
  // numerically on the threshold and crashing); 0.8 puts GNP/MB+HHi clearly above it. Consistent with
  // the data: MB+HHi Mki67 is high, so vismo-treated MB keeps proliferating (no single-cell arrest).

  Ptch1_transcription: => Ptch1_mRNA; k_Ptch1_basal + k_Ptch1_Gli*Gli_act^n_Gli_act/(K_Gli_Ptch^n_Gli_act + Gli_act^n_Gli_act);
  Ptch1_mRNA_degradation: Ptch1_mRNA => ; k_Ptch1_mRNA_deg*Ptch1_mRNA;
  Ptch1_translation: Ptch1_mRNA => Ptch1_mRNA + Ptch1_free; k_Ptch1_translation*Ptch1_mRNA;
  Ptch1_degradation: Ptch1_free => ; k_Ptch1_deg*Ptch1_free;
  SHH_Ptch_binding: Ptch1_free => SHH_Ptch; k_SHH_Ptch_bind*SHH*Ptch1_free;
  SHH_Ptch_release: SHH_Ptch => Ptch1_free; k_SHH_Ptch_release*SHH_Ptch;
  SHH_Ptch_degradation: SHH_Ptch => ; k_SHH_Ptch_deg*SHH_Ptch;
  Smo_activation: => Smo_active; k_Smo_act/(1 + Ptch1_free*Ptch1_copy_number/K_Ptch_Smo)*(1 - GDC0449);
  Smo_inactivation: Smo_active => ; k_Smo_inact*Smo_active;
  Gli_rep_to_act: Gli_rep => Gli_act; k_Gli_rep_to_act*Smo_active^2/(K_Smo_Gli_switch^2 + Smo_active^2)*Gli_rep;
  Gli_act_to_rep: Gli_act => Gli_rep; k_Gli_act_to_rep*(1 - Smo_active^2/(K_Smo_Gli_switch^2 + Smo_active^2))*Gli_act;
  Gli1_transcription: => Gli1_mRNA; Vmax_Gli1_tx*Gli_act^n_Gli_act/(K_Gli_act_Gli1^n_Gli_act + Gli_act^n_Gli_act)*(1 - Gli_rep^n_Gli_rep/(K_Gli_rep_Gli1^n_Gli_rep + Gli_rep^n_Gli_rep))*(Ptch1_copy_number + (1 - Ptch1_copy_number)*g_smo_Gli1_broken) + k_Gli1_auto*Gli1_epi;
  Gli1_epi_on:  => Gli1_epi; k_epi_on*(1 - Ptch1_copy_number)*Smo_active^n_Gli1_auto/(K_Gli1_auto^n_Gli1_auto + Smo_active^n_Gli1_auto)*(1 - Gli1_epi);
  Gli1_epi_off: Gli1_epi => ; k_epi_off*Gli1_epi;
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


# ---- Two-step Rb: mono- (CyclinD) and hyper- (CyclinE/A) phosphorylation ----
# Heldt has a single phospho-Rb (`pRb`); any kinase converts Rb->pRb AND releases E2f in one step,
# so CyclinD drives pRb to max within ~20 min of division (it conflates mono- and hyper-phospho-Rb).
# Biology (Narasimha 2014; Sanidas 2019): CyclinD-CDK4/6 MONO-phosphorylates Rb in early G1 (still
# binds/represses E2f); CyclinE-CDK2 HYPER-phosphorylates at the R-point, RELEASING E2f (commitment).
# Here: Rb -> Rbm (mono, CyclinD, size-gated) -> Rbh (hyper, CyclinE/A); E2f stays bound on Rb & Rbm,
# released only at the hyper step. The experimental phospho-Rb(Ser807/811) cycling marker maps to the
# hyper form, so `pRb := Rbh` (alias) -> G0 = pRb-negative now coincides with p27-positive. The mono
# state Rbm is the literal transient-G0 buffer (E2f-bound, growing) before the Skp2-p27 feedforward
# fires the hyper switch. Toggle: build_model_v44(with_two_step_rb=True). Re-uses existing rate
# constants (kPhRbCd mono; kPhRbCe/kPhRbCa hyper; kDpRb dephos; kAsRbE2f/kDsRbE2f binding).
def _apply_two_step_rb(m, with_growth):
    # Keep `pRb` as the HYPER form (its Heldt annotation is literally "hyperphosphorylated"), add only
    # `Rbm` (mono) + `RbmE2f` (mono.E2f). Reaction IDs with Heldt annotations are kept; the mono->hyper
    # and mono-E2f reactions are appended with new IDs. `pRb` = experimental phospho-Rb(Ser807/811)
    # cycling marker -> existing scripts that read pRb now get the committed/hyper form.
    # The cell-size gate sits on the CyclinD->p27 clearance (below); the mono-phosphorylation step is
    # left ungated (gating it stalls commitment and the cell over-grows). The transient G0 is the
    # short p27-high window after division before the Skp2-p27 feedforward fires.
    cdk_d = "kPhRbCd*Cd"
    p27_clear = "kDeP21Cd*Cd*commit_gate" if with_growth else "kDeP21Cd*Cd"
    m = m.replace(
        "species Rb in Cell, pRb in Cell, E2f in Cell, RbE2f in Cell, E1 in Cell;",
        "species Rb in Cell, Rbm in Cell, pRb in Cell, E2f in Cell, RbE2f in Cell, RbmE2f in Cell, E1 in Cell;")
    m = m.replace("tRb := Rb + pRb + RbE2f;",
                  "tRb := Rb + Rbm + pRb + RbE2f + RbmE2f;")
    m = m.replace("tE2f := E2f + RbE2f;", "tE2f := E2f + RbE2f + RbmE2f;")
    # mono-phosphorylation by CyclinD (keep IDs); add hyper-phosphorylation by CyclinE/A (releases E2f)
    m = m.replace(
        "Phosphorylation_of_Rb: Rb => pRb; Cell*((kPhRbCd*Cd + kPhRbCe*Ce + kPhRbCa*Ca)*Rb);",
        f"Phosphorylation_of_Rb: Rb => Rbm; Cell*({cdk_d}*Rb);\n"
        f"  Hyper_Phosphorylation_of_Rb: Rbm => pRb; Cell*((kPhRbCe*Ce + kPhRbCa*Ca)*Rbm);")
    m = m.replace(
        "Phosphorylation_Rb_in_Rb_E2F_complexes: RbE2f => pRb + E2f; Cell*((kPhRbCd*Cd + kPhRbCe*Ce + kPhRbCa*Ca)*RbE2f);",
        f"Phosphorylation_Rb_in_Rb_E2F_complexes: RbE2f => RbmE2f; Cell*({cdk_d}*RbE2f);\n"
        f"  Hyper_Phosphorylation_Rb_in_RbmE2f: RbmE2f => pRb + E2f; Cell*((kPhRbCe*Ce + kPhRbCa*Ca)*RbmE2f);")
    # dephosphorylation: hyper->mono->hypo (keep the Dephosphorylation_of_Rb ID for the mono step)
    m = m.replace(
        "Dephosphorylation_of_Rb: pRb => Rb; Cell*kDpRb*pRb;",
        "Dephosphorylation_of_Rb: Rbm => Rb; Cell*kDpRb*Rbm;\n"
        "  Dephosphorylation_of_pRb_hyper: pRb => Rbm; Cell*kDpRb*pRb;")
    # E2f binds mono-Rb too (still represses); E2f degradation in the mono complex
    m = m.replace(
        "Association_dissociation_of_Rb_and_E2F: Rb + E2f -> RbE2f; Cell*(kAsRbE2f*Rb*E2f - kDsRbE2f*RbE2f);",
        "Association_dissociation_of_Rb_and_E2F: Rb + E2f -> RbE2f; Cell*(kAsRbE2f*Rb*E2f - kDsRbE2f*RbE2f);\n"
        "  Association_dissociation_of_Rbm_and_E2F: Rbm + E2f -> RbmE2f; Cell*(kAsRbE2f*Rbm*E2f - kDsRbmE2f*RbmE2f);")
    m = m.replace(
        "Degradation_of_E2F_in_Rb_E2F_complexes: RbE2f => Rb; Cell*kDeE2f*RbE2f;",
        "Degradation_of_E2F_in_Rb_E2F_complexes: RbE2f => Rb; Cell*kDeE2f*RbE2f;\n"
        "  Degradation_of_E2F_in_RbmE2f: RbmE2f => Rbm; Cell*kDeE2f*RbmE2f;")
    m = m.replace("\n  pRb = 5;", "\n  pRb = 5;\n  Rbm = 0;\n  RbmE2f = 0;")
    # CyclinD-CDK4/6 clears/titrates p27 (canonical D-CDK4/6 -> p27 sequestration). With the two-step
    # Rb, CyclinD no longer releases E2f directly, so without this it can't escape the p27-CDK2 block
    # and the cell dead-locks in G0. This makes commitment fire when the cyclin D1/p27 ratio crosses
    # threshold (Fan-Meyer 2021) -- pRb(hyper) stays low until CyclinE/CDK2 is freed and fires it.
    m = m.replace("kDeP21 + kDeP21Cy*Skp2*(Ce + Ca) + kDeP21aRc*Cdt2*aRc",
                  f"kDeP21 + {p27_clear} + kDeP21Cy*Skp2*(Ce + Ca) + kDeP21aRc*Cdt2*aRc")
    m = m.replace("\n  kDpRb = 0.05;",
                  "\n  kDpRb = 0.05;\n  kDeP21Cd = 0.15;\n  kDsRbmE2f = 1.5;")
    return m


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


def build_model_v44(hu=None, with_ezh2=True, with_hh=True, with_growth=True,
                    with_skp2=True, with_two_step_rb=False, with_cd_sat=True,
                    with_h3k27_memory=False, params=None):
    """Build v44 = Heldt 2018 core + mitotic switch + HU->fork-speed coupling.

    with_ezh2=True (default): add the EZH2 epigenetic layer and make CyclinD (Cd) dynamic.
    with_hh=True (default; requires with_ezh2): drive Cd from the full Hedgehog/MYCN module
      (SHH, GDC0449, MYCN, Ptch1_copy_number inputs) with EZH2 repression of CyclinD1.
    with_hh=False: Cd is a placeholder (constant mitogen x EZH2 repression).
    with_ezh2=False: Heldt's constant Cd=0.65 (core engine only).
    with_growth=True (default): cell-growth-gated restriction point (mass-scaled CyclinE/A
      synthesis) so G1 length is the growth time, decoupled from CyclinD level.
    `hu` optionally overrides the default HU=0 in the string.
    """
    m = _load_heldt()
    if _DNA_RXN_OLD not in m:
        raise RuntimeError("Heldt DNA-synthesis reaction not found in expected form.")
    if "\nend" not in m:
        raise RuntimeError("Heldt model 'end' marker not found.")

    # 1. HU-scale the replication fork flux
    m = m.replace(_DNA_RXN_OLD, _DNA_RXN_NEW)
    # 1b. faster replication fork: S-phase ~3.5h (Heldt default 0.0093 gave ~10h, unrealistically long;
    #     ~5x brings S to the data DMSO S proportion ~15.7% of the cycle). G1/G0 fills the rest.
    m = m.replace("kSyDna = 0.0093;", "kSyDna = 0.044;")
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

    if with_two_step_rb:
        # split Rb->pRb into Rb->Rbm (mono, CyclinD, size-gated)->Rbh (hyper, CyclinE/A, releases E2f)
        m = _apply_two_step_rb(m, with_growth)

    if with_cd_sat and not with_two_step_rb:
        # CyclinD1 -> Rb drive SATURATES (CDK4/6 kinase activity is bounded): kPhRbCd*Cd is replaced by
        # kPhRbCd*Cd/(K_CdRb + Cd). This DECOUPLES CyclinD1 LEVEL (which tracks the RNA-seq, up to ~7x in
        # MB once the HH pathway is de-saturated) from the bounded cell-cycle DRIVE -- so high CyclinD1
        # no longer makes the engine stiff. (Linear kPhRbCd*Cd crashed the rescue conditions at MB levels.)
        # INK4 CDK4/6 inhibitors (p16 + p18) are COMPETITIVE -> they raise the CyclinD1 half-max for Rb
        # phosphorylation: kPhRbCd*Cd/(K_CdRb*(1 + p16 + p18) + Cd). Higher INK4 tone raises the CyclinD1
        # threshold to commit, but MORE CyclinD1 still overcomes it -> EZH2i (de-represses CyclinD1) CAN
        # rescue an INK4-braked arrest. DISTINCT from CDK4/6i (palbociclib), modeled as kPhRbCd=0 (a Vmax
        # block) which raising CyclinD1 CANNOT bypass -> NOT rescuable.
        #   p16 (Cdkn2a): H3K27me3-silenced in GNPs (p16=0), ~139x induced in MB (the dramatic, MB-specific
        #                 INK4) -> p16>0 in MB.
        #   p18 (Cdkn2c): the CONSTITUTIVE INK4 -- substantially expressed in GNPs (p18=0.4 baseline; known
        #                 GNP/MB INK4 in the literature, Shh-maintained) and ~3x higher in MB (p18=1.2).
        # K_CdRb=0.357 with the GNP p18=0.4 baseline gives the same GNP net half-max (0.357*1.4=0.5) as the
        # earlier p18-free K_CdRb=0.5 -> GNP behavior preserved; MB tone is p16+p18.
        if "kPhRbCd*Cd " not in m:
            raise RuntimeError("Rb-phosphorylation (kPhRbCd*Cd) not found in expected form.")
        # INK4 (p16/p18) AND CIP/KIP p27 (P21) competitively inhibit CDK4/6: p27 acts on CDK4/6 too,
        # not only CDK2 (w_p27 weights its CDK4/6 inhibition; 0 = legacy CDK2-only).
        m = m.replace("kPhRbCd*Cd ", "kPhRbCd*Cd/(K_CdRb*(1 + p16 + p18 + w_p27*P21) + Cd) ")   # both Rb reactions
        m = m.replace("\n  kDpRb = 0.05;", "\n  kDpRb = 0.05;\n  K_CdRb = 0.319;\n  p16 = 0.0;\n  p18 = 0.464;\n  w_p27 = 1.0;")  # K_CdRb, p18 wide-search baked

    if with_growth:
        # size gate on S-entry (origin firing) + mass growth + halving at division
        if _FIRE_OLD not in m:
            raise RuntimeError("Heldt origin-firing reaction not found in expected form.")
        m = m.replace(_FIRE_OLD, _FIRE_NEW)
        # size gate on COMMITMENT: CyclinD->Rb trigger waits for size -> transient growth-timed G0.
        # With two-step Rb, commit_gate is already baked into the mono step (_apply_two_step_rb).
        if not with_two_step_rb:
            cd_rb = "kPhRbCd*Cd/(K_CdRb*(1 + p16 + p18 + w_p27*P21) + Cd) +" if with_cd_sat else "kPhRbCd*Cd +"
            if cd_rb not in m:
                raise RuntimeError("Rb-phosphorylation (CyclinD term) not found in expected form.")
            m = m.replace(cd_rb, cd_rb[:-2] + "*commit_gate +")
        blocks += GROWTH_BLOCK
        blocks = blocks.replace("preMPF = 0, Cdc20 = 0",
                                "preMPF = 0, Cdc20 = 0, mass = mass/2")

    if with_skp2:
        # Skp2 becomes a dynamic E2F-target species (Cdh1-degraded) -> Skp2-p27 feedforward R-point.
        if ", Skp2," not in m:
            raise RuntimeError("Heldt Skp2 const declaration not found in expected form.")
        m = m.replace(", Skp2,", ",", 1)            # remove Skp2 from the const list
        m = m.replace("\n  Skp2 = 1;", "\n  Skp2 = 0.05;")   # low initial (G0 level)
        # raise CyclinD->Rb so the feedforward fires for GNP+SHH (Cd~0.5) but not GNP-SHH (Cd~0.28).
        # NB kPhRbCd is GLOBAL (it raises the MB threshold too, which breaks the MB+HHi+EZH2i rescue) -- to widen
        # the GNP-specific sub-threshold range, raise GNP's baseline CKI brake (p18 default) instead, see below.
        m = m.replace("kPhRbCd = 0.2;", "kPhRbCd = 0.35;")   # wide-search baked: raises GNP commitment threshold (wider sub-threshold range; rescue preserved by co-tuned CKIs)
        blocks += SKP2_BLOCK

    if with_two_step_rb:
        # division reset: dephosphorylate all Rb states back to hypo (the daughter is born
        # phospho-Rb-negative / E2f re-sequestered), replacing the single-pRb reset.
        blocks = blocks.replace(
            "Rb = Rb + pRb, pRb = 0",
            "Rb = Rb + Rbm + pRb, Rbm = 0, pRb = 0, RbE2f = RbE2f + RbmE2f, RbmE2f = 0")
        # two-step-Rb defaults: M_commit slightly > birth mass gives a short growth-timed p27-high G0
        # while keeping REGULAR (size-homeostatic) GNP cycling -- M_commit>=2.6 makes low-CyclinD1 GNP
        # cycle irregularly (period-2 size oscillation), so 2.2 trades G0 length for cycle regularity.
        blocks = blocks.replace("M_commit = 1.3;", "M_commit = 2.2;")
        blocks = blocks.replace("a25 = 0.02;", "a25 = 0.05;")
        blocks = blocks.replace("KG2 = 0.85;", "KG2 = 0.9;")

    m = m.replace("\nend", blocks + "\nend")

    if with_h3k27_memory and with_ezh2:
        # Explicit H3K27me3 mark at the CyclinD1 locus (epigenetic MEMORY). EZH2 deposits the mark
        # (deposition blocked by EZH2i), it is removed by demethylation (+ optional replication-dilution
        # via aRc). The MARK (H3K27_Cd), not EZH2 directly, represses CyclinD1 -> after EZH2i the
        # de-repression has a real timescale (~ln2/k_demeth_cd) instead of being instantaneous.
        # CALIBRATION: k_meth_cd = k_demeth_cd makes the mark TRACK EZH2 at steady state
        # (H3K27_Cd_ss = EZH2*(1-EZH2i), dilution off), so every settled-state validation target is
        # UNCHANGED -- only the transient de-repression kinetics differ. The mark is INHERITED through
        # division (not in the E_div reset) = the memory; replication-dilution (k_dil_cd>0, gated on the
        # S-phase fork signal aRc) is off by default (needs recalibration of K_EZH2_repression if used).
        h3k27 = (
            "\n  species H3K27_Cd in Cell; H3K27_Cd = 0.5;"
            "\n  k_meth_cd = 0.0010; k_demeth_cd = 0.0010; k_dil_cd = 0.010;   // CALIBRATED (JP): CyclinD1 de-represses"
            "\n  // clearly by ~24h, ~max by ~48h after EZH2i; BOTH active demethylation (k_demeth, breaks the"
            "\n  // arrested-cell deadlock) + replication-dilution (k_dil*aRc, accelerates once cycling resumes)."
            "\n  H3K27_methylation:   => H3K27_Cd; Cell*k_meth_cd*EZH2*(1 - EZH2i);"
            "\n  H3K27_demethylation: H3K27_Cd => ; Cell*k_demeth_cd*H3K27_Cd;"
            "\n  H3K27_dilution:      H3K27_Cd => ; Cell*k_dil_cd*aRc*H3K27_Cd;"
        )
        m = m.replace("\nend", h3k27 + "\nend")
        # repress CyclinD1 via the mark instead of EZH2 directly (HH transcription path)
        m = m.replace("K_EZH2_repression + EZH2*(1 - EZH2i)", "K_EZH2_repression + H3K27_Cd")

    if hu is not None:
        m = m.replace("HU = 0;", f"HU = {hu};")
    if params:
        m = _apply_overrides(m, params)
    return m


if __name__ == "__main__":
    m = build_model_v44()
    print("v44 model built:", len(m), "chars")
    for tok in ["MPF", "Chk1 := aRc", "Cdc25a", "Wee1a", "vfork :=", "kSyDna*vfork*aRc",
                "E_div: at (MPF > MPF_div)", "size_gate :=", "*size_gate*Rc", "mass = mass/2", "Skp2_synthesis:", "P21 = P21_div"]:
        assert tok in m, f"missing {tok}"
    print("mitotic switch + HU fork coupling present.")
    import tellurium as te
    import numpy as np
    rr = te.loada(m)
    res = rr.simulate(0, 12000, 24000, selections=["time", "Dna", "MPF"])
    dna = res["Dna"]
    divisions = int(np.sum((dna[:-1] > 0.9) & (dna[1:] < 0.1)))   # count replication resets
    print(f"sanity: {divisions} divisions over 12000 min, "
          f"Dna in [{dna.min():.2f},{dna.max():.2f}], MPF max {res['MPF'].max():.2f}")
