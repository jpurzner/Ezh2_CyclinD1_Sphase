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
  f_commit_carry = 0.5672456517307928;   // Feature C (Spencer carryover): fraction of commitment state (phospho-Rb, CyclinE, low-p27) inherited across division. 0 = legacy hard-reset (every G1 crashes E2f -> looks like arrest); >0 lets a committed daughter keep E2f in G1 (immediate re-entry) while uncommitted daughters still reset toward G0.
  E_div: at (MPF > MPF_div): Dna = 0, Rc = 1, pRc = 0, aRc = 0, iRc = 0, Rb = Rb + (1 - f_commit_carry)*pRb, pRb = f_commit_carry*pRb, P21 = P21_div - f_commit_carry*(P21_div - P21), CeP21 = 0, CaP21 = 0, Skp2 = 0.05, Ce = Ce_div + f_commit_carry*Ce, Ca = Ca_div, E1 = E1/2, MPF = 0, preMPF = 0, Cdc20 = 0 ;
"""

# ---- HU -> replication fork speed (dNTP depletion slows forks) ----
HU_BLOCK = """
  // ===== v44: HU -> fork speed (hydroxyurea depletes dNTPs -> slower replication) =====
  HU = 0;                 // hydroxyurea dose (0 = none; 1 ~ 10 uM experiment)
  vmin_fork = 0.1;        // residual fork speed at saturating HU (>0 -> S finite, no arrest)
  KmHU_fork = 0.2;        // HU IC50 (wide-search baked, was 0.35)
  hHU_fork = 3;           // Hill coefficient
  vfork := vmin_fork + (1 - vmin_fork)*KmHU_fork^hHU_fork/(KmHU_fork^hHU_fork + HU^hHU_fork);
  // HU also blocks S-ENTRY (origin firing): the ATR-CHK1 replication-stress checkpoint suppresses
  // new origin firing under dNTP depletion, so HU-arrested cells hold at 2N/G1 (few enter S, G2
  // depletes) while cells already in S crawl on (fork slowdown -> long S, EZH2 accumulates).
  // fire_gate_HU = 1 EXACTLY at HU=0 -> every HU=0 result is unchanged.
  KmHU_fire = 0.25;       // HU IC50 for the origin-firing (S-entry) block (calibrated to MB55 HU phase folds)
  hHU_fire = 2;           // Hill coefficient
  fire_gate_HU := KmHU_fire^hHU_fire/(KmHU_fire^hHU_fire + HU^hHU_fire);
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
  M_commit = 1.2843242130809425;        // critical cell size for COMMITMENT (G0->G1; Skp2-p27 feedforward fires) ->
                          // transient G0 ~20% (MB); cell grows in G0 (p27 high) until mass>=M_commit
  n_size = 6;             // steepness of the size gates
  // ---- graded G1 / cell-cycle LENGTHENING (JP 2026-07-20): the cycle length is the mass-doubling time (ln2/mu),
  // so subthreshold mitogen (low CyclinD1) SLOWS growth -> longer Tc, GRADUALLY (no transient G0). This delays
  // S-phase / lengthens the low-E2F interval so the slow EZH2 response follows. mu_min_frac=1 => mu_eff=mu
  // (NO-OP, default unchanged). Set mu_min_frac<1 to engage lengthening at low mitogen (arrest only as Cd->0). ----
  mu_min_frac = 1;        // 1 = off (legacy). <1 = growth slows to mu*mu_min_frac at low CyclinD1 (Tc lengthens up to /mu_min_frac)
  K_g1len = 1.0;          // CyclinD1 midpoint at which lengthening half-engages
  n_g1len = 4;            // steepness of the mitogen dependence
  // ---- CELL-TYPE cycle length via CKIs (JP 2026-07-24): GNP(P7) cycles FAST (~16h); MB is SLOWER (~22h) because
  // higher INK4 (p16+p18) + transient G0s. mu_eff is DIVIDED by (1 + k_mu_cki*(p16+p18)) so high-CKI MB grows/cycles
  // slower than low-CKI GNP. k_mu_cki=0 => OFF (no-op, cycle = ln2/mu for all types; default unchanged). To engage a
  // GNP-16h/MB-22h split: raise base mu and set k_mu_cki>0 (GNP p18=0.464 -> ~16h; MB p16+p18=1.86 -> ~22h).
  k_mu_cki = 0;           // 0 = off (legacy, cell-type-invariant cycle). >0 = INK4-slowed MB cycle
  mu_eff := mu/(1 + k_mu_cki*(p16 + p18))*(mu_min_frac + (1 - mu_min_frac)*Cd^n_g1len/(K_g1len^n_g1len + Cd^n_g1len));
  size_gate := mass^n_size/(M_size^n_size + mass^n_size);
  commit_gate := mass^n_size/(M_commit^n_size + mass^n_size);
  Growth: => mass; Cell*mu_eff*mass;
"""
_FIRE_OLD = ("Phosphorylation_priming_of_replication_complexes: Rc => pRc; "
             "Cell*((kPhRc*(Ce + Ca)^n/(jCy^n + (Ce + Ca)^n))*Rc);")
_FIRE_NEW = ("Phosphorylation_priming_of_replication_complexes: Rc => pRc; "
             "Cell*((kPhRc*(Ce + Ca)^n/(jCy^n + (Ce + Ca)^n))*size_gate*fire_gate_HU*Rc);")

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
  // 2026-07-20: MITOGEN-DOSE Skp2 induction. Data (JP scRNA + MB table): Skp2 is HIGHER in MB (~2.3x GNP-P7) and
  // declines with differentiation -- but E2F-only + Cdh1-degradation makes MB Skp2 LOW (MB high Cdh1/p27). Skp2 is
  // induced by proliferative/mitogen signalling (high in SHH-MB), so add a CyclinD1-dose synthesis term. kSySkp2_Cd=0
  // => NO-OP (default unchanged).
  kSySkp2_Cd = 0;          // mitogen(CyclinD1)-dose Skp2 induction: >0 raises Skp2 in high-CyclinD1 (MB) cells
  K_Skp2_Cd = 2.0;         // CyclinD1 half-saturation for the dose term
  Skp2_synthesis: => Skp2; Cell*(kSySkp2bas + kSySkp2*E2f + kSySkp2_Cd*Cd/(K_Skp2_Cd + Cd));
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
  // p27 SPLIT (2026-07-15 transient-G0 reframe, docs/transient_g0_synthesis.md): P21 is the FUNCTIONAL
  // CDK-inhibitory NUCLEAR pool -- the only one that couples to the cycle (inhibits CDK2/CDK4/6, drives G0).
  // P27_total is a widespread NON-functional readout (cytoplasmic/phospho/G1 p27 that rides along with cycling and
  // Atoh1+ cells) for IHC comparison ONLY -- it drives NOTHING. This is why total-p27 IHC cannot measure the G0
  // fraction (p27 is widespread even in cycling cells); the functional G0 marker is pRb-hypophosphorylation.
  P27_wide = 0.30;                    // widespread non-functional p27 baseline (IHC-visible, CDK-inert)
  P27_total := P21 + P27_wide;        // total p27 (IHC-comparable readout); NOT coupled to the cycle
  kEZbas = 0.00043040856433872047; kEZE2f = 0.021216027072906932; K_E2f_EZ = 0.5; Kez_cd = 7.830608487520997;   // EZH2-stability re-fit (was 0.00027/0.022): co-fit with faster kDeEZ against all 9 EZH2-coupled validation targets
  kEZbas_Cd = 0.0013508978696408793;      // mitogen-dose-scaled but CYCLE-FLAT baseline transcription: carries the MB/GNP dose
                          // WITHOUT a within-cycle swing -> decouples EZH2 mRNA phase gradient (transcript S/G0)
                          // from the mitogen-dose ratio. Repression is via the H3K27me3 mark (integrates EZH2),
                          // so raising the flat baseline does not over-repress CyclinD1.
  K_Ce_EZ = 0.5; K_Ca_EZ = 0.8; wCe = 0.582;       // CycE(S-onset)/CycA(S-G2) gate weights (wCe baked)
  kDeEZm = 0.0098; kTlEZ = 0.015098862630302894; kDeEZ = 0.0011685177801161655;  // EZH2 t1/2 ~14h (was 0.00015=77h, under-constrained): pinned by HU-arrest in-S boost + G0-withdrawal IF; kTlEZ co-tuned to hold EZH2 level (EZi-fold/MB-GNP)
  kDeEZ_C1 = 0;  // APC/C-Cdh1 (C1)-gated EZH2 degradation: EZH2 is a Cdh1 substrate -> degraded in G0/G1 (C1 active), STABILIZED in S/G2/committed & HU-arrest (C1 phospho-off) so elongated/arrested S ACCUMULATES EZH2 (exp-verified). Default 0 = legacy single-rate turnover.
  // EZH2 is a Rb-E2f target driven by CyclinD1-CDK4/6: synthesis is cycle-gated (E2f x CycE/CycA, peaks S/G2)
  // AND MITOGEN-DOSE dependent (the *Cd/(Kez_cd+Cd) factor). The Cd term reproduces the dose-dependent EZH2
  // increase over a WIDE rShh range (Fig 4H) and the MB/GNP=2.05x (Fig 4J), while SATURATING (Kez_cd=2) so the
  // ~7x CyclinD1 fold compresses to ~2x EZH2. Without it EZH2 only tracks the binary commitment (MB/GNP~1.1).
  EZH2_tx: => EZH2m; Cell*(kEZbas + (kEZbas_Cd + kEZE2f*E2f/(K_E2f_EZ + E2f)*(wCe*Ce/(K_Ce_EZ + Ce) + (1 - wCe)*Ca/(K_Ca_EZ + Ca)))*Cd/(Kez_cd + Cd));
  EZH2m_deg: EZH2m => ; Cell*kDeEZm*EZH2m;
  EZH2_tl: EZH2m => EZH2m + EZH2; Cell*kTlEZ*EZH2m;
  EZH2_deg: EZH2 => ; Cell*(kDeEZ + kDeEZ_C1*C1)*EZH2;
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
# set the Cd scale). SHH, HHi are boundary inputs; Ptch1_copy_number, MYCN_amplification params.
HH_MYCN_BLOCK = """
  // ===== v44: Hedgehog + MYCN -> CyclinD1 (drives Cd) =====
  // HH species initialized near the PROLIFERATING steady state so Cd is ~0.65 from t=0
  // (starting from low Gli/Cd lets the cell miss its restriction-point window -> false G0).
  species $SHH = 0.5, $HHi = 0.0;
  species SHH_Ptch = 0.0, Ptch1_free = 0.3, Ptch1_mRNA = 0.6, Smo_active = 0.8;
  species Gli_rep = 0.1, Gli_act = 0.5, Gli1_mRNA = 1.0, Gli1 = 1.4, Gli1_epi = 0.0001;
  species MYCN = 0.4, Cd_mRNA = 2.6, Cd in Cell;
  Cd = 0.70;
  // CyclinD2 (JP 2026-07-28): a SEPARATE, weakly-Hh-responsive D-cyclin. Data (Chahin RNA-seq): CCND2 is the
  // DOMINANT D-cyclin (~5x CCND1 in GNP) and drops only ~35% under vismodegib (vs CCND1 ~85%), so it is a large
  // buffered CDK4/6 floor. Modeled as basal + a small Gli-driven part + an MB developmental elevation (Cd2_expr),
  // no EZH2 repression. It adds to the CyclinD1->Rb drive with weight w_Cd2. Defaults 0 (neutral) until engaged.
  species Cd2 in Cell; Cd2 = 0.0;
  k_Cd2_bas = 0.0; k_Cd2_Gli = 0.0; K_Cd2_Gli = 0.457; k_Cd2_deg = 1.0; Cd2_expr = 0.0; w_Cd2 = 0.0;
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
  // CILIUM cell-cycle DUTY CYCLE (JP 2026-07-21; hh-signalling-timing.md): ~94% of EGL GNP progenitors are ciliated
  // but cilia RESORB pre-mitotically (G2/M) -> Shh reception is phase-gated (on G1/S, off through G2/M), restored
  // after division. cilium_gate multiplies the SHH->Ptch1 reception. w_cilium=0 => OFF (default; gate=1, no change).
  w_cilium = 0; K_cil = 0.6; n_cil = 4;
  cilium_gate := 1 - w_cilium*Cb^n_cil/(K_cil^n_cil + Cb^n_cil);   // ->1 in G1/S (ciliated), ->(1-w) at high CyclinB (G2/M, cilium resorbed)
  // HH/MYCN->CyclinD1 DE-SATURATED + refit WITH the Gli->Ptch1 negative feedback (v44_recalibrate_gli.py):
  // Gli1 MB/GNP 7.2x, CyclinD1 6.8x, vismo crashes MB CyclinD1 to 0.142 of MB (~1.0x a cycling GNP).
  // Gli supplies ~86% of MB CyclinD1; Mycn/basal the HHi-resistant residual. Used WITH the saturating
  // CyclinD1->Rb drive (with_cd_sat). The Gli->Ptch1 loop (k_Ptch1_basal/Gli, K_Gli_Ptch above) damps
  // GNP Gli (transient overshoot, adaptive) and is BROKEN in MB (low functional Ptch1 -> constitutive Gli).
  k_Smo_act = 1.573; k_Smo_inact = 1.2; K_Ptch_Smo = 0.1025;
  // DELAYED Ptch1 negative feedback (JP 2026-07-21): the real Gli->Ptch1 TRANSCRIPTIONAL feedback acts over HOURS,
  // not the quasi-static minutes the fast Ptch1_free uses. Ptch1_slow accumulates with Gli_act on an hours timescale
  // and adds to the Ptch1 inhibition of Smo -> Gli/CyclinD1 OVERSHOOTS then DROPS (adaptation). Gated by
  // Ptch1_copy_number (BROKEN in MB, like the fast loop). k_ptch1_slow_on=0 => OFF (default, validation-preserving).
  species Ptch1_slow; Ptch1_slow = 0.0;
  k_ptch1_slow_on = 0; k_ptch1_slow_off = 0.006; w_ptch1_slow = 4.0;
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
  // large Gli1 residual after vismo (MB_HHi Gli1 413 / MB 17881 = 2.3% at 24h, 41x GNP+vismo's 10)
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
  // MYCN_expr (JP 2026-07-27): SHH-MB is NOT MYCN-amplified (that is Group-3 MB); MYCN is ELEVATED because it
  // is a Hedgehog/Gli target. MYCN_expr is a MB-specific, DEVELOPMENTALLY Gli-set elevated-expression term
  // (additive to the basal), 0 in GNP. Data (Chahin RNA-seq): MYCN rises with chronic Hh dose (wt GNP 3077 ->
  // Ptch+/- 4251 -> MB 8808) but is ACUTELY vismodegib-RESISTANT in MB (8808->7862, -11%) vs GNP (-20%) -- i.e.
  // developmentally locked, acutely Gli-independent. An acute Gli1-driven self-activation LATCH cannot realize
  // this (the buffer + the GNP/MB discrimination are the same knob, and GNP's own 20% drop puts it too close to
  // any threshold; ~50-candidate search + fine threshold mapping, 2026-07-27), so the lock is represented as a
  // developmental state (this MB-specific term), not an acute switch. Replaces the MYCN_amplification multiplier.
  MYCN_expr = 0.0;
  // deprecated (failed mechanism, kept default-inert): the acute Gli1-Hill + self-activation latch.
  n_MYCN_Gli = 2; k_MYCN_auto = 0.0; K_MYCN_auto = 0.6; n_MYCN_auto = 4;
  k_Cd_tx_MYCN = 21.66; K_MYCN_Cd = 1.655; n_MYCN_Cd = 3.658;   // wide-search baked (was 35.22)
  K_EZH2_repression = 0.4846191484325215;   // EZH2-DIRECT CyclinD1 repression. NB: DEAD in the with_h3k27_chain DEFAULT (repression there is the PRC2-occupancy Hill f0_prc2+(1-f0)/(1+(PRC2/K_prc2)^n); tune EZH2i via K_prc2/a_rw_prc2). Only LIVE in the legacy w_ezdir (memory/lumped-mark) branches.
  k_Cd_translation = 0.801; k_Cd_deg = 1.0;   // wide-search baked (was 0.75). Cd protein scale: GNP (and MB+HHi == cycling-GNP level by
  // the data) must CLEANLY clear the cycling threshold. The desaturated Gli->Cd recalibration to
  // MB_HHi dropped GNP Cd toward the bistable knife-edge (0.4 hysteretic, 0.5/0.65 left MB+HHi
  // numerically on the threshold and crashing); 0.8 puts GNP/MB+HHi clearly above it. Consistent with
  // the data: MB+HHi Mki67 is high, so vismo-treated MB keeps proliferating (no single-cell arrest).

  Ptch1_transcription: => Ptch1_mRNA; k_Ptch1_basal + k_Ptch1_Gli*Gli_act^n_Gli_act/(K_Gli_Ptch^n_Gli_act + Gli_act^n_Gli_act);
  Ptch1_mRNA_degradation: Ptch1_mRNA => ; k_Ptch1_mRNA_deg*Ptch1_mRNA;
  Ptch1_translation: Ptch1_mRNA => Ptch1_mRNA + Ptch1_free; k_Ptch1_translation*Ptch1_mRNA;
  Ptch1_degradation: Ptch1_free => ; k_Ptch1_deg*Ptch1_free;
  SHH_Ptch_binding: Ptch1_free => SHH_Ptch; k_SHH_Ptch_bind*SHH*cilium_gate*Ptch1_free;
  SHH_Ptch_release: SHH_Ptch => Ptch1_free; k_SHH_Ptch_release*SHH_Ptch;
  SHH_Ptch_degradation: SHH_Ptch => ; k_SHH_Ptch_deg*SHH_Ptch;
  Smo_activation: => Smo_active; k_Smo_act/(1 + (Ptch1_free + w_ptch1_slow*Ptch1_slow)*Ptch1_copy_number/K_Ptch_Smo)*(1 - HHi);
  Smo_inactivation: Smo_active => ; k_Smo_inact*Smo_active;
  Ptch1_slow_on:  => Ptch1_slow; k_ptch1_slow_on*Gli_act^2*(1 - Ptch1_slow);   // slow (hours) delayed accumulation, STEEP in Gli (near-0 at baseline -> transient overshoot preserved)
  Ptch1_slow_off: Ptch1_slow => ; k_ptch1_slow_off*Ptch1_slow;
  Gli_rep_to_act: Gli_rep => Gli_act; k_Gli_rep_to_act*Smo_active^2/(K_Smo_Gli_switch^2 + Smo_active^2)*Gli_rep;
  Gli_act_to_rep: Gli_act => Gli_rep; k_Gli_act_to_rep*(1 - Smo_active^2/(K_Smo_Gli_switch^2 + Smo_active^2))*Gli_act;
  Gli1_transcription: => Gli1_mRNA; Vmax_Gli1_tx*Gli_act^n_Gli_act/(K_Gli_act_Gli1^n_Gli_act + Gli_act^n_Gli_act)*(1 - Gli_rep^n_Gli_rep/(K_Gli_rep_Gli1^n_Gli_rep + Gli_rep^n_Gli_rep))*(Ptch1_copy_number + (1 - Ptch1_copy_number)*g_smo_Gli1_broken) + k_Gli1_auto*Gli1_epi;
  Gli1_epi_on:  => Gli1_epi; k_epi_on*(1 - Ptch1_copy_number)*Smo_active^n_Gli1_auto/(K_Gli1_auto^n_Gli1_auto + Smo_active^n_Gli1_auto)*(1 - Gli1_epi);
  Gli1_epi_off: Gli1_epi => ; k_epi_off*Gli1_epi;
  Gli1_mRNA_degradation: Gli1_mRNA => ; k_Gli1_mRNA_deg*Gli1_mRNA;
  Gli1_translation: Gli1_mRNA => Gli1_mRNA + Gli1; k_Gli1_translation*Gli1_mRNA;
  Gli1_degradation: Gli1 => ; k_Gli1_deg*Gli1;
  MYCN_synthesis: => MYCN; (k_MYCN_synth_basal + MYCN_expr) + k_MYCN_synth_Gli*(Gli_act + Gli1)/(K_Gli_MYCN + Gli_act + Gli1);
  MYCN_degradation: MYCN => ; k_MYCN_deg*MYCN;
  CycD1_transcription: => Cd_mRNA; (k_Cd_tx_basal + k_Cd_tx_Gli_max*(Gli_act + Gli1)^n_Gli_act/(K_Gli_act_CycD^n_Gli_act + (Gli_act + Gli1)^n_Gli_act)*(K_Gli_rep_CycD^n_Gli_rep/(K_Gli_rep_CycD^n_Gli_rep + Gli_rep^n_Gli_rep)) + k_Cd_tx_MYCN*MYCN^n_MYCN_Cd/(K_MYCN_Cd^n_MYCN_Cd + MYCN^n_MYCN_Cd))*(K_EZH2_repression/(K_EZH2_repression + EZH2*(1 - EZH2i)));
  CycD1_mRNA_degradation: Cd_mRNA => ; k_Cd_mRNA_deg*Cd_mRNA;
  Cd_translation: Cd_mRNA => Cd_mRNA + Cd; k_Cd_translation*Cd_mRNA;
  CycD2_synthesis: => Cd2; (k_Cd2_bas + Cd2_expr + k_Cd2_Gli*(Gli_act + Gli1)/(K_Cd2_Gli + Gli_act + Gli1));
  CycD2_degradation: Cd2 => ; k_Cd2_deg*Cd2;
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
def _apply_two_step_rb(m, with_growth, decouple_commit=False):
    # Keep `pRb` as the HYPER form (its Heldt annotation is literally "hyperphosphorylated"), add only
    # `Rbm` (mono) + `RbmE2f` (mono.E2f). Reaction IDs with Heldt annotations are kept; the mono->hyper
    # and mono-E2f reactions are appended with new IDs. `pRb` = experimental phospho-Rb(Ser807/811)
    # cycling marker -> existing scripts that read pRb now get the committed/hyper form.
    # The cell-size gate sits on the CyclinD->p27 clearance (below); the mono-phosphorylation step is
    # left ungated (gating it stalls commitment and the cell over-grows). The transient G0 is the
    # short p27-high window after division before the Skp2-p27 feedforward fires.
    # CyclinD-CDK4/6 kinase, competitively braked by INK4 (p16/p18) + p27 (saturating, bounded Vmax) --
    # applied to the mono-phosphorylation step (Rb=>Rbm) AND the p27-clearance (both done by ACTIVE CDK4/6).
    # 2026-07-15 (transient-G0 reframe, docs/transient_g0_synthesis.md): the p27-clearance is now gated on CDK4/6
    # ACTIVITY (the numerator carries kPhRbCd, the same active-kinase rate that phosphorylates Rb), NOT merely on
    # CyclinD LEVEL. Consequences, matching JP's CDK4/6i data: (i) CDK4/6i (kPhRbCd=0) stops p27 clearance -> p27
    # ACCUMULATES onto CDK2 -> the reversible p27-high G0 FILLS (was: p27 stayed ~0 under CDK4/6i, wrong);
    # (ii) baseline G0 EMERGES from the CDK4/6-INK4 balance -- GNP (high mitogen, no INK4) clears p27 fast -> G0~0;
    # MB (high p16/p18 braking CDK4/6) clears slowly -> a modest baseline G0 -- so the MB G0 comes from real INK4
    # biology, NOT an inflated birth-p27. The w_p27*P21 denominator self-brake (mutual CDK4/6<->p27 antagonism) is
    # kept. G0 is now read off pRb-HYPOphosphorylation (the functional switch), p27 is the supporting readout.
    # CyclinD1 + CyclinD2 both feed CDK4/6 -> Rb (Cd2 weighted by w_Cd2; w_Cd2=0 => CyclinD1-only, default-neutral)
    cdk_d = "kPhRbCd*(Cd + w_Cd2*Cd2)/(K_CdRb*(1 + p16 + p18 + w_p27*P21) + (Cd + w_Cd2*Cd2))"
    _clear_brake = "/(K_CdRb*(1 + p16 + p18 + w_p27*P21) + (Cd + w_Cd2*Cd2))"
    # decouple_commit (JP 2026-07-26): drop the growth/size gate (commit_gate) from p27-clearance so the
    # G0->G1 commitment is governed PURELY by the CyclinD1/CDKi balance (cdk_d) + birth-p27 (Spencer/Cappell),
    # NOT by growth rate. Lets a fast 16h core cycle keep a strong CyclinD1-vs-CDKi commitment threshold, and
    # makes MB's transient G0 (longer AVERAGE cycle) come from its high CDKi, not from a slower engine.
    _commit_gated = with_growth and not decouple_commit
    p27_clear = f"kDeP21Cd*kPhRbCd*(Cd + w_Cd2*Cd2){_clear_brake}*commit_gate" if _commit_gated else f"kDeP21Cd*kPhRbCd*(Cd + w_Cd2*Cd2){_clear_brake}"
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
    _p27_target = "kDeP21 + kDeP21Cy*Skp2*(Ce + Ca) + kDeP21aRc*Cdt2*aRc"
    assert m.count(_p27_target) >= 1, "two-step p27-clearance splice: target string not found (silent no-op)"
    m = m.replace(_p27_target,
                  f"kDeP21 + {p27_clear} + kDeP21Cy*Skp2*(Ce + Ca) + kDeP21aRc*Cdt2*aRc")
    m = m.replace("\n  kDpRb = 0.05;",
                  "\n  kDpRb = 0.05;\n  kDeP21Cd = 0.15;\n  kDsRbmE2f = 1.5;"
                  "\n  K_CdRb = 0.319;\n  p16 = 0.0;\n  p18 = 0.464;\n  w_p27 = 1.0;")
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
                    with_skp2=True, with_two_step_rb=True, with_cd_sat=True,
                    with_h3k27_memory=False, with_h3k27_dilution=True, with_prc2=True,
                    with_h3k27_chain=True, with_mother_g2=False, with_ezh2_conc=False,
                    with_mitogen_tracker=False, with_diff_gene=False, decouple_commit=True,
                    mycn_autoreg=False, with_cdki_species=True, with_p21_pip_degron=True,
                    with_p27_optionB=False, with_cdk6_gli=False, with_cd_hyper_escape=False,
                    with_mark_amplifier=False, with_proximal_distal=False,
                    params=None):   # 2026-07-30 BAKED: dynamic CDKI. 2026-08-01 BAKED: with_p21_pip_degron (un-map — restore inherited PIP-degron machinery p27->p21; species-correct re-cal, 30/32). 2026-08-05: with_p27_optionB (EXPLORATION, default OFF; Fan-Meyer p27-inhibitory CDK4/6 — buffered p27 = INACTIVE Cd, not counted in the Rb drive; needs re-cal before default)
    """Build v44 = Heldt 2018 core + mitotic switch + HU->fork-speed coupling.

    with_ezh2=True (default): add the EZH2 epigenetic layer and make CyclinD (Cd) dynamic.
    with_hh=True (default; requires with_ezh2): drive Cd from the full Hedgehog/MYCN module
      (SHH, HHi, MYCN, Ptch1_copy_number inputs) with EZH2 repression of CyclinD1.
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
    m = m.replace("kSyDna = 0.0093;", "kSyDna = 0.04265121846938828;")
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
        # mycn_autoreg (DEPRECATED, default-off): the acute Gli1-driven MYCN self-activation LATCH. It cannot
        # realize the developmental Hh-lock (buffer strength and GNP/MB discrimination are the same knob; see the
        # MYCN_expr note). Superseded by the additive MYCN_expr elevated-expression term (now the default). Flag
        # kept only so callers/env do not break; it is intentionally inert.
        _ = mycn_autoreg
        # EZH2/EZH2m dilution at division. DILUTION convention (with_ezh2_conc=False, DEFAULT since 2026-07-19):
        # EZH2 is DILUTED at division (discrete /2) so that a cell's EZH2 tracks its DIVISION HISTORY -- this is
        # what makes an ELONGATED / HU-ARRESTED S ACCUMULATE EZH2 (exp-verified by JP): dividing cells reset it,
        # arrested cells don't. Paired with LOW baseline kDeEZ (committed/S cells STABLE) + APC/Cdh1(C1)-gated
        # kDeEZ_C1 (clears EZH2 in G0). The concentration convention (with_ezh2_conc=True, default 2026-07-13..07-19)
        # preserved EZH2 across division so arrest could NOT accumulate it -> the S-elongation EZH2 boost was flat
        # (1.009); reverted here. kTlEZ stays at its string default (no concentration-convention x0.78 rescale).
        if not with_ezh2_conc:
            blocks = blocks.replace("MPF = 0, preMPF = 0, Cdc20 = 0 ;",
                                    "MPF = 0, preMPF = 0, Cdc20 = 0, EZH2 = EZH2/2, EZH2m = EZH2m/2 ;")
        else:
            # concentration-convention re-fit (refit_ezh2_conc.py grid best): dropping the halving inflates the
            # EZH2 level ~47% -> restore it with kTlEZ x0.78 (kDeEZ unchanged). Reaches 26/28. The one new miss is
            # the CyclinD1 GNP+HHi/GNP fold (0.29 vs noisy target 0.16) -- it sat at its passing EDGE (0.22) even
            # WITH the halving, and trades off against EZH2 G0/cycling within (kTlEZ,kDeEZ) (2-param structural
            # ceiling). All EZH2 phase ratios + levels hold. PROMOTED to default (JP 2026-07-13); the lost target
            # is a noisy HHi fold accepted per the loosen-noisy-targets policy. with_ezh2_conc=False = legacy 27/28.
            blocks = blocks.replace("kTlEZ = 0.015098862630302894", "kTlEZ = 0.011771")

    if with_two_step_rb:
        # split Rb->pRb into Rb->Rbm (mono, CyclinD, size-gated)->Rbh (hyper, CyclinE/A, releases E2f)
        m = _apply_two_step_rb(m, with_growth, decouple_commit)

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
        # division reset WITH Spencer carryover (f_commit_carry): the daughter INHERITS a fraction of the
        # mono (Rbm) + hyper (pRb) Rb states -- committed daughters exit mitosis phospho-Rb-high (Moser
        # 2018: CDK2inc daughters are pRb-high right after anaphase), uncommitted daughters reset toward
        # hypo Rb (G0). NOT a full dephosphorylation. Conserves tRb and tE2f. Operates on `m` (the E_div
        # event lives in the model header, not `blocks` -- the old blocks.replace silently no-op'd once
        # f_commit_carry rewrote the reset). NB the CKI carried at birth (P21 = p27 pool) is the DOMINANT
        # CIP/KIP in GNP/MB (p27, not p21). The E_div event lives in blocks (MITOSIS_BLOCK).
        blocks = blocks.replace(
            "Rb = Rb + (1 - f_commit_carry)*pRb, pRb = f_commit_carry*pRb",
            "Rb = Rb + (1 - f_commit_carry)*(pRb + Rbm), pRb = f_commit_carry*pRb, Rbm = f_commit_carry*Rbm, "
            "RbE2f = RbE2f + (1 - f_commit_carry)*RbmE2f, RbmE2f = f_commit_carry*RbmE2f")
        # two-step-Rb defaults: M_commit slightly > birth mass gives a short growth-timed p27-high G0
        # while keeping REGULAR (size-homeostatic) GNP cycling. (Was a no-op: target string is now the
        # baked single-step M_commit; STARTING value -- the two-step needs its own re-optimization.)
        blocks = blocks.replace("M_commit = 1.2843242130809425;", "M_commit = 2.2;")

    if with_mother_g2:
        # ⚰️ GRAVESTONE -- SETTLED NEGATIVE RESULT (2026-07, do not revive; git preserves it). This mother-G2
        # p27 integrator does NOT reproduce the Overton bifurcation: the birth-p27 floor P21_div + g*Sg2 can
        # only RAISE p27, so it saturates ~1 (+ a frequency confound). See memory v44-mother-g2-negative-result.
        # Kept default-off and unused; the working model uses a per-type IDENTITY birth-p27 (P21_DIV_MB).
        # PROTOTYPE (exploratory, opt-in, NOT promoted): the daughter's birth p27 is set by the MOTHER's
        # G2 mitogen history (Spencer 2013 R1 window; Min 2020 CyclinD-translation integrator) instead of
        # being the fixed parameter P21_div. Sg2 integrates a p27-inducing signal during the mother's G2
        # (Dna~1, pre-mitosis), proportional to the mitogen DEFICIT (low Cd -> more accumulation). At
        # division the daughter's effective birth p27 = P21_div + g_moth*Sg2, and Sg2 resets to 0 for the
        # daughter's own G2. -> low mother-G2 mitogen -> high birth p27 -> daughter goes transient-G0 /
        # CDK2low; the fate is a threshold on an INHERITED condition, not a within-cycle race. p27 (P21
        # pool) is the dominant CIP/KIP in GNP/MB. Sanity-check only; needs calibration + an ensemble to
        # reproduce the Overton graded->binary bifurcation.
        mg2 = (
            "\n  species Sg2 in Cell; Sg2 = 0;   // mother-G2 p27 integrator (inherited birth-p27 signal)"
            "\n  k_sg2 = 0.003; K_mit_g2 = 0.35; g_moth = 1.2; k_sg2_dec = 0.0004; n_g2 = 8;"
            "\n  Sg2_accumulate: => Sg2; Cell*k_sg2*(Dna^n_g2/(0.9^n_g2 + Dna^n_g2))*(K_mit_g2/(K_mit_g2 + Cd));"
            "\n  Sg2_decay: Sg2 => ; Cell*k_sg2_dec*Sg2;"
        )
        blocks += mg2
        blocks = blocks.replace(
            "P21 = P21_div - f_commit_carry*(P21_div - P21)",
            "P21 = (P21_div + g_moth*Sg2) - f_commit_carry*((P21_div + g_moth*Sg2) - P21)")
        blocks = blocks.replace("preMPF = 0, Cdc20 = 0", "preMPF = 0, Cdc20 = 0, Sg2 = 0")

    if with_mitogen_tracker:
        # ⚰️ GRAVESTONE -- SHELVED, WRONG SIGN FOR MB (2026-07, do not revive; git preserves it). A mitogen-
        # DEFICIT tracker forces MB (high mitogen: Cd ~12 vs GNP ~1.7) to be born p27-LOW, i.e. ~10-30x BELOW
        # GNP -- the OPPOSITE of the data (MB p27 transcript HIGHER). MB p27 is MYCN/INK4/identity-driven, not a
        # mitogen deficit, so routing inherited fate through Mtr is the wrong causal model for the paper's cell
        # type. The working model uses the per-type IDENTITY constant P21_DIV_MB (=1.39). Kept default-off/unused.
        # PHASE 1 of the daughter-transfer plan (docs/daughter_transfer_plan.md). The daughter's birth p27 is
        # SET (not the fixed P21_div, and NOT added to it -- that was the mother_g2 FLOOR trap) by a low-pass
        # tracker Mtr of the mother's mitogen DEFICIT. Mtr is a relaxation LEVEL-tracker (steady state = the
        # signal -> frequency-clean, fixing the mother_g2 confound) and is NOT reset at division -> both sisters
        # inherit it (Spencer ~98% concordance). P27set is a LOW-floored INCREASING function of the deficit Mtr:
        # high deficit = low mitogen -> high Mtr -> high birth p27 -> transient G0; high mitogen -> low Mtr ->
        # birth p27 ~ P21_lo -> immediate/CDK2inc. fc_p27 (small, DECOUPLED from the pRb-marker f_commit_carry)
        # keeps the effective floor ~ P21_lo so high-mitogen daughters actually reach the immediate basin.
        # Driver = Cd (protein) [open question]; the Str/p21 stress arm is Phase 3 (not here). Params are the
        # plan's STARTING points -- a single-lineage re-pin (GNP settled P27set ~ 0.6 through the loop) comes
        # after we see the Phase-1 interaction. Default OFF.
        mtr = (
            "\n  species Mtr in Cell; Mtr = 0.5;   // mitogen-DEFICIT low-pass tracker (NOT reset at division)"
            "\n  k_M = 0.0025; K_m = 2.5; n_m = 3;   // 1/k_M ~ 6.7h ~ mother R-window; K_m centers GNP (settled"
            "\n  //                                     Cd ~ 1.9, NOT the 0.70 init) on the sigmoid shoulder -> P27set ~ 0.6"
            "\n  mdef := K_m^n_m/(K_m^n_m + Cd^n_m);   // deficit: HIGH at low Cd (low mitogen)"
            "\n  Mtr_track: => Mtr; Cell*k_M*mdef;"
            "\n  Mtr_relax: Mtr => ; Cell*k_M*Mtr;   // net dMtr/dt = k_M*(mdef - Mtr)"
            "\n  P21_lo = 0.06; P21_hi = 1.10; K_S = 0.60; h_S = 4; fc_p27 = 0.10;"
            "\n  P27set := P21_lo + (P21_hi - P21_lo)*Mtr^h_S/(K_S^h_S + Mtr^h_S);   // INCREASING in deficit Mtr"
        )
        blocks += mtr
        blocks = blocks.replace(
            "P21 = P21_div - f_commit_carry*(P21_div - P21)",
            "P21 = (1 - fc_p27)*P27set + fc_p27*P21")

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

    if with_h3k27_dilution and with_ezh2 and not with_h3k27_memory:
        # Mean-field AUM H3K27me3 module at the Ccnd1 locus (Purzner spec + literature review, June 2026).
        # NOTE: this is now the DEFAULT repression (promoted June 2026); pass with_h3k27_dilution=False for
        # the legacy direct-EZH2 term, or with_h3k27_memory=True for the memory variant (which takes priority).
        # Mk = H3K27me3 occupancy over the ~7 kb / ~35-nucleosome bivalent Ccnd1 domain in [0,1]. This is
        # the MEAN-FIELD REDUCTION of the Berry-Dean-Howard 2017 A/U/M per-nucleosome model: one state Mk
        # (the U<->M arm with read-write); the "active/A" (H3K27ac) state's antagonism is supplied not by a
        # separate acetyl species but by the TRANSCRIPTION -> PRC2 reciprocal arm below (lit-review S3.3:
        # nascent RNA evicts PRC2 / inhibits its HMTase in real time). Well-mixed mean-field is justified at
        # N~35 because PRC2 read-write recruitment is non-local (Dodd 2007; lit-review S3.9).
        #
        # Continuous methylation (PRC2 catalytic-rate term, lit-review S3.6):
        #   * AUTOCATALYTIC read-write k_w*EZH2*Mk  (EED reads me3 -> allosterically activates EZH2)
        #   * de-novo nucleation floor k0  (accessory-driven, so a halved locus can reseed)
        #   * both on unmethylated substrate (1-Mk), both blocked by EZH2i (PRC2 inhibition)
        #   * RECIPROCAL ARM (1 - g*tx-Hill): ongoing Ccnd1 transcription (proxied by Cd_mRNA = nascent
        #     output) evicts PRC2 -> a double-negative => positive feedback that makes the BIVALENT,
        #     sub-saturating M_ss EMERGE from the antagonism balance (rather than be imposed), sharpens the
        #     switch and adds hysteresis. g_mk is the (deliberately weak/tunable) eviction strength: the
        #     measured 5.07x MB/GNP CyclinD1 fold pins the feedback to the weak regime. g_mk=0 recovers the
        #     prior read-write-only module exactly.
        # minus demethylation/turnover (del*Mk; tau_restore ~ ln2/del ~ 10 h, the lit-review fast end).
        # Discrete: the mark is HALVED at EARLY S (Dna>0.05 event) = replicative dilution -- the emergent
        # cycle period sets the dilution frequency, so the dilution-vs-restoration race is the T_cc
        # phenotype (Jadhav 2020). Readout = LEAKY frequency-acting Hill R = f0 + (1-f0)/(1+(Mk/K)^n):
        # H3K27me3 throttles initiation/burst frequency with a residual floor f0 (loaded Ser5P Pol II ->
        # leaky firing; lit-review S3.1/S3.8) -- it DAMPENS, never locks out. NOTE: reformulates the
        # repression term, so the MB/GNP fold is re-derived via EZH2/transcription -> M_ss (separate
        # calibration). Default OFF; not for use together with with_h3k27_memory.
        # Shared H3K27me3 (Mk) dynamics: basal turnover + Gli->Jmjd3/Kdm6b active eraser + replicative dilution.
        shared = (
            "\n  species Mk in Cell; Mk = 0.20;   // H3K27me3 occupancy at Ccnd1 domain [0,1] (init derepressed)"
            "\n  del_mk = 0.0015; k_jmjd3_gli = 0.047832059241740936;   // basal H3K27me3 turnover + Gli->Jmjd3/Kdm6b ACTIVE eraser (Shi 2014 ncomms6425): WRITER=EZH2(cycle) vs ERASER=Gli(mitogen) race (optimize_prc2 2026-07-10)"
            "\n  K_tx_mk = 5.119476334025431; p_tx_mk = 3.2652610235558535;   // nascent-transcription -> PRC2 eviction Hill on Cd_mRNA"
            "\n  Mk_turnover: Mk => ; Cell*del_mk*Mk;"
            "\n  Mk_demeth_jmjd3: Mk => ; Cell*k_jmjd3_gli*Gli1*Mk;   // MB (high Gli) actively strips the mark -> ChIP MB<GNP"
            "\n  Mk_replicative_dilution: at (Dna > 0.05): Mk = 0.5*Mk;"
        )
        if with_prc2 and with_h3k27_chain:
            # SERIAL METHYLATION CHAIN (JP 2026-07-11): resolve the lumped Mk into me0->me1->me2->me3.
            # PRC2 catalyses all three forward steps; me0->me1 and me1->me2 are FAST, me2->me3 (=> Mk) is
            # the SLOW rate-limiting step (biochemistry: the final transition is ~an order slower). This
            # gives me3 an accumulation LAG and makes it a division-timing readout (fast cycling dilutes it
            # out before the slow step finishes) -- the mechanistic origin of the "persistence" behaviour.
            # Mk STILL denotes me3 (the repressive, heritable state) so PRC2 read-write (a_rw*Mk) and the
            # CyclinD1 Hill are unchanged. m0 = 1 - m1_me - m2_me - Mk (implicit unmethylated pool). Eraser
            # (Gli->Jmjd3/Kdm6b) demethylates stepwise me3->me2->me1 so high-Gli MB strips me3 (ChIP MB<GNP).
            # Replication halves ALL modified fractions (new histones = me0). Default OFF (validation-preserving).
            mk = (
                "\n  species m1_me in Cell; m1_me = 0.05;   // H3K27me1 fraction"
                "\n  species m2_me in Cell; m2_me = 0.15;   // H3K27me2 fraction"
                "\n  species Mk in Cell; Mk = 0.15;   // H3K27me3 (= me3): repressive/heritable state at Ccnd1"
                "\n  del_mk = 0.0016134242382228271; k_jmjd3_gli = 0.108;   // passive turnover + Gli->Jmjd3/Kdm6b eraser (optimize_chain 2026-07-11; eraser stronger to keep ChIP MB<GNP at the faster me2->me3)"
                "\n  K_tx_mk = 5.119476334025431; p_tx_mk = 3.2652610235558535;   // nascent-tx -> PRC2 eviction Hill"
                "\n  a0_prc2 = 0.0005319329066695013; a_rw_prc2 = 0.003731831639056823; g_prc2 = 0.030623647714812362;"
                "\n  K_prc2 = 0.0028684066713754613; n_prc2 = 3.982266551808909; f0_prc2 = 0.18030409984207782;"
                "\n  kme1 = 9.706031386488556; kme2 = 8.266881635062989; kme3 = 4.555295001018562;   // serial rates x PRC2: FAST me0->me1,me1->me2; SLOW me2->me3 (rate-limiting; optimize_chain 2026-07-11, me3 lag ~13h)"
                "\n  d_me = 0.0014787645946766597;   // precursor (me1,me2) passive turnover"
                # ---- 2026-07-21 four mechanism additions (JP-approved); all NO-OP at these string defaults ----
                "\n  n_rw = 1;   // read-write cooperativity exponent: a_rw_prc2*Mk^n_rw (1 = linear/legacy; >1 = cooperative EED spreading -> bistable-capable, gap #1)"
                "\n  k_h33 = 0;   // replication-INDEPENDENT H3.3/HIRA transcription-coupled me3 loss (k_h33*txHill*Mk -> me0): mark decay in NON-dividing cells (gap #3). 0 = off"
                "\n  dil_frac = 0.5;   // discrete S-onset dilution fraction (0.5 = legacy instantaneous halving); set 0 to use the continuous k_dil_S term instead"
                "\n  k_dil_S = 0;   // continuous replication dilution x DNA-synth flux (kSyDna*vfork*aRc): spreads dilution over S so an elongated S recovers concurrently (gap #2). 0 = off; ~0.693 => ~50% over a full S"
                # EZH2i is a SAM-competitive CATALYTIC inhibitor: it blocks the SET-domain methyltransferase
                # but leaves the PRC2 complex BOUND and the existing H3K27me3 in place. So (1-EZH2i) acts ONLY on
                # the writing activity (PRC2, below, drives the me-chain); the CyclinD1 repression is driven by
                # PRC2 OCCUPANCY (PRC2_rep, NO EZH2i factor) -> after EZH2i, repression falls only as the mark
                # decays (turnover + eraser + replicative dilution) = mark-decay-LIMITED, GRADUAL de-repression,
                # not the old instant collapse. At EZH2i=0, PRC2 == PRC2_rep, so baseline is unchanged.
                "\n  PRC2 := EZH2*(1 - EZH2i)*(a0_prc2 + a_rw_prc2*Mk^n_rw)*(1 - g_prc2*Cd_mRNA^p_tx_mk/(K_tx_mk^p_tx_mk + Cd_mRNA^p_tx_mk));   // CATALYTIC (writes me): EZH2i-blocked"
                "\n  PRC2_rep := EZH2*(a0_prc2 + a_rw_prc2*Mk^n_rw)*(1 - g_prc2*Cd_mRNA^p_tx_mk/(K_tx_mk^p_tx_mk + Cd_mRNA^p_tx_mk));   // OCCUPANCY (represses CyclinD1): complex stays bound under EZH2i"
                "\n  me0_to_me1: => m1_me; Cell*PRC2*kme1*(1 - m1_me - m2_me - Mk);   // FAST"
                "\n  me1_to_me2: m1_me => m2_me; Cell*PRC2*kme2*m1_me;                // FAST"
                "\n  me2_to_me3: m2_me => Mk; Cell*PRC2*kme3*m2_me;                    // SLOW (rate-limiting)"
                "\n  me3_demeth: Mk => m2_me; Cell*(del_mk + k_jmjd3_gli*Gli1)*Mk;    // eraser strips me3->me2"
                "\n  me2_demeth: m2_me => m1_me; Cell*(d_me + k_jmjd3_gli*Gli1)*m2_me;"
                "\n  me1_turnover: m1_me => ; Cell*d_me*m1_me;                         // -> me0 (implicit)"
                "\n  me3_h33_exchange: Mk => ; Cell*k_h33*Cd_mRNA^p_tx_mk/(K_tx_mk^p_tx_mk + Cd_mRNA^p_tx_mk)*Mk;   // H3.3/HIRA replication-INDEPENDENT, transcription-coupled histone exchange -> me0 (gap #3)"
                "\n  Mk_replicative_dilution: at (Dna > 0.05): m1_me = (1 - dil_frac)*m1_me, m2_me = (1 - dil_frac)*m2_me, Mk = (1 - dil_frac)*Mk;"
                "\n  me1_dil_S: m1_me => ; Cell*k_dil_S*kSyDna*vfork*aRc*m1_me;   // continuous S-phase (fork-gated) dilution -> me0 (gap #2)"
                "\n  me2_dil_S: m2_me => ; Cell*k_dil_S*kSyDna*vfork*aRc*m2_me;"
                "\n  me3_dil_S: Mk => ; Cell*k_dil_S*kSyDna*vfork*aRc*Mk;"
            )
            m = m.replace("\nend", mk + "\nend")
            m = m.replace("(K_EZH2_repression/(K_EZH2_repression + EZH2*(1 - EZH2i)))",
                          "(f0_prc2 + (1 - f0_prc2)/(1 + (PRC2_rep/K_prc2)^n_prc2))")
        elif with_prc2:
            # FORMAL PRC2 COMPLEX (JP 2026-07-10): CyclinD1 repression is carried by PRC2 OCCUPANCY at the locus,
            # NOT the mark level. PRC2 = EZH2 complex abundance x (accessory recruitment a0 + EED read-write a_rw*Mk),
            # minus nascent-transcription eviction. The SAME PRC2 both WRITES H3K27me3 and REPRESSES CyclinD1, so
            # (i) the mark is a read-write AMPLIFIER not an independent repressor; (ii) CyclinD1 stays responsive to
            # EZH2 OE/i (paper) because ALL repression is EZH2/PRC2-mediated; (iii) MB keeps occupancy via high EZH2
            # (5.07 dose fold) while Gli/Jmjd3 keeps the MARK low (ChIP). Replaces the phenomenological w_ezdir blend.
            mk = shared + (
                "\n  a0_prc2 = 0.00037763461927116913; a_rw_prc2 = 0.004297010112265872; g_prc2 = 0.04244445931486211;   // PRC2 recruitment: accessory (mark-indep, sequence/SUZ12) + H3K27me3 read-write (EED); nascent-tx eviction (optimize_prc2 2026-07-10)"
                "\n  n_rw = 1;   // read-write cooperativity exponent (a_rw_prc2*Mk^n_rw); declared here too so the non-chain PRC2 branch builds (else 'n_rw missing'); _ts_bake sets 2"
                "\n  K_prc2 = 0.0034813553293576824; n_prc2 = 3.385026804758643; f0_prc2 = 0.05054636367014487;   // CyclinD1 repression Hill on PRC2 OCCUPANCY + leaky floor (Pol II retained)"
                # EZH2i = catalytic inhibitor: (1-EZH2i) on the WRITING (PRC2) only; repression via PRC2_rep OCCUPANCY (no EZH2i)
                "\n  PRC2 := EZH2*(1 - EZH2i)*(a0_prc2 + a_rw_prc2*Mk^n_rw)*(1 - g_prc2*Cd_mRNA^p_tx_mk/(K_tx_mk^p_tx_mk + Cd_mRNA^p_tx_mk));   // CATALYTIC (writes me): EZH2i-blocked"
                "\n  PRC2_rep := EZH2*(a0_prc2 + a_rw_prc2*Mk^n_rw)*(1 - g_prc2*Cd_mRNA^p_tx_mk/(K_tx_mk^p_tx_mk + Cd_mRNA^p_tx_mk));   // OCCUPANCY (represses CyclinD1): complex stays bound under EZH2i"
                "\n  Mk_methylation: => Mk; Cell*PRC2*(1 - Mk);   // PRC2 writes H3K27me3 on unmethylated substrate"
            )
            m = m.replace("\nend", mk + "\nend")
            m = m.replace("(K_EZH2_repression/(K_EZH2_repression + EZH2*(1 - EZH2i)))",
                          "(f0_prc2 + (1 - f0_prc2)/(1 + (PRC2_rep/K_prc2)^n_prc2))")
        else:
            # LEGACY (baked 9fbc99c): w_ezdir blend of a mark-Hill and EZH2-direct dose repression.
            mk = shared + (
                "\n  k_w_mk = 0.003105363657385844; k0_mk = 0.0007839663325000536;   // read-write, de-novo floor"
                "\n  g_mk = 0.11134937424607448; w_ezdir = 0.9006939037085142;   // eviction + weight of EZH2-DIRECT vs mark-Hill repression"
                "\n  K_mk = 0.17832905051073167; n_mk = 2.8959220031025135; f0_mk = 0.22565611687374498;   // Ccnd1 repression Hill + leaky floor"
                "\n  Mk_methylation: => Mk; Cell*EZH2*(1 - EZH2i)*(k_w_mk*Mk + k0_mk)*(1 - Mk)*(1 - g_mk*Cd_mRNA^p_tx_mk/(K_tx_mk^p_tx_mk + Cd_mRNA^p_tx_mk));"
            )
            m = m.replace("\nend", mk + "\nend")
            m = m.replace("(K_EZH2_repression/(K_EZH2_repression + EZH2*(1 - EZH2i)))",
                          "((1 - w_ezdir)*(f0_mk + (1 - f0_mk)/(1 + (Mk/K_mk)^n_mk)) + w_ezdir*(K_EZH2_repression/(K_EZH2_repression + EZH2*(1 - EZH2i))))")
        # ---- POISED BIVALENT DIFFERENTIATION GENE (JP 2026-07-21; added AFTER the mark-block if/elif/else so it
        # never breaks that chain). Off by default (with_diff_gene=False) -> default model unchanged. A heavily-
        # H3K27me3-marked differentiation locus that is ALSO H3K4me3+ (bivalent / poised): repressed in the
        # progenitor but "ready to fire" the instant H3K27me3 falls. Written by the SAME EZH2/PRC2 writer and diluted
        # by the SAME S-phase replication as CyclinD1, but with NO feedback to the cell cycle (a passive SENTINEL /
        # readout -> validation-neutral). Writing is SLOW (restoration ~ Tc, the me3 de-novo lag) so the mark is
        # DIVISION-RATE-SENSITIVE: at fast Tc it dilutes below the poised H3K4me3 level and the gene PREMATURELY
        # ACTIVATES (Diff -> 1) -- the "too fast -> differentiation" event of the speed-limiter framework. Also fires
        # under EZH2i (writer loss). Requires the chain (Mk) + growth (Dna/aRc/vfork/kSyDna). ----
        if with_diff_gene and with_prc2 and with_h3k27_chain:
            dg = (
                "\n  species me1_diff in Cell; me1_diff = 0.05;   // H3K27me1 at the differentiation locus"
                "\n  species me2_diff in Cell; me2_diff = 0.12;   // H3K27me2"
                "\n  species Mk_diff in Cell;  Mk_diff = 0.85;    // H3K27me3 (= me3): the heavy repressive mark at the diff locus"
                "\n  species Diff in Cell;     Diff = 0.0;         // PREMATURE differentiation-GENE ACTIVATION readout (0 = repressed/poised)"
                "\n  H3K4_diff = 0.55;                             // H3K4me3 (poising ACTIVE mark) -- constitutively present at the bivalent locus"
                "\n  a0_diff = 0.0014; a_rw_diff = 0.006; n_rw_diff = 2;   // PRC2 recruitment at the diff locus (accessory + me3 read-write)"
                "\n  kme1_diff = 6.0; kme2_diff = 5.0; kme3_diff = 0.35;   // SERIAL rates x PRC2_diff: FAST me0->me1,me1->me2; SLOW me2->me3 (=> me3) = the de-novo LAG -> me3 is division-rate-sensitive"
                "\n  del_diff = 0.0005; d_me_diff = 0.0007;   // me3->me2 turnover + precursor (me1,me2) turnover"
                "\n  k_dil_diff = 0.693;   // continuous S-phase (fork-gated) replicative dilution (halves all me over a full S)"
                "\n  K_diff = 0.28; n_diff = 8; k_on_diff = 0.05; k_off_diff = 0.02;   // HIGH activation threshold: the gene fires ONLY when me3 is LARGELY LOST (Mk_diff < ~0.28) -> robustly repressed"
                "\n  k_erase_diff = 0.12; K_erase = 0.7; n_erase = 10;   // bivalent RESOLUTION: an active gene erases its own me3 (KDM6B) -> latches the activation (silent below Diff~0.7)"
                "\n  g_diff_exit = 0.97; K_diff_exit = 0.5; n_diff_exit = 4;   // DOWNSTREAM (G0/differentiation is complex): once the gene is ON the cell COMPLETES its final mitosis then exits to G0 (represses CyclinD1 re-commitment; a few residual divisions during the transition)"
                "\n  PRC2_diff := EZH2*(1 - EZH2i)*(a0_diff + a_rw_diff*Mk_diff^n_rw_diff);   // PRC2 writing at the diff locus (EZH2-driven, me3 read-write)"
                "\n  diff_drive := H3K4_diff*K_diff^n_diff/(K_diff^n_diff + Mk_diff^n_diff);   // high only when me3 is LOW (repression relieved; H3K4me3 poises the gene)"
                "\n  diff_exit := 1 - g_diff_exit*Diff^n_diff_exit/(K_diff_exit^n_diff_exit + Diff^n_diff_exit);   // 1 when Diff low; -> (1-g) when the gene is active (bias to G0)"
                "\n  me0d_to_me1d: => me1_diff; Cell*PRC2_diff*kme1_diff*(1 - me1_diff - me2_diff - Mk_diff);   // FAST"
                "\n  me1d_to_me2d: me1_diff => me2_diff; Cell*PRC2_diff*kme2_diff*me1_diff;                    // FAST"
                "\n  me2d_to_me3d: me2_diff => Mk_diff; Cell*PRC2_diff*kme3_diff*me2_diff;                     // SLOW (rate-limiting) = de-novo lag"
                "\n  me3d_demeth: Mk_diff => me2_diff; Cell*del_diff*Mk_diff;                                  // passive turnover me3->me2"
                "\n  me2d_demeth: me2_diff => me1_diff; Cell*d_me_diff*me2_diff;"
                "\n  me1d_turn:   me1_diff => ; Cell*d_me_diff*me1_diff;                                       // -> me0 (implicit)"
                "\n  me3d_erase: Mk_diff => me2_diff; Cell*k_erase_diff*Diff^n_erase/(K_erase^n_erase + Diff^n_erase)*Mk_diff;   // active-gene self-erasure (latch)"
                "\n  meX_dil_S: me1_diff => ; Cell*k_dil_diff*kSyDna*vfork*aRc*me1_diff;                       // S-phase dilution (fork-gated)"
                "\n  meY_dil_S: me2_diff => ; Cell*k_dil_diff*kSyDna*vfork*aRc*me2_diff;"
                "\n  meZ_dil_S: Mk_diff => ; Cell*k_dil_diff*kSyDna*vfork*aRc*Mk_diff;"
                "\n  Diff_activate: => Diff; Cell*k_on_diff*diff_drive*(1 - Diff);                             // PREMATURE ACTIVATION when me3 collapses"
                "\n  Diff_repress:  Diff => ; Cell*k_off_diff*Mk_diff*Diff;                                    // re-repressed if me3 recovers (before resolution)"
            )
            m = m.replace("\nend", dg + "\nend")
            # DOWNSTREAM coupling: a fired diff gene biases the cell to G0 (completes its final mitosis, then stops
            # re-committing) by multiplying the CyclinD1 transcription by diff_exit. NOT an abrupt mid-cycle arrest --
            # differentiation itself (G0-dependent) is complex and not modelled in detail. Diff low => diff_exit=1.
            m = m.replace("(f0_prc2 + (1 - f0_prc2)/(1 + (PRC2_rep/K_prc2)^n_prc2))",
                          "(f0_prc2 + (1 - f0_prc2)/(1 + (PRC2_rep/K_prc2)^n_prc2))*diff_exit")

    if with_two_step_rb and with_ezh2:
        # Baked two-step-Rb calibration (optimize_twostep_g0 2026-07-14 JOINT G0+EZH2 re-opt, 28/28 -- the
        # critical-review root-cause fix). This round genuinely LANDS the CyclinD1 MB/GNP fold at 5.6 (was
        # 7.1, passing only on the wide band) and RECOVERS a real MB p27-high G0 dwell (~12%) at birth-p27
        # 1.39 (was the inflated 1.8 crutch), while FIXING both standing failures (GNP+HHi de-repression;
        # MB HU G2). Mechanism: (i) p27 clearance re-gated by *effective* CDK4/6 activity (w_p27*P21 self-
        # brake -> mutual antagonism, not a saturated switch); (ii) EZH2 synthesis weighted toward the
        # MITOGEN-DOSE term (Kez_cd 3.2->9.2) so EZH2 stays high through the G0 dwell -> decouples the
        # fold-vs-G0 tension. TRADE (data-forced): EZH2i de-repression is 1.80 (was 2.47) -- strong EZH2i
        # (>=2.2) is EXCLUDED by the fold+GNP+HHi+ChIP (proven by sweep), consistent with the weak-feedback
        # regime the 5.07 fold already pins. Applied BRANCH-SPECIFICALLY (with_two_step_rb=False keeps the
        # single-step bake). User `params` still override. pRb(hyper) = G0 marker (Moser 2018); R-point
        # logic (Narasimha/Sanidas); p27 (P21 pool) the dominant CIP/KIP in GNP/MB.
        if with_cdki_species:
            # EXPLORATION (JP 2026-07-30): model the individually-expressed CDKIs as dynamic species.
            # STAGE 1 (INK4): p18 (Cdkn2c) + p19 (Cdkn2d) as protein species with turnover; drop p16 (silent in GNP)
            # from the CDK4/6 brake. Input params p18/p19 = expression-set target LEVELS; the *_prot species relax
            # to them (kDe turnover) -> steady state == target, but with realistic protein kinetics and a hook for
            # transcriptional (Gli/EZH2) synthesis regulation later. CIP/KIP (p21,p57) = Stage 2.
            m = m.replace(
                "\n  p18 = 0.464;",
                "\n  p18 = 0.464;\n  p19 = 0.36;   // Cdkn2c(p18)/Cdkn2d(p19) INK4 target levels (RNA-seq: p18 dominant, p19 secondary)"
                "\n  kDe_p18 = 0.003; kDe_p19 = 0.03;   // INK4 turnover: p18 STABLE (tau~4h); p19^INK4d SHORT-lived (tau~23min, S-peaking) (reviewer)"
                "\n  species p18_prot in Cell; p18_prot = 0.464;\n  species p19_prot in Cell; p19_prot = 0.36;")
            m = m.replace(
                "\nend",
                "\n  p18_synthesis: => p18_prot; Cell*kDe_p18*p18;"
                "\n  p18_degradation: p18_prot => ; Cell*kDe_p18*p18_prot;"
                "\n  p19_synthesis: => p19_prot; Cell*kDe_p19*p19;"
                "\n  p19_degradation: p19_prot => ; Cell*kDe_p19*p19_prot;\nend", 1)
            m = m.replace("1 + p16 + p18 + w_p27*P21", "1 + p18_prot + p19_prot + w_p27*P21")
            # PROVENANCE FIX (reviewer 2026-07-30): CRL4^Cdt2 (kDeP21aRc*Cdt2*aRc) is a p21 route (PIP-degron), NOT
            # p27 (no PIP box). The Heldt species is p21, relabeled p27 -> the term rode along. Strip it off p27 here
            # (p27 then persists further into S) and re-attach it to p21a below (aRc-gated S-phase clearance = the
            # branch-divergence term: origins fire -> p21 slammed to 0; arrest -> term is 0, so it does NOT feed the runaway).
            m = m.replace(" + kDeP21aRc*Cdt2*aRc", "")
            # STAGE 2 (CIP/KIP): p21 (Cdkn1a) as a dynamic species sequestering CyclinE/A-CDK2 (own complexes,
            # mirroring p27/kAsCyP21 kinetics) + adding to the CDK4/6 brake via w_p27*(P21 + p21a). Synthesis is
            # basal + p53-inducible (kSyp21aP53; P53~0 default). GNP synthesis kept LOW (p21 ~10% of p27 in GNP) so
            # the GNP<->MB dichotomy holds; MB elevates it (data 2.8x). PCNA/Rc arm skipped (p27 dominates it).
            m = m.replace(
                "\nend",
                "\n  species p21a in Cell; p21a = 0.05;   // Cdkn1a (p21) functional pool"
                "\n  species Cep21a in Cell; Cep21a = 0;\n  species Cap21a in Cell; Cap21a = 0;"
                "\n  kSyp21a = 0.00002; kSyp21aP53 = 0.0001; kDep21a = 0.02; w_cdk2_p21 = 0.3;   // p21 synth+turnover; w_cdk2_p21 = CDK2 potency vs p27 (<1: p21 weaker; reviewer)"
                "\n  // NB start p21a NEAR-ZERO so flag-on preserves the GNP<->MB dichotomy; Stage 3 raises kSyp21a to the"
                "\n  // data level (p21 ~10% of p27 GNP, 2.8x in MB) WITH compensating re-tune (avoid the Skp2-gated runaway)."
                "\n  Synthesis_of_p21a: => p21a; Cell*(kSyp21a + kSyp21aP53*P53);"
                "\n  Assoc_CycE_Cdk2_p21a: Ce + p21a -> Cep21a; Cell*(kAsCyP21*w_cdk2_p21*Ce*p21a - kDsCyP21*Cep21a);"
                "\n  Assoc_CycA_Cdk2_p21a: Ca + p21a -> Cap21a; Cell*(kAsCyP21*w_cdk2_p21*Ca*p21a - kDsCyP21*Cap21a);"
                "\n  Deg_free_p21a: p21a => ; Cell*(kDep21a + kDeP21Cy*Skp2*(Ce + Ca) + kDeP21aRc*Cdt2*aRc)*p21a;"
                "\n  Deg_p21a_in_CycE: Cep21a => Ce; Cell*(kDep21a + kDeP21Cy*Skp2*(Ce + Ca) + kDeP21aRc*Cdt2*aRc)*Cep21a;"
                "\n  Deg_p21a_in_CycA: Cap21a => Ca; Cell*(kDep21a + kDeP21Cy*Skp2*(Ce + Ca) + kDeP21aRc*Cdt2*aRc)*Cap21a;"
                "\n  Deg_CycE_in_Cep21a: Cep21a => p21a; Cell*((kDeCe + kDeCeCa*Ca)*Cep21a);"
                "\n  Deg_CycA_in_Cap21a: Cap21a => p21a; Cell*((kDeCa + kDeCaC1*C1)*Cap21a);\nend", 1)
            m = m.replace("w_p27*P21", "w_p27*(P21 + p21a)")
            m = m.replace("tCe := Ce + CeP21;", "tCe := Ce + CeP21 + Cep21a;")
            m = m.replace("tCa := Ca + CaP21;", "tCa := Ca + CaP21 + Cap21a;")
            m = m.replace("CeP21 = 0, CaP21 = 0,", "CeP21 = 0, CaP21 = 0, Cep21a = 0, Cap21a = 0,")
            # p57 (Cdkn1c): DE-EMPHASIZED (JP 2026-07-30) -- minimally expressed in proliferating GNP/MB; largely a
            # Sox2+ quiescent-cell CKI, not a functional player here. Species/reactions kept (structure) but carried
            # OFF (kSyp57a = 0 -> p57a decays to ~0, vestigial). p53-independent; no CDK4/6 brake role.
            m = m.replace(
                "\nend",
                "\n  species p57a in Cell; p57a = 0.0;"
                "\n  species Cep57a in Cell; Cep57a = 0;\n  species Cap57a in Cell; Cap57a = 0;"
                "\n  kSyp57a = 0.0; kDep57a = 0.02; w_cdk2_p57 = 0.3;   // p57 OFF (Sox2+/quiescent CKI; de-emphasized in proliferating cells)"
                "\n  Synthesis_of_p57a: => p57a; Cell*kSyp57a;"
                "\n  Assoc_CycE_Cdk2_p57a: Ce + p57a -> Cep57a; Cell*(kAsCyP21*w_cdk2_p57*Ce*p57a - kDsCyP21*Cep57a);"
                "\n  Assoc_CycA_Cdk2_p57a: Ca + p57a -> Cap57a; Cell*(kAsCyP21*w_cdk2_p57*Ca*p57a - kDsCyP21*Cap57a);"
                "\n  Deg_free_p57a: p57a => ; Cell*(kDep57a + kDeP21Cy*Skp2*(Ce + Ca))*p57a;"
                "\n  Deg_p57a_in_CycE: Cep57a => Ce; Cell*(kDep57a + kDeP21Cy*Skp2*(Ce + Ca))*Cep57a;"
                "\n  Deg_p57a_in_CycA: Cap57a => Ca; Cell*(kDep57a + kDeP21Cy*Skp2*(Ce + Ca))*Cap57a;"
                "\n  Deg_CycE_in_Cep57a: Cep57a => p57a; Cell*((kDeCe + kDeCeCa*Ca)*Cep57a);"
                "\n  Deg_CycA_in_Cap57a: Cap57a => p57a; Cell*((kDeCa + kDeCaC1*C1)*Cap57a);\nend", 1)
            m = m.replace("w_p27*(P21 + p21a)", "w_p27*(P21 + p21a + p57a)")
            m = m.replace("tCe := Ce + CeP21 + Cep21a;", "tCe := Ce + CeP21 + Cep21a + Cep57a;")
            m = m.replace("tCa := Ca + CaP21 + Cap21a;", "tCa := Ca + CaP21 + Cap21a + Cap57a;")
            m = m.replace("Cep21a = 0, Cap21a = 0,", "Cep21a = 0, Cap21a = 0, Cep57a = 0, Cap57a = 0,")
            # STAGE 2.5 (reviewers 2026-07-30, round 2): commitment swap + degradation corrections.
            # KPC: p27 clearance is Skp2-indep, mitogen-insensitive, and acts on CYTOPLASMIC Ser10-P p27 -> a CONSTANT
            #   kDeKPC on the FREE pool ONLY (NOT the nuclear cyclin/PCNA complexes, NOT CdP21). Drop the old CDK4/6-gated
            #   clearance from every pool; basal+KPC put p27 t1/2 ~1-2h (kDeKPC ~0.006, not 0.05 = 13min bug).
            m = m.replace(" + kDeP21Cd*kPhRbCd*(Cd + w_Cd2*Cd2)/(K_CdRb*(1 + p18_prot + p19_prot + w_p27*(P21 + p21a + p57a)) + (Cd + w_Cd2*Cd2))", "")
            m = m.replace("Degradation_of_free_p21: P21 => ; Cell*((kDeP21 + ",
                          "Degradation_of_free_p21: P21 => ; Cell*((kDeP21 + kDeKPC + ")
            # CdP21 redistribution (option A, STOICHIOMETRIC + saturating): p27 binds CyclinD1-CDK4/6 (Cd CONSUMED -> the
            #   buffer SATURATES at Cd abundance -> non-linear collapse as Cd falls); the buffered Cd is returned to the
            #   Rb-drive NUMERATOR (Cd+Cd2+CdP21) so total drive = Cd_total (option A: the p27-CDK4 trimer is part of the
            #   drive, not a bystander). Free p27 is REMOVED from the Rb DENOMINATOR (it can't brake CDK4/6 until bound ->
            #   acts on CDK2 only; fixes the free-inhibits/bound-doesn't double-hit). INK4 promotes release. CdP21 is NOT
            #   reset at division (CyclinD1 persists through mitosis -> no mass sink; buffer stays full into G1). Stage 3 recal.
            m = m.replace(
                "(Cd + w_Cd2*Cd2)/(K_CdRb*(1 + p18_prot + p19_prot + w_p27*(P21 + p21a + p57a)) + (Cd + w_Cd2*Cd2))",
                "(Cd + w_Cd2*Cd2 + CdP21)/(K_CdRb*(1 + p18_prot + p19_prot) + (Cd + w_Cd2*Cd2 + CdP21))")   # p21/p57 dropped from the CDK4/6 brake (poor CDK4/6 modulators; they act on CDK2 only). INK4 = the CDK4/6 brake; p27 via CdP21 stoichiometry.
            m = m.replace(
                "\nend",
                "\n  kDeKPC = 0.006;   // p27 KPC clearance: CONSTANT, FREE-pool only; basal+KPC t1/2 ~1-2h"
                "\n  kSeqCd = 0.05; kRelCd = 0.01; w_ink4 = 1.0;   // CyclinD1-p27 buffer (bimolecular, saturating); INK4 releases (PLACEHOLDER)"
                "\n  species CdP21 in Cell; CdP21 = 0;   // p27 buffered on CyclinD1-CDK4/6 (option A: counts in the Rb drive)"
                "\n  Buffer_p27_on_CyclinD: Cd + P21 -> CdP21; Cell*(kSeqCd*Cd*P21 - kRelCd*(1 + w_ink4*(p18_prot + p19_prot))*CdP21);"
                "\n  Deg_p27_in_CdP21: CdP21 => Cd; Cell*kDeP21*CdP21;   // buffered p27 turns over at BASAL only; Cd released\nend", 1)
            m = m.replace("tP21 := P21 + CeP21 + CaP21 + iPcna + iRc;", "tP21 := P21 + CeP21 + CaP21 + iPcna + iRc + CdP21;")
            # PCNA/Rc arm — mechanistic MOVE REVERTED (2026-07-30, round 3). Moving iPcna/iRc p27->p21a is
            # NUMERICALLY PATHOLOGICAL: iPcna/iRc are replication-coupled shared species that hold mass while the
            # near-zero p21a pool is dragged strongly NEGATIVE (p21a min ~-0.4, corrupting the sim to 15/32) — the
            # "mechanistic version for free" is not free in this shared-species framework. So p21 keeps its Cdt2/
            # S-phase clearance via the aRc PROXY (in p21a's free/complex degradation, from the Stage-2 swap), which
            # is robust. KNOWN LIMITATION for review: the base-Heldt PCNA/Rc binding stays nominally on p27 (p27's
            # Cdt2 degradation IS removed); a fully mechanistic p21-PCNA arm needs dedicated (non-shared) p21-PCNA
            # species, deferred. tP21 keeps iPcna/iRc (Heldt convention).
            m = m.replace("tP21 := P21 + CeP21 + CaP21 + iPcna + iRc + CdP21;",
                          "tP21 := P21 + CeP21 + CaP21 + iPcna + iRc + CdP21;")  # unchanged (documented above)

            if with_p21_pip_degron:
                # ===== Mechanistic p21-PCNA / CRL4^Cdt2 arm (UN-MAP; JP path A, 2026-08-01). DEFAULT-OFF; nothing baked. =====
                # The PIP-degron machinery is INHERITED, not missing: Heldt's "P21" (which v44 remapped to p27) IS p21 --
                # iPcna = p21.PCNA, iRc = p21.loaded-PCNA (RCi). p27 has no PIP box, so it must not bind PCNA. This block
                # RESTORES the machinery to p21: repoint the PCNA/Rc binding p27(P21)->p21a, put the CRL4^Cdt2 degron on
                # the p21a-bound iPcna/iRc (loaded PCNA), drop the free-pool aRc PROXY, and fix the moieties. It KEEPS the
                # p21->fork (aRc->iRc) edge needed for the Barr-2017 Cdt2-depletion phenotype + the bistable G1/S switch.
                # REQUIRES p21 at a meaningful abundance (raise kSyp21a at runtime): the round-3 negativity was the
                # near-zero pool x kAsPcP21=100 stiffness breaking the p21 moiety, NOT the un-map itself. Quantitative
                # knobs (kSyp21a, w_cdk2_p21, Cdt2) are LEFT at their current values -- set per-experiment / calibrate
                # downstream. NB p21 abundance is pinned from TRANSCRIPT fold (Cdkn1a ~2.8x MB/GNP); the emergent p21
                # PROTEIN is a PREDICTION (no p21 protein quants). See docs/p21_pcna_arm_design.md.
                # (a) repoint PCNA/Rc binding p27(P21) -> p21a
                m = m.replace("aPcna + P21 -> iPcna; Cell*(kAsPcP21*aPcna*P21 - kDsPcP21*iPcna)",
                              "aPcna + p21a -> iPcna; Cell*(kAsPcP21*aPcna*p21a - kDsPcP21*iPcna)")
                m = m.replace("aRc + P21 -> iRc; Cell*(kAsPcP21*aRc*P21 - kDsPcP21*iRc)",
                              "aRc + p21a -> iRc; Cell*(kAsPcP21*aRc*p21a - kDsPcP21*iRc)")
                m = m.replace("Nuclear_export_of_inactive_PCNA: iPcna => P21; Cell*kExPc*iPcna;",
                              "Nuclear_export_of_inactive_PCNA: iPcna => p21a; Cell*kExPc*iPcna;")
                # (b) CRL4^Cdt2 degron restored ON the p21a-bound (loaded-PCNA) complexes = mechanistic S-phase clearance
                m = m.replace("iPcna => aPcna; Cell*((kDeP21 + kDeP21Cy*Skp2*(Ce + Ca))*iPcna);",
                              "iPcna => aPcna; Cell*((kDeP21 + kDeP21Cy*Skp2*(Ce + Ca) + kDeP21aRc*Cdt2*aRc)*iPcna);")
                m = m.replace("iRc => aRc; Cell*((kDeP21 + kDeP21Cy*Skp2*(Ce + Ca))*iRc);",
                              "iRc => aRc; Cell*((kDeP21 + kDeP21Cy*Skp2*(Ce + Ca) + kDeP21aRc*Cdt2*aRc)*iRc);")
                # (c) drop the free/CDK2 aRc PROXY on p21a (degradation is now mechanistic, only on PCNA-bound p21)
                m = m.replace(" + kDeP21aRc*Cdt2*aRc)*p21a;", ")*p21a;")
                m = m.replace(" + kDeP21aRc*Cdt2*aRc)*Cep21a;", ")*Cep21a;")
                m = m.replace(" + kDeP21aRc*Cdt2*aRc)*Cap21a;", ")*Cap21a;")
                # (d) moieties: p27(P21) loses iPcna/iRc; p21(p21a) total gains them
                m = m.replace("tP21 := P21 + CeP21 + CaP21 + iPcna + iRc + CdP21;",
                              "tP21 := P21 + CeP21 + CaP21 + CdP21;\n  tp21a := p21a + Cep21a + Cap21a + iPcna + iRc;")
                # guards: fail LOUD if any splice silently no-op'd (a mismatch would strip clearance -> false neutrality)
                assert "aPcna + p21a -> iPcna" in m and "aRc + p21a -> iRc" in m, "p21-PCNA repoint no-op"
                assert m.count("kDeP21aRc*Cdt2*aRc") == 2, "Cdt2 degron not exactly on iPcna+iRc"
                assert "tp21a :=" in m and "tP21 := P21 + CeP21 + CaP21 + CdP21;" in m, "moiety readouts not updated"
                assert "kDeP21aRc*Cdt2*aRc)*p21a;" not in m, "free-p21a proxy not removed"

            if with_p27_optionB:
                # ===== Option B (Fan-Meyer 2021, PMID 34320337): p27 is an INHIBITOR at CyclinD-CDK4/6. EXPLORATION. =====
                # Baked option A treats the buffered p27.CyclinD-CDK4/6 trimer (CdP21) as ACTIVE -- it sits in the Rb-drive
                # NUMERATOR, so p27 sequestration is drive-neutral and p27 is nearly inert as a brake (free P21 ~1% of the
                # pool; verified phaseA A1). Option B makes the buffered complex INACTIVE: drop CdP21 from the active drive.
                # The buffer reaction (Cd + P21 -> CdP21) still CONSUMES Cd, so raising p27 now removes CyclinD1 from the
                # active pool -> p27 INHIBITS CDK4/6 via the CyclinD1/p27 ratio (Fan-Meyer), on top of free p27 -> CDK2.
                # This is the mechanism JP endorsed (p27 drives G0; ratio tunes G1). BREAKS the option-A 30/32 calibration
                # -> needs re-cal (kSeqCd/kRelCd/kDeKPC/w_ink4/kSyP21 + K_CdRb) before it can be a default.
                _optA = "(Cd + w_Cd2*Cd2 + CdP21)/(K_CdRb*(1 + p18_prot + p19_prot) + (Cd + w_Cd2*Cd2 + CdP21))"
                assert m.count(_optA) >= 1, "option-B: option-A drive term not found (structure changed?)"
                m = m.replace(_optA,
                              "(Cd + w_Cd2*Cd2)/(K_CdRb*(1 + p18_prot + p19_prot) + (Cd + w_Cd2*Cd2))")
                assert "+ CdP21)/(K_CdRb" not in m, "option-B: CdP21 still in the Rb drive after replace"

        if with_cdk6_gli:
            # ===== CDK6 = the 2nd arm of the CyclinD-CDK4/6 drive governor (JP 2026-08-07). EXPLORATION, default OFF. =====
            # CDK6 is Hh/GLI2-driven (Hedgehog signaling drives MB growth via CDK6, JCI 2017 PMID 29202464: GLI2 binds
            # the Cdk6 promoter) AND EZH2-marked (JP ChIP) -- the biggest MB drive fold in the bulk data (Cdk6 MB/GNP
            # ~17x, vismo ~0.4x). Modeled as an EZH2-repressible, Gli-driven multiplier 'cdk6' on the CyclinD-arm Rb-
            # phosphorylation Vmax (kPhRbCd), mirroring the CyclinD1 IFFL (Gli up + EZH2 represses; a broad mark tunes
            # kinetics/threshold, so CDK6 is marked AND highly expressed). Makes EZH2 a TWO-PRONGED governor
            # (CyclinD1 + CDK6) of the whole mitogen-drive axis: under vismo both arms fall (2-arm arrest); under EZH2i
            # both de-repress (stronger vismo-resistance). GNP normalized ~1; re-balance kPhRbCd/K_CdRb/w_ink4 to hold
            # ~30/32 (params via _ts_bake / user params). Applies to BOTH Rb-phos drive reactions (CyclinD arm only;
            # kPhRbCe/kPhRbCa = CyclinE/A-CDK2, untouched).
            # REVISED (JP 2026-08-07): CDK6 protein MB/GNP ~5.2x (TMT log2 2.37; RNA 17x); EZH2-cKO GNP de-represses Cdk6
            # +1.24 log2 / 2.39x (strongest of the drive genes -> CDK6 is a direct EZH2 target); vismo ~0.41x. CDK6 >>
            # CyclinD1, so a large UNBOUND CDK6 pool acts as an INK4 (p18) SINK: free CDK6 sops up p18 monomers, RELIEVING
            # the p18 brake on CyclinD-CDK4/6 -- this is how MB's very high CDK6 OVERWHELMS the 3x-elevated p18 CDKI
            # (substrate titration). effective p18 in the brake = w_p18*p18_prot*K_cdk6_sink/(K_cdk6_sink + cdk6). Vismo/
            # EZH2-cKO change cdk6 -> change the sink -> change the effective brake (vismo lowers cdk6 -> p18 brake RESTORED
            # -> arrest; a 2nd reason vismo needs the axis). p18 kept a real brake (w_p18): GNP (low cdk6) stays p18-braked
            # (Uziel: p18-KO GNPs cycle w/o Shh), residual p18 in MB -> the Atoh1+ transient-G0 subpopulation. Active
            # CyclinD-CDK6 drive stays CyclinD-limited (kinase not limiting) -> NO Vmax multiplier; CDK6 acts via the sink.
            # Both p18 AND p19 are INK4s (bind CDK4/6 monomers) -> both are titrated by the unbound CDK6 sink, and both
            # (higher in MB: p18 3.7x, p19 1.6x) contribute to the CDK4/6 brake / G0 threshold (JP 2026-08-07).
            assert m.count("1 + p18_prot + p19_prot") >= 1, "cdk6: INK4 brake term not found (needs default cdki-on)"
            m = m.replace("1 + p18_prot + p19_prot",
                          "1 + (w_p18*p18_prot + w_p19*p19_prot)*K_cdk6_sink/(K_cdk6_sink + cdk6)")
            # CDK6 is a DYNAMIC species with the SAME mark-memory IFFL as CyclinD1 ("responsive like cyclin d1", JP):
            # Gli-driven synthesis x PRC2/H3K27me3-mark repression, REUSING the shared PRC2_rep occupancy (the mark Mk is
            # written by EZH2, ERASED by Gli via k_jmjd3_gli*Gli1). In MB (high Gli) the shared mark is erased -> CDK6
            # de-repressed, and it PERSISTS (mark memory + a SLOW cdk6 reservoir, k_cdk6_deg) -> CDK6 stays ~7x GNP even
            # under vismo (the persistent baseline the acute-factor version couldn't give). CDK6-specific K_prc2_cdk6
            # (more sensitive than CyclinD1) -> bigger MB fold (~17x RNA / ~5x protein) + EZH2-cKO de-repression (~2.39x).
            # Requires the H3K27 chain (PRC2_rep). Free cdk6 still sinks p18 (INK4 titration, wired above).
            assert "PRC2_rep :=" in m, "cdk6 mark-memory needs the H3K27 chain (PRC2_rep); build with with_h3k27_chain=True"
            m = m.replace("\nend",
                "\n  k_cdk6_bas = 0.0002; k_cdk6_Gli = 0.0188; K_Gli_cdk6 = 0.35; n_Gli_cdk6 = 4;   // CDK6 Gli-driven synthesis (mark-memory IFFL like CyclinD1); Gli Hill 0.21->0.40 gives ~5x -> MB fold. Calibrated (analytic; PRC2_rep GNP 0.011/MB 0.020)."
                "\n  f0_cdk6 = 0.42; K_prc2_cdk6 = 0.0029; n_prc2_cdk6 = 4;   // CDK6 PRC2/H3K27me3-mark repression (reuses shared PRC2_rep). f0=0.42 floors GNP&MB repression -> EZH2-cKO (PRC2_rep->0) de-represses 1/0.42 = 2.4x (data 2.39x)."
                "\n  k_cdk6_deg = 0.001;   // SLOW -> cdk6 reservoir -> persistent vismo-resistant baseline (like Cd_mRNA)"
                "\n  K_cdk6_sink = 2.0; w_p18 = 1.0; w_p19 = 1.0;   // free CDK6 sinks BOTH INK4s p18+p19: eff INK4 = (w_p18*p18+w_p19*p19)*K/(K+cdk6). w_p18/w_p19 = INK4 brake strengths (JP: p18,p19 + p27 drive G0)."
                "\n  species cdk6 in Cell; cdk6 = 1.0;   // CDK6 level (mark-memory reservoir); persists under vismo"
                "\n  Cdk6_synthesis: => cdk6; Cell*(k_cdk6_bas + k_cdk6_Gli*(Gli_act + Gli1)^n_Gli_cdk6/(K_Gli_cdk6^n_Gli_cdk6 + (Gli_act + Gli1)^n_Gli_cdk6))*(f0_cdk6 + (1 - f0_cdk6)/(1 + (PRC2_rep/K_prc2_cdk6)^n_prc2_cdk6));"
                "\n  Cdk6_degradation: cdk6 => ; Cell*k_cdk6_deg*cdk6;\nend", 1)
            assert "K_cdk6_sink/(K_cdk6_sink + cdk6)" in m and "Cdk6_synthesis:" in m, "cdk6 mark-memory wiring failed"

        if with_mark_amplifier:
            # ===== JP 2026-08-09 H3K27me3 AMPLIFIER (EXPLORATION, default OFF). Requires with_cdk6_gli + the chain. =====
            # JP hypothesis: a slow / out-of-cycle cell DILUTES H3K27me3 LESS -> the mark BANKS -> represses CDK6 + CyclinD1
            # -> drive erodes over a couple of divisions -> the transient-G0 dwell deepens. The GOVERNOR IFFL
            # (f0_prc2.../f0_cdk6...) CANNOT carry this: it must stay SATURATED to give the EZH2i de-repression fold (2.2x)
            # + the MB/GNP CyclinD1 fold (5x), and a single Hill cannot be both saturated (big STATIC fold) AND high-slope
            # (big DYNAMIC response) at one operating point -- de-saturating it drops validation to 27-29/32 (EZH2i fold
            # -> 1.03, MB/GNP -> 3.3; see docs/mark_amplifier_2026-08-09.md). So the amplifier is a SEPARATE,
            # EXCURSION-GATED repression term mark_amp, multiplied onto BOTH CDK6 synthesis and CyclinD1 transcription, with
            # midpoint K_mamp set ABOVE the normal cycling PRC2_rep (~0.03) so mark_amp ~= 1 during normal cycling (the
            # governor + ALL 31/32 validation conditions are ~unchanged by construction) and only drops (adds repression)
            # once a slow cell banks the mark above the cycling range -> lowers the CDKI arrest threshold + compounding
            # deepening. Gated on PRC2_rep (occupancy) so EZH2-cKO (PRC2_rep->0) switches the amplifier OFF too. f_mamp =
            # floor (max extra repression), K_mamp = mark-excursion midpoint, n_mamp = steepness. Tune via params.
            assert "PRC2_rep :=" in m, "mark amplifier needs the H3K27 chain (PRC2_rep)"
            assert "Cdk6_synthesis: => cdk6; Cell*(" in m, "mark amplifier needs with_cdk6_gli (CDK6 synthesis)"
            assert m.count("CycD1_transcription: => Cd_mRNA; (") == 1, "mark amplifier: CyclinD1 transcription reaction not in expected form"
            m = m.replace("\nend",
                "\n  f_mamp = 0.35; K_mamp = 0.045; n_mamp = 10;   // H3K27me3 amplifier: EXTRA repression of CDK6+CyclinD1 that engages ONLY when the mark banks above the cycling range (slow/out-of-cycle cells). K_mamp > cycling PRC2_rep -> ~no-op at baseline (31/32 preserved); slow cells -> mark_amp -> f_mamp."
                "\n  mark_amp := f_mamp + (1 - f_mamp)/(1 + (PRC2_rep/K_mamp)^n_mamp);   // excursion-gated amplifier (JP 2026-08-09)\nend", 1)
            m = m.replace("Cdk6_synthesis: => cdk6; Cell*(", "Cdk6_synthesis: => cdk6; Cell*mark_amp*(")
            m = m.replace("CycD1_transcription: => Cd_mRNA; (", "CycD1_transcription: => Cd_mRNA; mark_amp*(")
            assert "mark_amp :=" in m and "Cell*mark_amp*(" in m and "Cd_mRNA; mark_amp*(" in m, "mark amplifier wiring failed"

        if with_proximal_distal:
            # ===== JP 2026-08-09 TWO-COMPARTMENT PRC2 / H3K27me3 (EXPLORATION, default OFF). Requires with_cdk6_gli + chain. =====
            # JP's molecular mechanism: CDK6 & CyclinD1 are heavily H3K27me3-marked yet highly EXPRESSED because the
            # PROXIMAL promoter/TSS is kept DEPLETED of mark+complex (Pol2 elongation clears it), while the broad DISTAL
            # region holds H3K27me3 as a RESERVOIR. Repression is the PRC2 COMPLEX (not the mark); the mark recruits/holds
            # the complex via EED read-write, and makes the locus BIDIRECTIONALLY responsive -- when Pol2 stops elongating
            # the reservoir rapidly re-loads the proximal promoter -> re-suppression. Two compartments:
            #   * DISTAL RESERVOIR = Mk (the existing me-chain: read-write substrate, cumulative replicative dilution,
            #     slow turnover) = the MEMORY. Unchanged.
            #   * PROXIMAL OCCUPANCY = P_prox (NEW [0,1] species) = the actual repressor of CDK6+CyclinD1. LOADED by the
            #     reservoir (a_P*Mk^n_anch) but ONLY when Pol2 is NOT elongating (1 - elong_gate); EVICTED by active
            #     elongation (elong_gate) + a small basal. elong_gate = E2f Hill (the proliferation/elongation state).
            # Behaviour: CYCLING (E2f high -> elong_gate~1) -> proximal cleared -> P_prox~0 -> NO extra repression -> the
            # governor + ALL 31/32 validation conditions are UNCHANGED by construction. SUSTAINED ARREST (E2f collapses ->
            # elong_gate->0) -> the PERSISTENT reservoir re-loads P_prox -> deep, reservoir-SUSTAINED repression that does
            # NOT collapse with EZH2 (fixing the PRC2_rep=EZH2xmark collapse) -> a transient-G0 dwell that deepens over
            # time and is REVERSIBLE (mitogen/cycling return -> elongation evicts P_prox -> re-entry). Slow P kinetics +
            # sharp elong_gate keep normal G1 E2f-dips from loading P_prox (only SUSTAINED arrest does). a_P/n_anch =
            # reservoir-sustain strength (SCENARIO knob: small->EZH2-gated/leaky 'B'; large->reservoir-sustained 'A');
            # K_el/n_el = elongation gate; k_onP/k_offP = proximal kinetics; f_P/K_P/n_P = extra repression strength.
            assert not with_mark_amplifier, "with_proximal_distal supersedes with_mark_amplifier; do not enable both (double-wrap)"
            assert "species Mk in Cell" in m, "proximal/distal needs the H3K27 chain reservoir (Mk)"
            assert "Cdk6_synthesis: => cdk6; Cell*(" in m, "proximal/distal needs with_cdk6_gli (CDK6 synthesis)"
            assert m.count("CycD1_transcription: => Cd_mRNA; (") == 1, "proximal/distal: CyclinD1 transcription not in expected form"
            m = m.replace("\nend",
                "\n  a_P = 0.10; n_anch = 4; K_el = 0.20; n_el = 6;   // proximal load = reservoir a_P*Mk^n_anch gated by elongation-OFF; elong_gate = E2f Hill (K_el,n_el). DEFAULT = transient scenario 'C'."
                "\n  k_onP = 0.05; k_offP = 0.1; basal_offP = 0.002;   // proximal PRC2 loading / eviction kinetics (slow enough that transient G1 E2f-dips do not load P_prox)"
                "\n  f_P = 0.80; K_P = 0.35; n_P = 4;   // EXTRA proximal repression of CDK6+CyclinD1: Hill on P_prox (f_P = floor = max extra repression). SCENARIO knob: f_P 0.80/a_P 0.10 = TRANSIENT dwell ~148h (30/32); f_P 0.30/a_P 0.50 = PERMANENT-LOCK (28/32, data-excluded like the chromatin latch); f_P 0.90/a_P 0.05 = LEAKY ~89h. Boundary sharp at a_P~0.11."
                "\n  species P_prox in Cell; P_prox = 0.0;   // proximal PRC2 occupancy [0,1] = reservoir-sustained repressor in arrest"
                "\n  elong_gate := E2f^n_el/(K_el^n_el + E2f^n_el);   // Pol2 elongation / proliferation state: ~1 cycling, ->0 in sustained arrest"
                "\n  w_cd_prox = 1.0; w_cdk6_prox = 1.0;   // per-ARM strength of the proximal repression (1 = full, default; 0 = that arm un-repressed). Lets the H3K27me3->CyclinD1 vs ->CDK6 repression be weakened independently (JP mechanism decomposition)."
                "\n  prox_amp := f_P + (1 - f_P)/(1 + (P_prox/K_P)^n_P);   // extra repression from proximal PRC2 occupancy (1 when P_prox~0)"
                "\n  prox_amp_cd := 1 - w_cd_prox*(1 - prox_amp);     // CyclinD1 arm (w_cd_prox scales repression; =prox_amp at w=1)"
                "\n  prox_amp_cdk6 := 1 - w_cdk6_prox*(1 - prox_amp); // CDK6 arm"
                "\n  P_load: => P_prox; Cell*k_onP*a_P*Mk^n_anch*(1 - elong_gate)*(1 - P_prox);   // reservoir re-loads proximal when Pol2 stops elongating"
                "\n  P_evict: P_prox => ; Cell*(k_offP*elong_gate + basal_offP)*P_prox;   // active elongation (Pol2) + basal evict the proximal complex\nend", 1)
            m = m.replace("Cdk6_synthesis: => cdk6; Cell*(", "Cdk6_synthesis: => cdk6; Cell*prox_amp_cdk6*(")
            m = m.replace("CycD1_transcription: => Cd_mRNA; (", "CycD1_transcription: => Cd_mRNA; prox_amp_cd*(")
            assert "P_load:" in m and "Cell*prox_amp_cdk6*(" in m and "Cd_mRNA; prox_amp_cd*(" in m, "proximal/distal wiring failed"

        if with_cd_hyper_escape:
            # ===== ESCAPE from the p27 / bistable-OFF G0 (JP 2026-08-07 v44 transient-G0 rebuild). EXPLORATION, default OFF. =====
            # Diagnostic (birthp27_probe / escape trajectory): a cell born with high p27 LOCKS PERMANENTLY -- NOT because
            # p27 stays high (free p27 clears to ~0, CDK2 is not p27-inhibited), but because E2f is trapped on MONO-
            # phosphorylated Rb (RbmE2f) and the two-step Rb lets ONLY CyclinE/A hyper-phosphorylate Rb. So accumulating
            # CyclinD1 (Cd~13) cannot complete Rb-P / release E2f -> no CyclinE bootstrap -> permanent arrest, never a
            # transient dwell. FIX (Yang 2020 eLife 44571 / Chung 2019 Mol Cell 31543423: CDK4/6 activity ALONE drives
            # commitment, even in CyclinE/A-quadruple-null MEFs): let high CyclinD-CDK4/6 contribute WEAKLY to Rb HYPER-
            # phosphorylation (w_cd_hyper * the CyclinD-CDK4/6 saturating activity). A high-p27 daughter then DWELLS in G0
            # while CyclinD1 accumulates, then RE-ENTERS once CyclinD-CDK4/6 crosses the hyper-P threshold -> a TRANSIENT
            # G0. w_cd_hyper kept small so the normal-cycle two-step bistability + Shh-dependence are preserved.
            # _dact = the CyclinD-CDK4/6 saturating activity, matched to the CURRENT drive/brake so the escape is
            # consistent with the integrated model: option-B drops CdP21 from the drive; with_cdk6_gli sinks the INK4
            # brake by unbound CDK6. (JP 2026-08-09: fixed from the option-A-only original for the integrated model.)
            _drive = "(Cd + w_Cd2*Cd2)" if with_p27_optionB else "(Cd + w_Cd2*Cd2 + CdP21)"
            _brake = "(1 + (w_p18*p18_prot + w_p19*p19_prot)*K_cdk6_sink/(K_cdk6_sink + cdk6))" if with_cdk6_gli else "(1 + p18_prot + p19_prot)"
            _dact = f"kPhRbCd*{_drive}/(K_CdRb*{_brake} + {_drive})"
            assert m.count("(kPhRbCe*Ce + kPhRbCa*Ca)*Rbm)")==1 and m.count("(kPhRbCe*Ce + kPhRbCa*Ca)*RbmE2f)")==1, "escape: hyper-P reactions not in expected form (needs default cdki-on/option-A)"
            m = m.replace("(kPhRbCe*Ce + kPhRbCa*Ca)*Rbm)", f"(kPhRbCe*Ce + kPhRbCa*Ca + w_cd_hyper*{_dact})*Rbm)")
            m = m.replace("(kPhRbCe*Ce + kPhRbCa*Ca)*RbmE2f)", f"(kPhRbCe*Ce + kPhRbCa*Ca + w_cd_hyper*{_dact})*RbmE2f)")
            m = m.replace("\nend", "\n  w_cd_hyper = 0.15;   // weak CyclinD-CDK4/6 -> Rb hyper-P: escape from the p27/OFF-lock G0 (transient-G0 re-entry; Yang2020 CDK4/6-alone commits). PLACEHOLDER.\nend", 1)
            assert "w_cd_hyper*" in m, "escape: w_cd_hyper not wired"

        _ts_bake = {
            # 2026-07-27 CELL-TYPE SPLIT BAKED (JP-approved): GNP 16-17h / MB ~23.5h CORE cycle, decoupled commitment.
            # GNP cycle = data (Nakashima 15.9h + Contestabile 16.25h); MB longer via the SAME k_mu_cki slowdown driven
            # by MB's higher INK4 (p16+p18) -- the cycle-length split EMERGES from the CDKi, not a separate knob. Paired
            # with decouple_commit=True (default flipped): the G0/commit R-point is gated on CyclinD1/CDKi + birth-p27
            # (growth-INDEPENDENT), which restores cKO Shh-dependence + the strong EZH2 governor at the faster cycle.
            # Params from recal_16h_v3 (cKO-hard-constrained search). 26/29 (fails: CyclinD1 GNP+HHi = the basal-cut cost
            # of cKO; MB HU S+G2 = the structural HU limit). See memory mb-celltype-transient-g0-parameterization.
            'mu': 0.000769, 'k_mu_cki': 0.3619,   # k_mu_cki re-fit for the baked dynamic-CDKI model (was 0.307)
            # 2026-07-30 DYNAMIC-CDKI BAKED (recal_cdki_clean_best): individual expressed CDKI species now default.
            # p27+INK4(p18/p19)+CdP21 carry CDK regulation; p21 low-activity, p57 off. 30/32 (fails HU pair), Shh-dep.
            'kSeqCd': 2.13243, 'kRelCd': 0.01647, 'kDeKPC': 0.02022, 'w_ink4': 3.6204, 'kSyP21': 0.00116,
            # 2026-07-28 TWO-CYCLIN (JP): engage CyclinD2 (separate, Hh-buffered D-cyclin; CCND2 dominant + drops
            # only ~35% under vismo). Cd2 = basal + small Gli-driven + MB elevation (Cd2_expr, set per-condition);
            # feeds the CDK4/6->Rb drive with w_Cd2. Calibrated: D2 MB/GNP 2.39, GNP+HHi 0.66 (exact). w_Cd2<=0.2
            # keeps Shh-gating (D2's Gli-independent drive stays below the R-point). See cyclind1-d2-two-cyclin-model.
            'k_Cd2_bas': 2.0, 'k_Cd2_Gli': 2.49, 'K_Cd2_Gli': 0.3, 'w_Cd2': 0.079448,   # w_Cd2: two-cyclin recal
            # 2026-07-19 Cdh1-EZH2 re-bake (JP-approved; optimize_cdh1_ezh2). Adds APC/Cdh1(C1)-GATED EZH2
            # DEGRADATION (kDeEZ_C1) on top of the DILUTION convention (with_ezh2_conc=False, default flipped) so an
            # ELONGATED / HU-ARRESTED S ACCUMULATES EZH2 (exp-verified by JP): HU EZH2-in-S boost 1.009->1.24, EZH2
            # protein G2/G0 1.04->1.86, EZH2 G0/cycling 0.44->0.62 -- all previously ~flat. Mechanism: LOW baseline
            # kDeEZ (committed/S/arrest cells STABLE, EZH2 banks over the long arrest) + kDeEZ_C1*C1 clears EZH2 in G0
            # (Cdh1 active). Biologically grounded: EZH2 is a bona fide APC/C-Cdh1 substrate. Eraser still OFF
            # (k_jmjd3_gli=0); fold via direct Gli/MYCN transcription; vismodegib 0.128; mark half-life ~2.8h preserved.
            # Still 27/28; lone miss is MB HU G2 (0.352), the SAME structural near-miss (see v44-hu-s-g2-structural-limit).
            # Cell-cycle machinery (M_commit, kPhRbCd, Rb, f_commit_carry, M_size, g_prc2) UNCHANGED. User `params` override.
            'M_commit': 2.476412980954672, 'kPhRbCd': 0.1403, 'kDeP21Cd': 0.6300278361307281,   # kPhRbCd, K_CdRb: two-cyclin recal (split was 0.1433 / 0.3284)
            'kDsRbmE2f': 1.8821064877621487, 'K_CdRb': 0.4032, 'f_commit_carry': 0.37450197206637914,
            'kSyDna': 0.046863027361587935, 'M_size': 3.140622748974385, 'kEZbas': 0.0006220377408563062,
            'kEZbas_Cd': 0.00031843672436329825, 'kEZE2f': 0.01226393864518667, 'Kez_cd': 9.17574845440131,
            'kDeEZ': 0.0002766221763090641, 'kDeEZ_C1': 0.0003858762128411371, 'K_E2f_EZ': 0.4497729841633607,
            # 2026-07-21 (JP-approved): re-enable a MASSIVELY NERFED Gli->Jmjd3/Kdm6b eraser (was 0 by default; the
            # optimize_chain string default 0.108) + three more-biological mark mechanisms, all added default-preserving
            # then engaged here: n_rw (cooperative read-write, gap #1), dil_frac=0/k_dil_S (continuous S-phase dilution
            # spread over S so an elongated S recovers concurrently, gap #2), k_h33 (H3.3/HIRA transcription-coupled
            # replication-INDEPENDENT me3 loss, gap #3). Values chosen mild to preserve the 28/29 validation.
            'k_jmjd3_gli': 0.005, 'n_rw': 2, 'dil_frac': 0.0, 'k_dil_S': 0.693, 'k_h33': 0.003,
            'a0_prc2': 0.0028749753507882796,
            'a_rw_prc2': 0.0022657109194688194, 'g_prc2': 0.0704913572494475, 'K_prc2': 0.006856946735209077,
            'n_prc2': 4.015744557919901, 'f0_prc2': 0.17482554666158287, 'del_mk': 0.00044000677598027846,
            'KmHU_fire': 0.2797289308652902,
            # 2026-07-20 Skp2 MITOGEN-DOSE induction (JP scRNA + MB table: Skp2 HIGHER in MB ~2.3x -- the model's
            # E2F-only + Cdh1-degradation made MB Skp2 LOW because MB has high Cdh1). Re-optimized -> Skp2 MB/GNP 2.31,
            # HU-G2 recovered to 0.357 (structural floor), 28/29. kSySkp2_Cd=0 would disable it.
            'kSySkp2_Cd': 0.24457065227870578, 'K_Skp2_Cd': 10.296077904941832,
            'k_Cd_tx_Gli_max': 0.17494, 'k_Cd_tx_MYCN': 0.062,   # 2026-08-01 p21-un-map re-cal (was 0.152118 / 0.031157): the un-map raises MB-specific EZH2-in-S -> more MB CyclinD1 repression -> MB/GNP fold drops ~7%. Restored via the MB-specific MYCN driver (spares GNP+HHi, unlike Gli) + a small Gli bump. Species-correct targets (EZH2m/Gli1_mRNA). 30/32.
            'k_Cd_tx_basal': 0.000529,   # two-cyclin recal: RAISED (D2 now floors proliferation) -> D1 GNP+HHi 0.055->0.10 passes. tx scaled 1/800 (with k_Cd_mRNA_deg 1/800) so the mRNA is a SLOW reservoir; THEN halved 2026-07-24: lowers the mitogen-INDEPENDENT CyclinD1 floor so that even FULL de-repression (genetic f0_prc2=1 / EZH2 cKO) stays Shh-DEPENDENT (no division at SHH=0, matches JP's cKO: transcript UP but no Shh-independent division + differentiates). Governed threshold only 0.16->0.18, de-repression fold preserved (~3.5x full / ~1.6x EZH2i), validation-neutral 28/29. Shh-dependence now comes from the LOW basal (CyclinD1 needs Gli-drive to commit), NOT from EZH2 repression -> EZH2 = threshold/level GOVERNOR, not the Shh on/off switch
            # 2026-07-20 graded G1/cell-cycle LENGTHENING (JP): subthreshold mitogen SLOWS growth -> Tc lengthens
            # GRADUALLY (no transient G0) -> longer low-E2F G1 -> slow EZH2 (S-phase) response = the buffering.
            # K_g1len (0.6) sits BELOW the GNP/MB CyclinD1, so it is LATENT at the validation conditions (still 27/28,
            # only GNP period 22.6->~23h) and engages only under mitogen withdrawal/perturbation (Tc up to ~2.5x at
            # deep subthreshold). mu_min_frac=1 would disable it. Free params (no withdrawal target yet) -- tune when
            # withdrawal/EdU-plateau data land.
            'mu_min_frac': 0.4, 'K_g1len': 0.6,
            # 2026-07-21 CyclinD1 = the Hh MEMORY CARRIER (Ho, Tsai & Stearns 2020, Curr Biol 30:2829; JP acks), and
            # the memory is the cumulative CyclinD1 TRANSCRIPT reservoir (JP mechanism). The PROTEIN is labile/fast
            # (k_Cd_deg=1.0, a gated readout that forms rapidly when there is lots of mRNA); the mRNA is the SLOW
            # RESERVOIR, driven UP by Gli and DOWN by EZH2 -- an INCOHERENT FEEDFORWARD (Gli raises both the transcript
            # and EZH2, so when mitogen falls EZH2's repression lifts and COMPENSATES, holding the transcript up
            # longer; the more EZH2, the more de-repression 'room'). Implemented by slowing the mRNA (k_Cd_mRNA_deg
            # 0.8->0.001 = ~11-12h reservoir) with the transcription (k_Cd_tx_*, above) scaled 1/800 to HOLD both mRNA
            # and protein LEVELS (so g_prc2/H3.3 terms reading Cd_mRNA are unaffected; fold ratios preserved).
            # Validation-neutral (28/29); after Hh-off the mRNA reservoir t1/2 ~13h -> cell COASTS ~1 more cycle
            # (the paper's previous-cycle restriction point). Persistence tunable via k_Cd_mRNA_deg.
            'k_Cd_mRNA_deg': 0.001}
        # tolerate flag combos where some params are absent (e.g. with_h3k27_memory has no chain params)
        m = _apply_overrides(m, {k: v for k, v in _ts_bake.items() if (k + ' = ') in m})
    if hu is not None:
        m = m.replace("HU = 0;", f"HU = {hu};")
    if params:
        m = _apply_overrides(m, params)
    return m


if __name__ == "__main__":
    m = build_model_v44()
    print("v44 model built:", len(m), "chars")
    for tok in ["MPF", "Chk1 := aRc", "Cdc25a", "Wee1a", "vfork :=", "kSyDna*vfork*aRc",
                "E_div: at (MPF > MPF_div)", "size_gate :=", "*commit_gate", "mass = mass/2", "Skp2_synthesis:", "P21 = P21_div"]:
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
