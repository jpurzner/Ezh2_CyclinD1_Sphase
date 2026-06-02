"""
Model v43 - Intra-S Checkpoint extension of the v42 MYCN-EZH2-CyclinD1-HH model.

This module does NOT modify the published v42 model. It imports the v42 Antimony
string and programmatically injects an algebraic hydroxyurea (HU) effect on S-phase
PROGRESSION FLUX. The Gerard-Goldbeter cell cycle core is left untouched except for a
single surgically added multiplicative factor phi(HU) on the CycA/CDK2 (Ma) activation
step -- the CDK2-Tyr15 dephosphorylation that drives S-phase progression.

Biological rationale
--------------------
HU inhibits ribonucleotide reductase -> dNTP pools collapse -> replication forks slow.
dNTP pools equilibrate in minutes (fast) relative to the multi-hour cycle, so the
effective fork velocity / replication competence is modeled as an INSTANTANEOUS
algebraic function of HU (quasi-steady-state; no dNTP ODE):

    phi(HU) = phi_min + (1 - phi_min) / (1 + (HU / K_HU_phi)^h_HU)

    phi in [phi_min, 1]:  1 = uninhibited, phi_min = residual progression at saturating
                          HU. phi_min > 0 GUARANTEES graded slowing, never full arrest
                          (matches the data: HU prolongs S, cells still divide).
    K_HU_phi : cellular IC50 for HU (model units; HU=1.0 ~ 10 uM in the experiment).
    h_HU     : Hill coefficient (~2-4): switch-like S block near threshold.

Coupling -- phase-gated cell-cycle CLOCK (NOT a single-flux scaling)
-------------------------------------------------------------------
The Gerard-Goldbeter relaxation oscillator's period is ROBUST to scaling any single
cyclin activation/degradation/synthesis flux by a constant factor: the bistable
Cdc25/Wee1 switches fire at the same time regardless of linear rate, so a multiplicative
phi on Ma_activation (or Mb_activation, Cdc20, ...) only rescales cyclin AMPLITUDE and
leaves the period and S-dwell flat (verified: simulations/diag_phi_levers.py). Only the
GLOBAL timescale eps moves the period. We therefore couple phi to a PHASE-GATED clock:

    sphase := Ma^n_S / (K_S_phi^n_S + Ma^n_S)        # ~1 while replicating (S/G2), ~0 in G1/M
    eps    := eps0 * (1 - sphase * (1 - phi(HU)))    # clock slowed toward eps0*phi in S only

So replication stress (low phi) slows the cell-cycle ENGINE specifically within S, which
prolongs the S/G2 dwell, lengthens the period, lowers the division rate -- graded, never
arrested (phi >= phi_min > 0) -- WITHOUT throttling any one switch relative to its antagonist
(so no free-CycB inflation / phase-detector artifact). This is the canonical intra-S /
ATR-CHK1 checkpoint: an unfinished, stressed S-phase holds/slows progression toward mitosis.
At HU=0, phi=1 -> eps=eps0 everywhere -> the model is numerically identical to v42 in EVERY
phase, so ALL HU=0 invariants (within-cycle EZH2 gradient, GDC reduction, mean level) are
preserved automatically; only the HU response is shaped by (K_S_phi, n_S, phi params).

EZH2 response is EMERGENT: there is NO HU term in the EZH2 equation. EZH2 rises under
HU purely because cells spend longer in the high-synthesis (Me+Ma-gated, S-phase)
window. This is a prediction of the model, not a fitted coupling.

INVARIANT: at HU = 0, phi = 1.0 exactly, and the model is numerically identical to v42
(plus the EZH2 rewire). v43 must reproduce the v42 validation when HU = 0.
"""

from src.build_model_v42_mycn import build_model_v42


# Default HU -> S-progression parameters (algebraic phi(HU); calibrated in Phase 3).
# These REPLACE the former ATR/CHK1/D_stress block and all CHK1->Cdc25/Cdc20 levers,
# which were abandoned: every lever either arrested the Gerard-Goldbeter oscillator
# (no phi_min floor -> CHK1 could drive activation toward 0 -> bistable collapse) or was
# a phase-detector artifact. The bounded phi(HU) with a phi_min floor is the principled
# fix -- a graded, self-limiting slowdown of S-progression that never fully arrests.
HU_PHI_DEFAULTS = dict(
    phi_min=0.25,       # residual replication competence / fork velocity at saturating HU.
                        # >0 GUARANTEES graded slowing, never arrest. Fit from residual
                        # EdU/BrdU incorporation (or residual doubling) at high HU.
    K_HU_phi=0.6,       # cellular IC50 for HU (model units; HU=1.0 ~ 10 uM experiment).
    h_HU=3.0,           # Hill coefficient (~2-4): switch-like S block near threshold.
    # --- S-phase clock gate (where phi bites): see module docstring "Coupling" ---
    K_S_phi=0.08,       # CycA(Ma) half-on for the replication/S window (sphase gate).
    n_S=8,              # steepness of the S gate (sharp on/off of the clock slowdown).
)
# Calibration note (simulations/diag_phi_*.py): these give an EMERGENT S-phase EZH2
# boost ~1.3x at HU=1.0 (the 10 uM anchor; target 1.31x, p=0.005) with period ~19.6h ->
# ~27h, divisions 13 -> ~8 over 250h (graded, never arrested), S-phase ~2.5-3x longer.
# The boost is purely a consequence of longer dwell in the high-synthesis S window; there
# is NO HU term in the EZH2 equation. Reaching exactly 1.31x ties the EZH2 response to a
# substantial S-prolongation -- a falsifiable model prediction. Dial phi_min DOWN (stronger
# S-arrest) for a larger boost, UP for a milder one.


# EZH2 transcription rewire (Phase 3 design decision: track time-in-S, not pRBpp spike).
#
# The published v42 gates EZH2 transcription on pRBpp, a CycE/CDK2-driven spike at the
# G1/S transition that is high only ~3% of the cycle. That makes EZH2 effectively flat
# (tau~50h) and structurally unable to (a) show the measured within-cycle gradient
# (G0<G1<S<G2) or (b) integrate S-phase DURATION into EZH2 protein. We therefore rewire
# transcription to track Ma = active CycA/CDK2, the operational S/G2 marker (high ~16% of
# the cycle), and speed protein turnover so EZH2 ramps within a cycle. This makes the
# plan's central hypothesis ("longer S -> more EZH2") mechanistically realisable while
# preserving v42's BETWEEN-condition EZH2 behaviour (mean cycling level, G0/cycling ratio,
# GDC reduction, EZH2i->CycD1). Calibrated in Phase 3.
#
# Set rewire_ezh2=False in build_model_v43() to recover the exact v42 EZH2 wiring.
#
# Phase-3 CALIBRATED values (used in simulations/validate_v43.py). The shape params
# (K_Ma_EZH2, K_Me_EZH2, w_Me, k_EZH2_deg, k_EZH2_translation) give the within-cycle EZH2
# gradient and preserve the v42 between-condition invariants. There is NO HU/stress term
# here: the HU -> EZH2 response is EMERGENT (longer S-dwell under phi(HU); see module 3).
# NOTE: the transcription floor/amplitude (k_EZH2_mRNA_synth_basal=0.02,
# k_EZH2_mRNA_synth_E2F=3.5) live in the inherited v42 string and are set at RUNTIME
# (see calibrate_ezh2_rewire.run); they are NOT part of this dict.
EZH2_REWIRE_DEFAULTS = dict(
    K_Ma_EZH2=0.35,          # half-saturation of CycA/CDK2 (Ma) for EZH2 transcription
    K_Me_EZH2=0.35,          # half-saturation of CycE/CDK2 (Me) for EZH2 transcription
    w_Me=0.65,               # weight of the Me (S-ONSET) arm in the EZH2 transcription gate.
                             # The gate is a weighted SUM of two Michaelis-Menten arms:
                             #   w_Me * Me/(K_Me+Me) + (1-w_Me) * Ma/(K_Ma+Ma)
                             # Ma (CycA/CDK2) peaks late (S/G2); on its own the EZH2 protein
                             # keeps decaying from the prior G2 peak through G1 into early S,
                             # so the within-cycle minimum lands IN S and the gradient is
                             # non-monotone (S<G1). Me (CycE/CDK2) marks the G1/S transition;
                             # the Me arm restarts EZH2 transcription at S-ONSET so protein is
                             # already rising by S -> restores monotone G0<G1<S<G2. BOTH arms
                             # collapse at the arrested/quiescent floor (Ma,Me ~0), so blocking
                             # the cycle (GDC0449/HHi, serum-starve) drops the gate toward 0 ->
                             # restores the v42 GDC-reduction invariant. w_Me=0 recovers the
                             # pure Ma gate; w_Me~0.65 balances the gradient against GDC.
    k_EZH2_deg=0.20,         # faster EZH2 turnover (tau~5h) so within-cycle gradient appears
    k_EZH2_translation=0.80, # raised from v42's 0.12 to hold mean cycling EZH2 ~0.49
)


_CHECKPOINT_BLOCK = """
// ===================================================================
// MODULE 3: INTRA-S HU EFFECT (algebraic, quasi-steady-state)
// ===================================================================
// HU inhibits ribonucleotide reductase -> dNTP depletion -> slower
// replication fork velocity. dNTPs equilibrate in minutes (fast) vs the
// multi-hour cycle, so fork competence phi is an INSTANTANEOUS algebraic
// function of HU (no dNTP ODE). phi in [phi_min, 1]; phi_min > 0 -> graded
// slowing, never full arrest. phi does NOT scale a single cyclin flux (the
// GG period is robust to that); instead it slows the cell-cycle CLOCK (eps)
// specifically within the replication window via the sphase gate below
// (see the build step that redefines eps). At HU=0, phi = 1.0 exactly ->
// eps = eps0 everywhere -> identical to v42 (+ EZH2 rewire) in every phase.

// Hydroxyurea dose (0 = none, 1 = saturating ~10 uM in experiment)
HU = 0.0;
phi_min = {phi_min};
K_HU_phi = {K_HU_phi};
h_HU = {h_HU};

// Effective replication competence / fork velocity (algebraic assignment)
phi_HU := phi_min + (1 - phi_min) / (1 + (HU / K_HU_phi)^h_HU);

// Replication / S-phase progress gate: ~1 while replicating (CycA/Ma high),
// ~0 in G1/M. The clock slowdown (eps redefinition) bites ONLY here.
K_S_phi = {K_S_phi};
n_S = {n_S};
sphase := Ma^n_S / (K_S_phi^n_S + Ma^n_S);

"""

def _apply_ezh2_rewire(model, rewire_params):
    """Rewire EZH2 transcription from the pRBpp spike to Ma (CycA/CDK2, S/G2 marker).

    Pure string surgery on the inherited v42 string; v42 itself is never modified.
    Each replacement asserts its target exists so a v42 structure change fails loudly.
    """
    p = dict(EZH2_REWIRE_DEFAULTS)
    unknown = set(rewire_params or {}) - set(EZH2_REWIRE_DEFAULTS)
    if unknown:
        raise ValueError(f"Unknown EZH2 rewire parameter(s): {sorted(unknown)}")
    p.update(rewire_params or {})

    # 1. Add the new half-saturation constant next to the (now-unused) pRBpp one.
    if "K_pRBpp_EZH2 = 0.03;" not in model:
        raise RuntimeError("K_pRBpp_EZH2 declaration not found in expected v42 form.")
    model = model.replace(
        "K_pRBpp_EZH2 = 0.03;",
        "K_pRBpp_EZH2 = 0.03;\nK_Ma_EZH2 = {K_Ma_EZH2};".format(**p),
    )

    # 2. Speed EZH2 protein turnover so it ramps within a cycle (tau ~10h vs v42's ~50h).
    if "k_EZH2_deg = 0.02;" not in model:
        raise RuntimeError("k_EZH2_deg declaration not found in expected v42 form.")
    model = model.replace("k_EZH2_deg = 0.02;",
                          "k_EZH2_deg = {k_EZH2_deg};".format(**p))

    # 3. Raise translation to hold mean cycling EZH2 ~0.49 despite faster degradation.
    if "k_EZH2_translation = 0.12;" not in model:
        raise RuntimeError("k_EZH2_translation declaration not found in expected v42 form.")
    model = model.replace("k_EZH2_translation = 0.12;",
                          "k_EZH2_translation = {k_EZH2_translation};".format(**p))

    # 4. Swap the transcription driver: pRBpp spike -> a weighted SUM of an S-ONSET arm
    #    (Me = CycE/CDK2) and an S/G2 arm (Ma = CycA/CDK2). The Me arm restarts EZH2
    #    transcription at the G1/S transition so protein is already rising by S (monotone
    #    G0<G1<S<G2); the Ma arm carries it through G2. Both arms collapse at the arrested
    #    floor, restoring the GDC-reduction invariant. There is NO HU/stress term: the HU ->
    #    EZH2 response is EMERGENT -- under HU, phi(HU) prolongs the high-CycA S window, so the
    #    Me+Ma gate stays active longer and EZH2 accumulates more (a model prediction, not a
    #    fitted coupling). At HU=0, phi=1 -> identical to the DMSO baseline.
    tx_old = ("EZH2_mRNA_synthesis: -> EZH2_mRNA; k_EZH2_mRNA_synth_basal + "
              "k_EZH2_mRNA_synth_E2F * E2F / (K_E2F_EZH2 + E2F) * "
              "pRBpp / (K_pRBpp_EZH2 + pRBpp);")
    tx_new = ("EZH2_mRNA_synthesis: -> EZH2_mRNA; k_EZH2_mRNA_synth_basal + "
              "k_EZH2_mRNA_synth_E2F * E2F / (K_E2F_EZH2 + E2F) * "
              "(w_Me * Me / (K_Me_EZH2 + Me) + "
              "(1 - w_Me) * Ma / (K_Ma_EZH2 + Ma));")
    if tx_old not in model:
        raise RuntimeError("EZH2_mRNA_synthesis reaction not found in expected v42 form.")
    model = model.replace(tx_old, tx_new)

    # 4b. Declare the new constants (next to the half-saturation constant).
    model = model.replace(
        "K_Ma_EZH2 = {K_Ma_EZH2};".format(**p),
        "K_Ma_EZH2 = {K_Ma_EZH2};\nK_Me_EZH2 = {K_Me_EZH2};\n"
        "w_Me = {w_Me};".format(**p),
    )

    return model


def build_model_v43(checkpoint_params=None, rewire_ezh2=True, ezh2_rewire_params=None):
    """Build v43 = v42 + intra-S checkpoint (+ optional EZH2->Ma transcription rewire).

    Parameters
    ----------
    checkpoint_params : dict or None
        Overrides for any of HU_PHI_DEFAULTS (phi_min, K_HU_phi, h_HU). Name kept as
        `checkpoint_params` for harness backward-compatibility.
    rewire_ezh2 : bool
        If True (default), rewire EZH2 transcription to track Ma (CycA/CDK2) so EZH2
        integrates S-phase duration. If False, keep the exact v42 pRBpp wiring (used
        for the negative-control comparison and to confirm non-destructiveness).
    ezh2_rewire_params : dict or None
        Overrides for any of EZH2_REWIRE_DEFAULTS (only used when rewire_ezh2=True).

    Returns
    -------
    str : Antimony model string.
    """
    params = dict(HU_PHI_DEFAULTS)
    if checkpoint_params:
        unknown = set(checkpoint_params) - set(HU_PHI_DEFAULTS)
        if unknown:
            raise ValueError(f"Unknown HU/phi parameter(s): {sorted(unknown)}")
        params.update(checkpoint_params)

    model = build_model_v42()
    if rewire_ezh2:
        model = _apply_ezh2_rewire(model, ezh2_rewire_params)
    checkpoint_block = _CHECKPOINT_BLOCK.format(**params)

    # --- 1. Inject checkpoint species/params/reactions before the INITIAL CONDITIONS marker ---
    init_marker = "// INITIAL CONDITIONS"
    if init_marker not in model:
        raise RuntimeError(
            "v42 structure changed: '// INITIAL CONDITIONS' marker not found. "
            "Checkpoint injection aborted to avoid corrupting the inherited core."
        )
    # Insert just before the comment banner line that precedes the marker.
    insert_at = model.rfind("// ===", 0, model.index(init_marker))
    model = model[:insert_at] + checkpoint_block + "\n" + model[insert_at:]

    # --- 2. Phase-gated clock coupling: redefine the global timescale eps so that the
    #        cell-cycle ENGINE is slowed toward eps0*phi(HU) specifically within the
    #        replication (S/G2) window (sphase ~ 1), and left at eps0 elsewhere. This is
    #        the ONLY HU coupling. A single-flux scaling (e.g. Ma_activation * phi) does
    #        NOT slow the GG oscillator -- its period is robust to cyclin amplitude (see
    #        module docstring) -- so we modulate the clock instead. eps is a global
    #        time-scaling PARAMETER, not a core kinetic mechanism; gating it by phi(HU)
    #        injects the checkpoint as an external regulator without altering the GG
    #        reaction network. At HU=0, phi=1 -> eps = eps0 in all phases (identical to v42).
    eps_old = "eps = 150;"
    eps_new = ("eps0 = 150;\n"
               "eps := eps0 * (1 - sphase * (1 - phi_HU));")
    if eps_old not in model:
        raise RuntimeError("eps = 150; declaration not found in expected v42 form.")
    model = model.replace(eps_old, eps_new)

    # Rename the model so SBML/Antimony id differs from v42.
    model = model.replace("model ezh2_cyclind1_v42", "model ezh2_cyclind1_v43_checkpoint")

    return model


if __name__ == "__main__":
    m = build_model_v43()
    print("v43 checkpoint model built:", len(m), "characters")
    # Sanity: confirm injected tokens are present and the old hacks are gone.
    for token in ["phi_HU", "phi_min", "K_HU_phi", "h_HU", "K_S_phi", "n_S",
                  "sphase := Ma^n_S", "eps := eps0 * (1 - sphase * (1 - phi_HU))",
                  "ezh2_cyclind1_v43_checkpoint"]:
        assert token in m, f"missing {token}"
    for gone in ["D_stress", "ATR_active", "CHK1_active", "Ki_CHK1", "k_EZH2_stress",
                 "Pa * phi_HU * eps", "eps = 150;"]:
        assert gone not in m, f"stale token still present: {gone}"
    print("phi(HU) phase-gated clock module present; single-flux/stress hacks removed.")
    print("Defaults:", HU_PHI_DEFAULTS)
