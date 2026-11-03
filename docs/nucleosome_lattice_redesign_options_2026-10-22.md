# Redesign options: nucleosome-resolved promoter module (2026-10-22, draft for discussion)

Replaces the mean-field `Mk` + `P_prox` two-compartment module with an explicit
nucleosome lattice, keeping the core/distal separation. Written against
`two_compartment_prc2_2026-08-09.md`, `stochastic_amplifier_2026-10-22.md`, the
`transcriptional-model` MODEL.md kernel, and the spatial-model literature in
`chipseq_norm/notes/prc2-boundary-modelling.md`.

---

## 0. The structural blocker, first

**v44 cannot be restarted mid-arrest.** Chunked integration dies at t~6250 min
(CV_CONV_FAILURE), reproducibly, and identically with
`with_proximal_distal=False`. No tolerance or integrator setting fixes it.

Any nucleosome lattice needs per-step state updates, which means restarts. So a
lattice **cannot** live inside the Antimony/roadrunner ODE, and co-simulation by
interleaved `rr.simulate()` calls is blocked by the same bug. This is the single
biggest architectural constraint and it should drive the choice below, not the
biology.

Two ways out:
- **(a) One-way coupling first.** Drive the lattice with *recorded* trajectories
  of EZH2, Gli1, nascent Ccnd1 and S-phase timing from the existing ODE. No
  restarts, no porting. Validates the lattice against ChIP geometry before the
  loop is closed. Cheap, and most of the payoff is here.
- **(b) Port the cell-cycle core out of Antimony** into a plain Python RHS
  (scipy `solve_ivp`, LSODA/Radau + event handling). Removes the restart bug
  permanently and makes closed-loop co-simulation trivial. Real cost, but the
  stiffness problem is already biting independently of this redesign.

Recommendation: do (a) now, treat (b) as the thing that unlocks closed loop.

## 1. What the discussions already contain (Gli1 -> KDM6B/JMJD3)

Already implemented, in `build_model_v44_heldt.py`:

```
del_mk = 0.0015;  k_jmjd3_gli = 0.0478          # shared-mark version
Mk_demeth_jmjd3: Mk => ;  Cell*k_jmjd3_gli*Gli1*Mk
# me-chain version:
me3_demeth: Mk => m2_me;  Cell*(del_mk + k_jmjd3_gli*Gli1)*Mk
me2_demeth: m2_me => m1_me; Cell*(d_me + k_jmjd3_gli*Gli1)*m2_me
```
with `k_jmjd3_gli` = 0.108 (chain fit), 0.0478 (shared fit), 0.005 (nerfed
default). Cited to Shi 2014 (ncomms6425). Framed as **WRITER = EZH2 (cycle-gated)
vs ERASER = Gli (mitogen-gated)** — a race, with the eraser mitogen-gated so it
strips the mark in high-Gli MB and gives ChIP MB < GNP. Supporting code:
`simulations/optimize_jmjd3_chip.py`, `simulations/fig_v44_gli_kdm6b_eraser.py`.

**The thing that makes a lattice worth building:** the ChIP repo already has
**GNP_P7_KDM6B** and **GNP_P7_EZH2ko_Gli1-Flag**. In a mean-field model those
tracks are unusable — the eraser is a scalar rate. In a lattice they become a
*measured* occupancy profile psi_KDM6B(x), which is exactly the term Newar/Jost
had to put in as a static open-loop guess (their stated limitation: "UTX is an
open-loop fixed input, not dynamic"). Here it would be measured in x AND dynamic
in t through Gli1. That is a genuine advance over the published models, and it
is the strongest argument for going nucleosome-resolved.

## 2. Geometry: how many nucleosomes, and from what

Lattice length n = N_core + 2*N_flank, nucleosome pitch ~190 bp.

- **N_core** = width of the H3K27me3-depleted core / 190. Measured per gene from
  the ChIP profile.
- **N_flank** = decay length of the flanking me3 domain / 190, taken per side
  because the real profile is asymmetric. `transcriptional-model` MODEL.md
  already fits this: **lambda_up ~ 248 bp, lambda_dn ~ 4.7 kb** (priors were
  lambda_1 0.5-1 kb PIC-scale, lambda_2 5-8 kb domain-scale). So ~1-2
  nucleosomes upstream, ~25 downstream. Strongly asymmetric — the lattice should
  not be symmetric about the TSS.
- Typical sizes: N_core ~ 5-10, N_flank ~ 25 downstream / ~2 upstream, so
  n ~ 35-40 nucleosomes. Small enough for exact SSA.

Note N_core ~ 5-10 lands squarely in the regime where the stochastic result from
`stochastic_amplifier_2026-10-22.md` matters (Nsites 4-12 changed the locked
fraction from 0 to 100%). The site count is not a nuisance parameter.

## 3. The options

### Option A — two-box with data-derived site counts (minimal)
Keep the ODE and the two compartments; replace the scalar occupancies with
binomial occupancy over N_core and N_flank sites taken from the ChIP geometry.
- **Gains:** `Nsites` stops being a free parameter; the noise axis that we just
  showed controls the lock boundary becomes data-derived.
- **Does not give:** boundary position, core width as an output, or any use of
  the KDM6B/Gli1 ChIP tracks.
- **Cost:** days. No architecture change.

### Option B — 1D nucleosome lattice, one-way coupled (RECOMMENDED)
Explicit lattice, each nucleosome carrying a K27 state from {u, ac, me1, me2,
me3} (Newar's 5-state, so p300 and PRC2 compete for the same residue), plus a
K4me3 flag in the core.
- **PRC2:** nucleation at the CpG island, propagation LOCAL (nearest-neighbour
  dominant) with EED read-write boost `a0 + a_rw*f(me3 in neighbourhood)`.
  Hodges & Crabtree is the constraint: bounded domains need local propagation and
  kappa = k_prop/k_turn in [1.0, 1.5], with both rates ~0.1-0.2 h^-1.
- **Antagonists, all LOCAL and non-spreading** (the Newar/Jost result — only
  PRC2 carries the nonlinearity):
  - p300 -> K27ac in the core, rate proportional to promoter activity;
  - **KDM6B/JMJD3 with rate `k_jmjd3_gli*Gli1*psi_KDM6B(x)`** — the existing
    kinetics, now with a measured spatial profile;
  - transcription-coupled H3.3 exchange over the body, proportional to nascent
    output (the existing `k_h33` term, spatially resolved).
- **Replication:** at S-phase each nucleosome independently retained with p=1/2.
  This replaces `Mk -> Mk/2` AND supplies the partition noise currently bolted on
  in `sim_inheritance_partition.py` — it comes free and correct.
- **Output to the cell cycle:** one scalar, proximal PRC2 occupancy computed from
  the core's me3, feeding the existing `prox_amp`. The interface does not change.
- **Gains:** boundary position and width become OUTPUTS; the depleted core
  emerges instead of being assumed; stochasticity and partition noise are native;
  the KDM6B/Gli1 tracks are usable; per-gene geometry makes it portable beyond
  Ccnd1/Cdk6.
- **Cost:** weeks. Needs the one-way coupling harness (0a), not the port.

### Option C — full 5-state + K36/K4 + 3D contacts (research-grade)
Adds the K36me2 gene-body domain with Alabert's measured state-dependent
antagonism (11.9-fold K27-on-K27me1/K36me0; 4.5-fold K36-on-K27me0/K36me2), and
a 1/|i-j| contact kernel for long-range PRC2.
- **Why probably not, yet:** Anink-Groenen's result is that 3D connectivity
  *destroys* the local boundary and converts the system to a global bistable
  switch — so adding the contact kernel may remove the very boundary we are
  trying to model. And Jost needed no K36 term to fit TSS-proximal profiles.
  Add K36 only if Option B underfits the downstream flank.

## 4. Parameterisation — two separable identification steps

This is the part that makes the redesign tractable, and it does not require
kinetics.

**Step 1: profile shape fixes the spatial parameters (static ChIP alone).**
- flank decay length -> kappa = k_prop/k_turn (Hodges & Crabtree);
- core width -> antagonist strength (p300 + KDM6B) relative to PRC2;
- up/down asymmetry -> separate upstream/downstream kernels, already fitted
  (248 bp vs 4.7 kb).
No expression data needed. This can be done against the existing
`chipseq_norm` TSS matrices today.

**Step 2: the cross-tissue H3K27me3 -> TPM relationship fixes the gate.**
The `transcriptional-model` MODEL.md states that gates A (Q*exp(-L)), B
(Q*(1-m3bar)) and C (competitive sigmoid) are **not decidable on bulk data** —
A and B are numerically indistinguishable and ranking tests are vacuous. That is
true for a cross-GENE scatter at one condition. It is NOT obviously true for the
same gene across many cell lines and tissues: that gives a **within-gene
dose-response**, many (me3, TPM) pairs per promoter, which is a far better
conditioned identification problem. Worth checking formally before assuming it
is still unidentifiable.

**The test the two-compartment structure actually needs.** With per-gene geometry
in hand, ask whether TPM is predicted by *core* me3 alone or by *total domain*
me3. That is the direct empirical test of "repression is local (proximal
promoter)" — the central claim the module rests on — and it cannot be done
without nucleosome-resolved geometry. If total domain me3 wins, the
two-compartment split is wrong and the whole redesign should stop there.

## 5. Open question / what is needed

I could not locate the cross-cell-line/tissue H3K27me3-vs-TPM pull. Found:
- `Medulloblastoma_projects/H3K27me3_fast_change/Scott_Bulk_Quant/` — 55-sample
  mouse bulk RNA-seq (mm39, STAR + featureCounts, DESeq2 tables), with a session
  note about a server agent updating a private repo;
- `R_projects/chipseq_norm` — 21,600 TSS +/-250 bp windows, H3K4me3 / H3K27ac /
  H3K27me3 across GNP P7/P12 and MB.

Neither is the many-cell-line/tissue resource described. **Where does it live,
and does it carry ChIP profiles (needed for Step 1 geometry) or only
summarised H3K27me3 levels (Step 2 only)?** The answer decides whether Step 1
can be done per-gene across tissues or only on the GNP/MB data.
