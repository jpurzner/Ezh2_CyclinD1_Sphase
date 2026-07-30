# Dynamic CDKI module — design spec for review (round 2)

**Status:** built, opt-in behind `build_model_v44(with_cdki_species=True)`; **default OFF** (default model unchanged at 30/32). Round-2 incorporates the 2026-07-30 review (KPC magnitude/pools, saturating CdP21, §3.1↔§3.3 fix, PCNA-arm move). Rates are **inherited** (Heldt) / **data-anchored** (RNA-seq tones) / **placeholder** (to fit). **Not yet calibrated** — this is the structure to finalize on the Stage-3 run.

Purpose: replace the lumped/mislabeled CDKI representation (static `p16`,`p18` params + one `P21` species standing in for all CIP/KIP + a mis-attributed p21 route) with the **individually expressed** CDKIs as dynamic species.

---

## 1. Scope

RNA-seq (P7 wt GNP → SHH-MB; anchored so Cdkn2c `p18` GNP = 0.464, factor ≈ 4.04e-4):

| gene | GNP | MB | MB/GNP | tone GNP→MB | modeled |
|---|---|---|---|---|---|
| **Cdkn2c (p18)** INK4 | 1149 | 4288 | 3.7× | 0.46→1.73 | **yes** (dominant INK4) |
| **Cdkn2d (p19)** INK4 | 890 | 1449 | 1.6× | 0.36→0.58 | **yes** |
| Cdkn2a (p16) | 3 | 417 | 139× | ~0→0.17 | no (silent GNP) |
| Cdkn2b (p15) | 25 | 154 | 6× | — | no (negligible) |
| **Cdkn1b (p27)** CIP/KIP | 9168 | 10893 | 1.2× | dominant (~85%) | **yes** (`P21`) |
| **Cdkn1a (p21)** CIP/KIP | 890 | 2511 | 2.8× | secondary, MB-up | **yes** (`p21a`) |
| Cdkn1c (p57) | 495 | 173 | 0.35× | Sox2+/quiescent CKI | **de-emphasized** — species kept but carried OFF (`kSyp57a=0`); minimally expressed in proliferating GNP/MB (JP 2026-07-30) |

---

## 2. Species (10 added)

`p18_prot`,`p19_prot` (INK4 proteins) · `p21a`,`p57a` (free CIP/KIP) · `Cep21a`,`Cap21a`,`Cep57a`,`Cap57a` (p21/p57 · CyclinE/A-CDK2) · `CdP21` (p27 buffered on CyclinD1-CDK4/6).
Re-assigned: `iPcna`,`iRc` (PCNA / inhibited-replication-complex) now carry **p21** (moved off p27). Inherited: `P21`(=p27), `CeP21`,`CaP21`, `Ce`,`Ca`.

---

## 3. Reaction network (exact rate laws)

### 3.1 INK4 arm + CDK4/6 → Rb drive
```
p18_synthesis:  => p18_prot;  kDe_p18*p18       (steady state p18_prot = p18; hook for Gli/EZH2 later)
p18_degradation: p18_prot => ; kDe_p18*p18_prot     [p19 identical]
```
Rb drive (both Rb-phos reactions) — **p16 dropped; free p27 NOT in the brake (acts on CDK2 only); buffered CyclinD `CdP21` in numerator + denominator (option A)**:
```
kPhRbCd*(Cd + w_Cd2*Cd2 + CdP21) / ( K_CdRb*(1 + p18_prot + p19_prot + w_p27*(p21a + p57a)) + (Cd + w_Cd2*Cd2 + CdP21) )
```
INK4 (p18/p19) is the CDK4/6 competitive brake. `w_p27*(p21a+p57a)` is a **phenomenological** minor term (p21/p57 have no explicit CyclinD complex; note p21 is a poor CDK4/6 activator). p27's CDK4/6 effect is stoichiometric via `CdP21`.

### 3.2 CIP/KIP → CyclinE/A-CDK2 sequestration
p27 keeps CeP21/CaP21. p21/p57 get their own complexes (same `kAsCyP21`,`kDsCyP21`); **p21 also carries the PCNA/Rc arm** (§3.4):
```
Assoc_CycE_Cdk2_p21a: Ce + p21a -> Cep21a; kAsCyP21*Ce*p21a - kDsCyP21*Cep21a     [+ CycA, + p57a]
Synthesis_of_p21a: => p21a; kSyp21a + kSyp21aP53*P53      (p57a: kSyp57a, no p53 arm)
```
Free Ce/Ca = active CDK2. Cyclin degradation releases the CKI (`Cep21a => p21a`) — the switch that liberates active CDK2 (verified for p27, mirrored here).

### 3.3 CdP21 — CyclinD1-p27 redistribution (option A, stoichiometric + saturating)
p27 binds CyclinD1-CDK4/6; **Cd is consumed** (buffer saturates at Cd abundance → collapses non-linearly as Cd falls). Buffered Cd is returned to the Rb-drive numerator (§3.1) so **total drive = Cd_total** (option A: the p27·CDK4 trimer is part of the drive). INK4 promotes release → p27 floods CDK2:
```
Buffer_p27_on_CyclinD: Cd + P21 -> CdP21; kSeqCd*Cd*P21 - kRelCd*(1 + w_ink4*(p18_prot + p19_prot))*CdP21
Deg_p27_in_CdP21:      CdP21 => Cd;       kDeP21*CdP21            (basal only; KPC can't reach it; Cd released)
```
**CdP21 is NOT reset at division** (CyclinD1 persists through mitosis → buffer full into G1, no mass sink). This is the **new CyclinD1-gated commitment + INK4-arrest lever** replacing the old CDK4/6-gated clearance.

### 3.4 Degradation routes (post-review)
**p27 (`P21`)** — basal + **KPC (constant, FREE pool only)** + SCF^Skp2. (KPC acts on cytoplasmic Ser10-P p27 → not on complexes/CdP21.)
```
free P21:            (kDeP21 + kDeKPC + kDeP21Cy*Skp2*(Ce+Ca)) · P21
CeP21, CaP21:        (kDeP21          + kDeP21Cy*Skp2*(Ce+Ca)) · [complex]
```
**p21 (`p21a`)** — free/CDK2-bound: basal + Skp2; **CRL4^Cdt2 only on PCNA/Rc-bound p21** (mechanistic):
```
free/Cep21a/Cap21a:  (kDep21a + kDeP21Cy*Skp2*(Ce+Ca)) · [pool]
iPcna, iRc:          (kDep21a + kDeP21Cy*Skp2*(Ce+Ca) + kDeP21aRc*Cdt2*aRc) · [pool]
PCNA/Rc binding:     aPcna + p21a -> iPcna ;  aRc + p21a -> iRc ;  iPcna => p21a (export)
```
**p57 (`p57a`)** — basal + SCF^Skp2 (2nd ligase FBL12, not modeled).

---

## 4. Reviewer corrections applied

1. **Cdt2 provenance + PCNA arm moved to p21.** CRL4^Cdt2 and the entire PCNA/Rc binding are p21's (PIP-degron); p27 has none. Both moved p27→p21a; Cdt2 now acts on PCNA-bound p21 (`iPcna`/`iRc`), retiring the free-p21a `aRc` proxy.
2. **KPC magnitude + pool.** `kDeKPC` 0.05→**0.006** (basal+KPC t½ ~1–2h, not 13 min) and restricted to the **free** pool only. This restored a usable p27 level (was ≈0).
3. **CdP21 stoichiometric + saturating.** `Cd + P21 ⇌ CdP21` (Cd consumed) with buffered Cd back in the drive numerator (option A). Fixes the non-saturating catalytic version → capacity-limited, switch-like exit.
4. **§3.1↔§3.3 double-hit fixed.** Free p27 removed from the Rb denominator (it can't brake CDK4/6 until bound). p27→CDK2 only; p27→CDK4/6 via CdP21 stoichiometry.
5. **No division reset for CdP21** (was a mass sink; CyclinD1 persists at mitosis).
6. **INK4 turnover split** — p18 stable (`kDe_p18`=0.003), p19^INK4d short (`kDe_p19`=0.03).
7. **Confirmed intact:** Skp2 destruction of complexed CKI releases active Ce/Ca (`CeP21=>Ce`, `Cep21a=>Ce`).

---

## 5. Parameters — status

| param | value | status |
|---|---|---|
| `p18`,`p19` (targets), `p21a`/`p57a`/`p27` tones | data | **data-anchored**; per-condition in `validate` (Stage 3) |
| `kDe_p18`,`kDe_p19` | 0.003 / 0.03 | estimate |
| `kSyp21a`,`kSyp57a`,`kSyp21aP53`,`kDep21a`,`kDep57a` | 2e-5 / 1e-5 / 1e-4 / 0.02 / 0.02 | **placeholder** |
| `kDeKPC` | 0.006 | placeholder (t½ ~1–2h) |
| `kSeqCd`,`kRelCd`,`w_ink4` | 0.05 / 0.01 / 1.0 | **placeholder** (CdP21 buffer) |
| `kAsCyP21`,`kDsCyP21`,`kAsPcP21`,`kDsPcP21`,`kExPc`,`kDeP21`,`kDeP21Cy`,`kDeP21aRc`,`Cdt2`,`kDeCe`,`kDeCa`,`kExPc`… | — | **inherited (Heldt)** |

Cycling confirmed (structure sound): with `kSeqCd≈0.5, kDeKPC≈0.02, kSyP21≈0.001`, GNP does 18 divisions; free p27→0. Calibration lever = the free-p27 balance (buffering + KPC clearance + synthesis).

---

## 6. Open items for review (deferred to calibration unless noted)

- **Free-p27 Skp2 over-degradation** — T187 phosphorylation needs p27 *in complex*; the model applies `kDeP21Cy·Skp2·(Ce+Ca)` to the free pool too (mass-action proxy). Known simplification, known direction. (Same for p21a/p57a.)
- **INK4 release is catalytic, not competitive** — no explicit INK4·CDK4/6 complex; `kRelCd·(1+w_ink4·INK4)` accelerates release without INK4 occupying a site (unbounded in INK4). Reviewer suggests an explicit complex (also makes the §3.1 minor term stoichiometric). **Not yet added.**
- **p57 has no regulatory hook** — `kSyp57a` is a bare constant, but its 0.35× GNP→MB drop is the Cdkn1c-H3K27me3 story; the synthesis term is the natural EZH2/mark hook. **Candidate for the mark coupling.**
- **p19 split is currently cosmetic** — with constant synthesis, `kDe_p19` changes only relaxation time, not steady state. Needs cycle-coupled synthesis if the S-phase peak matters.
- **p21a/p57a degeneracy** — they share binding constants + the Skp2 term; against aggregate CDK2 activity `kSyp21a`/`kSyp57a` may be near-unidentifiable (check on the sloppy direction before spending a fit).
- **MB 23.6h cycle time** is an unvalidated output — confirm against data (no formal MB cycle-length target yet).
- **CdP21 buffers CyclinD1 only** (not Cd2); free `Cd` (not total) is still used in `mu_eff`/EZH2-dose (minor while buffering is modest) — flag for consistency.

---

## 7. Stage 3 calibration — RESULT

250-candidate search over `{kSeqCd, kRelCd, kDeKPC, kSyP21, w_ink4, k_mu_cki}` (CDKI on, INK4 data tones p18 1.73 / p19 0.58 MB), hard-gated on GNP-SHH arrest + GNP+SHH cycles, cKO post-filtered.

**SUCCEEDED — 30/32 with all four constraints:** best = `kSeqCd 0.149, kRelCd 0.026, kDeKPC 0.0086, kSyP21 0.00158, w_ink4 2.86, k_mu_cki 0.259`. Verified: **30/32**, GNP-SHH arrest 0 div, **cKO@0 = 0.00**, GNP 16.9h / MB 23.0h (dichotomy), CyclinD1/D2 folds + EZH2 MB/GNP + MB 2N% + CDK4/6i±EZH2i all pass; fails = the structural HU pair only. The **CdP21 buffer is the Shh-gated commitment lever** exactly as designed (SHH 0.5: free p27 0.006 → commit; SHH 0: 0.29 → arrest).

**Both caveats RESOLVED (round 3):**
1. **p21a negativity — FIXED.** Root cause: the mechanistic PCNA-arm move — `iPcna`/`iRc` are replication-coupled *shared* species that can't host a near-zero p21 pool (p21a → −0.4, corrupting the sim). **Reverted** the move: p21 keeps its Cdt2/S-phase clearance via the `aRc` proxy, `iPcna`/`iRc` stay on p27 (Heldt) with p27's Cdt2 removed. p21a is now positive everywhere. A fully mechanistic p21-PCNA arm needs dedicated (non-shared) p21-PCNA species — **deferred**.
2. **p21/p57 cannot be engaged at data-transcript levels — a genuine finding.** Added CDK2-potency multipliers `w_cdk2_p21/p57` and removed p21/p57 from the CDK4/6 brake (poor CDK4/6 modulators), but across ~560 sims p21/p57 at ~0.49/0.24 **arrest the cell** even with p27 decomposed to near-zero and at low potency. This is the reviewer's degeneracy + the biological point that **mRNA abundance overstates the active nuclear CDK-inhibitor pool** for these heavily post-translationally-regulated proteins. **Resolution: p21/p57 remain expressed-but-low-activity species; the functional CDK regulation is p27 + INK4 + CdP21.**

**Clean structure RE-CALIBRATED: 30/32** (PCNA-move reverted, p21/p57 inert, p21a positive): `kSeqCd 2.13, kRelCd 0.0165, kDeKPC 0.0202, kSyP21 0.00116, w_ink4 3.62, k_mu_cki 0.362`. Verified: **30/32** (fails only the structural HU pair MB HU S + MB HU G2), cKO@0 = 0.00, GNP 17.6h / MB 26.1h dichotomy, Shh-dependence via the CdP21 buffer.

**Status: calibrated exemplar of the STRUCTURE (NOT baked).** Default remains 30/32 two-cyclin. Params in `scratchpad/recal_cdki_clean_best.json`; validate via env `CDKI_SPECIES=1 P18_MB=1.73 P19_MB=0.58` (NB CDKI-on validate is slow/stiff — use direct cycle tests for quick checks).

*Source: `src/build_model_v44_heldt.py`, `with_cdki_species` block (before `_ts_bake`).*
