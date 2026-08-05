# Dynamic CDKI module — design spec v2 (post-review)

> **⚠ STATUS UPDATE 2026-08-05 — parts of this spec are STALE. See `docs/cdki_optionB_progress_2026-08-05.md`.**
> - **Phase A (§5) is DONE and baked**, not "NOT YET APPLIED": the shipped model (`with_cdki_species=True`, now the
>   DEFAULT, 30/32) already has KPC on the free pool, stoichiometric+saturating CdP21, P21 out of the Rb denominator,
>   mass conservation, PIP-degron default, and release products. All Phase A stop-checks verified (`phaseA_verify.py`).
>   The P21≈0 bug is fixed (p27 pool nonzero, but ~99% buffered and INERT under option A).
> - **Phase B1 (option B, p27-inhibitory) is IMPLEMENTED** behind `with_p27_optionB` (default OFF; harness `P27_OPTIONB`).
>   Reaches **31/32** (beats option A's 30/32), Fan-Meyer-correct, and accommodates a biological p27 t½ (~1.5 h).
> - **Characterization result:** option B is *necessary-not-sufficient* — p27 stays inert far from the CyclinD1/p27
>   threshold, so it does NOT by itself tune GNP G1 or source MB G0 (both need near-threshold / population / DREAM-B5).
>   A partial-Shh dose-response target must be added before promoting option B (it shifts Shh dose-sensitivity).

**Supersedes:** the 2026-07-30 design spec.
**Status:** built, opt-in behind `build_model_v44(with_cdki_species=True)`; **default OFF** (default model unchanged at 30/32).
**File:** `src/build_model_v44_heldt.py`, `with_cdki_species` block (before `_ts_bake`).

**What changed in v2:**
1. A second structural review identified **four defects that were not in the first review** — a `kDeKPC` magnitude/pool error that explains the P21≈0 symptom, a non-saturating CdP21 buffer, a sign contradiction between §3.1 and §3.3, and an incomplete PCNA un-map. These are §5.
2. A **major literature correction**: the cyclin D1/p27 ratio is a *published, validated* decision variable for exit to quiescence at mitotic exit (Fan & Meyer 2021). The CdP21 module is therefore the modelled form of an established mechanism, not a speculative extension. This raises the priority of getting it right. See §6.1.
3. New §7 on **which CKI can and cannot source a transient G0**, and new §2.1 on **when RNA is and is not an acceptable proxy** — the latter changes how three of the seven CKI tones should be set.
4. **§7.3 rewritten from the primary papers** (Uziel 2005 *Genes Dev*; Uziel 2006 *Cell Cycle*), read in full rather than from abstracts. This **withdraws two arguments** from the v2 first draft and adds a timer-vs-rate distinction that constitutes a falsifiable test of the INK4 arm, plus a new constraint on the MB cycle time (§9).

---

## 1. Scope — which CDKIs, and why

RNA-seq (normalized counts, P7 wt GNP → SHH-MB; anchored so Cdkn2c `p18` GNP = 0.464 model units, factor ≈ 4.04e-4):

| gene | GNP | MB | MB/GNP | modeled? | RNA a valid proxy for protein? |
|---|---|---|---|---|---|
| **Cdkn2c (p18)** INK4 | 1149 | 4288 | 3.7× | **yes** | **No — see §2.1(b)** |
| **Cdkn2d (p19)** INK4 | 890 | 1449 | 1.6× | **yes** | partly; p19 is short-lived and S-peaking |
| Cdkn2a (p16) INK4 | 3 | 417 | 139× | no — silent in GNP | n/a |
| Cdkn2b (p15) INK4 | 25 | 154 | 6× | no — negligible | n/a |
| **Cdkn1b (p27)** CIP/KIP | 9168 | 10893 | 1.2× | **yes** (`P21`) | **No — see §2.1(a)** |
| **Cdkn1a (p21)** CIP/KIP | 890 | 2511 | 2.8× | **yes** (`p21a`) | yes (p53/transcriptionally driven) |
| **Cdkn1c (p57)** CIP/KIP | 495 | 173 | 0.35× | **yes** (`p57a`) | yes, but **compositionally confounded — §2.1(c)** |

### 2.1 RNA→protein validity, per gene

This section is new and it changes three parameter decisions.

**(a) p27 — RNA is uninformative.** *CDKN1B* is constitutively transcribed. Polyak's original characterisation found a single 2.5-kb transcript at similar levels in all tissues examined, unchanged between exponentially proliferating and contact-inhibited cells, and unchanged on release from contact inhibition. Hengst & Reed (PMID 8596954) showed directly that **p27 protein varies through the cycle and on growth arrest while p27 mRNA abundance remains unchanged**, the post-transcriptional change being achieved by translational control (polysome loading, higher in G0) plus altered half-life.

Consequences:
- The 1.2× MB/GNP ratio is the signature of a non-transcriptionally-regulated gene, not evidence that p27 is elevated in MB.
- 9168 counts in P7 GNP is what almost any tissue or cell line returns; it is not a distinguishing feature of the system.
- **Internal confirmation:** the outer EGL has *Cdkn1b* at 9168 normalized counts and no detectable p27 by IHC. The RNA–protein disconnect is already demonstrated in our own tissue.
- **External confirmation from this lineage:** Uziel 2005 (**PMID 16260494**) states directly that IHC "is mandatory for studies of p27^Kip1 expression, which is predominantly regulated post-transcriptionally."
- **Do not set `kSyP21` (p27 synthesis) from the transcript.** Fit it to protein/IHC behaviour.
- **Missing mechanism:** p27 synthesis is itself a regulated variable (translation rate rises in G0). The model has it as a constant. This is a p27-specific positive feedback — arrest → higher p27 translation → deeper arrest — that p21 does not have. If p27 is to be capable of generating bistability on its own, this is where it lives. **Not implemented; see §8.**

**(b) p18 — RNA is retained while protein is lost; the disconnect is proven in human MB.** In the granule lineage, *Ink4c* mRNA is confined to the EGL and extinguished in post-mitotic neurons, but p18^Ink4c is a stable protein detectable through the migratory phase (Uziel 2005/2006, **PMIDs 16260494 / 16479172**).

The decisive evidence is human. In 30 human MBs analysed by 100K SNP array, there were **no deletions at the *INK4C* (*CDKN2C*) locus and *INK4C* RNA was still expressed**, yet p18^INK4C protein was undetectable by IHC in 19% of a 73-tumour cohort — the authors conclude that alterations in p18^INK4C protein expression "reflect post-transcriptional regulatory mechanisms," and close by asking how p18^INK4C protein levels are controlled post-transcriptionally (Uziel 2005, **PMID 16260494**).

**Direction of the error:** a transcript-derived p18 tone will *over*estimate protein in MB. Our 3.7× rise is a transcript observation and is **not** contradicted by the human data, which also retain *INK4C* RNA. See §7.3 for the full review.

**(c) p57 — the bulk drop is compositionally confounded.** A 0.35× bulk decrease is exactly what is expected if p57 is confined to a rare quiescent subpopulation while the Atoh1+ bulk expands around it. A 2025 SHH-MB preprint (bioRxiv 2025.04.27.650711, *not peer-reviewed*) reports nuclear p57 enriched in Sox2+ over Sox2− cells (p = 0.00028, n = 667 nuclei) and p57 protein ranking Nestin+ > Sox2+ > Atoh1+, with a cytoplasmic shift in cycling Atoh1+ cells. **"p57 is lost in MB" and "p57 is the CKI of the quiescent MB compartment" are both consistent with the same bulk number.** Bulk RNA-seq cannot distinguish them, and it is the second reading that bears on transient G0.

---

## 2. Species added (9)

| species | meaning |
|---|---|
| `p18_prot`, `p19_prot` | INK4 proteins (Cdkn2c, Cdkn2d) |
| `p21a` | p21 (Cdkn1a) free functional pool |
| `p57a` | p57 (Cdkn1c) free functional pool |
| `Cep21a`, `Cap21a` | p21·CyclinE-CDK2, p21·CyclinA-CDK2 |
| `Cep57a`, `Cap57a` | p57·CyclinE-CDK2, p57·CyclinA-CDK2 |
| `CdP21` | **p27 buffered on CyclinD-CDK4/6** |

Inherited: `P21` (= p27), `CeP21`, `CaP21`, `Ce`, `Ca`, `iPcna`, `iRc` (the last two repointed to p21a behind `with_p21_pip_degron`).

> **Provenance warning.** `P21` is the Heldt p21 species relabelled p27. Every inherited kinetic constant on it (`kDeP21`, `kDeP21Cy`, `kAsCyP21`, `kDsCyP21`, `kDeP21aRc`) is a **p21** parameter. The CRL4^Cdt2 arm was one visible consequence and has been removed; the rest of the parameter set has not been re-derived for p27. p21 and p27 differ in half-life, degron architecture, and sign at CDK4/6 (§7.1). Treat all inherited `P21` rates as placeholders pending calibration.

---

## 3. Reaction network (current implementation)

### 3.1 INK4 arm → CDK4/6 brake
```
p18_synthesis:   => p18_prot;  Cell*kDe_p18*p18
p18_degradation: p18_prot => ; Cell*kDe_p18*p18_prot     # steady state p18_prot = p18
p19_synthesis:   => p19_prot;  Cell*kDe_p19*p19
p19_degradation: p19_prot => ; Cell*kDe_p19*p19_prot
```
Rb drive (both Rb-phos reactions), p16 dropped:
```
kPhRbCd*(Cd + w_Cd2*Cd2) / ( K_CdRb*(1 + p18_prot + p19_prot + w_p27*(P21 + p21a + p57a)) + (Cd + w_Cd2*Cd2) )
```
**⚠ This denominator contradicts §3.3 — see §5.3.**

### 3.2 CIP/KIP → CyclinE/A-CDK2 sequestration
```
Assoc_CycE_Cdk2_p21a: Ce + p21a -> Cep21a; kAsCyP21*Ce*p21a - kDsCyP21*Cep21a
Assoc_CycA_Cdk2_p21a: Ca + p21a -> Cap21a; kAsCyP21*Ca*p21a - kDsCyP21*Cap21a
Synthesis_of_p21a:    => p21a;  kSyp21a + kSyp21aP53*P53
Synthesis_of_p57a:    => p57a;  kSyp57a
Deg_CycE_in_Cep21a: Cep21a => p21a; (kDeCe + kDeCeCa*Ca)*Cep21a
```
**⚠ Release products for the four new complexes not verified — see §5.7.**

### 3.3 CdP21 — CyclinD-p27 buffer (option A as implemented)
```
Sequester_p27_on_CyclinD: P21   => CdP21; kSeqCd*(Cd + w_Cd2*Cd2)*P21
Release_p27_from_CyclinD: CdP21 => P21;   kRelCd*(1 + w_ink4*(p18_prot + p19_prot))*CdP21
Deg_p27_in_CdP21:         CdP21 => ;      (kDeP21 + kDeKPC)*CdP21
```
**⚠ Three defects here — §5.1, §5.2, §5.4.**

### 3.4 Degradation routes
**p27 (`P21`)** — basal + KPC + SCF^Skp2; Cdt2 arm removed:
```
(kDeP21 + kDeKPC + kDeP21Cy*Skp2*(Ce+Ca)) · [P21 | CeP21 | CaP21 | iPcna | iRc]
```
**p21 (`p21a`)** — basal + SCF^Skp2 + CRL4^Cdt2:
```
(kDep21a + kDeP21Cy*Skp2*(Ce+Ca) + kDeP21aRc*Cdt2*aRc) · [p21a | Cep21a | Cap21a]
```
**p57 (`p57a`)** — basal + SCF^Skp2:
```
(kDep57a + kDeP21Cy*Skp2*(Ce+Ca)) · [p57a | Cep57a | Cap57a]
```

---

## 4. Corrections from review round 1 (2026-07-30) — applied

1. **Cdt2 provenance swap.** p27 has no PIP degron and is not a CRL4^Cdt2 substrate (Havens & Walter 2011; Abbas & Dutta 2009). Removed from p27, attached to p21a. ✅
2. **KPC ungated.** Replaced `kDeP21Cd·[CDK4/6 activity]` with a constant `kDeKPC`. ✅ *(But see §5.1 — the magnitude and pool assignment are wrong.)*
3. **CdP21 redistribution takes over commitment.** ✅ *(But see §5.2, §5.3.)*
4. **INK4 turnover split.** p18 stable (`kDe_p18=0.003`) vs p19 short-lived (`kDe_p19=0.03`). ✅ *(But see §5.6 — cosmetic as implemented.)*
5. **Confirmed correct, no change:** Skp2 destruction of complexed p27 releases active Ce/Ca (`CeP21 => Ce`). ✅

---

## 5. Corrections from review round 2 — **NOT YET APPLIED**

Ordered by expected impact. §5.1 alone probably explains the P21≈0 / CdP21≈0 symptom.

### 5.1 `kDeKPC` magnitude and pool assignment — **the P21≈0 bug**

Units are per-minute (`kDeP21`=0.0025 → t½ = 277 min = 4.6 h ✓). So `kDeP21 + kDeKPC` = **0.0525/min → t½ = 13 minutes**. That is ~20× basal and roughly 5–10× too fast; p27 t½ is ~1–3 h in cycling cells and longer in G0. No synthesis rate can accumulate p27 against that.

- **Fix (magnitude):** `kDeKPC` ≈ 0.004–0.010, putting basal+KPC at t½ ≈ 1–2 h.
- **Fix (pools):** KPC acts on Ser10-phosphorylated p27 **after CRM1 export** — its substrate is cytoplasmic and by definition not in nuclear cyclin–CDK complexes. Apply `kDeKPC` to **free `P21` only**, not to `CeP21`, `CaP21`, `iPcna`, `iRc`, or `CdP21`.
- **Consequence for §3.3:** `Deg_p27_in_CdP21` currently gives the reservoir the same 13-min half-life as the free pool. The buffer must carry **basal only**. This is the differential-stability requirement — p27 in cyclin E/A–CDK2 is the Skp2 substrate (T187 written by the kinase it is bound to); p27 in cyclin D–CDK4/6 is the protected reservoir. If they share a rate there is nothing to redistribute.

### 5.2 CdP21 does not saturate

`kSeqCd*(Cd + w_Cd2*Cd2)*P21` with Cd catalytic-but-not-consumed gives a steady state of
`CdP21/P21 = kSeqCd·(Cd + w_Cd2·Cd2) / (kRelCd·(1 + w_ink4·INK4))` — a **fixed ratio, linear in Cd, with no capacity limit**. Buffered p27 can exceed the number of cyclin D–CDK4/6 complexes without bound.

The buffer's defining property is that it is capacity-limited: its size is set by cyclin D abundance, so it collapses non-linearly as cyclin D falls. That collapse is where a cooperative-looking exit comes from for free.

**Fix:** make it a real bimolecular reaction (`Cd + P21 ⇌ CdP21`, Cd consumed) and put `CdP21` back into the §3.1 numerator:
```
kPhRbCd*(Cd + w_Cd2*Cd2 + CdP21) / ( ... )
```
Total drive = free + buffered = Cd_total, so the Rb drive is untouched (that *is* option A), but mass is conserved on both sides and the buffer saturates. This is also the mechanistically honest version — the trimer is a real species with real activity (§7.1).

### 5.3 §3.1 and §3.3 contradict each other

The Rb denominator contains `w_p27*(P21 + p21a + p57a)` — **free** CIP/KIP braking cyclin D. §3.3 asserts **buffered** p27 does not brake cyclin D. That is backwards: free p27 cannot affect CDK4/6 until it binds.

Dynamically this creates a hidden double-hit. As p27 flows into CdP21, free `P21` falls → the §3.1 brake *releases* → sequestration raises CDK4/6 activity. Then INK4-driven release raises free `P21` → brakes CDK4/6 again *on top of* the direct INK4 term. This is exactly the double-counting option A was chosen to avoid.

**Fix:** with an explicit complex present, drop `P21` from the denominator entirely and let the CDK4/6 effect come from stoichiometry. If the phenomenological term is retained for p21a/p57a (which have no explicit cyclin D complex), at minimum remove `P21` from it. Note also that one shared `w_p27` treats all three CIP/KIPs as equipotent at CDK4/6, which is wrong (§7.1).

### 5.4 Division reset `CdP21 = 0` is a mass sink

Where does that p27 go? If deleted rather than returned to `P21`, mass is lost every division — and the reservoir is emptied at the start of every G1, precisely when p27 is highest and the buffer should be fullest. Cyclin D is not degraded at mitosis, so there is no biological basis for the reset. **Check whether `tP21` is conserved across division.**

### 5.5 Correction #1 was applied incompletely — the whole PCNA arm belongs to p21

The *degradation term* moved off p27, but p27 still carries `iPcna` and `iRc` in the base (non-flagged) path. Those are Heldt's **p21** mechanism. p27 has no PIP box, does not bind PCNA, and does not inhibit replication complexes. Functionally, p27 can currently inhibit origin firing via `iRc` — a p21 function attributed to the wrong protein.

The old §6 note ("Omitted for the minor CIP/KIPs: the PCNA/replication-complex binding arm (p27 dominates it)") is exactly inverted and should be deleted. The `with_p21_pip_degron` flag already does the right un-map; the correction is to **make it the default rather than an opt-in**, which also retires the `aRc`-proxy note.

### 5.6 The p19 fix is cosmetic as implemented

`p19_prot` relaxes to a constant target. With constant synthesis, `kDe_p19` = 0.03 and 0.003 give the **same steady state** — the faster rate only shortens relaxation time. If p19's S-phase peak matters, the *synthesis* term needs cell-cycle coupling. If it does not, the split buys nothing and can be reverted to reduce parameter count.

### 5.7 Verify release products for the four new complexes

§3.4 lists degradation applying to `[p21a | Cep21a | Cap21a]` without stating products. If `Cep21a =>` destroys the cyclin as well, the same switch that was verified for p27 (`CeP21 => Ce`) is broken for p21/p57.

### 5.8 INK4 release is catalytic, not competitive

`kRelCd*(1 + w_ink4*(p18_prot + p19_prot))*CdP21` accelerates release without INK4 being consumed or occupying a CDK4/6 site. The real mechanism is competitive: INK4 binds CDK4/6 monomers, reducing the sites available to hold p27. With no explicit INK4·CDK4/6 complex, release is unbounded in INK4 and the "conserved mass" claim holds for p27 only. Given INK4 is the intended MB lever, add the explicit complex; it also makes §3.1's brake stoichiometric rather than phenomenological.

### 5.9 Identifiability warning: p21a and p57a are near-degenerate

They share binding constants, share the Skp2 term, and differ only in basal rate, synthesis, and p21's Cdt2 arm. Against aggregate CDK2 activity they will sit on the same sloppy direction. Check this in the planned Sethna/Transtrum analysis **before** spending a fit on `kSyp21a` and `kSyp57a` separately.

---

## 6. Literature corrections to the module's rationale

### 6.1 **MAJOR — the cyclin D1/p27 ratio is a published decision variable**

Fan & Meyer, *Cell Reports* 2021 (**PMID 34320337**): the variable memory of local cell density experienced by a **mother** cell controls cyclin D1 and p27 levels in newborn **daughters**, directing them to proliferation or quiescence. High density rapidly suppresses ERK, lowering cyclin D1 and raising p27 in daughters. Cell density and mitogen signals compete by shifting the **cyclin D1/p27 ratio** across a **single sharp threshold**; small changes in the ratio are converted by an ultrasensitive response controlling CDK4/6 activity and Rb hyperphosphorylation. They further report **residual cell-to-cell variability in the ratio among cells that experienced the same local density**.

Extended by Meyer-lab follow-up (*Communications Biology* 2026, from bioRxiv 2024.10.11.617852): mitogen signalling, contact inhibition and YAP–TEAD all converge on the **nuclear cyclin D1/p27 ratio in early G1**, which dictates Rb phosphorylation and cycle entry.

**Implications for this module:**
- p27 **is** competent to carry a proliferation–quiescence decision at mitotic exit. The earlier working assumption that p27 lacks a stochastic cell-autonomous input was wrong — the input is the mother's ERK history, and there is intrinsic residual variance on top of it.
- The decision variable is a **ratio**, not a level. `CdP21` + cyclin D is the modelled form of that ratio. This makes §5.2 (stoichiometric buffer) and §5.3 (sign contradiction) higher-priority than previously stated: they are defects in the model's central mechanism, not in an accessory module.
- **Testable now, no bench work:** after applying §5.1–5.3, add cell-to-cell variance in inherited cyclin D1 at division and check whether a CDK2-low subpopulation emerges with `p21a` held at zero. If it does, the module reproduces Fan & Meyer without a p21/p53 arm — consistent with p21 being undetectable in MB IHC and OLIG2-repressed in the stem compartment.
- **Discriminating prediction:** the p21 route is ~98% sister-correlated (input inherited before division; Spencer 2013, **PMID 24075009**). A cyclin-D-partitioning mechanism should show *lower* sister correlation, since partitioning noise at cytokinesis is anti-correlated between sisters. Sister correlation alone separates the two routes in live imaging.

### 6.2 p57's second ligase is FBL12, not FBXW7

Corrected in v1; retained here for the record. p57 phosphorylated at its C-terminal threonine is degraded in late G1/S by SCF^FBL12 and SCF^Skp2. Skp2-null cells accumulate p57, implying FBL12 cannot compensate for Skp2 loss. The Skp2 degron is Thr310 (human), and *in vitro* ubiquitylation by recombinant SCF^Skp2 requires cyclin E–CDK2 — so the shared `kDeP21Cy·Skp2·(Ce+Ca)` term is structurally right for p57 (Kamura 2003).

### 6.3 p21 has a third route not implemented

APC/C–Cdc20 in early mitosis. Flagged for verification, still not implemented. Add only if mitotic p21 dynamics turn out to matter to the dichotomy.

---

## 7. What each CKI can and cannot contribute to a transient G0

This section is new. It is the biology that should drive Stage 3 priorities.

### 7.1 p27 — competent, but conditionally signed at CDK4/6

**Positive evidence.** Fan & Meyer 2021 (§6.1). Adenoviral p27 vs p21 in breast lines gives greater cytotoxicity and G1/S arrest at lower MOI — more arrest per molecule. In adipogenesis the polarity is the reverse of the Spencer framing: rising p27 with falling cyclin D1 lengthens G1 reversibly, while PPARG-induced p21 and p18 produce the *permanent* exit (PMID 41272295).

**Negative evidence, and its limits.** Spencer 2013 (PMID 24075009) measured p27 in the single-cell marker panel and it did not discriminate CDK2^inc from CDK2^low daughters. Overton 2014 (PMID 25267623) showed p21-null MCF10A lose the split with p27 still present. Both establish that **p27 is not the variable carrying the decision in those systems** — not that p27 cannot enforce a CDK2-low state. Neither is a matched-stoichiometry substitution test, and neither measured **free nuclear p27**, which is the functional pool.

**The mechanistic asymmetry that matters for the model.** p27 associates with CDK4 in two alternative conformations — closed/inactive or open/activating — switched by Y88/Y89 phosphorylation (Brk/Src/Abl family). Guiley et al. (PMID 31831640) showed tyrosine phosphorylation of p27 is required for CDK4 activation and that **the absence of the corresponding tyrosine makes p21 a poor activator**. p21 by contrast is a monotone stoichiometric inhibitor (Yang/Meyer 2017, **PMID 28869970**: stoichiometric inhibition of cyclin D1–CDK4 by p21 controls Rb/E2F ultrasensitively; Hill coefficient rising from 4 to 10 over a fivefold increase in p21).

- At **CDK4/6**, p27's sign is set by tyrosine-kinase tone, not by its own concentration.
- At **CDK2**, p27 is inhibitory regardless of Y-phosphorylation (it binds cyclin E and occludes the substrate site), though pY88-p27 becomes a CDK2-dependent ubiquitination target.

**Model consequence:** implement option A as the **pY88 limit** — the trimer is the *active* species, contributing to the Rb drive (§5.2). Release from CdP21 then does two things in one step: removes an activator from CDK4/6 *and* adds an inhibitor to CDK2. That is a sharper switch than release-is-neutral, and it is the mechanistically correct reading.

**Missing:** the translational positive feedback (§2.1a).

### 7.2 p57 — the best candidate for a *source* of bimodality

*CDKN1C* is **imprinted** — monoallelically expressed. A methylation event on the single active allele silences the gene completely and heritably in that lineage. That is a cell-autonomous, all-or-none, heritable binary switch, which is the input structure a stochastic transient G0 requires, and it sits naturally inside the H3K27me3/epigenetic machinery already in the model.

Predicted population structure: most cells silence p57 and cycle; a minority retain the active allele and park. That matches the rare quiescent compartment, generated with no p53, no replication stress, no p21. The Shh→Gli→Dnmt→*Cdkn1c* hypothesis therefore is not only an explanation for the bulk drop — it is a candidate **source of the bimodality itself**.

**Model consequence:** `kSyp57a` as a low uniform constant set from the bulk transcript models the *average* and is the wrong object. If p57 is the rare-cell CKI, it must be a **bimodal or heritable** variable. It is also the species most in need of a regulatory hook (currently a bare constant, while the INK4 synthesis terms have one reserved).

Background: p57 is the validated deep-dormancy CKI in HSCs and adult NSCs (Matsumoto 2011; Zou 2011; Furutachi 2013 — PMIDs unverified, see §10), with p27 compensating on p57 loss.

For contrast, the MB genetics put p27 on a different axis: *Ptc1*^+/−;*Kip1*^+/− and *Ptc1*^+/−;*Kip1*^−/− mice develop medulloblastoma rapidly and with high penetrance, all tumours lose the wild-type *Ptc1* allele, but wild-type p27 protein is invariably retained in *Ptc1*^+/−;*Kip1*^+/− tumours — i.e. p27 is **haploinsufficient**, not deleted (Ayrault et al. 2009, **PMID 19147535**). Tumours neutralise p27 post-translationally rather than losing it, which is consistent with the p27+ Atoh1+ staining.

### 7.3 p18^Ink4c — primary-literature review (Uziel 2005 *Genes Dev* **PMID 16260494**; Uziel 2006 *Cell Cycle* **PMID 16479172**)

> **This section supersedes the v2 first-draft §7.3, which was written from abstracts and got two things wrong.** Corrections are flagged inline as ⚠.

#### (i) Expression: p18 precedes p27 along the GNP trajectory

In situ hybridisation at P7 (Uziel 2005, Fig. 1A): *Ink4c* is expressed **within the EGL in a pattern broader than *Math1*(Atoh1)**, with neither gene expressed in the deeper EGL adjacent to the Purkinje layer. The authors' reading: *Ink4c* "is induced in **dividing cells** and transiently maintained as GNPs exit the cycle and differentiate." *Ink4c* mRNA is confined to the EGL and absent in post-mitotic granule neurons.

p27^Kip1 is the mirror image: up-regulated in the **premigratory zone in Ki67-negative cells** that have already exited, with BrdU/p27 double staining confirming expression **only in post-mitotic cells**, and persisting in mature IGL neurons throughout adult life.

Timing (qPCR, whole cerebellum and Percoll-purified GNPs): *Ink4c* and *Math1* are concordant at P7, but *Ink4c* **declines less precipitously than *Math1* and is maintained through P12**; in purified GNPs it is restricted to P7–P12 and extinguished thereafter, and is less robust in glia/Purkinje fractions.

**Order along the trajectory: p18 ON in cycling GNPs → cell-cycle exit → p27 ON.** The abstracts of both papers ("expressed as CGNPs exit the cell cycle") compress this and lose the part that matters. Note also that mouse p18^Ink4c immunostaining "has not been readily achieved with the antibodies presently available" — all mouse p18 expression data here are RNA.

**Model consequence — p18 is a timer, not an exit marker.** A CKI present *during* proliferation whose level sets *when* exit occurs is a division-linked accumulating variable, i.e. the same class of object as the H3K27me3 dilution clock already in the model. Roussel's own framing is that p18 "times" the exit of GNPs from the cycle. Note that the authors state plainly that "the factors that control the transient expression of *Ink4c* in the EGL remain unknown" — this is precisely the gap the reserved Gli/EZH2 hook on the INK4 synthesis term is aimed at.

#### (ii) Functional phenotype: extends the Shh-responsiveness *window*, not the rate

Cultured purified GNPs, ±Shh, 22-h BrdU pulses (Uziel 2005, Figs. 3–4):

| comparison | day 1, no Shh | day 3, no Shh |
|---|---|---|
| WT vs *Ink4c*−/−;*p53*−/− DKO, P7 & P10 | **indistinguishable** | most WT exited; significant DKO fraction still cycling |
| WT vs DKO, P12 | ~50% of WT had lost BrdU capacity | — |
| WT / *Ptc1*+/− / *Ptc1*+/−;*Ink4c*−/− (Fig. 4) | ~32% / ~37% / ~33% — **no difference** | ~2% / ~10% / ~19% |

Single mutants (*Ink4c*−/− or *p53*−/−) each gave half-maximal BrdU incorporation relative to DKO; the two contribute **additively and independently**. The authors' conclusion: loss of *Ink4c* or *p53* "independently extends the **temporal window of responsiveness** of GNPs to Shh."

**This is a hard test for the model.** p18 does not change instantaneous proliferation rate — the day-1 difference is zero and the day-3 difference is large. Any implementation in which p18 acts as an instantaneous CDK4/6 brake will produce a day-1 difference and fail. This favours the timer reading in (i) and argues against relying on §3.1's competitive-brake term alone.

**Second hard constraint:** DKO GNPs **remain Shh-dependent** — only a small fraction stayed in cycle after 3 d without mitogen. CKI loss must not confer mitogen independence in the model.

#### (iii) Genetics: haploinsufficiency, and CKIs are retained in tumours

- **With p53:** up to 25% of *Ink4c*−/−;*p53*−/− mice develop MB by 5 months; neither single mutant does. **Both** *Ink4c* alleles are required in this background.
- **Irradiation, 4 Gy at P5** (Table 1): *Ink4c*−/−;*p53*+/+ **0/20** · *Ink4c*+/−;*p53*+/− 2/10 · *Ink4c*+/+;*p53*−/− 13/19 · *Ink4c*+/−;*p53*−/− 15/17 · *Ink4c*−/−;*p53*−/− 18/24 · *Ink4c*−/−;*p53*−/FL;Nestin-Cre 4/4 · *Ink4c*+/+;*p53*−/FL;Nestin-Cre 2/5.
- **With Ptc1** (Fig. 5A): *Ptc1*+/−;*Ink4c*+/+ **2/27** · *Ptc1*+/−;*Ink4c*+/− **31/68** · *Ptc1*+/−;*Ink4c*−/− **10/27**. Incidence in +/− does not differ significantly from −/−.
- Percoll-purified Math1+ tumour cells: wild-type *Ptc1* undetectable in 4/4 MBs while *Ink4c* remained expressed. Note added in proof extends this to 30 tumours — 6/6, 11/14 and 10/10 across the three genotypes expressed no detectable *Ptc1*, while **all 20 tumours from *Ink4c*+/+ and +/− animals continued to express *Ink4c***. **Ink4c is haploinsufficient; Ptc1 is not.**
- No *p53* mutations in 6 sequenced tumours (3 per genotype), p53 IHC characteristically low → **this route does not require p53 inactivation**, mimicking human MB where *TP53* mutation is infrequent.
- Spectral karyotyping of *Ink4c*/*p53* DKO tumours: the only consistent anomaly was loss of one copy of chromosome 13, where *Ptc1* resides.

**Convergence with p27.** Ayrault 2009 (**PMID 19147535**) found the identical pattern for *Kip1*: haploinsufficient, wild-type p27 protein invariably retained in *Ptc1*+/−;*Kip1*+/− tumours. **Neither CKI is deleted in MB; both remain expressed.** → **Do not zero any CKI in the MB condition.** MB tolerates its CKIs rather than removing them, which is consistent with the p27+ Atoh1+ staining and points to post-translational neutralisation (§7.1) rather than loss.

#### (iv) A quiescence-relevant phenotype of *Ink4c* loss

Irradiated *Ink4c*-null mice with wild-type p53 developed **no** MB, but invariably showed ectopic "nests" in the cerebellar molecular layer: **Ki67-negative**, GABA(A) receptor α6-positive (differentiated granule neuron marker), ~3-fold more numerous than in irradiated WT and located more deeply. So *Ink4c* loss produces **arrested, mislocalised cells that persist**, not proliferating ones — a migration/persistence phenotype, not a proliferation phenotype.

#### (v) Human MB

- **Methylation-specific PCR, 23 tumours:** 1 fully methylated + 3 hemimethylated *INK4C* promoter (**4/23, 17%**). Absent from 9 normal cerebella and from 11 immortalised MB cell lines lacking the SHH-activation signature. *INK4A* unmethylated. The authors note 17% "eclipses that of *TP53* mutations."
- **IHC, separate 73-tumour cohort:** 17 (23%, grade 2) strong · 42 (58%, grade 1) intermediate · **14 (19%, grade 0) virtually no signal**. All tumours cyclin D1-positive.
- **100K SNP array, 30 samples: no deletions at *INK4C*, and *INK4C* RNA still expressed** → loss of p18 protein is post-transcriptional.
- Of 46 comprehensively analysed: 3 with *TP53* mutations and 7 with *PTC1*-pathway mutations all expressed p18^INK4C; 6/7 *PTC1*-mutant samples scored grade 1. n precludes significant correlation with *TP53* or *PTC1* status.

#### (vi) ⚠ Corrections to the v2 first draft

1. ⚠ **"Direction conflict with human data" — withdrawn.** Human MB *retains* **INK4C RNA** and loses the **protein**. There is no RNA-level conflict with our 3.7×. The real conflict is RNA-vs-protein, and it applies within human MB itself.
2. ⚠ **"Compositional confound (more differentiating cells)" — runs the wrong way.** *Ink4c* mRNA is a **cycling**-GNP gene confined to the EGL and absent from post-mitotic neurons. A larger differentiated fraction would *lower*, not raise, bulk *Cdkn2c*.
3. ⚠ **"N-Myc tension" — largely dissolves.** *Ink4c* appears in the Fig. 5B set of ~30 genes **commonly overexpressed across mouse MB models** relative to P5 cerebellum, alongside Gli, Gli2, Foxm1, Patched-2, Math1, cyclin D1, cyclin D2 and N-Myc. Our 3.7× **replicates** a published mouse finding. (N-Myc loss does de-repress *Ink4c* and *Kip1* in the anlage — Zindy 2006, **PMID 16864777** — but that repression is evidently overridden in established tumours.)
4. **Retained and strengthened:** the resolving experiment is still p18 IHC with Atoh1 and pRb — but the prior has flipped. On this literature, p18 protein should be found in **Atoh1+/cycling** cells, not in the Atoh1− islands. Note the mouse antibody limitation above; a human/MB-validated antibody (M168, Santa Cruz, as used in Uziel 2005) may be the practical route.

#### (vii) Model directives from these two papers

| finding | directive |
|---|---|
| p18 ON in cycling GNPs, before p27 | Do not treat p18 as a post-exit marker. Consider p18 as a division-linked accumulating timer coupled to the epigenetic clock. |
| Day-1 no difference, day-3 large difference | p18 must not act as an instantaneous rate brake. Validate against the two-timepoint BrdU shape, not a single proliferation index. |
| DKO GNPs remain Shh-dependent | CKI loss must not confer mitogen independence in any simulation. |
| *Ink4c* and *Kip1* both haploinsufficient, both retained/expressed in tumours | Never set a CKI tone to zero in the MB condition. |
| *Ink4c* loss → Ki67-negative ectopic nests | *Ink4c* loss produces persistent arrested cells, not extra cycling ones. |
| Human MB: RNA retained, protein lost in 19% | Transcript-derived p18 tone over-estimates protein. Fit p18 protein, don't anchor it. |
| "Factors controlling transient *Ink4c* expression remain unknown" | The reserved Gli/EZH2 hook on the INK4 synthesis term addresses a stated open question — worth framing as such. |
| MB S-phase fraction ≈ P5 GNP (§9) | Re-examine the MB-slower assumption. |

### 7.4 INK4s are amplifiers, not sources — and p18 does **not** mark G0

**Does p18 anticorrelate with Ki67 or pRb at single-cell level? No such literature exists, and the mechanism says it shouldn't.**

**(a) The one direct comparison across a proliferation→arrest transition.** In CD40L-stimulated splenic B cells undergoing plasma-cell differentiation, immunoblotting after saturating immunoprecipitation showed **sustained expression of p18 and p19^INK4d through cell-cycle activation (S phase) *and* the subsequent onset of G1 arrest**. In the same experiment p21^Cip1 was progressively and markedly reduced during cell-cycle entry and arrest, and p27 diminished on CD40L then was restored in arrested cells (Tourigny et al., *Immunity* 2002;17:179–189). When all four CKIs were compared across one exit, **p18 was the one that did not change** — the definition of a marker that cannot stratify.

**(b) p18 protein is non-periodic by design.** Thullberg, Bartek & Lukas (*Oncogene* 2000;19:2870–6, **PMID 10851091**) showed the periodic oscillation of p19^INK4d protein is set by ubiquitin/proteasome degradation tracking its mRNA, and that **within the INK4 family this regulatory mode is restricted to p19^INK4d**. See also Kumar et al. (*PNAS* 2018, **PMID 29531090**) on phosphorylation-induced unfolding of p19^INK4d through the human cycle. p18 is flat across the cycle.

**(c) In stem cells, p18 sits on the self-renewal axis, not the quiescence axis.** Yuan/Cheng (*Nat Cell Biol* 2004, **PMID 15122268**): in vivo self-renewing divisions of HSCs are increased in the absence of p18, and "different CKIs have highly distinct effects on the kinetics of stem cells, possibly because of their active position in the cell cycle." Yu/Cheng (*Blood* 2006;107:1200–6): p18 and p21 affect HSC exhaustion in **opposite** manners. *Nat Commun* 2015;6:7328: p18 is a more potent inhibitor of HSC self-renewal than p27, and single-cell analysis showed p18 deletion **favoured self-renewing division**. Meanwhile p57/p27 govern HSC quiescence (Matsumoto 2011; Zou 2011). Clean dissociation: **p57/p27 = quiescence; p18 = division fate.**

**(d) The sign is not guaranteed.** In mouse ES cells, p18 overexpression **accelerated** growth of ES cells and embryoid bodies, with increased CDK4 expression proposed to release CDK2 from p21/p27 inhibition (Li et al., *PLoS One* 2012;7:e45212). ⚠ **Caution for §3.3/§5.8:** the model assumes raising INK4 pushes CIP/KIP onto CDK2 and arrests. In at least one system, raising p18 did the opposite via compensatory CDK4. Check that the `w_ink4` sweep does not enter a regime where the real system would reverse sign.

**(e) A uniform INK4 rise arrests uniformly.** Palbociclib, the pharmacological equivalent, abolishes both CDK4/6 and CDK2 reporter signals after mitosis (Yang 2020, **PMID 32255427**). Pennycook & Barr (*Open Biol* 2021;11:210125) further showed palbociclib-mediated arrest can occur **in the absence of both p21 and p27** — CDK4/6-route arrest does not require the CIP/KIP arm at all. INK4s amplify variance in cyclin D; they cannot create a subpopulation from a uniform input.

**(f) Ki67 was the wrong denominator anyway.** Miller/Spencer (*Cell Rep* 2018;24:1105–1112.e5) found Ki67 accumulates only in S/G2/M and is degraded continuously in G1 **and** G0, so its level depends on how long a cell has been quiescent: spontaneous G0 cells did not have distinctively low Ki67, with significant overlap against G1, and forced quiescence shifted the distribution **unimodally** downward. **Phospho-Rb S807/811 was bimodal** under the same conditions. Score arrest off pRb, not Ki67 — which the current IHC already does.

**Net:** the §7.3 timer reading of p18 survives and is strengthened — flat, non-periodic, sustained through both proliferation and arrest is exactly an accumulating threshold variable. But p18 IHC will tell you *where p18 protein is*, not *which cells are arrested*.

### 7.4b ⚠ New problem: p19^INK4d is an E2F1 target

Carcagno et al. (*PLoS One* 2011;6:e21938, **PMID 21765927**): p19^INK4d is transcriptionally regulated by **E2F1** through two response elements in its promoter; mRNA and protein accumulate periodically "reminiscent of cyclins" — low in G1, highest in S, low again approaching the next G1. Blocking E2F1 binding at the p19 promoter **stimulated proliferation and increased the S-phase fraction**. Proposed architecture: once cyclin E–CDK2 takes over as the driving kinase, E2F1-induced p19 shuts down the CDK4/6 complexes to close G1 — a hand-off.

**Consequence for the model.** `p19_prot` currently relaxes to a static expression-set target with no cycle coupling, and the 1.6× MB rise is read as more brake. If p19 is E2F-driven, its elevation in MB is at least partly a **readout of proliferation**, not a cause of arrest — and the model is missing a real negative feedback loop (E2F → p19 ⊣ CDK4/6 ⊣ E2F). This is a larger issue than the cosmetic half-life split in §5.6. **Decision required before calibration:** either couple `p19_prot` synthesis to E2F and accept the feedback loop, or drop p19 as a separate species and fold it into p18. A static E2F target modelled as an independent brake will misattribute a proliferation readout to inhibition, and the error runs in the same direction as the MB tone.

⚠ **This is a confound *class*, not a one-off.** *Mir17hg* has the same problem (§7.5c) — four E2F sites in its promoter, and it tracks the proliferation profile across our timecourse. Any E2F-target gene whose MB elevation is being read as a regulatory input needs this check. Current members: **Cdkn2d (p19), Mir17hg.** Screen the rest of the matrix for E2F targets before treating any fold-change as causal.

### 7.5 p21 (Cdkn1a) in GNP and MB — the canonical route may be structurally unavailable

**(a) p21 protein is not a documented GNP species.** Roussel and colleagues state flatly that **"the only CDK inhibitors expressed in the EGL are p18^Ink4c and p27^Kip1"** (Uziel 2006, **PMID 16479172**, citing Miyazawa 2000 and Zindy 2003). Miyazawa (*J Neurosci* 2000;20:5756–63) found only p27 at significant levels in developing and adult cerebellum, and **no compensatory upregulation of p18, p21 or p57 in p27-null GNPs**. In human MB, the ~14-tumour IHC series found **p21 undetectable in every tumour**. Our bulk *Cdkn1a* is 890 (GNP) → 2511 (MB): transcript present, protein unreported.

**(b) Three independent mechanisms suppress p21 in this exact lineage.**

1. **miR-17~92 — the strongest and most lineage-specific.** Northcott et al. (*Cancer Res* 2009, **PMID 19351822**): the miR-17/92 polycistron is amplified in **6%** of pediatric MBs; its components are the **most highly up-regulated miRNAs in medulloblastoma** across 90 primary human tumours; expression is **highest in the SHH subgroup**; Shh treatment of CGNPs increases miR-17/92; and **N-myc — but not Gli1 — induces it**. Ectopic miR-17/92 synergised with Shh and **enabled CGNPs to proliferate in the absence of Shh**. Uziel et al. (*PNAS* 2009, **PMID 19196975**): of 26 miRNAs overexpressed in proliferating mouse GNPs and MBs relative to mature cerebellum, **9 belong to the miR-17~92 family**, regardless of founding genotype (*Ink4c*−/−;*Ptch1*+/− and *Ink4c*−/−;*p53*−/−). Independently, Dicer1 ablation in GNPs produced a **selective increase in Cdkn1a/p21**, with miR-17-5p highly expressed in developing cerebellum. The miR-17/20a/106b/93 seed family are canonical *CDKN1A* repressors (Ivanovska et al., *MCB* 2008;28:2167–74). Therapeutic proof-of-concept: 8-mer LNA-antimiRs against these seed families inhibit MB progression in two SHH-MB mouse models (*Cancer Res* 2013;73:7068).
   → **Shh → N-Myc → miR-17~92 ⊣ p21 is a post-transcriptional repression circuit operating in GNPs and elevated in SHH-MB.** This reconciles transcript-present with protein-absent, and predicts that a 2.8× *Cdkn1a* transcript rise in MB need not produce any protein rise.
2. **OLIG2 transcriptional repression** in the quiescent stem compartment (Ligon et al., *Neuron* 2007) — ChIP binding at the *p21* promoter, p21 mRNA up in Olig2-null neurospheres, and OLIG2 does **not** bind *p27* or *p57*. Lineage-restricted and p21-specific.
3. **p53 dampened upstream even when wild-type.** In the NeuroD2:SmoA1 model with wild-type *Trp53*, endogenous I2PP2A/SET suppresses p53 by promoting phospho-MDM2^S166 accumulation and faster p53 degradation; I2PP2A knockdown induced p21 **only** in p53-wild-type ONS76 and not in p53-mutant Daoy or UW228 (*Mol Cancer Res* 2019;17:186).

**(c) ⚠ RETRACTED: "in GNPs the p53 response to damage is apoptosis, not p21 arrest."** This was wrong for the SHH-MB system and is retracted on the primary data. Tamayo-Orrego et al. (*Cell Rep* 2016;14:2925–37, **PMID 26997276**) measured cleaved caspase-3 across all three stages in *Ptch1*^+/−^ mice: apoptosis rises from EGL to preneoplasia (≈1.5 → ≈4.6 per 10⁴ µm², p ≤ 0.05) and then **does not change** into advanced MB (≈4.9, n.s.). Apoptosis is not the barrier to progression. Their discussion states the inverse of my claim: **p53 preferentially regulates cell senescence over apoptosis in these tumours**, and human SHH-MB shows a paucity of apoptosis-pathway mutations. (They also note the *Ptch1*-dependence-receptor prediction failed — preneoplasia had *more* apoptosis, not less.)

→ **`kSyp21aP53*P53` stays. What it feeds is a senescence arm, not an apoptotic one.**

### 7.5b The senescence trajectory (Tamayo-Orrego 2016, **PMID 26997276**)

**Stage-resolved IHC, *Ptch1*^+/−^ mice, n ≥ 5 animals** (areal densities read off Fig. 1B/D/F/H/J):

| stage | p21⁺ /10⁴µm² | p16^Ink4a⁺ /10⁴µm² | Ki67⁺ /10³µm² | pH3⁺ /10⁴µm² | c-casp3⁺ /10⁴µm² |
|---|---|---|---|---|---|
| P7 *Ptch1*^+/−^ EGL | ~1.0 | ~0.3 | ~17 | ~6.1 | ~1.5 |
| P14 preneoplasia | ~10 | ~5.3 | ~8.2 | ~2.9 | ~4.6 |
| Advanced MB | ~1.7 | ~0.9 | ~14.2 | ~4.6 | ~4.9 |

p21 and p16: both stage comparisons p ≤ 0.01. Ki67: EGL vs preneo p ≤ 0.001, EGL vs Adv MB **n.s.**, preneo vs Adv MB p ≤ 0.001. Single-cell anticorrelation is shown for p16: **p16^Ink4a⁺ cells in preneoplastic lesions are negative for proliferation** (Fig. S1B).

**⚠ Quantitative caveat before fitting any of these.** These are areal densities, not fractions, and the Ki67 denominator differs (10³ vs 10⁴ µm²). Cross-scaling gives p21⁺ ≈ 12% and p16⁺ ≈ 6% of the Ki67 density in preneoplasia — a minority that cannot obviously account for the ~52% drop in average proliferation. Either cell density per area differs across stages, the markers underestimate the arrested fraction, or non-senescent cells are also slowed. **Fractions, not densities, are needed to use these as fit targets.**

**Where the heterogeneity is.** The heterogeneous p21/Ki67 labelling is in **preneoplasia**, not advanced MB. Advanced tumours return to near-EGL p21 levels, significantly below preneoplasia. Scattered p21⁺ cells in advanced MB sit at ~1.7× the EGL density — a residual population, not the dominant biology. This is the paper's thesis: advanced MB has *evaded* senescence.

**⚠ NEW HARD CONSTRAINT — proliferation is non-monotonic in Hh drive.** EGL high → preneoplasia low → advanced MB high again, on both Ki67 and pH3, with senescence markers as the exact mirror image. The dip occurs as Hh signalling *increases* (*Ptch1*^+/−^ → *Ptch1* LOH). **A model in which proliferation rate increases monotonically with Shh drive cannot reproduce this.** The escape is specific: removing the p53 or p16 arm restores proliferation at the same high Hh drive.
→ **Test:** sweep Hh drive with the p53/p16 arm intact and look for a proliferation dip; disable the arm and the dip should vanish.
→ **Mechanism they nominate:** Hh hyperactivation after *Ptch1* LOH raises **N-Myc and cyclin D1** as candidate mediators of the oncogenic stress. This supports modelling replication stress as a function of cyclin D1 rather than as a constant (cf. Mille et al., *Dev Cell* 2014, **PMID 25263791**: Shh induces DNA damage in GCPs, requiring Boc and cyclin D1).

**Preneoplasia is a third condition the model lacks.** The senescent state sits between the GNP and MB conditions and is transient at the population level — most lesions regress. GFP-control transplant lesions **frequently displayed regions of differentiated granule neurons populating the IGL**, while p53^R270C^ lesions showed no features of differentiation. The senescence-competent route ends in differentiation and clearance, **not reversible parking** — this is not transient G0.

**Genetic sequence and frequencies.**
- Advanced MB: *Ptch1* LOH in 95% (18/19); somatic *p53* mutations in **36% (7/19)**, mostly heterozygous missense in the DNA-binding domain, dominant-negative. All *p53*-mutant tumours also had *Ptch1* LOH.
- LCM on P14 preneoplasia: 10/14 (71%) had *Ptch1* LOH; **none had *p53* mutations**; p53 IHC showed no nuclear accumulation. → **Ptch1 LOH precedes p53 mutation.**
- *Olig1-Gnas* model: *p53* mutations in 38% (3/8). Human SHH-MB *TP53* mutation rate 13–21%, with similar mutation spectrum.
- p53 pathway is deranged beyond mutation: nuclear p53 accumulation in 50% (3/6) of *p53*-WT MBs, and 50% (7/14) of *p53*-WT tumours failed to transactivate Pai-1, p21 and Noxa after UV.
- Functional proof: p53^R270C^ lentivirus into *Ptch1*^+/lacZ^ GCPs → larger lesions, **fewer p16⁺ cells**, accelerated MB (survival p = 0.0003). Genetic version *Ptch1*^+/−^;*p53*^+/R172H^ (Math1-Cre) accelerated MB (p < 0.0001) and was **stronger than *Ptch1*^+/−^;*p53*^+/−^** — i.e. dominant-negative/GOF, not simple null.

**⚠ p16/*Cdkn2a* is no longer optional (revises §1).** *Cdkn2a* was dropped from the model as "silent in GNP" (3 → 417 counts). Defensible for a two-endpoint comparison, but it is one of the two escape routes: **60% (6/10)** of *p53*-WT advanced MBs have low p16 mRNA, with **promoter methylation in 5/6** of those — no copy-number change, no mutations. Human SHH-MB: **22%** show elevated *CDKN2A* regulatory-region methylation; **3.6% (10/266)** show CDKN2A loss, subgroup-specific (0% WNT, 0.6% G3, 0.3% G4). Together p53 mutation (36%) and *Cdkn2a* inactivation (38%) explain senescence evasion in **74%** of advanced tumours; 64% of human SHH-MBs carry mutations compromising the P53 or CDKN2A pathway.
→ The mechanism is **promoter methylation — the same mechanism as the *Cdkn1c* hypothesis, in the same tumours.** If the Gli→Dnmt arm is real it has two substrates, and *CDKN2A* has published human methylation data to check against.

**What this paper does NOT address.** No p27, no p57, no p18 were examined. It does not speak to the p27⁺ Atoh1⁺ population — whatever those cells are, they are not the p21/p16 senescent population described here, which is a *preneoplastic* phenomenon that advanced tumours have selected against.

### 7.5c The miR-17~92 axis — and why our own data break the Shh→N-Myc reading

**Who defined the axis.** Northcott et al. (*Cancer Res* 2009;69(8):3249–55, **PMID 19351822**) is the paper. The defining experiment is CGNP transduction with N-myc, stabilised N-myc, Gli1 or Gli2: **N-myc, but not Gli, drives miR-17/92 expression.** Plus amplification in 6% of pediatric MBs, components the most highly up-regulated miRNAs across 90 primary tumours, expression highest in the SHH subgroup, and ectopic miR-17/92 enabling CGNPs to proliferate **without Shh**. The collaboration paper is Uziel et al. (*PNAS* 2009, **PMID 19196975**). The p21 link itself is neither MB paper but Ivanovska et al. (*Mol Cell Biol* 2008;28(7):2167–74) on the miR-106b seed family and *CDKN1A*.

**Our polyA bulk RNA-seq (mature miRNAs not captured; host transcript only):**

| | E15 | P7 wt | P14 | P28 Ptch | MB | MB+vismo |
|---|---|---|---|---|---|---|
| *Mir17hg* | 318 | 354 | **134** | 848 | 666 | 802 |

| contrast | log2FC | fold | padj |
|---|---|---|---|
| MB vs P7 wt | +0.75 | 1.68× | **1.0e-03 ✓** |
| P28 Ptch vs P7 Ptch | +1.51 | 2.84× | 0.25 (ns, n=2) |
| MB vs P28 Ptch | −0.82 | 0.56× | NA |
| MB vismo vs veh | +0.35 | 1.28× | 0.58 (ns) |
| P7 wt vismo vs veh | +0.09 | 1.07× | 0.91 (ns) |

Mature *Mir17*, *Mir18*, *Mir19a*, *Mir20a*, *Mir19b-1*, *Mir92-1* = 0 in every group (as expected for polyA-seq), as are the miR-106b~25 and miR-106a~363 paralogs.

**⚠ *Mir17hg* is a direct E2F target, not only a MYC target.** Sylvestre et al. (*JBC* 2007;282:2135–43, **PMID 17135249**) mapped a 1 kb promoter fragment containing **four E2F-binding sites** and showed endogenous E2F1, E2F2 and E2F3 all bind it and activate transcription, while miR-20a reciprocally represses E2F2/E2F3 translation. O'Donnell et al. (*Nature* 2005;435:839–43, **PMID 15944709**) showed c-Myc activates the cluster while miR-17-5p/miR-20a repress E2F1 — c-Myc simultaneously activating E2F1 transcription and limiting its translation. Reviewed as a negative feedback loop, extended to miR-106b~25, in Mendell (*Cell* 2008;133:217–22). See also He et al. (*Nature* 2005;435:828–33, **PMID 15944707**) for the cluster as an oncogene.

→ **This explains the vismodegib result.** The host has two independent inputs: MYC/MYCN and E2F. Under vismodegib *Mycn* falls but the tumour is still cycling, so E2F keeps *Mir17hg* up. Northcott's axis was defined in **acutely Shh-stimulated CGNPs**; it does not follow that the same input dominates in an established tumour.

→ **Our timecourse is itself a test.** E15 318 → P7 354 → **P14 134** → P28 Ptch 848 → MB 666. The P14 trough is exactly where GNPs exit the cycle and differentiate. That is a *proliferation* profile, not an Shh profile. **Action:** correlate *Mir17hg* across all six groups against an E2F-target signature (Mcm2-7, Ccne1, Pcna, Mki67). Tight correlation ⇒ it is a proliferation readout, demonstrable in our own data without external ChIP.

**⚠ Same confound class as *Cdkn2d* (§7.4b).** Both *Mir17hg* and *Cdkn2d* are E2F targets whose elevation in MB may **report proliferation rather than constitute an independent regulatory input**. Two rows in the matrix now share this hazard; treat it as a class, not case by case.

**Three caveats on the measurement.**
1. **The host is a two-sided proxy.** After Drosha/DGCR8 cleavage the host transcript is destabilised, so steady-state *Mir17hg* reflects transcription *minus* processing rate. If processing efficiency rises in MB, host counts can fall while mature miRNAs rise. The 1.68× is a lower bound in one direction and potentially wrong-signed in the other.
2. **"Shh-independent" overstates the statistics.** ns at 1.28× with these group sizes is absence of evidence. Write "no detectable suppression by vismodegib" and state the effect-size bound. Same for P28 Ptch 2.84× at padj 0.25, n=2 — the more interesting number and the one we cannot claim.
3. **Which member matters.** miR-17/20a/106b/93 share the seed targeting *CDKN1A*; miR-19a/b carry most of the transforming activity via PTEN. The host cannot distinguish them and they support different stories. Flag small-RNA-seq with that framing, not "per-miRNA resolution" generally.

**⚠ To verify before asserting:** miR-17~92 → *CDKN1C*/p57 is **not** established to my knowledge. The *CDKN1A* link is solid (Ivanovska 2008); p57 repression is more associated with miR-221/222 and miR-25. Check before writing "represses p21 and p57" into any figure legend — p57 is the species we most care about.

**Heatmap placement.** Add as its own row, labelled **"Mir17hg (pri-miR-17~92 host)"**, *not* folded into the cell-cycle module: it is a host transcript rather than the effector, and its E2F-target status means placing it among the CKI regulators implies a causal role it may not have in these samples.

**Prediction this buys.** If *Mir17hg* stays elevated under vismodegib, **p21 protein should stay repressed under vismodegib** — the repression outlives the Shh signal. Testable by p21 IHC or western on vismo-treated tumours, and it discriminates "p21 is Shh-repressed" from "p21 is proliferation-repressed."

**(d) Model directives.**

| finding | directive |
|---|---|
| p21 not a documented EGL protein; undetectable in MB IHC | `p21a` ≈ 0 in the GNP condition is **mechanistically correct**, not a numerical convenience. The earlier runaway on raising p21a to "data levels" was avoided for the right reason. |
| Shh→N-Myc→miR-17~92 ⊣ p21 | Do not scale `kSyp21a` by the 2.8× transcript fold without a repression term. The current decision (pin by fold, let protein emerge) is defensible; document the mechanism. |
| p53 → senescence, not apoptosis, in SHH-MB | **Corrected.** Keep `kSyp21aP53*P53`; it feeds a senescence arm. Senescence is an absorbing state — model it as such, not as reversible G0. |
| Proliferation non-monotonic in Hh drive (§7.5b) | **Highest-value new test.** Sweep Hh drive; a dip should appear at high drive with the p53/p16 arm intact and vanish without it. |
| p21 near EGL baseline in advanced MB | Consistent with `p21a` ≈ 0 in both GNP *and* advanced-MB conditions. The 2.8× transcript rise is absorbed post-transcriptionally. |
| p16/*Cdkn2a* is one of two senescence-escape routes | Reinstate *Cdkn2a* if the preneoplasia→MB transition is ever modelled. Silencing is by promoter methylation — shared mechanism with the *Cdkn1c* hypothesis. |
| p21 competence differs by *TP53* status | If the MB condition is ever split, *TP53*-mutant vs wild-type is the biologically real axis (matches the WHO SHH-activated *TP53*-mutant / wild-type split). |

### 7.6 Ranking for Stage 3

| CKI | role in a transient G0 | priority |
|---|---|---|
| **p27** | **validated decision variable via the cyclin D1/p27 ratio** (§6.1) | **highest — fix CdP21 first** |
| **p57** | plausible *source* of bimodality via monoallelic silencing | **high — needs a regulatory hook and a bimodal representation** |
| p21 | canonical source, but triply suppressed in this lineage (§7.5) and possibly replaced by a p53→apoptosis arm | **low in GNP; conditional on *TP53* status in MB** |
| p18 | timer/division-fate variable, not a G0 marker (§7.3, §7.4) | medium — needs the timer form |
| p19 | E2F target; elevation may be a proliferation readout (§7.4b) | **decide: couple to E2F or drop** |

---

## 8. Parameters — status

| param | value | status | note |
|---|---|---|---|
| `p18`,`p19` (targets) | 0.464 / 0.36 (GNP) | data-anchored | per-condition tones set in `validate` |
| `kDe_p18`,`kDe_p19` | 0.003 / 0.03 | estimate | p19 split is cosmetic as implemented (§5.6) |
| `kSyp21a`,`kSyp57a` | 2e-5 / 1e-5 | placeholder | near-degenerate (§5.9); p57a needs a hook (§7.2) |
| `kSyp21aP53` | 1e-4 | placeholder | P53≈0 default |
| `kDep21a`,`kDep57a` | 0.02 | estimate | own basal turnover |
| `kDeKPC` | 0.05 | **WRONG — §5.1** | → 0.004–0.010, free pool only |
| `kSeqCd`,`kRelCd`,`w_ink4` | 0.05 / 0.01 / 1.0 | placeholder | re-derive after §5.2 makes the buffer stoichiometric |
| `kSyP21` (p27 synthesis) | inherited | **do not set from transcript — §2.1a** | missing translational feedback |
| `w_p27` | — | **remove `P21` from it — §5.3** | also not equipotent across CIP/KIPs |
| inherited Heldt rates | — | **all are p21 parameters** | see provenance warning, §2 |

---

## 9. Validation status & Stage 3 plan

- Default (flag off): **30/32** (unchanged).
- Flag on, placeholder rates: builds + integrates; GNP cycles 17.2 h / MB 23.6 h / GNP arrests at SHH=0 — but p27 over-cleared (`P21`≈0, `CdP21`≈0), so Shh-dependence is not riding on CdP21. **§5.1 is the likely sole cause.**

> **⚠ New constraint on the MB cycle time.** Uziel 2005 (**PMID 16260494**) reports FACS DNA-content analysis of *uncultured* cells: the S-phase fraction of *Ptc1*^+/−;*Ink4c*^+/− MBs was **21.4% and 21.5%**, closely matching P5 GNPs (**22%, 23%**) and not P10 GNPs (**12%, 11%**). For comparable Ts this implies MB Tc ≈ peak-GNP Tc, against the current model's MB-slower-by-37% output (23.6 h vs 17.2 h). The constraint is tighter than it first appears: whole-population S-fraction includes any G0 cells, so if MB contains a quiescent compartment the *cycling* MB cells must have an even higher S-fraction. Reconcile before treating 23.6 h as a qualitative pass. (Caveat: n = 2 tumours, mouse, and a *Ptc1*/*Ink4c* genotype rather than ours.)

**Recommended order (changed from v1):**

1. **§5.1** — fix `kDeKPC` magnitude and restrict to the free pool; exempt `CdP21`. Re-run. This alone tests whether the rest of the structure is sound.
2. **§5.2** — make CdP21 stoichiometric, `CdP21` into the drive numerator.
3. **§5.3** — remove `P21` from the Rb denominator; resolve the double-representation.
4. **§5.4** — audit the division reset for mass conservation.
5. **§5.5** — make `with_p21_pip_degron` the default.
6. **§6.1 test** — add inherited cyclin D1 variance at division, `p21a` = 0, look for a CDK2-low subpopulation.
7. **§5.7, §5.8, §5.6, §5.9** — release products, INK4 competition, p19 decision, degeneracy check.
8. *Then* calibrate: wire per-condition tones into `validate`; fit `kDeKPC, kSeqCd, kRelCd, w_ink4, kSyp21a, kSyp57a` (+ tones) so that (i) p27 sits at a realistic level, (ii) commitment routes through CdP21, (iii) validation recovers ~30/32, (iv) Shh-dependence, EZH2 cKO, and the GNP/MB dichotomy verify.

**Bench work that would constrain the model:**
- **pT187-p27 IHC** to resolve the p27+ Atoh1+ population. Standard p27 antibodies correlate inversely with MIB-1 (Spearman R = −0.55) whereas pThr187-p27 correlates *positively* (R = 0.88); proliferating cells stain only for pT187-p27 and are unreactive with conventional antibodies (PMC554775). If our p27+ Atoh1+ cells are positive with a conventional antibody, they carry unphosphorylated p27 — genuinely inhibitory.
- **p27 nuclear vs cytoplasmic scoring**, given Ser10/CRM1 export driven by the mTOR/GSK-3 axis operating in CGNPs and MB (Bhatia 2009).
- **p18 / Atoh1 / pRb triple** (§7.3vi) — note the prior has flipped: on the Uziel data, p18 should be in Atoh1+/cycling cells, not the Atoh1− islands. Mouse p18^Ink4c immunostaining is documented as difficult; Uziel 2005 used M168 (Santa Cruz) on human tissue microarrays.
- **Sox2 added to the p27/Atoh1 panel** — if the p27+ Atoh1+ cells are Sox2−, that is a second quiescent state absent from the published hierarchy.
- **Mine the Min & Spencer 2019 MCF10A RNA-seq** (CDK2^low vs cycling vs forced quiescence) for CDKN1A/CDKN1B levels — the direct answer to "how much p27 do the reference cells have."

---

## 10. References

**Convention:** PMIDs marked ✅ were read off a retrieved PubMed record and are also given **inline in the body text** at each point of use. Entries without ✅ give journal/volume/DOI only — the PMID was not independently verified, so none is quoted anywhere in this document rather than risk a wrong number. Those are findable by title.

**Verified inline PMIDs, quick list:** 24075009 (Spencer 2013) · 25267623 (Overton 2014) · 28317845 (Barr 2017) · 34320337 (Fan & Meyer 2021) · 28869970 (Yang 2017) · 32255427 (Yang 2020) · 27368103 (Cappell 2016) · 29875408 (Cappell 2018) · 8596954 (Hengst & Reed 1996) · 31831640 (Guiley 2019) · 16260494 (Uziel 2005) · 16479172 (Uziel 2006) · 16864777 (Zindy 2006) · 19147535 (Ayrault 2009) · 17363588 (Zindy 2007) · 24515438 (Cerqueira 2014) · 41272295 (CDK4/6-CDK2-ERK differentiation commitment) · 10851091 (Thullberg 2000) · 29531090 (Kumar 2018) · 21765927 (Carcagno 2011) · 15122268 (Yuan 2004) · 19351822 (Northcott 2009) · 19196975 (Uziel 2009) · 26997276 (Tamayo-Orrego 2016) · 25263791 (Mille 2014) · 17135249 (Sylvestre 2007) · 15944709 (O'Donnell 2005) · 15944707 (He 2005)

**Proliferation–quiescence decision / transient G0**
- Spencer SL, Cappell SD, Tsai FC, Overton KW, Wang CL, Meyer T. The proliferation-quiescence decision is controlled by a bifurcation in CDK2 activity at mitotic exit. *Cell* 2013;155(2):369–83. **PMID 24075009** ✅
- Overton KW, Spencer SL, Noderer WL, Meyer T, Wang CL. Basal p21 controls population heterogeneity in cycling and quiescent cell cycle states. *PNAS* 2014;111(41):E4386–93. **PMID 25267623** ✅
- Barr AR, Cooper S, Heldt FS, Butera F, Stoy H, Mansfeld J, Novák B, Bakal C. DNA damage during S-phase mediates the proliferation-quiescence decision in the subsequent G1 via p21 expression. *Nat Commun* 2017;8:14728. **PMID 28317845** ✅
- **Fan Y, Meyer T. Molecular control of cell density-mediated exit to quiescence. *Cell Rep* 2021;36(4):109436. PMID 34320337** ✅ — *the cyclin D1/p27 ratio paper; §6.1*
- Yang HW, Chung M, Kudo T, Meyer T. Competing memories of mitogen and p53 signalling control cell-cycle entry. *Nature* 2017;549(7672):404–8. **PMID 28869970** ✅
- Yang HW, Cappell SD, Jaimovich A, Liu C, Chung M, Daigh LH, et al. Stress-mediated exit to quiescence restricted by increasing persistence in CDK4/6 activation. *eLife* 2020;9:e44571. **PMID 32255427** ✅
- Arora M, Moser J, Phadke H, Basha AA, Spencer SL. Endogenous replication stress in mother cells leads to quiescence of daughter cells. *Cell Rep* 2017;19:1351–64.
- Min M, Spencer SL. Spontaneously slow-cycling subpopulations of human cells originate from activation of stress-response pathways. *PLoS Biol* 2019;17(3):e3000178.
- Cappell SD, Chung M, Jaimovich A, Spencer SL, Meyer T. Irreversible APC^Cdh1 inactivation underlies the point of no return for cell-cycle entry. *Cell* 2016;166(1):167–80. **PMID 27368103** ✅
- Cappell SD, Mark KG, Garbett D, Pack LR, Rape M, Meyer T. EMI1 switches from being a substrate to an inhibitor of APC/C^CDH1 to start the cell cycle. *Nature* 2018;558(7709):313–7. **PMID 29875408** ✅
- Heldt FS, Barr AR, Cooper S, Bakal C, Novák B. A comprehensive model for the proliferation–quiescence decision in response to endogenous DNA damage in human cells. *PNAS* 2018;115(10):2532–7. doi:10.1073/pnas.1715345115
- Yao G et al. — quiescence depth / Rb-E2F threshold: Kwon JS et al. *Cell Rep* 2017;20:3223–35; Fujimaki K et al. *PNAS* 2019;116:22624–34; Wang X et al. *Nat Commun* 2017;8:321.

**p27 biochemistry and regulation**
- Hengst L, Reed SI. Translational control of p27^Kip1 accumulation during the cell cycle. *Science* 1996;271(5257):1861–4. **PMID 8596954** ✅
- Polyak K et al. Cloning of p27^Kip1, a cyclin-dependent kinase inhibitor... *Cell* 1994;78:59–66.
- Pagano M et al. Role of the ubiquitin–proteasome pathway in regulating abundance of p27. *Science* 1995;269:682–5.
- **Guiley KZ, Stevenson JW, Lou K, et al. p27 allosterically activates cyclin-dependent kinase 4 and antagonizes palbociclib inhibition. *Science* 2019;366(6471):eaaw2106. PMID 31831640** ✅
- Grimmler M, Wang Y, Mund T, et al. Cdk-inhibitory activity and stability of p27^Kip1 are directly regulated by oncogenic tyrosine kinases. *Cell* 2007;128:269–80.
- Besson A, Hwang HC, Cicero S, et al. Discovery of an oncogenic activity in p27^Kip1 that causes stem cell expansion and a multiple tumor phenotype. *Genes Dev* 2007;21:1731–46.
- Besson A, Gurian-West M, Chen X, et al. A pathway in quiescent cells that controls p27^Kip1 stability, subcellular localization, and tumor suppression. *Genes Dev* 2006;20(1):47–64. (p27^S10A and p27^CK− knock-ins.)
- Kamura T, Hara T, Matsumoto M, et al. Cytoplasmic ubiquitin ligase KPC regulates proteolysis of p27^Kip1 at G1 phase. *Nat Cell Biol* 2004;6:1229–35.
- Sheaff RJ, Groudine M, Gordon M, Roberts JM, Clurman BE. Cyclin E-CDK2 is a regulator of p27^Kip1. *Genes Dev* 1997;11:1464–78.
- Malek NP, Sundberg H, McGrew S, et al. A mouse knock-in model exposes sequential proteolytic pathways that regulate p27^Kip1 in G1 and S phase. *Nature* 2001;413:323–7.
- p27-pThr187 IHC: "p27^Kip1 is expressed in proliferating cells in its form phosphorylated on threonine 187." *(PMC554775 — search by title for PMID.)*

**Skp2 / APC-Cdh1 / CRL4-Cdt2**
- Wei W, Ayad NG, Wan Y, Zhang GJ, Kirschner MW, Kaelin WG. Degradation of the SCF component Skp2 in cell-cycle phase G1 by the anaphase-promoting complex. *Nature* 2004;428:194–8.
- Bashir T, Dorrello NV, Amador V, Guardavaccaro D, Pagano M. Control of the SCF^Skp2–Cks1 ubiquitin ligase by the APC/C^Cdh1 ubiquitin ligase. *Nature* 2004;428:190–3.
- Kamura T, Hara T, Kotoshiba S, et al. Degradation of p57^Kip2 mediated by SCF^Skp2-dependent ubiquitylation. *PNAS* 2003;100(18):10231–6.
- Abbas T, Dutta A. p21 in cancer: intricate networks and multiple activities. *Nat Rev Cancer* 2009;9:400–14.
- Havens CG, Walter JC. Mechanism of CRL4^Cdt2, a PCNA-dependent E3 ubiquitin ligase. *Genes Dev* 2011;25:1568–82.

**INK4 biology and G0-marker validity (§7.4)**
- Thullberg M, Bartek J, Lukas J. Ubiquitin/proteasome-mediated degradation of p19^INK4d determines its periodic expression during the cell cycle. *Oncogene* 2000;19(24):2870–6. **PMID 10851091** ✅
- Kumar A, Gopalswamy M, Wolf A, Brockwell DJ, Hatzfeld M, Balbach J. Phosphorylation-induced unfolding regulates p19^INK4d during the human cell cycle. *PNAS* 2018;115(13):3344–9. **PMID 29531090** ✅
- Carcagno AL, Ogara MF, Sonzogni SV, et al. E2F1-mediated upregulation of p19^INK4d determines its periodic expression during cell cycle and regulates cellular proliferation. *PLoS One* 2011;6(7):e21938. **PMID 21765927** ✅
- Yuan Y, Shen H, Franklin DS, Scadden DT, Cheng T. In vivo self-renewing divisions of haematopoietic stem cells are increased in the absence of the early G1-phase inhibitor, p18^INK4C. *Nat Cell Biol* 2004;6(5):436–42. **PMID 15122268** ✅
- Tourigny MR, Ursini-Siegel J, Lee H, et al. CDK inhibitor p18^INK4c is required for the generation of functional plasma cells. *Immunity* 2002;17(2):179–89.
- Yu H, Yuan Y, Shen H, Cheng T. Hematopoietic stem cell exhaustion impacted by p18^INK4C and p21^Cip1/Waf1 in opposite manners. *Blood* 2006;107(3):1200–6.
- Xu Y, et al. Small-molecule inhibitors targeting INK4 protein p18^INK4C enhance ex vivo expansion of haematopoietic stem cells. *Nat Commun* 2015;6:7328.
- Li Y, Pal R, Sung L-Y, et al. An opposite effect of the CDK inhibitor p18^INK4c on embryonic stem cells compared with tumor and adult stem cells. *PLoS One* 2012;7(9):e45212.
- Miller I, Min M, Yang C, Tian C, Gookin S, Carter D, Spencer SL. Ki67 is a graded rather than a binary marker of proliferation versus quiescence. *Cell Rep* 2018;24(5):1105–1112.e5.
- Pennycook BR, Barr AR. Palbociclib-mediated cell cycle arrest can occur in the absence of the CDK inhibitors p21 and p27. *Open Biol* 2021;11(11):210125.

**p21 in the granule lineage and MB (§7.5)**
- Northcott PA, Fernandez-L A, Hagan JP, et al. The miR-17/92 polycistron is up-regulated in Sonic hedgehog-driven medulloblastomas and induced by N-myc in Sonic hedgehog-treated cerebellar neural precursors. *Cancer Res* 2009;69(8):3249–55. **PMID 19351822** ✅
- Uziel T, Karginov FV, Xie S, et al. The miR-17~92 cluster collaborates with the Sonic Hedgehog pathway in medulloblastoma. *PNAS* 2009;106(8):2812–7. **PMID 19196975** ✅
- Murphy BL, Obad S, Bihannic L, et al. Silencing of the miR-17~92 cluster family inhibits medulloblastoma progression. *Cancer Res* 2013;73(23):7068–78.
- Ivanovska I, Ball AS, Diaz RL, et al. MicroRNAs in the miR-106b family regulate p21/CDKN1A and promote cell cycle progression. *Mol Cell Biol* 2008;28(7):2167–74.
- Miyazawa K, Himi T, Garcia V, Yamagishi H, Sato S, Ishizaki Y. A role for p27/Kip1 in the control of cerebellar granule cell precursor proliferation. *J Neurosci* 2000;20(15):5756–63.
- Zindy F, Nilsson LM, Nguyen L, et al. Hemangiosarcomas, medulloblastomas, and other tumors in *Ink4c*/*p53*-null mice. *Cancer Res* 2003;63(17):5420–7.
- Dicer1 ablation in GNPs / CDKN1A de-repression: *Cerebellum* 2016, doi:10.1007/s12311-016-0821-x.
- I2PP2A/SET suppression of p53 in SHH-MB: *Mol Cancer Res* 2019;17(1):186.
- Jeffers JR, et al. Puma is an essential mediator of p53-dependent and -independent apoptotic pathways. *Cancer Cell* 2003;4(4):321–8. / Villunger A, et al. p53- and drug-induced apoptotic responses mediated by BH3-only proteins Puma and Noxa. *Science* 2003;302(5647):1036–8.

**Senescence in SHH-MB (§7.5b)**
- Tamayo-Orrego L, Wu C-L, Bouchard N, Khedher A, Swikert SM, Remke M, Skowron P, Taylor MD, Charron F. Evasion of cell senescence leads to medulloblastoma progression. *Cell Rep* 2016;14(12):2925–37. **PMID 26997276** ✅ (open access, CC BY-NC-ND)
- Tamayo-Orrego L, Swikert SM, Charron F. Evasion of cell senescence in SHH medulloblastoma. *Cell Cycle* 2016;15(16):2102–7. (Extra View companion.)
- Mille F, Tamayo-Orrego L, Lévesque M, et al. The Shh receptor Boc promotes progression of early medulloblastoma to advanced tumors. *Dev Cell* 2014;31(1):34–47. **PMID 25263791** ✅ (Shh induces DNA damage in GCPs, requiring Boc and CyclinD1.)
- Zhukova N, Ramaswamy V, Remke M, et al. Subgroup-specific prognostic implications of TP53 mutation in medulloblastoma. *J Clin Oncol* 2013;31(23):2927–35.
- Kool M, Jones DT, Jäger N, et al. Genome sequencing of SHH medulloblastoma predicts genotype-related response to smoothened inhibition. *Cancer Cell* 2014;25(3):393–405.
- Olive KP, Tuveson DA, Ruhe ZC, et al. Mutant p53 gain of function in two mouse models of Li-Fraumeni syndrome. *Cell* 2004;119(6):847–60. (Source of the p53^R172H knockin.)
- Thomas WD, Chen J, Gao YR, et al. Patched1 deletion increases N-Myc protein stability as a mechanism of medulloblastoma initiation and progression. *Oncogene* 2009;28(13):1605–15.
- Barnes EA, Kong M, Ollendorff V, Donoghue DJ. Patched1 interacts with cyclin B1 to regulate cell cycle progression. *EMBO J* 2001;20(9):2214–23.

**miR-17~92 and the E2F feedback loop (§7.5c)**
- Northcott PA, Fernandez-L A, Hagan JP, et al. The miR-17/92 polycistron is up-regulated in Sonic hedgehog-driven medulloblastomas and induced by N-myc in Sonic hedgehog-treated cerebellar neural precursors. *Cancer Res* 2009;69(8):3249–55. **PMID 19351822** ✅
- Sylvestre Y, De Guire V, Querido E, Mukhopadhyay UK, Bourdeau V, Major F, Ferbeyre G, Chartrand P. An E2F/miR-20a autoregulatory feedback loop. *J Biol Chem* 2007;282(4):2135–43. **PMID 17135249** ✅ (four E2F sites in the miR-17~92 promoter; E2F1–3 all bind and activate.)
- O'Donnell KA, Wentzel EA, Zeller KI, Dang CV, Mendell JT. c-Myc-regulated microRNAs modulate E2F1 expression. *Nature* 2005;435(7043):839–43. **PMID 15944709** ✅
- He L, Thomson JM, Hemann MT, et al. A microRNA polycistron as a potential human oncogene. *Nature* 2005;435(7043):828–33. **PMID 15944707** ✅
- Woods K, Thomson JM, Hammond SM. Direct regulation of an oncogenic micro-RNA cluster by E2F transcription factors. *J Biol Chem* 2007;282(4):2130–4.
- Petrocca F, Visone R, Onelli MR, et al. E2F1-regulated microRNAs impair TGFβ-dependent cell-cycle arrest and apoptosis in gastric cancer. *Cancer Cell* 2008;13(3):272–86.
- Mendell JT. miRiad roles for the miR-17-92 cluster in development and disease. *Cell* 2008;133(2):217–22.
- Ventura A, Young AG, Winslow MM, et al. Targeted deletion reveals essential and overlapping functions of the miR-17 through 92 family of miRNA clusters. *Cell* 2008;132(5):875–86.

**GNP / medulloblastoma**
- **Uziel T, Zindy F, Xie S, et al. The tumor suppressors Ink4c and p53 collaborate independently with Patched to suppress medulloblastoma formation. *Genes Dev* 2005;19(22):2656–67. PMID 16260494** ✅
- **Uziel T, Zindy F, Sherr CJ, Roussel MF. The CDK inhibitor p18^Ink4c is a tumor suppressor in medulloblastoma. *Cell Cycle* 2006;5(4):363–5. PMID 16479172** ✅
- Zindy F, Knoepfler PS, Xie S, Sherr CJ, Eisenman RN, Roussel MF. N-Myc and the cyclin-dependent kinase inhibitors p18^Ink4c and p27^Kip1 coordinately regulate cerebellar development. *PNAS* 2006;103(31):11579–83. **PMID 16864777** ✅
- Ayrault O, Zindy F, Rehg J, Sherr CJ, Roussel MF. Two tumor suppressors, p27^Kip1 and Patched-1, collaborate to prevent medulloblastoma. *Mol Cancer Res* 2009;7(1):33–40. **PMID 19147535** ✅
- Zindy F, Uziel T, Ayrault O, et al. Genetic alterations in mouse medulloblastomas and generation of tumors de novo from primary cerebellar granule neuron precursors. *Cancer Res* 2007;67(6):2676–84. **PMID 17363588** ✅
- Hatton BA et al. p27^Kip1 in the ND2:SmoA1 model — progression, not initiation. *Biomarker Research* 2013;1:14.
- Vanner RJ, Remke M, Gallo M, et al. Quiescent Sox2+ cells drive hierarchical growth and relapse in Sonic Hedgehog subgroup medulloblastoma. *Cancer Cell* 2014;26(1):33–47.
- Zhang L, He X, Liu X, et al. Single-cell transcriptomics in medulloblastoma reveals tumor-initiating progenitors and oncogenic cascades during tumorigenesis and relapse. *Cancer Cell* 2019;36(3):302–318.
- Ligon KL, Huillard E, Mehta S, et al. Olig2-regulated lineage-restricted pathway controls replication competence in neural stem cells and malignant glioma. *Neuron* 2007;53(4):503–17. (OLIG2 represses *CDKN1A*.)
- Bhatia B, Northcott PA, Hambardzumyan D, et al. Tuberous sclerosis complex suppression in cerebellar development and medulloblastoma. *Cancer Res* 2009;69(18):7224–34. (mTOR/GSK-3 → cytoplasmic p27 in CGNPs and MB.)
- Nakashima K, Umeshima H, Kengaku M. Cerebellar granule cells are predominantly generated by terminal symmetric divisions. *Dev Dyn* 2015;244(6):748–58. (~15.9 h M-to-M.)
- *Preprint, not peer-reviewed:* Dexamethasone-induced p57-mediated quiescence drives chemotherapy resistance in SHH medulloblastoma. bioRxiv 2025.04.27.650711.

**Cip/Kip genetics**
- Matsumoto A, Takeishi S, Kanie T, et al. p57 is required for quiescence and maintenance of adult hematopoietic stem cells. *Cell Stem Cell* 2011;9(3):262–71.
- Zou P, Yoshihara H, Hosokawa K, et al. p57^Kip2 and p27^Kip1 cooperate to maintain hematopoietic stem cell quiescence through interactions with Hsc70. *Cell Stem Cell* 2011;9(3):247–61.
- Furutachi S, Matsumoto A, Nakayama KI, Gotoh Y. p57 controls adult neural stem cell quiescence and modulates the pace of lifelong neurogenesis. *EMBO J* 2013;32(7):970–81.
- Susaki E et al. Cip/Kip triple-knockout MEFs. *Biochem Biophys Res Commun* 2012.
- Cerqueira A, Martín A, Symonds CE, et al. Genetic characterization of the role of the Cip/Kip family of proteins as cyclin-dependent kinase inhibitors and assembly factors. *Mol Cell Biol* 2014;34(8):1452–9. **PMID 24515438** ✅

**CDK4/6 inhibition**
- Pennycook BR, Barr AR. Palbociclib-mediated cell cycle arrest can occur in the absence of the CDK inhibitors p21 and p27. *Open Biol* 2021;11(11):210125.
- **Inactivation of CDK4/6, CDK2, and ERK in G1-phase triggers differentiation commitment. PMID 41272295** ✅ (adipogenesis; p27 → reversible G1 lengthening, PPARG-driven p21+p18 → permanent exit.)
