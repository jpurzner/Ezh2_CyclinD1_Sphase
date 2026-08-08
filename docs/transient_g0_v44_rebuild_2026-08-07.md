# v44 transient-G0 rebuild — escape mechanism + a biological reframe (2026-08-07)

JP chose the integrated v44 rebuild for biological accuracy. Building it produced (a) a real, working piece — the
**re-entry / escape mechanism** — and (b) a series of robust structural findings that **reframe what an MB transient-G0
cell is**. Nothing baked; default preserved at 30/32.

## Piece built: the escape / re-entry mechanism (`with_cd_hyper_escape`, default OFF)
The escape diagnostic showed a cell that falls into G0 **locks permanently** — not because p27 stays high (free p27
clears to ~0), but because E2f is trapped on **mono**-phosphorylated Rb (RbmE2f) and the two-step Rb lets **only**
CyclinE/A hyper-phosphorylate Rb. So accumulating CyclinD1 (Cd~13) cannot complete Rb-P / release E2f → no CyclinE
bootstrap → permanent arrest. Fix (Yang 2020 eLife 44571 / Chung 2019 Mol Cell: **CDK4/6 activity alone drives
commitment**, even in CyclinE/A-quadruple-null MEFs): let high CyclinD-CDK4/6 contribute weakly to Rb **hyper**-P
(`w_cd_hyper`). This converts the permanent lock into a re-enterable state — the necessary re-entry component of any
transient G0. Implemented, builds, default off (30/32 preserved). `simulations/escape_test.py`.

## The robust finding: a p27-brake *dwell* is structurally impossible in high-drive MB
Birth-p27 (`P21_div`) has **essentially zero effect on the MB cycle** across every combination tested — option A,
option B (31/32 params), slow p27 clearance (kDeKPC 0.02→0.0003), a coherent **CDK2-low** birth (Ce_div/Ca_div → 0.001)
with high p27 (3.5), and with the escape on. MB is **24.4 h in every case, no dwell**. Reason: MB's drive is so high it
re-commits fast from *any* birth state (the same drive-dominance / "superset" result as the stochastic-commitment and
option-B work). No birth-brake can hold a high-drive MB cell in G0; a brake strong enough to hold it instead **locks**
it (permanent). There is no transient middle for the high-drive cycling bulk.

## Reframe: an MB transient-G0 cell is a LOW-DRIVE subpopulation, not a p27-brake dwell
Because the high-drive bulk cannot dwell, the transient-G0 cells must be a **low-drive subpopulation** — low
CyclinD1/Gli (hence low CyclinD1/p27 ratio → p27-high / CDK2-low → G0), which **re-enters** when drive recovers (the
escape mechanism). This is biologically well-supported:
- **Vanner 2014 (Cancer Cell):** quiescent **Sox2⁺** cells drive SHH-MB growth and relapse — a low-proliferation,
  re-entry-competent subpopulation. That *is* the transient-G0 cell.
- **Fan-Meyer:** the decision variable is the **CyclinD1/p27 ratio** — a low-CyclinD1 cell has a low ratio → G0, even
  with normal p27. So "high p27" is satisfied via a low ratio (low CyclinD1), matching the data (p27 not elevated).
- **MB-specificity is the re-entry:** a low-drive MB cell re-enters (oncogenic drive + escape); a low-drive GNP cell
  **differentiates** (permanent exit) — the same low-drive state, opposite outcome. That is the clean "MB transient
  G0, GNP none."
- **This is the same drive axis as #2 (graded G1):** at low Shh/CyclinD1, G1 lengthened to ~50 h (a de-facto transient
  G0). #3 is the low-drive *tail* of that #2 distribution — one mechanism, not two.

So the faithful transient-G0 layer is: **a population with CyclinD1/Gli (drive) heterogeneity → the low-drive tail
sits in a long-G1/G0 dwell → re-enters via the escape (MB) or exits via differentiation (GNP).** It is a
*drive-heterogeneity* population layer, not a *birth-p27* mechanism.

## Where this needs JP's call (biology)
The rebuild has shifted the mechanism from "high-p27 daughter dwells on a p27 brake" (JP's original #3 framing) to
"low-CyclinD1 subpopulation dwells at low ratio and re-enters." Both are "transient G0 with high p27," but the second
is what the model (and Vanner/Fan-Meyer) support. Before building the drive-heterogeneity population layer +
differentiation-vs-re-entry split, confirm this is the intended biology.

## Status
- Built: escape/re-entry (`with_cd_hyper_escape`), default off, 30/32 preserved.
- Established: birth-p27 dwell impossible in high-drive MB (drive-dominance).
- Reframed: transient G0 = low-drive Sox2⁺-like subpopulation + escape (MB) vs differentiation (GNP), = the low-drive
  tail of #2.
- Next (pending JP confirm): drive-heterogeneity population layer + a differentiation/permanent-exit route for the
  low-drive GNP state, so low-drive MB re-enters (transient G0) while low-drive GNP exits.
Probes: `simulations/{birthp27_probe,escape_test,escape_test2}.py`.
