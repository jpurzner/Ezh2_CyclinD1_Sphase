# BioRender prompt — v44 model architecture (Supp. Fig. S8A)

Reference for us (not for BioRender): `simulations/fig_v44_wiring.png`. Mechanistically current except for
the EZH2i placement (see corrections). The prompt below (≤1000 words) is a node/edge spec; paste it
into BioRender's AI assistant. It is fully standalone; do not attach the PNG, the assistant cannot parse it.

---

## PROMPT

Build a publication-quality BioRender schematic (Supplementary Fig. S8A) of an ODE model of the cell cycle in cerebellar granule-neuron progenitors (GNP) and SHH-medulloblastoma (MB). Message: Cyclin D1 integrates a mitogenic drive (Hedgehog→Gli, MYCN) against an epigenetic brake (PRC2/EZH2 → H3K27me3 at the *Ccnd1* promoter); the cycle sets its own brake (EZH2 is an E2F target; replication dilutes the mark); vismodegib arrests MB, EZH2 inhibitor rescues.

**FORMAT.** Landscape 3:2, 180 mm wide, labels ≥7 pt at print, flat BioRender icons, white background, thin dark outlines, low-saturation module fills, sans-serif. ≤32 labelled nodes. No rate constants, variable names or model-internal notation (no k_, K_, f0, "Dna 0→1", "E_div", "Cd protein").

**LAYOUT.** Four zones left→right; a cell-cycle wheel fills the right third. Main axis, one straight horizontal line, thickest stroke: Gli/MYCN → *Ccnd1* promoter → *Ccnd1* mRNA → **CYCLIN D1** → CDK4/6 → Rb → E2F → wheel. Cyclin D1 is the largest node, centred at ~45 % width.

**ZONE 1 — MITOGENIC DRIVE** (green, top-left; membrane strip along the top edge).
Nodes: SHH (ligand), Ptch1 (12-pass receptor), Smo (7TM receptor), Gli activator (TF icon), Gli1 (small "pathway readout" tag), MYCN (TF icon), vismodegib (red pill).
Edges: SHH → Ptch1 (binds); Ptch1 ⊣ Smo; Smo → Gli; Gli → *Ccnd1* promoter (thick); Gli → MYCN; Gli → Ptch1 (curved dashed arc, "negative feedback, lost in MB"); Gli → Jmjd3/Kdm6b (into zone 2); MYCN → *Ccnd1* promoter, tag "mitogen-independent floor, amplified in MB"; vismodegib ⊣ Smo.
Badges (two mini-cells at the zone bottom): "GNP: Ptch1 intact, MYCN normal" / "MB: Ptch1 loss, MYCN ×3, INK4 high".

**ZONE 2 — *CCND1* LOCUS** (purple, centre-left; visual focus).
Draw a chromatin fibre of six nucleosomes: promoter (left three, decorated with H3K27me3 flags), gene body (right three), RNA Pol II sitting on the promoter, a wavy nascent transcript.
Nodes: PRC2 (EZH2 large, EED and SUZ12 small) bound over the promoter; H3K27me3 flags; Jmjd3/Kdm6b (demethylase, scissors motif); nascent *Ccnd1* RNA; mature *Ccnd1* mRNA (exits right); EZH2 protein node at the top of the zone ("stable protein, made in S/G2, integrates S-phase length; lost on cycle exit"); tazemetostat (red pill).
Edges: EZH2 → PRC2 (assembly); PRC2 → flags, solid, "writes me3 (me1→me2→me3, last step slow)"; flags → PRC2, dashed loop, "EED reads me3, recruits PRC2 (read-write)"; Jmjd3 ⊣ flags, "erases; Gli-induced, mitogen-gated"; nascent RNA ⊣ PRC2, "transcription evicts PRC2"; flags ⊣ Pol II, "leaky repression: dampens, never silences (Pol II stays)"; tazemetostat ⊣ the write arrow ONLY, tag "catalytic: blocks writing, PRC2 stays bound, mark decays ~24 h", never on the PRC2 node; dashed arrow from the wheel's S segment → flags, "replication halves the mark each S phase"; Pol II → *Ccnd1* mRNA → CYCLIN D1 (thick teal).
Side label: "writer = cell cycle (EZH2) vs eraser = mitogen (Gli): their race sets the mark; MB carries less mark than GNP".

**ZONE 3 — RESTRICTION POINT** (blue, centre-right).
Nodes: Cyclin D1–CDK4/6 complex; Rb as three states in a row (Rb → mono-P Rb → hyper-P Rb); E2F; Cyclin E/A–CDK2; p16/p18 (INK4, brake icon); p27 (brake icon); Skp2 (SCF ubiquitin-ligase icon); cell-growth glyph (small→large cell) with a gate symbol; palbociclib (red pill).
Edges: CYCLIN D1 → CDK4/6; CDK4/6 → Rb→mono-P, "primes (step 1)"; CDK2 → mono-P→hyper-P, "commits (step 2)"; hyper-P Rb → E2F, "released"; E2F → Cyclin E/A; E2F → Skp2; E2F → EZH2 (long dashed, routed along the top to zone 2, "E2F × Cyclin E/A drive EZH2 in S/G2"); Skp2 ⊣ p27, "degrades"; p27 ⊣ CDK2 and ⊣ CDK4/6; p16/p18 ⊣ CDK4/6, tag "raised in MB: the Cyclin D1 that cycles a GNP arrests an MB"; palbociclib ⊣ CDK4/6; growth glyph → gate on the E2F→wheel path, dotted, "growth dilutes Rb, licensing commitment and S entry". Bracket Skp2/p27/E2F: "Skp2–p27 bistable commitment switch".

**ZONE 4 — CELL-CYCLE WHEEL** (right; G0 amber side-lobe off G1, G1 light blue, S teal, G2 orange, M red).
On the wheel: G1/S: "origins fire (CDK2 + growth gate)". S: "DNA replication; fork speed sets S length"; hydroxyurea (red pill) ⊣ forks, "slow forks → long S → more EZH2"; CHK1 icon ⊣ the G2/M switch, "active forks hold mitosis". G2/M: Cyclin B–CDK1 with Cdc25 (→) and Wee1 (⊣) as a small bistable-switch glyph. M→G1: APC/C–Cdc20, "degrades cyclins"; division marker (two daughter cells). One division callout: "daughters born p27-high → transient G0 (MB born higher, more G0); mass halves; H3K27me3 mark and EZH2 concentration inherited". G0 lobe tag: "p27-high, phospho-Rb-negative, reversible".

**NUMBERED LOOPS** (circled ①②③; long dashed arcs routed around zones, never through boxes):
① E2F → EZH2 → PRC2 → H3K27me3 ⊣ *Ccnd1* → Cyclin D1 → E2F. "Proliferation writes its own brake: negative feedback, a buffer not a switch."
② Wheel S segment → mark halving → less repression → more Cyclin D1 → faster cycling. "Faster cycling dilutes the brake."
③ Gli/MYCN → transcription ⊣ PRC2 → less mark → more transcription. "Double-negative: an active locus stays poised; vismodegib arrest self-reinforces."

**INSET** (bottom-right; four mini-cells with tiny wheels, green check or red cross): "MB + vismodegib → arrest (Gli drive gone; MYCN floor plus intact brake sits below threshold)"; "+ EZH2i → rescue (brake released)"; "+ CDK4/6i → no rescue (downstream of the lesion)"; "GNP − SHH → exit; EZH2 lost with the cycle, so EZH2i cannot rescue".

**ARROW GRAMMAR** (legend bottom-left, five entries): solid arrow activation/conversion; flat bar inhibition; dashed curve feedback (hours); red flat bar drug; grey dotted growth gate. Main axis 3 pt, other edges 1.5 pt, feedback 1 pt.

**COLOUR.** Green Hedgehog/MYCN; purple PRC2/EZH2/mark; teal *Ccnd1* mRNA and Cyclin D1; blue restriction point; wheel phase colours as above; drugs red outline; growth grey.

**DO NOT DRAW.** Separate cell-mass or size-threshold boxes (the growth glyph covers them); a basal mark-turnover loop; commitment-carryover notes; me1/me2/me3 as separate nodes; multi-line paragraphs inside the figure; zone-crossing arrows other than the main axis, E2F→EZH2, S→mark and loops ①–③.

**DELIVER.** Editable BioRender file; SVG and PDF at 180 mm; a second export with loops ①–③ and the inset hidden (main-text variant).

---

## Corrections relative to the current matplotlib wireframe

- **EZH2i placement.** The current figure draws tazemetostat inhibiting the PRC2 complex node.
  Since commit `5b6c22d` EZH2i is a catalytic inhibitor: it blocks the write step only, PRC2 stays
  bound, and de-repression is paced by mark decay. The prompt puts the drug on the write arrow.
- **Label collisions to avoid.** "Gli→Jmjd3: mitogen-gated eraser" sits on the mRNA → Cyclin D1
  arrow; the legend covers CHK1; "translate (k_Cd_translation)" overlaps the Cyclin D1 box; the
  "S-PHASE ÷2" annotation crosses the growth panel.
- **Figure number.** `docs/FIGURES.md` lists the wiring diagram as Supp. Fig. 7A;
  `simulations/PAPER_MODEL_SECTIONS_v44.md` numbers it S8A. The prompt uses S8A.
