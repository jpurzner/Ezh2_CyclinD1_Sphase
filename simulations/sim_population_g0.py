"""Population layer: GNP vs MB transient-G0 fractions + a quiescence-onset-vs-mitogen curve + Spencer sister
concordance, from the DATA-GROUNDED distribution machinery (not a placeholder, not the shelved mitogen tracker).

Reuses the sim_inheritance_partition engine (per-lineage sim, Dna-reset division detection, per-cycle G0 dwell)
but runs it for BOTH cell types with the calibrated per-lineage draws (scale-free, medians only):
    k_Cd_translation = KTL0 * LogN(0, 0.633)   (CyclinD1 abundance CV 0.70, data)
    kTlEZ            = KTLEZ0 * LogN(0, 0.51)   (EZH2 abundance CV ~0.55, data)
    birth-p27 = MED[cond] * LogN(0, 0.32) INHERITED (CV 0.33, data; shared by sisters) * small partition noise
GNP birth-p27 median 0.6, MB 1.39 (P21_DIV_MB after the 2026-07-14 joint re-opt). MB adds its INK4/mitogen identity.
The birth-p27 cell heterogeneity is INHERITED (drawn once per lineage, shared by sisters) with only small
independent PARTITION noise per division -- this is what gives Spencer's ~98% sister concordance (see PART).

NB birth-p27 is a per-TYPE IDENTITY constant -- SHH (mitogen) does NOT set it. So panel B is a
quiescence-ONSET curve (arrest/pause fraction vs mitogen at a fixed birth-p27), NOT Overton's mother-mitogen-
history -> birth-p21 -> daughter-fate map. The GNP<MB contrast follows from the hand-set medians + shared
lognormal draws, i.e. population heterogeneity around an identity constant -- do not present it as an emergent
mitogen-graded Overton split (that would need a data-correct mitogen->birth-p27 link, which does not yet exist;
the mitogen tracker was backwards for MB and is shelved).

Three readouts:
  (A) snapshot G0 time-fraction, GNP vs MB (mean over cycling cells of the fraction of time in transient G0 --
      the flow/live-imaging-comparable metric, cf. Moser 2018 CDK2-reporter 1-31%). NOT the old per-lineage
      "ever-paused" classification, which over-reported (~57% for MB) by counting a whole lineage as CDK2low
      if ANY cycle dwelled.
  (B) quiescence-onset: snapshot G0 fraction vs mitogen (SHH), at fixed per-type birth-p27.
  (C) Spencer: sister concordance -- two daughters share the INHERITED birth-p27 + abundance draws, differ only
      by small independent partition noise; fraction that land in the SAME fate (~98% when the median is away
      from the fate threshold, which the inheritance now ensures).

Run:  ./venv/bin/python simulations/sim_population_g0.py [--n=120] [--workers=9] [--pilot] [--fresh]
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import multiprocessing as mp
import tellurium as te
from src.build_model_v44_heldt import build_model_v44

SEL = ['time', 'MPF', 'P21', 'aRc', 'Dna', 'Cd']
KTL0, KTLEZ0 = 0.801, 0.004
# CyclinD1 + p27 spreads are deliberately kept WIDE/flexible: their quantification is not trustworthy (EGL cell
# packing is too tight to count p27-high cells accurately, and culture is confounded by high differentiation
# rates), so we do NOT pin them to a precise measured CV. Wider p27/CyclinD1 spread also gives GNP a realistic
# small p27-high (differentiating) tail rather than an implausible 0% (JP: GNPs have some p27-high cells in the
# outer EGL, fewer than MB). CyclinD1 variation additionally softens the (otherwise sharp ~1.35) per-cell G0 threshold.
CD_SDLOG, P27_SDLOG, EZ_SDLOG = 0.70, 0.42, 0.51      # CyclinD1 ~CV0.85, p27 ~CV0.44 (flexible), EZH2 CV0.55
# Birth-p27 has TWO noise sources, and separating them is what makes Spencer concordance work:
#   * INHERITED cell-to-cell heterogeneity (CV 0.33 = P27_SDLOG): drawn once per lineage, SHARED by sisters
#     (a mother's daughters inherit ~the same p27 setpoint) -> sisters land in the SAME fate -> high concordance.
#   * PARTITION asymmetry (PART, small): independent per sister/division -> the ~2% of sisters that DISagree.
# (Previously birth-p27 was redrawn independently every division at sigma 0.30 with no inheritance, so sisters
#  were uncorrelated -> concordance collapsed to ~50% whenever the median sat near the fate threshold.)
PART = 0.08                                           # SMALL partition-asymmetry sigma (independent per sister;
                                                      # Spencer sisters are highly similar -> ~98% concordance -> small)
SETTLE, MEAS, CHUNK, NP = 5000.0, 12000.0, 60.0, 9
DWELL_CUT, P27_HI = 2.0, 0.1                          # transient if mean G0 dwell >= 2h; G0 = pre-S & p27>0.1
# cell-type identities (per-condition params) + birth-p27 median
# Birth-p27 population medians. The DETERMINISTIC validation setpoints are 0.6 (GNP) / 1.39 (MB); here the GNP
# median is raised toward the differentiating tail (still < the ~1.35 pause threshold, so the median GNP cell still
# cycles -- GNP validation is unaffected) so the population shows the p27-high GNP minority JP observes. Exact values
# are intentionally soft (see the CV note above); the load-bearing claim is qualitative: GNP < MB, both nonzero.
GNP = dict(params=dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464, kSyP21=0.002), med=0.9)
MB  = dict(params=dict(Ptch1_copy_number=0.1, MYCN_amplification=2.8, p16=0.306, p18=1.553, kSyP21=0.002), med=1.39)

_RR = None
def _init(_):
    global _RR
    _RR = te.loada(build_model_v44())
    _RR.integrator.setValue('relative_tolerance', 1e-6)
    try: _RR.integrator.setValue('maximum_num_steps', 300000); _RR.integrator.setValue('maximum_time_step', 5.0)
    except Exception: pass


def _run_lineage(shh, params, ktl, ktlez, bp, part, seed):
    """One lineage. bp = the INHERITED birth-p27 setpoint (already carries the CV-0.33 cell heterogeneity,
    shared by sisters). Returns (dwells_hours_per_cycle, g0_time_fraction): the per-cycle G0 dwell hours
    AND the fraction of complete-cycle time spent in transient G0 (the flow/live-imaging-comparable snapshot
    metric, cf. Moser 2018 CDK2-reporter). Arrest -> (empty, nan)."""
    rng = np.random.default_rng(seed)
    _RR.reset()
    _RR['SHH'] = shh
    for k, v in params.items(): _RR[k] = v
    _RR['k_Cd_translation'] = float(ktl); _RR['kTlEZ'] = float(ktlez)
    _RR['P21_div'] = float(bp)
    draw = lambda: bp * np.exp(rng.normal(-part * part / 2.0, part)) if part > 0 else bp
    T, P21, ARC, DNA = [], [], [], []
    t = 0.0; TEND = SETTLE + MEAS
    while t < TEND:
        if part > 0: _RR['P21_div'] = float(draw())         # SMALL partition redraw around the inherited bp
        r = None
        for atol in (1e-9, 1e-8, 1e-7):
            try:
                _RR.integrator.setValue('absolute_tolerance', atol)
                r = _RR.simulate(t, t + CHUNK, NP, selections=SEL); break
            except Exception:
                r = None
        if r is None: break
        if t + CHUNK > SETTLE:
            for A, s in [(T, 'time'), (P21, 'P21'), (ARC, 'aRc'), (DNA, 'Dna')]:
                A.append(r[s][1:])
        t += CHUNK
    if not T: return (np.array([]), np.nan)
    T = np.concatenate(T); P21 = np.concatenate(P21); ARC = np.concatenate(ARC); DNA = np.concatenate(DNA)
    dt = np.median(np.diff(T))
    div = np.where((DNA[:-1] > 0.9) & (DNA[1:] < 0.1))[0]   # canonical division = Dna reset
    if len(div) < 2: return (np.array([]), np.nan)
    inS = (ARC > 0.05) & (DNA < 0.98); inG2 = DNA >= 0.98; preS = ~inS & ~inG2
    g0 = preS & (P21 > P27_HI)
    dwells = np.array([float(g0[a:b].sum() * dt / 60.0) for a, b in zip(div[:-1], div[1:])])
    g0_frac = float(g0[div[0]:div[-1]].mean())             # fraction of complete-cycle TIME in G0 (snapshot)
    return (dwells, g0_frac)


def _work(task):
    i, shh, tag, params, ktl, ktlez, bp, rep = task
    # partition noise is rep-dependent (sisters differ); the inherited bp is shared (rep-independent).
    # DETERMINISTIC seed (no Python hash() -- that is randomized across runs, making the figure irreproducible).
    seed = (int(i) * 100003 + int(round(shh * 1000)) * 131 + (1 if tag == 'MB' else 0) * 17 + int(rep) * 7) & 0xFFFFFFFF
    dwells, g0_frac = _run_lineage(shh, params, ktl, ktlez, bp, PART, seed)
    if len(dwells) == 0:
        return (i, shh, tag, rep, 'arrest', np.nan, 0)
    cls = 'transient' if np.mean(dwells) >= DWELL_CUT else 'immediate'
    return (i, shh, tag, rep, cls, g0_frac, len(dwells))    # r[5] = snapshot G0 time-fraction (Moser-comparable)


def main():
    def arg(f, d):
        for a in sys.argv:
            if a.startswith(f + '='): return a.split('=', 1)[1]
        return d
    N = int(arg('--n', 120)); WORKERS = int(arg('--workers', max(1, mp.cpu_count() - 1)))
    if '--pilot' in sys.argv: N = 24
    SHH_SWEEP = [0.2, 0.3, 0.4, 0.5, 0.7]
    rng = np.random.default_rng(11)
    ktl_i = np.clip(KTL0 * np.exp(rng.normal(0, CD_SDLOG, N)), 0.15, 5.0)     # shared per-lineage draws (paired
    ktlez_i = np.clip(KTLEZ0 * np.exp(rng.normal(0, EZ_SDLOG, N)), 0.0008, 0.02)  # GNP/MB use the SAME abundance cell)
    p27z_i = rng.normal(0, 1, N)                           # shared INHERITED birth-p27 z (CV 0.33) -- sisters inherit it
    bp = lambda med, i: float(med * np.exp(P27_SDLOG * p27z_i[i]))   # inherited birth-p27 setpoint for lineage i
    tasks = []
    for tag, cond in (('GNP', GNP), ('MB', MB)):
        for shh in SHH_SWEEP:                              # Overton: mitogen sweep
            for i in range(N):
                tasks.append((i, shh, tag, cond['params'], float(ktl_i[i]), float(ktlez_i[i]), bp(cond['med'], i), 0))
        for i in range(N):                                 # Spencer: rep=1 at SHH=0.5 (paired sister -- shares bp, diff partition noise)
            tasks.append((i, 0.5, tag, cond['params'], float(ktl_i[i]), float(ktlez_i[i]), bp(cond['med'], i), 1))
    print(f"population G0: N={N}/type x {len(SHH_SWEEP)} SHH + sisters on {WORKERS} workers ({len(tasks)} runs) ...", flush=True)
    with mp.get_context('spawn').Pool(WORKERS, initializer=_init, initargs=(None,)) as pool:
        res = pool.map(_work, tasks, chunksize=4)

    def cell(tag, shh, rep):   # -> list of rows for a condition
        return [r for r in res if r[2] == tag and abs(r[1] - shh) < 1e-6 and r[3] == rep]
    # HEADLINE metric = the flow/live-imaging-comparable SNAPSHOT G0 time-fraction (mean over cycling cells of the
    # fraction of time each spends in transient G0), directly comparable to Moser 2018 (CDK2-reporter, 1-31%). This
    # REPLACES the old per-lineage "ever-paused" classification, which counted a whole lineage as CDK2low if ANY
    # cycle dwelled and so massively over-reported (~57% for MB).
    def snap(rows):
        v = [r[5] for r in rows if not np.isnan(r[5])]
        return float(np.mean(v)) if v else np.nan
    # per-lineage fate TYPE (for the composition bars + Spencer concordance): transient pausers (re-enter, Moser)
    # vs immediate vs permanent arrest -- arrest is NOT lumped into CDK2low (a category error).
    cdk2low = lambda rows: np.mean([r[4] == 'transient' for r in rows]) if rows else np.nan
    arrestf = lambda rows: np.mean([r[4] == 'arrest' for r in rows]) if rows else np.nan

    print("\n(A) POPULATION G0 (SHH=0.5) -- snapshot G0 time-fraction (flow/Moser-comparable) + fate composition:")
    for tag in ('GNP', 'MB'):
        rows = cell(tag, 0.5, 0)
        comp = {c: np.mean([r[4] == c for r in rows]) for c in ('immediate', 'transient', 'arrest')}
        print(f"  {tag}: snapshot CDK2low {snap(rows)*100:.0f}% | fate: immediate {comp['immediate']*100:.0f}% / "
              f"transient {comp['transient']*100:.0f}% / arrest {comp['arrest']*100:.0f}%")

    print("\n(B) QUIESCENCE-ONSET -- snapshot G0 fraction vs mitogen (SHH), fixed per-type birth-p27:")
    print(f"  {'SHH':>5} | " + ' '.join(f'{s:>5}' for s in SHH_SWEEP))
    for tag in ('GNP', 'MB'):
        print(f"  {tag:>5} | " + ' '.join(f'{snap(cell(tag,s,0))*100:4.0f}%' for s in SHH_SWEEP))

    print("\n(C) SPENCER sister concordance (same lineage draws + inherited identity, independent partition noise):")
    conc = {}
    for tag in ('GNP', 'MB'):
        a = {r[0]: r[4] for r in cell(tag, 0.5, 0)}; b = {r[0]: r[4] for r in cell(tag, 0.5, 1)}
        both = [i for i in a if i in b]
        conc[tag] = float(np.mean([a[i] == b[i] for i in both])) if both else np.nan
        print(f"  {tag}: {conc[tag]*100:.0f}% same-fate sisters (n={len(both)} pairs; target ~98%)")

    # ---- figure: fig_v44_population_g0 (paper) ----
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    fig.suptitle("Population layer: GNP vs medulloblastoma transient-G0 -- quiescence-onset vs mitogen + Spencer sister concordance",
                 fontsize=13, fontweight='bold', y=1.0)
    fig.text(0.5, 0.905, f"N={N}/type. p27 & CyclinD1 spreads kept wide (EGL/culture G0 quantification is untrustworthy) -> GNP keeps a small p27-high "
             "minority, not 0%; birth-p27 heterogeneity inherited by sisters. Snapshot G0 fraction: absolute scale uncertain, GNP < MB is the claim.",
             ha='center', fontsize=8, color='#555', style='italic')
    COL = {'GNP': '#1b9e77', 'MB': '#762A83'}
    for j, tag in enumerate(('GNP', 'MB')):                                  # (A) snapshot G0 fraction at SHH=0.5
        v = snap(cell(tag, 0.5, 0)) * 100
        ax[0].bar(j, v, color=COL[tag], edgecolor='k')
        ax[0].text(j, v + 0.6, f'{v:.0f}%', ha='center', fontsize=12, fontweight='bold')
    ax[0].axhspan(1, 31, color='#f5b041', alpha=0.12); ax[0].text(1.45, 31, 'Moser 1-31%', fontsize=7.5, va='bottom', ha='right', color='#a0680a')
    ax[0].set_xticks([0, 1]); ax[0].set_xticklabels(['GNP', 'MB']); ax[0].set_ylim(0, 38); ax[0].set_ylabel('snapshot CDK2low G0 time-fraction (%)')
    ax[0].set_title('(A) Snapshot G0 fraction (SHH=0.5)\nsmall in GNP, larger in MB (both quiescent minorities)', fontsize=10.5, fontweight='bold')
    for tag in ('GNP', 'MB'):                                               # (B) quiescence onset vs mitogen
        y = [snap(cell(tag, s, 0)) * 100 for s in SHH_SWEEP]
        ax[1].plot(SHH_SWEEP, y, '-o', color=COL[tag], lw=2, ms=5, label=tag)
    ax[1].axhspan(1, 31, color='#f5b041', alpha=0.12); ax[1].text(SHH_SWEEP[-1], 31, 'Moser 1-31%', fontsize=7.5, va='bottom', ha='right', color='#a0680a')
    ax[1].set_xlabel('mitogen (SHH)'); ax[1].set_ylabel('snapshot CDK2low G0 time-fraction (%)'); ax[1].set_ylim(-2, 45)
    ax[1].set_title('(B) Quiescence onset vs mitogen (SHH)', fontsize=10.5, fontweight='bold'); ax[1].legend(fontsize=9); ax[1].grid(alpha=0.25)
    for j, tag in enumerate(('GNP', 'MB')):                                 # (C) Spencer sister concordance
        ax[2].bar(j, conc[tag] * 100, color=COL[tag], edgecolor='k')
        ax[2].text(j, conc[tag] * 100 + 1, f'{conc[tag]*100:.0f}%', ha='center', fontsize=10, fontweight='bold')
    ax[2].axhline(98, color='k', ls='--', alpha=0.5); ax[2].text(1.45, 98, 'Spencer ~98%', fontsize=8, va='bottom', ha='right')
    ax[2].set_xticks([0, 1]); ax[2].set_xticklabels(['GNP', 'MB']); ax[2].set_ylim(0, 105); ax[2].set_ylabel('same-fate sisters (%)')
    ax[2].set_title('(C) Spencer sister concordance', fontsize=10.5, fontweight='bold')
    plt.tight_layout()
    out = os.path.join(os.path.dirname(__file__), 'fig_v44_population_g0')
    plt.savefig(out + '.png', dpi=170); plt.savefig(out + '.pdf'); plt.close()
    print("\nSaved: fig_v44_population_g0.png / .pdf")


if __name__ == '__main__':
    main()
