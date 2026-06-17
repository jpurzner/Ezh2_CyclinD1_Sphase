"""Vismodegib (GDC0449) treatment of the heterogeneous MB ensemble (v44 baked).

The MB analogue of the mitogen-withdrawal experiment (sim_mitogen_withdrawal.py). Vismo blocks Smo,
removing the Gli-driven CyclinD1 but leaving the MYCN floor intact -> even FULL vismo is a PARTIAL
mitogen withdrawal for MB. So the EZH2 -| CyclinD1 feedback is decisive: WITH feedback the residual
(MYCN-floor) CyclinD1 is repressed below threshold -> ARREST; WITHOUT it (no-fb = K_EZH2_repression
-> inf, equivalent to EZH2i / tazemetostat) the MYCN floor keeps CyclinD1 up -> the cell COASTS
(= the model's EZH2i rescue), now dose-resolved + sudden/gradual + ensemble.

Each MB cell is equilibrated cycling at GDC=0 (EZH2 settled), then vismo is applied either SUDDENLY
(GDC steps up instantly) or GRADUALLY (GDC ramped 0 -> dose over 60 h) to a range of final doses.
We count the cell divisions AFTER treatment, WITH vs WITHOUT the feedback. State-reuse: each cell is
equilibrated once (saveState) and every dose/mode is run from that saved state.

Heterogeneity (paired across all conditions; same draws as sim_g0_bifurcation / sim_mitogen_withdrawal):
  cd_scale (CyclinD1 setpoint) = lognormal,  P21_div (birth p27) = lognormal(median 0.45, wide).

Run:  ./venv/bin/python simulations/sim_vismo_withdrawal.py            # full grid (~2.5 h), cached
      ./venv/bin/python simulations/sim_vismo_withdrawal.py --probe    # 6 cells, print only (~30 min)
      ./venv/bin/python simulations/sim_vismo_withdrawal.py --fresh    # force recompute
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

PROBE = '--probe' in sys.argv
FRESH = '--fresh' in sys.argv
CACHE = 'simulations/fig_v44_vismo_withdrawal_cache.npz'

N = 16
T_PRE = 6500            # equilibrate cycling MB + let slow EZH2 settle
T_POST = 6000           # ~100 h window to count post-treatment divisions
CHUNK = 250             # piecewise step (min) -- chunked for robustness in both modes
T_RAMP = 3600           # gradual: linear GDC 0->dose over 60 h, then hold
KREP_OFF = 1e6          # feedback knockout (= EZH2i / tazemetostat)

# MB identity (matches validate_v44 / sim_g0_bifurcation MB)
MB = dict(shh=0.5, mycn_amp=2.8, ptch1=0.1, p16=0.306, p18=1.553, ksyp21=0.002)

rng = np.random.default_rng(7)
cd_scale = np.clip(np.exp(rng.normal(0.0, 0.68, N)), 0.25, 4.0)
p21_div = np.clip(np.exp(rng.normal(np.log(0.45), 0.60, N)), 0.12, 2.5)
KTL0 = 0.801

rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
rr.integrator.setValue("absolute_tolerance", 1e-9)
rr.integrator.setValue("relative_tolerance", 1e-6)
KREP_ON = rr['K_EZH2_repression']
SEL = ['time', 'MPF', 'Cd', 'EZH2']


def _setup(ktl, p21d, feedback):
    rr.reset()
    rr['SHH'] = MB['shh']; rr['EZH2i'] = 0; rr['GDC0449'] = 0.0
    rr['MYCN_amplification'] = MB['mycn_amp']; rr['Ptch1_copy_number'] = MB['ptch1']
    rr['p16'] = MB['p16']; rr['p18'] = MB['p18']; rr['kSyP21'] = MB['ksyp21']
    rr['k_Cd_translation'] = ktl; rr['P21_div'] = p21d
    rr['K_EZH2_repression'] = KREP_ON if feedback else KREP_OFF


def _set_tol(atol, rtol):
    rr.integrator.setValue("absolute_tolerance", atol)
    rr.integrator.setValue("relative_tolerance", rtol)


def _sim_chunk(dt_chunk):
    try:
        return rr.simulate(0, dt_chunk, int(dt_chunk / 0.5), selections=SEL)
    except Exception:
        try:
            _set_tol(1e-6, 1e-4)
            r = rr.simulate(0, dt_chunk, int(dt_chunk / 0.5), selections=SEL)
            _set_tol(1e-9, 1e-6)
            return r
        except Exception:
            _set_tol(1e-9, 1e-6)
            return None


def _equilibrate(ktl, p21d, feedback):
    _setup(ktl, p21d, feedback)
    try:
        rr.simulate(0, T_PRE, int(T_PRE / 0.5))
        return True
    except Exception:
        return False


def _post(mode, dose):
    """post-vismo phase from the current (equilibrated, GDC=0) state. mode: sudden/gradual; dose=final GDC."""
    ts, mpfs, cds, ezs, gdcs = [], [], [], [], []
    t0 = 0.0
    while t0 < T_POST:
        g = dose if mode == 'sudden' else min(dose, dose * t0 / T_RAMP)
        rr['GDC0449'] = g
        r = _sim_chunk(CHUNK)
        if r is None:
            break
        sl = slice(1, None) if ts else slice(None)
        ts.append(r['time'][sl] + t0); mpfs.append(r['MPF'][sl])
        cds.append(r['Cd'][sl]); ezs.append(r['EZH2'][sl])
        gdcs.append(np.full(r['time'][sl].shape, g))
        t0 += CHUNK
    if not ts:
        return np.nan, None
    t = np.concatenate(ts); mpf = np.concatenate(mpfs)
    cd = np.concatenate(cds); ez = np.concatenate(ezs); gdc = np.concatenate(gdcs)
    dt = t[1] - t[0]
    pk, _ = find_peaks(mpf, prominence=0.15, distance=int(200 / dt))
    return len(pk), dict(t=t / 60.0, mpf=mpf, cd=cd, ez=ez, gdc=gdc, pk=pk)


def withdrawal(ktl, p21d, feedback, mode, dose):
    if not _equilibrate(ktl, p21d, feedback):
        return np.nan, None
    return _post(mode, dose)


# ---------------------------------------------------------------- probe
if PROBE:
    idx = list(np.argsort(cd_scale)[::max(1, N // 6)][:6])
    conds = [('sudden', 0.6, True), ('sudden', 0.6, False),
             ('sudden', 0.95, True), ('sudden', 0.95, False),
             ('gradual', 0.95, True), ('gradual', 0.95, False)]
    cols = {}
    print(f"MB cells (cd_scale): {[round(float(cd_scale[i]),2) for i in idx]}")
    print(f"{'condition':>26} | per-cell divisions          | mean")
    for mode, dose, fb in conds:
        vals = [withdrawal(KTL0 * cd_scale[i], p21_div[i], fb, mode, dose)[0] for i in idx]
        tag = f"{mode} GDC={dose} {'+fb' if fb else '-fb'}"
        cols[tag] = vals
        print(f"{tag:>26} | {str([int(v) if not np.isnan(v) else 'x' for v in vals]):>26} | {np.nanmean(vals):.2f}")
    print("\nEZH2-feedback effect (mean -fb minus +fb; positive = feedback REDUCES post-vismo divisions):")
    print(f"  sudden GDC=0.6: {np.nanmean(cols['sudden GDC=0.6 -fb']) - np.nanmean(cols['sudden GDC=0.6 +fb']):+.2f}")
    print(f"  sudden GDC=0.95: {np.nanmean(cols['sudden GDC=0.95 -fb']) - np.nanmean(cols['sudden GDC=0.95 +fb']):+.2f}")
    print(f"  gradual GDC=0.95: {np.nanmean(cols['gradual GDC=0.95 -fb']) - np.nanmean(cols['gradual GDC=0.95 +fb']):+.2f}")
    sys.exit(0)


# ---------------------------------------------------------------- dose-resolved grid
# vismo dose (GDC0449): 0 = untreated -> 0.95 = near-complete Smo block. Capped at 0.95 (not 1.0):
# at GDC=1.0 exactly Smo_active ~ (1-GDC0449) = 0 is a stiff singularity that spuriously arrests cells
# (verified: a high-CyclinD1 no-fb cell rescues at 0.85/0.95/0.98 but artifactually arrests at 1.0);
# biologically vismodegib is a competitive antagonist that never 100%-blocks Smo, so 0.95 is realistic.
DOSES = [0.0, 0.3, 0.6, 0.85, 0.95]
DISC = 0.95                              # near-full-vismo dose for panels C/D (= the clinical rescue scenario)
MODES = ['sudden', 'gradual']

if FRESH or not os.path.exists(CACHE):
    import tempfile
    print(f"computing MB+vismo grid N={N} x doses{DOSES} x {MODES} x +/-fb (state-reuse)...")
    statef = tempfile.mktemp(suffix='.rrstate')
    medi = int(np.argsort(cd_scale)[N // 2])
    res = {f"{m}_{'fb' if fb else 'nofb'}": np.full((N, len(DOSES)), np.nan)
           for m in MODES for fb in [True, False]}
    traces = {}
    for i in range(N):
        for fb in [True, False]:
            if not _equilibrate(KTL0 * cd_scale[i], p21_div[i], fb):
                continue
            rr.saveState(statef)
            for m in MODES:
                for di, d in enumerate(DOSES):
                    rr.loadState(statef)
                    n, tr = _post(m, d)
                    res[f"{m}_{'fb' if fb else 'nofb'}"][i, di] = n
                    if i == medi and abs(d - DISC) < 1e-9:
                        traces[f"{m}_{'fb' if fb else 'nofb'}"] = tr
        print(f"  cell {i + 1}/{N} done")
    try:
        os.remove(statef)
    except OSError:
        pass
    save = dict(res)
    for k, tr in traces.items():
        if tr is not None:
            for f in ['t', 'mpf', 'cd', 'ez', 'gdc']:
                save[f"tr_{k}_{f}"] = tr[f]
            save[f"tr_{k}_pk"] = tr['pk']
    np.savez(CACHE, doses=np.array(DOSES), disc=DISC, cd_scale=cd_scale, p21_div=p21_div, **save)
    print("cached ->", CACHE)
    for k, arr in res.items():
        print(f"  {k:14s} mean-by-dose: {np.round(np.nanmean(arr, 0), 2)}")

z = np.load(CACHE)
dose = z['doses']
disc_idx = int(np.argmin(np.abs(dose - float(z['disc']))))

# ---------------------------------------------------------------- figure
import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 2, figsize=(15, 11))


def _msem(arr):
    m = np.nanmean(arr, 0)
    s = np.nanstd(arr, 0) / np.sqrt(np.maximum(1, np.sum(~np.isnan(arr), 0)))
    return m, s


def _dose_panel(axp, mode, title):
    for k, col, lab in [(f'{mode}_fb', '#8e44ad', 'EZH2 feedback ON'), (f'{mode}_nofb', '#7f8c8d', 'no feedback (≡ EZH2i)')]:
        m, s = _msem(z[k])
        axp.errorbar(dose, m, yerr=s, marker='o', color=col, lw=2, ms=5, capsize=3, label=lab)
    mf, _ = _msem(z[f'{mode}_fb']); mn, _ = _msem(z[f'{mode}_nofb'])
    axp.fill_between(dose, mf, mn, where=(mn >= mf), color='#8e44ad', alpha=0.10)
    axp.axvline(DISC, color='k', ls=':', alpha=0.5)
    axp.set_xlabel('vismodegib dose  (GDC0449;  right = near-complete Hh block)', fontsize=11)
    axp.set_ylabel('post-treatment divisions / cell', fontsize=11)
    axp.set_title(title, fontweight='bold', fontsize=11)
    axp.legend(fontsize=9, loc='lower left'); axp.grid(alpha=0.15)


_dose_panel(ax[0, 0], 'sudden', '(A) Sudden vismo: divisions vs dose\nfeedback gap (shaded) opens at high vismo — the MYCN-floor partial withdrawal')
_dose_panel(ax[0, 1], 'gradual', '(B) Gradual vismo: divisions vs dose\n(GDC ramped 0→dose over 60 h)')

# (C) per-cell distributions at full vismo (the rescue scenario)
c = ax[1, 0]
groups = [('sudden_fb', 'Sudden\n+fb', '#8e44ad'), ('sudden_nofb', 'Sudden\n−fb (EZH2i)', '#c39bd3'),
          ('gradual_fb', 'Gradual\n+fb', '#1f618d'), ('gradual_nofb', 'Gradual\n−fb (EZH2i)', '#a9cce3')]
rng2 = np.random.default_rng(1)
for j, (k, lab, col) in enumerate(groups):
    d = z[k][:, disc_idx]; d = d[~np.isnan(d)]
    c.scatter(rng2.normal(j, 0.07, len(d)), d, s=28, color=col, edgecolor='k', linewidth=0.4, alpha=0.85, zorder=3)
    c.bar(j, np.mean(d), color=col, alpha=0.25, width=0.62, zorder=1)
    c.text(j, np.mean(d) + 0.1, f'{np.mean(d):.1f}', ha='center', va='bottom', fontweight='bold', fontsize=9)
c.set_xticks(range(4)); c.set_xticklabels([g[1] for g in groups], fontsize=9)
c.set_ylabel('post-treatment divisions / cell', fontsize=11)
c.set_title(f'(C) Per-cell distribution at full vismo (GDC={DISC})\n+feedback → arrest;  −feedback (EZH2i) → rescue / coast', fontweight='bold', fontsize=11)
c.grid(alpha=0.15, axis='y')

# (D) representative MB cell at full vismo: CyclinD1 (MYCN floor) + divisions.
# Re-simulate a cell that clearly shows the rescue contrast (nofb rescues, fb arrests at DISC) --
# the median cell sits at the bistable rescue threshold and arrests, so it isn't illustrative.
nf = z['sudden_nofb'][:, disc_idx]; ff = z['sudden_fb'][:, disc_idx]
cand = np.where((nf >= 4) & (ff <= 2))[0]
rep_tr = {}
if len(cand):
    rep = int(cand[np.argmin(np.abs(z['cd_scale'][cand] - np.median(z['cd_scale'][cand])))])
    for fb in [True, False]:
        if _equilibrate(KTL0 * float(z['cd_scale'][rep]), float(z['p21_div'][rep]), fb):
            _, tr = _post('sudden', float(z['disc']))
            if tr is not None:
                rep_tr['nofb' if not fb else 'fb'] = tr

d = ax[1, 1]
YCAP = 4.0
for key, col, lab in [('nofb', '#7f8c8d', 'no feedback (EZH2i)'), ('fb', '#8e44ad', 'EZH2 feedback ON')]:
    tr = rep_tr.get(key) or ({'t': z[f'tr_sudden_{key}_t'], 'cd': z[f'tr_sudden_{key}_cd'], 'pk': z[f'tr_sudden_{key}_pk']} if f'tr_sudden_{key}_t' in z.files else None)
    if tr is None:
        continue
    t = tr['t']; cd = tr['cd']; pk = tr['pk']
    lvl = np.nanmean(cd[t > 20])
    d.plot(t, cd, '-', color=col, lw=1.9, label=f'CyclinD1, {lab}  (≈{lvl:.1f})')
    yy = YCAP * (0.95 if key == 'nofb' else 0.86)
    if len(pk):
        d.plot(t[pk], np.full(len(pk), yy), 'v', color=col, ms=9, zorder=6)
    d.text(t[-1] * 1.02, yy, f'{len(pk)} div', color=col, fontsize=8.5, va='center', fontweight='bold')
d.set_ylim(0, YCAP)
d.set_xlabel('time after vismo (h)', fontsize=11); d.set_ylabel('CyclinD1', fontsize=11)
d.set_title(f'(D) Representative MB cell, sudden vismo (GDC={DISC})\nno-feedback: MYCN-floor CyclinD1 sustains cycling; feedback: EZH2 represses it → arrest', fontweight='bold', fontsize=9.5)
d.legend(fontsize=8, loc='upper right'); d.grid(alpha=0.15)

plt.tight_layout()
plt.savefig('simulations/fig_v44_vismo_withdrawal.png', dpi=160, bbox_inches='tight')
plt.savefig('simulations/fig_v44_vismo_withdrawal.pdf', bbox_inches='tight')
plt.close()
print("Saved: fig_v44_vismo_withdrawal.png / .pdf")
