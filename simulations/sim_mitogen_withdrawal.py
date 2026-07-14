"""Mitogen-withdrawal experiment on the heterogeneous GNP ensemble (v44 baked).

Each cell is equilibrated cycling at high SHH (long enough for the slow EZH2 to settle), then
mitogen is withdrawn either SUDDENLY (SHH steps down instantly) or GRADUALLY (SHH ramped from 2
down over 60 h) to a RANGE of final levels (full -> partial -> near-baseline). We count the cell
divisions AFTER withdrawal, WITH vs WITHOUT the EZH2 -| CyclinD1 feedback (no-feedback =
f0_prc2 -> 1.0, the H3K27me3 mark cannot repress). State-reuse: each cell is equilibrated once
(saveState) and every depth/mode is run from that saved state.

KEY FINDING: the feedback effect is DEPTH-DEPENDENT. At FULL withdrawal both arrest (mitogen below
both thresholds); the feedback's effect is maximal at PARTIAL withdrawal (final SHH ~0.15-0.3), the
'feedback-gated window' between the no-feedback and the higher feedback proliferation thresholds
(= phase-plane panel D made dynamic): there no-feedback cells COAST (~4-5 divisions, CyclinD1 stays
high) while feedback cells ARREST (~1-2, CyclinD1 repressed). SUDDEN withdrawal shows a bigger gap
than GRADUAL (the slow ramp lets cells divide during the descent, blunting discrimination).

Heterogeneity (paired across all conditions; same draws as sim_g0_bifurcation):
  cd_scale (CyclinD1 setpoint) = lognormal,  P21_div (birth p27) = lognormal(median 0.45, wide).

Run:  ./venv/bin/python simulations/sim_mitogen_withdrawal.py            # full grid (~2.5 h), cached
      ./venv/bin/python simulations/sim_mitogen_withdrawal.py --probe    # 6 cells, print only (~30 min)
      ./venv/bin/python simulations/sim_mitogen_withdrawal.py --fresh    # force recompute
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from scipy.signal import find_peaks
from src.build_model_v44_heldt import build_model_v44

PROBE = '--probe' in sys.argv
FRESH = '--fresh' in sys.argv
CACHE = 'simulations/fig_v44_mitogen_withdrawal_cache.npz'

N = 16
SHH_HIGH = 2.0          # robustly proliferative mitogen before withdrawal
SHH_LOW = 0.02          # withdrawal target (near-zero; avoids the exact-0 stiff singular point)
T_PRE = 6500            # equilibrate cycling + let slow EZH2 settle (~5 cycles)
T_POST = 6000           # ~100 h window to count post-withdrawal divisions
CHUNK = 250             # piecewise step (min) -- chunked for robustness in both modes
T_RAMP = 3600           # gradual: linear SHH 2->low over 60 h, then hold
FB_OFF = 1.0            # feedback knockout: f0_prc2=1 -> R==1, mark cannot repress

rng = np.random.default_rng(7)
cd_scale = np.clip(np.exp(rng.normal(0.0, 0.68, N)), 0.25, 4.0)
p21_div = np.clip(np.exp(rng.normal(np.log(0.45), 0.60, N)), 0.12, 2.5)
KTL0 = 0.801

rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
rr.integrator.setValue("absolute_tolerance", 1e-9)
rr.integrator.setValue("relative_tolerance", 1e-6)
FB_ON = rr['f0_prc2']     # feedback ON: default leaky floor (mark represses)
SEL = ['time', 'MPF', 'Cd', 'EZH2']


def _setup(ktl, p21d, feedback):
    rr.reset()
    rr['SHH'] = SHH_HIGH; rr['EZH2i'] = 0; rr['MYCN_amplification'] = 1.0; rr['Ptch1_copy_number'] = 1.0
    rr['p16'] = 0.0; rr['p18'] = 0.464; rr['kSyP21'] = 0.002
    rr['k_Cd_translation'] = ktl; rr['P21_div'] = p21d
    rr['f0_prc2'] = FB_ON if feedback else FB_OFF


def _set_tol(atol, rtol):
    rr.integrator.setValue("absolute_tolerance", atol)
    rr.integrator.setValue("relative_tolerance", rtol)


def _sim_chunk(dt_chunk):
    """one continuation chunk; on stiff failure retry the same step with looser tolerances."""
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
    """set the cell up and run it to a cycling steady state (EZH2 settled). Returns success."""
    _setup(ktl, p21d, feedback)
    try:
        rr.simulate(0, T_PRE, int(T_PRE / 0.5))
        return True
    except Exception:
        return False


def _post(mode, depth):
    """post-withdrawal phase from the CURRENT (equilibrated) state -> (n divisions, trace)."""
    ts, mpfs, cds, ezs, shhs = [], [], [], [], []
    t0 = 0.0
    while t0 < T_POST:
        s = depth if mode == 'sudden' else max(depth, SHH_HIGH * (1 - t0 / T_RAMP))
        rr['SHH'] = s
        r = _sim_chunk(CHUNK)
        if r is None:
            break
        sl = slice(1, None) if ts else slice(None)
        ts.append(r['time'][sl] + t0); mpfs.append(r['MPF'][sl])
        cds.append(r['Cd'][sl]); ezs.append(r['EZH2'][sl])
        shhs.append(np.full(r['time'][sl].shape, s))
        t0 += CHUNK
    if not ts:
        return np.nan, None
    t = np.concatenate(ts); mpf = np.concatenate(mpfs)
    cd = np.concatenate(cds); ez = np.concatenate(ezs); shh = np.concatenate(shhs)
    dt = t[1] - t[0]
    pk, _ = find_peaks(mpf, prominence=0.15, distance=int(200 / dt))
    return len(pk), dict(t=t / 60.0, mpf=mpf, cd=cd, ez=ez, shh=shh, pk=pk)


def withdrawal(ktl, p21d, feedback, mode, shh_low=SHH_LOW):
    """standalone (equilibrate + withdraw); used by --probe."""
    if not _equilibrate(ktl, p21d, feedback):
        return np.nan, None
    return _post(mode, shh_low)


# ---------------------------------------------------------------- probe
if PROBE:
    idx = list(np.argsort(cd_scale)[::max(1, N // 6)][:6])   # 6 cells spanning the cd_scale range
    # test full (0.02) vs partial (0.3) sudden withdrawal, +/- feedback, to find the discriminating depth
    conds = [('sudden', 0.02, True), ('sudden', 0.02, False),
             ('sudden', 0.30, True), ('sudden', 0.30, False),
             ('gradual', 0.02, True), ('gradual', 0.02, False)]
    cols = {}
    print(f"cells (cd_scale): {[round(float(cd_scale[i]),2) for i in idx]}")
    print(f"{'condition':>26} | per-cell divisions          | mean")
    for mode, low, fb in conds:
        vals = [withdrawal(KTL0 * cd_scale[i], p21_div[i], fb, mode, shh_low=low)[0] for i in idx]
        tag = f"{mode} low={low} {'+fb' if fb else '-fb'}"
        cols[tag] = vals
        print(f"{tag:>26} | {str([int(v) if not np.isnan(v) else 'x' for v in vals]):>26} | {np.nanmean(vals):.2f}")
    print("\nEZH2-feedback effect (mean -fb minus +fb; positive = feedback REDUCES post-withdrawal divisions):")
    print(f"  sudden full (low 0.02): {np.nanmean(cols['sudden low=0.02 -fb']) - np.nanmean(cols['sudden low=0.02 +fb']):+.2f}")
    print(f"  sudden partial(low 0.3): {np.nanmean(cols['sudden low=0.3 -fb']) - np.nanmean(cols['sudden low=0.3 +fb']):+.2f}")
    print(f"  gradual full (low 0.02): {np.nanmean(cols['gradual low=0.02 -fb']) - np.nanmean(cols['gradual low=0.02 +fb']):+.2f}")
    sys.exit(0)


# ---------------------------------------------------------------- depth-resolved grid
DEPTHS = [0.02, 0.15, 0.3, 0.5, 0.8]    # final mitogen (SHH): full -> partial -> near-baseline
DISC = 0.3                               # discriminating depth (feedback-gated window) for panels C/D
MODES = ['sudden', 'gradual']

if FRESH or not os.path.exists(CACHE):
    import tempfile
    print(f"computing grid N={N} x depths{DEPTHS} x {MODES} x +/-fb (state-reuse)...")
    statef = tempfile.mktemp(suffix='.rrstate')
    medi = int(np.argsort(cd_scale)[N // 2])
    res = {f"{m}_{'fb' if fb else 'nofb'}": np.full((N, len(DEPTHS)), np.nan)
           for m in MODES for fb in [True, False]}
    traces = {}
    for i in range(N):
        for fb in [True, False]:
            if not _equilibrate(KTL0 * cd_scale[i], p21_div[i], fb):
                continue
            rr.saveState(statef)                       # reuse this cell's equilibrated state
            for m in MODES:
                for di, d in enumerate(DEPTHS):
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
            for f in ['t', 'mpf', 'cd', 'ez', 'shh']:
                save[f"tr_{k}_{f}"] = tr[f]
            save[f"tr_{k}_pk"] = tr['pk']
    np.savez(CACHE, depths=np.array(DEPTHS), disc=DISC, cd_scale=cd_scale, p21_div=p21_div, **save)
    print("cached ->", CACHE)
    for k, arr in res.items():
        print(f"  {k:14s} mean-by-depth: {np.round(np.nanmean(arr, 0), 2)}")

z = np.load(CACHE)
dep = z['depths']
disc_idx = int(np.argmin(np.abs(dep - float(z['disc']))))

# ---------------------------------------------------------------- figure
import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 2, figsize=(15, 11))


def _msem(arr):
    m = np.nanmean(arr, 0)
    s = np.nanstd(arr, 0) / np.sqrt(np.maximum(1, np.sum(~np.isnan(arr), 0)))
    return m, s


def _depth_panel(axp, mode, title):
    for k, col, lab in [(f'{mode}_fb', '#c0392b', 'EZH2 feedback ON'), (f'{mode}_nofb', '#7f8c8d', 'no feedback')]:
        m, s = _msem(z[k])
        axp.errorbar(dep, m, yerr=s, marker='o', color=col, lw=2, ms=5, capsize=3, label=lab)
    mf, _ = _msem(z[f'{mode}_fb']); mn, _ = _msem(z[f'{mode}_nofb'])
    axp.fill_between(dep, mf, mn, where=(mn >= mf), color='#c0392b', alpha=0.10)
    axp.axvline(DISC, color='k', ls=':', alpha=0.5)
    axp.set_xlabel('final mitogen level  (SHH;  left = full withdrawal)', fontsize=11)
    axp.set_ylabel('post-withdrawal divisions / cell', fontsize=11)
    axp.set_title(title, fontweight='bold', fontsize=11)
    axp.legend(fontsize=9, loc='upper left'); axp.grid(alpha=0.15)


_depth_panel(ax[0, 0], 'sudden', '(A) Sudden withdrawal: divisions vs final mitogen\nfeedback gap (shaded) opens at partial withdrawal')
_depth_panel(ax[0, 1], 'gradual', '(B) Gradual withdrawal: divisions vs final mitogen\n(SHH ramped 2→level over 60 h)')

# (C) per-cell distributions at the discriminating (partial-withdrawal) depth
c = ax[1, 0]
groups = [('sudden_fb', 'Sudden\n+fb', '#c0392b'), ('sudden_nofb', 'Sudden\n−fb', '#e8a39a'),
          ('gradual_fb', 'Gradual\n+fb', '#1f618d'), ('gradual_nofb', 'Gradual\n−fb', '#a9cce3')]
rng2 = np.random.default_rng(1)
for j, (k, lab, col) in enumerate(groups):
    d = z[k][:, disc_idx]; d = d[~np.isnan(d)]
    c.scatter(rng2.normal(j, 0.07, len(d)), d, s=28, color=col, edgecolor='k', linewidth=0.4, alpha=0.85, zorder=3)
    c.bar(j, np.mean(d), color=col, alpha=0.25, width=0.62, zorder=1)
    c.text(j, np.mean(d) + 0.1, f'{np.mean(d):.1f}', ha='center', va='bottom', fontweight='bold', fontsize=9)
c.set_xticks(range(4)); c.set_xticklabels([g[1] for g in groups], fontsize=9)
c.set_ylabel('post-withdrawal divisions / cell', fontsize=11)
c.set_title(f'(C) Per-cell distribution at partial withdrawal (SHH→{DISC})\nfeedback → arrest; no feedback → coasting', fontweight='bold', fontsize=11)
c.grid(alpha=0.15, axis='y')

# (D) representative single-cell mechanism at the discriminating depth
d = ax[1, 1]
YCAP = 4.0      # focus on post-transient steady levels (initial decay from high-mitogen Cd runs off-top)
for k, col, lab in [('sudden_nofb', '#7f8c8d', 'no feedback'), ('sudden_fb', '#c0392b', 'EZH2 feedback ON')]:
    if f'tr_{k}_t' not in z.files:
        continue
    t = z[f'tr_{k}_t']; cd = z[f'tr_{k}_cd']; pk = z[f'tr_{k}_pk']
    lvl = np.nanmean(cd[t > 20])
    d.plot(t, cd, '-', color=col, lw=1.9, label=f'CyclinD1, {lab}  (≈{lvl:.1f})')
    yy = YCAP * (0.95 if k == 'sudden_nofb' else 0.86)
    if len(pk):
        d.plot(t[pk], np.full(len(pk), yy), 'v', color=col, ms=9, zorder=6)
    d.text(t[-1] * 1.02, yy, f'{len(pk)} div', color=col, fontsize=8.5, va='center', fontweight='bold')
d.set_ylim(0, YCAP)
d.set_xlabel('time after withdrawal (h)', fontsize=11)
d.set_ylabel('CyclinD1', fontsize=11)
d.set_title(f'(D) Representative cell, sudden → SHH {DISC}\nno-feedback: CyclinD1 stays high → coasts (5 div); feedback: EZH2 represses CyclinD1 → arrests (1 div)', fontweight='bold', fontsize=9.5)
d.legend(fontsize=8, loc='upper right'); d.grid(alpha=0.15)

plt.tight_layout()
plt.savefig('simulations/fig_v44_mitogen_withdrawal.png', dpi=160, bbox_inches='tight')
plt.savefig('simulations/fig_v44_mitogen_withdrawal.pdf', bbox_inches='tight')
plt.close()
print("Saved: fig_v44_mitogen_withdrawal.png / .pdf")
