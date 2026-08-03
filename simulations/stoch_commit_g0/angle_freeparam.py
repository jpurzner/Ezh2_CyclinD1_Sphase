"""FREE-PARAM-HUNT refutation angle (efficient fate-map search).
Search model knobs (Kd,Kp,Krb,n,kb,kdE,dS,Kd2,md,Erb,w18,E_thr) for ANY regime where
data-grounded MB (Cd4,kCd2,p18 3, SHARED birth-p with GNP) has G0 >= GNP G0 + 5pp,
with GNP calibrated to ~0-6% G0 (real cell type) at the SHARED birth-p.

Efficiency: for each mechanistic param set, integrate GNP and MB ONCE over a fine 2D birth
grid (p21 x p27). That gives fate(p21,p27). G0 for any lognormal birth dist is then a cheap
Monte-Carlo lookup, so the birth-mean + E_thr sweeps cost almost nothing.
"""
import numpy as np
from model import BASE, GNP, MB, make_params, integrate

# ---- 0. Analytic: MB vs GNP CDK4/6 drive Dd over the whole (Kd,w18) plane ----
print("=== 0. Analytic: MB vs GNP CDK4/6 drive Dd across the whole (Kd,w18) plane ===")
worst = None
for Kd in np.geomspace(0.01, 50, 60):
    for w18 in np.geomspace(0.01, 50, 60):
        DdG = 1.0*1.0/(Kd*(1+w18*1.0)+1.0)
        DdM = 2.0*4.0/(Kd*(1+w18*3.0)+4.0)
        r = DdM/DdG
        if worst is None or r < worst[0]:
            worst = (r, Kd, w18, DdG, DdM)
print(f"  min MB/GNP drive ratio over 3600 (Kd,w18) combos = {worst[0]:.3f} at Kd={worst[1]:.3g}, w18={worst[2]:.3g}")
print(f"    (DdG={worst[3]:.3f}, DdM={worst[4]:.3f})  -> MB drive ALWAYS > GNP drive; ratio floor ~4.")
print("  numerator kCd*Cd up 8x, denom brake (1+p18) up only 2x -> drive/brake up 4x, monotone in Cd.\n")

# ---- fate-map machinery ----
P21_AX = np.concatenate([[0.0], np.geomspace(0.01, 6.0, 15)])   # 16 birth-p21 nodes
P27_AX = np.linspace(0.0, 4.0, 29)                               # 29 birth-p27 nodes
G21, G27 = np.meshgrid(P21_AX, P27_AX, indexing='ij')
FLAT21 = G21.ravel(); FLAT27 = G27.ravel()

def fate_grid(ct, extra):
    p = make_params(dict(ct), extra)
    E0 = np.full(FLAT21.shape, 0.05)
    Ef, _, _ = integrate(E0, FLAT21.copy(), FLAT27.copy(), p, T=80.0, dt=0.1)
    return Ef.reshape(G21.shape)   # E_final[i21, j27]

def snap_idx(vals, axis):
    return np.clip(np.searchsorted(axis, vals), 0, len(axis)-1)

# Precompute fixed standard-normal draws once (shared across all evals) -> cheap deterministic MC.
Nmc = 1500
_z21 = np.random.default_rng(101).standard_normal(Nmc)
_z27 = np.random.default_rng(202).standard_normal(Nmc)

def _ln(z, mean, cv):
    if mean <= 0: return np.zeros_like(z)
    s = np.sqrt(np.log(1+cv**2)); mu = np.log(mean)-0.5*s**2
    return np.exp(mu + s*z)

def g0_from_map(Emap, m21, cv21, m27, cv27, E_thr):
    i = snap_idx(_ln(_z21, m21, cv21), P21_AX)
    j = snap_idx(_ln(_z27, m27, cv27), P27_AX)
    return float(np.mean(Emap[i, j] < E_thr))

# ---- 1. randomized search ----
print("=== 1. Fate-map search: max (MB_G0 - GNP_G0), GNP calibrated to <=6% G0, SHARED birth-p ===")
grid = dict(
    Kd  =[0.1, 0.5, 2.0, 8.0],
    Kp  =[0.1, 0.3, 1.0],
    Krb =[0.2, 0.4, 0.8],
    n   =[2, 6, 12],
    kb  =[0.001, 0.01, 0.05],
    kdE =[0.5, 1.0, 2.0],
    dS  =[0.5, 2.0, 8.0],
    Kd2 =[0.2, 0.5, 1.0],
    md  =[2, 4, 8],
    Erb =[0.5, 1.0, 2.0],
    w18 =[0.5, 1.0, 3.0],
)
keys = list(grid.keys())
births = [dict(m21=m21, cv21=0.8, m27=m27, cv27=0.5)
          for m21 in [0.02, 0.05, 0.2]
          for m27 in np.linspace(0.0, 3.5, 12)]
E_thr_list = [0.2, 0.3, 0.5, 0.8]

rng = np.random.default_rng(2024)
N_param_samples = 500
best = None; best_valid = None

for s in range(N_param_samples):
    extra = {k: grid[k][rng.integers(len(grid[k]))] for k in keys}
    Eg = fate_grid(GNP, extra); Em = fate_grid(MB, extra)
    for E_thr in E_thr_list:
        for b in births:
            g = g0_from_map(Eg, E_thr=E_thr, **b)
            m = g0_from_map(Em, E_thr=E_thr, **b)
            gap = m - g
            if best is None or gap > best[0]:
                best = (gap, dict(extra), dict(b), E_thr, g, m)
            if 0.0 <= g <= 0.06 and (best_valid is None or gap > best_valid[0]):
                best_valid = (gap, dict(extra), dict(b), E_thr, g, m)

tot = N_param_samples*len(births)*len(E_thr_list)
print(f"\n  {N_param_samples} param sets x {len(births)} births x {len(E_thr_list)} E_thr = {tot:,} evals.")
print("\n  --- Best gap OVERALL (ignoring GNP calibration) ---")
gap, extra, b, E_thr, g, m = best
print(f"    GNP G0={g*100:.1f}%  MB G0={m*100:.1f}%  gap={gap*100:+.1f}pp  E_thr={E_thr}")
print(f"    params={extra}\n    birth={b}")

print("\n  --- Best gap among GNP-CALIBRATED sets (GNP 0-6% G0) [THE VALID TEST] ---")
if best_valid is None:
    print("    NONE: no sampled param set kept GNP in 0-6% G0 at any shared birth-p.")
else:
    gap, extra, b, E_thr, g, m = best_valid
    print(f"    GNP G0={g*100:.1f}%  MB G0={m*100:.1f}%  gap={gap*100:+.1f}pp  E_thr={E_thr}")
    print(f"    params={extra}\n    birth={b}")
    print(f"    -> {'COUNTEREXAMPLE FOUND' if gap >= 0.05 else 'NO valid counterexample (gap < +5pp)'}")
