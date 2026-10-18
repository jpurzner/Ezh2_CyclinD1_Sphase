"""Read fd_probe.json and report step-size convergence of the central-difference elasticities."""
import json, os, numpy as np

D = os.path.dirname(os.path.abspath(__file__))
S = json.load(open(os.path.join(D, 'fd_probe.json')))
C, STEPS, base = S['central'], S['steps'], S['base']

BIN = {'GNP+SHH cycles (>0 div)', 'GNP-SHH arrest (0 div)', 'GNP+HHi arrest (0 div)',
       'GNP serum-starve arrest', 'MB cycles (>0 div)', 'MB+CDK4/6i arrest (0)',
       'MB+CDK4/6i+EZH2i no rescue'}

print('Central-difference elasticity  dlog(target)/dlog(param)  vs relative step h')
print('A smooth output plateaus as h shrinks; a non-smooth one wanders.\n')

for p in C:
    hs = sorted([h for h in C[p]], key=float, reverse=True)
    tg = [t for t in base if t not in BIN and abs(base[t]) > 1e-9]
    # rank targets by |elasticity| at the largest step
    h0 = hs[0]
    rank = sorted(tg, key=lambda t: -abs(C[p][h0].get(t, 0)))[:4]
    print(f'--- {p} ---')
    print(f"   {'target':32s} " + ' '.join(f'{float(h):>9.3f}' for h in hs))
    for t in rank:
        vals = [C[p][h].get(t, np.nan) for h in hs]
        print(f'   {t[:32]:32s} ' + ' '.join(f'{v:9.3f}' for v in vals))
    # convergence metric: spread of the 3 smallest steps relative to their mean
    small = hs[-3:]
    rows = []
    for t in tg:
        v = np.array([C[p][h].get(t, np.nan) for h in small], float)
        if np.nanmax(np.abs(v)) < 1e-6:
            continue
        rows.append(np.nanstd(v) / max(np.nanmean(np.abs(v)), 1e-12))
    if rows:
        print(f'   -> across the 3 smallest steps (h={small}): median rel. spread '
              f'{np.median(rows):.2f}, 90th pct {np.percentile(rows,90):.2f}\n')
    else:
        print('   -> no resolvable sensitivity\n')

# global: which step size gives the most self-consistent Jacobian?
print('\n=== which step size is self-consistent? ===')
print('For each h, compare its elasticity vector to the one at the next smaller h (cosine sim).')
for p in C:
    hs = sorted([h for h in C[p]], key=float, reverse=True)
    tg = [t for t in base if t not in BIN and abs(base[t]) > 1e-9]
    sims = []
    for a, b in zip(hs, hs[1:]):
        va = np.array([C[p][a].get(t, 0.0) for t in tg])
        vb = np.array([C[p][b].get(t, 0.0) for t in tg])
        n = np.linalg.norm(va) * np.linalg.norm(vb)
        sims.append(float(va @ vb / n) if n > 1e-14 else np.nan)
    print(f'  {p:20s} ' + ' '.join(f'{float(a):.3f}->{float(b):.3f}:{s:6.3f}'
                                   for (a, b), s in zip(zip(hs, hs[1:]), sims)))
