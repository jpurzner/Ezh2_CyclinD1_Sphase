"""Targeted: in regimes where GNP sits at a genuine small-but-nonzero G0 (real cell type),
what does data-grounded MB do at the SAME shared birth-p? Report the closest MB gets to GNP."""
import numpy as np
from model import GNP, MB, run_ensemble

print("Regimes where GNP is calibrated to a small nonzero G0; MB uses SAME birth-p (Cd4,kCd2,p18=3).")
print(" also reported: p18 needed for MB to MATCH GNP, and the best MB-GNP gap found.\n")
best_gap=-9
for extra in [None, dict(Kd=2.0,n=2,dS=8.0), dict(Krb=0.8,kdE=0.5), dict(Kp=1.0,Kd=0.1)]:
    for m27 in np.linspace(0.05, 1.2, 12):
        g=run_ensemble(dict(GNP),N=1200,extra=extra,m21=0.05,cv21=0.8,m27=m27,cv27=0.5,seed=3)['g0_frac']
        if 0.01 <= g <= 0.06:
            m=run_ensemble(dict(MB),N=1200,extra=extra,m21=0.05,cv21=0.8,m27=m27,cv27=0.5,seed=3)['g0_frac']
            gap=m-g
            if gap>best_gap: best_gap=gap; best=(extra,m27,g,m)
            print(f"  extra={extra}  m27={m27:.2f}: GNP={g*100:4.1f}%  MB={m*100:4.1f}%  gap={gap*100:+.1f}pp")
print(f"\nBEST valid gap (GNP in 1-6%): {best_gap*100:+.1f}pp  at {best}")

# how much p18 does MB need just to EQUAL the calibrated GNP ~6%?
print("\np18 required for MB to reach GNP-like ~6% G0 (Cd4,kCd2, default wiring, m27=0.30):")
for p18 in [3,8,12,16,20,25]:
    m=run_ensemble(dict(Cd=4.0,kCd=2.0,p18=p18),N=1200,m21=0.05,cv21=0.8,m27=0.30,cv27=0.5,seed=3)['g0_frac']
    print(f"  p18={p18:2d}x (measured=3x): MB G0={m*100:4.1f}%")
