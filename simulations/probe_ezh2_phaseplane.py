"""Quick validation of the EZH2-CyclinD1 phase-plane construction (before the full figure).

Cd-nullcline (analytic, fast var): Cd = A * K_rep/(K_rep+EZH2),  A = ktl*D/(k_Cd_deg*k_Cd_mRNA_deg)
  D = mitogen drive = Cd_mRNA production with repression removed (read from steady Gli/MYCN).
EZH2-nullcline (numeric, slow var): set EZH2i=1 (feedback OFF), sweep ktl -> (mean Cd, mean EZH2).
Real orbit: full model feedback ON -> project trajectory onto (Cd, EZH2).
Prints the operating point + ranges so we can sanity-check shapes/intersection.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
rr.integrator.setValue("absolute_tolerance", 1e-9)
rr.integrator.setValue("relative_tolerance", 1e-6)

# constants
C = {k: rr[k] for k in ['k_Cd_translation', 'k_Cd_deg', 'k_Cd_mRNA_deg', 'K_EZH2_repression',
                         'k_Cd_tx_basal', 'k_Cd_tx_Gli_max', 'K_Gli_act_CycD', 'K_Gli_rep_CycD',
                         'n_Gli_act', 'n_Gli_rep', 'k_Cd_tx_MYCN', 'K_MYCN_Cd', 'n_MYCN_Cd']}
print("constants:", {k: round(v, 4) for k, v in C.items()})


def drive_D(ga, g1, grep, mycn):
    """Cd_mRNA production with EZH2 repression = 1 (the mitogen drive D)."""
    gli = ga + g1
    gli_term = C['k_Cd_tx_Gli_max'] * gli**C['n_Gli_act'] / (C['K_Gli_act_CycD']**C['n_Gli_act'] + gli**C['n_Gli_act']) \
        * (C['K_Gli_rep_CycD']**C['n_Gli_rep'] / (C['K_Gli_rep_CycD']**C['n_Gli_rep'] + grep**C['n_Gli_rep']))
    mycn_term = C['k_Cd_tx_MYCN'] * mycn**C['n_MYCN_Cd'] / (C['K_MYCN_Cd']**C['n_MYCN_Cd'] + mycn**C['n_MYCN_Cd'])
    return C['k_Cd_tx_basal'] + gli_term + mycn_term


def run(shh, ezh2i, mycn_amp, ptch1, ktl=None, p16=0.0, p18=0.464, ksyp21=0.002, T=14000, settle=8000):
    rr.reset()
    rr['SHH'] = shh; rr['EZH2i'] = ezh2i; rr['MYCN_amplification'] = mycn_amp
    rr['Ptch1_copy_number'] = ptch1; rr['p16'] = p16; rr['p18'] = p18; rr['kSyP21'] = ksyp21
    if ktl is not None:
        rr['k_Cd_translation'] = ktl
    r = rr.simulate(0, T, int(T/0.5), selections=['time', 'Cd', 'EZH2', 'Gli_act', 'Gli1', 'Gli_rep', 'MYCN', 'Cd_mRNA'])
    m = r['time'] >= settle
    return {k: r[k][m] for k in r.colnames}


GNP = dict(shh=0.5, ezh2i=0, mycn_amp=1.0, ptch1=1.0)

# --- operating point (feedback ON) ---
op = run(**GNP)
cd_op, ez_op = op['Cd'].mean(), op['EZH2'].mean()
D = drive_D(op['Gli_act'].mean(), op['Gli1'].mean(), op['Gli_rep'].mean(), op['MYCN'].mean())
A = C['k_Cd_translation'] * D / (C['k_Cd_deg'] * C['k_Cd_mRNA_deg'])
print(f"\nGNP operating point (feedback ON): Cd={cd_op:.3f}  EZH2={ez_op:.3f}")
print(f"  orbit ranges: Cd[{op['Cd'].min():.3f},{op['Cd'].max():.3f}]  EZH2[{op['EZH2'].min():.3f},{op['EZH2'].max():.3f}]")
print(f"  drive D={D:.3f}  A(=ktl*D/degs)={A:.3f}")
print(f"  Cd-nullcline at EZH2_op: {A*C['K_EZH2_repression']/(C['K_EZH2_repression']+ez_op):.3f}  (should ~= Cd_op {cd_op:.3f})")

# --- EZH2-nullcline: feedback OFF, sweep ktl ---
print("\nEZH2-nullcline (feedback OFF, sweep k_Cd_translation):")
print("   ktl   ->  meanCd   meanEZH2")
for ktl in [0.1, 0.25, 0.45, 0.7, 0.801, 1.2, 2.0, 3.5]:
    o = run(shh=0.5, ezh2i=1, mycn_amp=1.0, ptch1=1.0, ktl=ktl)
    print(f"  {ktl:5.3f}  ->  {o['Cd'].mean():6.3f}   {o['EZH2'].mean():6.3f}")
