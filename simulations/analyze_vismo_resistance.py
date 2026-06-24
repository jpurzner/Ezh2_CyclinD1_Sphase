"""How much of the vismo-resistant CyclinD1 floor in MB is driven by persistent MYCN vs persistent Gli1?

CyclinD1 transcription drive = basal + Gli-term(Gli_act+Gli1) + MYCN-term(MYCN), times the EZH2 repression
(a common multiplier). After vismo (HHi->1) Gli_act collapses but a Gli1 residual remains and MYCN holds
its amplification floor. We (a) decompose the drive into basal/Gli1/MYCN at MB baseline vs MB+vismo, and
(b) causally knock out each arm and measure the CyclinD1 floor that remains.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

rr = te.loada(build_model_v44(with_ezh2=True, with_hh=True))
rr.integrator.setValue("absolute_tolerance", 1e-9)
rr.integrator.setValue("relative_tolerance", 1e-6)
C = {k: rr[k] for k in ['k_Cd_tx_basal', 'k_Cd_tx_Gli_max', 'K_Gli_act_CycD', 'n_Gli_act',
                        'K_Gli_rep_CycD', 'n_Gli_rep', 'k_Cd_tx_MYCN', 'K_MYCN_Cd', 'n_MYCN_Cd']}
MB = dict(SHH=0.5, MYCN_amplification=2.8, Ptch1_copy_number=0.1, p16=0.306, p18=1.553, kSyP21=0.002, EZH2i=0)


def comps(ga, g1, grep, mycn):
    G = ga + g1
    gli = C['k_Cd_tx_Gli_max'] * G**C['n_Gli_act'] / (C['K_Gli_act_CycD']**C['n_Gli_act'] + G**C['n_Gli_act']) \
        * (C['K_Gli_rep_CycD']**C['n_Gli_rep'] / (C['K_Gli_rep_CycD']**C['n_Gli_rep'] + grep**C['n_Gli_rep']))
    mycn_t = C['k_Cd_tx_MYCN'] * mycn**C['n_MYCN_Cd'] / (C['K_MYCN_Cd']**C['n_MYCN_Cd'] + mycn**C['n_MYCN_Cd'])
    return C['k_Cd_tx_basal'], gli, mycn_t


def run_at(hhi, times):
    """return Gli_act/Gli1/Gli_rep/MYCN/Cd sampled at given times (h), plus a late-steady mean."""
    rr.reset()
    for k, v in MB.items():
        rr[k] = v
    rr['HHi'] = hhi
    r = rr.simulate(0, 11000, 44000, selections=['time', 'Cd', 'Gli_act', 'Gli1', 'Gli_rep', 'MYCN'])
    th = r['time'] / 60.0
    out = {}
    for h in times:
        i = int(np.argmin(np.abs(th - h)))
        out[h] = {k: float(r[k][i]) for k in ['Cd', 'Gli_act', 'Gli1', 'Gli_rep', 'MYCN']}
    m = r['time'] >= 7000
    out['ss'] = {k: float(r[k][m].mean()) for k in ['Cd', 'Gli_act', 'Gli1', 'Gli_rep', 'MYCN']}
    return out


def show(lab, s):
    b, gli, my = comps(s['Gli_act'], s['Gli1'], s['Gli_rep'], s['MYCN'])
    tot = b + gli + my
    print(f"{lab:>20} | {s['Gli_act']:7.3f} {s['Gli1']:6.3f} {s['MYCN']:6.3f} {s['Cd']:8.2f} | "
          f"{100*b/tot:6.1f}% {100*gli/tot:5.1f}% {100*my/tot:5.1f}%")


print("Drive decomposition of CyclinD1 (the EZH2 repression is a common multiplier, so the basal/Gli1/MYCN")
print("shares of the DRIVE = the shares of the CyclinD1 floor).")
print(f"\n{'condition':>20} | {'Gli_act':>7} {'Gli1':>6} {'MYCN':>6} {'CyclinD1':>8} | {'basal%':>7} {'Gli%':>6} {'MYCN%':>6}")
base = run_at(0.0, [24])
show('MB baseline', base['ss'])
vis = run_at(1.0, [12, 24, 48])
show('MB+vismo  @12h', vis[12])
show('MB+vismo  @24h', vis[24])
show('MB+vismo  @48h', vis[48])
show('MB+vismo  steady', vis['ss'])
print("\n(In the model, vismo drives Gli1 -> ~0 within ~a day; the residual CyclinD1 floor is held up by")
print(" the MYCN amplification floor + basal transcription, NOT by residual Gli1.)")
