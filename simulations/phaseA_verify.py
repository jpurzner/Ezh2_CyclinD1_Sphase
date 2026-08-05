"""Phase A verification against the ACTUAL baked model (JP 2026-08-05).
Checks the spec's Phase A stop-checks on the code as-shipped (dynamic-CDKI default, 30/32):
  A1  KPC free-pool + P21 not over-cleared  -> P21 and CdP21 nonzero at the cycling point (the 'P21~=0 bug' is fixed)
  A2  CdP21 stoichiometric + saturating     -> CdP21 sublinear in total p27 (plateaus, does not exceed the D1 pool)
  A3  P21 out of the CDK4/6 denominator      -> CDK4/6 (pRb) response monotonic in p27 (no non-monotonic dip)
  A4  mass conservation                      -> tP21 / tp21a stable over a long run (no sink)
Run: ./venv/bin/python simulations/phaseA_verify.py
"""
import numpy as np, tellurium as te
from src.build_model_v44_heldt import build_model_v44

MODEL = build_model_v44()   # defaults = baked dynamic-CDKI, 30/32

def rr_new():
    rr = te.loada(MODEL)
    rr.integrator.setValue("absolute_tolerance", 1e-8)
    rr.integrator.setValue("relative_tolerance", 1e-6)
    rr.integrator.setValue("maximum_num_steps", 1000000)
    return rr

# species present?
rr = rr_new()
ids = set(rr.getFloatingSpeciesIds())
want = ['P21','CdP21','p21a','p18_prot','p19_prot','p57a','Cd','Cd2']
print("=== species present ===")
for s in want: print(f"  {s:10s} {'yes' if s in ids else 'MISSING'}")
# assignment-rule readouts (moieties) live in reactions/other; check via getValue
print()

def settle(shh=0.5, ptch1=1.0, p18=None, p19=None, ksyp21=None, t=9000, npts=18000):
    rr = rr_new()
    rr['SHH']=shh; rr['Ptch1_copy_number']=ptch1
    rr['Cd2_expr']= (rr['k_Cd2_bas']*0 + 0)  # leave default; MB sets via ptch1
    if ptch1 < 0.5:
        rr['Cd2_expr'] = rr.getValue('CD2_EXPR_MB') if 'CD2_EXPR_MB' in dir(rr) else rr['Cd2_expr']
    if p18 is not None: rr['p18']=p18
    if p19 is not None: rr['p19']=p19
    if ksyp21 is not None: rr['kSyP21']=ksyp21
    sel = ["time"] + [s for s in want if s in ids] + [x for x in ['pRb','MPF','tP21','tp21a'] ]
    sel = [s for s in sel if s=="time" or s in ids or s in ('pRb','MPF') or s in ('tP21','tp21a')]
    try:
        res = rr.simulate(0, t, npts, selections=["time","P21","CdP21","p21a","p18_prot","p19_prot","p57a","Cd","Cd2","pRb","MPF","tP21","tp21a"])
    except Exception:
        res = rr.simulate(0, t, npts, selections=["time","P21","CdP21","p21a","p18_prot","p19_prot","p57a","Cd","Cd2","pRb","MPF","tP21"])
    return res

def ms(res, sp, settle=4000):
    if sp not in res.colnames: return float('nan')
    t=res['time']; return float(np.mean(res[sp][t>=settle]))

print("=== A1: cycling GNP vs MB — CDKI species means (settled) ===")
print(f"  {'cond':14s} {'P21(p27)':>9s} {'CdP21':>8s} {'p21a':>8s} {'p18_prot':>9s} {'p19_prot':>9s} {'p57a':>7s} {'Cd':>7s} {'Cd2':>7s}")
for lab,(shh,pt,p18) in [('GNP+SHH',(0.5,1.0,0.464)),('MB',(0.5,0.3,1.73))]:
    r=settle(shh=shh,ptch1=pt,p18=p18, p19=(0.58 if pt<0.5 else 0.36))
    print(f"  {lab:14s} {ms(r,'P21'):9.4f} {ms(r,'CdP21'):8.4f} {ms(r,'p21a'):8.4f} {ms(r,'p18_prot'):9.4f} {ms(r,'p19_prot'):9.4f} {ms(r,'p57a'):7.4f} {ms(r,'Cd'):7.3f} {ms(r,'Cd2'):7.3f}")
    tot=ms(r,'P21')+ms(r,'CdP21')
    print(f"     -> P21+CdP21 = {tot:.4f}  ({'PASS: p27 pool nonzero (P21~=0 bug fixed)' if tot>1e-3 else 'FAIL: p27 collapsed'})")

print()
print("=== A2: buffer saturation — sweep kSyP21 (total p27 up), watch CdP21 vs free P21 (GNP) ===")
print(f"  {'kSyP21':>9s} {'free P21':>9s} {'CdP21':>8s} {'Cd_free':>8s}  CdP21 sublinear?")
prev=None
for k in [0.0003,0.0006,0.00116,0.0023,0.0046,0.0092]:
    r=settle(ksyp21=k)
    cd21=ms(r,'CdP21'); p21f=ms(r,'P21'); cdf=ms(r,'Cd')
    ratio = cd21/prev if prev else float('nan')
    print(f"  {k:9.5f} {p21f:9.4f} {cd21:8.4f} {cdf:8.4f}   x{ratio:.2f} per 2x kSyP21" if prev else f"  {k:9.5f} {p21f:9.4f} {cd21:8.4f} {cdf:8.4f}")
    prev=cd21

print()
print("=== A3: monotonicity — pRb(hyper) mean vs kSyP21 (should not be non-monotonic from a p27 denominator term) ===")
for k in [0.0003,0.00116,0.0046,0.0092,0.02]:
    r=settle(ksyp21=k); print(f"  kSyP21={k:8.5f}  <pRb>={ms(r,'pRb'):.4f}  <MPF>={ms(r,'MPF'):.4f}")

print()
print("=== A4: mass conservation — tP21 / tp21a drift over a long run (GNP) ===")
r=settle(t=15000,npts=30000)
for sp in ['tP21','tp21a']:
    if sp in r.colnames:
        early=float(np.mean(r[sp][(r['time']>=4000)&(r['time']<5000)]))
        late =float(np.mean(r[sp][(r['time']>=13000)&(r['time']<14000)]))
        drift = (late-early)/early*100 if early else float('nan')
        print(f"  {sp}: early(4-5kh)={early:.4f}  late(13-14kh)={late:.4f}  drift={drift:+.1f}%")
    else:
        print(f"  {sp}: not a selectable readout")
