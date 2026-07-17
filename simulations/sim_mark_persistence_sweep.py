"""Overnight sweep: can we make H3K27me3 persist LONGER (slower decay) while keeping MB<GNP (ChIP)?
Sweep transcription-eviction g_prc2 (up, to replace the Gli eraser for MB<GNP), Gli eraser k_jmjd3_gli (down, to
slow decay), and mark turnover del_mk (down, to slow decay). For each: MB/GNP steady-state mark ratio (ChIP target
<1, ideally ~0.5), MB & GNP mark half-life after EZH2i (hours, want LONGER), and the MB/GNP CyclinD1 fold (want the
5.07x MB>GNP intact)."""
import os, sys
os.environ.setdefault('TWO_STEP_RB','1'); os.environ.setdefault('H3K27_CHAIN','1')
sys.path.insert(0,'/Users/jpurzner/Dropbox/Q_research/py_projects/Ezh2_CyclinD1_Sphase')
import numpy as np, multiprocessing as mp, tellurium as te
from src.build_model_v44_heldt import build_model_v44
UPH=62.8
GNP=dict(Ptch1_copy_number=1.0, MYCN_amplification=1.0, p16=0.0, p18=0.464)
MB =dict(Ptch1_copy_number=0.1, MYCN_amplification=2.8, p16=0.306, p18=1.553)
_M=None
def _init(_):
    global _M; _M=build_model_v44()
def measure(cond, over):
    rr=te.loada(_M); rr['SHH']=0.5
    for k,v in {**cond,**over}.items(): rr[k]=v
    rr.reset(); rr['SHH']=0.5
    for k,v in {**cond,**over}.items(): rr[k]=v
    try:
        r=rr.simulate(0,9000,3000,selections=['time','Mk','Cd'])
    except Exception:
        return np.nan, np.nan, np.nan
    m=r['time']>=6000
    mk=float(np.mean(np.asarray(r['Mk'])[m])); cd=float(np.mean(np.asarray(r['Cd'])[m]))
    rr['EZH2i']=1.0
    try:
        r2=rr.simulate(9000,9000+4000,2500,selections=['time','Mk'])
    except Exception:
        return mk, cd, np.nan
    mk2=np.asarray(r2['Mk']); t2=np.asarray(r2['time'])-9000
    idx=np.where(mk2<=mk2[0]/2)[0]
    thalf=(t2[idx[0]]/UPH) if len(idx) else (4000/UPH)   # censored at 63h
    return mk, cd, thalf
def _one(a):
    g,kje,dm=a
    over={'g_prc2':g,'k_jmjd3_gli':kje,'del_mk':dm}
    mkg,cdg,hg=measure(GNP,over); mkm,cdm,hm=measure(MB,over)
    ratio=mkm/mkg if mkg>1e-6 else np.nan
    fold=cdm/cdg if cdg>1e-6 else np.nan
    return (g,kje,dm,ratio,hm,hg,fold)
if __name__=='__main__':
    grid=[(g,kje,dm) for g in (0.07,0.25,0.5,0.9) for kje in (0.0,0.03,0.08,0.22) for dm in (0.0006,0.0023)]
    with mp.Pool(8, initializer=_init, initargs=(None,)) as pool:
        res=pool.map(_one, grid)
    print(f"{'g_prc2':>7} {'k_jmjd3':>8} {'del_mk':>7} | {'MB/GNP':>7} {'MBhalf_h':>9} {'GNPhalf_h':>9} {'CdFold':>7}")
    for g,kje,dm,ratio,hm,hg,fold in sorted(res, key=lambda x:(x[3] if not np.isnan(x[3]) else 9)):
        flag=' <-- MB<GNP + longer' if (ratio<0.9 and hm>6) else ''
        print(f"{g:>7.2f} {kje:>8.2f} {dm:>7.4f} | {ratio:>7.2f} {hm:>9.1f} {hg:>9.1f} {fold:>7.2f}{flag}")
