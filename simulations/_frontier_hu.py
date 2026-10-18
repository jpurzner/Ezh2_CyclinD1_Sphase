"""Definitive HU-S vs HU-G2 Pareto frontier. Fixes the optimizer-best NON-HU knobs (eraser off, vismodegib
fixed, fold ok) and finely sweeps the HU knobs (KmHU_fire x vmin_fork x kSyDna). Records MB HU S fold and
MB HU G2 fold at each point and flags any where BOTH pass (HU_S >= 0.816 AND HU_G2 <= 0.322). If none pass
both, the anti-correlation is a real model-structure limit (single fork-speed + BrdU threshold), not a
parameter miss. Parallel (ProcessPool). Run: ./venv/bin/python simulations/_frontier_hu.py"""
import os, json, subprocess, itertools
from concurrent.futures import ProcessPoolExecutor, as_completed
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(ROOT, 'venv', 'bin', 'python')
VAL = os.path.join(ROOT, 'simulations', 'validate_v44.py')
# optimizer-best NON-HU knobs (27/28; vismodegib 0.125, fold 4.53, eraser off)
BASE = {"kEZbas": 0.0007532662046698461, "kEZbas_Cd": 0.0003300678454462545, "kEZE2f": 0.005263697336374053,
        "Kez_cd": 4.7622239862331845, "K_prc2": 0.00487639924098308, "a_rw_prc2": 0.002,
        "a0_prc2": 0.00055375281726531, "f0_prc2": 0.1650540828720158, "n_prc2": 3.2059901282096046,
        "del_mk": 0.00042598321696347383, "kDeP21Cd": 0.7133316541581528, "k_Cd_tx_Gli_max": 44.213004058727,
        "k_Cd_tx_MYCN": 12.635298820913174, "k_Cd_tx_basal": 0.46574115766809326, "k_jmjd3_gli": 0.0}
HU_S_PASS, HU_G2_PASS = 0.816, 0.322   # validate bands: 1.36*0.6 and 0.23*1.4
GRID = list(itertools.product(
    [0.20, 0.30, 0.40, 0.55, 0.75],      # KmHU_fire (firing block)
    [0.08, 0.10, 0.13, 0.16, 0.20],      # vmin_fork (fork speed under HU; extend above default 0.1)
    [0.045, 0.055]))                     # kSyDna


def one(args):
    kmf, vmin, ksy = args
    p = {**BASE, 'KmHU_fire': kmf, 'vmin_fork': vmin, 'kSyDna': ksy}
    rp = os.path.join(ROOT, 'simulations', f'_fr_{kmf}_{vmin}_{ksy}.json')
    env = dict(os.environ, TWO_STEP_RB='1', H3K27_CHAIN='1', VALIDATE_PARAMS=json.dumps(p), VALIDATE_RESULTS=rp)
    try:
        subprocess.run([PY, VAL], env=env, cwd=ROOT, capture_output=True, timeout=600)
        out = json.load(open(rp)); os.remove(rp)
    except Exception as e:
        return (kmf, vmin, ksy, None, None, None)
    ch = {c['name']: c for c in out['checks']}
    hus = ch.get('MB HU S fold', {}).get('actual', float('nan'))
    hug2 = ch.get('MB HU G2 fold', {}).get('actual', float('nan'))
    mbs = ch.get('MB S count% (flow; BrdU Ts~3h)', {}).get('actual', float('nan'))
    return (kmf, vmin, ksy, hus, hug2, mbs)


if __name__ == '__main__':
    print(f"Pareto sweep: {len(GRID)} points. PASS-both = HU_S>={HU_S_PASS} AND HU_G2<={HU_G2_PASS}")
    print(f"{'KmHU':>5} {'vmin':>5} {'kSy':>6} | {'HU_S':>6} {'HU_G2':>6} {'MB_S%':>6} | both?")
    rows = []
    with ProcessPoolExecutor(max_workers=6) as ex:
        futs = {ex.submit(one, g): g for g in GRID}
        for f in as_completed(futs):
            rows.append(f.result())
    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    both = []
    for kmf, vmin, ksy, hus, hug2, mbs in rows:
        if hus is None:
            print(f"{kmf:>5.2f} {vmin:>5.2f} {ksy:>6.3f} | ERROR"); continue
        ok = (hus >= HU_S_PASS) and (hug2 <= HU_G2_PASS)
        if ok: both.append((kmf, vmin, ksy, hus, hug2))
        print(f"{kmf:>5.2f} {vmin:>5.2f} {ksy:>6.3f} | {hus:>6.3f} {hug2:>6.3f} {mbs:>6.2f} | {'YES <==' if ok else ''}")
    print(f"\n{'='*50}\nPoints passing BOTH HU-S and HU-G2: {len(both)}")
    if both:
        for b in both: print("  ", b)
    else:
        print("  NONE -> HU-S/HU-G2 anti-correlation is a MODEL-STRUCTURE limit (single fork-speed + BrdU threshold).")
    # closest joint miss (min over max-normalized violation)
    def viol(r):
        _, _, _, hus, hug2, _ = r
        if hus is None: return 9e9
        return max(0, (HU_S_PASS - hus) / HU_S_PASS) + max(0, (hug2 - HU_G2_PASS) / HU_G2_PASS)
    best = min(rows, key=viol)
    print(f"closest joint point: KmHU={best[0]} vmin={best[1]} kSy={best[2]} -> HU_S={best[3]:.3f} HU_G2={best[4]:.3f} (total viol {viol(best):.3f})")
