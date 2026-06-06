"""Aggregate + rank the parallel rescue searches. Run when the searches are done (or anytime):
  ./venv/bin/python simulations/v44_rescue_search_report.py

Reads all simulations/v44_rescue_search_*.jsonl, keeps guard-passing rows, ranks by objective
(distance to MB 100 / HHi 25 / RES 75 / CDK 0), and reports the best param sets + the region that works.
"""
import glob, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

rows = []
for f in sorted(glob.glob(os.path.join(os.path.dirname(__file__), "v44_rescue_search_*.jsonl"))):
    for ln in open(f):
        try:
            r = json.loads(ln)
        except Exception:
            continue
        if 'obj' in r and r.get('guard_ok'):
            rows.append(r)

if not rows:
    print("no guard-passing rows yet"); raise SystemExit

rows.sort(key=lambda r: r['obj'])
total = sum(1 for f in glob.glob(os.path.join(os.path.dirname(__file__), "v44_rescue_search_*.jsonl"))
            for _ in open(f))
print(f"{total} total evals, {len(rows)} guard-passing.  Target: MB 100 / HHi 25 / RES 75 / CDK 0\n")
print(f"{'rank':>4} {'obj':>7} {'MB':>5} {'HHi':>5} {'RES':>5} {'CDK':>4} | {'p16':>5} {'sigma':>6} {'kez':>5} {'ktl0':>5}")
for i, r in enumerate(rows[:20]):
    print(f"{i+1:>4} {r['obj']:>7.2f} {r['MB']:>5.0f} {r['HHi']:>5.0f} {r['RES']:>5.0f} {r['CDK']:>4.0f} | "
          f"{r['p16']:>5.2f} {r['sigma']:>6.2f} {r['kez']:>5.2f} {r['ktl0']:>5.2f}")

# region that works: mean +/- range over the top 15
top = rows[:15]
print("\nparameter region of the top-15 (median [min, max]):")
import statistics as st
for k in ('p16', 'sigma', 'kez', 'ktl0'):
    v = [r[k] for r in top]
    print(f"  {k:6}: {st.median(v):.3f}  [{min(v):.3f}, {max(v):.3f}]")

b = rows[0]
print("\nBEST:")
print(json.dumps({k: b[k] for k in ('p16', 'sigma', 'kez', 'ktl0', 'MB', 'HHi', 'RES', 'CDK', 'obj')}, indent=0))
print("\nTo apply: set MB p16 + K_EZH2_repression + k_Cd_translation in the builder, sigma in"
      " fig_v44_fig5_population.py (cd_scale), then re-validate (and add crash-exclusion to the figure).")
