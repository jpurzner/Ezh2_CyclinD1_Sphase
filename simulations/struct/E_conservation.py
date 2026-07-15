"""STRUCTURE ANALYSIS — step E: conservation laws & true dimension.

Left null space of the stoichiometry matrix N (species x reactions) = the conserved moieties (fixed pools).
These are the WITHIN-CYCLE-EPOCH invariants: the E_div / Mk-dilution events reshuffle them at division.
Verified numerically (y.x constant between events). Output: the conserved-pool table + a dimension ledger.

Run: ./venv/bin/python simulations/struct/E_conservation.py
"""
import os, sys, json
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1')
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import tellurium as te
from fractions import Fraction
import sympy as sp
from src.build_model_v44_heldt import build_model_v44

OUT = os.path.dirname(os.path.abspath(__file__))
rr = te.loada(build_model_v44())
rr.conservedMoietyAnalysis = False
N = rr.getFullStoichiometryMatrix()
species = list(N.rownames)
Narr = np.array(N)
n_sp = len(species)
rank = np.linalg.matrix_rank(Narr)
n_cons = n_sp - rank
print(f"species: {n_sp}  reactions: {Narr.shape[1]}  rank(N): {rank}  conservation laws: {n_cons}")

# integer left-null basis via sympy (interpretable moieties): null space of N^T
Nsym = sp.Matrix([[sp.Rational(str(round(float(v), 6))).limit_denominator(1000) for v in row] for row in Narr])
left_null = Nsym.T.nullspace()   # vectors y in R^n_sp with y^T N = 0
print(f"sympy left-null basis vectors: {len(left_null)}\n")

# simulate once to numerically verify each conservation between events
r = rr.simulate(0, 3000, 3000, selections=['time'] + species)
traj = {s: r[s] for s in species}

def clean_int(vec):
    # scale a rational vector to smallest integer vector
    denoms = [Fraction(str(v)).limit_denominator(1000).denominator for v in vec]
    from math import gcd
    L = 1
    for d in denoms: L = L * d // gcd(L, d)
    iv = [int(round(float(v) * L)) for v in vec]
    g = 0
    for x in iv: g = gcd(g, abs(x))
    if g > 1: iv = [x // g for x in iv]
    # sign: majority positive
    if sum(1 for x in iv if x < 0) > sum(1 for x in iv if x > 0): iv = [-x for x in iv]
    return iv

moieties = []
for vec in left_null:
    iv = clean_int([vec[i] for i in range(n_sp)])
    terms = [(species[i], iv[i]) for i in range(n_sp) if iv[i] != 0]
    if not terms: continue
    # numeric total over the run (should be ~constant between events; report mean + drift)
    tot = np.zeros(len(r['time']))
    for s, c in terms: tot = tot + c * traj[s]
    expr = " + ".join(f"{c}*{s}" if c != 1 else s for s, c in terms)
    moieties.append(dict(species=[s for s, _ in terms], coeffs=[c for _, c in terms],
                         expr=expr, mean=float(np.mean(tot)),
                         cv=float(np.std(tot) / (abs(np.mean(tot)) + 1e-12))))

moieties.sort(key=lambda m: len(m['species']), reverse=True)
print("=== CONSERVED MOIETIES (within-epoch invariants) ===")
for m in moieties:
    print(f"  [{len(m['species'])} sp]  {m['expr']}")
    print(f"           total~{m['mean']:.3g}  (CV over run incl. events {m['cv']:.2%})")
print(f"\nTRUE DIMENSION: {n_sp} species - {n_cons} conservation laws = {rank} independent dof")

json.dump(dict(n_species=n_sp, rank=int(rank), n_conservation=int(n_cons), moieties=moieties),
          open(os.path.join(OUT, 'E_conservation.json'), 'w'), indent=2)
print("wrote E_conservation.json")
