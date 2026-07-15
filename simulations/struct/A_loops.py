"""STRUCTURE ANALYSIS — step A: signed influence graph + feedback loops.

Numeric species-Jacobian at several limit-cycle operating points -> a SIGNED influence digraph (edge x_j->x_i
signed by sign of d(dx_i/dt)/dx_j). The numeric Jacobian AUTOMATICALLY threads the := assignment rules (Cdc25a,
Wee1a, Chk1, commit_gate, PRC2, ...) because perturbing a species perturbs everything downstream of it -- so the
gotcha the plan warned about (dropped links) does NOT bite the numeric approach; algebraic intermediates are
compiled into direct species->species edges. We then: (1) find the strongly-connected FEEDBACK CORE vs the
feed-forward periphery, (2) enumerate short cycles + sign/classify them, (3) verify the hypothesized L1-L5/N1-N5
loops, (4) render a signed-Jacobian heatmap. Event edges (E_div, Mk-dilution) are outside the Jacobian -> noted.

Run: ./venv/bin/python simulations/struct/A_loops.py
"""
import os, sys, json
os.environ.setdefault('TWO_STEP_RB', '1'); os.environ.setdefault('H3K27_CHAIN', '1')
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np
import tellurium as te
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from src.build_model_v44_heldt import build_model_v44

OUT = os.path.dirname(os.path.abspath(__file__))
MODULES = {
    'Rb-E2f': ['Rb', 'Rbm', 'pRb', 'E2f', 'RbE2f', 'RbmE2f'],
    'CDK2-p27': ['Ce', 'Ca', 'CeP21', 'CaP21', 'P21'],
    'Skp2': ['Skp2'],
    'Replication': ['aPcna', 'iPcna', 'Rc', 'pRc', 'aRc', 'iRc', 'Dna'],
    'p53-damage': ['P53', 'Dam', 'Pr'],
    'MPF-mitotic': ['MPF', 'preMPF', 'Cdc20'],
    'Emi1-Cdh1': ['E1', 'C1', 'pC1', 'E1C1'],
    'Mass': ['mass'],
    'EZH2': ['EZH2m', 'EZH2'],
    'Hedgehog': ['SHH_Ptch', 'Ptch1_free', 'Ptch1_mRNA', 'Smo_active', 'Gli_rep', 'Gli_act',
                 'Gli1_mRNA', 'Gli1', 'Gli1_epi', 'MYCN'],
    'CyclinD1': ['Cd_mRNA', 'Cd'],
    'H3K27chain': ['m1_me', 'm2_me', 'Mk'],
}
MOD_OF = {s: m for m, ss in MODULES.items() for s in ss}

rr = te.loada(build_model_v44())
sp = rr.getFloatingSpeciesIds()
idx = {s: i for i, s in enumerate(sp)}

# --- operating points across one settled cycle ---
rr.reset()
rr.simulate(0, 3000, 3000)            # settle
# grab Jacobians at several points spanning ~one cycle
Js = []; phases = []
t = 3000.0
for k in range(8):
    rr.simulate(t, t + 175, 20)       # advance ~175 min
    t += 175
    try:
        J = np.array(rr.getFullJacobian())
        Js.append(J)
        phases.append(dict(Dna=float(rr['Dna']), MPF=float(rr['MPF']), P21=float(rr['P21']), aRc=float(rr['aRc'])))
    except Exception as e:
        pass
print(f"collected {len(Js)} Jacobians across the cycle")

# --- signed adjacency: edge j->i if |J[i,j]| exceeds a relative threshold at ANY phase; flag sign flips ---
n = len(sp)
EPS = 1e-6
sign = np.zeros((n, n), dtype=int)     # sign of j->i (i=row target, j=col source)
flip = np.zeros((n, n), dtype=bool)
Jmax = np.max([np.abs(J) for J in Js], axis=0)
scale = np.percentile(Jmax[Jmax > 0], 50)   # median nonzero magnitude for relative thresholding
for i in range(n):
    for j in range(n):
        if i == j:
            continue
        vals = [J[i, j] for J in Js if abs(J[i, j]) > EPS]
        if not vals:
            continue
        pos = sum(1 for v in vals if v > 0); neg = sum(1 for v in vals if v < 0)
        if pos and neg:
            flip[i, j] = True
        sign[i, j] = 1 if pos >= neg else -1

n_edges = int(np.sum(sign != 0))
print(f"signed edges (j->i): {n_edges}  sign-flipping edges: {int(np.sum(flip))}")

# --- build signed DiGraph (drop self-loops for cycle analysis; record them separately) ---
G = nx.DiGraph()
G.add_nodes_from(sp)
selfreg = {}
for i in range(n):
    for j in range(n):
        if i == j or sign[i, j] == 0:
            continue
        G.add_edge(sp[j], sp[i], sign=int(sign[i, j]))
for i in range(n):
    d = 0
    # diagonal self-effect sign (net autoregulation) from mean J
    dv = np.mean([J[i, i] for J in Js])
    selfreg[sp[i]] = int(np.sign(dv)) if abs(dv) > EPS else 0

# --- strongly connected components: the feedback core vs feed-forward periphery ---
sccs = [c for c in nx.strongly_connected_components(G) if len(c) > 1]
sccs.sort(key=len, reverse=True)
core = sccs[0] if sccs else set()
periphery = set(sp) - core
print(f"\nFEEDBACK CORE (largest SCC): {len(core)} species")
print("  modules in core:", sorted(set(MOD_OF.get(s, '?') for s in core)))
print(f"FEED-FORWARD periphery: {len(periphery)} species -> {sorted(periphery)}")
print("  periphery modules:", sorted(set(MOD_OF.get(s, '?') for s in periphery)))

# --- short cycles in the core, signed + classified ---
Gc = G.subgraph(core).copy()
cyc_pos = 0; cyc_neg = 0; sample = []
try:
    for cyc in nx.simple_cycles(Gc, length_bound=4):
        s = 1
        for a, b in zip(cyc, cyc[1:] + cyc[:1]):
            s *= G[a][b]['sign']
        if s > 0: cyc_pos += 1
        else: cyc_neg += 1
        if len(sample) < 40:
            sample.append(dict(cycle=cyc, sign='+' if s > 0 else '-',
                               mods=sorted(set(MOD_OF.get(x, '?') for x in cyc))))
except Exception as e:
    print("cycle enum note:", repr(e)[:100])
print(f"\nshort cycles (len<=4) in core: {cyc_pos} POSITIVE (switch/memory) + {cyc_neg} NEGATIVE (clock/homeostat)")

# --- verify hypothesized loops exist (edge presence + sign of a representative link) ---
def edge(a, b):
    return G[a][b]['sign'] if G.has_edge(a, b) else 0
known = {
    'L1 Rb-E2f (E2f<->Ce)': [('E2f', 'Ce'), ('Ce', 'pRb')],
    'L2 Skp2-p27 (Skp2 -| P21)': [('Skp2', 'P21'), ('P21', 'Ce')],
    'L4 Wee1/Cdc25-MPF (self)': [('MPF', 'MPF')],  # self via algebraic Cdc25a/Wee1a
    'N1 MPF-Cdc20': [('MPF', 'Cdc20'), ('Cdc20', 'MPF')],
    'N2 EZH2 -| CyclinD1': [('Cd', 'EZH2m'), ('EZH2', 'Cd_mRNA')],
    'N5 mass homeostat (mass->Cd? via commit)': [('mass', 'Rbm')],
}
print("\n=== hypothesized-loop edge check (sign of representative links) ===")
loopchk = {}
for name, edges in known.items():
    signs = []
    for a, b in edges:
        if a == b:
            signs.append(('self', selfreg.get(a, 0)))
        else:
            signs.append((f'{a}->{b}', edge(a, b)))
    loopchk[name] = signs
    print(f"  {name}: " + ", ".join(f'{e}={s:+d}' for e, s in signs))

# --- signed-Jacobian heatmap (ordered by module) ---
order = [s for m in MODULES for s in MODULES[m] if s in idx]
oi = [idx[s] for s in order]
Jm = np.mean(Js, axis=0)
H = np.sign(Jm[np.ix_(oi, oi)]) * np.log10(np.abs(Jm[np.ix_(oi, oi)]) + 1e-9)
fig, ax = plt.subplots(figsize=(13, 12))
vmax = np.percentile(np.abs(H[H != 0]), 95)
im = ax.imshow(H, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='equal')
ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=90, fontsize=6)
ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=6)
ax.set_title('Signed Jacobian d(dx_i/dt)/dx_j  (row i = affected, col j = source)\nred = activates, blue = inhibits; blocks = modules', fontsize=11)
# module boundaries
b = 0
for m in MODULES:
    k = sum(1 for s in MODULES[m] if s in idx)
    if k == 0: continue
    ax.add_patch(plt.Rectangle((b - .5, b - .5), k, k, fill=False, ec='k', lw=1.4))
    ax.text(b + k / 2 - .5, -1.5, m, rotation=45, fontsize=7, ha='left', va='bottom')
    b += k
plt.colorbar(im, fraction=0.046, pad=0.04, label='sign x log10|J|')
plt.tight_layout()
fig.savefig(os.path.join(OUT, 'A_signed_jacobian.png'), dpi=140, bbox_inches='tight')
print("wrote A_signed_jacobian.png")

json.dump(dict(n_species=n, n_edges=n_edges, n_signflip=int(np.sum(flip)),
               core=sorted(core), core_modules=sorted(set(MOD_OF.get(s, '?') for s in core)),
               periphery=sorted(periphery), periphery_modules=sorted(set(MOD_OF.get(s, '?') for s in periphery)),
               cycles_pos=cyc_pos, cycles_neg=cyc_neg, cycle_sample=sample, loop_checks=loopchk,
               selfreg={k: v for k, v in selfreg.items() if v != 0}),
          open(os.path.join(OUT, 'A_loops.json'), 'w'), indent=2)
print("wrote A_loops.json")
