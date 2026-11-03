"""Per-gene lattice: geometry taken from the measured ChIP profile, simulated
as a POPULATION of cells, averaged to compare against bulk ChIP.

The point of the population is that bulk ChIP cannot distinguish "every cell
has an intermediate level here" from "some cells are marked and some are not".
At these promoters H3K27me3 and H3K4me3 overlap in bulk, and the overlap may be
true bivalency in each cell or a mixture of two populations. Only a simulated
ensemble can produce the second, so the comparison to bulk has to be made on
the ENSEMBLE MEAN, not on one trajectory.

What is data and what is model
------------------------------
INPUT  (measured, per site): H3K4me3 and H3K27ac -> the "active" occupancy that
       antagonises PRC2; the NDR position -> where histone exchange peaks.
OUTPUT (predicted): the H3K27me3 profile.
So this is not a fit of me3 to itself. The non-trivial content is whether the
measured active profile, plus PRC2 kinetics taken from the literature, puts the
me3 boundary in the right PLACE with the right SHARPNESS and AMPLITUDE -- and
whether single-cell boundaries have to be sharp or graded to average correctly.
"""
import numpy as np
from dataclasses import dataclass, field, replace

from particle_lattice import ParticleLattice, P, U, ME1, ME2, ME3, AC, HOUR


@dataclass
class GeneGeom:
    """Lattice geometry read off a gene's binned ChIP profile."""
    name: str
    x: np.ndarray            # site centre, bp from TSS, oriented 5'->3'
    active: np.ndarray       # 0-1, measured K4me3/K27ac occupancy
    exch: np.ndarray         # 0-1, histone exchange weight (NDR-peaked)
    nucl: np.ndarray         # 0-1, PRC2 nucleation weight
    target_me3: np.ndarray   # measured me3, same sites (the thing to predict)
    pitch: float = 200.0
    gli1: np.ndarray = None  # 0-1, measured Gli1 occupancy -> where KDM6B goes

    @property
    def n(self):
        return len(self.x)

    @classmethod
    def from_bins(cls, path, tss, strand, pitch=200.0, half=10000.0,
                  active_marks=("H3K4me3", "H3K27ac"), nucleation="uniform",
                  nucl_sites=(), nucl_sigma=600.0, exch_sigma=450.0):
        """`nucleation` selects where PRC2 is recruited, which is NOT K4me3 --
        K4me3 is the antagonist, and seeding PRC2 there inverts the geometry.
        Options: "uniform" (no bias), "sites" (Gaussian bumps at `nucl_sites`,
        bp from TSS, 5'->3' oriented -- use neighbouring promoters), or
        "antisite" (bias AWAY from the active zone).
        """
        rows = [l.rstrip("\n").split("\t") for l in open(path)]
        hdr, dat = rows[0], rows[1:]
        col = {c: hdr.index(c) for c in hdr}
        mid = np.array([(int(r[1]) + int(r[2])) / 2 for r in dat]) - tss
        if strand == "-":
            mid = -mid
        order = np.argsort(mid)
        mid = mid[order]
        get = lambda c: np.array([float(r[col[c]]) for r in dat])[order]

        edges = np.arange(-half, half + pitch, pitch)
        ctr = 0.5 * (edges[:-1] + edges[1:])
        idx = np.digitize(mid, edges) - 1
        ok = (idx >= 0) & (idx < len(ctr))

        def agg(v):
            out = np.zeros(len(ctr)); cnt = np.zeros(len(ctr))
            np.add.at(out, idx[ok], v[ok]); np.add.at(cnt, idx[ok], 1)
            return out / np.clip(cnt, 1, None)

        def unit(v):
            lo, hi = np.percentile(v, 1), np.percentile(v, 99)
            return np.clip((v - lo) / (hi - lo), 0, 1) if hi > lo else v * 0

        act = np.max(np.vstack([unit(agg(get(m))) for m in active_marks]), axis=0)
        if nucleation == "sites" and len(nucl_sites):
            nuc = np.zeros(len(ctr))
            for s0 in nucl_sites:
                nuc = np.maximum(nuc, np.exp(-0.5 * ((ctr - s0) / nucl_sigma) ** 2))
        elif nucleation == "antisite":
            nuc = np.clip(1.0 - act, 0, 1)
        else:
            nuc = np.ones(len(ctr))
        me3 = agg(get("H3K27me3"))
        exch = np.exp(-0.5 * ((ctr - 150.0) / exch_sigma) ** 2)
        gli = unit(agg(get("Gli1"))) if "Gli1" in col else np.zeros(len(ctr))
        return cls(name=path.split("/")[-1].split("_")[0], x=ctr, active=act,
                   exch=exch / exch.max(), nucl=nuc, target_me3=me3, pitch=pitch,
                   gli1=gli)


class GeneLattice(ParticleLattice):
    """ParticleLattice whose geometry comes from a gene's measured profile."""

    def __init__(self, geom: GeneGeom, par=None, seed=0, init_me3=0.8,
                 active_scale=1.0):
        self.g = geom
        self.p = par or P()
        self.x = geom.x
        self.n = geom.n
        self.nucl = geom.nucl
        self.exch = geom.exch
        self.in_core = geom.active > 0.5
        self.rng = np.random.default_rng(seed)
        self.state = np.zeros(self.n, dtype=np.int8)
        if init_me3 > 0:
            self.state[self.rng.random(self.n) < init_me3] = ME3
        # the measured active profile IS the K4me3 field, scaled per cell
        self.k4 = np.clip(geom.active * active_scale, 0, 1)

        d = np.abs(np.subtract.outer(np.arange(self.n), np.arange(self.n))).astype(float)
        K = np.exp(-(d - 1) / self.p.allo_range_nuc)
        K[d == 0] = 0.0
        K[d > 6] = 0.0
        self.K = K / np.clip(K.sum(1, keepdims=True), 1e-9, None)
        w = 1.0 + self.p.nucleation_gain * self.nucl
        self.w_prc2 = w / w.sum()
        wc = self.in_core.astype(float) + 1e-6
        self.w_core = wc / wc.sum()


def population(geom, n_cells=200, hours=600, dt=30.0, par=None, cycle_h=24.0,
               init_me3=0.8, burn_frac=0.5, seed0=0,
               prc2_cv=0.4, active_cv=0.25, record_every_h=4.0):
    """Simulate `n_cells` independent cells and return per-cell and mean profiles.

    Cell-to-cell variety comes from three places: the stochastic dynamics, a
    lognormal spread in PRC2 abundance (`prc2_cv`), and a spread in how strongly
    the active machinery is engaged (`active_cv`). Without the latter two every
    cell converges to the same boundary and the ensemble mean is just a single
    trajectory with less noise -- which cannot reproduce a graded bulk profile.
    """
    base = par or P()
    rng = np.random.default_rng(seed0)
    finals, traces = [], []
    for c in range(n_cells):
        npr = base.n_prc2 * float(np.exp(rng.normal(0, prc2_cv) - prc2_cv**2 / 2))
        asc = float(np.clip(rng.normal(1.0, active_cv), 0.0, None))
        lat = GeneLattice(geom, replace(base, n_prc2=npr), seed=seed0 + 1000 + c,
                          init_me3=init_me3, active_scale=asc)
        o = lat.run(hours, dt, cycle_h=cycle_h, record_every_h=record_every_h)
        prof = o["profile"]
        k0 = int(len(prof) * burn_frac)
        finals.append(prof[k0:].mean(0))        # this cell's time-averaged profile
        traces.append(prof)
    finals = np.array(finals)
    return dict(x=geom.x, per_cell=finals, mean=finals.mean(0),
                sd=finals.std(0), t=o["t"], traces=np.array(traces),
                target=geom.target_me3, name=geom.name)
