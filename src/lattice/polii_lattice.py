"""Lattice with explicit RNA Pol II, mutually antagonistic with PRC2.

What this adds over `gene_lattice`
----------------------------------
Previously transcription entered only as a constant `tx` scaling a fixed,
NDR-shaped exchange field. Histone turnover was therefore an input, not a
consequence, and the antagonism ran one way: the active marks suppressed PRC2
catalysis, but nothing suppressed transcription in return.

Here Pol II is explicit:

  LOADING    Poisson at the promoter, with the rate repressed by promoter me3

                 k_load_eff = k_load / (1 + (m_prom / Km)^h)

  ELONGATION each polymerase advances v_elong sites per second

  EVICTION   each site a polymerase crosses loses its mark with probability
             p_evict -- so exchange is now produced BY transcription

That closes the loop: me3 --| Pol II --| (eviction) --| me3. A double negative
is positive feedback on me3, and combined with the existing active --| PRC2
term it is a mutual-antagonism toggle.

Why that matters: read-write alone cannot be bistable, because the writing
flux W(0) = beta > 0 means a straight loss line crosses it exactly once (see
docs/lattice_model_2026-10-25.pdf). Pol II-driven loss is NOT a straight line
-- it is steep while Pol II is loading and flattens once me3 shuts Pol II off
-- so W can cross it three times. Mean field says this needs a Hill
coefficient h >= 2 and a Pol II-driven loss of order 10x the basal writing
rate; h = 1 is monostable at every coupling strength tested.
"""
import numpy as np
from dataclasses import dataclass, replace

from gene_lattice import GeneLattice, GeneGeom
from particle_lattice import P, U, ME1, ME2, ME3, AC, HOUR


@dataclass
class PolII:
    k_load: float = 0.0020      # /s attempts at the promoter when fully de-repressed
    Km_me3: float = 0.25        # promoter me3 giving half-maximal repression
    h_me3: float = 3.0          # cooperativity of PRC2 repression of loading
    v_elong: float = 0.010      # sites/s  (~2 kb/min at 200 bp sites)
    p_evict: float = 0.55       # chance a passing polymerase evicts the mark
    prom_half: int = 5          # sites either side of the TSS that count as promoter
    basal_exch: float = 0.02 / HOUR   # Pol II-independent turnover, everywhere


class PolIILattice(GeneLattice):
    def __init__(self, geom, par=None, pol=None, seed=0, init_me3=0.8,
                 active_scale=1.0, tss_site=None):
        super().__init__(geom, par, seed=seed, init_me3=init_me3,
                         active_scale=active_scale)
        self.pol = pol or PolII()
        self.tss = int(np.argmin(np.abs(geom.x))) if tss_site is None else tss_site
        lo = max(0, self.tss - self.pol.prom_half)
        hi = min(self.n, self.tss + self.pol.prom_half + 1)
        self.prom = slice(lo, hi)
        self.pols = []              # positions (float site index) of elongating Pol II
        self.loaded = 0             # cumulative, for reporting
        self._occ = np.zeros(self.n)   # time-averaged Pol II occupancy
        self.t = 0.0
        # KDM6B is recruited by Gli1 (they co-IP), so the eraser goes where
        # Gli1 is bound -- the proximal promoter and the enhancers -- rather
        # than being spread over the core as in the parent class.
        g = geom.gli1 if getattr(geom, "gli1", None) is not None else None
        self.gli1 = np.zeros(self.n) if g is None else np.asarray(g, float)
        self.gli1_on = 1.0          # set to 0 to withdraw the mitogen
        self._set_kdm6_weights()
        self.init_times = []        # when each polymerase initiated (s)
        self.prom_trace = []        # (t, promoter me3) sampled each step

    # ------------------------------------------------------------------
    def _set_kdm6_weights(self):
        w = self.gli1 * self.gli1_on
        if w.sum() <= 0:                      # Gli1 withdrawn: no targeting
            self.w_kdm6 = None
        else:
            self.w_kdm6 = w / w.sum()

    def withdraw_gli1(self):
        """Mitogen removal: Gli1 leaves, so KDM6B is no longer recruited."""
        self.gli1_on = 0.0
        self._set_kdm6_weights()

    def promoter_me3(self):
        return float((self.state[self.prom] == ME3).mean())

    def load_rate(self):
        p = self.pol
        return p.k_load / (1.0 + (self.promoter_me3() / p.Km_me3) ** p.h_me3)

    def step(self, dt):
        p = self.pol
        # The parent fires KDM6B on `self.w_core`. Point that at the Gli1 sites
        # for the duration of the parent call, and disable it when Gli1 is gone.
        _saved_core, _saved_n = self.w_core, self.p.n_kdm6
        if self.w_kdm6 is None:
            self.p = replace(self.p, n_kdm6=0.0)
        else:
            self.w_core = self.w_kdm6
        # --- PRC2 / p300 / KDM6 / HDAC / basal demethylation, unchanged ----
        # (the parent also applies its NDR exchange term; we neutralise that
        #  by zeroing geom.exch in the caller so turnover comes from Pol II)
        super().step(dt)
        self.w_core = _saved_core
        if self.p.n_kdm6 != _saved_n:
            self.p = replace(self.p, n_kdm6=_saved_n)
        self.t += dt
        self.prom_trace.append(self.promoter_me3())

        # --- Pol II loading, repressed by promoter me3 ---------------------
        n_new = self.rng.poisson(self.load_rate() * dt)
        for _ in range(n_new):
            self.pols.append(float(self.tss))
            # spread initiations uniformly within the step so inter-event
            # intervals are not quantised to dt, which would fake burstiness
            self.init_times.append(self.t + float(self.rng.random()) * dt)
        self.loaded += n_new

        # --- elongation and mark eviction ----------------------------------
        if self.pols:
            adv = p.v_elong * dt
            keep = []
            for q in self.pols:
                a, b = q, q + adv
                i0, i1 = int(np.floor(a)), int(np.floor(b))
                for i in range(max(i0, 0), min(i1 + 1, self.n)):
                    if i >= i0 and self.rng.random() < p.p_evict:
                        self.state[i] = U
                if b < self.n:
                    keep.append(b)
                    # credit every site traversed this step, in proportion to
                    # the time spent there -- crediting only the END site
                    # biases occupancy by the phase of the step.
                    lo_s, hi_s = int(np.floor(a)), int(np.floor(b))
                    if hi_s == lo_s:
                        self._occ[min(lo_s, self.n - 1)] += dt
                    else:
                        span = b - a
                        for i in range(max(lo_s, 0), min(hi_s + 1, self.n)):
                            seg = min(b, i + 1) - max(a, i)
                            if seg > 0:
                                self._occ[i] += dt * seg / span
            self.pols = keep

        # --- Pol II-independent turnover, flat across the lattice ----------
        r = self.rng.random(self.n)
        self.state[(self.state != U) & (r < p.basal_exch * dt)] = U

    def run(self, hours, dt=30.0, cycle_h=None, record_every_h=1.0):
        out = super().run(hours, dt, cycle_h=cycle_h, record_every_h=record_every_h)
        out["pol_occ"] = self._occ / max(hours * HOUR, 1e-9)
        out["loaded"] = self.loaded
        out["init_times"] = np.array(self.init_times)
        out["prom_trace"] = np.array(self.prom_trace)
        return out


def population(geom, n_cells=60, hours=300, dt=30.0, par=None, pol=None,
               cycle_h=24.0, init_me3=0.8, burn_frac=0.5, seed0=0,
               prc2_cv=0.4, active_cv=0.25, record_every_h=4.0):
    """Ensemble with Pol II. `geom.exch` is zeroed so that ALL turnover comes
    from transcription plus the flat basal term -- otherwise the NDR field
    would still dominate and Pol II would be decorative."""
    g = GeneGeom(**{**geom.__dict__, "exch": np.zeros_like(geom.exch)})
    base = par or P()
    pbase = pol or PolII()
    rng = np.random.default_rng(seed0)
    finals, occs, loads = [], [], []
    for c in range(n_cells):
        npr = base.n_prc2 * float(np.exp(rng.normal(0, prc2_cv) - prc2_cv**2 / 2))
        asc = float(np.clip(rng.normal(1.0, active_cv), 0.0, None))
        lat = PolIILattice(g, replace(base, n_prc2=npr), pbase,
                           seed=seed0 + 2000 + c, init_me3=init_me3,
                           active_scale=asc)
        o = lat.run(hours, dt, cycle_h=cycle_h, record_every_h=record_every_h)
        prof = o["profile"]
        finals.append(prof[int(len(prof) * burn_frac):].mean(0))
        occs.append(o["pol_occ"]); loads.append(o["loaded"])
    finals = np.array(finals); occs = np.array(occs)
    return dict(x=g.x, per_cell=finals, mean=finals.mean(0), sd=finals.std(0),
                pol_occ=occs.mean(0), pol_occ_per_cell=occs,
                loaded=np.array(loads), target=g.target_me3, name=g.name)
