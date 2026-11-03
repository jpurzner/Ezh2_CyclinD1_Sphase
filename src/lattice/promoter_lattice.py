"""Nucleosome-resolved promoter model (H3K27 lattice) — ground-truth simulator.

Built OUTSIDE Antimony on purpose: Antimony has no arrays, cannot express
per-nucleosome random replication partitioning, and v44 cannot be restarted
mid-arrest (see docs/stochastic_amplifier_2026-10-22.md). This module is the
reference implementation; an Antimony-expressible surrogate is fitted to it
separately (see fit_lattice_surrogate.py).

STATE. One lattice of `n` nucleosomes at `pitch` bp spacing, asymmetric about
the TSS. Each nucleosome carries ONE H3K27 state, mutually exclusive, following
Newar/Jost (PLoS Comput Biol 2022):

    0 = u (unmodified)   1 = ac   2 = me1   3 = me2   4 = me3

K4me3 is a static per-position input (core only), not a dynamic species.

ANTAGONISM. PRC2 is the only enzyme with read-write feedback; every antagonist
is local and non-spreading, which is the Newar/Jost result:
  * p300 -> ac in the core, gated by promoter activity
  * KDM6B/JMJD3 -> me3->me2->me1, rate k_kdm6b*Gli1*psi_kdm[i]   (Gli1-driven,
    with a MEASURED spatial profile -- GNP_P7_KDM6B ChIP)
  * transcription-coupled H3.3 exchange -> u over the gene body
  * local active marks (ac + K4me3) additionally inhibit PRC2 catalysis
    (H3K4me3 occupies the EED allosteric pocket; Finogenova 2020, NSMB 2024)

SPREADING is local, per Hodges & Crabtree (PNAS 2012): a geometric kernel with
nearest-neighbour probability `p_local`. Their result is that bounded domains
do not form at all when marking is spatially uniform, and that the controlling
number is kappa = k_prop/k_turn, optimal in [1.0, 1.5].

REPLICATION. At S phase each nucleosome is INDEPENDENTLY reset to u with
probability `dil_frac` (new histone). This is the discrete truth that
`Mk -> 0.5*Mk` approximates, and it supplies partition noise for free.

Time unit: HOURS (Hodges & Crabtree fit k_prop, k_turn ~ 0.1-0.2 /h; Newar
k_ac = 1.03 /h, gamma_ac = 0.6 /h, gamma_me = 1.5 /h).
"""
from dataclasses import dataclass, field
import numpy as np

U, AC, ME1, ME2, ME3 = 0, 1, 2, 3, 4
NSTATES = 5


# --------------------------------------------------------------------------- geometry
@dataclass
class Geometry:
    """Promoter geometry in bp, converted to a nucleosome lattice.

    Defaults are the transcriptional-model fit: lambda_up ~ 248 bp,
    lambda_dn ~ 4.7 kb -- strongly asymmetric, so the lattice is NOT centred.
    """
    pitch: float = 190.0          # bp per nucleosome
    up_bp: float = 1500.0         # lattice extent upstream of TSS
    dn_bp: float = 8000.0         # lattice extent downstream
    core_up_bp: float = 400.0     # H3K27me3-depleted core, upstream edge
    core_dn_bp: float = 700.0     # depleted core, downstream edge
    cgi_sigma_bp: float = 300.0   # CpG-island nucleation width
    body_start_bp: float = 300.0  # gene body (H3.3 exchange) starts here

    def build(self):
        x = np.arange(-self.up_bp, self.dn_bp + self.pitch, self.pitch)
        n = len(x)
        in_core = (x >= -self.core_up_bp) & (x <= self.core_dn_bp)
        nucleation = np.exp(-0.5 * (x / self.cgi_sigma_bp) ** 2)
        body = 1.0 / (1.0 + np.exp(-(x - self.body_start_bp) / 200.0))
        return x, n, in_core, nucleation, body


# --------------------------------------------------------------------------- parameters
@dataclass
class Params:
    # --- PRC2 writing (per hour) --------------------------------------------
    k_me1: float = 1.2            # u -> me1   (fast)
    k_me2: float = 0.8            # me1 -> me2 (fast)
    k_me3: float = 0.25           # me2 -> me3 (RATE LIMITING, Newar/Jost)
    a0_prc2: float = 0.10         # genome-wide basal PRC2 floor (everywhere)
    a_rw_prc2: float = 1.35       # EED read-write boost from local me3
    n_rw: float = 2.0             # read-write cooperativity
    p_local: float = 0.75         # geometric kernel: P(nearest neighbour)
    kernel: str = 'local'         # 'local' (geometric, H&C) | 'contact' (1/|i-j|^a, 3D)
    contact_exp: float = 1.0      # exponent for the 'contact' kernel (Jost uses 1)
    nucleation_gain: float = 0.55 # CpG-island nucleation, ADDITIVE on top of a0

    # --- erasers / antagonists (per hour) -----------------------------------
    del_mk: float = 0.06          # basal demethylation (all me steps)
    k_kdm6b: float = 0.50         # KDM6B/JMJD3 max rate (x Gli1 x psi_kdm)
    k_ac: float = 1.03            # p300 u -> ac            (Newar fit)
    g_ac: float = 0.60            # HDAC ac -> u, uniform   (Newar fit)
    k_h33: float = 0.15           # transcription-coupled H3.3 exchange -> u
    w_inh_active: float = 0.85    # local ac/K4me3 inhibition of PRC2 catalysis

    # --- replication ---------------------------------------------------------
    dil_frac: float = 0.5         # P(nucleosome replaced by naive histone) at S

    # --- spatial profiles of the antagonists --------------------------------
    p300_core_only: bool = True   # p300 restricted to the depleted core
    kdm_core_width: float = 1200.0  # psi_KDM6B gaussian sd (bp) about the TSS


@dataclass
class Inputs:
    """One-way coupling: trajectories recorded from the v44 cell-cycle model."""
    EZH2: float = 1.0        # PRC2 catalytic scale
    Gli1: float = 1.0        # drives KDM6B
    activity: float = 1.0    # promoter activity -> p300 acetylation
    tx: float = 1.0          # nascent output -> H3.3 exchange
    EZH2i: float = 0.0       # EZH2 inhibitor occupancy [0,1]


def kappa_eff(par: Params, inp: Inputs, local_me3: float = 0.6, body: float = 1.0):
    """Hodges & Crabtree's controlling ratio, evaluated in the FLANK.

        kappa = (me3 propagation rate) / (me3 turnover rate)

    They find bounded domains only for kappa in roughly [1.0, 1.5]; below that
    the domain collapses, above it the mark fills the locus through a
    second-order phase transition. Evaluating it here makes that constraint
    explicit instead of an emergent accident of the parameter set.
    """
    P = inp.EZH2 * (1 - inp.EZH2i) * (par.a0_prc2 + par.a_rw_prc2 * local_me3 ** par.n_rw)
    gain = P * par.k_me3
    loss = par.del_mk + par.k_h33 * inp.tx * body
    return gain / loss


def a_rw_for_kappa(kappa: float, par: Params, inp: Inputs,
                   local_me3: float = 0.6, body: float = 1.0):
    """Solve a_rw_prc2 so the flank sits at a target kappa."""
    loss = par.del_mk + par.k_h33 * inp.tx * body
    P_need = kappa * loss / par.k_me3
    base = inp.EZH2 * (1 - inp.EZH2i)
    return (P_need / base - par.a0_prc2) / local_me3 ** par.n_rw


# --------------------------------------------------------------------------- model
class PromoterLattice:
    def __init__(self, geom: Geometry = None, par: Params = None, seed: int = 0,
                 init_me3: float = 0.0):
        self.geom = geom or Geometry()
        self.par = par or Params()
        self.x, self.n, self.in_core, self.nucl, self.body = self.geom.build()
        self.rng = np.random.default_rng(seed)
        p = self.par
        # geometric spreading kernel, normalised, local-dominant
        d = np.abs(np.subtract.outer(np.arange(self.n), np.arange(self.n))).astype(float)
        with np.errstate(divide='ignore', invalid='ignore'):
            if p.kernel == 'local':
                K = p.p_local * (1 - p.p_local) ** (d - 1)
            elif p.kernel == 'contact':
                K = 1.0 / np.power(d, p.contact_exp)      # polymer contact prob
            else:
                raise ValueError(f'unknown kernel {p.kernel!r}')
        K[d == 0] = 0.0
        self.K = K / K.sum(1, keepdims=True)
        # antagonist spatial profiles
        self.psi_p300 = self.in_core.astype(float) if p.p300_core_only else np.ones(self.n)
        self.psi_kdm = np.exp(-0.5 * (self.x / p.kdm_core_width) ** 2)
        self.k4 = self.in_core.astype(float)      # static K4me3 in the core
        self.state = np.zeros(self.n, dtype=np.int8)
        if init_me3 > 0:
            self.state[self.rng.random(self.n) < init_me3] = ME3

    # -------------------------------------------------------------- helpers
    def me3(self):
        return self.state == ME3

    def prc2_activity(self, inp: Inputs):
        """Per-site PRC2 catalytic rate: nucleation + local read-write, minus
        inhibition by local active marks."""
        p = self.par
        m3 = self.me3().astype(float)
        local_me3 = self.K @ m3
        active = np.clip(self.k4 + (self.state == AC).astype(float), 0, 1)
        local_active = self.K @ active
        boost = p.a0_prc2 + p.nucleation_gain * self.nucl + p.a_rw_prc2 * local_me3 ** p.n_rw
        return inp.EZH2 * (1 - inp.EZH2i) * boost * (1 - p.w_inh_active * local_active)

    def rate_matrix(self, inp: Inputs):
        """Per-site transition rates, shape (n, NSTATES, NSTATES)."""
        p, n = self.par, self.n
        P = np.clip(self.prc2_activity(inp), 0, None)
        erase = p.del_mk + p.k_kdm6b * inp.Gli1 * self.psi_kdm
        exch = p.k_h33 * inp.tx * self.body
        R = np.zeros((n, NSTATES, NSTATES))
        R[:, U,   ME1] = P * p.k_me1
        R[:, ME1, ME2] = P * p.k_me2
        R[:, ME2, ME3] = P * p.k_me3
        R[:, ME3, ME2] = erase
        R[:, ME2, ME1] = erase
        R[:, ME1, U]   = p.del_mk
        R[:, U,   AC]  = p.k_ac * inp.activity * self.psi_p300
        R[:, AC,  U]   = p.g_ac
        for s in (AC, ME1, ME2, ME3):            # H3.3 exchange resets to u
            R[:, s, U] += exch
        return R

    # -------------------------------------------------------------- stepping
    MAX_P = 0.15          # cap on per-substep transition probability

    def step(self, dt: float, inp: Inputs):
        """One tau-leap step, auto-subdivided so no substep exceeds MAX_P."""
        R = self.rate_matrix(inp)
        pmax = (R[np.arange(self.n), self.state].sum(1) * dt).max()
        nsub = int(np.ceil(pmax / self.MAX_P)) if pmax > self.MAX_P else 1
        for _ in range(nsub):
            self._substep(dt / nsub, inp)

    def _substep(self, dt: float, inp: Inputs):
        R = self.rate_matrix(inp)
        cur = R[np.arange(self.n), self.state]            # (n, NSTATES)
        ptrans = cur * dt
        tot = ptrans.sum(1)
        u = self.rng.random(self.n)
        fire = u < tot
        if fire.any():
            idx = np.where(fire)[0]
            pr = ptrans[idx] / tot[idx, None]
            cdf = np.cumsum(pr, 1)
            draws = self.rng.random(len(idx))
            self.state[idx] = (draws[:, None] > cdf).sum(1).astype(np.int8)

    def replicate(self):
        """S phase: each nucleosome independently replaced by a naive histone."""
        m = self.rng.random(self.n) < self.par.dil_frac
        self.state[m] = U

    def run(self, hours: float, dt: float, inp: Inputs,
            cycle_h: float = None, record_every: float = None):
        """Integrate, optionally with a replication event every `cycle_h`."""
        nsteps = int(round(hours / dt))
        rec_n = max(1, int(round((record_every or dt) / dt)))
        t, prof, core, times = 0.0, [], [], []
        next_S = cycle_h if cycle_h else np.inf
        for i in range(nsteps):
            self.step(dt, inp)
            t += dt
            if t >= next_S:
                self.replicate()
                next_S += cycle_h
            if i % rec_n == 0:
                prof.append(self.me3().astype(float))
                core.append(float(self.me3()[self.in_core].mean())
                            if self.in_core.any() else np.nan)
                times.append(t)
        return dict(t=np.array(times), profile=np.array(prof),
                    core_me3=np.array(core), x=self.x)


def ensemble(n_cells=40, hours=400.0, dt=0.05, inp=None, geom=None, par=None,
             cycle_h=None, seed0=0, burn_frac=0.5, init_me3=0.0):
    """Mean me3 profile over independent cells — the ChIP-comparable observable."""
    inp = inp or Inputs()
    profs, cores = [], []
    for c in range(n_cells):
        lat = PromoterLattice(geom, par, seed=seed0 + c, init_me3=init_me3)
        out = lat.run(hours, dt, inp, cycle_h=cycle_h, record_every=1.0)
        k = int(len(out['t']) * burn_frac)
        profs.append(out['profile'][k:].mean(0))
        cores.append(np.nanmean(out['core_me3'][k:]))
    lat0 = PromoterLattice(geom, par, seed=seed0, init_me3=init_me3)
    return dict(x=lat0.x, in_core=lat0.in_core,
                profile=np.mean(profs, 0), profile_sd=np.std(profs, 0),
                core_me3=float(np.nanmean(cores)), core_me3_sd=float(np.nanstd(cores)))
