"""Particle-based H3K27 promoter lattice.

Replaces the continuous "PRC2 activity field" of promoter_lattice.py with
DISCRETE ENZYME MOLECULES that bind, dwell, catalyse and unbind. Written against
the six-part literature census (docs/lit_*.md, docs/parameters_consolidated_2026-10-23.md).

WHY PARTICLES. Measured abundance puts every writer and eraser at ~1 or fewer
molecules per target promoter domain:
    PRC2   ~1-4 per domain        p300  ~1-7 per domain
    KDM6B  ~0.3 per domain        HDAC  1 per 18-59 NUCLEOSOMES
So PRC2, KDM6 and p300 are modelled as particles. HDAC is the one enzyme dense
enough to remain a field, and is kept as one.

BINDING (all grade A, Tatavosian 2018 mESC Suppl. Table 2):
    tau_search    210 s     tau_stable 10.0 s     tau_transient 0.9 s
    total bound fraction 23%, STABLE bound fraction 4.6%
  -> 4.6/23 = 20% of binding events enter the stable mode.
  Internal check: 10/(10+210) = 4.5% duty cycle vs 4.6% measured. Independent.

CATALYSIS DURING BOTH MODES. The census showed that if only the stable 4.6%
catalyses, the required kcat (0.6-2.1 /min) exceeds every measured nucleosomal
value (max 0.474 /min). Allowing the 0.9 s transient population to catalyse
brings the requirement to 0.12-0.42 /min, inside the measured range. That is the
central modelling claim here and it is falsifiable.

READ-WRITE ACTS ON CATALYSIS, NOT BINDING. H3K27me3 changes kcat and leaves K_M
and binding affinity unaffected (Margueron 2009; null result for binding across
Wang 2017 / Cookis 2025 / Finogenova 2020). Allosteric stimulation 7-12x, and
measured to act over 1-4 nucleosomes in cis with no trans component
(Ge/Muir JACS 2019: 8.8x proximal, 3.2x distal).

NUCLEATION IS AN ACCESSORY-SUBUNIT PROPERTY. PRC2 core shows no positional
preference for a CpG island at 2 nm accuracy (Heenan 2020), so the nucleation
term is attributed to PCL/MTF2, and enters as a bias on WHERE binding happens
and on the probability of entering the stable mode.

TURNOVER is the measured 0.0044 /h (Zheng 2014, me3 t1/2 158 h), not the 0.12 /h
the field model used. In vivo loss is dominated by replicative dilution
(Reveron-Gomez 2018; Jadhav 2020; Alabert 2020 rejects demethylation by BIC).

Time unit: SECONDS internally (binding is a seconds-scale process), reported in hours.
"""
from dataclasses import dataclass, field
import numpy as np

U, AC, ME1, ME2, ME3 = 0, 1, 2, 3, 4
HOUR = 3600.0


@dataclass
class Geom:
    pitch: float = 190.0
    up_bp: float = 3000.0
    dn_bp: float = 8000.0
    core_up_bp: float = 400.0
    core_dn_bp: float = 700.0
    cgi_sigma_bp: float = 300.0
    exch_centre_bp: float = 150.0    # NDR sits just downstream of the TSS
    exch_sigma_bp: float = 450.0     # histone turnover is PROMOTER-LOCAL

    def build(self):
        x = np.arange(-self.up_bp, self.dn_bp + self.pitch, self.pitch)
        in_core = (x >= -self.core_up_bp) & (x <= self.core_dn_bp)
        nucl = np.exp(-0.5 * (x / self.cgi_sigma_bp) ** 2)
        # Exchange is a LOCAL peak at the nucleosome-depleted region, not a
        # gene-body sigmoid. Kraushaar 2013 / Deaton 2016: promoters are the
        # highest-turnover class; PcG domains are among the LOWEST. The earlier
        # 8 kb body profile demanded 18 methylation steps/h from a domain whose
        # measured enzyme complement can supply ~1.4.
        exch = np.exp(-0.5 * ((x - self.exch_centre_bp) / self.exch_sigma_bp) ** 2)
        return x, len(x), in_core, nucl, exch


@dataclass
class P:
    # --- PRC2 binding (Tatavosian 2018, grade A) --------------------------
    tau_search: float = 210.0       # s, unbound search time
    tau_stable: float = 10.0        # s
    tau_transient: float = 0.9      # s
    p_stable: float = 0.20          # fraction of events entering the stable mode
    # COMMITTED DEFAULT (2026-10-23). targets = 10,000 measured from H3K27me3
    # epic2 peaks (TSS with a peak within +/-5 kb; EL_GNP1 10,069 / EL_Tumor1
    # 11,537 -- the two samples that pass histone QC). PRC2 = 152,334 (Itzhak).
    # 152,334 / 10,000 = 15.2. See the note on why the HIGH abundance estimate.
    n_prc2: float = 15.0            # PRC2 molecules serving this domain

    # --- PRC2 catalysis ----------------------------------------------------
    # COMMITTED: 0.47 /min, the fastest measured nucleosomal value (PALI1 PIR,
    # 100 mM KCl). Chosen because the model is not viable below ~0.3; this is an
    # assumption that the model is right, not an independent measurement.
    kcat: float = 0.47              # per MINUTE while bound (measured 0.001-0.474)
    allo_max: float = 9.0           # allosteric stimulation, 7-12x (Margueron; Lee)
    allo_range_nuc: float = 2.5     # 1-4 nucleosomes (Ge/Muir, in cis, no trans)
    # In vivo ratio from Zheng 2014 Table 1B (me0->me1 3.58, me1->me2 0.41,
    # me2->me3 0.54 x10^-2 /h) = 8.7 : 1.0 : 1.3. The 10:4:1 peptide ratio is
    # not measured step-resolved on nucleosomes anywhere.
    step_rel: tuple = (8.7, 1.0, 1.3)    # u->me1 : me1->me2 : me2->me3
    w_inh_active: float = 0.90      # K4me3 / K27ac block of CATALYSIS (not binding)
    nucleation_gain: float = 3.0    # PCL-mediated bias on where PRC2 binds

    # --- erasers / antagonists --------------------------------------------
    del_mk: float = 0.0044 / HOUR   # /s  measured me3 loss (Zheng t1/2 158 h)
    n_kdm6: float = 0.3             # KDM6 molecules per domain (859 copies/cell)
    kdm6_kcat: float = 0.70         # /min, peptide upper bound (Walport 2016)
    kdm6_tau: float = 10.0          # s, NO measurement exists - assumed PRC2-like
    n_p300: float = 2.0             # p300 molecules per domain
    p300_kcat: float = 3.1          # /min, H3K18 proxy (no H3K27 kcat exists)
    p300_tau: float = 26.0          # s, FRAP (Ferrie 2023)
    hdac_rate: float = 0.093 / 60   # /s  CoREST on nucleosomal H3K27ac (Wang 2020)
    k_exch: float = 0.60 / HOUR     # /s  peak H3.3 exchange AT THE NDR
    tx: float = 1.0                 # transcriptional activity, scales exchange

    dil_frac: float = 0.5


class ParticleLattice:
    def __init__(self, geom=None, par=None, seed=0, init_me3=0.8):
        self.g = geom or Geom()
        self.p = par or P()
        self.x, self.n, self.in_core, self.nucl, self.exch = self.g.build()
        self.rng = np.random.default_rng(seed)
        self.state = np.zeros(self.n, dtype=np.int8)
        if init_me3 > 0:
            self.state[self.rng.random(self.n) < init_me3] = ME3
        self.k4 = self.in_core.astype(float)
        # allosteric sensing kernel: exponential over allo_range_nuc nucleosomes
        d = np.abs(np.subtract.outer(np.arange(self.n), np.arange(self.n))).astype(float)
        K = np.exp(-(d - 1) / self.p.allo_range_nuc)
        K[d == 0] = 0.0
        K[d > 6] = 0.0                      # measured: no stimulation beyond ~4
        self.K = K / np.clip(K.sum(1, keepdims=True), 1e-9, None)
        # where enzymes bind: PCL-biased for PRC2, core-biased for p300/KDM6
        w = 1.0 + self.p.nucleation_gain * self.nucl
        self.w_prc2 = w / w.sum()
        wc = self.in_core.astype(float) + 1e-6
        self.w_core = wc / wc.sum()

    # ---------------------------------------------------------------- helpers
    def _allo(self):
        return 1.0 + (self.p.allo_max - 1.0) * (self.K @ (self.state == ME3).astype(float))

    def _inhib(self):
        act = np.clip(self.k4 + (self.state == AC), 0, 1)
        return np.clip(1.0 - self.p.w_inh_active * (self.K @ act), 0.0, 1.0)

    def _prc2_fire(self, i, dwell):
        """One bound PRC2 at site i for `dwell` seconds. Returns True if it wrote."""
        s = self.state[i]
        if s == AC or s == ME3:
            return False
        step = {U: 0, ME1: 1, ME2: 2}[s]
        # kcat refers to the RATE-LIMITING step. Published nucleosomal kcat is
        # an aggregate H3K27-methylation turnover number, so normalising to the
        # FASTEST step (as I first did) penalises the limiting step ~9x twice over.
        k = (self.p.kcat / 60.0) * self.p.step_rel[step] / min(self.p.step_rel)
        k *= self._allo()[i] * self._inhib()[i]
        if self.rng.random() < 1.0 - np.exp(-k * dwell):
            self.state[i] = {U: ME1, ME1: ME2, ME2: ME3}[s]
            return True
        return False

    def _eraser_fire(self, i, dwell, kcat_min):
        s = self.state[i]
        if s not in (ME1, ME2, ME3):
            return False
        if self.rng.random() < 1.0 - np.exp(-(kcat_min / 60.0) * dwell):
            self.state[i] = {ME3: ME2, ME2: ME1, ME1: U}[s]
            return True
        return False

    def _p300_fire(self, i, dwell):
        if self.state[i] != U:
            return False
        if self.rng.random() < 1.0 - np.exp(-(self.p.p300_kcat / 60.0) * dwell):
            self.state[i] = AC
            return True
        return False

    # ---------------------------------------------------------------- stepping
    def _enzyme_events(self, dt, n_mol, tau_bound, weights, fire, p_stable=None,
                       tau_transient=None):
        """Poisson-many binding events in dt, each at a weighted-random site."""
        if n_mol <= 0:
            return 0
        cycle = self.p.tau_search + tau_bound
        n_ev = self.rng.poisson(n_mol * dt / cycle)
        fired = 0
        for _ in range(n_ev):
            i = int(self.rng.choice(self.n, p=weights))
            if p_stable is None:
                d = self.rng.exponential(tau_bound)
            else:
                stable = self.rng.random() < p_stable * (0.3 + 0.7 * self.nucl[i])
                d = self.rng.exponential(self.p.tau_stable if stable else tau_transient)
            fired += bool(fire(i, d))
        return n_ev

    def step(self, dt):
        p = self.p
        # --- discrete enzymes ------------------------------------------------
        tau_mean = p.p_stable * p.tau_stable + (1 - p.p_stable) * p.tau_transient
        self._enzyme_events(dt, p.n_prc2, tau_mean, self.w_prc2,
                            self._prc2_fire, p_stable=p.p_stable,
                            tau_transient=p.tau_transient)
        self._enzyme_events(dt, p.n_kdm6, p.kdm6_tau, self.w_core,
                            lambda i, d: self._eraser_fire(i, d, p.kdm6_kcat))
        self._enzyme_events(dt, p.n_p300, p.p300_tau, self.w_core, self._p300_fire)
        # --- fields: HDAC, basal demethylation, transcription-coupled exchange
        r = self.rng.random(self.n)
        isac = self.state == AC
        self.state[isac & (r < p.hdac_rate * dt)] = U
        # NOTE: iterate LOW->HIGH so a nucleosome demoted this step is not
        # re-tested and cascaded. Forward order let a single draw take
        # me3 -> me2 -> me1 -> u in one step, tripling the effective loss rate.
        r = self.rng.random(self.n)
        for s, tgt in ((ME1, U), (ME2, ME1), (ME3, ME2)):
            m = (self.state == s) & (r < p.del_mk * dt)
            self.state[m] = tgt
        r = self.rng.random(self.n)
        self.state[(self.state != U) & (r < p.k_exch * p.tx * self.exch * dt)] = U

    def replicate(self):
        self.state[self.rng.random(self.n) < self.p.dil_frac] = U

    def run(self, hours, dt=30.0, cycle_h=None, record_every_h=1.0):
        nst = int(hours * HOUR / dt)
        rec = max(1, int(record_every_h * HOUR / dt))
        nxt = cycle_h * HOUR if cycle_h else np.inf
        t = 0.0
        prof, times = [], []
        for i in range(nst):
            self.step(dt); t += dt
            if t >= nxt:
                self.replicate(); nxt += cycle_h * HOUR
            if i % rec == 0:
                prof.append((self.state == ME3).astype(float)); times.append(t / HOUR)
        return dict(t=np.array(times), profile=np.array(prof), x=self.x)


def ensemble(n_cells=12, hours=1500, dt=30.0, geom=None, par=None, cycle_h=None,
             init_me3=0.8, burn=0.7, seed0=0):
    profs = []
    for c in range(n_cells):
        lat = ParticleLattice(geom, par, seed=seed0 + c, init_me3=init_me3)
        o = lat.run(hours, dt, cycle_h=cycle_h)
        profs.append(o['profile'][int(len(o['t']) * burn):].mean(0))
    l0 = ParticleLattice(geom, par, seed=seed0, init_me3=init_me3)
    pf = np.mean(profs, 0)
    return dict(x=l0.x, in_core=l0.in_core, profile=pf, sd=np.std(profs, 0))


def metrics(x, pf, in_core):
    core = float(pf[in_core].mean())
    plat = float(np.percentile(pf[x > 2000], 75))
    out = dict(core=core, plateau=plat, depletion=plat - core)
    if plat - core > 0.05:
        half = 0.5 * (core + plat)
        xs, ps = x[x > 0], pf[x > 0]
        ab = np.where(ps > half)[0]
        out['edge_bp'] = float(xs[ab[0]]) if len(ab) else float('nan')
    else:
        out['edge_bp'] = float('nan')
    return out
