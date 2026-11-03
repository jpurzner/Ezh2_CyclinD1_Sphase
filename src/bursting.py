"""Canonical CDKI transcriptional-bursting driver (JP 2026-11-02).

The burst SPECIFICATION is baked into the model (`with_cdki_bursting=True` declares burst_duty, burst_tau,
burst_a_on/off, burst_a_on/off_p18). This module is the single implementation that APPLIES it, so every sim
and figure inherits one canonical bursting behaviour instead of re-declaring a telegraph.

Why the switching lives here and not in the antimony: a deterministic ODE cannot generate independent
per-cell stochastic switching, and SSA is infeasible for this 65-species event-driven model. So the driver
steps a two-state (telegraph) process per gene and writes the multiplier onto the model's `*_tx` knobs,
which scale ONLY transcription -> mRNA -> protein (no leak into mu_eff).

Established behaviour this encodes (see docs/FIGURES.md + memory):
  - bursts are INDEPENDENT per gene (rho=0): correlating them destroys MB-vs-GNP specificity.
  - p21 (p21a_tx) is the dominant transient-G0 driver; p18 (p18_tx) is a milder secondary.

Typical use:

    import tellurium as te
    from src.build_model_v44_heldt import build_model_v44
    from src.bursting import CdkiBurstDriver

    rr = te.loada(build_model_v44())
    drv = CdkiBurstDriver(rr, seed=11)
    for each cell:
        rr.reset(); ...set per-cell params...
        drv.reset()                      # fresh telegraph state for this cell
        while t < T:
            drv.step(DT)                 # updates p21a_tx / p18_tx on the model
            rr.simulate(t, t + DT, 2, selections=SEL)
            t += DT

Set `p21_only=True` to hold p18 static (the attribution control).
"""
import numpy as np

# model parameter name -> (spec param for ON, spec param for OFF)
GENE_KNOBS = {
    'p21a_tx': ('burst_a_on', 'burst_a_off'),
    'p18_tx': ('burst_a_on_p18', 'burst_a_off_p18'),
}
_SPEC = ('burst_duty', 'burst_tau', 'burst_a_on', 'burst_a_off', 'burst_a_on_p18', 'burst_a_off_p18')


class CdkiBurstDriver:
    """Independent two-state telegraph on each CDKI transcription knob, configured from the model."""

    def __init__(self, rr, seed=None, genes=('p21a_tx', 'p18_tx'), p21_only=False, **overrides):
        self.rr = rr
        self.rng = np.random.default_rng(seed)
        self.genes = ['p21a_tx'] if p21_only else list(genes)
        self._static = [g for g in genes if g not in self.genes]   # held at 1.0 (not bursting)
        cfg = {}
        for p in _SPEC:
            try:
                cfg[p] = float(rr[p])
            except Exception as e:                                  # model built without with_cdki_bursting
                raise RuntimeError(
                    f"model has no '{p}': rebuild with build_model_v44(with_cdki_bursting=True) "
                    f"(it is the default), or pass the value explicitly as an override."
                ) from e
        cfg.update(overrides)
        self.duty = cfg['burst_duty']; self.tau = cfg['burst_tau']
        self.amp = {'p21a_tx': (cfg['burst_a_on'], cfg['burst_a_off']),
                    'p18_tx': (cfg['burst_a_on_p18'], cfg['burst_a_off_p18'])}
        self.k_on = self.duty / self.tau
        self.k_off = (1.0 - self.duty) / self.tau
        self.reset()

    def reset(self):
        """New telegraph state — call at the start of each cell so realizations are independent."""
        self.on = {g: bool(self.rng.random() < self.duty) for g in self.genes}
        self._apply()
        for g in self._static:
            try: self.rr[g] = 1.0
            except Exception: pass

    def _apply(self):
        for g in self.genes:
            a_on, a_off = self.amp[g]
            self.rr[g] = a_on if self.on[g] else a_off

    def step(self, dt):
        """Advance each gene's telegraph by dt (min) and write the multipliers onto the model."""
        for g in self.genes:
            if self.on[g]:
                if self.rng.random() < self.k_off * dt: self.on[g] = False
            elif self.rng.random() < self.k_on * dt:
                self.on[g] = True
        self._apply()

    @property
    def state(self):
        return dict(self.on)


def simulate_bursting(rr, t0, t1, dt, driver, selections=None):
    """Run rr from t0->t1 in dt chunks with the telegraph applied each chunk.

    Returns {selection: np.array} of the per-chunk endpoint values, or None if the integrator failed
    (retries with progressively smaller internal steps / looser tolerance first, as deep CDKI arrest is
    stiff). Does NOT reset the driver — call driver.reset() yourself per cell.
    """
    sel = list(selections) if selections else None
    out = {k: [] for k in (sel or [])}
    t = t0
    while t < t1:
        driver.step(dt)
        r = None
        for pts, rtol, atol in [(2, 1e-6, 1e-8), (12, 1e-5, 1e-7), (40, 1e-4, 1e-6), (120, 1e-3, 1e-5)]:
            try:
                rr.integrator.setValue('relative_tolerance', rtol)
                rr.integrator.setValue('absolute_tolerance', atol)
                r = rr.simulate(t, t + dt, pts, selections=sel)
                break
            except Exception:
                continue
        if r is None:
            return None
        if sel:
            for k in sel:
                out[k].append(float(r[k][-1]))
        t += dt
    return {k: np.array(v) for k, v in out.items()} if sel else {}
