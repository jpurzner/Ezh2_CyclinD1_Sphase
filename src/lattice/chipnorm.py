"""Python port of the chipseq.norm histogram cutoff + long-tail trimmer.

Faithful to R/binarize.R (`binarize_counts`) and R/tail_trim.R
(`trim_upper_tail`) from the chipseq_norm package. Ported because the lattice
dataset is npz/python and the Frontenac node has no R.

Two substitutions, both validated in validate_chipnorm.py:
  * mixtools::normalmixEM  ->  a hand-written 1-D two-component EM with free
    variances (below). Same algorithm and same initialisation as the R
    (means at the 25th/75th percentile, sigma = sd/2, lambda = 0.6/0.4).
    Written out rather than taken from sklearn so that identical code runs
    locally (to validate against the real R) and on the cluster.
  * stats::loess           ->  tricube-weighted LOCAL LINEAR regression
    evaluated on a rank grid. This is loess(degree=1); the grid keeps it O(n)
    instead of O(n^2) so it runs on 10^6+ values.

WHY THIS MATTERS for the lattice dataset: its shipped `me3_bin` is saturated
(core plateau = 1.000 in 11/20 ENCODE mouse tissues; index `cut_me3` as low as
-2.035) and its `state` layer depends on ac/k4 fits that failed in many samples
(genome-wide Active fraction spans 0.001-0.577). Re-deriving the cutoffs with
this method replaces both.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm


# --------------------------------------------------------------- loess(deg=1)
def _local_linear(x, y, xout, span):
    """Tricube-weighted local linear fit — loess with degree=1."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    n = len(x); k = max(3, int(np.ceil(span * n)))
    out = np.empty(len(xout))
    for i, x0 in enumerate(xout):
        d = np.abs(x - x0)
        if k < n:
            h = np.partition(d, k - 1)[k - 1]
        else:
            h = d.max()
        if h <= 0:
            out[i] = y[np.argmin(d)]; continue
        u = np.clip(d / h, 0, 1)
        w = (1 - u ** 3) ** 3
        sw = w.sum()
        if sw <= 0:
            out[i] = y[np.argmin(d)]; continue
        xm = (w * x).sum() / sw; ym = (w * y).sum() / sw
        sxx = (w * (x - xm) ** 2).sum()
        b = (w * (x - xm) * (y - ym)).sum() / sxx if sxx > 0 else 0.0
        out[i] = ym + b * (x0 - xm)
    return out


# --------------------------------------------------------------- tail trimmer
def _local_linear_sorted(x, y, xout, span, degree=2):
    """Tricube local polynomial fit where x is sorted and unit-spaced.

    `degree=2` matches R's `loess()` default, which is what chip_tools calls.
    Degree 1 reads a convex knee as straighter than it is and pushes the
    crossing ~30% further down the rank axis.
    """
    n = len(x); k = max(degree + 2, int(np.ceil(span * n)))
    half = k // 2
    out = np.empty(len(xout))
    pos = np.searchsorted(x, xout)
    for i, (x0, c) in enumerate(zip(xout, pos)):
        a = max(0, c - half); b = min(n, c + half + 1)
        xs = x[a:b]; ys = y[a:b]
        d = np.abs(xs - x0); h = d.max()
        if h <= 0:
            out[i] = ys[np.argmin(d)]; continue
        w = (1 - np.clip(d / h, 0, 1) ** 3) ** 3
        if w.sum() <= 0:
            out[i] = ys[np.argmin(d)]; continue
        dx = xs - x0
        A = np.vander(dx, degree + 1, increasing=True)   # [1, dx, dx^2, ...]
        sw = np.sqrt(w)
        try:
            coef, *_ = np.linalg.lstsq(A * sw[:, None], ys * sw, rcond=None)
            out[i] = coef[0]                              # value at dx = 0
        except np.linalg.LinAlgError:
            out[i] = float(np.average(ys, weights=w))
    return out


N_REF = 20000   # rows in the features x samples matrix the R default was tuned on

# Slope thresholds as K value-units per unit RANK FRACTION; the per-rank
# threshold is K/n. K_MINMAX reproduces the R default (3e-4 per rank at
# n = N_REF on a [0,1] scale).
#
# K_ROBUST is set on the MAD scale against three shapes: a long-tailed
# unimodal, a bimodal with a planted tail, and a bimodal with NO tail. The
# last one is the trap -- its rank curve flattens in the VALLEY between the
# modes, so too permissive a threshold fires there and "trims" a third of the
# data (K = 6 cuts 28.8% off that case). K = 25 clears it with margin, trims
# 2-4% on all three, and on the R reference fixture lands at 0.398 against
# R's 0.405. The answer is insensitive to K between 12 and 50.
K_MINMAX = 6.0
K_ROBUST = 25.0


def background_scale(x, bins=512, smooth=5):
    """Anchor on the BACKGROUND mode and its left half-width.

    Returns (z, centre, scale, ok). `ok` is False when the tallest mode sits
    above the median, which means the marked population outnumbers the
    background and "background mode" is no longer meaningful; the caller then
    falls back to `robust_scale`.

    Why not the median and MAD: both drift as the marked fraction grows,
    because marked features are part of the sample they summarise. On a
    synthetic background N(0,1) plus marked N(4,1.5), the MAD goes 1.02 ->
    2.84 as marked goes 2% -> 50%, so z compresses and the bound fraction
    falls just because there is more signal. Taking the modal bin and
    measuring spread only BELOW it holds at 0.95-1.05 across that whole range,
    because nothing to the right of the mode can move either statistic.
    """
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    h, e = np.histogram(x, bins=bins)
    c = 0.5 * (e[:-1] + e[1:])
    if smooth > 1:
        h = np.convolve(h, np.ones(smooth) / smooth, mode="same")
    centre = float(c[int(np.argmax(h))])
    lo = x[x < centre]
    if lo.size < 50:
        return robust_scale(x) + (False,)
    scale = float(np.median(np.abs(lo - centre)) * 1.4826)
    if scale <= 0:
        q25, q75 = np.percentile(x, [25, 75])
        scale = float((q75 - q25) / 1.349) or float(np.std(x)) or 1.0
    ok = centre <= float(np.median(x))
    return (x - centre) / scale, centre, scale, ok


def robust_scale(x):
    """Center on the median, scale by the MAD. Returns (z, centre, scale).

    This replaces `wins_norm` as the scale that `chip_segment` works on.
    Winsorising clipped the top 0.1% and then min-max scaled, which anchored
    the scale to the 99.9th percentile -- a statistic that MOVES with the
    number of marked features. A sample with more H3K27me3 got a higher
    hi-clip, a wider span, a compressed distribution and therefore a LOWER
    bound fraction, which is why bound fraction came out anti-correlated with
    the marked-gene flag (r = -0.42 over 46 samples).

    The median and MAD are set by the unmarked bulk instead, so adding marked
    features barely moves them, and no data is clipped away.
    """
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    centre = float(np.median(x))
    scale = float(np.median(np.abs(x - centre)) * 1.4826)
    if scale <= 0:                      # discrete / heavily tied data
        q25, q75 = np.percentile(x, [25, 75])
        scale = float((q75 - q25) / 1.349)
    if scale <= 0:
        scale = float(np.std(x)) or 1.0
    return (x - centre) / scale, centre, scale


def trim_upper_tail(values, derivative_threshold="auto", span_value=0.01,
                    span_deriv=0.1, max_rank_frac=1.0, ngrid=1500,
                    rescale=True):
    """Rank-slope upper-tail cutoff. Returns (cutoff_value, crossing_rank).

    Mirrors R/tail_trim.R: rank descending (rank 1 = largest), optionally keep
    only the top `max_rank_frac`, smooth value~rank, take the first derivative
    walking from high rank to low rank, smooth it, then report the first rank
    (from rank 1 upward) at which the smoothed slope falls to or below
    `derivative_threshold`.

    `derivative_threshold="auto"` uses ``3e-4 * N_REF / n``. The R's own help
    warns that 3e-4 "is on the scale of the smoothed derivative, which depends
    on how many rows you have". The slope is per unit RANK, so it falls as 1/n:
    a fixed 3e-4 that works on ~20k features is 1-2 orders of magnitude too
    permissive on a few million per-site values, every slope sits below it, and
    the crossing fires at rank 1 -- reporting that the tail starts at the top of
    the distribution and trimming nothing. Scaling by N_REF/n makes the
    threshold a slope per unit RANK FRACTION, which is sample-size invariant
    because value-vs-rank-fraction is just the quantile function. Pass a float
    to reproduce the R exactly at a given n.

    TWO SCALE CORRECTIONS vs a naive port, both needed for the default
    threshold to mean what it means in the R:
      * the slope is per UNIT RANK (the R diffs consecutive ranks); evaluating
        on a grid requires dividing by the grid spacing.
      * `derivative_threshold=3e-4` is documented as being on a [0,1] value
        scale. With `rescale=True` the values are min-max scaled before
        smoothing and the cutoff is mapped back, so the default transfers to
        log-ratio or rlog data.
    """
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if v.size < 50:
        return np.nan, np.nan
    order = np.argsort(-v, kind="mergesort")
    sv = v[order]                      # descending; index 0 = rank 1
    ranks = np.arange(1, sv.size + 1)
    if max_rank_frac < 1.0:
        cap = int(np.ceil(sv.size * max_rank_frac))
        sv = sv[:cap]; ranks = ranks[:cap]
    if sv.size < 50:
        return np.nan, np.nan
    if derivative_threshold == "auto":
        K = K_ROBUST if rescale == "robust" else K_MINMAX
        derivative_threshold = K / sv.size
    if rescale == "robust":
        work, _c, _s = robust_scale(sv)
        back = lambda u: u * _s + _c
    elif rescale and float(sv.max()) > float(sv.min()):
        lo_v, hi_v = float(sv.min()), float(sv.max())
        work = (sv - lo_v) / (hi_v - lo_v)
        back = lambda u: u * (hi_v - lo_v) + lo_v
    else:
        work = sv
        back = lambda u: u
    grid = np.unique(np.linspace(ranks[0], ranks[-1], min(ngrid, sv.size)).astype(int))
    # ranks are sorted and unit-spaced, so the k nearest neighbours of a grid
    # point are a CONTIGUOUS window. That turns the loess from O(n) per grid
    # point (np.partition over all n) into O(k), which is what makes this
    # usable on the 10^6-10^7 sites in a lattice sample.
    sm = _local_linear_sorted(ranks, work, grid, span_value)
    # derivative taken from HIGH rank to LOW rank, as the R does after order(-Rank)
    g_rev, s_rev = grid[::-1], sm[::-1]
    # |step|: walking down the rank axis the spacing is negative, and the R
    # does not divide at all, so its slopes are positive. Dividing by the
    # signed step flips every sign and the crossing fires at rank 1.
    step = np.abs(np.diff(g_rev))
    step = np.where(step == 0, np.nan, step)
    d1 = np.concatenate([[np.nan], np.diff(s_rev) / step])   # per unit rank, positive
    ok = np.isfinite(d1)
    d1s = np.full_like(d1, np.nan)
    gi = g_rev[ok][::-1]; di = d1[ok][::-1]      # ascending for the sorted fit
    d1s[ok] = _local_linear_sorted(gi, di, gi, span_deriv)[::-1]
    # back to ascending rank, then first crossing from rank 1
    g_asc, d_asc, s_asc = g_rev[::-1], d1s[::-1], s_rev[::-1]
    hit = np.where(np.isfinite(d_asc) & (d_asc <= derivative_threshold))[0]
    if hit.size == 0:
        return np.nan, np.nan
    j = hit[0]
    return float(back(s_asc[j])), float(g_asc[j])


# --------------------------------------------------------------- GMM cutoff
def _em_2comp(x, mu, sigma, lam, maxit=10000, tol=1e-10):
    """1-D two-component Gaussian mixture EM, free variances (arbvar=TRUE)."""
    mu = np.array(mu, float); sigma = np.array(sigma, float); lam = np.array(lam, float)
    prev = -np.inf
    for _ in range(maxit):
        d = np.stack([lam[k] * norm.pdf(x, mu[k], sigma[k]) for k in range(2)])
        tot = d.sum(0)
        tot = np.where(tot <= 0, 1e-300, tot)
        r = d / tot
        ll = float(np.log(tot).sum())
        nk = r.sum(1)
        nk = np.where(nk <= 0, 1e-300, nk)
        lam = nk / x.size
        mu = (r * x).sum(1) / nk
        sigma = np.sqrt(np.maximum((r * (x - mu[:, None]) ** 2).sum(1) / nk, 1e-12))
        if abs(ll - prev) < tol:
            return mu, sigma, lam, ll, True
        prev = ll
    return mu, sigma, lam, prev, False


def binarize_counts(values, mu=None, random_state=0, subsample=200000):
    """2-component Gaussian mixture; cutoff where P(lower component | x) = 0.5.

    Returns dict with cutoff_50, cutoff_75, mu_lower, mu_upper, sigma, lambda_,
    separation, converged. Mirrors R/binarize.R.
    """
    x = np.asarray(values, float)
    x = x[np.isfinite(x)]
    if x.size < 100:
        return dict(converged=False, cutoff_50=np.nan, cutoff_75=np.nan)
    fit_x = x
    if subsample and x.size > subsample:
        rng = np.random.default_rng(random_state)
        fit_x = rng.choice(x, subsample, replace=False)
    if mu is None:
        mu = np.quantile(fit_x, [0.25, 0.75])
    sd = float(np.std(fit_x, ddof=1))
    m, s, w, ll, conv = _em_2comp(fit_x, mu, [sd / 2, sd / 2], [0.6, 0.4])
    lo = int(np.argmin(m))

    def post_lower(xx):
        a = w[lo] * norm.pdf(xx, m[lo], s[lo])
        b = sum(w[i] * norm.pdf(xx, m[i], s[i]) for i in range(2))
        return a / b if b > 0 else 0.0

    def solve(p):
        f = lambda xx: p - post_lower(xx)
        lo_b, hi_b = m[lo], max(float(x.max()), m[1 - lo] + 6 * s[1 - lo])
        try:
            if f(lo_b) * f(hi_b) > 0:
                return float(np.mean(m))
            return float(brentq(f, lo_b, hi_b, xtol=1e-12))
        except Exception:
            return float(np.mean(m))

    return dict(converged=bool(conv), loglik=ll,
                cutoff_50=solve(0.5), cutoff_75=solve(0.75),
                mu_lower=float(m[lo]), mu_upper=float(m[1 - lo]),
                sigma=[float(s[lo]), float(s[1 - lo])],
                lambda_=[float(w[lo]), float(w[1 - lo])],
                separation=float((m[1 - lo] - m[lo]) / s.mean()))


def binarize_with_tail_removal(values, derivative_threshold=3e-4,
                               max_rank_frac=0.01, **kw):
    """Long-tail marks: trim the upper shelf, fit the GMM on the remainder.

    Returns (result_dict, tail_cutoff). Values above `tail_cutoff` are the
    third ('tail') state, as in categorize_marks().
    """
    v = np.asarray(values, float); v = v[np.isfinite(v)]
    tail, _ = trim_upper_tail(v, derivative_threshold=derivative_threshold,
                              max_rank_frac=max_rank_frac)
    core = v[v <= tail] if np.isfinite(tail) else v
    res = binarize_counts(core, **kw)
    res["tail_cutoff"] = float(tail) if np.isfinite(tail) else np.nan
    res["frac_tail"] = float((v > tail).mean()) if np.isfinite(tail) else 0.0
    return res


# ------------------------------------------------------------------ chip_tools
# Ported from github.com/jpurzner/chip_tools (R/soften_zero.R, R/split_counts.R).
# These are NOT in chipseq_norm and are needed for zero-inflated data such as
# log2tpm, where ~26% of genes are exactly zero and the spike dominates the fit.

def soften_zero(counts, scale_max=1.5, breaks=100, random_state=0):
    """Thin a dominant zero bin so it stops swamping a mixture fit.

    Faithful port of chiptools::soften_zero: histogram the data, take the
    tallest bin OTHER than the first, and randomly drop exact zeros until at
    most `scale_max` x that height remain. Dropped values become NaN.
    """
    x = np.asarray(counts, float).copy()
    fin = np.isfinite(x)
    if fin.sum() < 2:
        return x
    hist, _ = np.histogram(x[fin], bins=breaks)
    if hist.size < 2:
        return x
    max_other = hist[1:].max()
    zidx = np.where(x == 0)[0]
    n_remove = len(zidx) - (max_other * scale_max)
    if n_remove <= 0:
        return x
    rng = np.random.default_rng(random_state)
    drop = rng.choice(zidx, int(np.floor(n_remove)), replace=False)
    x[drop] = np.nan
    return x


def split_counts(values, min_k=2, max_k=4, random_state=0, subsample=200000):
    """Fit k = min_k..max_k component mixtures, keep the lowest-BIC model.

    Port of chiptools::split_counts. Matters because forcing k=2 on a
    distribution that is not bimodal returns a cutoff anyway — which is a
    likely source of the implausible shipped cutoffs in the lattice index
    (`cut_me3` as low as -2.035, `cut_k4 = inf`).
    """
    x = np.asarray(values, float); x = x[np.isfinite(x)]
    if x.size < 100:
        return dict(converged=False, best_k=np.nan, cutoff_50=np.nan)
    fit_x = x
    if subsample and x.size > subsample:
        fit_x = np.random.default_rng(random_state).choice(x, subsample, replace=False)
    best = None
    for k in range(min_k, max_k + 1):
        qs = np.quantile(fit_x, np.linspace(0.1, 0.9, k))
        sd = float(np.std(fit_x, ddof=1))
        mu, sig, lam = np.array(qs), np.full(k, sd / 2), np.full(k, 1.0 / k)
        prev = -np.inf; conv = False
        for _ in range(2000):
            d = np.stack([lam[j] * norm.pdf(fit_x, mu[j], sig[j]) for j in range(k)])
            tot = np.where(d.sum(0) <= 0, 1e-300, d.sum(0))
            r = d / tot; ll = float(np.log(tot).sum())
            nk = np.maximum(r.sum(1), 1e-300)
            lam = nk / fit_x.size
            mu = (r * fit_x).sum(1) / nk
            sig = np.sqrt(np.maximum((r * (fit_x - mu[:, None]) ** 2).sum(1) / nk, 1e-12))
            if abs(ll - prev) < 1e-10:
                conv = True; break
            prev = ll
        npar = 3 * k - 1
        bic = -2 * prev + npar * np.log(fit_x.size)
        if best is None or bic < best["bic"]:
            best = dict(bic=bic, k=k, mu=mu.copy(), sigma=sig.copy(),
                        lam=lam.copy(), loglik=prev, converged=conv)
    o = np.argsort(best["mu"])
    m, s, w = best["mu"][o], best["sigma"][o], best["lam"][o]

    def post_lowest(xx):
        a = w[0] * norm.pdf(xx, m[0], s[0])
        b = sum(w[i] * norm.pdf(xx, m[i], s[i]) for i in range(len(m)))
        return a / b if b > 0 else 0.0

    f = lambda xx: 0.5 - post_lowest(xx)
    try:
        cut = float(brentq(f, m[0], max(float(x.max()), m[-1] + 6 * s[-1]), xtol=1e-12))
    except Exception:
        cut = float(0.5 * (m[0] + m[1])) if len(m) > 1 else np.nan
    return dict(converged=bool(best["converged"]), best_k=int(best["k"]),
                bic=float(best["bic"]), cutoff_50=cut,
                mu=m.tolist(), sigma=s.tolist(), lambda_=w.tolist())


def wins_norm(x, low=500, high=20, normalize=True):
    """DEPRECATED -- kept only for R parity tests; `chip_segment` no longer
    uses it. See `robust_scale` for why clipping-and-min-max-scaling made
    bound fractions incomparable across samples.

    Winsorise by POINT COUNT then unit-scale — chiptools::wins_norm_histones.

    Clips to the `low`-th smallest and `high`-th largest value, then maps to
    [0,1]. This is the step that makes every downstream default meaningful:
    `trim_upper_tail`'s 3e-4 threshold, `binarize_counts`'s mu = c(0.3, 0.65)
    and its `exclude_extreme` window are all defined on this scale.

    `low`/`high` are point counts tuned for ~20k genes. Pass fractions of n
    via low_frac/high_frac for a lattice with millions of sites.
    """
    v = np.asarray(x, float); fin = np.isfinite(v)
    if fin.sum() < 10:
        return v
    s = np.sort(v[fin])
    lo = s[min(low, s.size - 1)] if low > 0 else s[0]
    hi = s[-min(high, s.size)] if high > 0 else s[-1]
    w = np.clip(v, lo, hi)
    if not normalize:
        return w
    rng = np.nanmax(w) - np.nanmin(w)
    return np.zeros_like(w) if rng == 0 else (w - np.nanmin(w)) / rng


def background_null_cutoff(values, alpha=0.01):
    """Cutoff from a background NULL rather than a two-component fit.

    Takes the modal bin as the background centre, measures spread from the
    left side only (nothing right of the mode can contaminate it), mirrors
    that into a Gaussian null and cuts at its upper `alpha` quantile.

    Use this instead of `binarize_counts` whenever the marked fraction may be
    small. A two-component mixture needs both components to be visible; when
    marked features are 2-10% of the data the fit splits the BACKGROUND in
    two and puts the cutoff inside it. On a synthetic background N(0,1) plus
    marked N(4,1.5) the GMM calls 38-41% of features positive at a true 2-10%,
    while this returns 3.5/6.0/8.9% against a true 2/5/10% and holds its
    cutoff at 2.1-2.5 across marked fractions from 2% to 50%.

    Assumes the background is roughly symmetric about its mode. On a strongly
    skewed background it is approximate -- check against the histogram.
    """
    v = np.asarray(values, float); v = v[np.isfinite(v)]
    if v.size < 200:
        return dict(cutoff=np.nan, centre=np.nan, scale=np.nan, ok=False)
    _z, centre, sc, ok = background_scale(v)
    cut = float(centre + norm.ppf(1.0 - alpha) * sc)
    return dict(cutoff=cut, centre=float(centre), scale=float(sc), ok=bool(ok),
                frac_pos=float(np.mean(v > cut)))


def chip_segment(values, derivative_threshold="auto", max_rank_frac=1.0,
                 scale="background", method="background_null", alpha=0.01,
                 min_tail_gap=1.0, diagnostics=True, nbins=400):
    """Segment a signal distribution into unbound / bound / tail.

    robust-scale -> trim_upper_tail -> two-component GMM on the remainder.
    Cutoffs come back on the ORIGINAL scale.

    NO WINSORISATION. The previous version clipped the bottom 2.5% and top
    0.1% of points and min-max scaled what was left, which anchored the scale
    to the 99.9th percentile. That statistic rises with the number of marked
    features, so a sample with more signal got a compressed distribution and a
    lower bound fraction -- bound fraction came out ANTI-correlated with the
    marked-gene flag (r = -0.42 across 46 samples). `robust_scale` anchors on
    the median and MAD, which the unmarked bulk sets, and discards nothing.

    `scale="none"` fits on the raw values (for data already on a sensible
    scale); the threshold is then in raw units per unit rank fraction.
    """
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if v.size < 200:
        return dict(converged=False, bound_cutoff=np.nan, tail_cutoff=np.nan)

    if scale == "background":
        z, centre, sc, ok = background_scale(v)
        if not ok:                      # marked mode outnumbers background
            z, centre, sc = robust_scale(v)
    elif scale == "robust":
        z, centre, sc = robust_scale(v)
    elif scale == "none":
        z, centre, sc = v, 0.0, 1.0
    else:
        raise ValueError("scale must be 'background', 'robust' or 'none'")
    back = lambda u: u * sc + centre

    # z is already on the MAD scale, so pass the robust K explicitly: calling
    # with rescale=False would otherwise make "auto" pick the min-max K, which
    # is 4x too permissive here and fires in the valley of a bimodal fit.
    thr = (K_ROBUST / z.size) if derivative_threshold == "auto" \
        else derivative_threshold
    tail_z, rank = trim_upper_tail(z, derivative_threshold=thr,
                                   max_rank_frac=max_rank_frac, rescale=False)
    core = z[z < tail_z] if np.isfinite(tail_z) else z

    # Initialise the two components from the core's own quantiles. The old
    # mu = c(0.5, max - 0.2) and the 0.1 < x < 1 window were both coordinates
    # on the winsorised [0,1] scale and have no meaning here.
    if method == "background_null":
        nul = background_null_cutoff(core, alpha=alpha)
        bound_z = nul["cutoff"]
        fit = dict(converged=nul["ok"], cutoff_50=bound_z,
                   mu_lower=nul["centre"], mu_upper=np.nan,
                   sigma=[nul["scale"], np.nan], lambda_=[np.nan, np.nan],
                   separation=np.nan)
    elif method == "gmm":
        fit = binarize_counts(core) if core.size >= 200 else dict(
            converged=False, cutoff_50=np.nan)
        bound_z = fit.get("cutoff_50", np.nan)
    else:
        raise ValueError("method must be 'background_null' or 'gmm'")
    tail_cut = float(back(tail_z)) if np.isfinite(tail_z) else np.nan
    bound_cut = float(back(bound_z)) if np.isfinite(bound_z) else np.nan
    # When marked features are a small minority, the knee the trimmer finds IS
    # the marked population, not an extreme shelf beyond it, and the two
    # cutoffs land on top of each other. Require at least `min_tail_gap` MADs
    # between them: on a synthetic sweep the gap is 0.1-0.35 MAD when there is
    # no real shelf and 3.5-4.5 MAD when there is.
    tail_collapsed = (np.isfinite(tail_z) and np.isfinite(bound_z)
                      and (tail_z - bound_z) < min_tail_gap)
    if tail_collapsed:
        tail_cut = np.nan
    out = dict(
        bound_cutoff=bound_cut, tail_cutoff=tail_cut,
        bound_cutoff_z=float(bound_z), tail_cutoff_z=float(tail_z),
        crossing_rank=float(rank), centre=float(centre), scale=float(sc),
        separation=fit.get("separation", np.nan),
        converged=bool(fit.get("converged", False)),
        mu=[fit.get("mu_lower", np.nan), fit.get("mu_upper", np.nan)],
        sigma=fit.get("sigma", [np.nan, np.nan]),
        lambda_=fit.get("lambda_", [np.nan, np.nan]),
        tail_collapsed=bool(tail_collapsed),
        frac_tail=float(np.mean(v > tail_cut)) if np.isfinite(tail_cut) else 0.0,
        frac_pos=float(np.mean(v > bound_cut)) if np.isfinite(bound_cut) else np.nan,
        frac_bound=(float(np.mean((v > bound_cut) & (v <= tail_cut)))
                    if np.isfinite(tail_cut) and np.isfinite(bound_cut)
                    else float(np.mean(v > bound_cut))
                    if np.isfinite(bound_cut) else np.nan))
    if diagnostics:
        out["hist"] = segmentation_hist(v, out, nbins=nbins)
    return out


def segmentation_hist(values, seg, nbins=400, lo_q=0.0005, hi_q=0.9995):
    """Everything needed to redraw a `chip_segment` call, without the values.

    Kept small (a few kB) so it can be carried in a results file and plotted
    later -- the point being that a cutoff should never be reported as a bare
    number. Range is clipped to [lo_q, hi_q] for display only; `n_above` says
    how many points fall past the right edge.
    """
    v = np.asarray(values, float); v = v[np.isfinite(v)]
    lo = float(np.quantile(v, lo_q)); hi = float(np.quantile(v, hi_q))
    if not (hi > lo):
        lo, hi = float(v.min()), float(v.max()) or 1.0
    counts, edges = np.histogram(v, bins=nbins, range=(lo, hi))
    return dict(counts=counts.astype(int).tolist(),
                lo=lo, hi=hi, n=int(v.size),
                n_above=int((v > hi).sum()), n_below=int((v < lo).sum()),
                bound_cutoff=seg.get("bound_cutoff"),
                tail_cutoff=seg.get("tail_cutoff"),
                tail_collapsed=bool(seg.get("tail_collapsed", False)),
                centre=seg.get("centre"), scale=seg.get("scale"),
                frac_pos=seg.get("frac_pos"), frac_tail=seg.get("frac_tail"))
