"""Burst statistics from a Pol II initiation train.

The question this serves: how much transcriptional bursting survives when
H3K27me3 is depleted from the central region, and how it falls as me3
encroaches. Initiation events are the observable; everything here is computed
from their times.
"""
import numpy as np


def burst_stats(times, total_s, gap_factor=0.25, bin_h=1.0):
    """Statistics of one cell's initiation train.

    A burst is a run of initiations whose gaps are shorter than
    `gap_factor` x the mean interval. That is a threshold, not a fit, so it is
    only meaningful compared ACROSS conditions computed the same way.

    Returns rate (per hour), CV^2 of intervals (1 = Poisson, >1 = bursty),
    Fano factor of hourly counts (1 = Poisson), burst frequency per hour, mean
    burst size, and the fraction of initiations that fall in a burst.
    """
    HOUR = 3600.0
    t = np.sort(np.asarray(times, float))
    n = t.size
    out = dict(n=n, rate_h=n / (total_s / HOUR) if total_s > 0 else np.nan,
               cv2=np.nan, fano=np.nan, burst_freq_h=np.nan,
               burst_size=np.nan, frac_in_burst=np.nan)
    if n < 6:
        return out
    iv = np.diff(t)
    mu = iv.mean()
    if mu > 0:
        out["cv2"] = float(iv.var() / mu ** 2)
    nb = max(1, int(total_s / (bin_h * HOUR)))
    counts, _ = np.histogram(t, bins=nb, range=(0, total_s))
    if counts.mean() > 0:
        out["fano"] = float(counts.var() / counts.mean())
    thr = gap_factor * mu
    brk = np.where(iv > thr)[0]
    starts = np.concatenate([[0], brk + 1])
    ends = np.concatenate([brk, [n - 1]])
    sizes = ends - starts + 1
    real = sizes[sizes >= 2]
    out["burst_freq_h"] = float(len(sizes) / (total_s / HOUR))
    out["burst_size"] = float(sizes.mean())
    out["frac_in_burst"] = float(real.sum() / n) if n else np.nan
    return out


def encroachment(profile, x, core_bp=1000.0):
    """How far me3 has come into the central region.

    Returns the mean me3 over +/- core_bp of the TSS -- 0 is a fully cleared
    core, 1 is fully encroached.
    """
    return float(np.asarray(profile)[np.abs(np.asarray(x)) <= core_bp].mean())
