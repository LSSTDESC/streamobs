#!/usr/bin/env python
"""Shared Roman F158 size-envelope star classifier.

ONE source of truth for the single-band F158 size-envelope star classifier used by
the Roman DC2 selection-function products. Both
``scripts/roman/create_streamobs_files_hlwas.py`` (selection-function generator) and
``scripts/roman/build_roman_galaxy_misclass.py`` (galaxy-misclassification curve)
import :func:`build_env_classifier` from here so the classifier cannot drift between
them.

Method (see ``docs/source/selection_function_methodology.md`` for the science):
A detection is a star iff ``lower(mag) < size_sb < upper(mag)``, where
``size_sb = sqrt(lambda1) * 3600"`` is the windowed F158 semi-major axis from the
per-band second moments (``x2/y2/xywin_world_H158``). Working in ``L = log10(size)``,
the band is symmetric about the per-magnitude stellar locus ``L0(mag)`` with a
half-width ``Delta(mag)`` (dex) tuned per magnitude bin to the COMPLETE purity
target (0.875), capped at the bright end by the stellar log-size scatter
(``N_SIG * sigma_L``), single-peaked, PCHIP-splined, and frozen faintward of
``ENV_FREEZE``. The upper bound additionally flares blueward of ``ENV_UP_KNEE`` to
retain bright, slightly-resolved stars.

This module has no I/O and no side effects: it only fits the classifier on a
provided matched det->truth catalog and returns callables + the fitted locus/Delta
artifacts (so the generator can also draw the stellar-locus curve).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

# --------------------------------------------------------------------------- #
# Classifier hyper-parameters (the single authoritative copy).
# --------------------------------------------------------------------------- #
BAND = "H158"  # mock column name == Roman F158
FLAG_CUT = 1  # flags < 1, i.e. flags == 0
ENV_PURITY = 0.875  # COMPLETE target (DES Y6 0<=EXT_XGB<=1)
ENV_N_SIG = 4.0  # bright-end cap: half-width <= N_SIG * stellar log-size scatter
ENV_MAX_D = 0.6  # absolute ceiling on the half-width (dex), safety only
ENV_FREEZE = 24.0  # hold Delta constant (fixed dex offset) faintward of this mag
ENV_UP_KNEE, ENV_UP_BRIGHT, ENV_UP_BRIGHT_VAL = (
    23.0,
    18.0,
    0.15,
)  # upper-bound bright flare (arcsec @ mag)
SN_GATE = 1.0857 / 5.0  # magerr at S/N=5; gate the envelope fit to S/N>5 detections


def size_sb(df) -> np.ndarray:
    """Single-band F158 semi-major axis (arcsec) from windowed second moments.

    ``size_sb = sqrt(lambda1) * 3600"``, where ``lambda1`` is the larger eigenvalue
    of the windowed second-moment matrix (``x2/y2/xywin_world_H158``).
    """
    x2 = np.asarray(df["x2win_world_H158"])
    y2 = np.asarray(df["y2win_world_H158"])
    xy = np.asarray(df["xywin_world_H158"])
    half = 0.5 * (x2 + y2)
    root = np.sqrt(np.clip((0.5 * (x2 - y2)) ** 2 + xy**2, 0, None))
    return np.sqrt(np.clip(half + root, 0, None)) * 3600.0


@dataclass
class EnvelopeClassifier:
    """Fitted F158 size-envelope classifier + the artifacts the docs/plots need.

    Attributes
    ----------
    classify : callable
        ``classify(df) -> bool array`` (True inside the stellar size band).
    env_upper_size, env_lower_size : callable
        Upper/lower size boundary (arcsec) vs F158 mag.
    Dfun : callable
        Half-width ``Delta(mag)`` (dex).
    locus_lin, L0_at, sigL_at : callable
        Stellar locus (linear size, log10 size) and robust log-size scatter vs mag.
    locus_mag, locus_logsize : np.ndarray
        Bin midpoints and the fitted log10 stellar-locus size (for plotting the locus).
    """

    classify: Callable
    env_upper_size: Callable
    env_lower_size: Callable
    Dfun: Callable
    locus_lin: Callable
    L0_at: Callable
    sigL_at: Callable
    locus_mag: np.ndarray
    locus_logsize: np.ndarray


def build_env_classifier(cat: pd.DataFrame) -> EnvelopeClassifier:
    """Fit the F158 size envelope on clean true stars/galaxies in ``cat``.

    ``cat`` must carry (at least) ``matched``, ``flags``, ``truth_gal_star``,
    ``mag_auto_H158``, ``magerr_auto_H158`` and the three ``*win_world_H158`` second
    moments. ``size_sb`` is computed internally if absent. The catalog is not
    mutated.

    Returns an :class:`EnvelopeClassifier`.
    """
    cat = cat.copy()
    if "size_sb" not in cat.columns:
        cat["size_sb"] = size_sb(cat)
    base = (
        cat["matched"]
        & (cat["flags"] < FLAG_CUT)
        & np.isfinite(cat["size_sb"])
        & (cat["size_sb"] > 0)
        & np.isfinite(cat[f"mag_auto_{BAND}"])
        & (cat[f"magerr_auto_{BAND}"] < SN_GATE)
    )

    # stellar locus L0(mag) and robust log-size scatter sigma_L(mag) from clean true stars
    fit = cat.loc[base & (cat["truth_gal_star"] == 1)]
    Ls = np.log10(fit["size_sb"].values)
    ms = fit[f"mag_auto_{BAND}"].values
    lb = np.arange(17, 27.51, 0.25)
    lm = 0.5 * (lb[1:] + lb[:-1])
    mu = np.full(lm.size, np.nan)
    sg = np.full(lm.size, np.nan)
    sbin = np.digitize(ms, lb) - 1
    for i in range(lm.size):
        v = Ls[sbin == i]
        v = v[np.isfinite(v)]
        if v.size >= 30:
            mu[i] = np.median(v)
            sg[i] = 1.4826 * np.median(np.abs(v - mu[i])) + 1e-9
    sm = (
        lambda a: pd.Series(a)
        .rolling(3, center=True, min_periods=1)
        .median()
        .ffill()
        .bfill()
        .to_numpy()
    )
    mu, sg = sm(mu), sm(sg)
    L0_at = lambda m: np.interp(m, lm, mu)
    sigL_at = lambda m: np.interp(m, lm, sg)
    locus_lin = lambda m: 10.0 ** L0_at(m)

    # per-mag purity-target half-width Delta(mag): single-peaked, PCHIP-splined, frozen >24
    eb = np.arange(18, 28.01, 0.25)
    emid = 0.5 * (eb[1:] + eb[:-1])
    ebi = lambda m: np.digitize(m, eb) - 1
    fit_all = cat.loc[base]
    m_all = fit_all[f"mag_auto_{BAND}"].values
    d_all = np.abs(np.log10(fit_all["size_sb"].values) - L0_at(m_all))
    lab_all = fit_all["truth_gal_star"].values == 1

    def fit_delta(X, keep=0.02):
        mb = ebi(m_all)
        Dc = np.full(emid.size, np.nan)
        for i in range(emid.size):
            mm = mb == i
            ns = lab_all[mm].sum()
            if mm.sum() < 200 or ns < 10:
                continue
            ds = d_all[mm]
            ls = lab_all[mm]
            qs = np.unique(np.nanquantile(ds, np.linspace(0, 1, 200)))
            best = qs[0]
            for t in qs[::-1]:  # loose (wide) -> tight
                s = ds < t
                n = s.sum()
                if n == 0:
                    continue
                if (s & ls).sum() / n >= X and (s & ls).sum() / ns >= keep:
                    best = t
                    break
            Dc[i] = best
        valid = np.isfinite(Dc)
        xb = emid[valid]
        cap = np.minimum(
            ENV_N_SIG * sigL_at(xb), ENV_MAX_D
        )  # bright end tracks stellar scatter
        Db = np.minimum(np.clip(Dc[valid], 0, None), cap)
        Db = (
            pd.Series(Db)
            .rolling(5, center=True, min_periods=1)
            .median()
            .to_numpy()
            .copy()
        )
        pk = int(np.argmax(Db))  # single-peaked: rise then only tighten
        Db[: pk + 1] = np.maximum.accumulate(Db[: pk + 1])
        Db[pk:] = np.minimum.accumulate(Db[pk:])
        pch = PchipInterpolator(xb, Db)
        hi = min(ENV_FREEZE, xb[-1])
        return lambda m: pch(np.clip(np.asarray(m, float), xb[0], hi))

    Dfun = fit_delta(ENV_PURITY)

    def env_upper_size(m):
        """Upper size boundary (arcsec): symmetric band, flaring blueward of the knee."""
        m = np.asarray(m, float)
        base_sz = locus_lin(m) * 10.0 ** Dfun(m)
        u_knee = float(locus_lin(ENV_UP_KNEE) * 10.0 ** float(Dfun(ENV_UP_KNEE)))
        frac = np.clip((ENV_UP_KNEE - m) / (ENV_UP_KNEE - ENV_UP_BRIGHT), 0.0, 1.0)
        ramp = u_knee + (ENV_UP_BRIGHT_VAL - u_knee) * frac
        return np.where(m < ENV_UP_KNEE, ramp, base_sz)

    def env_lower_size(m):
        """Lower size boundary (arcsec): tight, symmetric in dex about the locus."""
        m = np.asarray(m, float)
        return locus_lin(m) * 10.0 ** (-Dfun(m))

    def classify_star(df):
        """Envelope star classifier: True where lower(mag) < size_sb < upper(mag)."""
        sz = (
            df["size_sb"].values
            if "size_sb" in getattr(df, "columns", [])
            else size_sb(df)
        )
        m = df[f"mag_auto_{BAND}"].values
        with np.errstate(invalid="ignore"):
            ok = np.isfinite(sz) & (sz > 0) & np.isfinite(m)
            return ok & (sz > env_lower_size(m)) & (sz < env_upper_size(m))

    return EnvelopeClassifier(
        classify=classify_star,
        env_upper_size=env_upper_size,
        env_lower_size=env_lower_size,
        Dfun=Dfun,
        locus_lin=locus_lin,
        L0_at=L0_at,
        sigL_at=sigL_at,
        locus_mag=lm,
        locus_logsize=mu,
    )
