#!/usr/bin/env python
"""Measure the LSST<->Roman star/galaxy CLASSIFICATION correlation in DC2 (task T3, Part A).

WHY
---
StreamObs draws an independent uniform per object *per survey* when deciding
"selected in survey k?" (``StreamInjector.inject`` spawns a child RNG per survey,
``detect_flag`` draws its own uniforms), and the joint analysis channel ANDs the two
per-survey flags (``fig4_metricC2_ingredients._cmd_cols``, ``which == "joint"``).  The
joint selection probability is therefore ``p_lsst * p_roman`` by construction.

Physically, whether a compact galaxy is mistaken for a star is driven by its angular
size relative to each PSF, so a galaxy compact enough to fool Roman (0.11" PSF FWHM
in H158-ish terms) is almost certainly compact enough to fool LSST (~0.7" seeing).
That coupling is close to COMONOTONIC: ``P(both) ~ min(p_lsst, p_roman)``.  Roman-Rubin
DC2 shares one input truth catalog, so the joint rate is directly MEASURABLE -- that is
what this script does.

WHAT IS AND IS NOT MEASURED HERE
--------------------------------
The matched table is built from Roman DETECTIONS positionally matched to LSST OBJECTS
(1"), so every row is an object DETECTED BY BOTH surveys.  This script therefore
measures the coupling of the CLASSIFICATION step conditional on double detection:

    p_lsst   = P(LSST calls it point-like  | detected by both)
    p_roman  = P(Roman calls it point-like | detected by both)
    p_joint  = P(both call it point-like   | detected by both)

That is the dominant term for the contaminant question, because the shipped
misclassification curves are themselves conditional on detection
(``lsst_dc2_galaxy_misclass_cutr.csv``: "COMPACT true galaxies ... that are DETECTED
(S/N>5, clean) yet classified as point sources"; the Roman curve's denominator is
likewise the detected+matched galaxies), and the injector multiplies that curve by the
detection efficiency (``Survey.get_gal_misclassification_detection`` ->
``mis(dm) * det(dm)``).  The remaining DETECTION-side coupling needs a truth-centric
join (objects LSST detects but Roman misses are absent from this table by construction)
and is NOT covered here; see the "detection coupling" note at the end of the printout.

THE HETEROGENEITY TRAP (why the naive ratio is not the right null)
-----------------------------------------------------------------
Within any single Delta-m bin the objects are not identical: a galaxy that is bright in
F158 is usually also bright in r, so the per-object probabilities p_l(i), p_r(i) are
themselves correlated across the bin.  Even under the code's INDEPENDENT draws,

    E[1(both)] = mean_i p_l(i) p_r(i)  =  mean(p_l) mean(p_r) + Cov_i(p_l, p_r)  >  p_lsst p_roman

so ``p_joint / (p_lsst * p_roman) > 1`` is expected even when the code is right.  The
primary null here is therefore the CELL-WEIGHTED independent prediction: the sample is
partitioned into fine 2-D cells of (Delta-m_lsst, Delta-m_roman), independence is
imposed within each cell, and the predictions are re-aggregated,

    pred_indep  = sum_c n_c * p_l,c * p_r,c / n
    pred_comono = sum_c n_c * min(p_l,c, p_r,c) / n

which removes the bulk of the induced covariance.  Both the naive and cell-weighted
versions are reported; ``r_indep``/``r_comono`` use the cell-weighted ones.

INTERPRETATION RULE (from the handoff)
--------------------------------------
* ``r_indep  ~ 1`` -> the current independent draws are correct; Parts B/C become a
  robustness footnote.
* ``r_comono ~ 1`` -> adopt the shared-uniform coupling and rerun the joint grid.
* in between   -> fit the Gaussian-copula rho per source type that reproduces
  ``p_joint(Delta m)`` (``rho`` is solved for here, per source type and per bin).

INPUTS
------
* ``data/surveys/roman_dc2/roman_lsst_matched.parquet`` -- one row per Roman detection
  positionally matched to an LSST object (built by ``build_roman_galaxy_misclass.py``).
* ``data/surveys/roman_dc2/roman_dc2_det_truth.parquet`` -- to fit the Roman F158
  size-envelope classifier exactly as the shipped products do (cached; see ``--refit``).
* ``data/surveys/roman_dc2/cosmodc2_galaxy_size_true.parquet`` -- cosmoDC2 true sizes
  for the compact-galaxy cut, joined on ``lsst_true_cosmodc2_id``.

CLASSIFIERS (identical to the shipped selection functions)
* Roman : F158 size-envelope classifier, ``roman_star_classifier.build_env_classifier``
  (``ENV_PURITY = 0.875``, ``size_sb = sqrt(lambda1)*3600"`` from the windowed H158
  second moments).
* LSST  : ``extendedness < 0.5`` (``build_lsst_dc2_galaxy_misclass.EXT_CUT``).

Run (streamobs env):
    /astro/store/shiren/conda-envs/stream_team/envs/streamobs/bin/python \
        scripts/roman/measure_joint_misclassification.py [--size-source cosmodc2_true] [--refit]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import healpy as hp
import numpy as np
import pandas as pd
import pyarrow.dataset as pads
import pyarrow.parquet as pq
from scipy.optimize import brentq
from scipy.stats import norm

REPO = Path(__file__).resolve().parents[2]
ROMAN_DIR = REPO / "data" / "surveys" / "roman_dc2"
LSST_DIR = REPO / "data" / "surveys" / "lsst_dc2"
ART_DIR = REPO / "artifacts"

MATCHED = ROMAN_DIR / "roman_lsst_matched.parquet"
DET_TRUTH = ROMAN_DIR / "roman_dc2_det_truth.parquet"
COSMODC2_SIZES = ROMAN_DIR / "cosmodc2_galaxy_size_true.parquet"
ROMAN_MAGLIM = ROMAN_DIR / "roman_dc2_maglim_f158_nside1024.fits.gz"
LSST_MAGLIM = LSST_DIR / "lsst_dc2_maglim_r_nside1024.fits.gz"
CLF_CACHE = ROMAN_DIR / "_cache_env_classifier_bounds.npz"

OUT_CSV = ART_DIR / "joint_misclassification_correlation.csv"
OUT_PNG = ART_DIR / "joint_misclassification_correlation.png"
OUT_META = ART_DIR / "joint_misclassification_correlation.meta.json"

# --- conventions copied from the two shipped product builders ---------------- #
FLAG_CUT = 1  # Roman SExtractor flags == 0            (build_roman_galaxy_misclass)
EXT_CUT = 0.5  # LSST extendedness < 0.5 -> point source (build_lsst_dc2_galaxy_misclass)
TRUTH_GAL, TRUTH_STAR = 1, 2  # LSST truth_type
GAL_SIZE_MAX = 0.3  # arcsec: "compact" galaxy (both products)
JOINT_F158_CUT = 25.5  # the joint channel's Roman cut (fig4_metricC2_ingredients._CUTS)

# Delta-m binning for the reported curves (relative to each survey's reference maglim)
DM_BINS = np.arange(-9.0, 2.001, 0.5)
DM_MID = 0.5 * (DM_BINS[1:] + DM_BINS[:-1])
N_MIN = 200  # minimum objects per reported bin
CELL = 0.25  # mag: 2-D cell size for the heterogeneity-corrected null

MATCHED_COLS = [
    "roman_true_gal_star",
    "roman_true_mag_H158",
    "roman_true_ind",
    "roman_obs_flags",
    "roman_obs_mag_auto_H158",
    "roman_obs_magerr_auto_H158",
    "roman_obs_x2win_world_H158",
    "roman_obs_y2win_world_H158",
    "roman_obs_xywin_world_H158",
    "roman_obs_matched",
    "roman_obs_match_sep_arcsec",
    "lsst_obs_extendedness",
    "lsst_obs_clean",
    "lsst_true_truth_type",
    "lsst_true_mag_r",
    "lsst_true_matched",
    "lsst_true_cosmodc2_id",
    "lsst_obs_ra",
    "lsst_obs_dec",
]

sys.path.insert(0, str(Path(__file__).resolve().parent))
from roman_star_classifier import build_env_classifier, size_sb  # noqa: E402


# --------------------------------------------------------------------------- #
# statistics helpers
# --------------------------------------------------------------------------- #
def wilson(k, n, z=1.0):
    """Wilson score interval (default z=1 -> 68% CI). Returns (lo, hi) arrays."""
    k = np.asarray(k, float)
    n = np.asarray(n, float)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = k / n
        d = 1.0 + z**2 / n
        c = (p + z**2 / (2 * n)) / d
        h = z * np.sqrt(np.clip(p * (1 - p) / n + z**2 / (4 * n**2), 0, None)) / d
    return np.where(n > 0, c - h, np.nan), np.where(n > 0, c + h, np.nan)


_GL_X, _GL_W = np.polynomial.legendre.leggauss(64)


def bvn_cdf(z1, z2, rho):
    """Bivariate normal CDF Phi2(z1, z2; rho), vectorized over z1/z2 for a scalar rho.

    Plackett's identity: Phi2(z1,z2;rho) = Phi(z1)Phi(z2) + int_0^rho phi2(z1,z2;t) dt,
    evaluated with 64-node Gauss-Legendre quadrature.  Accurate to ~1e-10 for |rho|<0.99
    and far faster than calling a multivariate-normal CDF per element.
    """
    z1 = np.asarray(z1, float)
    z2 = np.asarray(z2, float)
    base = norm.cdf(z1) * norm.cdf(z2)
    if rho == 0.0:
        return base
    t = 0.5 * rho * (_GL_X + 1.0)  # nodes on [0, rho]
    w = 0.5 * rho * _GL_W
    om = 1.0 - t**2
    pdf = np.exp(
        -(z1[..., None] ** 2 - 2 * t * z1[..., None] * z2[..., None] + z2[..., None] ** 2)
        / (2 * om)
    ) / (2 * np.pi * np.sqrt(om))
    return base + (pdf * w).sum(axis=-1)


def gauss_copula_joint(p1, p2, rho):
    """P(U1<=p1, U2<=p2) under a Gaussian copula with correlation rho (vectorized)."""
    p1 = np.asarray(p1, float)
    p2 = np.asarray(p2, float)
    deg = (p1 <= 0) | (p2 <= 0) | (p1 >= 1) | (p2 >= 1)
    safe1 = np.where(deg, 0.5, p1)
    safe2 = np.where(deg, 0.5, p2)
    if rho >= 0.999:
        out = np.minimum(p1, p2)
    elif rho <= -0.999:
        out = np.clip(p1 + p2 - 1.0, 0, None)
    else:
        out = bvn_cdf(norm.ppf(safe1), norm.ppf(safe2), rho)
    # degenerate marginals: the copula value is fixed by the marginals alone
    return np.where(deg, np.minimum(p1, p2) * ((p1 > 0) & (p2 > 0)), out)


def fit_rho(w, p_l, p_r, p_joint_target):
    """Gaussian-copula rho reproducing the cell-weighted joint rate. NaN if unbracketed.

    ``w``/``p_l``/``p_r`` are per-cell weights and marginals; the model prediction is the
    weight-averaged copula joint probability.
    """

    def model(rho):
        return float(np.sum(w * gauss_copula_joint(p_l, p_r, rho)))

    if not np.isfinite(p_joint_target):
        return np.nan
    lo, hi = model(-0.99), model(0.99)
    if not (lo <= p_joint_target <= hi):
        return np.nan
    try:
        return float(brentq(lambda r: model(r) - p_joint_target, -0.99, 0.99, xtol=1e-4))
    except Exception:
        return np.nan


def cell_predictions(dm_l, dm_r, star_l, star_r, cell=CELL, fit_copula=True):
    """Heterogeneity-corrected independent / comonotonic nulls (+ copula rho).

    Partitions on (Delta-m_lsst, Delta-m_roman) cells of size ``cell`` mag, imposes the
    coupling inside each cell, and re-aggregates by cell occupancy.  Cells with a single
    object carry no information about within-cell independence but are still weighted
    (their p's are 0/1); that is the usual bias-variance trade-off of this estimator, so
    the naive (single-cell) numbers are reported alongside.

    Both Delta-m arrays must be finite. A non-finite entry cannot be assigned to a
    cell, and silently dropping it would compute the nulls on a different sample than
    the measured ``p_joint`` it is compared against — so it is rejected loudly. Filter
    the sample (e.g. require a finite truth magnitude in BOTH bands) before calling.
    """
    n = len(dm_l)
    if n == 0:
        return dict(pred_indep=np.nan, pred_comono=np.nan, rho=np.nan, n_cells=0)
    bad = ~(np.isfinite(dm_l) & np.isfinite(dm_r))
    if bad.any():
        raise ValueError(
            f"cell_predictions got {int(bad.sum())}/{n} objects with a non-finite "
            "Delta-m; require a finite truth magnitude in both bands when building "
            "the sample (dropping them here would bias the nulls relative to the "
            "measured joint rate)."
        )
    il = np.floor(dm_l / cell).astype(np.int64)
    ir = np.floor(dm_r / cell).astype(np.int64)
    # dense 2-D bin index (no hashing, so no collisions)
    il -= il.min()
    ir -= ir.min()
    key = il * (ir.max() + 1) + ir
    n_c = np.bincount(key).astype(float)
    k_l = np.bincount(key, weights=star_l.astype(float))
    k_r = np.bincount(key, weights=star_r.astype(float))
    occ = n_c > 0
    n_c, k_l, k_r = n_c[occ], k_l[occ], k_r[occ]
    w = n_c / n
    p_l, p_r = k_l / n_c, k_r / n_c
    out = dict(
        pred_indep=float(np.sum(w * p_l * p_r)),
        pred_comono=float(np.sum(w * np.minimum(p_l, p_r))),
        n_cells=int(occ.sum()),
        rho=np.nan,
    )
    if fit_copula:
        out["rho"] = fit_rho(w, p_l, p_r, float((star_l & star_r).mean()))
    return out


# --------------------------------------------------------------------------- #
# inputs
# --------------------------------------------------------------------------- #
def local_maglim(path, ra, dec):
    """Per-object maglim from a band's HEALPix depth map (the injector's own lookup).

    ``_inject_one_survey`` evaluates the selection curves at
    ``survey.get_maglim(band, pixel=hp.ang2pix(nside_of_map, ra, dec))``, i.e. at each
    object's LOCAL depth -- not at the map median. Mirrors that, returning NaN on
    uncovered pixels.
    """
    m = hp.read_map(str(path))
    nside = hp.get_nside(m)
    pix = hp.ang2pix(nside, np.asarray(ra, float), np.asarray(dec, float), lonlat=True)
    out = np.asarray(m, float)[pix]
    return np.where(np.isfinite(out) & (out > 0) & (out != hp.UNSEEN), out, np.nan)


def ref_maglim(path):
    m = hp.read_map(str(path))
    good = (m != hp.UNSEEN) & np.isfinite(m) & (m > 0)
    return float(np.median(m[good])), int(good.sum())


CLF_MAG_GRID = np.arange(14.0, 30.0 + 1e-9, 0.005)


def get_classifier(refit=False):
    """Return ``classify(mag, size) -> bool`` for the F158 size-envelope classifier.

    The envelope is fit on the FULL Roman det->truth catalog -- exactly the Stage-B fit
    of ``build_roman_galaxy_misclass.py``, so this cannot drift from the shipped
    selection functions.  The fitted classifier holds closures (not picklable), so the
    cache stores the two size boundaries sampled on a 0.005-mag grid; the boundaries are
    smooth (PCHIP splines of a monotone half-width), so linear interpolation on that grid
    is exact to well below the size resolution.
    """
    if CLF_CACHE.exists() and not refit:
        print(f"loading cached envelope bounds {CLF_CACHE.name} (--refit to rebuild)")
        z = np.load(CLF_CACHE)
        mg, lo, hi = z["mag"], z["lower"], z["upper"]
    else:
        print("fitting the F158 size-envelope classifier on the full Roman catalog ...")
        t0 = time.time()
        dset = pads.dataset(str(DET_TRUTH), format="parquet")
        fit_cols = [
            "matched",
            "flags",
            "mag_auto_H158",
            "magerr_auto_H158",
            "truth_gal_star",
            "x2win_world_H158",
            "y2win_world_H158",
            "xywin_world_H158",
        ]
        fit_cat = dset.to_table(
            columns=fit_cols, filter=pads.field("matched") == True  # noqa: E712
        ).to_pandas()
        clf = build_env_classifier(fit_cat)
        del fit_cat
        print(
            f"  Delta@21={float(clf.Dfun(21)):.3f} Delta@24={float(clf.Dfun(24)):.3f} dex "
            f"({time.time() - t0:.0f}s, {len(CLF_MAG_GRID)}-point boundary grid cached)"
        )
        mg = CLF_MAG_GRID
        lo = np.asarray(clf.env_lower_size(mg), float)
        hi = np.asarray(clf.env_upper_size(mg), float)
        np.savez_compressed(CLF_CACHE, mag=mg, lower=lo, upper=hi)

    def classify(mag, size):
        mag = np.asarray(mag, float)
        size = np.asarray(size, float)
        with np.errstate(invalid="ignore"):
            ok = np.isfinite(size) & (size > 0) & np.isfinite(mag)
            return (
                ok
                & (size > np.interp(mag, mg, lo))
                & (size < np.interp(mag, mg, hi))
            )

    return classify


def load_matched():
    print(f"reading {MATCHED.name} ({len(MATCHED_COLS)} columns) ...")
    df = pq.read_table(MATCHED, columns=MATCHED_COLS).to_pandas()
    print(f"  {len(df):,} rows")
    return df


def join_true_sizes(df):
    sz = pd.read_parquet(COSMODC2_SIZES, columns=["cosmodc2_id", "size_true"])
    lut = pd.Series(sz["size_true"].values, index=sz["cosmodc2_id"].astype("i8").values)
    ids = pd.to_numeric(df["lsst_true_cosmodc2_id"], errors="coerce")
    out = ids.map(lut).to_numpy(dtype=float)
    print(f"  cosmoDC2 size_true coverage: {np.isfinite(out).mean():.2%}")
    return out


# --------------------------------------------------------------------------- #
# the measurement
# --------------------------------------------------------------------------- #
def measure(sub, label, size_desc, axis="roman", conditioning="scalar"):
    """Per-Delta-m-bin rates + coupling nulls for one source type.

    ``axis`` selects which survey's Delta-m the curve is binned on ("roman" -> F158,
    "lsst" -> r).  The cell-weighted nulls always partition on BOTH axes, so they are
    unaffected by this choice; only the reported binning changes.

    ``conditioning`` selects WHICH Delta-m feeds the cells:

    * ``"scalar"`` -- Delta-m from each band's median (reference) maglim, the
      convention the shipped curves are keyed on.
    * ``"local"`` -- Delta-m from each object's own per-pixel maglim, which is what
      the injector actually evaluates its curves at
      (``_inject_one_survey`` -> ``survey.get_maglim(band, pixel=...)``).

    This matters for the verdict. Under "scalar", any spatial depth variation that is
    correlated BETWEEN the two surveys sits inside the cells as unmodelled shared
    heterogeneity and inflates ``r_indep`` -- so the scalar ``r_indep`` is an UPPER
    BOUND on the code's error, not a point estimate. Conditioning on local maglim
    removes that term and gives the like-for-like comparison against what the code
    does. The reported binning stays on the scalar Delta-m either way, so rows from
    the two conditionings are directly comparable bin by bin.
    """
    suffix = "" if conditioning == "scalar" else "_local"
    dm_l = sub["dm_lsst" + suffix].to_numpy()
    dm_r = sub["dm_roman" + suffix].to_numpy()
    dm_axis = (
        sub["dm_roman"].to_numpy() if axis == "roman" else sub["dm_lsst"].to_numpy()
    )
    sl = sub["star_lsst"].to_numpy()
    sr = sub["star_roman"].to_numpy()
    sj = sl & sr

    rows = []

    def _row(name, mask, dm_bin=np.nan):
        n = int(mask.sum())
        if n < N_MIN:
            return
        kl, kr, kj = int(sl[mask].sum()), int(sr[mask].sum()), int(sj[mask].sum())
        pl, pr, pj = kl / n, kr / n, kj / n
        cp = cell_predictions(dm_l[mask], dm_r[mask], sl[mask], sr[mask])
        lo_l, hi_l = wilson(kl, n)
        lo_r, hi_r = wilson(kr, n)
        lo_j, hi_j = wilson(kj, n)
        rows.append(
            dict(
                source_type=label,
                selection=name,
                dm_axis=axis,
                conditioning=conditioning,
                dm_bin=dm_bin,
                n=n,
                p_lsst=pl,
                p_lsst_lo=float(lo_l),
                p_lsst_hi=float(hi_l),
                p_roman=pr,
                p_roman_lo=float(lo_r),
                p_roman_hi=float(hi_r),
                p_joint=pj,
                p_joint_lo=float(lo_j),
                p_joint_hi=float(hi_j),
                pred_indep=cp["pred_indep"],
                pred_comono=cp["pred_comono"],
                pred_indep_naive=pl * pr,
                pred_comono_naive=min(pl, pr),
                r_indep=pj / cp["pred_indep"] if cp["pred_indep"] > 0 else np.nan,
                r_comono=pj / cp["pred_comono"] if cp["pred_comono"] > 0 else np.nan,
                r_indep_naive=pj / (pl * pr) if pl * pr > 0 else np.nan,
                r_comono_naive=pj / min(pl, pr) if min(pl, pr) > 0 else np.nan,
                rho_copula=cp["rho"],
                n_cells=cp["n_cells"],
                size_source=size_desc,
            )
        )

    for lo, hi, mid in zip(DM_BINS[:-1], DM_BINS[1:], DM_MID):
        _row("dm_bin", (dm_axis >= lo) & (dm_axis < hi), dm_bin=float(mid))
    if axis == "roman":  # axis-independent aggregates: emit once
        _row("ALL", np.ones(len(sub), bool))
        _row("F158_lt_25.5", sub["mag_roman_true"].to_numpy() < JOINT_F158_CUT)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# validation: do our marginals reproduce the SHIPPED selection-function curves?
# --------------------------------------------------------------------------- #
PROD_MAG_BINS = np.arange(15.0, 29.01, 0.25)  # both product builders use this grid
PROD_MAG_MID = 0.5 * (PROD_MAG_BINS[1:] + PROD_MAG_BINS[:-1])


def _binned_rate(mag, flag, n_min=20):
    n = np.histogram(mag, PROD_MAG_BINS)[0]
    k = np.histogram(mag[flag], PROD_MAG_BINS)[0]
    with np.errstate(invalid="ignore"):
        p = np.where(n >= n_min, k / np.maximum(n, 1), np.nan)
    return n, p


def validate_marginals(mag_roman, star_roman, compact_meas, mag_lsst, star_lsst, compact_true):
    """Reproduce each shipped misclassification curve from THIS sample and compare.

    The marginals must land on the shipped curves, otherwise the measured joint rate is
    being compared against the wrong nulls.  Exact agreement is not expected:

    * Roman: the shipped curve uses the same merged table and the measured-size compact
      cut, but without this script's extra LSST-side quality cuts (clean, finite
      extendedness, LSST truth matched), so small differences are normal.
    * LSST: the shipped curve is built over the FULL DC2 object/truth skims with an
      explicit psFlux S/N>5 cut, whereas this sample is restricted to the Roman
      footprint overlap and has no S/N column -- so the amplitude can differ by a few
      per cent; the SHAPE is the thing to check.
    """
    print("\n----- marginal validation vs the shipped curves -----")
    for name, csv, mag, flag, comp, magcol in [
        (
            "Roman F158 (measured-size compact cut)",
            ROMAN_DIR / "roman_galaxy_misclass_cutf158.csv",
            mag_roman,
            star_roman,
            compact_meas,
            "mag_F158",
        ),
        (
            "LSST r (cosmoDC2 true-size compact cut)",
            LSST_DIR / "lsst_dc2_galaxy_misclass_cutr.csv",
            mag_lsst,
            star_lsst,
            compact_true,
            "mag_r",
        ),
    ]:
        ship = pd.read_csv(csv, comment="#", header=None, names=[magcol, "delta_mag", "eff"])
        n, p = _binned_rate(mag[comp], flag[comp])
        mine = pd.DataFrame({magcol: PROD_MAG_MID, "n": n, "p_here": p}).dropna()
        cmp_ = ship.merge(mine, on=magcol, how="inner")
        cmp_ = cmp_[(cmp_.eff > 1e-3) & (cmp_.n >= 200)]
        if not len(cmp_):
            print(f"  {name}: no overlapping bins to compare")
            continue
        ratio = cmp_.p_here / cmp_.eff
        print(
            f"  {name}: {len(cmp_)} bins with eff>1e-3, "
            f"p_here/p_shipped median={np.median(ratio):.3f} "
            f"[p16={np.percentile(ratio, 16):.3f}, p84={np.percentile(ratio, 84):.3f}]"
        )
        show = cmp_.iloc[np.linspace(0, len(cmp_) - 1, min(8, len(cmp_))).astype(int)]
        for _, r in show.iterrows():
            print(
                f"      {magcol}={r[magcol]:.2f}  shipped={r.eff:.5f}  "
                f"here={r.p_here:.5f}  ratio={r.p_here / r.eff:.3f}  (N={int(r.n):,})"
            )


def cell_scan(sub, label, cells=(1.0, 0.5, 0.25, 0.1, 0.05)):
    """r_indep vs the conditioning cell size -- is the residual coupling or heterogeneity?

    ``r_indep`` compares the measured joint rate against independence imposed inside
    (Delta-m_lsst, Delta-m_roman) cells.  Any variable the cells do not resolve (size,
    surface brightness, blending) leaves within-cell covariance behind and inflates
    r_indep.  If r_indep falls steadily toward 1 as the cells shrink, the residual is
    unresolved heterogeneity, not a genuine coupling of the two draws.  If it plateaus
    above 1, the plateau is the real coupling.

    NOTE the code's own prediction uses per-object probabilities that depend on Delta-m
    ALONE, so the 0.25-mag row is the number that judges the code; the finer rows only
    diagnose WHERE the residual comes from.
    """
    print(f"\n----- cell-size scan ({label}) -----")
    print(f"  {'cell':>6} {'n_cells':>8} {'pred_indep':>12} {'r_indep':>9} {'r_comono':>9}")
    dm_l = sub["dm_lsst"].to_numpy()
    dm_r = sub["dm_roman"].to_numpy()
    sl = sub["star_lsst"].to_numpy()
    sr = sub["star_roman"].to_numpy()
    pj = float((sl & sr).mean())
    for c in cells:
        cp = cell_predictions(dm_l, dm_r, sl, sr, cell=c, fit_copula=False)
        ri = pj / cp["pred_indep"] if cp["pred_indep"] > 0 else np.nan
        rc = pj / cp["pred_comono"] if cp["pred_comono"] > 0 else np.nan
        print(
            f"  {c:6.2f} {cp['n_cells']:8d} {cp['pred_indep']:12.6g} "
            f"{ri:9.3f} {rc:9.3f}"
        )
    print(f"  (measured p_joint = {pj:.6g}, N = {len(sub):,})")


def make_figure(tab, meta):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.4), sharex=True)
    for j, (stype, title) in enumerate(
        [("galaxies", 'compact true galaxies (size$<$0.3")'), ("stars", "true stars")]
    ):
        t = tab[
            (tab.source_type == stype)
            & (tab.selection == "dm_bin")
            & (tab.dm_axis == "roman")
            & (tab.conditioning == "scalar")
        ].sort_values("dm_bin")
        ax, axr = axes[0, j], axes[1, j]
        if len(t):
            x = t.dm_bin.values
            ax.errorbar(
                x,
                t.p_lsst,
                yerr=[t.p_lsst - t.p_lsst_lo, t.p_lsst_hi - t.p_lsst],
                fmt="o-",
                ms=3.5,
                color="C0",
                label=r"$p_{\rm LSST}$ (extendedness$<$0.5)",
            )
            ax.errorbar(
                x,
                t.p_roman,
                yerr=[t.p_roman - t.p_roman_lo, t.p_roman_hi - t.p_roman],
                fmt="s-",
                ms=3.5,
                color="C3",
                label=r"$p_{\rm Roman}$ (F158 size envelope)",
            )
            ax.errorbar(
                x,
                t.p_joint,
                yerr=[t.p_joint - t.p_joint_lo, t.p_joint_hi - t.p_joint],
                fmt="D",
                ms=5,
                color="k",
                label=r"$p_{\rm joint}$ measured",
            )
            ax.plot(
                x, t.pred_indep, "--", color="tab:green", lw=2, label="independent null"
            )
            ax.plot(
                x,
                t.pred_comono,
                ":",
                color="tab:purple",
                lw=2.4,
                label=r"comonotonic $\min(p_1,p_2)$",
            )
            axr.axhline(1.0, color="0.6", lw=1)
            axr.plot(x, t.r_indep, "--o", color="tab:green", ms=3.5, label="vs independent")
            axr.plot(
                x, t.r_comono, ":s", color="tab:purple", ms=3.5, label="vs comonotonic"
            )
            axr.plot(
                x,
                t.r_indep_naive,
                "-",
                color="tab:green",
                alpha=0.35,
                lw=1.2,
                label="vs independent (naive)",
            )
        ax.set_title(f"{title}   (N={int(t.n.sum()):,} in bins)", fontsize=11)
        ax.set_ylabel("P(classified point-like)")
        ax.set_yscale("log")
        ax.set_ylim(1e-4, 1.6)
        ax.legend(fontsize=8, loc="lower right")
        ax.grid(alpha=0.25)
        axr.set_xlabel(r"$\Delta m = m^{\rm true}_{\rm F158} - {\rm maglim_{F158}}$")
        axr.set_ylabel(r"$p_{\rm joint}\,/\,$prediction")
        axr.set_yscale("log")
        axr.set_ylim(0.05, 60)
        axr.legend(fontsize=8, loc="upper left")
        axr.grid(alpha=0.25)
    fig.suptitle(
        "LSST$\\times$Roman classification coupling in DC2 (conditional on detection by both)\n"
        f"size source: {meta['size_source']}   |   "
        f"maglim$_{{F158}}$={meta['ref_maglim_f158']:.2f}, maglim$_r$={meta['ref_maglim_r']:.2f}",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT_PNG, dpi=140)
    print(f"wrote {OUT_PNG.relative_to(REPO)}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--size-source",
        default="cosmodc2_true",
        choices=["cosmodc2_true", "measured_roman"],
        help="compact-galaxy size: cosmoDC2 true size (default, matches the LSST "
        "product) or the interim measured Roman F158 size_sb",
    )
    ap.add_argument("--refit", action="store_true", help="refit the Roman classifier")
    ap.add_argument(
        "--cell-scan",
        action="store_true",
        help="also report r_indep vs the 2-D cell size, to separate a real coupling "
        "from leftover within-cell heterogeneity (see cell_scan)",
    )
    args = ap.parse_args()
    ART_DIR.mkdir(exist_ok=True)

    ml_f158, npix_f158 = ref_maglim(ROMAN_MAGLIM)
    ml_r, npix_r = ref_maglim(LSST_MAGLIM)
    print(f"ref maglim F158 = {ml_f158:.4f} ({npix_f158} px);  r = {ml_r:.4f} ({npix_r} px)")

    classify = get_classifier(refit=args.refit)
    df = load_matched()

    # ---- Roman measured size + per-survey classification --------------------------
    # size_sb() only indexes the three moment columns, so hand it a dict (no df copy)
    df["size_sb"] = size_sb(
        {
            "x2win_world_H158": df["roman_obs_x2win_world_H158"].to_numpy(),
            "y2win_world_H158": df["roman_obs_y2win_world_H158"].to_numpy(),
            "xywin_world_H158": df["roman_obs_xywin_world_H158"].to_numpy(),
        }
    )
    df["star_roman"] = classify(
        df["roman_obs_mag_auto_H158"].to_numpy(), df["size_sb"].to_numpy()
    )
    ext = df["lsst_obs_extendedness"].to_numpy()
    df["star_lsst"] = np.isfinite(ext) & (ext < EXT_CUT)

    # ---- shared quality base ------------------------------------------------------
    base = (
        df["roman_obs_matched"].astype(bool).to_numpy()  # drop spurious detections
        & (df["roman_obs_flags"].to_numpy() < FLAG_CUT)
        & df["lsst_obs_clean"].astype(bool).to_numpy()
        & np.isfinite(df["lsst_obs_extendedness"].to_numpy())
        & df["lsst_true_matched"].astype(bool).to_numpy()
        & np.isfinite(df["roman_true_mag_H158"].to_numpy())
        & np.isfinite(df["lsst_true_mag_r"].to_numpy())
        & np.isfinite(df["size_sb"].to_numpy())
    )
    print(f"quality base: {int(base.sum()):,} / {len(df):,}")


    df["dm_roman"] = df["roman_true_mag_H158"].to_numpy() - ml_f158
    df["dm_lsst"] = df["lsst_true_mag_r"].to_numpy() - ml_r
    df["mag_roman_true"] = df["roman_true_mag_H158"].to_numpy()

    # LOCAL per-pixel depth -- what the injector actually evaluates its curves at.
    # Positions come from the LSST object (the Roman detection is within 1", far
    # inside either depth map's pixel scale).
    _ra = df["lsst_obs_ra"].to_numpy()
    _dec = df["lsst_obs_dec"].to_numpy()
    ml_f158_loc = local_maglim(ROMAN_MAGLIM, _ra, _dec)
    ml_r_loc = local_maglim(LSST_MAGLIM, _ra, _dec)
    df["dm_roman_local"] = df["roman_true_mag_H158"].to_numpy() - ml_f158_loc
    df["dm_lsst_local"] = df["lsst_true_mag_r"].to_numpy() - ml_r_loc
    print(
        "local maglim: F158 "
        f"p16/50/84 = {np.nanpercentile(ml_f158_loc, [16, 50, 84]).round(3)} "
        f"(scalar {ml_f158:.3f}); r "
        f"p16/50/84 = {np.nanpercentile(ml_r_loc, [16, 50, 84]).round(3)} "
        f"(scalar {ml_r:.3f})"
    )
    _both = np.isfinite(ml_f158_loc) & np.isfinite(ml_r_loc)
    if _both.sum() > 10:
        print(
            "  spatial depth correlation between surveys: "
            f"r = {np.corrcoef(ml_f158_loc[_both], ml_r_loc[_both])[0, 1]:+.4f} "
            "(this is the term the scalar conditioning leaves in the cells)"
        )

    # The local-maglim pass needs a covered pixel in BOTH depth maps. Require it in the
    # shared base so the scalar and local passes run on the SAME objects -- otherwise
    # the two r_indep values would differ partly through sample choice, not conditioning.
    _finite_local = np.isfinite(df["dm_roman_local"].to_numpy()) & np.isfinite(
        df["dm_lsst_local"].to_numpy()
    )
    n_lost = int((base & ~_finite_local).sum())
    base = base & _finite_local
    print(
        f"  dropped {n_lost:,} on pixels uncovered by one of the two depth maps "
        f"-> {int(base.sum()):,} in the common base"
    )

    # ---- truth labels: require the two truth catalogs to AGREE --------------------
    rgs = df["roman_true_gal_star"].to_numpy()
    ltt = df["lsst_true_truth_type"].to_numpy()
    is_gal = (rgs == 0) & (ltt == TRUTH_GAL)
    is_star = (rgs == 1) & (ltt == TRUTH_STAR)
    print(
        f"truth-agreeing galaxies: {int((base & is_gal).sum()):,}   "
        f"stars: {int((base & is_star).sum()):,}   "
        f"disagreements dropped: {int((base & ~is_gal & ~is_star).sum()):,}"
    )

    # ---- compact-galaxy size (both definitions; one is chosen as primary) ---------
    size_true = join_true_sizes(df)
    size_meas = df["size_sb"].to_numpy()
    compact_true = np.isfinite(size_true) & (size_true < GAL_SIZE_MAX)
    compact_meas = np.isfinite(size_meas) & (size_meas < GAL_SIZE_MAX)
    if args.size_source == "cosmodc2_true":
        compact = compact_true
        size_desc = f'cosmoDC2 size_true < {GAL_SIZE_MAX}" (joined on lsst_true_cosmodc2_id)'
    else:
        compact = compact_meas
        size_desc = (
            f'measured Roman F158 size_sb < {GAL_SIZE_MAX}" (INTERIM; conflates with PSF)'
        )

    # ---- marginal validation against the shipped curves --------------------------
    gal_base = base & is_gal
    validate_marginals(
        df["roman_true_mag_H158"].to_numpy()[gal_base],
        df["star_roman"].to_numpy()[gal_base],
        compact_meas[gal_base],
        df["lsst_true_mag_r"].to_numpy()[gal_base],
        df["star_lsst"].to_numpy()[gal_base],
        compact_true[gal_base],
    )

    keep = [
        "dm_roman",
        "dm_lsst",
        "dm_roman_local",
        "dm_lsst_local",
        "star_lsst",
        "star_roman",
        "mag_roman_true",
        "roman_true_ind",
    ]
    gal = df.loc[base & is_gal & compact, keep].drop_duplicates("roman_true_ind")
    star = df.loc[base & is_star, keep].drop_duplicates("roman_true_ind")
    print(f"compact true galaxies: {len(gal):,}   true stars: {len(star):,}")

    t0 = time.time()
    tab = pd.concat(
        [
            measure(sub, lab, size_desc, axis=ax, conditioning=cond)
            for lab, sub in [("galaxies", gal), ("stars", star)]
            for ax in ("roman", "lsst")
            for cond in ("scalar", "local")
        ],
        ignore_index=True,
    )
    print(f"binned measurement: {time.time() - t0:.0f}s")

    if args.cell_scan:
        cell_scan(gal, "compact galaxies, all")
        cell_scan(
            gal[gal["mag_roman_true"] < JOINT_F158_CUT], "compact galaxies, F158<25.5"
        )
        cell_scan(star, "stars, all")

    tab.to_csv(OUT_CSV, index=False, float_format="%.6g")
    print(f"wrote {OUT_CSV.relative_to(REPO)} ({len(tab)} rows)")

    meta = dict(
        script=str(Path(__file__).relative_to(REPO)),
        matched_table=str(MATCHED),
        n_matched_rows=int(len(df)),
        n_compact_galaxies=int(len(gal)),
        n_stars=int(len(star)),
        size_source=size_desc,
        ref_maglim_f158=ml_f158,
        ref_maglim_r=ml_r,
        roman_classifier="F158 size envelope (roman_star_classifier.build_env_classifier, ENV_PURITY=0.875)",
        lsst_classifier=f"extendedness < {EXT_CUT}",
        conditioning="both surveys DETECT the object (matched table is det<->obj 1'' matched)",
        detection_coupling_measured=False,
        cell_size_mag=CELL,
        joint_f158_cut=JOINT_F158_CUT,
    )
    OUT_META.write_text(json.dumps(meta, indent=2))
    make_figure(tab, meta)

    # ---- printed summary ----------------------------------------------------------
    pd.set_option("display.width", 200, "display.max_columns", 50)
    for stype in ("galaxies", "stars"):
        t = tab[tab.source_type == stype]
        print(f"\n===== {stype} =====")
        cols = [
            "selection",
            "dm_axis",
            "dm_bin",
            "n",
            "p_lsst",
            "p_roman",
            "p_joint",
            "pred_indep",
            "pred_comono",
            "r_indep",
            "r_comono",
            "rho_copula",
        ]
        print(t[cols].to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    # ---- the verdict table: scalar (upper bound) vs local (like-for-like) --------
    print("\n===== conditioning comparison: r_indep vs what the code actually does =====")
    print(
        "  scalar = Delta-m from the median maglim (what the shipped curves are keyed\n"
        "  on; leaves correlated spatial depth inside the cells -> UPPER BOUND).\n"
        "  local  = Delta-m from each object's own pixel depth (what the injector\n"
        "  evaluates -> like-for-like)."
    )
    piv = tab[tab.selection.isin(["ALL", "F158_lt_25.5"])]
    for stype in ("galaxies", "stars"):
        for sel in ("ALL", "F158_lt_25.5"):
            row = piv[(piv.source_type == stype) & (piv.selection == sel)]
            got = {r.conditioning: r for r in row.itertuples()}
            if "scalar" in got and "local" in got:
                sc, lo = got["scalar"], got["local"]
                print(
                    f"  {stype:<9} {sel:<13} r_indep: scalar={sc.r_indep:.3f} -> "
                    f"local={lo.r_indep:.3f}   r_comono: {sc.r_comono:.3f} -> "
                    f"{lo.r_comono:.3f}   rho: {sc.rho_copula:.3f} -> {lo.rho_copula:.3f}"
                )

    print(
        "\nNOTE: this measures the CLASSIFICATION coupling conditional on detection by "
        "both surveys.\nThe DETECTION coupling (objects one survey detects and the other "
        "misses) is not measurable\nfrom this det-centric table and needs a truth-centric "
        "join -- separate step."
    )


if __name__ == "__main__":
    main()
