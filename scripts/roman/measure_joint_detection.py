#!/usr/bin/env python
"""Measure the LSST<->Roman DETECTION correlation in DC2 (task T3, Part B).

WHY
---
This is the companion to ``measure_joint_misclassification.py`` (read that script
first -- this one imports its statistics helpers and mirrors its output schema).
That script measured the CLASSIFICATION coupling conditional on an object being
DETECTED BY BOTH surveys, because its input (``roman_lsst_matched.parquet``) is
built from Roman detections positionally matched to LSST objects -- by
construction it cannot contain an object that only one survey detects.

This script fills the missing piece: the DETECTION coupling. It builds a
TRUTH-CENTRIC sample (every true source in the Roman x LSST DC2 overlap, whether
or not either survey detected it) and measures

    p_lsst   = P(LSST detects)
    p_roman  = P(Roman detects)
    p_joint  = P(both detect)

and the same for DETECTION x CLASSIFICATION (the quantity the injector actually
multiplies, ``Survey.get_gal_misclassification_detection`` -> ``mis(dm)*det(dm)``,
``streamobs/surveys.py:543-556``), since a joint contaminant/complete-sample
estimate needs the full per-object selection probability, not detection alone.

Like the classification script, the injector draws an INDEPENDENT uniform per
survey (``StreamInjector`` spawns a child RNG per survey; ``detect_flag`` at
``streamobs/observed.py:1237`` draws its own uniform) and the joint analysis
channel ANDs the two flags (``fig4_metricC2_ingredients._cmd_cols``, ``which ==
"joint"``), so the joint selection probability is ``p_lsst * p_roman`` by
construction. Whether that is right is an empirical question about how strongly
correlated the two surveys' detections are -- e.g. a source dropped by LSST
because of a chance blend with a foreground galaxy is not obviously any less
likely to be detected by Roman, so DETECTION coupling need not be as strong as
the CLASSIFICATION coupling (which is driven by the same physical quantity,
angular size, in both surveys). That is exactly what this script measures.

WHAT IS AND IS NOT MEASURED HERE
---------------------------------
Every row in the sample is a TRUE source (star or compact galaxy) that exists in
BOTH surveys' input truth catalogs (Roman-Rubin DC2 shares one underlying
cosmoDC2/star truth, see the cross-match section below), regardless of whether
either survey detected it. This is therefore the DETECTION-coupling analog of
the CLASSIFICATION-coupling script's ``p_lsst``/``p_roman``/``p_joint``, plus the
combined DETECTION x CLASSIFICATION probability that is what the injector
actually draws against. It is NOT a remeasurement of the classification-only
coupling (already done in ``measure_joint_misclassification.py``) -- the
classification flag used here is unconditional on detection (``P(detected AND
classified point-like)``), not ``P(classified point-like | detected)``.

THE HETEROGENEITY TRAP -- SAME AS THE CLASSIFICATION SCRIPT
-------------------------------------------------------------
Within a Delta-m bin, objects are not identical: a galaxy bright in F158 is
usually also bright in r, so per-object detection probabilities p_l(i), p_r(i)
are themselves correlated across the bin, which inflates the naive product
``p_lsst * p_roman`` relative to the correct within-bin-independent null even
when the code's independent draws are exactly right. The primary null used here
is therefore ``cell_predictions()``'s cell-weighted independent/comonotonic
prediction (2-D cells of 0.25 mag in (Delta-m_lsst, Delta-m_roman)), imported
verbatim from ``measure_joint_misclassification.py`` -- see that module's
docstring for the derivation. Both cell-weighted and naive ratios are reported;
``r_indep``/``r_comono`` are the cell-weighted ones.

INTERPRETATION RULE (same as the handoff for Part A)
------------------------------------------------------
* ``r_indep  ~ 1`` -> the current independent draws are correct for detection.
* ``r_comono ~ 1`` -> detection should be drawn from a shared uniform.
* in between   -> ``rho_copula`` (fitted per source type/bin) is the actionable number.

REUSED ESTIMATORS (imported, not reimplemented)
-------------------------------------------------
``wilson``, ``cell_predictions``, ``fit_rho``, ``gauss_copula_joint``, ``bvn_cdf``,
``get_classifier``, ``size_sb``, ``ref_maglim`` are imported directly from
``measure_joint_misclassification.py`` (added to ``sys.path``). ``get_classifier``
loads the SAME cached F158 size-envelope boundary grid the classification script
built (``data/surveys/roman_dc2/_cache_env_classifier_bounds.npz``), so the two
scripts cannot use different Roman star/galaxy classifiers.

INPUTS AND EXACT CONVENTIONS COPIED (source of each)
-------------------------------------------------------
* Roman truth, per tile: ``truth/dc2_index_<tile>.fits.gz`` under the roman_mock
  data root (see ``build_roman_dc2_det_truth.py``). Columns used: ``ind, ra, dec,
  gal_star, mag_H158`` (``mag_H158 == 0.0`` -> NaN, "no flux in band" convention
  from that script). The tile universe (1039 real tiles) is discovered the same
  way as that script's ``discover_tiles()``: a non-empty ``dc2_index_*`` truth
  file with a matching ``dc2_det_*`` detection file.
* Roman DETECTED flag: a truth object (keyed on ``(tile, truth_ind)``, since
  ``truth_ind`` repeats across tiles -- see ``roman_dc2_det_truth.README.md``) is
  Roman-detected iff it appears as ``truth_ind`` on a ``matched == True`` row of
  ``roman_dc2_det_truth.parquet``. S/N>5 and the 1" det->truth match are already
  baked into that table's construction (``build_roman_dc2_det_truth.py``). TWO
  parallel definitions are carried all the way through the measurement (all four
  appear in the output CSV's ``quantity`` column):
    - ``detection`` / ``detection_classification``: bare ``matched`` (the
      handoff's definition for this task -- NO ``flags`` cut).
    - ``detection_flagcut`` / ``detection_classification_flagcut``: ALSO
      requires ``flags < FLAG_CUT`` (``flags == 0``), matching what the SHIPPED
      ``roman_stellar_efficiency_cutf158.csv`` / ``build_roman_galaxy_misclass.py``
      apply (added per coordinator request, to check whether the flags cut moves
      ``r_indep``/``r_comono``/``rho`` -- it removes ~14% of matched detections,
      see the printed summary). Where a truth object has more than one matched
      detection (shredding), the NEAREST one PASSING THE RELEVANT MASK is kept
      (the flags cut is applied before picking the representative detection, not
      after -- ``dedup_nearest_lookup``), mirroring
      ``create_streamobs_files_hlwas.py``'s / ``build_roman_galaxy_misclass.py``'s
      own order of operations.
* Roman classification: F158 size-envelope classifier via ``get_classifier()``
  (``roman_star_classifier.build_env_classifier``), applied to the kept
  detection's ``mag_auto_H158`` and ``size_sb(x2/y2/xywin_world_H158)`` --
  identical inputs to ``measure_joint_misclassification.py``.
* LSST truth, per tract: ``dc2_run2.2i_truth_merged_summary_skim_tract_<tract>.fits``.
  Columns: ``ra, dec, truth_type`` (1=galaxy, 2=star, 3=SN, SNe dropped), ``id,
  cosmodc2_id, mag_r``.
* LSST DETECTED + classified flag: copied EXACTLY from
  ``scripts/lsst/build_lsst_dc2_galaxy_misclass.py`` /
  ``scripts/lsst/create_streamobs_files_lsst_dc2.py``: match every LSST OBJECT to
  its nearest TRUTH source (``dc2_object_run2.2i_dr6_skim_tract_<tract>.fits``,
  object -> truth direction, NOT truth -> object), ``MATCH_RADIUS = 1.0"``,
  ``SNR_DEPTH = 5`` on ``psFlux_r/psFluxErr_r``, require ``clean``; where several
  objects map to the same truth (deblends), keep the nearest
  (``drop_duplicates("truth_id")`` after sorting on separation). A truth object is
  LSST-detected iff it survives that cut; LSST-classified iff its kept object has
  ``extendedness < EXT_CUT = 0.5``.
* Compact-galaxy size: ``cosmodc2_galaxy_size_true.parquet`` (``cosmodc2_id,
  size_true``), joined on the LSST truth ``cosmodc2_id`` -- identical join to
  ``measure_joint_misclassification.py``'s ``join_true_sizes``.
  ``GAL_SIZE_MAX = 0.3"``.
* Delta-m convention: ``delta_m = truth_mag - REF_MAGLIM``, ``REF_MAGLIM`` = the
  median of the survey's maglim map over covered pixels, via the imported
  ``ref_maglim()`` (``roman_dc2_maglim_f158_nside1024.fits.gz``,
  ``lsst_dc2_maglim_r_nside1024.fits.gz``) -- same maps, same convention as Part A.
* Truth-label agreement: Roman ``gal_star == 1`` & LSST ``truth_type == 2`` for
  stars; Roman ``gal_star == 0`` & LSST ``truth_type == 1`` for galaxies.
  Disagreements are dropped and counted (imported ``TRUTH_GAL``/``TRUTH_STAR``,
  ``EXT_CUT``, ``GAL_SIZE_MAX``, ``JOINT_F158_CUT`` constants -- same values as
  Part A, module-level import).

CROSS-SURVEY TRUTH MATCH (the piece Part A did not need)
------------------------------------------------------------
Roman-Rubin DC2 shares one input truth catalog with the wider LSST DC2 run2.2i
sim (Troxel et al. 2023), so a genuine counterpart should coincide to well under
1". Matching is nearest-neighbor by exact great-circle separation (haversine),
using a scipy ``cKDTree`` on an equirectangular projection ONLY to build the
spatial index (a local flat-sky approximation good to << 1" over the several-
degree field; the separation actually reported and cut on is the exact haversine
value, not the projected one). The separation histogram is printed and the match
radius is chosen from it (see ``main()``); a generous LSST-tract bbox around the
Roman footprint discovers the 15 overlapping LSST tracts automatically.

MANDATORY CAVEATS FOR THE READER
------------------------------------
* This sample is restricted to the observed OVERLAP of the two truth catalogs
  after the cross-match radius cut -- Roman truth objects with no LSST-truth
  counterpart (or vice versa) within the match radius are dropped, which also
  drops the rare cases where the truth catalogs themselves disagree (e.g. one
  side lacks a source the other has). The effective area and object counts are
  printed and stored in the ``.meta.json``.
* The bare ``matched`` Roman-detected definition (``detection`` /
  ``detection_classification`` in the ``quantity`` column) is DIAGNOSTIC ONLY.
  ``create_streamobs_files_hlwas.py`` (section 2: "... and passed the paper's
  `flags == 0` cut") and ``build_roman_galaxy_misclass.py`` (``clean =
  merged["roman_obs_flags"].values < FLAG_CUT``) both apply ``flags < FLAG_CUT``
  to define Roman detection/classification in the SHIPPED selection-function
  products, so the ``*_flagcut`` quantities are AUTHORITATIVE -- use those, not
  the bare ones, when quoting a number. The two are NOT interchangeable: for
  STARS the flags cut moves ``r_indep`` by ~1-2% (it removes a similar fraction
  of detections regardless of classification outcome, so it largely cancels in
  the ratio), but for GALAXIES it is NOT negligible -- e.g. compact galaxies,
  ``detection_classification``, ``F158<25.5``: ``r_indep`` goes from 0.9995
  (bare) to 1.52 (flagcut), because ``flags<FLAG_CUT`` removes a much larger
  share of the (rare) misclassified-as-point-source detections than of
  detections overall. The bare numbers are kept in the output only to show that
  this movement is itself part of the finding.
* No S/N or depth cut is applied beyond what is already baked into each survey's
  own detection catalogs -- this measures detection AS THE INJECTOR MODELS IT,
  not some idealized common depth.
* Delta-m here is measured from the SCALAR median maglim (``ref_maglim()``),
  while the injector evaluates its efficiency curves at each object's LOCAL
  per-pixel maglim, which raised a concern that correlated inter-survey depth
  variation could sit inside this script's Delta-m cells as unmodelled shared
  heterogeneity and inflate ``r_indep``. CHECKED, and found NEGLIGIBLE: the
  coordinator measured the inter-survey spatial depth correlation directly
  (each object's own per-pixel maglim from both depth maps, same maps and
  ~ the same 17 sq deg footprint as this sample) and found ``r = -0.040``
  (i.e. no meaningful shared-depth term -- depth is uniform to ~0.06 mag,
  p16-p84, across this footprint in both bands) and confirmed the galaxy
  ``r_indep`` is unchanged (1.273 -> 1.273) between the scalar- and local-maglim
  conditioning. The ``r_indep``/``r_comono`` values reported here are therefore
  POINT ESTIMATES of the injector's independence error, not upper bounds.
* The compact-GALAXY sample requires a cosmoDC2 ``size_true`` join, which only
  covers a fraction of the cross-matched galaxy footprint (this is the ONLY
  compactness measure available for a galaxy neither survey detects, so a
  measured-size proxy cannot substitute here). ``report_size_coverage()`` prints
  and stores (``compact_galaxy_size_coverage`` in the ``.meta.json``) the covered
  vs. total RA/Dec range, effective healpix area, and LSST-/Roman-detected rates
  on the covered subset vs. the full cross-matched galaxy sample, so the reader
  can judge whether the compact-galaxy result generalizes beyond the covered
  region or should be scoped to it. (Found during review: the FIRST version of
  this join used ``data/surveys/roman_dc2/cosmodc2_galaxy_size_true.parquet``,
  which turned out to be DETECTION-biased -- see the note above
  ``join_true_sizes_full`` -- and was replaced by a join against the full,
  unfiltered per-healpix cosmoDC2 skims; the representativeness check now
  PASSES, see ``.meta.json``.)
* The paper's joint background does NOT restrict to compact galaxies
  (``fig4_metricC2_ingredients.py``: ``gal_idx = np.where(~is_star)[0]`` is
  every galaxy, fed through ``inject_galaxies`` -> the same
  ``get_gal_misclassification_detection`` used for compact ones). An
  ``all_galaxies`` source type (identical machinery, no ``size_true<0.3"`` cut)
  is therefore reported alongside ``galaxies`` -- ``all_galaxies`` is the
  quantity that actually feeds the paper's contaminant density; ``galaxies``
  (compact-only) is kept for comparison, since compactness is plausibly itself
  correlated with the detection coupling.
* ``create_streamobs_files_hlwas.py`` notes that ~26% of Roman truth STARS have
  a duplicate entry at the exact same position (different ``ind``, only the
  J129 band populated) and collapses truth stars to unique POSITIONS before
  computing the shipped stellar-efficiency curve (a star counts as detected if
  ANY member entry matched). This script does NOT do that collapse, and its
  finite-Delta-m filter would drop a star if the Roman<->LSST cross-match lands
  on the J129-only duplicate (H158 NaN) rather than the full-band one.
  CHECKED with ``check_star_position_collapse()`` (run on the full Roman truth
  catalog, independent of the LSST cross-match) -- two different answers at two
  different levels, both reported in ``star_duplicate_position_check`` in the
  ``.meta.json``:
    - On the RAW Roman truth catalog (230,820 star rows, 15.4% at a duplicate
      position -- close to, not identical to, the ~26% figure, likely because
      that number was measured on the full HLWAS footprint and this is the
      ~17 sq deg DC2 overlap only), collapsing by position is NOT negligible:
      ``p_roman`` (bare) rises from 0.806 to 0.952 and (flags-cut) from 0.759
      to 0.897 -- an ~18% relative shift if you average per ROW instead of per
      unique star.
    - But THIS SCRIPT's actual analysis sample is not "all Roman truth star
      rows" -- it is one Roman truth row per LSST truth star, chosen by nearest-
      neighbor cross-match, so the raw-catalog duplication rate does not
      directly transfer. The marginal validation against
      ``roman_stellar_efficiency_cutf158.csv`` (itself built WITH the
      position-collapse) lands at ratio ~1.006 (median) for this script's
      actual star sample, which is the more direct evidence for this script's
      numbers specifically. The two checks are not fully reconciled here (that
      would need tracing exactly which member of a duplicate pair the
      cross-match selects, not done for time); readers should treat the star
      ``r_indep``/``r_comono`` values as validated empirically to ~1%, without
      a first-principles explanation for why the raw-catalog duplication effect
      does not visibly propagate into them.

Run (streamobs env):
    /astro/store/shiren/conda-envs/stream_team/envs/streamobs/bin/python \
        scripts/roman/measure_joint_detection.py [--refresh] [--tracts N] [--cross-radius X]
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time
from pathlib import Path

import fitsio
import healpy as hp
import numpy as np
import pandas as pd
import pyarrow.dataset as pads
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parents[2]
ROMAN_DIR = REPO / "data" / "surveys" / "roman_dc2"
LSST_DIR = REPO / "data" / "surveys" / "lsst_dc2"
ART_DIR = REPO / "artifacts"

ROMAN_MOCK_ROOT = Path("/astro/store/shire/stream_team/stream_finding/data/roman_mock")
ROMAN_TRUTH_DIR = ROMAN_MOCK_ROOT / "truth"
ROMAN_DET_DIR = ROMAN_MOCK_ROOT / "det"
LSST_RAW_DIR = Path("/astro/store/shire/stream_team/stream_finding/data/lsst_dc2")

DET_TRUTH = ROMAN_DIR / "roman_dc2_det_truth.parquet"
COSMODC2_SIZES = ROMAN_DIR / "cosmodc2_galaxy_size_true.parquet"
ROMAN_MAGLIM = ROMAN_DIR / "roman_dc2_maglim_f158_nside1024.fits.gz"
LSST_MAGLIM = LSST_DIR / "lsst_dc2_maglim_r_nside1024.fits.gz"

CACHE_ROMAN_TRUTH = ROMAN_DIR / "_cache_roman_truth_all.parquet"
CACHE_ROMAN_DETLOOKUP = ROMAN_DIR / "_cache_roman_detected_raw.parquet"
CACHE_LSST_TRUTH = LSST_DIR / "_cache_lsst_truth_detflag_overlap.parquet"

OUT_CSV = ART_DIR / "joint_detection_correlation.csv"
OUT_PNG = ART_DIR / "joint_detection_correlation.png"
OUT_META = ART_DIR / "joint_detection_correlation.meta.json"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from measure_joint_misclassification import (  # noqa: E402
    CELL,
    DM_BINS,
    DM_MID,
    EXT_CUT,
    FLAG_CUT,
    GAL_SIZE_MAX,
    JOINT_F158_CUT,
    N_MIN,
    TRUTH_GAL,
    TRUTH_STAR,
    bvn_cdf,  # noqa: F401  (re-exported for anyone importing this module)
    cell_predictions,
    fit_rho,  # noqa: F401
    gauss_copula_joint,  # noqa: F401
    get_classifier,
    ref_maglim,
    size_sb,
    wilson,
)

# --- conventions copied from the LSST product builders (see module docstring) - #
LSST_MATCH_RADIUS = 1.0  # arcsec, LSST object -> nearest truth
SNR_DEPTH = 5  # LSST psFlux_r / psFluxErr_r cut

# --- Roman truth "no flux" sentinel (build_roman_dc2_det_truth.py convention) -- #
ROMAN_MAG_NODATA = 0.0

ARCSEC_PER_DEG = 3600.0
ROMAN_FOOTPRINT = dict(ra_lo=50.5, ra_hi=56.5, dec_lo=-42.3, dec_hi=-37.6)  # generous


# --------------------------------------------------------------------------- #
# geometry helpers
# --------------------------------------------------------------------------- #
def haversine_arcsec(ra1, dec1, ra2, dec2):
    """Exact great-circle separation (arcsec), vectorized."""
    ra1r, dec1r, ra2r, dec2r = (np.radians(np.asarray(a, float)) for a in (ra1, dec1, ra2, dec2))
    dra = ra2r - ra1r
    ddec = dec2r - dec1r
    a = np.sin(ddec / 2) ** 2 + np.cos(dec1r) * np.cos(dec2r) * np.sin(dra / 2) ** 2
    return 2.0 * np.arcsin(np.sqrt(np.clip(a, 0, 1))) * (180.0 / np.pi) * ARCSEC_PER_DEG


def nearest_match(ra_q, dec_q, ra_ref, dec_ref, tangent=None):
    """Nearest ref point for each query point.

    The ``cKDTree`` is built on a local equirectangular projection (flat-sky,
    good to << 1" of distortion over a several-degree field) purely to make the
    nearest-neighbor SEARCH fast; the separation returned is the EXACT
    great-circle (haversine) distance for the pair the tree finds, so reported
    separations carry no projection error.
    """
    ra_ref = np.asarray(ra_ref, float)
    dec_ref = np.asarray(dec_ref, float)
    ra_q = np.asarray(ra_q, float)
    dec_q = np.asarray(dec_q, float)
    if tangent is None:
        tangent = (float(np.median(ra_ref)), float(np.median(dec_ref)))
    ra0, dec0 = tangent
    cosd = np.cos(np.radians(dec0))

    def _proj(ra, dec):
        return np.column_stack(
            [(ra - ra0) * cosd * ARCSEC_PER_DEG, (dec - dec0) * ARCSEC_PER_DEG]
        )

    tree = cKDTree(_proj(ra_ref, dec_ref))
    _, idx = tree.query(_proj(ra_q, dec_q), k=1, workers=-1)
    sep = haversine_arcsec(ra_q, dec_q, ra_ref[idx], dec_ref[idx])
    return idx, sep


def _native(a):
    """FITS arrays are big-endian; fitsio/pandas want native byte order."""
    if a.dtype.byteorder in (">", "="):
        try:
            return a.astype(a.dtype.newbyteorder("="))
        except Exception:
            return a
    return a


# --------------------------------------------------------------------------- #
# Roman truth-centric catalog: all truth objects + Roman-detected/classified flag
# --------------------------------------------------------------------------- #
def discover_roman_tiles():
    """Same recipe as build_roman_dc2_det_truth.py's discover_tiles(): a non-empty
    per-tile truth file with a matching detection file."""
    tiles = []
    for p in sorted(glob.glob(str(ROMAN_TRUTH_DIR / "dc2_index_*.fits.gz"))):
        name = Path(p).name[len("dc2_index_") : -len(".fits.gz")]
        parts = name.split("_")
        if len(parts) != 2:
            continue
        try:
            float(parts[0])
            float(parts[1])
        except ValueError:
            continue
        if os.path.getsize(p) <= 1000:
            continue
        if not (ROMAN_DET_DIR / f"dc2_det_{name}.fits.gz").exists():
            continue
        tiles.append(name)
    return sorted(tiles)


def load_all_roman_truth(refresh=False, limit=0):
    if CACHE_ROMAN_TRUTH.exists() and not refresh and not limit:
        print(f"loading cached Roman truth {CACHE_ROMAN_TRUTH.name}")
        return pd.read_parquet(CACHE_ROMAN_TRUTH)
    tiles = discover_roman_tiles()
    if limit:
        tiles = tiles[:limit]
    print(f"reading {len(tiles)} Roman truth tiles ...")
    t0 = time.time()
    frames = []
    for name in tiles:
        p = ROMAN_TRUTH_DIR / f"dc2_index_{name}.fits.gz"
        d = fitsio.read(str(p), columns=["ind", "ra", "dec", "gal_star", "mag_H158"])
        mag = _native(d["mag_H158"]).astype("f8")
        mag = np.where(mag == ROMAN_MAG_NODATA, np.nan, mag)
        frames.append(
            pd.DataFrame(
                {
                    "tile": name,
                    "ind": _native(d["ind"]).astype("i8"),
                    "ra": _native(d["ra"]).astype("f8"),
                    "dec": _native(d["dec"]).astype("f8"),
                    "gal_star": _native(d["gal_star"]).astype("i2"),
                    "mag_H158": mag,
                }
            )
        )
    out = pd.concat(frames, ignore_index=True)
    print(f"  {len(out):,} Roman truth rows across {len(tiles)} tiles ({time.time()-t0:.0f}s)")
    if not limit:
        out.to_parquet(CACHE_ROMAN_TRUTH)
    return out


def load_roman_detected_raw(refresh=False):
    """One row per Roman DETECTION matched to truth (S/N>5 and the 1"
    det->truth match are already baked into ``roman_dc2_det_truth.parquet``).
    NOT deduplicated -- ``dedup_nearest_lookup`` does that, optionally after a
    ``flags`` cut, since which detection is "nearest" can depend on whether
    flagged detections are allowed (see ``build_roman_truth_catalog``)."""
    if CACHE_ROMAN_DETLOOKUP.exists() and not refresh:
        print(f"loading cached Roman matched-detections {CACHE_ROMAN_DETLOOKUP.name}")
        return pd.read_parquet(CACHE_ROMAN_DETLOOKUP)
    print(f"reading {DET_TRUTH.name} (matched==True rows only) ...")
    t0 = time.time()
    dset = pads.dataset(str(DET_TRUTH), format="parquet")
    cols = [
        "tile",
        "truth_ind",
        "match_sep_arcsec",
        "flags",
        "mag_auto_H158",
        "magerr_auto_H158",
        "x2win_world_H158",
        "y2win_world_H158",
        "xywin_world_H158",
    ]
    df = dset.to_table(columns=cols, filter=pads.field("matched") == True).to_pandas()  # noqa: E712
    print(f"  {len(df):,} matched detection rows ({time.time()-t0:.0f}s)")
    df["truth_ind"] = np.round(df["truth_ind"].to_numpy()).astype("i8")
    df.to_parquet(CACHE_ROMAN_DETLOOKUP)
    return df


def dedup_nearest_lookup(df, flag_ok_mask=None):
    """One row per (tile, truth_ind): the NEAREST matched detection, optionally
    restricted to a quality mask first (so "nearest" is computed among only the
    detections that pass the mask -- e.g. ``flags < FLAG_CUT`` -- matching how
    ``create_streamobs_files_hlwas.py`` / ``build_roman_galaxy_misclass.py``
    apply the flags cut BEFORE picking the representative detection, not after)."""
    d = df if flag_ok_mask is None else df.loc[flag_ok_mask]
    return (
        d.sort_values("match_sep_arcsec")
        .drop_duplicates(["tile", "truth_ind"], keep="first")
        .reset_index(drop=True)
    )


def build_roman_truth_catalog(refresh=False, limit=0):
    """All Roman truth objects + two parallel detected/classified flag pairs:

    * ``roman_detected`` / ``roman_star_class`` -- bare ``matched`` (the
      handoff's definition for this task; NO ``flags`` cut).
    * ``roman_detected_flagcut`` / ``roman_star_class_flagcut`` -- ALSO
      requires ``flags < FLAG_CUT`` (the cut the shipped Roman selection-
      function products apply), per the coordinator's request to check whether
      it moves ``r_indep``/``r_comono``/``rho`` materially.
    """
    truth = load_all_roman_truth(refresh=refresh, limit=limit)
    raw = load_roman_detected_raw(refresh=refresh)
    flag_ok = raw["flags"].to_numpy() < FLAG_CUT
    print(
        f"Roman matched detections: {len(raw):,} total, {int(flag_ok.sum()):,} pass "
        f"flags<{FLAG_CUT} ({flag_ok.mean():.1%}) -- {len(raw) - int(flag_ok.sum()):,} removed"
    )

    classify = get_classifier(refit=False)

    def _attach(det, suffix):
        det = det.rename(columns={"truth_ind": "ind"})
        cols = ["tile", "ind", "mag_auto_H158", "x2win_world_H158", "y2win_world_H158", "xywin_world_H158"]
        m = truth.merge(det[cols], on=["tile", "ind"], how="left", indicator=True)
        detected = (m["_merge"] == "both").to_numpy()
        sz = size_sb(
            {
                "x2win_world_H158": m["x2win_world_H158"].to_numpy(),
                "y2win_world_H158": m["y2win_world_H158"].to_numpy(),
                "xywin_world_H158": m["xywin_world_H158"].to_numpy(),
            }
        )
        star_class = np.zeros(len(m), bool)
        star_class[detected] = classify(m.loc[detected, "mag_auto_H158"].to_numpy(), sz[detected])
        print(
            f"  [{suffix}] {int(detected.sum()):,} / {len(m):,} detected "
            f"({detected.mean():.2%}); {int(star_class.sum()):,} classified point-like"
        )
        return detected, star_class

    merged = truth.copy()
    merged["roman_detected"], merged["roman_star_class"] = _attach(
        dedup_nearest_lookup(raw), "bare matched"
    )
    merged["roman_detected_flagcut"], merged["roman_star_class_flagcut"] = _attach(
        dedup_nearest_lookup(raw, flag_ok), f"flags<{FLAG_CUT}"
    )
    return merged[
        [
            "tile", "ind", "ra", "dec", "gal_star", "mag_H158",
            "roman_detected", "roman_star_class",
            "roman_detected_flagcut", "roman_star_class_flagcut",
        ]
    ]


# --------------------------------------------------------------------------- #
# LSST truth-centric catalog: all truth objects + LSST-detected/classified flag
# --------------------------------------------------------------------------- #
def discover_overlap_tracts(pad_deg=0.3):
    """LSST tracts whose truth-file footprint overlaps the (padded) Roman bbox."""
    fp = ROMAN_FOOTPRINT
    ra_lo, ra_hi = fp["ra_lo"] - pad_deg, fp["ra_hi"] + pad_deg
    dec_lo, dec_hi = fp["dec_lo"] - pad_deg, fp["dec_hi"] + pad_deg
    tracts = []
    for f in sorted(glob.glob(str(LSST_RAW_DIR / "dc2_run2.2i_truth_merged_summary_skim_tract_*.fits"))):
        tract = Path(f).name.split("_")[-1].split(".")[0]
        d = fitsio.read(f, columns=["ra", "dec"])
        r_lo, r_hi = _native(d["ra"]).min(), _native(d["ra"]).max()
        de_lo, de_hi = _native(d["dec"]).min(), _native(d["dec"]).max()
        if not (r_hi < ra_lo or r_lo > ra_hi or de_hi < dec_lo or de_lo > dec_hi):
            tracts.append(tract)
    return tracts


def process_lsst_tract(tract):
    """Object->truth match (LSST convention, see module docstring); returns one
    row per truth object (galaxy or star) with lsst_detected / lsst_star_class."""
    tf = LSST_RAW_DIR / f"dc2_run2.2i_truth_merged_summary_skim_tract_{tract}.fits"
    of = LSST_RAW_DIR / f"dc2_object_run2.2i_dr6_skim_tract_{tract}.fits"
    tru = fitsio.read(
        str(tf), columns=["ra", "dec", "truth_type", "id", "cosmodc2_id", "cosmodc2_hp", "mag_r"]
    )
    obj = fitsio.read(str(of), columns=["ra", "dec", "psFlux_r", "psFluxErr_r", "extendedness", "clean"])

    tru_type = _native(tru["truth_type"]).astype("i4")
    sel_type = np.isin(tru_type, [TRUTH_GAL, TRUTH_STAR])
    tru_ra = _native(tru["ra"])[sel_type]
    tru_dec = _native(tru["dec"])[sel_type]
    tru_id = _native(tru["id"]).astype("i8")[sel_type]
    tru_cid = _native(tru["cosmodc2_id"]).astype("i8")[sel_type]
    tru_chp = _native(tru["cosmodc2_hp"]).astype("i8")[sel_type]
    tru_type = tru_type[sel_type]
    tru_mag_r = _native(tru["mag_r"]).astype("f8")[sel_type]

    obj_ra = _native(obj["ra"])
    obj_dec = _native(obj["dec"])
    idx, sep = nearest_match(obj_ra, obj_dec, tru_ra, tru_dec)
    matched = sep < LSST_MATCH_RADIUS
    sn = _native(obj["psFlux_r"]).astype("f8") / _native(obj["psFluxErr_r"]).astype("f8")
    clean = _native(obj["clean"]).astype(bool)
    with np.errstate(invalid="ignore"):
        det_ok = matched & (sn > SNR_DEPTH) & clean

    ext = _native(obj["extendedness"]).astype("f8")
    cand = pd.DataFrame(
        {"truth_id": tru_id[idx[det_ok]], "sep": sep[det_ok], "extendedness": ext[det_ok]}
    )
    cand = cand.sort_values("sep").drop_duplicates("truth_id", keep="first")
    detected_id = set(cand["truth_id"].to_numpy())
    classified_id = set(cand.loc[cand["extendedness"] < EXT_CUT, "truth_id"].to_numpy())

    out = pd.DataFrame(
        {
            "id": tru_id,
            "ra": tru_ra,
            "dec": tru_dec,
            "truth_type": tru_type,
            "cosmodc2_id": tru_cid,
            "cosmodc2_hp": tru_chp,
            "mag_r": tru_mag_r,
            "tract": tract,
        }
    )
    out["lsst_detected"] = out["id"].isin(detected_id)
    out["lsst_star_class"] = out["id"].isin(classified_id)
    return out


def build_lsst_truth_catalog(refresh=False, n_tracts=0):
    if CACHE_LSST_TRUTH.exists() and not refresh and not n_tracts:
        print(f"loading cached LSST truth {CACHE_LSST_TRUTH.name}")
        return pd.read_parquet(CACHE_LSST_TRUTH)
    tracts = discover_overlap_tracts()
    if n_tracts:
        tracts = tracts[:n_tracts]
    print(f"{len(tracts)} LSST tracts overlap the Roman footprint: {tracts}")
    frames = []
    t0 = time.time()
    for tract in tracts:
        tt0 = time.time()
        df = process_lsst_tract(tract)
        frames.append(df)
        print(
            f"  tract {tract}: {len(df):,} truth (gal/star), "
            f"{df.lsst_detected.mean():.1%} detected ({time.time()-tt0:.0f}s)"
        )
    out = pd.concat(frames, ignore_index=True)
    before = len(out)
    out = out.drop_duplicates("id").reset_index(drop=True)
    print(f"LSST truth: {before:,} rows across tracts -> {len(out):,} unique ({time.time()-t0:.0f}s)")
    if not n_tracts:
        out.to_parquet(CACHE_LSST_TRUTH)
    return out


# --------------------------------------------------------------------------- #
# cross-survey truth match + compact-size join
# --------------------------------------------------------------------------- #
# NOTE (found during review): ``data/surveys/roman_dc2/cosmodc2_galaxy_size_true.parquet``
# (``COSMODC2_SIZES`` below, ~6.18M rows) is NOT a spatial/random subset of the
# cosmoDC2 truth catalog -- it was built by joining sizes only for galaxies that
# were already DETECTED/matched (the same ``needed_ids`` filtering
# ``load_cosmodc2_sizes`` does in ``scripts/lsst/build_lsst_dc2_galaxy_misclass.py``,
# whose ``needed_ids`` come from a det-centric matched sample; its row count,
# ~6.18M, closely tracks ``measure_joint_misclassification.py``'s det<->det
# matched-table size, 6.33M). Joining a truth-centric sample against it makes
# "has a size" almost synonymous with "was detected" -- exactly the selection
# effect this script exists to avoid. ``join_true_sizes`` below (unused,
# kept only so the failure mode is documented) is the naive join that produces
# that bias; ``join_true_sizes_full`` is the fix: it reads sizes for EVERY
# galaxy in the needed cosmoDC2 healpix (hp32) pixels directly from the raw
# per-pixel skims (``cosmoDC2_v1.1.4_image_summary_skim_hp32_<hp>.fits``,
# unfiltered by detection), so coverage no longer depends on being detected.
def join_true_sizes(cosmodc2_id):  # noqa: D401 -- kept for the record, NOT used (see note above)
    sz = pd.read_parquet(COSMODC2_SIZES, columns=["cosmodc2_id", "size_true"])
    lut = pd.Series(sz["size_true"].to_numpy(), index=sz["cosmodc2_id"].astype("i8").to_numpy())
    out = pd.Series(np.asarray(cosmodc2_id, "i8")).map(lut).to_numpy(dtype=float)
    print(f"  [DETECTION-BIASED, do not use for the compact-galaxy sample] cosmoDC2 size_true coverage: {np.isfinite(out).mean():.2%}")
    return out


CACHE_COSMODC2_SIZE_FULL = LSST_DIR / "_cache_cosmodc2_size_true_full_hp.parquet"


def load_cosmodc2_sizes_for_hp(hp_list, refresh=False):
    """Sizes for EVERY galaxy in the given cosmoDC2 hp32 pixels (no detection
    filter), read directly from the raw per-pixel skims. Cached across pixels
    already fetched by a previous run; only missing pixels are (re-)read."""
    hp_list = sorted(int(h) for h in hp_list)
    cached = pd.DataFrame(columns=["cosmodc2_hp", "galaxy_id", "size_true"])
    if CACHE_COSMODC2_SIZE_FULL.exists() and not refresh:
        cached = pd.read_parquet(CACHE_COSMODC2_SIZE_FULL)
    have = set(cached["cosmodc2_hp"].unique().tolist())
    missing = [h for h in hp_list if h not in have]
    if missing:
        print(f"  reading {len(missing)} cosmoDC2 hp32 skim(s) not yet cached: {missing}")
        frames = [cached] if len(cached) else []
        for hp in missing:
            f = LSST_RAW_DIR / f"cosmoDC2_v1.1.4_image_summary_skim_hp32_{hp}.fits"
            if not f.exists():
                print(f"    hp {hp}: skim file missing -- skipped")
                continue
            d = fitsio.read(str(f), columns=["galaxy_id", "size_true"])
            df = pd.DataFrame(
                {
                    "cosmodc2_hp": hp,
                    "galaxy_id": _native(d["galaxy_id"]).astype("i8"),
                    "size_true": _native(d["size_true"]).astype("f4"),
                }
            )
            frames.append(df)
            print(f"    hp {hp}: {len(df):,} galaxies (ALL, unfiltered by detection)")
        cached = pd.concat(frames, ignore_index=True)
        cached.to_parquet(CACHE_COSMODC2_SIZE_FULL)
    return cached[cached["cosmodc2_hp"].isin(hp_list)]


def join_true_sizes_full(cosmodc2_id, cosmodc2_hp, refresh=False):
    """cosmoDC2 ``size_true`` for every galaxy, keyed by ``galaxy_id ==
    cosmodc2_id`` -- NOT detection-biased (see the module-level note above).
    Memory note: avoids building a Python ``set``/dict of the (up to tens of
    millions of) ids; the lookup is a ``pandas.Series`` (hash-indexed under the
    hood) used with ``.map()``, matching the coordinator's guidance."""
    hp_arr = np.asarray(cosmodc2_hp, "i8")
    needed_hp = sorted(int(h) for h in np.unique(hp_arr) if h >= 0)
    print(f"cosmoDC2 size join (full, unbiased): {len(needed_hp)} unique hp32 pixel(s) needed: {needed_hp}")
    sz = load_cosmodc2_sizes_for_hp(needed_hp, refresh=refresh)
    lut = pd.Series(sz["size_true"].to_numpy(), index=sz["galaxy_id"].to_numpy())
    lut = lut[~lut.index.duplicated(keep="first")]
    out = pd.Series(np.asarray(cosmodc2_id, "i8")).map(lut).to_numpy(dtype=float)
    print(f"  cosmoDC2 size_true coverage (full hp32 skim join): {np.isfinite(out).mean():.2%}")
    return out


def cross_match(lsst_truth, roman_truth, radius_arcsec):
    idx, sep = nearest_match(
        lsst_truth["ra"].to_numpy(), lsst_truth["dec"].to_numpy(),
        roman_truth["ra"].to_numpy(), roman_truth["dec"].to_numpy(),
    )
    print("\ncross-match separation percentiles (LSST truth -> nearest Roman truth):")
    pcts = [10, 25, 50, 75, 90, 95, 99]
    vals = np.percentile(sep, pcts)
    for p, v in zip(pcts, vals):
        print(f"  p{p:>2d} = {v:.4f}\"")
    n_within = {r: int((sep < r).sum()) for r in (0.05, 0.1, 0.2, 0.3, 0.5, 1.0)}
    print(f"  N(sep < r): {n_within}  (of {len(sep):,} LSST truth objects)")

    keep = sep < radius_arcsec
    sample = lsst_truth.loc[keep].copy()
    sample["roman_sep_arcsec"] = sep[keep]
    ridx = idx[keep]
    sample["roman_gal_star"] = roman_truth["gal_star"].to_numpy()[ridx]
    sample["roman_detected"] = roman_truth["roman_detected"].to_numpy()[ridx]
    sample["roman_star_class"] = roman_truth["roman_star_class"].to_numpy()[ridx]
    sample["roman_detected_flagcut"] = roman_truth["roman_detected_flagcut"].to_numpy()[ridx]
    sample["roman_star_class_flagcut"] = roman_truth["roman_star_class_flagcut"].to_numpy()[ridx]
    sample["mag_H158_true"] = roman_truth["mag_H158"].to_numpy()[ridx]
    print(
        f"cross-match radius = {radius_arcsec}\" -> {len(sample):,} / {len(lsst_truth):,} "
        f"LSST truth objects matched to Roman truth ({len(sample)/len(lsst_truth):.1%})"
    )
    return sample


def check_star_position_collapse(roman_truth):
    """Sanity check (coordinator request). ``create_streamobs_files_hlwas.py``
    notes that ~26% of Roman truth STARS have a duplicate entry at the exact
    same ``(ra, dec)`` under a different ``ind``, carrying only the J129 flux
    (other bands, INCLUDING H158, are ``0.0`` -> NaN); the shipped stellar-
    efficiency product therefore collapses truth stars to unique POSITIONS and
    counts a star detected if ANY member entry matched. This script does NOT do
    that collapse, and its finite-Delta-m filter (see ``main()``) drops a row
    whose picked entry has NaN H158 -- i.e. it WOULD drop a star from the
    sample if the Roman<->LSST cross-match happens to land on the J129-only
    duplicate rather than the full-band one. This checks it directly on the
    FULL Roman truth catalog (independent of the LSST cross-match, so it also
    covers duplicates that never entered the final sample) -- NOTE: this
    raw-catalog check is NOT the same population as this script's actual star
    sample (one Roman row per cross-matched LSST truth star), and DOES show a
    non-negligible ~18% relative shift in aggregate p_roman on collapsing;
    see the module docstring's caveats for how that reconciles (or doesn't,
    within the time budget of this task) with the ~1% marginal-validation
    agreement of this script's actual star sample against
    ``roman_stellar_efficiency_cutf158.csv``."""
    stars = roman_truth[roman_truth["gal_star"] == 1]
    n_uncollapsed = len(stars)
    p_bare_uncollapsed = float(stars["roman_detected"].mean())
    p_flagcut_uncollapsed = float(stars["roman_detected_flagcut"].mean())
    coll = stars.groupby(["ra", "dec"], sort=False).agg(
        roman_detected=("roman_detected", "any"),
        roman_detected_flagcut=("roman_detected_flagcut", "any"),
    )
    n_collapsed = len(coll)
    p_bare_collapsed = float(coll["roman_detected"].mean())
    p_flagcut_collapsed = float(coll["roman_detected_flagcut"].mean())
    print("\n----- star duplicate-position collapse sanity check (coordinator request) -----")
    print(
        f"  {n_uncollapsed:,} truth-star ROWS -> {n_collapsed:,} unique POSITIONS "
        f"({(1 - n_collapsed / n_uncollapsed):.1%} are same-position duplicates, "
        f"cf. the ~26% figure in create_streamobs_files_hlwas.py)"
    )
    print(
        f"  p_roman (bare)        : uncollapsed={p_bare_uncollapsed:.4f}  "
        f"collapsed={p_bare_collapsed:.4f}  ratio={p_bare_collapsed / p_bare_uncollapsed:.4f}"
    )
    print(
        f"  p_roman (flags<{FLAG_CUT})   : uncollapsed={p_flagcut_uncollapsed:.4f}  "
        f"collapsed={p_flagcut_collapsed:.4f}  ratio={p_flagcut_collapsed / p_flagcut_uncollapsed:.4f}"
    )
    return dict(
        n_uncollapsed=n_uncollapsed,
        n_collapsed=n_collapsed,
        duplicate_position_frac=1 - n_collapsed / n_uncollapsed,
        p_roman_bare_uncollapsed=p_bare_uncollapsed,
        p_roman_bare_collapsed=p_bare_collapsed,
        p_roman_flagcut_uncollapsed=p_flagcut_uncollapsed,
        p_roman_flagcut_collapsed=p_flagcut_collapsed,
    )


def report_size_coverage(gal_all, size_true, nside=1024):
    """Diagnose whether the cosmoDC2 size-covered galaxy subset is a
    depth-representative sub-sample of the full cross-matched galaxy sample.

    ``size_true`` is the ONLY compactness measure available for a galaxy that
    NEITHER survey detects (there is no measured size to fall back on), so
    restricting the compact-galaxy sample to size-covered objects is
    unavoidable -- this only CHECKS whether that restriction is safe. If the
    detected fractions on the covered subset track the full cross-matched
    sample to a few percent, the size-covered footprint is not a biased
    (e.g. shallower/deeper) sub-region and the galaxy result generalizes; if
    not, the compact-galaxy conclusion should be scoped to the covered region.
    """
    covered = np.isfinite(size_true)
    n_tot, n_cov = len(gal_all), int(covered.sum())
    print("\n----- cosmoDC2 size-coverage diagnostic (compact-galaxy sample only) -----")
    print(f"  size_true coverage: {n_cov:,} / {n_tot:,} cross-matched true galaxies ({n_cov/n_tot:.2%})")
    ra, dec = gal_all["ra"].to_numpy(), gal_all["dec"].to_numpy()
    for label, m in [("covered", covered), ("all cross-matched", np.ones(n_tot, bool))]:
        print(f"  [{label:>18s}] RA {ra[m].min():.3f}-{ra[m].max():.3f}  Dec {dec[m].min():.3f}-{dec[m].max():.3f}")
    pixarea = hp.nside2pixarea(nside, degrees=True)
    pix_cov = np.unique(hp.ang2pix(nside, ra[covered], dec[covered], lonlat=True))
    pix_all = np.unique(hp.ang2pix(nside, ra, dec, lonlat=True))
    area_cov, area_all = len(pix_cov) * pixarea, len(pix_all) * pixarea
    print(
        f"  effective area (nside={nside}): covered={area_cov:.3f} sq deg ({len(pix_cov)} px)  "
        f"vs all cross-matched={area_all:.3f} sq deg ({len(pix_all)} px)  "
        f"-> covered/all = {area_cov/area_all:.2%}"
    )
    lsst_det = gal_all["lsst_detected"].to_numpy()
    roman_det = gal_all["roman_detected"].to_numpy()
    out = dict(
        n_total=n_tot,
        n_covered=n_cov,
        coverage_frac=n_cov / n_tot,
        area_covered_sqdeg=float(area_cov),
        area_all_sqdeg=float(area_all),
        lsst_detected_rate_covered=float(lsst_det[covered].mean()),
        lsst_detected_rate_all=float(lsst_det.mean()),
        roman_detected_rate_covered=float(roman_det[covered].mean()),
        roman_detected_rate_all=float(roman_det.mean()),
    )
    print(
        f"  lsst_detected  rate: covered={out['lsst_detected_rate_covered']:.4f}  "
        f"all={out['lsst_detected_rate_all']:.4f}  "
        f"ratio={out['lsst_detected_rate_covered']/out['lsst_detected_rate_all']:.3f}"
    )
    print(
        f"  roman_detected rate: covered={out['roman_detected_rate_covered']:.4f}  "
        f"all={out['roman_detected_rate_all']:.4f}  "
        f"ratio={out['roman_detected_rate_covered']/out['roman_detected_rate_all']:.3f}"
    )
    lsst_dev = abs(out["lsst_detected_rate_covered"] - out["lsst_detected_rate_all"])
    roman_dev = abs(out["roman_detected_rate_covered"] - out["roman_detected_rate_all"])
    out["lsst_rate_abs_deviation"] = float(lsst_dev)
    out["roman_rate_abs_deviation"] = float(roman_dev)
    passed = (lsst_dev < 0.03) and (roman_dev < 0.03)
    out["representativeness_check_passed"] = bool(passed)
    verdict = "PASS" if passed else ("FAIL (>10pp)" if max(lsst_dev, roman_dev) > 0.10 else "MARGINAL")
    print(
        f"  representativeness check: LSST |Delta|={lsst_dev:.4f}  Roman |Delta|={roman_dev:.4f}  "
        f"-> {verdict} (pass threshold: both < 0.03 absolute)"
    )
    return out


# --------------------------------------------------------------------------- #
# the measurement (mirrors measure_joint_misclassification.measure())
# --------------------------------------------------------------------------- #
def measure(sub, label, quantity, flag_l_col, flag_r_col, size_desc, axis="roman"):
    dm_l = sub["dm_lsst"].to_numpy()
    dm_r = sub["dm_roman"].to_numpy()
    dm_axis = dm_r if axis == "roman" else dm_l
    sl = sub[flag_l_col].to_numpy()
    sr = sub[flag_r_col].to_numpy()
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
                quantity=quantity,
                selection=name,
                dm_axis=axis,
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
    if axis == "roman":
        _row("ALL", np.ones(len(sub), bool))
        _row("F158_lt_25.5", sub["mag_H158_true"].to_numpy() < JOINT_F158_CUT)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# validation: reproduce the SHIPPED stellar-efficiency curves (stars only)
# --------------------------------------------------------------------------- #
PROD_MAG_BINS = np.arange(15.0, 29.01, 0.25)
PROD_MAG_MID = 0.5 * (PROD_MAG_BINS[1:] + PROD_MAG_BINS[:-1])


def _binned_rate(mag, flag, n_min=20):
    n = np.histogram(mag, PROD_MAG_BINS)[0]
    k = np.histogram(mag[flag], PROD_MAG_BINS)[0]
    with np.errstate(invalid="ignore"):
        p = np.where(n >= n_min, k / np.maximum(n, 1), np.nan)
    return n, p


def _compare(name, csv, mag, flag_det, flag_detcls, magcol, det_col, detcls_col):
    ship = pd.read_csv(csv, comment="#", header=None)
    ncols = ship.shape[1]
    names = [magcol, "delta_mag", "detection_eff", "classification_eff", "classification_detection_eff"][:ncols]
    ship.columns = names
    for quantity, flag, ship_col in [
        ("detection", flag_det, "detection_eff"),
        ("detection_classification", flag_detcls, "classification_detection_eff"),
    ]:
        if ship_col not in ship.columns:
            continue
        n, p = _binned_rate(mag, flag)
        mine = pd.DataFrame({magcol: PROD_MAG_MID, "n": n, "p_here": p}).dropna()
        cmp_ = ship.merge(mine, on=magcol, how="inner")
        cmp_ = cmp_[(cmp_[ship_col] > 1e-3) & (cmp_.n >= 200)]
        if not len(cmp_):
            print(f"  {name} [{quantity}]: no overlapping well-populated bins")
            continue
        ratio = cmp_.p_here / cmp_[ship_col]
        print(
            f"  {name} [{quantity}]: {len(cmp_)} bins, "
            f"p_here/p_shipped median={np.median(ratio):.3f} "
            f"[p16={np.percentile(ratio,16):.3f}, p84={np.percentile(ratio,84):.3f}]"
        )
        show = cmp_.iloc[np.linspace(0, len(cmp_) - 1, min(6, len(cmp_))).astype(int)]
        for _, r in show.iterrows():
            print(
                f"      {magcol}={r[magcol]:.2f}  shipped={r[ship_col]:.5f}  "
                f"here={r.p_here:.5f}  ratio={r.p_here / r[ship_col]:.3f}  (N={int(r.n):,})"
            )


def validate_marginals(star_sample):
    """Reproduce the shipped stellar-efficiency curves from this truth-centric
    sample. Two expected sources of mismatch (documented in the module
    docstring): (1) this script's Roman-detected flag omits the shipped
    product's ``flags==0`` cut, so the Roman detection efficiency here should
    run a bit ABOVE the shipped curve; (2) this sample is restricted to the
    Roman/LSST cross-match overlap (not the full LSST or Roman footprint), so
    per-bin counts are smaller and noisier than the shipped curves'."""
    print("\n----- marginal validation vs the shipped stellar-efficiency curves (stars) -----")
    _compare(
        "LSST r stellar efficiency",
        LSST_DIR / "lsst_dc2_stellar_efficiency_cutr.csv",
        star_sample["mag_r"].to_numpy(),
        star_sample["lsst_detected"].to_numpy(),
        (star_sample["lsst_detected"] & star_sample["lsst_star_class"]).to_numpy(),
        "mag_r",
        "detection_eff",
        "classification_detection_eff",
    )
    _compare(
        "Roman F158 stellar efficiency",
        ROMAN_DIR / "roman_stellar_efficiency_cutf158.csv",
        star_sample["mag_H158_true"].to_numpy(),
        star_sample["roman_detected"].to_numpy(),
        (star_sample["roman_detected"] & star_sample["roman_star_class"]).to_numpy(),
        "mag_f158",
        "detection_eff",
        "classification_detection_eff",
    )
    # bright-star sanity check
    for name, mag, flag in [
        ("LSST", star_sample["mag_r"].to_numpy(), star_sample["lsst_detected"].to_numpy()),
        ("Roman", star_sample["mag_H158_true"].to_numpy(), star_sample["roman_detected"].to_numpy()),
    ]:
        bright = mag < np.nanpercentile(mag, 5)
        faint = mag > np.nanpercentile(mag, 97)
        print(
            f"  {name} bright-star (< p5 mag) detected frac = {flag[bright].mean():.3f}  "
            f"(N={int(bright.sum())}); faint-star (> p97 mag) detected frac = "
            f"{flag[faint].mean():.3f} (N={int(faint.sum())})"
        )


# --------------------------------------------------------------------------- #
# figure
# --------------------------------------------------------------------------- #
def make_figure(tab, meta):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.4), sharex=True)
    # AUTHORITATIVE quantity (Roman flags<FLAG_CUT applied, matches the shipped
    # products -- see the module docstring); "all_galaxies" (no compact cut) is
    # what the paper's joint background actually uses.
    for j, (stype, title) in enumerate(
        [("all_galaxies", "true galaxies (no compact cut)"), ("stars", "true stars")]
    ):
        t = tab[
            (tab.source_type == stype)
            & (tab.selection == "dm_bin")
            & (tab.dm_axis == "roman")
            & (tab.quantity == "detection_classification_flagcut")
        ].sort_values("dm_bin")
        ax, axr = axes[0, j], axes[1, j]
        if len(t):
            x = t.dm_bin.values
            ax.errorbar(
                x, t.p_lsst, yerr=[t.p_lsst - t.p_lsst_lo, t.p_lsst_hi - t.p_lsst],
                fmt="o-", ms=3.5, color="C0", label=r"$p_{\rm LSST}$ (det$\times$cls)",
            )
            ax.errorbar(
                x, t.p_roman, yerr=[t.p_roman - t.p_roman_lo, t.p_roman_hi - t.p_roman],
                fmt="s-", ms=3.5, color="C3", label=r"$p_{\rm Roman}$ (det$\times$cls)",
            )
            ax.errorbar(
                x, t.p_joint, yerr=[t.p_joint - t.p_joint_lo, t.p_joint_hi - t.p_joint],
                fmt="D", ms=5, color="k", label=r"$p_{\rm joint}$ measured",
            )
            ax.plot(x, t.pred_indep, "--", color="tab:green", lw=2, label="independent null")
            ax.plot(x, t.pred_comono, ":", color="tab:purple", lw=2.4, label=r"comonotonic $\min(p_1,p_2)$")
            axr.axhline(1.0, color="0.6", lw=1)
            axr.plot(x, t.r_indep, "--o", color="tab:green", ms=3.5, label="vs independent")
            axr.plot(x, t.r_comono, ":s", color="tab:purple", ms=3.5, label="vs comonotonic")
        ax.set_title(f"{title}   (N={int(t.n.sum()):,} in bins)", fontsize=11)
        ax.set_ylabel("P(detected AND classified point-like)")
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
        "LSST$\\times$Roman DETECTION$\\times$classification coupling in DC2 (truth-centric, "
        "Roman flags$<$FLAG\\_CUT applied)\n"
        f"maglim$_{{F158}}$={meta['ref_maglim_f158']:.2f}, maglim$_r$={meta['ref_maglim_r']:.2f}   |   "
        f"cross-match radius={meta['cross_match_radius_arcsec']}\"",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT_PNG, dpi=140)
    print(f"wrote {OUT_PNG.relative_to(REPO)}")


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refresh", action="store_true", help="rebuild the Roman/LSST truth caches")
    ap.add_argument("--tracts", type=int, default=0, help="limit to first N LSST tracts (smoke test)")
    ap.add_argument("--roman-tiles", type=int, default=0, help="limit to first N Roman tiles (smoke test)")
    ap.add_argument("--cross-radius", type=float, default=0.2, help="Roman<->LSST truth match radius (arcsec)")
    args = ap.parse_args()
    ART_DIR.mkdir(exist_ok=True)

    ml_f158, npix_f158 = ref_maglim(ROMAN_MAGLIM)
    ml_r, npix_r = ref_maglim(LSST_MAGLIM)
    print(f"ref maglim F158 = {ml_f158:.4f} ({npix_f158} px);  r = {ml_r:.4f} ({npix_r} px)")

    roman_truth = build_roman_truth_catalog(refresh=args.refresh, limit=args.roman_tiles)
    star_dup_check = check_star_position_collapse(roman_truth)
    lsst_truth = build_lsst_truth_catalog(refresh=args.refresh, n_tracts=args.tracts)

    sample = cross_match(lsst_truth, roman_truth, args.cross_radius)

    sample["dm_roman"] = sample["mag_H158_true"].to_numpy() - ml_f158
    sample["dm_lsst"] = sample["mag_r"].to_numpy() - ml_r

    # a small number of Roman truth entries are duplicate stamps carrying only
    # ONE band's flux (see create_streamobs_files_hlwas.py's "gotcha" note on
    # ~26% of truth stars having a same-position duplicate with other bands ==
    # 0.0 -> NaN); drop non-finite delta_m rows before binning/cell_predictions,
    # which otherwise chokes on NaN -> int64 cast (cell_predictions has no NaN
    # guard, matching measure_joint_misclassification.py's own base-quality mask).
    # Done BEFORE any other boolean-mask derived from `sample` so all downstream
    # masks stay aligned to the (possibly shorter) filtered index.
    finite_dm = np.isfinite(sample["dm_roman"].to_numpy()) & np.isfinite(sample["dm_lsst"].to_numpy())
    n_dropped = int((~finite_dm).sum())
    if n_dropped:
        print(f"dropping {n_dropped:,} / {len(sample):,} rows with non-finite delta_m (NaN true mag)")
    sample = sample.loc[finite_dm].reset_index(drop=True)

    # ---- truth labels: require the two truth catalogs to AGREE ------------------
    is_gal = (sample["roman_gal_star"].to_numpy() == 0) & (sample["truth_type"].to_numpy() == TRUTH_GAL)
    is_star = (sample["roman_gal_star"].to_numpy() == 1) & (sample["truth_type"].to_numpy() == TRUTH_STAR)
    n_disagree = int((~is_gal & ~is_star).sum())
    print(
        f"truth-agreeing galaxies: {int(is_gal.sum()):,}   stars: {int(is_star.sum()):,}   "
        f"disagreements dropped: {n_disagree:,}"
    )

    # ---- compact-galaxy size: FULL cosmoDC2 hp32 skims, NOT the detection-
    # biased cosmodc2_galaxy_size_true.parquet (see the note above
    # join_true_sizes_full -- that parquet's ~6.18M rows track the det<->det
    # matched-table size almost exactly, so joining against it makes "has a
    # size" nearly synonymous with "was detected", defeating the point of a
    # truth-centric sample). -----------------------------------------------
    print("joining cosmoDC2 true sizes (full hp32 skims, unbiased by detection) ...")
    gal_all = sample.loc[is_gal].reset_index(drop=True)
    size_true = join_true_sizes_full(
        gal_all["cosmodc2_id"].to_numpy(), gal_all["cosmodc2_hp"].to_numpy(), refresh=args.refresh
    )
    compact = np.isfinite(size_true) & (size_true < GAL_SIZE_MAX)
    size_desc = f'cosmoDC2 size_true < {GAL_SIZE_MAX}" (joined on cosmodc2_id via the full, unfiltered hp32 skims)'

    # coordinator-requested diagnostic (now on the UNBIASED size join): is the
    # size-covered galaxy subset depth-representative of the full cross-matched
    # galaxy sample? PASS criterion: LSST/Roman detected rates agree to a few %.
    coverage_report = report_size_coverage(gal_all, size_true)

    keep_cols = [
        "dm_roman", "dm_lsst", "mag_r", "mag_H158_true",
        "lsst_detected", "lsst_star_class",
        "roman_detected", "roman_star_class",
        "roman_detected_flagcut", "roman_star_class_flagcut",
        "id",
    ]
    gal = gal_all.loc[compact, keep_cols].drop_duplicates("id").copy()
    all_gal = gal_all.loc[:, keep_cols].drop_duplicates("id").copy()
    star = sample.loc[is_star, keep_cols].drop_duplicates("id").copy()
    for sub in (gal, all_gal, star):
        sub["lsst_detcls"] = sub["lsst_detected"] & sub["lsst_star_class"]
        sub["roman_detcls"] = sub["roman_detected"] & sub["roman_star_class"]
        sub["roman_detcls_flagcut"] = sub["roman_detected_flagcut"] & sub["roman_star_class_flagcut"]
    print(
        f"compact true galaxies (post cross-match): {len(gal):,}   "
        f"ALL true galaxies (no compact cut, coordinator request -- this is what "
        f"fig4_metricC2_ingredients.py's ``gal_idx = np.where(~is_star)[0]`` actually "
        f"injects): {len(all_gal):,}   true stars: {len(star):,}"
    )

    # ---- validation ---------------------------------------------------------------
    validate_marginals(star)

    # ---- measurement ----------------------------------------------------------------
    # Four quantities per source type: bare detection / detection+classification
    # (the handoff's definition), and their Roman-``flags<FLAG_CUT`` counterparts
    # (coordinator's request to check whether the shipped-product flags cut
    # moves r_indep/r_comono/rho). LSST's flag definition is unchanged in all four
    # -- the flags cut is a Roman-only SExtractor quality flag.
    t0 = time.time()
    all_gal_size_desc = "N/A -- ALL true galaxies, no compact-size cut applied (coordinator request: the paper's joint background is NOT restricted to compact galaxies)"
    tab = pd.concat(
        [
            measure(sub, lab, quantity, lcol, rcol, sdesc, axis=ax)
            for lab, sub, sdesc in [
                ("galaxies", gal, size_desc),
                ("all_galaxies", all_gal, all_gal_size_desc),
                ("stars", star, size_desc),
            ]
            for quantity, lcol, rcol in [
                ("detection", "lsst_detected", "roman_detected"),
                ("detection_classification", "lsst_detcls", "roman_detcls"),
                ("detection_flagcut", "lsst_detected", "roman_detected_flagcut"),
                ("detection_classification_flagcut", "lsst_detcls", "roman_detcls_flagcut"),
            ]
            for ax in ("roman", "lsst")
        ],
        ignore_index=True,
    )
    print(f"binned measurement: {time.time()-t0:.0f}s")

    tab.to_csv(OUT_CSV, index=False, float_format="%.6g")
    print(f"wrote {OUT_CSV.relative_to(REPO)} ({len(tab)} rows)")

    meta = dict(
        script=str(Path(__file__).relative_to(REPO)),
        n_lsst_truth_overlap_candidates=int(len(lsst_truth)),
        n_cross_matched=int(len(sample)),
        cross_match_radius_arcsec=args.cross_radius,
        n_compact_galaxies=int(len(gal)),
        n_all_galaxies=int(len(all_gal)),
        n_stars=int(len(star)),
        n_truth_disagreements_dropped=n_disagree,
        size_source=size_desc,
        ref_maglim_f158=ml_f158,
        ref_maglim_r=ml_r,
        roman_classifier="F158 size envelope (roman_star_classifier.build_env_classifier, ENV_PURITY=0.875)",
        lsst_classifier=f"extendedness < {EXT_CUT}",
        roman_detected_definition="matched==True in roman_dc2_det_truth.parquet (bare, NO flags==0 cut) for the 'detection'/'detection_classification' quantities [DIAGNOSTIC ONLY]; matched==True AND flags<FLAG_CUT (flags==0) for the '*_flagcut' quantities [AUTHORITATIVE -- create_streamobs_files_hlwas.py:726 and build_roman_galaxy_misclass.py both apply flags==0 to define Roman detection/classification in the shipped selection-function products, so the '*_flagcut' quantities are the ones that match the code being validated; the bare-matched quantities are kept only to show how much the cut matters, which for galaxies is NOT negligible -- see 'flagcut_moves_r_indep_materially_for_galaxies']",
        roman_flag_cut=FLAG_CUT,
        authoritative_quantities=["detection_flagcut", "detection_classification_flagcut"],
        flagcut_moves_r_indep_materially_for_galaxies="YES for galaxies (e.g. detection_classification, F158<25.5: r_indep 0.9995 [bare] -> 1.52 [flagcut], because flags<1 drops Roman p_roman by a factor ~2.4 while p_joint drops only ~1.4x); NO for stars (e.g. detection_classification, ALL: r_indep 1.028 -> 1.042, a 1.4% shift) -- see the printed summary / CSV for the full breakdown.",
        lsst_detected_definition=f"nearest LSST object within {LSST_MATCH_RADIUS}\" with psFlux_r/psFluxErr_r>{SNR_DEPTH} and clean (unchanged across all four quantities -- flags is a Roman-only SExtractor QA flag)",
        conditioning="truth-centric: ALL truth objects in the Roman x LSST cross-match overlap, whether or not either survey detected them",
        cell_size_mag=CELL,
        joint_f158_cut=JOINT_F158_CUT,
        overlap_tracts=discover_overlap_tracts() if not args.tracts else None,
        compact_galaxy_size_coverage=coverage_report,
        star_duplicate_position_check=star_dup_check,
    )
    OUT_META.write_text(json.dumps(meta, indent=2, default=str))
    make_figure(tab, meta)

    # ---- printed summary --------------------------------------------------------
    pd.set_option("display.width", 220, "display.max_columns", 50)
    for stype in ("galaxies", "all_galaxies", "stars"):
        for quantity in ("detection", "detection_classification", "detection_flagcut", "detection_classification_flagcut"):
            t = tab[(tab.source_type == stype) & (tab.quantity == quantity)]
            print(f"\n===== {stype} :: {quantity} =====")
            cols = [
                "selection", "dm_axis", "dm_bin", "n",
                "p_lsst", "p_roman", "p_joint",
                "pred_indep", "pred_comono", "r_indep", "r_comono", "rho_copula",
            ]
            print(t[cols].to_string(index=False, float_format=lambda v: f"{v:.4g}"))

    # ---- sample precision (coordinator request): N + Wilson CI on p_joint -----------
    print("\n----- sample precision: p_joint N and Wilson 68% CI (ALL, F158<25.5) -----")
    prec = tab[(tab.dm_axis == "roman") & (tab.selection.isin(["ALL", "F158_lt_25.5"]))]
    for stype in ("galaxies", "all_galaxies", "stars"):
        for _, r in prec[prec.source_type == stype].iterrows():
            print(
                f"  [{stype:>13s} | {r.quantity:>32s} | {r.selection:>12s}] N={int(r.n):>10,}  "
                f"p_joint={r.p_joint:.5f}  [{r.p_joint_lo:.5f}, {r.p_joint_hi:.5f}]  "
                f"r_indep={r.r_indep:.3f}  r_comono={r.r_comono:.3f}  rho={r.rho_copula:.3f}"
            )


if __name__ == "__main__":
    main()
