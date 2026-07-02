#!/usr/bin/env python
"""Build a positionally-matched Roman<->LSST merged table and a galaxy
MISCLASSIFICATION efficiency curve for the Roman DC2 selection function (D2).

Adapted from ``lsst_dc2_scratch/create_contam_model.py`` (do NOT edit that scratch
file). That script positionally matched (LSST object, Roman det, Roman truth)
triplets and built a ``contam_cat_roman`` sample (Roman ``class_star > 0.5`` &
truth ``gal_star < 0.5`` = true galaxies Roman calls stars). Here we:

  1. Use the Roman det<->truth match that is already baked into
     ``roman_dc2_det_truth.parquet`` (``matched``, ``truth_*`` columns, 1" match).
  2. Positionally match Roman det <-> LSST object within 1" (the existing recipe).
  3. SAVE the merged table = LSST_true + roman_true + lsst_obs + roman_obs
     (one row per matched object), namespaced by origin, to
     ``roman_lsst_matched.parquet``.
  4. Compute the galaxy-MISCLASSIFICATION efficiency curve: among true galaxies
     (Roman truth ``gal_star == 0``) restricted to COMPACT galaxies
     (``size < GAL_SIZE_MAX = 0.3"``), the fraction classified as stars by the
     Roman size-envelope classifier, vs magnitude. Written to
     ``roman_galaxy_misclass_cutf158.csv``.

RE-RUN ARCHITECTURE (critical): the galaxy SIZE used for the compact cut comes
from a SINGLE swappable input -- the top-level ``SIZE_SOURCE`` switch:

  * ``"measured_roman"`` (INTERIM, used now): measured Roman F158 size
    ``size_sb = sqrt(lambda1) * 3600"`` from the per-band windowed second moments
    (``x2/y2/xywin_world_H158``) -- the SAME size the existing classifier uses.
    NOTE: measured size conflates with the PSF (~0.27" FWHM in H158) and is an
    interim proxy.
  * ``"cosmodc2_true"`` (FUTURE, not yet available): a cosmoDC2 ``size_true`` file
    joined by LSST ``cosmodc2_id``. See ``load_true_galaxy_sizes`` below for the
    documented hook/stub.

Switching ``SIZE_SOURCE`` + re-running regenerates everything with the new size
cut, NO other code changes.

Run (streamobs conda env):
    python scripts/roman/build_roman_galaxy_misclass.py                   # measured_roman (default)
    python scripts/roman/build_roman_galaxy_misclass.py --size-source cosmodc2_true
    python scripts/roman/build_roman_galaxy_misclass.py --force           # re-do the LSST matching cache
"""

import argparse
import sys
from glob import glob
from pathlib import Path

import fitsio
import numpy as np
import pandas as pd
import pyarrow.dataset as pads
import healpy as hp
from astropy.coordinates import SkyCoord
import astropy.units as u

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
REPO = Path(__file__).resolve().parents[2]
ROMAN_DIR = REPO / "data" / "surveys" / "roman_dc2"
ROMAN_PARQUET = ROMAN_DIR / "roman_dc2_det_truth.parquet"
LSST_DIR = Path("/astro/store/shire/stream_team/stream_finding/data/lsst_dc2")
MAGLIM_MAP = ROMAN_DIR / "roman_dc2_maglim_f158_nside1024.fits.gz"

MERGED_OUT = ROMAN_DIR / "roman_lsst_matched.parquet"
MISCLASS_CSV = ROMAN_DIR / "roman_galaxy_misclass_cutf158.csv"
CACHE_OUT = ROMAN_DIR / "roman_lsst_matched.parquet"   # the merged table doubles as the cache

BAND = "H158"                       # Roman F158
MATCH_RADIUS_ROMAN_LSST = 1.0       # arcsec (existing recipe in create_contam_model)
GAL_SIZE_MAX = 0.3                  # arcsec: "compact" galaxy cut
FLAG_CUT = 1                        # flags == 0 (paper-exact, matches HLWAS products)

# magnitude binning for the efficiency curve (same convention as the HLWAS script)
MAG_BINS = np.arange(15.0, 29.01, 0.25)
MAG_MID = 0.5 * (MAG_BINS[1:] + MAG_BINS[:-1])

# --------------------------------------------------------------------------- #
# size source switch (see module docstring)  --  the ONE knob to flip
# --------------------------------------------------------------------------- #
SIZE_SOURCE = "measured_roman"      # "measured_roman" (interim) | "cosmodc2_true" (future)

# FUTURE hook: path to a cosmoDC2 true-size file, joined by LSST cosmodc2_id.
# Not yet available -- populate when the file lands, then run with
# --size-source cosmodc2_true.
COSMODC2_SIZE_FILE = ROMAN_DIR / "cosmodc2_galaxy_size_true.parquet"
COSMODC2_SIZE_ID_COL = "cosmodc2_id"      # join key in that file (matches LSST truth)
COSMODC2_SIZE_VALUE_COL = "size_true"     # arcsec; semi-major axis half-light or similar


# Roman columns we actually need from the 8.2 GB parquet (read ONLY these).
ROMAN_COLS = [
    "tile",
    "alphawin_j2000", "deltawin_j2000",
    "mag_auto_H158", "magerr_auto_H158", "flux_auto_H158", "fluxerr_auto_H158",
    "class_star", "awin_world", "flags", "det_sn",
    "x2win_world_H158", "y2win_world_H158", "xywin_world_H158",
    "matched", "match_sep_arcsec",
    "truth_ind", "truth_ra", "truth_dec", "truth_gal_star",
    "truth_mag_Y106", "truth_mag_J129", "truth_mag_H158", "truth_mag_F184",
    "truth_bb_mag",
]

# LSST object (measured) columns.
LSST_OBJ_COLS = [
    "objectId", "ra", "dec",
    "mag_g", "mag_r", "mag_i", "mag_z", "mag_y",
    "mag_g_cModel", "mag_r_cModel", "mag_i_cModel", "mag_z_cModel", "mag_y_cModel",
    "extendedness", "blendedness", "clean",
]
# LSST truth columns.
LSST_TRUTH_COLS = [
    "id", "ra", "dec", "truth_type", "cosmodc2_id",
    "mag_u", "mag_g", "mag_r", "mag_i", "mag_z", "mag_y", "redshift",
]


# --------------------------------------------------------------------------- #
# Roman F158 size envelope classifier: imported from the shared module
# scripts/roman/roman_star_classifier.py (the SINGLE source of truth, also imported
# by create_streamobs_files_hlwas.py so the two cannot drift).
# --------------------------------------------------------------------------- #
sys.path.insert(0, str(Path(__file__).resolve().parent))
from roman_star_classifier import size_sb, build_env_classifier  # noqa: E402


# --------------------------------------------------------------------------- #
# Galaxy size source (the swappable input)
# --------------------------------------------------------------------------- #
def load_true_galaxy_sizes():
    """FUTURE hook: load cosmoDC2 true galaxy sizes, joined by LSST cosmodc2_id.

    Returns a DataFrame with columns [cosmodc2_id, size_true(arcsec)] or None if
    the file is not yet available. When SIZE_SOURCE == "cosmodc2_true", the
    merged table's ``lsst_true_cosmodc2_id`` is joined to this to get the true
    size used for the compact-galaxy cut. This is the ONLY change needed to
    re-run with true sizes -- no other code edits.
    """
    if not COSMODC2_SIZE_FILE.exists():
        return None
    df = pd.read_parquet(COSMODC2_SIZE_FILE,
                         columns=[COSMODC2_SIZE_ID_COL, COSMODC2_SIZE_VALUE_COL])
    return df.rename(columns={COSMODC2_SIZE_ID_COL: "cosmodc2_id",
                              COSMODC2_SIZE_VALUE_COL: "size_true"})


def resolve_galaxy_size(merged, size_source):
    """Return the galaxy size (arcsec) used for the compact cut, per SIZE_SOURCE.

    measured_roman : measured Roman F158 size_sb (interim; conflates with PSF).
    cosmodc2_true  : cosmoDC2 true size joined by lsst_true_cosmodc2_id.
    """
    if size_source == "measured_roman":
        sz = size_sb(merged.rename(columns={
            "roman_obs_x2win_world_H158": "x2win_world_H158",
            "roman_obs_y2win_world_H158": "y2win_world_H158",
            "roman_obs_xywin_world_H158": "xywin_world_H158",
        }))
        return np.asarray(sz), "measured Roman F158 size_sb = sqrt(lambda1)*3600\" (INTERIM proxy; PSF~0.27\")"

    if size_source == "cosmodc2_true":
        true_sizes = load_true_galaxy_sizes()
        if true_sizes is None:
            raise FileNotFoundError(
                f"SIZE_SOURCE='cosmodc2_true' but {COSMODC2_SIZE_FILE} is not available. "
                "Provide the cosmoDC2 true-size file (see load_true_galaxy_sizes) "
                "or run with --size-source measured_roman."
            )
        key = merged["lsst_true_cosmodc2_id"].values
        lut = dict(zip(true_sizes["cosmodc2_id"].values, true_sizes["size_true"].values))
        sz = np.array([lut.get(k, np.nan) for k in key], dtype=float)
        return sz, f"cosmoDC2 true size joined by cosmodc2_id from {COSMODC2_SIZE_FILE.name}"

    raise ValueError(f"unknown SIZE_SOURCE={size_source!r}")


# --------------------------------------------------------------------------- #
# I/O + matching
# --------------------------------------------------------------------------- #
def read_roman_bbox(dset, ra_lo, ra_hi, dec_lo, dec_hi):
    """Read the Roman det->truth rows inside an RA/Dec bbox (memory-aware).

    Reads ONLY ROMAN_COLS and only the rows whose detection lands in the bbox
    (a small margin is added by the caller). Returns a DataFrame.
    """
    f = pads.field
    filt = ((f("alphawin_j2000") >= ra_lo) & (f("alphawin_j2000") <= ra_hi)
            & (f("deltawin_j2000") >= dec_lo) & (f("deltawin_j2000") <= dec_hi))
    tbl = dset.to_table(columns=ROMAN_COLS, filter=filt)
    return tbl.to_pandas()


def load_lsst_tract(obj_file, truth_file):
    """Load one LSST tract: object (measured) + truth, matched truth->object 1".

    Returns (lsst_obs_df, lsst_true_df) row-aligned: lsst_true is the nearest
    truth within 1" of each object (NaN/empty where none).
    """
    def _native(a):
        # FITS arrays are big-endian; convert to native byte order for pandas.
        if a.dtype.byteorder == ">":
            return a.byteswap().view(a.dtype.newbyteorder())
        return a
    obj = fitsio.read(obj_file, columns=LSST_OBJ_COLS)
    tru = fitsio.read(truth_file, columns=LSST_TRUTH_COLS)
    obj_df = pd.DataFrame({c: _native(obj[c]) for c in LSST_OBJ_COLS})
    tru_df = pd.DataFrame({c: _native(tru[c]) for c in LSST_TRUTH_COLS})
    return obj_df, tru_df


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--size-source", default=SIZE_SOURCE,
                        choices=["measured_roman", "cosmodc2_true"],
                        help="galaxy size source for the compact-galaxy cut")
    parser.add_argument("--force", action="store_true",
                        help="recompute the LSST<->Roman match even if the merged table exists")
    args = parser.parse_args()
    size_source = args.size_source

    # REF_MAGLIM = median of the DC2 F158 maglim map (covered pixels only).
    mlm = hp.read_map(str(MAGLIM_MAP))
    covered = (mlm != hp.UNSEEN) & np.isfinite(mlm) & (mlm > 0)
    REF_MAGLIM = float(np.median(mlm[covered]))
    print(f"REF_MAGLIM (DC2 F158 maglim-map median) = {REF_MAGLIM:.4f}  "
          f"({int(covered.sum())} covered pixels)")
    print(f"SIZE_SOURCE = {size_source}")

    # ----------------------------------------------------------------------- #
    # Stage A: per-tract positional match (Roman det<->LSST object 1"), build merged
    # ----------------------------------------------------------------------- #
    if MERGED_OUT.exists() and not args.force:
        print(f"loading cached merged table {MERGED_OUT}")
        merged = pd.read_parquet(MERGED_OUT)
        n_roman_seen = n_roman_truth = n_lsst_obj = None
    else:
        dset = pads.dataset(str(ROMAN_PARQUET), format="parquet")
        obj_files = sorted(glob(str(LSST_DIR / "dc2_object_run2.2i_dr6_skim_tract_*.fits")))
        truth_files = {Path(f).name.split("_")[-1].split(".")[0]: f
                       for f in glob(str(LSST_DIR / "dc2_run2.2i_truth_merged_summary_skim_tract_*.fits"))}

        margin = 30.0 / 3600.0   # 30" tract margin to catch edge matches
        merged_parts = []
        n_roman_seen = n_roman_truth = n_lsst_obj = 0

        for obj_file in obj_files:
            tract = Path(obj_file).name.split("_")[-1].split(".")[0]
            tf = truth_files.get(tract)
            if tf is None:
                print(f"  tract {tract}: no truth file, skipping")
                continue

            lsst_obs, lsst_true = load_lsst_tract(obj_file, tf)
            n_lsst_obj += len(lsst_obs)

            ra_lo, ra_hi = lsst_obs["ra"].min() - margin, lsst_obs["ra"].max() + margin
            dec_lo, dec_hi = lsst_obs["dec"].min() - margin, lsst_obs["dec"].max() + margin
            roman = read_roman_bbox(dset, ra_lo, ra_hi, dec_lo, dec_hi)
            n_roman_seen += len(roman)
            if len(roman) == 0 or len(lsst_obs) == 0:
                print(f"  tract {tract}: roman={len(roman)} lsst_obj={len(lsst_obs)} -> 0 matches")
                continue

            # (a) Roman det<->Roman truth is ALREADY baked in (matched/truth_*).
            #     Keep only Roman detections that matched a truth object (1") so the
            #     merged row carries roman_true_*.
            roman = roman[roman["matched"].astype(bool)].reset_index(drop=True)
            n_roman_truth += len(roman)
            if len(roman) == 0:
                continue

            # (b) Roman det <-> LSST object, 1" (existing recipe): for each LSST
            #     object find nearest Roman detection.
            c_roman = SkyCoord(ra=roman["alphawin_j2000"].values * u.deg,
                               dec=roman["deltawin_j2000"].values * u.deg)
            c_obj = SkyCoord(ra=lsst_obs["ra"].values * u.deg,
                             dec=lsst_obs["dec"].values * u.deg)
            idx, d2d, _ = c_obj.match_to_catalog_sky(c_roman)
            sel = d2d.arcsec < MATCH_RADIUS_ROMAN_LSST

            lsst_obs_m = lsst_obs[sel].reset_index(drop=True)
            roman_m = roman.iloc[idx[sel]].reset_index(drop=True)

            # (c) LSST object <-> LSST truth, 1": attach the true type/cosmodc2_id.
            c_true = SkyCoord(ra=lsst_true["ra"].values * u.deg,
                              dec=lsst_true["dec"].values * u.deg)
            c_objm = SkyCoord(ra=lsst_obs_m["ra"].values * u.deg,
                              dec=lsst_obs_m["dec"].values * u.deg)
            tidx, td2d, _ = c_objm.match_to_catalog_sky(c_true)
            true_ok = td2d.arcsec < MATCH_RADIUS_ROMAN_LSST
            lsst_true_m = lsst_true.iloc[tidx].reset_index(drop=True)
            # null out LSST-truth fields where no truth within 1"
            for col in lsst_true_m.columns:
                if lsst_true_m[col].dtype.kind in "fc":
                    lsst_true_m.loc[~true_ok, col] = np.nan

            # assemble one namespaced row per matched object
            part = pd.DataFrame()
            for c in lsst_true_m.columns:
                part[f"lsst_true_{c}"] = lsst_true_m[c].values
            part["lsst_true_match_sep_arcsec"] = td2d.arcsec
            part["lsst_true_matched"] = true_ok
            for c in lsst_obs_m.columns:
                part[f"lsst_obs_{c}"] = lsst_obs_m[c].values
            # split roman into truth (already-matched) vs obs columns
            roman_true_cols = [c for c in roman_m.columns if c.startswith("truth_")]
            roman_obs_cols = [c for c in roman_m.columns if not c.startswith("truth_")]
            for c in roman_true_cols:
                part[f"roman_true_{c[len('truth_'):]}"] = roman_m[c].values
            for c in roman_obs_cols:
                part[f"roman_obs_{c}"] = roman_m[c].values
            part["roman_lsst_sep_arcsec"] = d2d.arcsec[sel]
            part["lsst_tract"] = tract
            merged_parts.append(part)
            print(f"  tract {tract}: roman_bbox={len(roman)+0:,} lsst_obj={len(lsst_obs):,} "
                  f"-> {len(part):,} roman-lsst matches "
                  f"({int(true_ok.sum()):,} with lsst-truth)")

        merged = pd.concat(merged_parts, ignore_index=True)
        # An LSST object can match the same Roman det; drop duplicate Roman dets
        # keeping the closest LSST object (one row per Roman detection).
        before = len(merged)
        merged = (merged.sort_values("roman_lsst_sep_arcsec")
                  .drop_duplicates("roman_true_ind", keep="first")
                  .reset_index(drop=True))
        print(f"deduplicated on Roman truth_ind: {before:,} -> {len(merged):,} rows")

        merged.to_parquet(MERGED_OUT, index=False)
        print(f"saved merged table -> {MERGED_OUT} ({len(merged):,} rows, "
              f"{MERGED_OUT.stat().st_size/1e6:.1f} MB)")

    # ----------------------------------------------------------------------- #
    # Stage B: build the env-star classifier from the FULL Roman catalog
    # ----------------------------------------------------------------------- #
    print("fitting the F158 size-envelope classifier on the full Roman catalog ...")
    dset = pads.dataset(str(ROMAN_PARQUET), format="parquet")
    fit_cols = ["matched", "flags", "mag_auto_H158", "magerr_auto_H158",
                "truth_gal_star", "x2win_world_H158", "y2win_world_H158", "xywin_world_H158"]
    fit_cat = dset.to_table(columns=fit_cols,
                            filter=pads.field("matched") == True).to_pandas()
    clf = build_env_classifier(fit_cat)
    classify_star = clf.classify
    print(f"  size-envelope classifier: Delta@21={float(clf.Dfun(21)):.3f} "
          f"Delta@24={float(clf.Dfun(24)):.3f} dex (frozen)")
    del fit_cat

    # ----------------------------------------------------------------------- #
    # Stage C: galaxy misclassification curve (true galaxies, compact, vs mag)
    # ----------------------------------------------------------------------- #
    # True galaxies per Roman truth (gal_star == 0), clean flags, detected (S/N>5
    # is already in the catalog by construction). Dedup on Roman truth_ind.
    is_gal = merged["roman_true_gal_star"].values == 0
    clean = merged["roman_obs_flags"].values < FLAG_CUT
    gal = merged[is_gal & clean].copy()
    gal = gal.sort_values("roman_obs_match_sep_arcsec").drop_duplicates("roman_true_ind")
    print(f"true galaxies (gal_star==0, clean, matched): {len(gal):,}")

    # size used for the compact cut -- the swappable input
    gal_size, size_desc = resolve_galaxy_size(gal, size_source)
    print(f"galaxy size source = {size_desc}")
    finite_sz = np.isfinite(gal_size)
    if finite_sz.sum():
        p = np.percentile(gal_size[finite_sz], [10, 50, 90]).round(3)
        print(f"galaxy size p10/50/90 = {p} arcsec")

    # the size-envelope classifier needs roman_obs_* columns under the bare names
    gal_for_clf = gal.rename(columns={
        "roman_obs_mag_auto_H158": "mag_auto_H158",
        "roman_obs_x2win_world_H158": "x2win_world_H158",
        "roman_obs_y2win_world_H158": "y2win_world_H158",
        "roman_obs_xywin_world_H158": "xywin_world_H158",
    })
    gal_for_clf["size_sb"] = size_sb(gal_for_clf)
    env_star = classify_star(gal_for_clf)

    compact = finite_sz & (gal_size < GAL_SIZE_MAX)
    mag_true = gal["roman_true_mag_H158"].values
    have_mag = np.isfinite(mag_true)
    sel = compact & have_mag
    print(f"compact (size<{GAL_SIZE_MAX}\") true galaxies with truth mag: {int(sel.sum()):,}")

    mg = mag_true[sel]
    fs = env_star[sel]
    n_gal = np.histogram(mg, MAG_BINS)[0]
    n_false = np.histogram(mg[fs], MAG_BINS)[0]
    with np.errstate(invalid="ignore"):
        misclass_eff = np.where(n_gal >= 20, n_false / np.maximum(n_gal, 1), np.nan)

    # ----------------------------------------------------------------------- #
    # Stage D: write the misclassification CSV
    #   columns: delta_mag, mag_F158, classification_eff
    #   mag_F158 = REF_MAGLIM - delta_mag  =>  delta_mag = REF_MAGLIM - mag_F158
    # ----------------------------------------------------------------------- #
    tab = pd.DataFrame({
        "delta_mag": REF_MAGLIM - MAG_MID,
        "mag_F158": MAG_MID,
        "classification_eff": misclass_eff,
    })
    tab = tab[n_gal >= 20].copy().fillna(0.0)

    header = (
        f"Roman DC2 galaxy MISCLASSIFICATION efficiency curve (D2)\n"
        f"fraction of COMPACT true galaxies (Roman truth gal_star==0, size<{GAL_SIZE_MAX}\")\n"
        f"classified as stars by the Roman F158 size-envelope classifier, vs F158 mag.\n"
        f"REF_MAGLIM = {REF_MAGLIM:.4f}  (median of DC2 F158 maglim map "
        f"{MAGLIM_MAP.name});  mag_F158 = REF_MAGLIM - delta_mag\n"
        f"SIZE_SOURCE = {size_source}  ({size_desc})\n"
        f"delta_mag,mag_F158,classification_eff"
    )
    np.savetxt(MISCLASS_CSV, tab.values, delimiter=",", header=header, fmt="%.6f")
    print(f"\nwrote {MISCLASS_CSV} ({len(tab)} rows)")

    # report the curve
    print("\ngalaxy misclassification curve (mag_F158, delta_mag, eff, N_compact_gal):")
    full = pd.DataFrame({"mag_F158": MAG_MID, "delta_mag": REF_MAGLIM - MAG_MID,
                         "eff": misclass_eff, "n_gal": n_gal})
    for _, r in full[full.n_gal >= 20].iterrows():
        print(f"  F158={r.mag_F158:6.3f}  delta={r.delta_mag:+6.3f}  "
              f"eff={r.eff:6.4f}  N={int(r.n_gal):>6}")

    print("\nMatch-stage row counts:")
    if n_roman_seen is not None:
        print(f"  Roman dets read (in LSST bboxes): {n_roman_seen:,}")
        print(f"  Roman dets with truth match (1\"): {n_roman_truth:,}")
        print(f"  LSST objects loaded:              {n_lsst_obj:,}")
    print(f"  merged table rows (Roman det<->LSST obj 1\"): {len(merged):,}")
    print(f"  merged table file: {MERGED_OUT} ({MERGED_OUT.stat().st_size/1e6:.1f} MB)")
    print(f"  misclassification CSV: {MISCLASS_CSV}")


if __name__ == "__main__":
    main()
