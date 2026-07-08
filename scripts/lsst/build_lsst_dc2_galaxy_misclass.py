#!/usr/bin/env python
# coding: utf-8
"""LSST DC2 galaxy-misclassification curve.

Fraction of COMPACT true galaxies (truth_type==1, cosmoDC2 ``size_true`` < 0.3")
that are DETECTED (S/N>5, clean) yet classified as point sources
(``extendedness`` < 0.5) — i.e. galaxies that leak into a star sample — as a
function of r magnitude.  LSST analog of
``scripts/roman/build_roman_galaxy_misclass.py``, adapted for LSST: the classifier
is the ``extendedness`` flag (not a size envelope), the match is a direct
single-survey object->truth match, and true sizes come from the cosmoDC2
``size_true`` per-healpix skims joined by ``cosmodc2_id``.

Output (mirrors the Roman product; not yet wired into lsst_dc2.yaml — consumption
by the injector is deferred, same as Roman):
  data/surveys/lsst_dc2/lsst_dc2_galaxy_misclass_cutr.csv
    columns: mag_r,delta_mag,missclassification_eff   (delta_mag = mag_r - maglim_r)

Run with the streamobs env:
  conda activate streamobs
  python scripts/lsst/build_lsst_dc2_galaxy_misclass.py [--tracts N] [--refresh]
"""

import argparse
import re
from glob import glob
from pathlib import Path

import astropy.units as u
import healpy as hp
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits

REPO = Path(__file__).resolve().parents[2]
DATA_DIR = Path("/astro/store/shire/stream_team/stream_finding/data/lsst_dc2")
OUT_DIR = REPO / "data/surveys/lsst_dc2"
MAGLIM_MAP = OUT_DIR / "lsst_dc2_maglim_r_nside1024.fits.gz"
CACHE_GAL = OUT_DIR / "_cache_matched_galaxies.parquet"
CACHE_SIZE = OUT_DIR / "_cache_cosmodc2_sizes.parquet"
MISCLASS_CSV = OUT_DIR / "lsst_dc2_galaxy_misclass_cutr.csv"

BAND = "r"
SNR_DEPTH = 5
EXT_CUT = 0.5  # extendedness < cut -> classified as point source (star)
TRUTH_GAL = 1  # truth_type: 1=galaxy, 2=star, 3=SN
MATCH_RADIUS = 1.0  # arcsec
GAL_SIZE_MAX = 0.3  # arcsec: "compact" galaxy cut (matches Roman)
NSIDE = 1024
MAG_BINS = np.arange(15.0, 29.0 + 1e-6, 0.25)
MAG_MID = 0.5 * (MAG_BINS[1:] + MAG_BINS[:-1])
EFF_DELTA_MIN = -11.0  # drop curve rows with delta_mag < this (bright/saturation cut)

OBJ_COLS = [
    "ra",
    "dec",
    f"mag_{BAND}",
    f"psFlux_{BAND}",
    f"psFluxErr_{BAND}",
    "extendedness",
    "clean",
]
TRU_COLS = [
    "ra",
    "dec",
    "truth_type",
    "id",
    "cosmodc2_id",
    "cosmodc2_hp",
    f"mag_{BAND}",
]


def load_cols(path, cols):
    with fits.open(path) as f:
        d = f[1].data
        out = {}
        for c in cols:
            a = np.asarray(d[c])
            if a.dtype.byteorder == ">":
                a = a.astype(a.dtype.newbyteorder("="))
            out[c] = a
        return pd.DataFrame(out)


def build_matched_galaxies(n_tracts=0):
    """Per tract: match objects -> nearest truth (<1"); keep DETECTED (S/N>5, clean)
    matched true-galaxy object rows, one per true galaxy (nearest)."""
    frames = []
    obj_files = sorted(glob(str(DATA_DIR / "dc2_object_run2.2i_dr6_skim_tract_*.fits")))
    if n_tracts:
        obj_files = obj_files[:n_tracts]
    for of in obj_files:
        tract = re.search(r"tract_(\d+)\.fits", of).group(1)
        tf = DATA_DIR / f"dc2_run2.2i_truth_merged_summary_skim_tract_{tract}.fits"
        if not tf.exists():
            print(f"tract {tract}: truth skim missing -- skipped")
            continue
        obj = load_cols(of, OBJ_COLS)
        tru = load_cols(tf, TRU_COLS)
        c_obj = SkyCoord(ra=obj["ra"].values * u.deg, dec=obj["dec"].values * u.deg)
        c_tru = SkyCoord(ra=tru["ra"].values * u.deg, dec=tru["dec"].values * u.deg)
        idx, d2d, _ = c_obj.match_to_catalog_sky(c_tru)
        obj["match_sep"] = d2d.arcsec
        obj["matched"] = d2d.arcsec < MATCH_RADIUS
        obj["truth_type"] = tru["truth_type"].values[idx]
        obj["truth_id"] = tru["id"].values[idx]
        obj["cosmodc2_id"] = tru["cosmodc2_id"].values[idx]
        obj["cosmodc2_hp"] = tru["cosmodc2_hp"].values[idx]
        obj[f"truth_mag_{BAND}"] = tru[f"mag_{BAND}"].values[idx]
        sn = obj[f"psFlux_{BAND}"] / obj[f"psFluxErr_{BAND}"]
        det_ok = (sn > SNR_DEPTH) & obj["clean"].astype(bool)
        g = obj[obj["matched"] & (obj["truth_type"] == TRUTH_GAL) & det_ok].copy()
        # one detection per true galaxy (nearest), like the Roman builder
        g = g.sort_values("match_sep").drop_duplicates("truth_id")
        frames.append(
            g[
                [
                    "truth_id",
                    "cosmodc2_id",
                    "cosmodc2_hp",
                    f"truth_mag_{BAND}",
                    f"mag_{BAND}",
                    "extendedness",
                ]
            ]
        )
        print(f"tract {tract}: {len(g):,} detected matched true galaxies")
    return pd.concat(frames, ignore_index=True).drop_duplicates("truth_id")


def load_cosmodc2_sizes(needed_ids, needed_hps):
    """Join cosmoDC2 ``size_true`` (arcsec) for the needed cosmodc2_ids, reading
    only the per-healpix skims covering the galaxies' cosmodc2_hp pixels."""
    need = set(np.asarray(needed_ids, dtype="i8"))
    frames = []
    for hpix in sorted(set(int(h) for h in needed_hps if np.isfinite(h))):
        f = DATA_DIR / f"cosmoDC2_v1.1.4_image_summary_skim_hp32_{hpix}.fits"
        if not f.exists():
            print(f"  cosmoDC2 hp {hpix}: skim missing -- skipped")
            continue
        d = load_cols(f, ["galaxy_id", "size_true"])
        d = d[d["galaxy_id"].astype("i8").isin(need)]
        frames.append(d)
        print(f"  cosmoDC2 hp {hpix}: {len(d):,} matched galaxy sizes")
    if not frames:
        return pd.DataFrame(columns=["cosmodc2_id", "size_true"])
    out = pd.concat(frames, ignore_index=True).drop_duplicates("galaxy_id")
    return out.rename(columns={"galaxy_id": "cosmodc2_id"})


def main(n_tracts=0, refresh=False):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    mlm = hp.read_map(MAGLIM_MAP, verbose=False)
    MAGLIM_REF = float(np.median(mlm[mlm != hp.UNSEEN]))
    print(f"reference maglim (r map median) = {MAGLIM_REF:.3f}")

    # ---- matched galaxies (cached) --------------------------------------
    if CACHE_GAL.exists() and not refresh and not n_tracts:
        print(f"loading cached galaxy match ({CACHE_GAL.name}; --refresh to rebuild)")
        gal = pd.read_parquet(CACHE_GAL)
    else:
        gal = build_matched_galaxies(n_tracts)
        if not n_tracts:
            gal.to_parquet(CACHE_GAL)
    print(f"{len(gal):,} detected matched true galaxies")

    # ---- cosmoDC2 true sizes (cached) -----------------------------------
    if CACHE_SIZE.exists() and not refresh and not n_tracts:
        print(f"loading cached cosmoDC2 sizes ({CACHE_SIZE.name})")
        sizes = pd.read_parquet(CACHE_SIZE)
    else:
        print("joining cosmoDC2 size_true:")
        sizes = load_cosmodc2_sizes(
            gal["cosmodc2_id"].values, gal["cosmodc2_hp"].values
        )
        if not n_tracts:
            sizes.to_parquet(CACHE_SIZE)
    lut = dict(zip(sizes["cosmodc2_id"].astype("i8"), sizes["size_true"]))
    gal["size_true"] = gal["cosmodc2_id"].astype("i8").map(lut)
    cov = gal["size_true"].notna()
    print(
        f"size coverage: {cov.mean():.1%} of matched galaxies "
        f"({int(cov.sum()):,}/{len(gal):,})"
    )
    if cov.sum():
        p = np.percentile(gal.loc[cov, "size_true"], [10, 50, 90]).round(3)
        print(f"size_true p10/50/90 = {p} arcsec")

    # ---- misclassification curve for COMPACT galaxies -------------------
    compact = gal["size_true"] < GAL_SIZE_MAX
    mag_true = gal[f"truth_mag_{BAND}"].values
    sel = compact.values & np.isfinite(mag_true)
    print(f'compact (size<{GAL_SIZE_MAX}") detected true galaxies: {int(sel.sum()):,}')

    mg = mag_true[sel]
    is_star = gal["extendedness"].values[sel] < EXT_CUT  # classified as point source
    n_gal = np.histogram(mg, MAG_BINS)[0]
    n_false = np.histogram(mg[is_star], MAG_BINS)[0]
    with np.errstate(invalid="ignore"):
        misclass_eff = np.where(n_gal >= 20, n_false / np.maximum(n_gal, 1), np.nan)

    tab = pd.DataFrame(
        {
            "mag_r": MAG_MID,
            "delta_mag": MAG_MID - MAGLIM_REF,
            "missclassification_eff": misclass_eff,
        }
    )
    tab = tab[n_gal >= 20].copy().fillna(0.0)
    # bright cut, matching the stellar efficiency curve convention
    _bright = tab["delta_mag"] < EFF_DELTA_MIN
    tab = tab[~_bright].copy()
    print(
        f"  dropped {int(_bright.sum())} bins with delta_mag < {EFF_DELTA_MIN} (bright cut)"
    )
    header = (
        "LSST DC2 galaxy MISCLASSIFICATION efficiency curve\n"
        f'fraction of COMPACT true galaxies (truth_type==1, cosmoDC2 size_true<{GAL_SIZE_MAX}")\n'
        "that are DETECTED (S/N>5, clean) yet classified as point sources "
        "(extendedness<0.5), vs r mag.\n"
        f"reference maglim (median of {MAGLIM_MAP.name}) = {MAGLIM_REF:.4f}; "
        "delta_mag = mag_r - maglim.\n"
        "mag_r,delta_mag,missclassification_eff"
    )
    np.savetxt(MISCLASS_CSV, tab.values, delimiter=",", header=header, fmt="%.6f")
    print(f"\nwrote {MISCLASS_CSV.relative_to(REPO)} ({len(tab)} rows)")
    print("\ncurve (mag_r, delta_mag, missclass_eff, N_compact):")
    full = pd.DataFrame(
        {
            "mag_r": MAG_MID,
            "delta_mag": MAG_MID - MAGLIM_REF,
            "eff": misclass_eff,
            "n": n_gal,
        }
    )
    for _, r in full[full.n >= 20].iterrows():
        print(
            f"  r={r.mag_r:6.3f}  delta={r.delta_mag:+6.3f}  eff={r.eff:6.4f}  N={int(r.n):>6}"
        )


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--tracts",
        type=int,
        default=0,
        help="limit to first N tracts (0=all; bypasses cache)",
    )
    ap.add_argument(
        "--refresh", action="store_true", help="rebuild the match + size caches"
    )
    args = ap.parse_args()
    main(n_tracts=args.tracts, refresh=args.refresh)
