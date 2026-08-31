#!/usr/bin/env python
"""Validate the DES Y6 EXT_XGB star selection against an external deep truth sample.

Why
---
The DES stellar `classification_eff` shipped by streamobs is derived from Balrog
via a surrogate for `EXT_XGB` plus a per-magnitude deconvolution
(``scripts/des/build_des_xgb_surrogate.py``).  That chain has one residual
assumption: the confusion terms ``a`` and ``b`` are measured on the real
catalogue's *mixed* star+galaxy population, while what is actually needed is the
same quantity conditioned on true stars.  The assumption is second-order, but it
should be checked rather than asserted.

This script performs the check the DES team performed themselves.  Bechtol et
al. 2025 (arXiv:2501.05739) Sec. IV.2 validate their faint-domain classifiers by
matching Y6 Gold at 0.5" against deep data in the SXDF/XMM-LSS field.  We use
the SPLASH-SXDF catalogue (Mehta et al. 2018, ApJS 235, 36), which covers the
same field with far deeper multi-wavelength photometry (HSC + Suprime-Cam +
VIDEO + IRAC) and ships its own star flag, so it can act as ground truth at DES
depth.

What it measures
----------------
Directly on the REAL catalogue, in bins of magnitude:

  * stellar completeness of the `0 <= EXT_XGB <= 1` selection
    = P(selected | true star)
  * stellar contamination  = fraction of the selection that are not true stars

These are exactly the quantities Table A.3 of the Gold paper reports integrated
over two magnitude ranges, so the integrated numbers are an independent check on
this script itself:

    0 <= EXT_XGB <= 1 : 98.0% eff / 4.0% contam over 17.5 <= i <= 22.5
                        94.3% eff / 12.5% contam over 16.5 <= i <= 23.5
    EXT_XGB == 0      : 92.1% / 1.0%   and   79.3% / 1.5%

and the per-magnitude curve is what the Balrog-derived `classification_eff`
should reproduce once detection is divided out.

Result (2026-08-31, SXDF, 0.5" match, 151,580 SPLASH-classified matches)
------------------------------------------------------------------------
**Completeness reproduces the paper**, which is the quantity streamobs ships as
`classification_eff`:

    selection          range        this work    paper (Table A.3)
    0 <= EXT_XGB <= 1  17.5-22.5      0.988          0.980
    0 <= EXT_XGB <= 1  16.5-23.5      0.964          0.943
    EXT_XGB == 0       17.5-22.5      0.947          0.921
    EXT_XGB == 0       16.5-23.5      0.811          0.793

**Contamination does not, and should not be trusted from this comparison**
(0.223 vs 0.040 for EXT_XGB <= 1 over 17.5-22.5).  The cause is a property of
the truth sample, not of the products.  SPLASH's `STAR_FLAG` is a *pure but
incomplete* star selector:

  * flag == 1 sits in a tight point-source locus -- median HSC FLUX_RADIUS 3.07
    (p16-p84 3.00-3.21) at 17.5 < i < 19, versus 6.76 for flag == 0.  So the
    star sample is clean.
  * but the DES-selected objects that SPLASH calls "galaxy" have median
    FLUX_RADIUS 3.351 with p16 = 3.082 -- more than half of them lie ON the
    stellar locus.  They are stars that failed SPLASH's (multi-band, IRAC-
    dependent) star criterion, not galaxies DES misclassified.

A pure-but-incomplete star label gives an unbiased P(selected | star) and an
inflated P(not star | selected).  Completeness is therefore validated and
contamination is not measurable here; the same limitation means this field
cannot validate the galaxy-misclassification product either.  Doing so would
need a *complete* galaxy label -- e.g. the HSC PDR3 concentration
(i_psfflux_mag - i_cmodel_mag) the DES team used for their own Fig. 3.

Caveats this cannot escape
--------------------------
SXDF is one deep field of a few deg^2.  It constrains the *shape* of the
classification efficiency versus magnitude, not per-region variation across the
5000 deg^2 DES footprint, and its seeing/depth are not representative of the
median survey tile.  Treat a disagreement as a flag to investigate, not as a
correction factor to apply blindly.
"""

from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np
import pandas as pd

BAD = -9.0e8
MAG_CAP = 37.0

# SPLASH-SXDF STAR_FLAG encoding, measured on the released catalogue:
#   1   -> star          (6,309 rows)
#   0   -> galaxy        (390,447)
#   -99 -> UNCLASSIFIED  (772,302 -- 66% of the catalogue)
# -99 is *not* "galaxy".  Folding it into the galaxy class would inflate the
# apparent contamination of the stellar selection enormously, since two thirds of
# the field is unclassified.  It is excluded from both the numerator and the
# denominator, exactly as KNN_CLASS == 0 is on the Balrog side.
SPLASH_STAR = 1
SPLASH_GALAXY = 0
SPLASH_UNCLASSIFIED = -99


def load_splash(path: str, radius_pad: float = 0.2):
    """Load the SPLASH-SXDF catalogue, keeping only what the match needs."""
    from astropy.io import fits

    with fits.open(path, memmap=True) as hdul:
        tab = hdul[1]
        n_expected = tab.header["NAXIS2"]
        names = {c.name.upper(): c.name for c in tab.columns}
        need = ["RA", "DEC", "STAR_FLAG"]
        missing = [c for c in need if c not in names]
        if missing:
            raise SystemExit(f"SPLASH is missing {missing}; columns are {sorted(names)[:40]}")
        # FITS is big-endian; pandas/numpy hashing needs native byte order.
        def native(a):
            a = np.asarray(a)
            return a.astype(a.dtype.newbyteorder("=")) if a.dtype.byteorder == ">" else a

        cols = {c: native(tab.data[names[c]]) for c in need}
        for opt in ("ID", "MAG_AUTO_hsc_i", "MAG_AUTO_hsc_g", "CHI_STAR", "ZPHOT"):
            key = opt.upper()
            if key in names:
                cols[opt] = native(tab.data[names[key]])
    d = pd.DataFrame(cols)
    print(f"  SPLASH: {len(d):,} rows (header says {n_expected:,})")
    if len(d) != n_expected:
        raise SystemExit("SPLASH row count disagrees with its header -- truncated file?")
    print(f"  STAR_FLAG values: {pd.Series(d['STAR_FLAG']).value_counts().to_dict()}")
    return d


def load_des_in_field(gold_dir: str, ra0, ra1, dec0, dec1):
    """Read only the des_y6_gold HATS partitions overlapping the field."""
    import pyarrow.parquet as pq

    cols = ["RA", "DEC", "EXT_XGB", "MAG_AUTO_I", "MAG_AUTO_G", "FLAGS_FOREGROUND"]
    frames = []
    for f in sorted(glob.glob(f"{gold_dir}/**/*.parquet", recursive=True)):
        try:
            pf = pq.ParquetFile(f)
            # cheap pre-filter on the partition's own RA/DEC statistics
            md = pf.metadata
            keep = False
            for rg in range(md.num_row_groups):
                st_ra = md.row_group(rg).column(
                    pf.schema_arrow.names.index("RA")).statistics
                st_dec = md.row_group(rg).column(
                    pf.schema_arrow.names.index("DEC")).statistics
                if st_ra is None or st_dec is None:
                    keep = True
                    break
                if (st_ra.max >= ra0 and st_ra.min <= ra1
                        and st_dec.max >= dec0 and st_dec.min <= dec1):
                    keep = True
                    break
            if not keep:
                continue
            t = pf.read(columns=cols).to_pandas()
            t = t[(t.RA >= ra0) & (t.RA <= ra1) & (t.DEC >= dec0) & (t.DEC <= dec1)]
            if len(t):
                frames.append(t)
        except Exception as exc:
            print(f"    skip {Path(f).name}: {type(exc).__name__}")
    if not frames:
        raise SystemExit("no des_y6_gold rows in the requested field")
    d = pd.concat(frames, ignore_index=True)
    for c in ("MAG_AUTO_I", "MAG_AUTO_G"):
        d.loc[(d[c] < BAD) | (d[c] > 90), c] = np.nan
    print(f"  DES: {len(d):,} rows in field")
    return d


def crossmatch(des: pd.DataFrame, splash: pd.DataFrame, radius_arcsec: float):
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    c_des = SkyCoord(des["RA"].to_numpy() * u.deg, des["DEC"].to_numpy() * u.deg)
    c_sp = SkyCoord(splash["RA"].to_numpy() * u.deg, splash["DEC"].to_numpy() * u.deg)
    idx, sep, _ = c_des.match_to_catalog_sky(c_sp)
    ok = sep.arcsec < radius_arcsec
    print(f"  matched {ok.sum():,}/{len(des):,} DES objects within "
          f"{radius_arcsec}\" ({ok.mean():.3f})")
    out = des[ok].copy().reset_index(drop=True)
    out["splash_star"] = splash["STAR_FLAG"].to_numpy()[idx[ok]]
    out["sep"] = sep.arcsec[ok]
    return out


def curves(m: pd.DataFrame, ext_max: int, magcol: str, bins: np.ndarray):
    """Completeness and contamination of 0 <= EXT_XGB <= ext_max vs magnitude."""
    sel = (m["EXT_XGB"] >= 0) & (m["EXT_XGB"] <= ext_max)
    star = m["splash_star"] == SPLASH_STAR
    mag = m[magcol].to_numpy()
    rows = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        b = (mag >= lo) & (mag < hi)
        n_star = int((b & star).sum())
        n_sel = int((b & sel).sum())
        rows.append({
            "mag": 0.5 * (lo + hi),
            "n": int(b.sum()),
            "n_star": n_star,
            "n_selected": n_sel,
            "completeness": float((b & sel & star).sum() / n_star) if n_star else np.nan,
            "contamination": float((b & sel & ~star).sum() / n_sel) if n_sel else np.nan,
        })
    return pd.DataFrame(rows)


def integrated(m: pd.DataFrame, ext_max: int, magcol: str, lo: float, hi: float):
    sel = (m["EXT_XGB"] >= 0) & (m["EXT_XGB"] <= ext_max)
    star = m["splash_star"] == SPLASH_STAR
    b = (m[magcol] >= lo) & (m[magcol] < hi)
    n_star, n_sel = int((b & star).sum()), int((b & sel).sum())
    return {
        "range": f"{lo}-{hi}",
        "n_star": n_star,
        "n_selected": n_sel,
        "efficiency": float((b & sel & star).sum() / n_star) if n_star else np.nan,
        "contamination": float((b & sel & ~star).sum() / n_sel) if n_sel else np.nan,
    }


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument(
        "--splash",
        default="/astro/store/shire/pferguso/projects/delve/stream_obs_update/"
                "splash_sxdf/SPLASH_SXDF_Mehta+_v1.6.fits")
    ap.add_argument(
        "--gold",
        default="/astro/store/shire/hats/catalogs/des/des_y6_gold/des_y6_gold/dataset")
    ap.add_argument("--radius", type=float, default=0.5, help="match radius, arcsec")
    ap.add_argument("--magcol", default="MAG_AUTO_I")
    ap.add_argument("--out", default="artifacts/des_y6")
    args = ap.parse_args()

    print("loading SPLASH-SXDF:")
    sp = load_splash(args.splash)
    ra0, ra1 = float(sp["RA"].min()), float(sp["RA"].max())
    dec0, dec1 = float(sp["DEC"].min()), float(sp["DEC"].max())
    print(f"  field: RA {ra0:.3f}-{ra1:.3f}, Dec {dec0:.3f}-{dec1:.3f}")

    print("\nloading des_y6_gold in field:")
    des = load_des_in_field(args.gold, ra0, ra1, dec0, dec1)

    print("\ncrossmatching:")
    m = crossmatch(des, sp, args.radius)
    m = m[np.isfinite(m[args.magcol])]
    n_before = len(m)
    m = m[m["splash_star"].isin([SPLASH_STAR, SPLASH_GALAXY])].reset_index(drop=True)
    print(f"  dropped {n_before - len(m):,} SPLASH-unclassified (STAR_FLAG=-99) "
          f"matches; they are neither star nor galaxy truth")
    print(f"  usable with finite {args.magcol}: {len(m):,}")
    print(f"  usable (SPLASH-classified): {len(m):,}")
    print(f"  SPLASH-truth star fraction: "
          f"{(m['splash_star'] == SPLASH_STAR).mean():.4f}")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    result = {"radius_arcsec": args.radius, "magcol": args.magcol,
              "n_matched": int(len(m))}

    for ext_max, label in ((1, "complete"), (0, "pure")):
        bins = np.arange(17.0, 25.0 + 1e-9, 0.5)
        c = curves(m, ext_max, args.magcol, bins)
        c.to_csv(out / f"des_y6_splash_validation_ext{ext_max}.csv", index=False)
        print(f"\n=== 0 <= EXT_XGB <= {ext_max} ({label}) ===")
        print(c[c["n_star"] > 20].to_string(index=False,
                                            float_format=lambda v: f"{v:.4f}"))
        result[f"ext_le_{ext_max}"] = {
            "paper_ranges": [integrated(m, ext_max, args.magcol, 17.5, 22.5),
                             integrated(m, ext_max, args.magcol, 16.5, 23.5)],
        }

    print("\n=== against Bechtol et al. Table A.3 ===")
    print("  selection        range        this work            paper")
    ref = {("ext_le_1", "17.5-22.5"): (0.980, 0.040),
           ("ext_le_1", "16.5-23.5"): (0.943, 0.125),
           ("ext_le_0", "17.5-22.5"): (0.921, 0.010),
           ("ext_le_0", "16.5-23.5"): (0.793, 0.015)}
    for key in ("ext_le_1", "ext_le_0"):
        for got in result[key]["paper_ranges"]:
            exp = ref.get((key, got["range"]))
            print(f"  {key:10s} {got['range']:>11s}  "
                  f"eff {got['efficiency']:.3f} cont {got['contamination']:.3f}   "
                  f"eff {exp[0]:.3f} cont {exp[1]:.3f}")

    (out / "des_y6_splash_validation.json").write_text(json.dumps(result, indent=2))
    print(f"\nwrote {out}/des_y6_splash_validation*.csv/.json")


if __name__ == "__main__":
    main()
