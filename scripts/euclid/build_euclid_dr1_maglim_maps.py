#!/usr/bin/env python
"""Build the euclid/dr1 release: DR1 footprint at the fiducial (Q1) depth.

Euclid DR1 covers ~2100 deg^2 of the Wide Survey. No DR1 catalogue is available
here, but its *input coverage map* is (``cgv_map_dr1input_o13.fits.gz`` from
cosmos.esa.int: a binary nside-8192 NESTED HEALPix mask, 2108.5 deg^2), and DR1
is observed to the same nominal Wide depth as the Q1 deep-field pass. So the
release is the DR1 footprint filled with a *uniform* depth per band -- the
anchored Q1 median from ``euclid_q1_audit.json`` (26.2 / 24.5 / 24.5 / 24.5 for
VIS / Y / J / H) -- and it carries the Q1-derived selection-function tables by
symlink, exactly as lsst_dp2 carries the lsst_dc2 tables and the Roman HLWAS
tiers carry roman_dc2's.

Footprint resolution. The binary mask is degraded to the nside-128 grid every
streamobs release ships at by taking the covered *fraction* of each coarse
pixel (this is what ``make_cache.py`` in the dr1_coverage project already
caches at nside 2048, which is read when present). A coarse pixel counts as
covered when that fraction is >= --min-frac (default 0.5), which preserves the
total area (2108 -> ~2100 deg^2) rather than inflating the footprint by 10-20%
as an any-subpixel rule would at this resolution. Pass --min-frac 0 for the
any-subpixel convention.

Usage
-----
  python scripts/euclid/build_euclid_dr1_maglim_maps.py
  python scripts/euclid/build_euclid_dr1_maglim_maps.py --min-frac 0 --figdir docs/source/_static/euclid_dr1
"""

from __future__ import annotations

import argparse
import json
import os
import pathlib

import healpy as hp
import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[2]
Q1_DIR = REPO / "data/surveys/euclid_q1"
OUT = REPO / "data/surveys/euclid_dr1"
TAG = "euclid_dr1"
NSIDE = 128
COVERAGE_DIR = pathlib.Path("/astro/store/shire/pferguso/projects/euclid/dr1_coverage")
CACHE = COVERAGE_DIR / "cache/euclid_dr1_frac_nside2048.fits"
SOURCE = COVERAGE_DIR / "cgv_map_dr1input_o13.fits.gz"
BANDS = ["VIS", "Y", "J", "H"]
# the Q1-derived tables every Euclid release shares
TABLES = [
    "euclid_q1_stellar_efficiency_cutvis.csv",
    "euclid_q1_galaxy_misclass_cutvis.csv",
    "euclid_q1_photoerror_vis.csv",
    "euclid_q1_photoerror_vis_catalog.csv",
    "euclid_q1_photoerror_vis_nocut.csv",
    "euclid_q1_photoerror_vis_catalog_nocut.csv",
]


def coverage_fraction(nside_out):
    """Covered fraction per nside_out pixel (RING), from the cache or the raw mask."""
    if CACHE.exists():
        frac = hp.read_map(str(CACHE), nest=True, dtype=np.float32)
        nside_in = hp.npix2nside(frac.size)
        print(f"read {CACHE.name} (nside {nside_in}, NESTED, fractional)")
    else:
        from astropy.io import fits

        print(f"streaming {SOURCE.name} (nside 8192 binary mask) -- a few minutes")
        with fits.open(SOURCE, memmap=False) as h:
            d = h[1].data["T"]
            nside_in = 8192
            frac = np.empty(hp.nside2npix(nside_in), dtype=np.float32)
            pos = 0
            for i in range(0, d.shape[0], 8192):
                blk = np.asarray(d[i : i + 8192], dtype=np.float32).ravel()
                frac[pos : pos + blk.size] = blk
                pos += blk.size
    fac = (nside_in // nside_out) ** 2
    # NESTED ordering makes degrading a reshape + mean over child pixels
    frac_out = frac.reshape(-1, fac).mean(axis=1)
    return hp.reorder(frac_out, n2r=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-frac", type=float, default=0.5,
                    help="coarse pixel is covered when its covered fraction is >= this (default 0.5)")
    ap.add_argument("--figdir", default=None, help="write a footprint figure here")
    args = ap.parse_args()

    audit = json.loads((Q1_DIR / "euclid_q1_audit.json").read_text())
    depth = audit["m5_anchored"]
    OUT.mkdir(parents=True, exist_ok=True)

    frac = coverage_fraction(NSIDE)
    pixarea = hp.nside2pixarea(NSIDE, degrees=True)
    print(f"DR1 covered area from fractions: {frac.sum() * pixarea:.1f} deg^2")
    covered = frac >= max(args.min_frac, 1e-9)
    print(f"covered pixels at nside {NSIDE} with frac >= {args.min_frac}: {covered.sum()} = {covered.sum() * pixarea:.1f} deg^2")

    for b in BANDS:
        m = np.where(covered, depth[b], hp.UNSEEN).astype(np.float32)
        fn = OUT / f"{TAG}_maglim_{b.lower()}_nside{NSIDE}.fits.gz"
        hp.write_map(str(fn), m, overwrite=True, dtype=np.float32)
        print(f"  wrote {fn.name}: uniform {depth[b]:.3f}")

    for csv in TABLES:
        src, dst = Q1_DIR / csv, OUT / csv
        if not src.exists():
            raise FileNotFoundError(f"{src} missing -- run build_euclid_selection_function.py first")
        if os.path.islink(dst):
            if os.readlink(dst) != os.path.relpath(src, OUT):
                dst.unlink()
                os.symlink(os.path.relpath(src, OUT), dst)
            print(f"  linked {csv}")
        elif dst.exists():
            raise FileExistsError(f"{dst} exists and is not a symlink -- resolve manually")
        else:
            os.symlink(os.path.relpath(src, OUT), dst)
            print(f"  symlinked {csv} -> ../euclid_q1/{csv}")

    dr1_audit = {
        "tag": TAG,
        "source_coverage_map": str(SOURCE),
        "coverage_nside_native": 8192,
        "nside": NSIDE,
        "min_frac": args.min_frac,
        "area_deg2_fractional": float(frac.sum() * pixarea),
        "area_deg2_shipped": float(covered.sum() * pixarea),
        "depth_uniform": depth,
        "depth_source": "euclid_q1_audit.json m5_anchored (Wide-survey 5 sigma point-source depths)",
        "tables_from": "euclid_q1",
    }
    (OUT / f"{TAG}_audit.json").write_text(json.dumps(dr1_audit, indent=2))
    print(f"  wrote {TAG}_audit.json")

    if args.figdir:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        figdir = pathlib.Path(args.figdir)
        figdir.mkdir(parents=True, exist_ok=True)
        fig = plt.figure(figsize=(11, 5.2))
        shown = np.where(covered, frac, np.nan)
        hp.mollview(shown, fig=fig.number, coord=["C"], title="", cbar=False, flip="astro",
                    badcolor="white", cmap="viridis", min=0, max=1, notext=True)
        hp.graticule(dpar=30, dmer=60, alpha=0.3)
        plt.title(
            f"Euclid DR1 input coverage on the nside-{NSIDE} grid: {covered.sum() * pixarea:.0f} deg$^2$ "
            f"shipped at uniform depth $I_E$ = {depth['VIS']:.1f}, $Y_E J_E H_E$ = {depth['Y']:.1f}\n"
            f"(colour = covered fraction of the coarse pixel; pixels below {args.min_frac:g} are dropped)",
            fontsize=10,
        )
        fig.savefig(figdir / f"{TAG}_footprint.png", dpi=140, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {figdir / f'{TAG}_footprint.png'}")


if __name__ == "__main__":
    main()
