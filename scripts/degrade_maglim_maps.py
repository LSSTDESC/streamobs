#!/usr/bin/env python
"""Degrade shipped maglim maps to a coarser nside, in place.

For releases whose depth maps come from a reducer that can be re-run at the
target resolution (des/yr6, delve/dr3_gold), prefer re-running it: the reducer
degrades the *input* map, then applies the anchor shift, then masks outliers at
the target nside, which is not the same operation as degrading the finished
product.

This script is for the releases where re-running is impractical -- lsst/dc2 and
the Roman releases, whose generators need 100+ GB of derivation inputs. Their
shipped maps are already final (scaled, anchored where applicable, masked), so
degrading them is a pure spatial average and the only thing that changes is
resolution.

The average is partial-coverage aware: a coarse pixel takes the mean of its
*valid* subpixels and is only invalid if it has none. That is the same
reduction HealSparse's degrade() applies and the same one the Balrog reducer
uses, so every release ends up treated identically.

Upsampling is refused. A map cannot be made finer than its source, and a
filename claiming otherwise is exactly the defect that once shipped the DES
maps labelled nside1024 while holding nside-512 data.

Usage
-----
  python scripts/degrade_maglim_maps.py --nside 128 lsst_dc2 roman_dc2
  python scripts/degrade_maglim_maps.py --nside 128 --dry-run lsst_dc2
"""

import argparse
import pathlib
import re

import healpy as hp
import healsparse as hsp
import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[1]
DATA = REPO / "data/surveys"


def degrade_one(path, nside_out, dry_run=False):
    """Degrade one maglim map to nside_out. Returns (old, new) stats."""
    name = path.name
    m = re.search(r"nside_?(\d+)", name)
    if not m:
        raise SystemExit(f"{name}: filename does not state an nside")
    claimed = int(m.group(1))

    if path.suffix == ".hsp":
        smap = hsp.HealSparseMap.read(str(path))
    else:
        dense = hp.read_map(str(path), dtype=np.float64)
        good = np.isfinite(dense) & (dense > 0) & (dense != hp.UNSEEN)
        clean = np.where(good, dense, hp.UNSEEN)
        # healpy maps are RING; HealSparseMap(healpix_map=...) assumes NEST by
        # default, so it MUST be told otherwise. Getting this wrong scatters the
        # footprint across the sky while leaving the pixel count and the median
        # depth almost unchanged -- it does not look wrong in a summary.
        smap = hsp.HealSparseMap(
            healpix_map=clean, nside_coverage=32, sentinel=hp.UNSEEN, nest=False
        )

    nside_in = smap.nside_sparse
    if nside_in != claimed:
        print(f"    NOTE {name}: filename says {claimed}, map is {nside_in}")
    if nside_in < nside_out:
        raise SystemExit(
            f"{name}: source is nside {nside_in}, refusing to upsample to "
            f"{nside_out}. Regenerate from a finer source instead."
        )

    v0 = smap[smap.valid_pixels]
    v0 = v0[np.isfinite(v0)]
    area0 = smap.valid_pixels.size * hp.nside2pixarea(nside_in, degrees=True)

    out = smap if nside_in == nside_out else smap.degrade(nside_out, reduction="mean")
    v1 = out[out.valid_pixels]
    v1 = v1[np.isfinite(v1)]
    area1 = out.valid_pixels.size * hp.nside2pixarea(nside_out, degrees=True)

    # The footprint must stay put. A pixel-ordering mistake relocates it while
    # barely changing the pixel count or the median, so it does not show up in
    # the summary statistics. Test it exactly: every valid input pixel has one
    # parent at the coarser resolution, so the output's valid set must be
    # precisely the set of those parents. (A centroid comparison is not usable
    # here -- the HLWAS wide footprint wraps RA 0->360, so its circular mean is
    # meaningless.)
    if nside_in != nside_out:
        levels = int(np.log2(nside_in // nside_out))
        expected = np.unique(smap.valid_pixels // (4**levels))
        got = np.unique(out.valid_pixels)
        if not np.array_equal(expected, got):
            raise SystemExit(
                f"{name}: degraded footprint is not the parent set of the input "
                f"({got.size:,} pixels vs {expected.size:,} expected, "
                f"{np.intersect1d(expected, got).size:,} shared). This is a "
                "pixel-ordering bug, not a resolution effect. Refusing to write."
            )

    new_name = re.sub(
        r"nside_?\d+",
        lambda mm: mm.group(0).replace(str(claimed), str(nside_out)),
        name,
    )
    dst = path.parent / new_name

    print(f"    {name}")
    print(
        f"      nside {nside_in} -> {nside_out} | median "
        f"{np.median(v0):.3f} -> {np.median(v1):.3f} "
        f"({np.median(v1) - np.median(v0):+.3f}) | area "
        f"{area0:,.0f} -> {area1:,.0f} deg^2"
    )
    if dry_run:
        print(f"      would write {new_name}")
        return

    if path.suffix == ".hsp":
        out.write(str(dst), clobber=True)
    else:
        dense_out = out.generate_healpix_map(nside=nside_out, nest=False)
        dense_out = np.where(np.isfinite(dense_out), dense_out, hp.UNSEEN)
        hp.write_map(str(dst), dense_out, overwrite=True, dtype=np.float64)
    print(f"      wrote {new_name}  ({dst.stat().st_size / 1024:.0f} KB)")
    if dst != path:
        path.unlink()
        print(f"      removed {name}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "releases", nargs="+", help="release directory names under data/surveys/"
    )
    ap.add_argument("--nside", type=int, default=128, help="target nside (default 128)")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    for rel in args.releases:
        d = DATA / rel
        if not d.is_dir():
            raise SystemExit(f"no such release directory: {d}")
        maps = sorted(
            p
            for p in d.glob("*maglim*")
            if p.suffix in (".gz", ".hsp") and "5_sig" not in p.name
        )
        print(f"  {rel}: {len(maps)} map(s)")
        for p in maps:
            degrade_one(p, args.nside, args.dry_run)
    print(
        "\nUpdate the maglim_map_* filenames in the matching "
        "config/surveys/*.yaml, then rebuild the figures and the archive."
    )


if __name__ == "__main__":
    main()
