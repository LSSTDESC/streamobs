#!/usr/bin/env python
"""Stage the LSST DP2 depth maps into data/surveys/lsst_dp2/ at a target nside.

DP2 has no injection catalogue, so unlike des/yr6 and delve/dr3_gold there is no
reducer step and nothing to truth-anchor against: the maps ship on their own
native 5-sigma scale, exactly as lsst_dc2 does. All this script does is degrade
the source HealSparse maps to the resolution the release ships at and write them
under a name that states that resolution.

Degrading uses reduction="mean" -- the same reduction MaglimMap.to_healpix
applies in the Balrog reducer, so the two releases treat resolution identically.

Upsampling is refused. A map cannot be made finer than its source, and a
filename claiming otherwise is exactly the defect that shipped the DES maps
labelled nside1024 while holding nside-512 data.

Usage
-----
  python scripts/lsst/build_dp2_maglim_maps.py --nside 128
  python scripts/lsst/build_dp2_maglim_maps.py --nside 128 --src /path/to/maps
"""

import argparse
import pathlib

import healpy as hp
import healsparse as hsp
import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[2]
OUT = REPO / "data/surveys/lsst_dp2"
DEFAULT_SRC = pathlib.Path(
    "/astro/users/pferguso/stream_team/stream_finding/data/rubin_dp2/depth_maps"
)
STEM = "dp2_deepCoadd_psf_maglim_consolidated_map_weighted_mean"
BANDS = ("g", "r")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nside", type=int, default=128,
                    help="nside to ship at (default 128)")
    ap.add_argument("--src", default=str(DEFAULT_SRC),
                    help="directory holding the source .hsp maps")
    ap.add_argument("--src-nside", type=int, default=512,
                    help="nside of the source filenames (default 512)")
    args = ap.parse_args()

    src = pathlib.Path(args.src)
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"source {src}\ntarget {OUT}  (nside {args.nside})\n")

    for b in BANDS:
        f = src / f"{STEM}_{b}_nside_{args.src_nside}.hsp"
        if not f.exists():
            raise SystemExit(f"missing source map: {f}")
        m = hsp.HealSparseMap.read(str(f))
        if m.nside_sparse < args.nside:
            raise SystemExit(
                f"band {b}: source is nside {m.nside_sparse}, cannot upsample to "
                f"{args.nside}. Fetch a finer map instead of relabelling a coarse one."
            )

        v0 = m[m.valid_pixels]
        v0 = v0[np.isfinite(v0)]
        area0 = m.valid_pixels.size * hp.nside2pixarea(m.nside_sparse, degrees=True)

        out = m if m.nside_sparse == args.nside else m.degrade(args.nside,
                                                               reduction="mean")
        v1 = out[out.valid_pixels]
        v1 = v1[np.isfinite(v1)]
        area1 = out.valid_pixels.size * hp.nside2pixarea(out.nside_sparse, degrees=True)

        dst = OUT / f"{STEM}_{b}_nside_{args.nside}.hsp"
        out.write(str(dst), clobber=True)

        print(f"  {b}: nside {m.nside_sparse} -> {out.nside_sparse}")
        print(f"     median {np.median(v0):.3f} -> {np.median(v1):.3f} "
              f"({np.median(v1) - np.median(v0):+.3f})")
        print(f"     area   {area0:.0f} -> {area1:.0f} deg^2")
        print(f"     wrote  {dst.name}  ({dst.stat().st_size / 1024:.0f} KB)")

    # drop any stale maps at other resolutions so the release never carries two
    for old in sorted(OUT.glob(f"{STEM}_*_nside_*.hsp")):
        if f"_nside_{args.nside}.hsp" not in old.name:
            old.unlink()
            print(f"  removed stale {old.name}")

    print(f"\ndone -> {OUT}\nUpdate maglim_map_g / maglim_map_r in "
          "config/surveys/lsst_dp2.yaml to match.")


if __name__ == "__main__":
    main()
