#!/usr/bin/env python
"""Apply the tracked bright-end cut to the galaxy-misclassification curves.

The misclassification rate is the fraction of *detected true galaxies* that the
survey's classifier calls a point source. Brightward of some delta_mag the
true-galaxy counts per bin get small and the rate stops being a measurement:
it swings between 0 and tens of percent from one bin to the next. Where that
sets in is survey-specific, so the threshold lives in each pipeline's tracked
corrections YAML rather than being hardcoded here.

Unlike the photo-error curves, the misclassification curves had no ``*_raw.csv``
provenance. This writes one the first time it runs, from the current file, so
the uncut measurement is preserved and the cut stays reversible and re-runnable.

Each pipeline's YAML supplies the rule under a ``<band>_misclass`` key:

    g_misclass:
      cut_bright:
        delta_mag_min: -4.0

The three generators that produce these curves
(scripts/des/balrog_selection_function.py for DES and DELVE,
scripts/lsst/build_lsst_dc2_galaxy_misclass.py,
scripts/roman/build_roman_galaxy_misclass.py) each write the *uncut* curve, so
re-run this afterwards if any of them is regenerated. The DES reducer applies
the cut inline as well, so DES and DELVE stay correct on a plain re-run.

Usage
-----
  python scripts/apply_misclass_cuts.py            # all releases
  python scripts/apply_misclass_cuts.py des_yr6    # one release
  python scripts/apply_misclass_cuts.py --dry-run
"""

import argparse
import pathlib

import numpy as np
import pandas as pd
import yaml

REPO = pathlib.Path(__file__).resolve().parents[1]
DATA = REPO / "data/surveys"

# release -> (misclass filename, corrections YAML, curve id)
TARGETS = {
    "des_yr6": ("des_yr6_galaxy_misclass_cutg.csv",
                "scripts/des/des_photoerror_corrections.yaml", "g_misclass"),
    "delve_dr3_gold": ("delve_dr3_gold_galaxy_misclass_cutg.csv",
                       "scripts/des/delve_photoerror_corrections.yaml", "g_misclass"),
    "lsst_dc2": ("lsst_dc2_galaxy_misclass_cutr.csv",
                 "scripts/lsst/lsst_photoerror_corrections.yaml", "r_misclass"),
    "roman_dc2": ("roman_galaxy_misclass_cutf158.csv",
                  "scripts/roman/roman_photoerror_corrections.yaml", "f158_misclass"),
}


def read_curve(path):
    """Read a misclass CSV under either header convention.

    Returns (DataFrame, header_text, is_commented) so the writer can put the
    header back the way the release already had it.
    """
    comment_hdr = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                comment_hdr = line.lstrip("#").strip()
            else:
                break
    if comment_hdr:
        names = [c.strip() for c in comment_hdr.split(",")]
        return pd.read_csv(path, comment="#", names=names), comment_hdr, True
    df = pd.read_csv(path)
    return df, ",".join(df.columns), False


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("releases", nargs="*", default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    wanted = args.releases or list(TARGETS)
    for rel in wanted:
        if rel not in TARGETS:
            raise SystemExit(f"unknown release {rel!r}; known: {', '.join(TARGETS)}")
        fname, yaml_rel, curve_id = TARGETS[rel]
        path = DATA / rel / fname
        raw = path.with_name(path.stem + "_raw.csv")
        cfg = REPO / yaml_rel

        rules = (yaml.safe_load(open(cfg)) or {}).get(curve_id) or {}
        cut = (rules.get("cut_bright") or {}).get("delta_mag_min")
        if cut is None:
            print(f"  {rel}: no cut_bright rule for {curve_id} in {yaml_rel} — skipped")
            continue
        if not path.exists():
            print(f"  {rel}: {fname} not found — skipped")
            continue

        # first run: the current file IS the uncut measurement, so keep it
        if not raw.exists():
            if args.dry_run:
                print(f"  {rel}: would save provenance -> {raw.name}")
            else:
                raw.write_bytes(path.read_bytes())
                print(f"  {rel}: saved provenance -> {raw.name}")

        src = raw if raw.exists() else path
        df, header, commented = read_curve(src)
        keep = df["delta_mag"] >= cut
        out = df[keep].reset_index(drop=True)
        print(f"  {rel}: cut delta_mag < {cut} — {len(df)} -> {len(out)} rows "
              f"(dropped {int((~keep).sum())})")
        if args.dry_run:
            continue
        # keep whichever header convention this release already used: DES and
        # DELVE write a plain header row, LSST and Roman a '#'-commented one
        np.savetxt(path, out.values, delimiter=",", header=header,
                   comments="# " if commented else "", fmt="%.6f")

    if not args.dry_run:
        print("\nRebuild the archive and the survey figures afterwards.")


if __name__ == "__main__":
    main()
