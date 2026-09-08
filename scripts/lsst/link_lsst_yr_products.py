#!/usr/bin/env python
"""Symlink the LSST DC2 selection-function tables into the dependent LSST releases.

The extrapolated year releases (lsst_yr1..lsst_yr5) and the DP2 release use the
DC2-derived efficiency / misclassification / photo-error tables and differ only
in their depth maps — RubinSim baseline maglim maps per year for lsst_yr*, and
DP2's own measured 5-sigma deep-coadd maps for lsst_dp2. Same convention as the
Roman HLWAS tiers, whose tables are symlinked from roman_dc2 by
scripts/roman/build_hlwas_maglim_maps.py.

Note that lsst_dp2 is real commissioning data carrying a simulation-derived
selection function, which is a weaker assumption than for the year releases --
see config/surveys/lsst_dp2.yaml.

Symlinks keep a single copy of each table on dev machines;
bin/build_data_archive.py dereferences them when staging data.zip, so
downloaded data packages contain real per-survey copies.

Run from the repo root (any env with Python 3):
  python scripts/lsst/link_lsst_yr_products.py
"""

import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DC2_DIR = REPO / "data/surveys/lsst_dc2"

TABLES = [
    "lsst_dc2_stellar_efficiency_cutr.csv",
    "lsst_dc2_galaxy_misclass_cutr.csv",
    "lsst_dc2_photoerror_r.csv",
    "lsst_dc2_photoerror_r_catalog.csv",
    # forced-photometry (non-reference band) curves -- required by every tier
    "lsst_dc2_photoerror_r_nocut.csv",
    "lsst_dc2_photoerror_r_catalog_nocut.csv",
]
TARGETS = [f"lsst_yr{y}" for y in (1, 2, 3, 4, 5)] + ["lsst_dp2"]


def main():
    for target in TARGETS:
        yr_dir = REPO / f"data/surveys/{target}"
        yr_dir.mkdir(parents=True, exist_ok=True)
        print(f"{target}:")
        for csv in TABLES:
            src = DC2_DIR / csv
            dst = yr_dir / csv
            if not src.exists():
                raise FileNotFoundError(
                    f"{src} missing — run scripts/lsst/create_streamobs_files_lsst_dc2.py "
                    "and build_lsst_dc2_galaxy_misclass.py first"
                )
            if os.path.islink(dst):
                print(f"    already linked: {csv}")
            elif dst.exists():
                raise FileExistsError(f"{dst} exists and is not a symlink — resolve manually")
            else:
                # relative link: survives moving/renaming the repo checkout
                os.symlink(os.path.relpath(src, yr_dir), dst)
                print(f"    symlinked: {csv} -> ../lsst_dc2/{csv}")


if __name__ == "__main__":
    main()
