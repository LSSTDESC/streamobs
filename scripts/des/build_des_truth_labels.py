#!/usr/bin/env python
"""Build the DES Y6 Balrog truth star/galaxy label table.

The DES Y6 Balrog files carry **no truth star/galaxy label**.  Anbajagane et al.
2025 (arXiv:2501.05683, Sec. 3.1 footnote 14 and Sec. 4.3) define injected stars
via the *colour-based* classifier of the parent Y3 deep-field catalogue (Hartley
& Choi et al. 2021) -- explicitly "using photometric colors (and not
morphology)".  That classifier lives in the deep-field catalogue's ``KNN_CLASS``
column, which is joinable to Balrog because the injected ``id`` **is** the
deep-field ``COADD_OBJECT_ID``.

This script downloads (or reads) the four Y3 deep-field photometry parquets and
emits one compact lookup table for the reducer's ``--truth-labels`` argument.

Source
------
https://desdr-server.ncsa.illinois.edu/despublic/y3a2_files/deepfields/Y3_DEEP_FIELDS_PHOTOM/
  Y3_DEEP_FIELDS_PHOTOM-{COSMOS,SN-C3,SN-E2,SN-X3}-0000.parquet   (~5.3 GB total)

KNN_CLASS encoding (determined empirically -- the released JSON sidecars carry
only Oracle type info, no value dictionary)
------------------------------------------------------------------------------
Measured on bright (18 < BDF_MAG_I < 22) deep-field objects, where the deep
photometry makes morphology unambiguous:

    class    n      med BDF_T   frac|T|<0.02   med(PSF-BDF)   med J-Ks
    0      19821     1.324         0.149          0.778         0.150
    1      20543     0.891         0.015          0.629         0.505
    2       4411    -0.004         0.924         -0.007        -0.206

  * ``KNN_CLASS == 2`` -> **STAR**.  Point-like (median BDF_T ~ 0, concentration
    PSF-BDF ~ 0) and on the stellar locus in J-Ks.
  * ``KNN_CLASS == 1`` -> **GALAXY**.  Resolved, red in J-Ks.
  * ``KNN_CLASS == 0`` -> **UNCLASSIFIED**, *not* galaxy.  Only 21% of these have
    finite BDF_MAG_J/KS at all, i.e. the colour classifier had no NIR to run on.
    They are 22.5% of injections (every injection has ``in_VHS_footprint == 1``,
    so this is "too faint for usable NIR photometry", not "outside NIR
    coverage").  **They must be dropped from BOTH the numerator and the
    denominator** of the stellar efficiency -- counting them as galaxies would
    bias the galaxy-misclassification curve, and counting them as stars would
    bias everything.
  * ``KNN_CLASS == 3`` -> 2 objects in the whole catalogue; ignored.

Why not a morphological proxy
-----------------------------
A ``|bdf_T| < 0.02`` cut on the injected deep-field morphology (the fallback a
previous iteration used) selects 7.64% of injections but is only **~37% pure**:
per 5M injections it picks up ~141k true stars against ~140k galaxies and ~102k
unclassified.  The colour labels are not a refinement here, they are the
difference between a stellar selection function and a mostly-galaxy one.

Outputs
-------
    <out>/des_y6_deepfield_truth.parquet    ID, KNN_CLASS, TILENAME, BDF_MAG_I
                                            (2,826,988 rows, ~34 MB)
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

BASE_URL = (
    "https://desdr-server.ncsa.illinois.edu/despublic/y3a2_files/"
    "deepfields/Y3_DEEP_FIELDS_PHOTOM"
)
FIELDS = ("COSMOS", "SN-C3", "SN-E2", "SN-X3")
COLUMNS = ["ID", "KNN_CLASS", "TILENAME", "BDF_MAG_I"]

# Semantic constants -- see the module docstring for the evidence.
KNN_UNCLASSIFIED = 0
KNN_GALAXY = 1
KNN_STAR = 2


def field_path(data_dir: Path, field: str) -> Path:
    return data_dir / f"Y3_DEEP_FIELDS_PHOTOM-{field}-0000.parquet"


def fetch(data_dir: Path) -> None:
    """Download any missing deep-field parquet (resumable)."""
    data_dir.mkdir(parents=True, exist_ok=True)
    for field in FIELDS:
        dest = field_path(data_dir, field)
        if dest.exists():
            print(f"  have {dest.name} ({dest.stat().st_size / 1e9:.2f} GB)")
            continue
        url = f"{BASE_URL}/{dest.name}"
        print(f"  fetching {dest.name} ...")
        subprocess.run(
            ["curl", "-L", "--fail", "--retry", "5", "--retry-delay", "10",
             "-C", "-", "-o", str(dest), url],
            check=True,
        )


def build(data_dir: Path, out_dir: Path) -> pd.DataFrame:
    frames = []
    for field in FIELDS:
        path = field_path(data_dir, field)
        tab = pq.ParquetFile(path).read(columns=COLUMNS).to_pandas()
        print(f"  {field:8s} {len(tab):>9,} rows")
        frames.append(tab)
    deep = pd.concat(frames, ignore_index=True)

    n_dup = int(deep["ID"].duplicated().sum())
    if n_dup:
        # The four fields are disjoint; a duplicate would mean the join key is
        # not unique and the reindex-based lookup would silently pick one.
        raise SystemExit(f"deep-field ID is not unique: {n_dup} duplicates")

    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "des_y6_deepfield_truth.parquet"
    deep.to_parquet(out, index=False)

    counts = deep["KNN_CLASS"].value_counts().sort_index()
    print(f"\n  total {len(deep):,} rows, {deep['ID'].nunique():,} unique IDs")
    for cls, n in counts.items():
        label = {
            KNN_UNCLASSIFIED: "unclassified",
            KNN_GALAXY: "galaxy",
            KNN_STAR: "STAR",
        }.get(int(cls), "other")
        print(f"    KNN_CLASS={int(cls)} ({label:12s}) {n:>9,}  {n / len(deep):.4f}")
    print(f"\n  wrote {out}  ({out.stat().st_size / 1e6:.1f} MB)")
    return deep


def audit_join(deep: pd.DataFrame, balrog: Path, n_rows: int) -> None:
    """Cross-check the join and the class encoding against the Balrog file."""
    import h5py

    lut = pd.Series(deep["KNN_CLASS"].to_numpy(), index=deep["ID"].to_numpy())
    with h5py.File(balrog, "r") as fh:
        ids = fh["id"][:n_rows]
        bdf_t = fh["bdf_T"][:n_rows]
        vhs = fh["in_VHS_footprint"][:n_rows]

    knn = lut.reindex(ids).to_numpy()
    matched = np.isfinite(knn)
    print(f"\n  join audit on {len(ids):,} injections")
    print(f"    matched to deep field : {matched.mean():.4f}")
    print(f"    in_VHS_footprint == 1 : {(vhs == 1).mean():.4f}")

    for cls in (KNN_UNCLASSIFIED, KNN_GALAXY, KNN_STAR):
        sel = knn == cls
        if not sel.any():
            continue
        print(
            f"    class {cls}: frac={sel.mean():.4f}  "
            f"frac|bdf_T|<0.02={(np.abs(bdf_t[sel]) < 0.02).mean():.4f}"
        )

    classified = np.isin(knn, [KNN_GALAXY, KNN_STAR])
    print(f"    star fraction among classified: "
          f"{(knn[classified] == KNN_STAR).mean():.4f}")

    # Purity of the morphological proxy this table replaces.
    compact = np.abs(bdf_t) < 0.02
    if compact.any():
        print(f"    |bdf_T|<0.02 selects {compact.mean():.4f} of injections, "
              f"purity for true stars = {(knn[compact] == KNN_STAR).mean():.4f}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--data-dir",
        default="/astro/store/shire/pferguso/projects/delve/stream_obs_update/y3_deep_fields",
        help="where the deep-field parquets live (downloaded here if missing)",
    )
    ap.add_argument("--out", default=None, help="output dir (default: --data-dir)")
    ap.add_argument("--no-fetch", action="store_true", help="fail instead of downloading")
    ap.add_argument(
        "--audit",
        default="/astro/store/shire/pferguso/des_y6_balrog/fiducial_injected_sof.hdf5",
        help="Balrog injected file to cross-check the join against ('' to skip)",
    )
    ap.add_argument("--audit-rows", type=int, default=5_000_000)
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    out_dir = Path(args.out) if args.out else data_dir

    if not args.no_fetch:
        print("deep-field catalogues:")
        fetch(data_dir)
    missing = [f for f in FIELDS if not field_path(data_dir, f).exists()]
    if missing:
        raise SystemExit(f"missing deep-field parquets: {missing}")

    print("\nbuilding truth table:")
    deep = build(data_dir, out_dir)

    if args.audit and os.path.exists(args.audit):
        audit_join(deep, Path(args.audit), args.audit_rows)


if __name__ == "__main__":
    sys.exit(main())
