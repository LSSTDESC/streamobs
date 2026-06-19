#!/usr/bin/env python
"""
Stage streamobs *runtime* data products and zip them for upload.

This is the build-side counterpart to ``bin/download_data.py``: it assembles the
small per-survey product files that ``Survey.load`` needs (maglim maps,
completeness / photo-error tables, the shared ebv map, ...) into a clean tree at
``archive/data/`` and writes ``archive/data.zip`` ready to upload (e.g. Zenodo).
Update ``BASE_DATA_URL`` in ``download_data.py`` to the new record afterward.

Why this exists (the symlink trap): on a dev machine ``data/surveys/roman_dc2``
is a *symlink* to the ~13 GB Roman mock, and the HLWAS tier CSVs are symlinks to
the roman_dc2 copies. Zipping those directly would either capture broken links or
pull in gigabytes of derivation source. This script therefore:
  * **dereferences** symlinks (copies the real file contents), and
  * **excludes** large derivation source/intermediate files (parquets, raw mock
    detections/truth, external skims, provenance), keeping only runtime products.

Usage:
    python bin/build_data_archive.py            # stage + zip into archive/
    python bin/build_data_archive.py --no-zip   # stage only
    python bin/build_data_archive.py --list     # dry-run: show what would be included
"""

import argparse
import fnmatch
import os
import shutil
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
OUT = REPO / "archive"
STAGE = OUT / "data"
ZIP_PATH = OUT / "data.zip"

# Directories (path relative to data/) skipped entirely — derivation source /
# intermediates, not needed to load or inject a survey at runtime.
EXCLUDE_DIRS = {
    "surveys/lsst_dc2",  # external LSST DC2 skims (contamination derivation)
    "surveys/roman_dc2/det",  # raw Roman mock SExtractor detections
    "surveys/roman_dc2/truth",  # raw Roman mock truth tiles
    "surveys/roman_hlwas",  # bare placeholder dir: raw HLWAS exptime maps (build
    # *inputs*); the runtime products are the derived
    # maglim maps in roman_hlwas_{wide,medium,all}/
}
# Directory basename globs skipped wherever they appear.
EXCLUDE_DIR_GLOBS = ["*_tiles", "*_tiles_*", "__pycache__"]

# File basename globs skipped — large / derivation / provenance / non-runtime.
EXCLUDE_FILE_GLOBS = [
    "*.parquet",  # det_truth, lsst_matched, cosmodc2_size, truth_stars
    "andy_*.fits",  # reference matched catalogs (4.8 GB)
    "cosmoDC2_*",  # cosmoDC2 size skims (derivation input)
    "dc2_object_*",  # LSST DC2 object skims
    "dc2_run2.2i_truth_*",  # LSST DC2 truth skims
    "*_raw.csv",  # photo-error provenance (raw, pre-afterburner)
    "roman_galaxy_misclass_*.csv",  # analysis output (injector does not consume it)
    "map_HLWAS-*",  # raw HLWAS exposure-time maps (build inputs)
    "*_rough_maglim*",  # intermediate rough maglim maps
    "*.README.md",
    ".DS_Store",
]
# Safety net: warn + skip any single file larger than this that slipped past the
# globs above (all genuine runtime products are < a few MB).
MAX_FILE_MB = 50.0


def _excluded_dir(rel: str) -> bool:
    rel = rel.replace(os.sep, "/")
    if rel in EXCLUDE_DIRS:
        return True
    base = rel.rsplit("/", 1)[-1]
    return any(fnmatch.fnmatch(base, g) for g in EXCLUDE_DIR_GLOBS)


def _excluded_file(name: str) -> bool:
    return any(fnmatch.fnmatch(name, g) for g in EXCLUDE_FILE_GLOBS)


def collect():
    """Walk data/ (following symlinks), return (kept, skipped) lists of (relpath, bytes)."""
    kept, skipped = [], []
    for root, dirs, files in os.walk(DATA, followlinks=True):
        rel_root = os.path.relpath(root, DATA)
        rel_root = "" if rel_root == "." else rel_root
        # prune excluded directories in place (don't descend — avoids the huge dirs)
        dirs[:] = [
            d for d in sorted(dirs) if not _excluded_dir(f"{rel_root}/{d}".lstrip("/"))
        ]
        for fname in sorted(files):
            rel = f"{rel_root}/{fname}".lstrip("/")
            src = Path(root) / fname
            try:
                size = src.stat().st_size  # stat() follows symlinks
            except OSError:
                continue
            if _excluded_file(fname):
                skipped.append((rel, size, "pattern"))
            elif size > MAX_FILE_MB * 1024 * 1024:
                skipped.append((rel, size, f">{MAX_FILE_MB:.0f}MB"))
            else:
                kept.append((rel, size))
    return kept, skipped


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--list", action="store_true", help="dry run: list included/skipped, no copy"
    )
    ap.add_argument(
        "--no-zip", action="store_true", help="stage into archive/data/ but don't zip"
    )
    args = ap.parse_args()

    kept, skipped = collect()
    kept_mb = sum(s for _, s in kept) / 1024 / 1024

    print("=" * 78)
    print(f"Runtime data products under {DATA} (symlinks dereferenced)")
    print("=" * 78)
    for rel, size in kept:
        print(f"  + {rel}  ({size/1024:.0f} KB)")
    print(f"\n  {len(kept)} files, {kept_mb:.1f} MB total")
    if skipped:
        print(f"\n  excluded {len(skipped)} file(s) (derivation/large/provenance):")
        for rel, size, why in skipped:
            print(f"    - {rel}  ({size/1024/1024:.0f} MB, {why})")

    if args.list:
        return 0

    # Stage real copies into archive/data/
    if STAGE.exists():
        shutil.rmtree(STAGE)
    for rel, _ in kept:
        dst = STAGE / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(DATA / rel, dst)  # copy2 follows symlinks -> real file
    print(f"\nStaged {len(kept)} files -> {STAGE}")

    if args.no_zip:
        return 0

    # Zip with a top-level data/ entry so download_data.py's extractall(repo_root) works.
    with zipfile.ZipFile(ZIP_PATH, "w", zipfile.ZIP_DEFLATED) as zf:
        for rel, _ in kept:
            zf.write(STAGE / rel, arcname=f"data/{rel}")
    print(f"Wrote {ZIP_PATH}  ({ZIP_PATH.stat().st_size/1024/1024:.1f} MB)")
    print(
        "\nNext: upload archive/data.zip and update BASE_DATA_URL in bin/download_data.py"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
