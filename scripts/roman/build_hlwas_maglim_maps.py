"""
Build exposure-time-scaled magnitude-limit maps for the Roman HLWAS tiers.

Reads the HLWAS F158 exposure-time HEALPix sparse maps (.hsp) from the
roman_notebooks auxiliary data directory, applies a quasi-depth recipe to
produce maglim maps in the same truth-anchored S/N=5 convention as the DC2
F158 maglim map, and writes the results to data/surveys/roman_hlwas_<tier>/.

Quasi-depth recipe (user-confirmed):
    depth(pix) = median_reported_depth + 1.25 * log10( t(pix) / t_median )

where t_median is the median exposure time over the tier's footprint.  The
resulting maps are expressed in the same convention as the DC2 F158 map
(truth-anchored S/N=5, keyed to the DC2 median of 26.38 mag) so that the
DC2-derived completeness and photometric-error tables apply unchanged.

Input exposure-time maps:
    ~/software/roman_notebooks/notebooks/footprint_visualization/aux_data/
    map_HLWAS-wide_F158.hsp, map_HLWAS-medium_F158.hsp, map_HLWAS-all_F158.hsp

    Source: https://github.com/spacetelescope/roman_notebooks/tree/main/
            notebooks/footprint_visualization/aux_data

Median reported 5σ point-source depths used as reference (from STScI HLWAS
community-survey page,
https://roman-docs.stsci.edu/roman-community-defined-surveys/high-latitude-wide-area-survey):
    wide   = 26.2 AB (F158)
    medium = 26.4 AB (F158)
    all    = 26.2 AB (wide-tier reference; the 'all' map median t equals the
             wide median because the wide tier dominates the pixel count)

Output:
    data/surveys/roman_hlwas_wide/roman_hlwas_wide_maglim_f158_nside1024.fits.gz
    data/surveys/roman_hlwas_medium/roman_hlwas_medium_maglim_f158_nside1024.fits.gz
    data/surveys/roman_hlwas_all/roman_hlwas_all_maglim_f158_nside1024.fits.gz

Usage:
    conda activate streamobs
    python scripts/roman/build_hlwas_maglim_maps.py
"""

import os

import healpy as hp
import healsparse as hsp
import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
AUX_DATA_DIR = os.path.expanduser(
    "~/software/roman_notebooks/notebooks/footprint_visualization/aux_data"
)

# Tier configuration: (hsp filename, STScI 5-sigma reported depth, description)
TIER_CONFIG = {
    "wide": {
        "hsp_file": "map_HLWAS-wide_F158.hsp",
        "reported_depth": 26.2,
        "depth_note": "STScI HLWAS page, wide tier F158 5-sigma depth",
    },
    "medium": {
        "hsp_file": "map_HLWAS-medium_F158.hsp",
        "reported_depth": 26.4,
        "depth_note": "STScI HLWAS page, medium tier F158 5-sigma depth",
    },
    "all": {
        "hsp_file": "map_HLWAS-all_F158.hsp",
        "reported_depth": 26.2,
        "depth_note": (
            "Wide-tier reference depth; 'all' map median t equals the wide median "
            "because the wide tier dominates pixel count (deep+ultradeep are <0.5% "
            "of pixels)"
        ),
    },
}

# Output nside: match DC2 maglim maps
NSIDE_OUT = 1024

# Map dtype: match DC2 maglim maps (big-endian float32)
DTYPE_OUT = np.float32


def load_exptime_map(hsp_path: str) -> np.ndarray:
    """Read a healsparse exposure-time map and return a full RING HEALPix array.

    Unseen pixels are returned as NaN.  The output nside matches the sparse
    map's nside_sparse (typically 4096).

    Parameters
    ----------
    hsp_path : str
        Path to the .hsp file.

    Returns
    -------
    np.ndarray
        Full HEALPix array (RING ordering, float64) with NaN for unseen pixels.
    """
    hsp_map = hsp.HealSparseMap.read(hsp_path)
    nside = hsp_map.nside_sparse
    arr = hsp_map.generate_healpix_map(nside=nside, nest=False)
    # Replace healpy UNSEEN sentinel with NaN
    arr = np.where(arr == hp.UNSEEN, np.nan, arr)
    return arr


def build_maglim_map(
    exptime_arr: np.ndarray,
    reported_depth: float,
    nside_out: int = NSIDE_OUT,
) -> tuple[np.ndarray, float, float]:
    """Apply the quasi-depth recipe and degrade to output nside.

    Recipe:
        depth(pix) = reported_depth + 1.25 * log10( t(pix) / t_median )

    Pixels with t <= 0 or NaN are set to NaN in the output.

    Parameters
    ----------
    exptime_arr : np.ndarray
        Full-sky exposure-time map at the input nside (RING ordering).
    reported_depth : float
        Median 5-sigma reported depth (AB mag) for this tier.
    nside_out : int
        Output nside for the maglim map.

    Returns
    -------
    maglim_out : np.ndarray
        Maglim map at nside_out (RING ordering, float32).  Unobserved pixels
        are NaN.
    t_median : float
        Median exposure time over valid footprint pixels (seconds).
    maglim_median : float
        Median maglim over valid output pixels.
    """
    # Compute t_median over valid (finite, positive) pixels at full resolution
    valid_mask = np.isfinite(exptime_arr) & (exptime_arr > 0)
    t_vals = exptime_arr[valid_mask]
    t_median = float(np.median(t_vals))

    # Apply recipe at full resolution (before degrading) to preserve spatial structure
    maglim_full = np.full_like(exptime_arr, np.nan)
    maglim_full[valid_mask] = reported_depth + 1.25 * np.log10(t_vals / t_median)

    # Degrade to output nside using mean (ud_grade fills unseen with mean of valid
    # sub-pixels; pixels with no valid sub-pixels remain at the fill value).
    # We set NaN pixels to the healpy UNSEEN sentinel before ud_grade and restore
    # NaN afterward, so healpy does not include them in the mean.
    fill_val = hp.UNSEEN
    tmp = np.where(np.isfinite(maglim_full), maglim_full, fill_val)
    tmp_deg = hp.ud_grade(tmp, nside_out=nside_out, order_in="RING", order_out="RING")
    maglim_out = np.where(tmp_deg == fill_val, np.nan, tmp_deg).astype(DTYPE_OUT)

    # Measured output median
    valid_out = np.isfinite(maglim_out)
    maglim_median = float(np.median(maglim_out[valid_out]))

    return maglim_out, t_median, maglim_median


def write_maglim_map(maglim_map: np.ndarray, out_path: str) -> None:
    """Write a maglim map to a gzipped FITS file (healpy RING convention).

    Unobserved pixels (NaN) are stored as healpy UNSEEN so the file is
    portable with standard healpy read_map.

    Parameters
    ----------
    maglim_map : np.ndarray
        Full-sky HEALPix array (RING ordering).
    out_path : str
        Output file path (should end in .fits.gz).
    """
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    # Replace NaN with healpy UNSEEN for on-disk storage
    map_out = np.where(np.isfinite(maglim_map), maglim_map, hp.UNSEEN).astype(DTYPE_OUT)
    hp.write_map(out_path, map_out, overwrite=True)
    print(f"    Written: {out_path}")


def symlink_dc2_tables(tier: str) -> None:
    """Create symlinks from the HLWAS tier data dir to the DC2 CSV files.

    The completeness and photometric-error tables are derived from the DC2 mock
    and shared across all HLWAS tiers.  The survey loader searches for files in
    data/surveys/roman_hlwas_<tier>/ first; symlinks keep a single copy of the
    files under roman_dc2/ while making them visible to the loader.

    Parameters
    ----------
    tier : str
        Tier name (e.g. 'wide', 'medium', 'all').
    """
    dc2_dir = os.path.join(REPO_ROOT, "data", "surveys", "roman_dc2")
    tier_dir = os.path.join(REPO_ROOT, "data", "surveys", f"roman_hlwas_{tier}")
    os.makedirs(tier_dir, exist_ok=True)

    csv_files = [
        "roman_stellar_efficiency_cutf158.csv",
        "roman_photoerror_f158_catalog.csv",
        "roman_photoerror_f158.csv",
    ]
    for csv in csv_files:
        src = os.path.join(dc2_dir, csv)
        dst = os.path.join(tier_dir, csv)
        if os.path.islink(dst):
            pass  # already linked
        elif os.path.exists(dst):
            pass  # real file already there, leave it
        else:
            os.symlink(src, dst)
            print(f"    Symlinked: {csv} -> roman_dc2/{csv}")


def main() -> None:
    print("=" * 70)
    print("Building Roman HLWAS maglim maps")
    print("=" * 70)
    print(f"Output nside: {NSIDE_OUT}")
    print(f"Aux data dir: {AUX_DATA_DIR}")
    print()

    for tier, cfg in TIER_CONFIG.items():
        print(f"--- Tier: {tier} ---")
        hsp_path = os.path.join(AUX_DATA_DIR, cfg["hsp_file"])
        reported_depth = cfg["reported_depth"]

        if not os.path.exists(hsp_path):
            raise FileNotFoundError(
                f"Exposure-time map not found: {hsp_path}\n"
                "Clone https://github.com/spacetelescope/roman_notebooks and ensure "
                "aux_data/ is present at "
                "~/software/roman_notebooks/notebooks/footprint_visualization/aux_data/"
            )

        print(f"  Reading: {cfg['hsp_file']}")
        exptime_arr = load_exptime_map(hsp_path)
        nside_in = hp.get_nside(exptime_arr)
        n_valid_in = int(np.sum(np.isfinite(exptime_arr) & (exptime_arr > 0)))
        print(f"  Input nside: {nside_in}, valid pixels: {n_valid_in:,}")
        print(f"  Reported depth (reference): {reported_depth:.1f} AB")
        print(f"  Depth note: {cfg['depth_note']}")

        maglim_out, t_median, maglim_median = build_maglim_map(
            exptime_arr, reported_depth, nside_out=NSIDE_OUT
        )

        n_valid_out = int(np.sum(np.isfinite(maglim_out)))
        print(f"  t_median over footprint: {t_median:.1f} s")
        print(f"  Output nside: {NSIDE_OUT}, valid pixels: {n_valid_out:,}")
        print(f"  >>> Measured maglim_median (F158, {tier}): {maglim_median:.4f}")

        out_dir = os.path.join(REPO_ROOT, "data", "surveys", f"roman_hlwas_{tier}")
        out_name = f"roman_hlwas_{tier}_maglim_f158_nside{NSIDE_OUT}.fits.gz"
        out_path = os.path.join(out_dir, out_name)
        write_maglim_map(maglim_out, out_path)

        print("  Symlinking DC2 selection-function tables...")
        symlink_dc2_tables(tier)
        print()

    print("=" * 70)
    print("Done.")
    print("=" * 70)


if __name__ == "__main__":
    main()
