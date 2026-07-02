"""
Build exposure-time-scaled magnitude-limit maps for the Roman HLWAS tiers.

Reads the HLWAS exposure-time HEALPix sparse maps (.hsp) from the
roman_notebooks auxiliary data directory, applies a quasi-depth recipe to
produce maglim maps anchored to the DC2 truth-anchored depth scale per band,
and writes the results to data/surveys/roman_hlwas_<tier>/.

Quasi-depth recipe (Option B — DC2 truth-anchored, exposure-ratio scaling):
    depth(pix) = DC2_REF_DEPTH_band + 1.25 * log10( t(pix) / DC2_REF_EXPTIME )

where:
  DC2_REF_DEPTH_band  = median (over valid pixels) of the DC2 maglim map for
                        that band (truth-anchored S/N=5, read from
                        roman_dc2_maglim_<band_lower>_nside1024.fits.gz at
                        runtime — one reference depth per band).
  DC2_REF_EXPTIME     = 770.0 s  (DC2 HLIS reference per-pixel exposure =
                        5.5 dithers × 140 s/exposure, Troxel et al. 2023
                        arXiv:2209.06829, Sec. 3.1).
  t(pix)              = per-pixel exposure time in seconds, read directly from
                        the .hsp maps (already in seconds; quantized in ~107.5 s
                        single-exposure units).

This anchors all HLWAS (tier, band) combinations identically to the DC2 truth
scale via the exposure-time ratio, without mixing ETC vintage references across
tiers.  The resulting maps are expressed in the same truth-anchored S/N=5
convention as the DC2 maps so that DC2-derived completeness and photometric-
error tables apply unchanged via delta_mag.

Processed (tier, band) combinations:
    wide/F158, medium/F158, all/F158, all/F106

Maps are produced only where BOTH a .hsp map (map_HLWAS-<tier>_<BAND>.hsp)
and a DC2 reference map (roman_dc2_maglim_<band_lower>_nside1024.fits.gz)
exist.  No F129 maps are available for HLWAS; F106 only for the 'all' tier.

Input exposure-time maps:
    ~/software/roman_notebooks/notebooks/footprint_visualization/aux_data/
    map_HLWAS-wide_F158.hsp, map_HLWAS-medium_F158.hsp,
    map_HLWAS-all_F158.hsp, map_HLWAS-all_F106.hsp

    Source: https://github.com/spacetelescope/roman_notebooks/tree/main/
            notebooks/footprint_visualization/aux_data

Output:
    data/surveys/roman_hlwas_wide/roman_hlwas_wide_maglim_f158_nside1024.fits.gz
    data/surveys/roman_hlwas_medium/roman_hlwas_medium_maglim_f158_nside1024.fits.gz
    data/surveys/roman_hlwas_all/roman_hlwas_all_maglim_f158_nside1024.fits.gz
    data/surveys/roman_hlwas_all/roman_hlwas_all_maglim_f106_nside1024.fits.gz

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

# DC2 reference exposure time: 5.5 dithers x 140 s/exposure (Troxel et al. 2023,
# arXiv:2209.06829, Sec. 3.1).
DC2_REF_EXPTIME = 770.0  # seconds

# (tier, band) combinations to process.
# Each entry: (tier, band_upper, note)
# The .hsp file is map_HLWAS-<tier>_<BAND>.hsp (uppercase band in filename).
# The DC2 reference map is roman_dc2_maglim_<band_lower>_nside1024.fits.gz.
TIER_BAND_CONFIGS = [
    (
        "wide",
        "F158",
        "HLWAS wide-tier F158 exposure-time map",
    ),
    (
        "medium",
        "F158",
        (
            "HLWAS medium-tier F158 exposure-time map. "
            "Single-band F158 wide ≈ medium because medium's extra depth "
            "vs wide is in additional HLWAS bands, not F158 alone."
        ),
    ),
    (
        "all",
        "F158",
        (
            "HLWAS all-tiers stacked F158 exposure-time map "
            "(wide + medium + deep + ultra-deep). "
            "Deep/ultradeep pixels (<0.5% of the footprint) yield a deeper tail; "
            "map median stays near the wide/medium F158 level."
        ),
    ),
    (
        "all",
        "F106",
        (
            "HLWAS all-tiers stacked F106 exposure-time map "
            "(wide + medium + deep + ultra-deep). "
            "F106 has a higher extinction coefficient than F158 but otherwise "
            "follows the same DC2-anchored quasi-depth recipe."
        ),
    ),
]

# Output nside: match DC2 maglim maps
NSIDE_OUT = 1024

# Map dtype: match DC2 maglim maps (big-endian float32)
DTYPE_OUT = np.float32


def get_dc2_ref_depth(band: str) -> float:
    """Read the DC2 maglim map for *band* and return the median over valid pixels.

    Parameters
    ----------
    band : str
        Band name (case-insensitive, e.g. 'F158' or 'f158').

    Returns
    -------
    float
        Median maglim over valid pixels.
    """
    dc2_path = os.path.join(
        REPO_ROOT,
        "data",
        "surveys",
        "roman_dc2",
        f"roman_dc2_maglim_{band.lower()}_nside1024.fits.gz",
    )
    if not os.path.exists(dc2_path):
        raise FileNotFoundError(
            f"DC2 reference maglim map not found: {dc2_path}"
        )
    dc2_map = hp.read_map(dc2_path, nest=False, verbose=False)
    valid = dc2_map > hp.UNSEEN + 1
    return float(np.median(dc2_map[valid]))


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
    dc2_ref_depth: float,
    nside_out: int = NSIDE_OUT,
) -> tuple[np.ndarray, float, float]:
    """Apply the Option B quasi-depth recipe and degrade to output nside.

    Recipe (DC2 truth-anchored, exposure-ratio scaling):
        depth(pix) = dc2_ref_depth + 1.25 * log10( t(pix) / DC2_REF_EXPTIME )

    Pixels with t <= 0 or NaN are set to NaN in the output.

    Parameters
    ----------
    exptime_arr : np.ndarray
        Full-sky exposure-time map at the input nside (RING ordering).
        Values are in seconds, used directly.
    dc2_ref_depth : float
        DC2 truth-anchored median depth for the band (AB mag, S/N=5).
    nside_out : int
        Output nside for the maglim map.

    Returns
    -------
    maglim_out : np.ndarray
        Maglim map at nside_out (RING ordering, float32).  Unobserved pixels
        are NaN.
    t_median : float
        Median exposure time over valid footprint pixels (seconds, diagnostic only).
    maglim_median : float
        Median maglim over valid output pixels.
    """
    # Identify valid pixels at full resolution
    valid_mask = np.isfinite(exptime_arr) & (exptime_arr > 0)
    t_vals = exptime_arr[valid_mask]

    # t_median is a diagnostic — does NOT enter the depth formula
    t_median = float(np.median(t_vals))

    # Apply Option B recipe at full resolution (before degrading) to preserve
    # spatial structure.  t_vals are already in seconds — use directly.
    maglim_full = np.full_like(exptime_arr, np.nan)
    maglim_full[valid_mask] = dc2_ref_depth + 1.25 * np.log10(t_vals / DC2_REF_EXPTIME)

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
    print("Building Roman HLWAS maglim maps (Option B: DC2 truth-anchored)")
    print("=" * 70)
    print(f"Output nside: {NSIDE_OUT}")
    print(f"Aux data dir: {AUX_DATA_DIR}")
    print(f"DC2_REF_EXPTIME: {DC2_REF_EXPTIME:.1f} s")
    print(f"Recipe: depth(pix) = DC2_REF_DEPTH_band + 1.25 * log10( t(pix) / {DC2_REF_EXPTIME:.1f} )")
    print()

    processed_tiers = set()

    for tier, band, note in TIER_BAND_CONFIGS:
        print(f"--- Tier: {tier}, Band: {band} ---")

        # Per-band DC2 reference depth (read at runtime)
        dc2_ref_depth = get_dc2_ref_depth(band)
        print(f"  DC2_REF_DEPTH_{band} (read from DC2 {band} map): {dc2_ref_depth:.6f} AB")

        hsp_filename = f"map_HLWAS-{tier}_{band}.hsp"
        hsp_path = os.path.join(AUX_DATA_DIR, hsp_filename)

        if not os.path.exists(hsp_path):
            raise FileNotFoundError(
                f"Exposure-time map not found: {hsp_path}\n"
                "Clone https://github.com/spacetelescope/roman_notebooks and ensure "
                "aux_data/ is present at "
                "~/software/roman_notebooks/notebooks/footprint_visualization/aux_data/"
            )

        print(f"  Reading: {hsp_filename}")
        exptime_arr = load_exptime_map(hsp_path)
        nside_in = hp.get_nside(exptime_arr)
        n_valid_in = int(np.sum(np.isfinite(exptime_arr) & (exptime_arr > 0)))
        print(f"  Input nside: {nside_in}, valid pixels: {n_valid_in:,}")
        print(f"  Note: {note}")

        maglim_out, t_median, maglim_median = build_maglim_map(
            exptime_arr, dc2_ref_depth, nside_out=NSIDE_OUT
        )

        n_valid_out = int(np.sum(np.isfinite(maglim_out)))
        print(f"  t_median over footprint (diagnostic): {t_median:.1f} s")
        print(f"  Output nside: {NSIDE_OUT}, valid pixels: {n_valid_out:,}")
        print(f"  >>> Measured maglim_median ({band}, {tier}): {maglim_median:.4f}")

        out_dir = os.path.join(REPO_ROOT, "data", "surveys", f"roman_hlwas_{tier}")
        out_name = f"roman_hlwas_{tier}_maglim_{band.lower()}_nside{NSIDE_OUT}.fits.gz"
        out_path = os.path.join(out_dir, out_name)
        write_maglim_map(maglim_out, out_path)

        if tier not in processed_tiers:
            print("  Symlinking DC2 selection-function tables...")
            symlink_dc2_tables(tier)
            processed_tiers.add(tier)
        print()

    print("=" * 70)
    print("Done.")
    print("=" * 70)


if __name__ == "__main__":
    main()
