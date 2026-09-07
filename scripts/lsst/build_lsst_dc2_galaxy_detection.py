#!/usr/bin/env python
"""LSST DC2 galaxy DETECTION efficiency curve  P(galaxy detected | delta_mag_r).

WHY THIS PRODUCT EXISTS
-----------------------
The joint Roman x Rubin channel is being redefined as

    detected by BOTH surveys  AND  classified point-like by ROMAN

i.e. LSST contributes photometry (the g-F158 colour) and depth, but not morphology.
Expressing that needs the LSST-side probability that a galaxy is *detected*, with no
classification term. No such product existed:

* ``lsst_dc2_galaxy_misclass_cutr.csv`` is P(classified point-like | DETECTED), conditional
  on detection -- it cannot supply the detection term.
* ``lsst_dc2_stellar_efficiency_cutr.csv`` is for STARS and is not a substitute: measured
  over the Roman overlap at F158 < 25.5, galaxies are detected at 0.619 against 0.855 for
  stars (see ``artifacts/joint_detection_correlation.csv``).
* ``Survey.get_gal_misclassification_detection`` returns misclassification x detection
  combined, so the detection term cannot be recovered from the shipped curves either.

Definition (mirrors ``build_lsst_dc2_galaxy_misclass.py`` exactly except for the quantity):
denominator = COMPACT true galaxies (``truth_type == 1``, cosmoDC2 ``size_true < 0.3"``) in
the footprint, detected or not; numerator = those with an LSST object within 1" passing
``clean`` and ``psFlux_r/psFluxErr_r > 5``. Binned on TRUE ``mag_r``;
``delta_mag = mag_r - maglim_r`` with ``maglim_r`` the median of the r-band depth map, the
convention every LSST product here uses.

SIZE CUT: the compact cut is the same one the misclassification curve applies, so both
factors of the contaminant probability refer to the same population. Extended galaxies --
which neither survey would ever classify as a star -- do not dilute the denominator. An
all-galaxy column is emitted for reference; it differs by only ~1-2% because 95.9% of true
galaxies are already compact.

Rows are filtered on the COMPACT count, never the all-galaxy count. There are no compact
galaxies brighter than r ~ 19 (bright galaxies are nearby and large), so filtering on the
all-galaxy count would keep bright bins in which the compact efficiency is undefined, and
zero-filling those would tell the injector that bright galaxies are never detected.

FOOTPRINT CAVEAT
----------------
Built from the truth-centric cache produced by ``scripts/roman/measure_joint_detection.py``,
which covers the 18 LSST tracts overlapping the Roman footprint (~17 deg^2), NOT all 79 DC2
tracts like the misclassification product. That is the region the joint channel actually
operates in, so it is the right footprint for this purpose, but the two products therefore
do not share a footprint. Use ``--all-tracts`` to rebuild over the full DC2 skims (slower;
rebuilds the cache from scratch).

Run (streamobs env):
    /astro/store/shiren/conda-envs/stream_team/envs/streamobs/bin/python \
        scripts/lsst/build_lsst_dc2_galaxy_detection.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import healpy as hp
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OUT_DIR = REPO / "data" / "surveys" / "lsst_dc2"
MAGLIM_MAP = OUT_DIR / "lsst_dc2_maglim_r_nside1024.fits.gz"
CACHE_TRUTH = OUT_DIR / "_cache_lsst_truth_detflag_overlap.parquet"
OUT_CSV = OUT_DIR / "lsst_dc2_galaxy_detection_cutr.csv"

BAND = "r"
TRUTH_GAL = 1
GAL_SIZE_MAX = 0.3
MAG_BINS = np.arange(15.0, 29.0 + 1e-6, 0.25)
MAG_MID = 0.5 * (MAG_BINS[1:] + MAG_BINS[:-1])
N_MIN = 20
EFF_DELTA_MIN = -11.0  # bright/saturation cut, matching the other LSST curves

sys.path.insert(0, str(REPO / "scripts" / "roman"))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--all-tracts",
        action="store_true",
        help="rebuild the truth catalogue over all 79 DC2 tracts instead of reusing the "
        "18-tract Roman-overlap cache (slow)",
    )
    args = ap.parse_args()

    mlm = hp.read_map(str(MAGLIM_MAP))
    good = (mlm != hp.UNSEEN) & np.isfinite(mlm) & (mlm > 0)
    MAGLIM_REF = float(np.median(mlm[good]))
    print(f"reference maglim (r map median) = {MAGLIM_REF:.4f}")

    if args.all_tracts:
        import measure_joint_detection as D  # noqa: E402

        print("rebuilding LSST truth over all tracts (this is slow) ...")
        tru = D.build_lsst_truth_catalog(refresh=True)
        footprint = "all DC2 tracts"
    else:
        if not CACHE_TRUTH.exists():
            raise SystemExit(
                f"{CACHE_TRUTH} not found -- run scripts/roman/measure_joint_detection.py "
                "first, or pass --all-tracts."
            )
        tru = pd.read_parquet(
            CACHE_TRUTH,
            columns=["truth_type", "mag_r", "lsst_detected", "cosmodc2_id", "cosmodc2_hp",
                     "tract"],
        )
        footprint = f"{tru['tract'].nunique()} LSST tracts overlapping the Roman footprint"
    print(f"footprint: {footprint};  {len(tru):,} truth rows")

    gal = tru[(tru["truth_type"].to_numpy() == TRUTH_GAL)
              & np.isfinite(tru["mag_r"].to_numpy())]
    print(f"true galaxies with a finite truth mag_r: {len(gal):,}")

    mag = gal["mag_r"].to_numpy()
    det = gal["lsst_detected"].to_numpy().astype(bool)

    # compact subset, for the reference column only
    import measure_joint_detection as D  # noqa: E402

    size = D.join_true_sizes_full(
        gal["cosmodc2_id"].to_numpy(), gal["cosmodc2_hp"].to_numpy()
    )
    compact = np.isfinite(size) & (size < GAL_SIZE_MAX)
    print(f"compact fraction: {compact.mean():.4f}")

    def curve(mask):
        n = np.histogram(mag[mask], MAG_BINS)[0]
        k = np.histogram(mag[mask & det], MAG_BINS)[0]
        with np.errstate(invalid="ignore"):
            p = np.where(n >= N_MIN, k / np.maximum(n, 1), np.nan)
        return n, p

    n_all, p_all = curve(np.ones(len(gal), bool))
    n_cmp, p_cmp = curve(compact)

    tab = pd.DataFrame(
        {
            "mag_r": MAG_MID,
            "delta_mag": MAG_MID - MAGLIM_REF,
            "detection_eff": p_cmp,          # PRIMARY: compact, matching the misclass curve
            "detection_eff_allgal": p_all,   # reference only
            "n_compact": n_cmp,
        }
    )
    # Filter on the COMPACT count, not the all-galaxy count. There are no compact galaxies
    # brighter than r ~ 19 (bright galaxies are nearby and large), so filtering on n_all
    # would keep bright bins where the compact efficiency is undefined -- and filling those
    # with 0.0 would tell the injector bright galaxies are never detected. compact is a
    # subset of all, so every surviving row has both quantities well defined.
    tab = tab[n_cmp >= N_MIN].copy()
    assert tab[["detection_eff", "detection_eff_allgal"]].notna().all().all(), (
        "NaN survived the n_compact filter -- do not zero-fill, investigate"
    )
    bright = tab["delta_mag"] < EFF_DELTA_MIN
    tab = tab[~bright].copy()
    print(f"dropped {int(bright.sum())} bins with delta_mag < {EFF_DELTA_MIN} (bright cut)")

    header = (
        "LSST DC2 galaxy DETECTION efficiency curve\n"
        f"fraction of COMPACT true galaxies (truth_type==1, cosmoDC2 size_true<{GAL_SIZE_MAX}\"),\n"
        "detected or not, that ARE detected: an LSST object within 1\" passing clean and\n"
        "psFlux_r/psFluxErr_r > 5. NO classification term -- this is the detection factor\n"
        "alone, for a joint selection that takes morphology from Roman only.\n"
        "Same compact population as lsst_dc2_galaxy_misclass_cutr.csv, so the two factors\n"
        "refer to the same objects; extended galaxies, which neither survey would call a\n"
        "star, do not dilute the denominator.\n"
        f"footprint: {footprint}\n"
        f"reference maglim (median of {MAGLIM_MAP.name}) = {MAGLIM_REF:.4f}; "
        "delta_mag = mag_r - maglim.\n"
        "detection_eff_allgal drops the size cut (reference only; differs ~1-2%).\n"
        f"rows filtered on n_compact >= {N_MIN}; no zero-filling.\n"
        "mag_r,delta_mag,detection_eff,detection_eff_allgal,n_compact"
    )
    np.savetxt(
        OUT_CSV,
        tab[["mag_r", "delta_mag", "detection_eff", "detection_eff_allgal", "n_compact"]].values,
        delimiter=",",
        header=header,
        fmt=["%.6f", "%.6f", "%.6f", "%.6f", "%d"],
    )
    print(f"\nwrote {OUT_CSV.relative_to(REPO)} ({len(tab)} rows)")
    show = tab.iloc[np.linspace(0, len(tab) - 1, min(14, len(tab))).astype(int)]
    print(show.to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
