"""Script to build precomputed CMD histogram resources for the light background method.

Run once per survey (or survey combination) on a cluster.  The output parquet
files are not tracked by git — distribute them via the project data repository
(e.g. Zenodo) and place them under ``data/background/`` before using
``Background(..., method='light')``.

Storage paths follow the canonical convention:
  data/background/{dir_name}/{source_type}_{bands_str}.parquet

where {dir_name} is the underscore-joined sorted survey names and {bands_str}
is the sorted concatenated band names, e.g.:
  - Single survey LSST g×r :  data/background/lsst/stars_gr.parquet
  - LSST g × Roman F158   :  data/background/lsst_roman/stars_gF158.parquet
"""

import numpy as np
import pandas as pd

from streamobs.background import BackgroundResourceBuilder, BackgroundStorage
from streamobs.surveys import Survey

# ---------------------------------------------------------------------------
# Load truth catalogs (edit paths for your environment)
# ---------------------------------------------------------------------------

star_path = "/pbs/home/m/mpelissi/notebooks/Detection/dc2_manip/catalogs/true_stars_cone_r4.0deg_wpositions.parquet"
df_stars = pd.read_parquet(star_path)

gals_path = "/pbs/home/m/mpelissi/notebooks/Detection/dc2_manip/catalogs/true_gals_cone_r1.0_sizecut0.3.parquet"
df_galaxies = pd.read_parquet(gals_path)

# ---------------------------------------------------------------------------
# Single-survey example — LSST g × r
#
# Load the survey whose photometric model (completeness curves, error model,
# extinction coefficients) should be applied to the truth catalog.  The
# magnitude limits are overridden grid-point by grid-point; the survey release
# only affects the selection-function shape.
#
# Required truth columns:
#   lsst_g_true, lsst_r_true          ← CMD bands
#   lsst_r_true                        ← completeness_band (already present)
# ---------------------------------------------------------------------------

survey_lsst = Survey.load("lsst", release="yr4")

builder = BackgroundResourceBuilder(surveys=survey_lsst)

# Build CMD grids for stars (optimised mag range/bins for stellar locus)
builder.build(
    catalog_stars=df_stars,
    bands=("g", "r"),
    maglim_min=23.5,       # lower end of the magnitude limit grid
    maglim_max=27.5,       # upper end
    maglim_step=0.25,      # step size between grid points
    max_delta=1.0,         # discard pairs with |maglim_b2 - maglim_b1| >= max_delta
    n_bins_color=125,
    n_bins_mag=125,
    color_range=(-0.5, 2.0),
    mag_range=(16.0, 28.0),
    area_ref_deg2=np.pi * 4**2,   # sky area of the truth catalog in deg²
    source_type="stars",
)

# Build CMD grids for galaxies (separate catalog, different mag/color ranges)
builder.build(
    catalog_galaxies=df_galaxies,
    bands=("g", "r"),
    maglim_min=23.5,
    maglim_max=27.5,
    maglim_step=0.5,
    max_delta=1.0,
    n_bins_color=80,
    n_bins_mag=80,
    color_range=(-1.0, 2.0),
    mag_range=(20.0, 29.0),
    area_ref_deg2=np.pi * 1**2,   # sky area of the galaxy truth catalog
    source_type="galaxies",
)

# Save — canonical path: data/background/lsst/stars_gr.parquet, galaxies_gr.parquet
storage = BackgroundStorage(base_path="../../data/background", survey_name="lsst")
builder.save(storage, source_type="both")

# ---------------------------------------------------------------------------
# Multi-survey example — LSST g × Roman F158
#
# Pass two Survey instances (one per CMD band).  The builder sorts them into
# canonical order automatically; storage lands under the canonical directory.
#
# Required truth columns:
#   lsst_g_true                        ← LSST CMD band (color axis)
#   lsst_r_true                        ← LSST completeness_band (auto-included)
#   roman_F158_true                    ← Roman CMD band (magnitude axis)
#   roman_F158_true                    ← Roman completeness_band (already present)
#
# Uncomment and edit paths to run.
# ---------------------------------------------------------------------------

# survey_roman = Survey.load("roman", release="dc2")
#
# builder_multi = BackgroundResourceBuilder(surveys=[survey_lsst, survey_roman])
# builder_multi.build(
#     catalog_stars=df_stars,     # must contain lsst_g_true, lsst_r_true, roman_F158_true
#     bands=("g", "F158"),        # lsst covers g, roman covers F158
#     maglim_min=23.5,
#     maglim_max=27.5,
#     maglim_step=0.25,
#     max_delta=1.0,
#     n_bins_color=125,
#     n_bins_mag=125,
#     color_range=(-0.5, 2.5),
#     mag_range=(16.0, 28.0),
#     area_ref_deg2=np.pi * 4**2,
#     source_type="stars",
# )
#
# # Canonical path: data/background/lsst_roman/stars_gF158.parquet
# storage_multi = BackgroundStorage(
#     base_path="../../data/background", survey_name="lsst_roman"
# )
# builder_multi.save(storage_multi, source_type="stars")

# Note: resources are not tracked by git.  Distribute via the data repository
# (e.g. Zenodo) and download before using Background(..., method='light').
