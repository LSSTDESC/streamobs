import pandas as pd
import numpy as np
from streamobs.background import BackgroundResourceBuilder, BackgroundStorage

star_path = "/pbs/home/m/mpelissi/notebooks/Detection/dc2_manip/catalogs/true_stars_cone_r4.0deg_wpositions.parquet"
df_stars    = pd.read_parquet(star_path)

gals_path = "/pbs/home/m/mpelissi/notebooks/Detection/dc2_manip/catalogs/true_gals_cone_r1.0_sizecut0.3.parquet"
df_galaxies = pd.read_parquet(gals_path)

# To be changed depending on the survey you want to build the background for
# (e.g. 'lsst', 'euclid', 'roman', etc.)
# The results won't depends on the release. It is just a way to load survey's
# properties (completeness, errors etc), but the magnitudes limits will be set by the user.
builder = BackgroundResourceBuilder(survey_name='lsst', release='yr4')

# Build the background resources for stars and galaxies

# Optimized range for stars and galaxies can be different
builder.build(
    catalog_stars=df_stars,
    bands=('g', 'r'),
    maglim_min=23.5,    # lower end of the magnitude limit grid
    maglim_max=27.5,    # upper end
    #maglim_min=24.5,    # lower end of the magnitude limit grid
    #maglim_max=24.7,    # upper end
    maglim_step=0.2,    # step size between grid points
    max_delta=1.0,      # discard pairs with |maglim_b2 - maglim_b1| >= max_delta
    n_bins_color=125,
    n_bins_mag=125,
    color_range=(-0.5, 2.0),
    mag_range=(16.0, 28.0),
    area_ref_deg2=np.pi*(4)**2,   # sky area of the truth catalog in deg²
    source_type='stars',
)

builder.build(
    catalog_galaxies=df_galaxies,
    bands=('g', 'r'),
    maglim_min=23.5,    # lower end of the magnitude limit grid
    maglim_max=27.5,    # upper end
    #maglim_min=24.5,    # lower end of the magnitude limit grid
    #maglim_max=24.7,    # upper end
    maglim_step=0.2,    # step size between grid points
    max_delta=1.0,      # discard pairs with |maglim_b2 - maglim_b1| >= max_delta
    n_bins_color=80,
    n_bins_mag=80,
    color_range=(-1., 2.0),
    mag_range=(20.0, 29.0),
    area_ref_deg2=np.pi*(1)**2,   # sky area of the truth catalog in deg²
    source_type='galaxies',
)

# Save the result to disc
storage = BackgroundStorage(base_path='../../data/background', survey_name='lsst')
builder.save(storage, source_type='both')

# Note that the ressources are not tracked by git. Thus, if needed, they must be
# added to the data repository (e.g. Zenodo) and downloaded before running the
# code that uses them.
