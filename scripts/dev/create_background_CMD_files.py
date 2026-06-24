import pandas as pd
from streamobs.background import BackgroundResourceBuilder, BackgroundStorage


df_stars    = pd.read_parquet('/path/to/true_stars.parquet')
df_galaxies = pd.read_parquet('/path/to/true_galaxies.parquet')

# To be changed depending on the survey you want to build the background for
# (e.g. 'lsst', 'euclid', 'roman', etc.)
# The results won't depends on the release. It is just a way to load survey's
# properties (completeness, errors etc), but the magnitudes limits will be set by the user.
builder = BackgroundResourceBuilder(survey_name='lsst', release='yr4')

# Build the background resources for stars and galaxies
builder.build(
    catalog_stars=df_stars,
    catalog_galaxies=df_galaxies,
    bands=('g', 'r'),
    maglim_min=23.5,    # lower end of the magnitude limit grid
    maglim_max=27.0,    # upper end
    maglim_step=0.1,    # step size between grid points
    max_delta=1.0,      # discard pairs with |maglim_b2 - maglim_b1| >= max_delta
    n_bins_color=50,
    n_bins_mag=50,
    color_range=(-1, 2),
    mag_range=(15, 30),
    area_ref_deg2=1.0,   # sky area of the truth catalog in deg²
    source_type='both',
)

# Save the result to disc
storage = BackgroundStorage(base_path='../data/background', survey_name='lsst')
builder.save(storage, source_type='both')

# Note that the ressources are not tracked by git. Thus, if needed, they must be
# added to the data repository (e.g. Zenodo) and downloaded before running the
# code that uses them.
