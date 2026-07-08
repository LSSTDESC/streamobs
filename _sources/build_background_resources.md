# Building background resources

The [light background method](background_light.md) reads precomputed CMD
histograms from parquet files in `data/background/`. These files are not tracked
by git and must be built once per survey (or survey combination) by the package
developer, then included in the data repository.

## Storage path convention

One parquet file per `(source_type, bands)` combination, located at:

```
data/background/{dir_name}/{source_type}_{bands_str}.parquet
```

The directory name and band string are derived from the **canonical sort** of
`(survey.name, band)` pairs (alphabetical), so any input order produces the
same path:

| Survey combination | `bands` parameter | Canonical path |
|---|---|---|
| LSST only | `('g', 'r')` | `lsst/stars_gr.parquet` |
| LSST + Roman | `('g', 'F158')` | `lsst_roman/stars_gF158.parquet` |
| Roman + LSST (reversed) | `('F158', 'g')` | `lsst_roman/stars_gF158.parquet` |

## What you need

- A DataFrame of **true** (pre-observation) positions and magnitudes for each
  population. These are not part of `streamobs`.
- Required truth columns per survey:
  - `<namespace>_<band>_true` for each CMD band (the two bands passed to `build`)
  - `<namespace>_<completeness_band>_true` for each survey's `completeness_band`
    (auto-included in the detection flag; must be present even if not a CMD band)
- Positions (`ra`, `dec`) in the catalog are optional. If absent, they are sampled
  uniformly over a sky area of `area_ref_deg2`.

## CMD histogram format

Each grid point stores a 2-D histogram of detected objects in
`(color, mag_ref_band)` space, **normalized by the reference sky area**
(counts per deg²). At generation time, the density is scaled by the pixel
area to obtain the expected count for each HEALPix pixel.

## Building the grid — single survey

```python
import pandas as pd
from streamobs.surveys import Survey
from streamobs.background import BackgroundResourceBuilder, BackgroundStorage

df_stars    = pd.read_parquet('/path/to/true_stars.parquet')
df_galaxies = pd.read_parquet('/path/to/true_galaxies.parquet')

survey = Survey.load('lsst', release='yr4')
builder = BackgroundResourceBuilder(surveys=survey)

# Build CMD grids for stars and galaxies in separate calls
# (catalogs and mag/color ranges may differ between populations)
builder.build(
    catalog_stars=df_stars,
    bands=('g', 'r'),
    maglim_min=23.5,    # lower end of the magnitude limit grid
    maglim_max=27.5,    # upper end
    maglim_step=0.25,   # step size between grid points
    max_delta=1.0,      # discard pairs with |maglim_b2 - maglim_b1| >= max_delta
    n_bins_color=125,
    n_bins_mag=125,
    color_range=(-0.5, 2.0),
    mag_range=(16.0, 28.0),
    area_ref_deg2=50.3,  # sky area of the truth catalog in deg²
    source_type='stars',
)

builder.build(
    catalog_galaxies=df_galaxies,
    bands=('g', 'r'),
    maglim_min=23.5,
    maglim_max=27.5,
    maglim_step=0.5,
    max_delta=1.0,
    n_bins_color=80,
    n_bins_mag=80,
    color_range=(-1.0, 2.0),
    mag_range=(20.0, 29.0),
    area_ref_deg2=3.1,
    source_type='galaxies',
)

storage = BackgroundStorage(base_path='data/background', survey_name='lsst')
builder.save(storage, source_type='both')
# → data/background/lsst/stars_gr.parquet
# → data/background/lsst/galaxies_gr.parquet
```

## Building the grid — multi-survey

Pass a list of two surveys (one per CMD band). The builder resolves the
canonical order and writes to the matching directory automatically.

```python
from streamobs.surveys import Survey
from streamobs.background import BackgroundResourceBuilder, BackgroundStorage

lsst  = Survey.load('lsst',  release='yr4')
roman = Survey.load('roman', release='dc2')

# Truth catalog must include:
#   lsst_g_true    ← LSST CMD band (color axis)
#   lsst_r_true    ← LSST completeness_band (auto-added even if not in CMD bands)
#   roman_F158_true ← Roman CMD band (magnitude axis) + Roman completeness_band
df_stars = pd.read_parquet('/path/to/true_stars_lsst_roman.parquet')

builder = BackgroundResourceBuilder(surveys=[lsst, roman])
builder.build(
    catalog_stars=df_stars,
    bands=('g', 'F158'),    # lsst covers g, roman covers F158
    maglim_min=23.5,
    maglim_max=27.5,
    maglim_step=0.25,
    max_delta=1.0,
    n_bins_color=125,
    n_bins_mag=125,
    color_range=(-0.5, 2.5),
    mag_range=(16.0, 28.0),
    area_ref_deg2=50.3,
    source_type='stars',
)

storage = BackgroundStorage(base_path='data/background', survey_name='lsst_roman')
builder.save(storage, source_type='stars')
# → data/background/lsst_roman/stars_gF158.parquet
```

## Grid size guidance

| Parameter | Typical value | Effect |
|---|---|---|
| `maglim_step` | 0.25 mag | Smaller → more accurate interpolation, longer build time |
| `max_delta` | 1.0 mag | Keeps only grid points near the diagonal |
| `n_bins_color`, `n_bins_mag` | 125 | Resolution of each CMD histogram |
| `area_ref_deg2` | actual catalog area | Must match the truth catalog's sky coverage; used to normalize CMD to counts/deg² |

The build time scales as `O(N_pairs × N_catalog)`. For a 0.25 mag step grid
over [23.5, 27.5] with `max_delta=1.0` there are roughly 130 pairs per source
type.

## Note on resources

Resources are not tracked by git. Distribute them via the data repository
(e.g. Zenodo) and place them under `data/background/` before using
`Background(..., method='light')`.
