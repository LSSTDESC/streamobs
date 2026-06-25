# Building background resources

The [light background method](background_light.md) reads precomputed CMD
histograms from parquet files in `data/background/`. These files are not tracked
by git and must be built once per survey by the package developer, and then
included in the data base.

## File format

One parquet file per `(source_type, bands)` combination (e.g. `stars_gr.parquet`), located at:

```
data/background/{survey_name}/{source_type}_{bands}.parquet
```

Each row corresponds to one `(maglim_b2, maglim_b1)` grid point (`b2 = bands[1]` = reference band, `b1 = bands[0]` = color band) and stores the full 2-D CMD histogram in compressed form.

## What you need

- A DataFrame of **true** (pre-observation) positions and magnitudes for each population. These are not part of `streamobs`.
- Required columns: `ra`, `dec`, `<survey>_<band>_true` for each band (e.g. `lsst_g_true`, `lsst_r_true`).
- A loaded `Survey` object for the target survey/release.

## Building the grid

```python
import pandas as pd
from streamobs.background import BackgroundResourceBuilder, BackgroundStorage

df_stars    = pd.read_parquet('/path/to/true_stars.parquet')
df_galaxies = pd.read_parquet('/path/to/true_galaxies.parquet')

builder = BackgroundResourceBuilder(survey_name='lsst', release='yr4')
builder.build(
    catalog_stars=df_stars,
    catalog_galaxies=df_galaxies,
    bands=('g', 'r'),
    maglim_min=23.5,    # lower end of the magnitude limit grid
    maglim_max=27.5,    # upper end
    maglim_step=0.2,    # step size between grid points
    max_delta=1.0,      # discard pairs with |maglim_b2 - maglim_b1| >= max_delta
    n_bins_color=125,
    n_bins_mag=125,
    color_range=(-0.5, 2.0),
    mag_range=(16.0, 28.0),
    area_ref_deg2=1.0,   # sky area of the truth catalog in deg²
    source_type='both',
)
```

Each `(maglim_b2, maglim_b1)` pair within `max_delta` is injected independently into a uniform copy of the survey (no dust, constant magnitude limits). The result is a 2-D histogram of detected objects in `(color, mag_ref)` space.

## Saving to disk

```python
storage = BackgroundStorage(base_path='data/background', survey_name='lsst')
builder.save(storage, source_type='both')
```

This writes `data/background/lsst/stars_gr.parquet` and `data/background/lsst/galaxies_gr.parquet`.

## Grid size guidance

| Parameter | Typical value | Effect |
|---|---|---|
| `maglim_step` | 0.1 mag | Smaller → more accurate interpolation, longer build time |
| `max_delta` | 1.0 mag | Keeps only grid points near the diagonal |
| `n_bins_color`, `n_bins_mag` | 50 | Resolution of each CMD histogram |
| `area_ref_deg2` | 1 deg² | Larger → lower Poisson noise in each histogram cell |

The build time scales as `O(N_pairs × N_catalog)`. For a 0.5 mag step grid over [23.5, 27.0] with `max_delta=1.0` there are roughly 50 pairs per source type.
