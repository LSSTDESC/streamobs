# Building background resources (developer guide)

The **light background method** reads precomputed color–magnitude diagram
(CMD) histograms from parquet files stored in `data/background/`.  These
files are not tracked by git and must be built locally (or provided by the
package maintainers).

## What the resources are

For each combination of ``(source_type, bands)``, one parquet file is
stored with a long-format table:

```
maglim_r | maglim_g | color_center | mag_center | count | n_ref | area_ref_deg2
```

Each row represents one cell of a 2-D CMD histogram at a specific magnitude
limit pair.  The generator looks up the nearest ``(maglim_r, maglim_g)``
entry at runtime and samples from it.

## How to build resources

### 1. Prepare input catalogs

You need a DataFrame of **true** (pre-observation) positions and magnitudes
for stars and/or galaxies.  These catalogs are not part of `streamobs` and
must be supplied by the developer.

```python
import pandas as pd
df_stars = pd.read_parquet('/path/to/true_stars.parquet')
df_galaxies = pd.read_parquet('/path/to/true_galaxies.parquet')
```

### 2. Run the builder

```python
from streamobs.background import BackgroundResourceBuilder, BackgroundStorage

builder = BackgroundResourceBuilder(survey_name='lsst', release='yr5')
builder.build(
    catalog_stars=df_stars,
    catalog_galaxies=df_galaxies,
    bands=('g', 'r'),
    maglim_ref_values=[25.0, 25.5, 26.0, 26.5, 27.0],
    delta_range=(-1.0, 1.0),
    delta_step=0.1,
    source_type='both',
)
```

### 3. Save to disk

```python
storage = BackgroundStorage(survey_name='lsst', release='yr5')
builder.save(storage)
```

Files are written to `data/background/lsst_yr5/` — one parquet per
``(source_type, bands)`` combination, e.g. `stars_gr.parquet`.

## File naming convention

`BackgroundStorage.get_path(source_type, bands)` returns the full path:

```
data/background/{survey_name}_{release}/{source_type}_{bands_str}.parquet
```

Example: `data/background/lsst_yr5/galaxies_gr.parquet`.
