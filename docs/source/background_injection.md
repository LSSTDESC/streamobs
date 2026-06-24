# Injection method

The injection method passes a user-supplied truth catalog through the full survey pipeline — the same pipeline used for stream injection. It applies photometric errors, and the survey selection function object by object.

**Advantages**
- Full per-object photometric errors in all simulated bands.
- Easy to tune: swap the truth catalog or the survey to test different assumptions.
- Magnitude errors and per-band noise are propagated exactly as for stream stars.

**Limitation** — Requires truth catalogs (true positions and magnitudes) for each population. These are not bundled with `streamobs`.

## Usage

```python
from streamobs.background import Background

bg = Background(
    survey,
    method='injection',
    source_type='both',          # 'stars', 'galaxies', or 'both'
    bands=('g', 'r'),
    catalog_stars=df_true_stars,
    catalog_galaxies=df_true_galaxies,
)
catalog = bg.generate(phi1_limits=(-20, 20), phi2_limits=(-2, 2))
```

The output is a DataFrame with survey-namespaced magnitude and flag columns (e.g. `lsst_yr4_r_obs`, `lsst_yr4_flag_observed`), identical in format to a stream injection output.
