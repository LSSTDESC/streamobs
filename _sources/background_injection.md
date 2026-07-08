# Injection method

The injection method passes a user-supplied truth catalog through the full survey pipeline — the same pipeline used for stream injection. It applies photometric errors and the survey selection function object by object.

**Advantages**
- Full per-object photometric errors in all simulated bands.
- Easy to tune: swap the truth catalog or the survey to test different assumptions.
- Magnitude errors and per-band noise are propagated exactly as for stream stars.

**Limitation** — Requires truth catalogs (true positions and magnitudes) for each population. These are not bundled with `streamobs`.

## Single-survey usage

```python
from streamobs.background import Background

bg = Background(
    survey,
    method='injection',
    source_type='both',           # 'stars', 'galaxies', or 'both'
    bands=('g', 'r'),
    catalog_stars=df_true_stars,
    catalog_galaxies=df_true_galaxies,
)
catalog = bg.generate(phi1_limits=(-20, 20), phi2_limits=(-2, 2))
```

The output is a DataFrame with survey-namespace-prefixed magnitude and flag columns (e.g. `lsst_yr4_r_obs`, `lsst_yr4_flag_observed`), identical in format to a stream injection output.

## Multi-survey usage

Pass a list of two surveys to use different photometric models for each band (e.g. LSST for the color band, Roman for the reference band).

```python
from streamobs.surveys import Survey
from streamobs.background import Background

lsst  = Survey.load('lsst',  release='yr4')
roman = Survey.load('roman', release='dc2')

bg = Background(
    surveys=[lsst, roman],
    method='injection',
    source_type='stars',
    bands=('g', 'F158'),          # lsst covers g, roman covers F158
    catalog_stars=df_true_stars,  # must contain lsst_g_true, lsst_r_true,
                                  # roman_F158_true (+ completeness bands)
)
catalog = bg.generate(phi1_limits=(-20, 20), phi2_limits=(-2, 2))
```

The output has separate columns for each survey namespace (e.g. `lsst_yr4_g_obs`, `roman_dc2_F158_obs`, `lsst_yr4_flag_observed`, `roman_dc2_flag_observed`).

### Required truth columns

The truth catalog must contain a `<namespace>_<band>_true` column for every injected band **and** for each survey's `completeness_band` (if not already one of the injected bands). For example, if LSST's `completeness_band` is `r`, but you only inject `g`, you still need `lsst_r_true` so that detection flags can be computed correctly.
