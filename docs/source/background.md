# Generating background

`streamobs` provides two complementary ways to generate a background catalog
(stars and/or galaxies) for a stream realization.

| Method | Speed | Requirements |
|--------|-------|--------------|
| `'light'` | Fast | Precomputed CMD resource files |
| `'full'` | Slow | True background catalogs |

## Quick start — light method

```python
from streamobs.surveys import Survey
from streamobs.background import Background

survey = Survey.load('lsst', release='yr5')
bg = Background(survey, source_type='both', method='light')

catalog = bg.generate(
    phi1_limits=(-20, 20),
    phi2_limits=(-2, 2),
    gc_frame=frame,          # gala GreatCircleICRSFrame
)
```

## Stars only or galaxies only

```python
bg_stars = Background(survey, source_type='stars', method='light')
bg_gals  = Background(survey, source_type='galaxies', method='light')
```

## Full injection method

Requires a DataFrame of true background objects for each population:

```python
bg = Background(
    survey,
    method='full',
    source_type='both',
    catalog_stars=df_true_stars,
    catalog_galaxies=df_true_galaxies,
)
catalog = bg.generate(phi1_limits=(-20, 20), phi2_limits=(-2, 2))
```

See [build_background_resources.md](build_background_resources.md) to learn
how to build or rebuild the resource files used by the light method.
