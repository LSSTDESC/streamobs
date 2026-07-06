# Background generation

`streamobs` models two background populations: **stars** (field stars passing stellar selection) and **misclassified galaxies** (galaxies passing a star–galaxy separator). Both can be generated together or independently via `source_type`.

Two methods are available:

| Method | Speed | What you need | Page |
|--------|-------|---------------|------|
| `'light'` *(default)* | Fast | Precomputed CMD resource files *(included with other `StreamObs` data)* | [background_light](background_light.md) |
| `'injection'` | Slow | True background catalogs | [background_injection](background_injection.md) |

## Quick start — single survey

```python
from streamobs.surveys import Survey
from streamobs.background import Background

survey = Survey.load('lsst', release='yr4')
bg = Background(survey, source_type='both', method='light', bands=('g', 'r'))

catalog, meta = bg.generate(
    phi1_limits=(-20, 20),
    phi2_limits=(-2, 2),
    gc_frame=frame,   # gala GreatCircleICRSFrame
)
```

Use `source_type='stars'` or `source_type='galaxies'` to restrict to one population.

## Quick start — multi-survey

Two different surveys can cover the two CMD bands (one band each).  Pass them as a list; `streamobs` sorts them into canonical order internally.

```python
from streamobs.surveys import Survey
from streamobs.background import Background

lsst  = Survey.load('lsst',  release='yr4')
roman = Survey.load('roman', release='dc2')

bg = Background(
    surveys=[lsst, roman],
    source_type='stars',
    method='light',
    bands=('g', 'F158'),   # lsst covers g, roman covers F158
)

catalog, meta = bg.generate(
    phi1_limits=(-20, 20),
    phi2_limits=(-2, 2),
    gc_frame=frame,
)
```

The resource files must be built for each survey combination before use.
See [build_background_resources](build_background_resources.md) for details.
