# Background generation

`streamobs` models two background populations: **stars** (field stars passing stellar selection) and **misclassified galaxies** (galaxies passing a star–galaxy separator). Both can be generated together or independently via `source_type`.

Two methods are available:

| Method | Speed | What you need | Page |
|--------|-------|---------------|------|
| `'light'` *(default)* | Fast | Precomputed CMD resource files *(included with other `StreamObs` data)* | [background_light](background_light.md) |
| `'injection'` | Slow | True background catalogs | [background_injection](background_injection.md) |

## Quick start

```python
from streamobs.surveys import Survey
from streamobs.background import Background

survey = Survey.load('lsst', release='yr4')
bg = Background(survey, source_type='both', method='light')   # default

catalog, meta = bg.generate(
    phi1_limits=(-20, 20),
    phi2_limits=(-2, 2),
    gc_frame=frame,   # gala GreatCircleICRSFrame
)
```

Use `source_type='stars'` or `source_type='galaxies'` to restrict to one population.
