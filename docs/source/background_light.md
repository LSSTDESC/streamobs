# Light method

The light method generates background catalogs directly from precomputed color–magnitude diagram (CMD) grids without re-running the injection pipeline. It is the default method and is recommended for most use cases.

## How it works

For each HEALPix pixel in the requested sky region:

1. **Effective magnitude limit** — The observed magnitude limit combined with
   dust absorption to obtain the effective magnitude limit:

   $$m_{\mathrm{eff}} = m_{\mathrm{lim}}(\mathrm{pixel}) - A_{\mathrm{band}} \cdot E(B-V)(\mathrm{pixel})$$

   This maps the real per-pixel depth to the nearest pre-built grid point.
   When two surveys cover the two bands, each survey's depth and dust map
   are queried independently for its own band.

2. **Bilinear CMD interpolation** — The CMD histogram grid is stored at a discrete 2-D lattice of `(maglim_b2, maglim_b1)` pairs (reference band × color band). The generator bilinearly interpolates the four surrounding grid points to obtain a CMD at the exact effective magnitude limits of the pixel.

3. **Poisson sampling** — The CMD histogram is stored as a **density** (counts per deg²). The expected number of objects for a pixel is obtained by multiplying by the pixel area (in deg²) and drawing from a Poisson distribution. Objects are then sampled from the interpolated CMD histogram.

4. **Position sampling** — Positions are drawn uniformly within the pixel and converted to great-circle coordinates `(phi1, phi2)`.

## Output

The output is a tuple `(catalog, meta)`:

- **`catalog`** — DataFrame with per-survey magnitude columns (e.g. `lsst_yr4_g_obs`, `roman_dc2_F158_obs`) and a `source_type` column. Column names follow the `<namespace>_<band>_obs` convention used throughout `streamobs`.
- **`meta`** — dict with the following keys:

| Key | Content |
|-----|---------|
| `nside` | HEALPix resolution used for position sampling |
| `namespaces` | List of survey namespaces in canonical band order, e.g. `['lsst_yr4', 'roman_dc2']` |
| `bands` | List of band names in canonical order, e.g. `['g', 'F158']` |
| `band1` | First (color) band name |
| `band2` | Second (reference/magnitude) band name |
| `color_edges` | Color bin edges from the CMD grid |
| `mag_edges` | Magnitude bin edges from the CMD grid |

**What the light method gives you**

- Spatial distribution correlated with survey depth and dust.
- Magnitude distribution consistent with the survey selection function.

**What it does not give you** (use the [injection method](background_injection.md) for these)

- Per-object magnitude errors or noise.
- More than two photometric bands.

## Advantages and limitations

| | |
|---|---|
| **Fast** | No injection pipeline per pixel — CMD lookup + sampling only. |
| **No truth catalogs needed** | Resource files are provided by the survey developer. |
| **Accounts for depth variation** | Effective maglim is computed per pixel. |
| **Accounts for dust** | Extinction is folded into the effective magnitude limit. |
| **No magnitude error columns** | Errors are not in the output, but their effect is included in the CMD distribution. |
| **Two bands only** | The CMD is 2-D; additional bands require full injection. |
| **Multi-survey** | The two bands can come from different surveys (one per band). |

## Single-survey usage

```python
from streamobs.surveys import Survey
from streamobs.background import Background

survey = Survey.load('lsst', release='yr4')
bg = Background(survey, source_type='both', method='light', bands=('g', 'r'))

catalog, meta = bg.generate(
    phi1_limits=(-20, 20),
    phi2_limits=(-2, 2),
    gc_frame=frame,
    nside=4096,    # HEALPix resolution; auto-capped to maglim map resolution
)
```

## Multi-survey usage

Pass a list of two surveys — one per CMD band.  The canonical band order is determined automatically.

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
    nside=4096,
)
# meta['namespaces'] == ['lsst_yr4', 'roman_dc2']
# catalog columns: lsst_yr4_g_obs, roman_dc2_F158_obs, source_type
```

The resource files must exist for the requested survey combination. See [build_background_resources](build_background_resources.md) for how to build them.
