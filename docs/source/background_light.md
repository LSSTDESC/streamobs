# Light method

The light method generates background catalogs directly from precomputed color–magnitude diagram (CMD) grids without re-running the injection pipeline. It is the default method and is recommended for most use cases.

## How it works

For each HEALPix pixel in the requested sky region:

1. **Effective magnitude limit** — The observed magnitude limit combined with
   dust absorption to obtain the effective magnitude limit:

   $$m_{\mathrm{eff}} = m_{\mathrm{lim}}(\mathrm{pixel}) - A_{\mathrm{band}} \cdot E(B-V)(\mathrm{pixel})$$

   This maps the real per-pixel depth to the nearest pre-built grid point.

2. **Bilinear CMD interpolation** — The CMD histogram grid is stored at a discrete 2-D lattice of `(maglim_b2, maglim_b1)` pairs (reference band × color band). The generator bilinearly interpolates the four surrounding grid points to obtain a CMD at the exact effective magnitude limit of the pixel.

3. **Poisson sampling** — The expected number of objects is scaled from the reference area to the pixel area via a Poisson draw. Objects are then sampled from the interpolated CMD histogram.

4. **Position sampling** — Positions are drawn uniformly within the pixel and converted to great-circle coordinates `(phi1, phi2)`.

## Output

The output is a tuple `(catalog, meta)`:

- **`catalog`** — DataFrame with columns `phi1`, `phi2`, `mag_<band_ref>`, `mag_<band_color>`, `source_type`.
- **`meta`** — dict with `nside`, `color_edges`, `mag_edges`, `band1`, `band2`.

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
| **No magnitude errors columns** | Errors are not in the output, but their effect is included in the CMD distribution |
| **Two bands only** | The CMD is 2-D; additional bands require full injection. |

## Usage

```python
from streamobs.background import Background

bg = Background(survey, source_type='both', method='light', bands=('g', 'r'))
catalog, meta = bg.generate(
    phi1_limits=(-20, 20),
    phi2_limits=(-2, 2),
    gc_frame=frame,
    nside=4096,    # HEALPix resolution; auto-capped to maglim map resolution
)
```

The resource files must exist for the requested survey. See [doc](build_background_resources.md) for how to build them.
