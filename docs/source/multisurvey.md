# Injecting one or many surveys

`streamobs.observed.StreamInjector` injects observational effects (photometric
errors, observed magnitudes, detection flags) into a stream catalog. **The same
class handles a single survey or several at once** — there is no separate
multi-survey class. When several surveys are injected together, every survey
describes the *same physical stars*, so you can build catalogs that carry both
Rubin/LSST and Roman photometry for one stream.

## Constructing an injector

Pass one survey, or several:

```python
from streamobs.observed import StreamInjector

# One survey (loaded by name; `release` etc. forwarded to Survey.load).
# Its namespace is "{name}_{release}", here "lsst_dc2".
inj = StreamInjector("lsst", release="dc2")

# Several surveys as a list of specs — each spec is a survey name, a Survey, or a
# {"survey": ..., "release": ...} dict. The namespace is derived from each loaded
# Survey ("lsst_dc2", "roman_dc2"), NOT from any key you supply.
inj = StreamInjector([
    {"survey": "lsst", "release": "dc2"},
    {"survey": "roman", "release": "dc2"},
])

# A list of plain names loads each release-less, so the namespace is the bare
# name ("lsst", "roman").
inj = StreamInjector(["lsst", "roman"])
```

- The namespace (the column prefix) is always the survey's own
  `{name}_{release}` ({attr}`streamobs.surveys.Survey.namespace`). A `{key: spec}`
  dict is also accepted, but the keys are containers only — the namespace is
  re-derived from each `Survey`, not taken from the key.
- `primary` selects the survey whose footprint drives the shared sky placement
  (defaults to the first survey).
- `inj.surveys` is the `{namespace: Survey}` mapping; `inj.primary` (alias
  `inj.survey`) is the primary `Survey`, and `inj.primary_namespace` its namespace
  string.

## Injecting

```python
# Single survey: `bands` is the shorthand (defaults to ['r', 'g']).
out = inj.inject(df, bands=["r", "g"], stream_config=cfg, seed=42)

# Several surveys: give the bands per survey as a {namespace: [bands]} dict,
# keyed by each survey's namespace ({name}_{release}).
out = inj.inject(
    df,
    bands={"lsst_dc2": ["r", "g"], "roman_dc2": ["F106", "F158"]},
    stream_config=cfg,
    seed=42,
)
```

A plain list is rejected for a multi-survey injector (it is ambiguous), and a
`bands` dict referencing an unknown namespace raises `ValueError`.

`df` may already contain `ra`/`dec` or `phi1`/`phi2`, may be a fully empty frame
of length *N*, or any subset — anything missing is sampled from `stream_config`
(see *Completing a catalog* below). The output carries shared
`ra`/`dec` plus, **per survey**, the namespaced columns described in
[Output column convention](column_convention.md):

```
ra, dec,
lsst_dc2_r_true,  lsst_dc2_r_obs,  lsst_dc2_r_err,  lsst_dc2_g_true, ..., lsst_dc2_flag_observed,
roman_dc2_F106_true, roman_dc2_F106_obs, ..., roman_dc2_flag_observed
```

Useful `inject` keyword arguments:

| kwarg | meaning |
|---|---|
| `seed` | reproducibility (per-survey RNGs are spawned from it, so results are independent of survey order) |
| `dist` | distance modulus used directly (scalar or per-row vector) instead of sampling one — see below |
| `detection_mag_cut` | non-reference bands to apply the explicit SNR ≥ 5 cut to (see *S/N cut ownership* below) |
| `perfect_galstarsep` | also emit a `<survey>_flag_perfect_galstarsep` flag (detection only, no classification losses) |
| `dust_correction` | apply extinction correction to observed magnitudes (default `True`) |
| `mask_type`, `gc_frame` | forwarded to the `phi1`/`phi2` → `ra`/`dec` placement |

## The same physical star across surveys

For a multi-survey injection the isochrone draws **one set of initial masses**
(exactly `nstars`) and interpolates *those same masses* into every survey's
bands. So a star's LSST and Roman magnitudes describe the same object — the
true magnitudes are physically consistent and tightly correlated across surveys
rather than drawn independently.

This requires a **multi-survey isochrone** in the stream config: a top-level
`surveys:` mapping sharing one stellar population, e.g.

```yaml
stream:
  # ... density / track / distance_modulus ...
  isochrone:
    name: Marigo2017      # shared population
    age: 12.0
    z: 0.0006
    surveys:
      lsst_dc2:  {survey: lsst,  band_1: g,    band_2: r}
      roman_dc2: {survey: roman, band_1: F106, band_2: F158}
```

```{important}
Each `surveys:` **key is the column namespace** the isochrone produces
(`<key>_<band>_true`), and it must match the injecting survey's namespace
`{name}_{release}` — otherwise the true-magnitude columns the model emits won't
line up with the columns the injector looks for. Here the inner `survey:` is the
*ugali* filter set (no release), while the key carries the release. (In the
single-survey flat form the namespace is derived as `{survey}_{release}` for you;
in the multi-survey form you spell it out as the key.)
```

A single-survey isochrone (the flat `survey`/`band_1`/`band_2` form, optionally
with `release:`) is just the one-survey case of the same machinery and produces
`<namespace>_<band>_true` identically.

A complete, runnable example — the surveys, per-survey bands, the multi-survey
isochrone, and the shared stream geometry — is provided as a *scene* config in
[`config/scenes/roman_rubin_demo.yaml`](https://github.com/LSSTDESC/streamobs/blob/main/config/scenes/roman_rubin_demo.yaml):

```python
import yaml
from streamobs.observed import StreamInjector

scene = yaml.safe_load(open("config/scenes/roman_rubin_demo.yaml"))
inj = StreamInjector(scene["surveys"])              # {"lsst": "lsst", "roman": "roman"}
cat = inj.inject(
    df, bands=scene["survey_bands"],                # {"lsst": [...], "roman": [...]}
    stream_config=scene["stream"], seed=42,
)
```

```{note}
**Roman bands are converted Vega→AB automatically.** `ugali` returns Roman
isochrone magnitudes in Vega while the catalogs are AB, so `IsochroneModel`
applies a fixed per-band offset (`streamobs.model.ROMAN_VEGA_TO_AB`) to every
Roman band unconditionally. Non-Roman bands pass through unchanged; there is no
config flag.
```

```{note}
**`nstars` means exactly N stars.** `StreamModel.sample(size)` / the isochrone
draw return *exactly* that many stars (a fixed mass set), not a random-length IMF
realization. This is required so the same masses can be shared across surveys.
```

## Completing a catalog

`StreamInjector.complete_data(...)` is the public "fill in the rest from the
config" helper — the same completion `inject` runs internally, exposed so you can
build or inspect a completed catalog **without** injecting noise. It fills
`ra`/`dec` (converting from `phi1`/`phi2` if needed) and the per-survey
`<survey>_<band>_true` columns, **preserving anything already present**:

```python
# Partial input -> filled from the config; existing columns are kept.
full = inj.complete_data(df, bands=["r", "g"], stream_config=cfg, seed=1)
```

### Supplying a distance directly

Apparent magnitudes need a distance modulus. Normally it comes from the config's
`distance_modulus` model (which needs `phi1`), but you can pass `dist` directly —
a scalar (broadcast to all rows) or a per-row vector — to fill magnitudes without
a distance model or `phi1`:

```python
out = inj.inject(df, bands=["r", "g"], stream_config=cfg, dist=16.8)     # scalar
out = inj.inject(df, bands=["r", "g"], stream_config=cfg, dist=dist_arr)  # per-row
```

When given, `dist` overrides the configured distance model. Only rows that are
missing a `dist` value are set.

### Existing values are never overwritten

Completion fills only the **missing** rows of each column. If you supply one band
and request another, the supplied band is left untouched and only the missing one
is filled (newly-filled cells still come from one shared mass draw, so they are
mutually colour-consistent).

### Stellar masses (the `mass` column)

When an isochrone is configured, the shared **initial masses** drawn for the stars
are surfaced as a single un-namespaced `mass` column (one mass per star, shared by
all surveys — the same physical star). You can also go the other way and supply
your own masses: pass a fully-populated `mass` column in the input catalog and the
isochrone uses *those* masses instead of drawing fresh ones, so the sampled
magnitudes reproduce your simulation's exact stars. At the model level
{meth}`streamobs.model.IsochroneModel.sample` accepts a `masses=`
array and returns the masses it used. The mass grid resolution is controlled by
`IsochroneModel._MASS_STEPS` (default 4000) and a per-call `mass_steps=` override.

## S/N cut ownership

The reference band (`survey.completeness_band`, e.g. LSST `r`) is special: the
survey's **selection-function curves are estimated on true stars detected at
SNR ≥ 5** in that band, so the cut is already baked into both the completeness
and detection-efficiency curves. The injector therefore does **not** re-apply a
SNR cut to the reference band (doing so would double-count it); the explicit
`detection_mag_cut` loop applies SNR ≥ 5 only to the *other* injected bands (its
default is all injected bands except the reference band). Net effect: a star must
have SNR ≥ 5 in every injected band to be flagged observed, with the reference
band's cut owned entirely by the selection-function curves.

## See also

- [Output column convention](column_convention.md) — the `<survey>_<band>_…`
  scheme and the sample-vs-catalog error split.
- {class}`streamobs.observed.StreamInjector` — full API.
- {class}`streamobs.model.StreamModel` — the stream/isochrone config.
