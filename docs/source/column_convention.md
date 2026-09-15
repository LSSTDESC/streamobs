# Output column convention

Injected catalogs use a single, uniform naming scheme for every column the
injector produces, and it is **always survey-namespaced** — even for a single
survey.

## The scheme

Every column is prefixed with the survey's **namespace**, which is
`{name}_{release}` (e.g. `lsst_yr5`, `roman_dc2`), or just `{name}` when the
survey was loaded without a release. The namespace is always derived from the
`Survey` itself ({attr}`streamobs.surveys.Survey.namespace`). **Observed**,
**error**, and **flag** columns carry the full namespace (release included), so
the same survey at two releases produces distinct, non-colliding observations.
**True-magnitude columns are the exception**: because a star's intrinsic
(noiseless) magnitude does not depend on which survey/release observed it, they
are keyed on the survey **name only** — the release is dropped (e.g.
`roman_F158_true`, shared across `roman_dc2` and `roman_hlwas_*`). When you
construct a multi-survey injector from
a `{key: spec}` dict the keys are containers only — they do **not** become the
namespace (it is re-derived from each loaded `Survey`).

For a survey with namespace `<survey>` (= `{name}_{release}`), survey name
`<name>`, and a photometric band `<band>`:

| Column | Meaning |
|---|---|
| `<name>_<band>_true` | True (noiseless) apparent magnitude — keyed on the survey **name only** (release-independent) |
| `<survey>_<band>_obs` | Observed (noisy) magnitude; `"BAD_MAG"` where the noisy flux went negative |
| `<survey>_<band>_err` | Reported photometric error (the *catalog* error; see below) |
| `<survey>_flag_observed` | Detection flag — `True` if detected **and** classified as a star |
| `<survey>_flag_perfect_galstarsep` | Optional flag assuming perfect star/galaxy separation (detection only); emitted only when `perfect_galstarsep=True` |

Plus the shared, un-namespaced sky coordinates `ra`, `dec`.

Examples (LSST loaded with `release="yr5"`, Roman with `release="dc2"`):
`lsst_yr5_r_obs`, `lsst_yr5_g_err`, `roman_F158_true` (true mag — name only), `lsst_yr5_flag_observed`.

These names are produced by the helpers in `streamobs.columns`
(`true_col`, `obs_col`, `err_col`, `flag_col`, `perfect_flag_col`), which take a
`band` and a `survey` namespace.

```{important}
This convention intentionally **drops** the historical `mag_<band>` /
`mag_<band>_obs` / `magerr_<band>` names and the un-namespaced single-survey form
(`r_obs`, `flag_observed`, …). It is **not backward compatible** with catalogs or
downstream readers expecting those columns — everything is now namespaced by
survey, even when only one survey is injected.
```

## Four error curves: catalog vs. sample, cut vs. no-cut

Each survey carries **four** photometric-error curves, all functions of
`delta_mag = mag − maglim`. They split along two independent axes: *which
quantity* the curve describes (catalog vs. sample), and *which population* it
was measured on (the detected sample, vs. the same sample without the
reference-band S/N cut).

The first axis — catalog vs. sample:

- **Catalog error** — the survey's *reported* error (e.g. SExtractor `magerr`).
  This is what is written to `<survey>_<band>_err`, and it drives the S/N
  detection cut.
- **Sample error** — an optional curve giving the *true scatter* of
  observed − true magnitudes. This is what is used to **draw** the observed
  magnitude (`<survey>_<band>_obs`).

The split exists because, for the Roman DC2 catalogs, the true photometric
scatter is ≈ 2× the reported error — so the noise you draw and the error you
report are genuinely different. When a survey has no sample curve loaded, the
sample draw transparently falls back to the catalog curve, so the two are
identical and outputs match the single-curve behaviour.

The second axis — cut vs. `_nocut`:

- **The cut pair** (`*_photoerror_<band>.csv`, `*_catalog.csv`) is measured on
  the **detected** population and applies to the **reference band** — the band
  whose own detection the measurement was conditioned on
  (`completeness_band` in the survey config).
- **The `_nocut` pair** (`*_nocut.csv`, `*_catalog_nocut.csv`) is measured
  **without** the reference-band S/N cut and applies to **every other band**.
  Photometry in those bands is *forced* at the reference band's position, so it
  is not conditioned on its own detection and must not inherit that
  conditioning.

The two pairs agree brightward of the depth and diverge only faintward, where
the S/N cut truncates the detected sample: its measured scatter turns over and
falls, while the `_nocut` curve keeps rising. Using the cut curve for forced
photometry therefore *understates* faint-band errors — by up to 0.58 dex for
`des/yr6` and 0.38 dex for `delve/dr3_gold`.

`Survey.get_photo_error(band=...)` selects the right pair for you. A release
that ships no `_nocut` curves **raises** for any non-reference band rather than
silently applying the detected-population curve.

Select between them via
{meth}`streamobs.surveys.Survey.get_photo_error` with `kind="catalog"` (default,
reproduces the reported error) or `kind="sample"` (true scatter).

## See also

- [Injecting one or many surveys](multisurvey.md) — how these columns are
  produced and the multi-survey "same physical star" guarantee.
- `streamobs.columns` — the column-name helper functions.
- {meth}`streamobs.surveys.Survey.get_photo_error` — the two error curves.
