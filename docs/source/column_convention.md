# Output column convention

Injected catalogs use a single, uniform naming scheme for every column the
injector produces, and it is **always survey-namespaced** — even for a single
survey.

## The scheme

For a survey with namespace `<survey>` (the survey's name, or the key you gave it
when constructing a multi-survey `StreamInjector`) and a photometric band
`<band>`:

| Column | Meaning |
|---|---|
| `<survey>_<band>_true` | True (noiseless) apparent magnitude |
| `<survey>_<band>_obs` | Observed (noisy) magnitude; `"BAD_MAG"` where the noisy flux went negative |
| `<survey>_<band>_err` | Reported photometric error (the *catalog* error; see below) |
| `<survey>_flag_observed` | Detection flag — `True` if detected **and** classified as a star |
| `<survey>_flag_perfect_galstarsep` | Optional flag assuming perfect star/galaxy separation (detection only); emitted only when `perfect_galstarsep=True` |

Plus the shared, un-namespaced sky coordinates `ra`, `dec`.

Examples: `lsst_r_obs`, `lsst_g_err`, `roman_F158_true`, `lsst_flag_observed`.

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

## Two error curves: catalog vs. sample

Each survey carries **two** photometric-error curves, both functions of
`delta_mag = mag − maglim`:

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

Select between them via
{meth}`streamobs.surveys.Survey.get_photo_error` with `kind="catalog"` (default,
reproduces the reported error) or `kind="sample"` (true scatter).

## See also

- [Injecting one or many surveys](multisurvey.md) — how these columns are
  produced and the multi-survey "same physical star" guarantee.
- `streamobs.columns` — the column-name helper functions.
- {meth}`streamobs.surveys.Survey.get_photo_error` — the two error curves.
