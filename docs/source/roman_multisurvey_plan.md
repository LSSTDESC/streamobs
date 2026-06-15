# Multi-survey (Roman + Rubin) injection — design & roadmap

```{warning}
This is a **living design document** for in-progress work, not a description of
shipped behaviour. It records the agreed design, the phased rollout, and — most
importantly — the **behaviour changes that need discussion** before we rely on
them. Sections marked *(future)* are designed but not yet implemented.
```

## Motivation

`streamobs` was built for **single-survey** stream injection (LSST). We want to
inject mock streams that carry **both Roman and Rubin/LSST photometric columns**
in one catalog, where each band draws its photometric errors and detection
probability from its **own** survey magnitude-limit (maglim) HEALPix maps.

A working prototype of this idea exists outside the package (the
`rubin_roman_object_classification` proposal code: survey-indexed dictionaries,
`{survey}_{band}` columns, and a "sample masses once, interpolate per survey"
photometry step). The goal here is to bring that capability **into** `streamobs`
properly, reusing the package's existing {class}`~streamobs.surveys.Survey` /
maglim-map machinery (which the prototype lacked), rather than maintaining a
parallel codebase.

Two science realities drive specific features:

- The Roman DC2 characterization (see *Roman DC2 Survey Files*) found the **true
  photometric scatter is ~2× the reported catalog error**. So the error model
  needs a *sample* model (true scatter, used to draw observed magnitudes)
  distinct from a *catalog* model (the reported `magerr`).
- We need a **realistic background** (field stars + galaxies, with star/galaxy
  misclassification) **without** loading the heavy DC2 mock catalogs at runtime
  — using lightweight binned colour–magnitude distributions vs. maglim, and
  optionally TRILEGAL for the stellar component. *(Designed below; future work.)*

## Locked design decisions

1. **Multi-survey via a `MultiSurveyInjector` orchestrator** that holds
   `{survey_name: Survey}` and delegates per-survey work to the existing,
   unchanged `Survey` per-band API. Not a composite `Survey`; not manual
   re-runs.
2. **Column convention.** Multi-survey outputs use `{survey}_{band}_true`,
   `{survey}_{band}_obs`, and `{survey}_{band}_err` (e.g. `roman_f158_obs`,
   `lsst_r_err`). The single-survey `StreamInjector` **keeps emitting the legacy
   names** `mag_{band}_obs` / `magerr_{band}` / `flag_observed` unchanged.
3. **Sample-vs-catalog error split, backward-compatible by default.** `Survey`
   holds two error curves, both functions of `delta_mag = mag − maglim`:
   - `log_photo_error_catalog` — the survey's **reported** error curve (the
     existing `photoerror_*.csv`). Written as `magerr` and drives the S/N cut.
     Always present; loaded from the `log_photo_error_catalog` key or the legacy
     `log_photo_error` key.
   - `log_photo_error_sample` — an **optional second** curve giving the **true
     scatter** of observed−true magnitudes; drives the noise draw. Config key
     `log_photo_error_sample`.

   If no sample curve is supplied, the noise draw falls back to the catalog
   curve ⇒ existing outputs unchanged. Separating them is opt-in: supply the
   second CSV. (Earlier drafts proposed a scalar `report_error_factor` /
   inflation curve to *derive* one from the other; superseded by two independent
   curves, which is more general and matches the Roman DC2 products that measure
   the true scatter directly.)
4. **Sample stellar masses once, interpolate per survey.** This is a correctness
   requirement (see below) and the engine of the whole feature.
5. **Backward compatibility is by column-schema + public API**, not internal
   coupling. Downstream consumers read pre-generated CSVs by column name
   (`mag_g_obs`/`mag_r_obs`/`flag_observed`) and use the
   `StreamInjector(survey).inject(df, bands=[...])` API. Both are preserved;
   multi-survey support is purely additive.

## Behaviour changes for discussion

```{important}
These are deliberate changes from current behaviour. They are flagged here for
review **before** we depend on them — they are not silently adopted.
```

### S/N detection cut now applies to *all* injected bands *(Phase 2, adopted)*

Previously the injector applied its SNR ≥ 5 cut to a single hard-coded band
(`detection_mag_cut=["g"]`). Now the default is **every band passed to
`inject()`** — a star must have SNR ≥ 5 in *all* injected bands to be flagged
observed. For the default LSST `bands=["r", "g"]` this is stricter than before
(it adds the r-band SNR requirement on top of g), so detection counts drop
relative to the old default. Rationale: it generalizes cleanly to any band set
and makes "observed" mean "detected in everything you asked for". Callers can
restore any prior behaviour by passing `detection_mag_cut=[...]` explicitly
(e.g. `["g"]` for the old LSST default). **Adopted, but flagged for review.**

### `nstars` becomes "exactly N stars" (was an emergent IMF count)

Today {meth}`~streamobs.model.IsochroneModel.sample` converts `nstars` into a
total stellar mass and lets `ugali`'s `iso.simulate()` return a *random-length*
IMF realization — so the number of stars returned is stochastic and generally
**≠ `nstars`**.

The multi-survey requirement (the *same physical star* must get consistent Roman
**and** Rubin magnitudes) forces us instead to draw a fixed set of initial
masses once and interpolate each survey's magnitudes from them:

```python
init_mass, mass_pdf, *_ = base_iso.sample(mass_min=mass_min, mass_steps=mass_steps)
sampled_masses = rng.choice(init_mass[sel], size=nstars, p=imf_pdf)   # exactly nstars
# then, per survey:  mag_band = np.interp(sampled_masses, init_mass, mag_band) + dist_mod
```

This returns *exactly* `nstars` stars. It is almost certainly the right
semantics for injection (you control N, and it is shared across surveys), but it
changes what `StreamModel.sample(size)` returns and could shift normalization in
any analysis that relied on the old emergent count. **Open for discussion.**

### Where Roman Vega→AB conversion lives

`ugali` appears to deliver Roman isochrone magnitudes in **Vega**, while our
catalogs are **AB**; the prototype corrects this after `isochrone_factory` with a
`ROMAN_ZEROPOINTS[...]["diff"]` table.

**Current plan (do both, in order):**

1. **This branch:** implement the correction in `streamobs` as a single isolated
   shim — `IsochroneModel._apply_vega_to_ab`, gated by a per-survey `vega_to_ab`
   config flag — so the feature works now without blocking on an upstream
   release. It is structured so that if `ugali` later returns AB natively we just
   set `vega_to_ab: false` and the shim becomes a no-op / is removed.
2. **Off the critical path:** pursue native AB support (or an explicit
   photometric-system flag) in `ugali` itself, since `ugali` owns the isochrone
   photometric system and every downstream user would benefit.

*(If we decide to commit to the `ugali`-native route up front, the shim in
Phase 3 is dropped.)*

## Phased rollout

Phases 1–4 are the current branch; Phases 5–6 are designed here but deferred.

### Phase 1 — Sample vs. catalog error models on `Survey` ✅ *(implemented)*

Two independent error curves, both vs `delta_mag = mag − maglim`:

- ✅ `Survey`: replaced the single `log_photo_error` field with
  `log_photo_error_catalog` (reported error, the base curve) and
  `log_photo_error_sample` (optional true-scatter curve). Kept a read/write
  `log_photo_error` property aliasing the **catalog** model (the legacy field's
  meaning — back-compat for existing tests and any code that sets it).
- ✅ `Survey.get_photo_error(band, mag, maglim, kind="catalog")`: `kind` selects
  the curve via `_resolve_log_photo_error`. `kind="catalog"` returns the reported
  error; `kind="sample"` returns the true-scatter curve, falling back to the
  catalog curve when no sample curve is loaded. Default `kind="catalog"`
  reproduces today's numbers exactly.
- ✅ `SurveyFactory._load_survey_data`: loads the catalog curve from
  `log_photo_error_catalog` (or the legacy `log_photo_error`) key, and an
  optional sample curve from the `log_photo_error_sample` key. No factor /
  inflation logic — both curves are read directly via `set_photo_error`.
- ✅ `StreamInjector.inject`: draws noise with the **sample** error
  (`kind="sample"`); writes the **catalog** error as `magerr` and runs the S/N
  cut on it (`kind="catalog"`).

**Backward compatibility verified (94 passing tests from the test branch):** with
only the legacy `log_photo_error` curve and no `log_photo_error_sample`, the
sample draw falls back to the catalog curve, so the noise draw, `magerr`, and the
S/N cut are bit-for-bit identical to the previous single-curve behaviour. To opt
in, add `log_photo_error_sample: <file>` (the measured true-scatter curve) to a
survey's `survey_files`.

### Phase 2 — De-hardcode the injector to arbitrary bands ✅ *(implemented)*

- ✅ New `streamobs/columns.py` with `true_col` / `obs_col` / `err_col` helpers
  `(band, survey=None)` and `flag_col(survey=None)` (`survey=None` ⇒ legacy
  names `mag_<band>` / `mag_<band>_obs` / `magerr_<band>` / `flag_observed`;
  a survey name ⇒ the `<survey>_<band>_…` multi-survey convention for Phase 4).
- ✅ `observed.py`: removed the `bands in {"r","g"}` hard block; the true-mag
  read, the observed/err columns, the valid-flux check (now ANDs over every
  injected band), the S/N cut, the per-survey detection flag, and the stored
  flag all route through the `columns.py` helpers and the survey's
  `completeness_band`. Existing per-band nside handling already supports Roman
  nside=1024 vs. LSST nside=128.
- ⚠️ The S/N-cut default changed from the hard-coded `["g"]` to **all injected
  bands** — see *Behaviour changes for discussion* above.

**Validated:** single-band (`bands=["r"]`) and arbitrary band sets now inject
without the old hard block; legacy column names and the `inject(df, bands=[...])`
API are unchanged; the test-branch suite stays green (94 passing; the lone
`des_yr6` photo-error failure is pre-existing and unrelated).

### Phase 3 — Multi-band / multi-survey `IsochroneModel` ✅ *(implemented)*

- ✅ `IsochroneModel.create_isochrone` accepts either today's single-survey
  config or a multi-survey form (`surveys: {name: {survey, band_1, band_2,
  vega_to_ab}}` plus shared `name`/`age`/`z`/... at top level), building one
  `ugali` isochrone per survey (`self.isos`, `self.survey_bands`).
- ✅ New `sample_masses(...)` draws the initial masses *once* from the primary
  isochrone's IMF (exactly `nstars`), and `sample_multisurvey(...)` interpolates
  those shared masses into every survey's bands → `{(survey, band):
  apparent_mag}` (same physical star, consistent across surveys). The legacy
  `sample(nstars, dm)` still returns the `(mag_band_1, mag_band_2)` tuple (the
  primary survey's two bands) so existing callers are unchanged.
- ✅ `_apply_vega_to_ab(survey, band, mag)` shim — gated by a per-survey
  `vega_to_ab: {band: offset}` mapping (`AB = Vega + offset`); a no-op when
  absent. Applied in both the single- and multi-survey paths.
- ✅ `StreamModel.sample`/`complete_catalog` derive their magnitude columns from
  the isochrone via `_iso_mag_columns()` / `_sample_iso_mags()`: legacy
  `mag_g`/`mag_r` for a single-survey isochrone, `<survey>_<band>_true` for a
  multi-survey one. Naming routes through `columns.true_col`.

**Note on the chosen API:** the plan originally named the dict-returning method
`sample` and a tuple `sample_legacy`; to avoid `sample()` changing return *type*
by config (a foot-gun for existing callers), the implementation keeps `sample()`
as the `(mag_1, mag_2)` tuple and adds `sample_multisurvey()` for the dict.
`StreamModel` dispatches on `isochrone.multi_survey`.

**Validated:** single-survey model tests unchanged (test branch green, 94
passing); a two-isochrone multi-survey config produces consistent shared-mass
magnitudes, the Vega→AB offset shifts only the configured band, and
`StreamModel.sample` emits the `<survey>_<band>_true` columns.

### Phase 4 — `MultiSurveyInjector`

- Refactor the per-band body of `StreamInjector.inject` into a shared
  `_inject_one_survey(...)` helper.
- `MultiSurveyInjector(surveys).inject(data, survey_bands, ...)`:
  (1) one shared sky placement; (2) one shared true-magnitude fill (masses
  sampled once); (3) a per-survey loop writing `{survey}_{band}_obs/_err` and
  `{survey}_flag_observed` using each survey's own `completeness_band` and maglim
  maps. Per-survey RNG via `rng.spawn(...)` for order-independent reproducibility.
- A *scene* config (`config/scenes/roman_rubin_demo.yaml`) lists the surveys,
  bands, per-survey isochrones, and shared stream geometry.

### Phase 5 *(future)* — Lightweight background + galaxy misclassification

- New `streamobs/background.py`: a `CMDDistribution` (binned colour–magnitude
  distribution vs. maglim, stored raw so one file serves any isochrone) and a
  `BackgroundGenerator` that, per HEALPix pixel, looks up the local maglim,
  selects the nearest-maglim CMD slice, scales counts linearly by pixel area,
  Poisson-draws the count, samples multi-band magnitudes, places objects
  uniformly within the pixel, and optionally applies a matched filter.
- Independent, pluggable `StellarBackground` / `GalaxyBackground` models.
- Galaxy misclassification: a new `Survey.get_galaxy_misclassification` curve and
  an `is_galaxy`-aware `detect_flag`, so misclassified galaxies leak into the
  stellar sample (and `perfect_galstarsep=True` ⇒ no leakage). A build script
  reads DC2 truth+det **once, offline** to emit the lightweight files; the
  runtime never loads DC2.

### Phase 6 *(future)* — TRILEGAL + docs

- `TrilegalStellarBackground` implementing the same interface, **lazy** (reads a
  user-provided TRILEGAL table or a `fetch` callable; never imported at module
  import). Documentation polish for the `{survey}_{band}` outputs, the error
  split, and scenes.

## Development fixtures (this branch)

Survey data files (maglim maps, efficiency / photo-error CSVs) are **not**
committed — they are downloaded into the git-ignored `data/surveys/<survey>/`.
To keep this branch self-contained and testable without the real Roman/Rubin
data, we add **dummy surveys**: small committed configs
(`config/surveys/*_dummy.yaml`) plus a generator that synthesizes tiny HEALPix
maglim maps and CSV tables at test time. The real `roman_dc2.yaml` is **not**
recreated here, so it can land cleanly when the `roman_hlwas` work merges.

## Key files

| File | Phase | Role |
|---|---|---|
| `streamobs/surveys.py` | 1 | sample/catalog error fields, `get_photo_error(kind=)`, loader |
| `streamobs/columns.py` | 2 | NEW — column-name helpers |
| `streamobs/observed.py` | 1,2,4 | error wiring; de-hardcode bands; `_inject_one_survey`; `MultiSurveyInjector` |
| `streamobs/model.py` | 3 | multi-band `IsochroneModel`; band-generalized `StreamModel` |
| `streamobs/multisurvey.py` | 4 | NEW — orchestrator (or a class in `observed.py`) |
| `config/scenes/roman_rubin_demo.yaml` | 4 | NEW — multi-survey scene |
| `streamobs/background.py` | 5 *(future)* | NEW — lightweight background |
