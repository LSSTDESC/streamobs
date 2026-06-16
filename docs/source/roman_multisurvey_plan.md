# Multi-survey (Roman + Rubin) injection — design & roadmap

```{warning}
This is a **living design document** for in-progress work, not a description of
shipped behaviour. It records the agreed design, the phased rollout, and — most
importantly — the **behaviour changes that need discussion** before we rely on
them. Sections marked *(future)* are designed but not yet implemented.
```

```{important}
**Remove this file before merging the `roman_multisurvey` branch.** It is a
working design/roadmap doc, not user documentation. Before merge, migrate the
durable content — the column convention, the sample/catalog error split, the
multi-survey `StreamInjector` usage, and the Vega→AB handling — into the proper
docs pages (and the API docstrings), then delete `roman_multisurvey_plan.md`.
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

1. **One `StreamInjector` for one *or* many surveys.** A single injector holds
   `{survey_name: Survey}` (size 1 for the single-survey case) and delegates
   per-survey work to the existing, unchanged `Survey` per-band API via a shared
   `_inject_one_survey` helper. Not a composite `Survey`; not manual re-runs;
   and **no separate multi-survey class** — `StreamInjector` accepts a survey
   name/`Survey`, a `{name: spec}` dict, or a list of specs.
2. **Column convention (always survey-namespaced `{survey}_{band}_true/_obs/_err`).**
   A single naming scheme is used everywhere — `{survey}_{band}_true` /
   `{survey}_{band}_obs` / `{survey}_{band}_err` / `{survey}_flag_observed`
   (e.g. `roman_f158_obs`, `lsst_r_err`) — and it is **always namespaced by
   survey, even for a single survey** (a one-survey injection of LSST emits
   `lsst_r_obs`, not `r_obs`). There is no longer an un-namespaced single-survey
   form. **This intentionally drops the historical `mag_{band}` /
   `mag_{band}_obs` / `magerr_{band}` names — it is *not* backward compatible**
   with catalogs or downstream readers expecting those columns (decision made
   deliberately to keep one convention; see decision 5).
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
5. **Public API preserved; column schema deliberately changed.** The
   `StreamInjector(survey).inject(df, bands=[...])` single-survey call still
   works (the same `StreamInjector` also takes a dict/list for several surveys).
   The **column names are not** preserved — we adopt the always-namespaced
   `{survey}_{band}_true/_obs/_err` scheme (decision 2) *instead of* the
   historical `mag_g_obs`/`mag_r_obs`/`magerr_g` names. Downstream consumers that
   read those columns must be updated. This is an accepted break in exchange for
   one convention across single- and multi-survey output.

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

### `nstars` becomes "exactly N stars" (was an emergent IMF count) *(adopted — agreed)*

```{important}
**Adopted (Phase 4) and agreed.** {meth}`~streamobs.model.IsochroneModel.sample`
now draws *exactly* `nstars` (the shared-mass path described below), for **both**
single-survey and multi-survey isochrones. The single-survey
`ugali.simulate()` emergent-count path has been removed. This was reviewed and
**accepted as the intended semantics** (you control N, and it is shared across
surveys); it is no longer open for discussion.
```

Previously {meth}`~streamobs.model.IsochroneModel.sample` converted `nstars`
into a total stellar mass and let `ugali`'s `iso.simulate()` return a
*random-length* IMF realization — so the number of stars returned was stochastic
and generally **≠ `nstars`**.

The multi-survey requirement (the *same physical star* must get consistent Roman
**and** Rubin magnitudes) forces us instead to draw a fixed set of initial
masses once and interpolate each survey's magnitudes from them:

```python
init_mass, mass_pdf, *_ = base_iso.sample(mass_min=mass_min, mass_steps=mass_steps)
sampled_masses = rng.choice(init_mass[sel], size=nstars, p=imf_pdf)   # exactly nstars
# then, per survey:  mag_band = np.interp(sampled_masses, init_mass, mag_band) + dist_mod
```

This returns *exactly* `nstars` stars. It is the right semantics for injection
(you control N, and it is shared across surveys). It changes what
`StreamModel.sample(size)` returns relative to the old emergent count, but this
was reviewed and **agreed** — no further discussion needed.

### Roman Vega→AB conversion is automatic and unconditional

`ugali` delivers Roman isochrone magnitudes in **Vega**, while our catalogs are
**AB**. `IsochroneModel` corrects this **unconditionally** for every Roman band
using the module-level table `streamobs.model.ROMAN_VEGA_TO_AB` (AB = Vega +
offset). The per-band offsets are the mode of the by-chip Roman zeropoints
(`Roman_zeropoints_20240301.ecsv`), the same values the
`rubin_roman_object_classification` prototype used:

| band | F062 | F087 | F106 | F129 | F146 | F158 | F184 | F213 |
|---|---|---|---|---|---|---|---|---|
| AB − Vega | 0.153 | 0.481 | 0.660 | 1.051 | 1.164 | 1.315 | 1.556 | 1.837 |

There is **no config flag** — Roman bands are always converted, non-Roman bands
pass through unchanged. The code carries a `TODO` flagging that this conversion
ideally belongs in `ugali` itself (so isochrones are returned natively in AB);
when that lands the table can be removed.

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
  `(band, survey=None)` and `flag_col(survey=None)`. Injected catalogs are
  **always** survey-namespaced (`<survey>_<band>_…` / `<survey>_flag_observed`);
  `survey=None` is retained only as a low-level fallback that the injector itself
  never uses.
- ✅ `observed.py`: removed the `bands in {"r","g"}` hard block; the true-mag
  read, the observed/err columns, the valid-flux check (now ANDs over every
  injected band), the S/N cut, the per-survey detection flag, and the stored
  flag all route through the `columns.py` helpers and the survey's
  `completeness_band`. Existing per-band nside handling already supports Roman
  nside=1024 vs. LSST nside=128.
- ⚠️ The S/N-cut default changed from the hard-coded `["g"]` to **all injected
  bands** — see *Behaviour changes for discussion* above.

**Validated:** single-band (`bands=["r"]`) and arbitrary band sets inject without
the old hard block; the `inject(df, bands=[...])` API is unchanged (output now
namespaced); `tests/test_observed.py` + `tests/test_model.py` green.

### Phase 3 — Multi-band / multi-survey `IsochroneModel` ✅ *(implemented)*

- ✅ `IsochroneModel.create_isochrone` accepts either today's single-survey
  config or a multi-survey form (`surveys: {name: {survey, band_1, band_2}}`
  plus shared `name`/`age`/`z`/... at top level), building one `ugali` isochrone
  per survey (`self.isos`, `self.survey_bands`).
- ✅ New `sample_masses(...)` draws the initial masses *once* from the primary
  isochrone's IMF (exactly `nstars`), and `sample_multisurvey(...)` interpolates
  those shared masses into every survey's bands → `{(survey, band):
  apparent_mag}` (same physical star, consistent across surveys). The legacy
  `sample(nstars, dm)` still returns the `(mag_band_1, mag_band_2)` tuple (the
  primary survey's two bands) so existing callers are unchanged.
- ✅ `_to_ab(band, mag)` converts Roman bands Vega→AB **unconditionally** using
  the `ROMAN_VEGA_TO_AB` table (no config flag; non-Roman bands pass through).
  Applied in the shared `sample_multisurvey` path. See the Vega→AB section above.
- ✅ `StreamModel.sample`/`complete_catalog` derive their magnitude columns from
  the isochrone via `_iso_mag_columns()` / `_sample_iso_mags()`, which **always**
  emit `<survey>_<band>_true` (a single-survey isochrone simply has one survey;
  `IsochroneModel` tracks `surveys`/`survey_bands` in both config forms). Naming
  routes through `columns.true_col`.

**Note on the chosen API:** the plan originally named the dict-returning method
`sample` and a tuple `sample_legacy`; to avoid `sample()` changing return *type*
by config (a foot-gun for existing callers), the implementation keeps
`IsochroneModel.sample()` as the `(mag_1, mag_2)` tuple and adds
`sample_multisurvey()` for the dict. `StreamModel` always goes through
`sample_multisurvey()` (a single-survey isochrone is just the one-survey case),
so the emitted columns are uniformly `<survey>_<band>_true`.

**Validated:** model tests green; a two-isochrone multi-survey config produces
consistent shared-mass magnitudes, Roman bands are converted Vega→AB by the
fixed `ROMAN_VEGA_TO_AB` offsets, and `StreamModel.sample` emits the
`<survey>_<band>_true` columns.

### Phase 4 — one `StreamInjector` for one *or* many surveys ✅ *(implemented)*

- ✅ The per-band body of injection lives in a shared
  `StreamInjector._inject_one_survey(data, bands, survey, survey_namespace, ...)`
  helper. It assumes positions and true magnitudes are already present and only
  does the observed/err draw, detection flags, and S/N cut, routing every column
  name through `columns.py`. `survey`/`survey_namespace` are passed explicitly so
  the same method serves every survey.
- ✅ `StreamInjector` accepts **one survey or several** (a name/`Survey`, a
  `{namespace: spec}` dict, or a list); `__init__` normalizes to
  `self.surveys = {namespace: Survey}` with a `survey` property pointing at the
  `primary`. `inject(data, survey_bands=None, bands=None, stream_config=...)`:
  (1) one shared sky placement; (2) one shared true-magnitude fill via the
  isochrone (masses sampled once → every survey's `<survey>_<band>_true`);
  (3) a per-survey loop calling `_inject_one_survey` with `survey_namespace=<name>`,
  writing `<survey>_<band>_obs/_err` and `<survey>_flag_observed` from that
  survey's own `completeness_band` and maglim maps. Per-survey RNG via
  `rng.spawn(...)` for order-independent reproducibility. `columns.perfect_flag_col`
  namespaces the optional `perfect_galstarsep` flag. The separate
  `MultiSurveyInjector` class has been **removed**.
- ✅ A *scene* config (`config/scenes/roman_rubin_demo.yaml`) lists the surveys,
  per-survey bands, the multi-survey isochrone, and shared stream geometry.
- ✅ `notebooks/multisurvey_phases_demo.ipynb` walks through Phases 1–4 end to
  end (executes against `StubSurvey` + real `ugali` isochrones).

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
| `streamobs/columns.py` | 2 | NEW — column-name helpers (always namespaced) |
| `streamobs/observed.py` | 1,2,4 | error wiring; de-hardcode bands; unified one-or-many-survey `StreamInjector` (`_inject_one_survey` + `_complete_shared`) |
| `streamobs/model.py` | 3 | multi-band `IsochroneModel`; always-namespaced `StreamModel` |
| `config/scenes/roman_rubin_demo.yaml` | 4 | NEW — multi-survey scene |
| `streamobs/background.py` | 5 *(future)* | NEW — lightweight background |
