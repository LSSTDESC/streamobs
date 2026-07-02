# Roman HLWAS Survey Files

This page is the **data sheet** for the three real-footprint Roman High Latitude Wide
Area Survey (HLWAS) releases: `hlwas_wide`, `hlwas_medium`, and `hlwas_all`. All three
reuse the DC2-derived selection-function tables (:doc:`roman_dc2`); the only
tier-specific product is the F158 magnitude-limit map, built from the official HLWAS
exposure-time maps via the exposure-scaled quasi-depth recipe (**Option B**)
documented in :doc:`selection_function_methodology`. Validation figures are embedded
in that methodology page's *Validation & audits* section.

## Survey tiers

The HLWAS comprises four nested tiers; streamobs supports three:

| Release | Footprint | F158 measured map median | Valid pixels (nside=1024) | Area |
|---|---|---|---|---|
| `hlwas_wide`   | Wide tier only | 26.2842 AB | ~1,028,570 | ~3372 deg² |
| `hlwas_medium` | Medium tier only | 26.2894 AB | ~879,052 | ~2882 deg² |
| `hlwas_all`    | All tiers stacked (wide + medium + deep + ultra-deep) | 26.2894 AB | ~1,792,907 | ~5878 deg² |

F158 wide ≈ medium because the medium tier's extra depth relative to wide is in
*other* HLWAS bands, not F158 — both tiers have the same median F158 exposure time
(~645 s) and therefore nearly identical F158 maglim medians under the recipe. The deep
and ultra-deep tiers contribute <0.5% of `hlwas_all` pixels, so the stacked median
stays at the wide/medium level (their deeper pixels appear locally in the map but do
not shift footprint-wide statistics).

## Exposure-time maps

Exposure-time maps (F158, healsparse `.hsp`, nside=4096) are sourced from the
[spacetelescope/roman_notebooks](https://github.com/spacetelescope/roman_notebooks/tree/main/notebooks/footprint_visualization/aux_data)
repository (`map_HLWAS-wide_F158.hsp`, `map_HLWAS-medium_F158.hsp`,
`map_HLWAS-all_F158.hsp`). Values are total on-sky integration time in seconds
(accounting for the planned dither/pass strategy) and are used directly — no unit
conversion or per-tier normalisation. The F158 exposure time is quantised in ~107.5 s
single-exposure units.

## Quasi-depth maps (Option B)

The maglim maps apply the exposure-scaled quasi-depth recipe (derivation and
rationale in :doc:`selection_function_methodology`):

```
depth(pix) = DC2_REF_DEPTH + 1.25 * log10( t(pix) / DC2_REF_EXPTIME )
```

- `DC2_REF_DEPTH` ≈ 26.375 AB — median of the DC2 F158 truth-anchored maglim map
  (`roman_dc2_maglim_f158_nside1024.fits.gz`), read from the file at runtime (not
  hardcoded).
- `DC2_REF_EXPTIME` = 770.0 s — the DC2 HLIS reference per-pixel exposure time
  (5.5 dithers × 140 s; Troxel et al. 2023, Sec. 3.1).

Option B anchors all tiers to the same DC2 truth-anchored reference (not to per-tier
STScI ETC depths of mixed vintage), so the maps land exactly on the DC2-relative
`delta_mag` convention the shared tables need, and inter-tier comparisons are
self-consistent. The HLWAS medians (~26.28–26.29) sit slightly below the DC2 reference
because the typical HLWAS exposure (~645 s) is shorter than the DC2 reference (770 s):
`26.375 + 1.25 × log10(645/770) ≈ 26.28`.

### Measured map medians

| Release | t_median (s) | DC2_REF_DEPTH | DC2_REF_EXPTIME | measured median (AB) |
|---|---|---|---|---|
| `hlwas_wide`   | 645.1 | 26.375 | 770.0 | 26.2842 |
| `hlwas_medium` | 645.1 | 26.375 | 770.0 | 26.2894 |
| `hlwas_all`    | 645.1 | 26.375 | 770.0 | 26.2894 |

Maps are written at nside=1024 (RING, float32) to match the DC2 maps, by
`scripts/roman/build_hlwas_maglim_maps.py`, to `data/surveys/roman_hlwas_<tier>/`
(gitignored). For reference, the
[STScI community-defined HLWAS median 5σ point-source depths](https://roman-docs.stsci.edu/roman-community-defined-surveys/high-latitude-wide-area-survey)
are F158 26.2 (wide) and 26.4 (medium); these differ from F158-only because the
medium tier's reported depth includes additional bands.

## Selection-function tables

The completeness and photo-error tables are **identical** to `roman_dc2` and live
under `data/surveys/roman_dc2/`; the build script symlinks them into each tier's data
directory so the loader finds them at the default path:

| File | Description |
|---|---|
| `roman_stellar_efficiency_cutf158.csv` | F158 detection + classification efficiency vs delta_mag |
| `roman_photoerror_f158_catalog.csv` | Median reported magerr vs delta_mag (S/N cut) |
| `roman_photoerror_f158.csv` | Truth-based scatter of (obs − true) vs delta_mag (noise draw) |

See :doc:`roman_dc2` for these products and :doc:`selection_function_methodology` for
how they are derived.

## Configuration

Each tier is a YAML in `config/surveys/`: `roman_hlwas_wide.yaml`,
`roman_hlwas_medium.yaml`, `roman_hlwas_all.yaml` (releases `hlwas_wide`,
`hlwas_medium`, `hlwas_all`). All use F158 only (F106/F129 HLWAS maps are not yet
built; F184 is excluded per the DC2 documentation). Extinction coefficients and the
saturation threshold are copied from `roman_dc2.yaml` (same instrument). The
`delta_saturation` (= saturation − map median) is keyed to each tier's measured
median:

| Release | map_median (AB) | delta_saturation |
|---|---|---|
| `hlwas_wide`   | 26.2842 | −9.2842 |
| `hlwas_medium` | 26.2894 | −9.2894 |
| `hlwas_all`    | 26.2894 | −9.2894 |

The column namespace per release is `{name}_{release}`: `roman_hlwas_wide`,
`roman_hlwas_medium`, `roman_hlwas_all`.

## Using the surveys in streamobs

```python
from streamobs.surveys import Survey

survey_wide = Survey.load('roman', release='hlwas_wide')
maglim = survey_wide.get_maglim('F158', pixel=123456)
eff    = survey_wide.get_completeness('F158', mag=25.0, maglim=maglim)
err    = survey_wide.get_photo_error('F158', magnitude=25.0, maglim=maglim)
```

## Regenerating the maps

```bash
conda activate streamobs
python scripts/roman/build_hlwas_maglim_maps.py
```

This reads the healsparse exposure-time maps, reads `DC2_REF_DEPTH` from the DC2 F158
maglim map at runtime, applies the Option B recipe, writes the three nside=1024 maglim
maps, and symlinks the DC2 CSV files into each tier's data directory.

## Caveats

- The quasi-depth recipe assumes photon-noise-limited exposures and a uniform PSF;
  real HLWAS depths vary with detector gaps, background, and read noise in the
  shortest exposures.
- The DC2-derived completeness/photo-error tables characterise the deep DC2 reference
  simulation; their use for HLWAS relies on the selection function depending on
  magnitude only through `delta_mag = mag − maglim`.
- The `hlwas_all` map includes deep/ultra-deep pixels (t up to ~23,600 s; local maglim
  up to ~29.5 AB) where the shallow DC2 selection function may be unreliable.
- ETC-derived or truth-anchored depth maps for the real HLWAS are a future upgrade
  that will supersede this quasi-depth approach.
