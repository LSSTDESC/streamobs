# streamobs survey data products

Runtime data for `streamobs.surveys.Survey.load()`: per-release depth maps and
selection-function tables. These files are **not tracked by git** — they are
distributed as a single `data.zip` and fetched with:

```bash
python bin/download_data.py
```

Everything here is a *runtime product*. The large derivation inputs (Balrog
injection catalogs, DC2 truth skims, HealSparse depth-map sources, deep-field
parquet) are deliberately excluded — see `bin/build_data_archive.py` for the
exclusion rules.

## The product contract

Every release is keyed to

```
delta_mag = mag_true - maglim(pixel)
```

where `maglim` is the **5σ** point-source depth at the source's HEALPix pixel and
`mag_true` is the **reddened** true apparent magnitude. A release supplies:

| File | Columns | Drives |
|---|---|---|
| `*_maglim_<band>_nside<N>.fits.gz` | HEALPix RING map | local depth per pixel |
| `*_stellar_efficiency_cut<band>.csv` | `mag_<band>, delta_mag, detection_eff, classification_eff, classification_detection_eff` | detection + star-classification probability |
| `*_photoerror_<band>.csv` (**sample**) | `delta_mag, log_mag_err` | the magnitude **noise draw** (truth scatter of obs − true) |
| `*_photoerror_<band>_catalog.csv` (**catalog**) | `delta_mag, log_mag_err` | the reported `magerr`, which drives the S/N cut |
| `*_photoerror_<band>_nocut.csv`, `*_catalog_nocut.csv` | `delta_mag, log_mag_err` | the same two curves for **forced-photometry** bands |
| `*_galaxy_misclass_cut<band>.csv` | `mag_<band>, delta_mag, missclassification_eff` | stellar contamination from misclassified galaxies |
| `*_audit.json` | — | provenance: row counts, anchors, inflation factor |

Conventions that matter if you are reading these tables directly:

- **The S/N > 5 cut is already baked into the efficiency curves.** Do not
  re-apply it.
- `detection_eff = det/all`, `classification_eff = cls/det`,
  `classification_detection_eff = cls/all`, with the denominator being all
  injected true stars. Bins with fewer than 20 stars are dropped.
- `detection_eff` and `classification_detection_eff` are clamped to zero for
  `delta_mag > 1`.
- The two photo-error curves are **not** interchangeable. The sample curve is
  the truth scatter of (observed − true) and is what you add as noise; the
  catalog curve is the median reported error and is what a survey would print.
  The first exceeds the second by the release's error-inflation factor.
- **Only the reference band's photometry is conditioned on detection.** Every
  other band is measured at the reference band's position, so it is *forced*
  photometry and must use the `_nocut` curves, which are measured without the
  reference-band S/N cut. `Survey.get_photo_error(band=...)` selects the right
  pair; it raises rather than silently applying the detected-population curve.
  The two pairs agree brightward of the depth and diverge faintward, where the
  S/N cut truncates the detected sample and its measured scatter turns over
  instead of continuing to rise.
- The misclassification column is spelled `missclassification_eff` (two s's).
  That is the key `SurveyFactory` reads; a differently-spelled column loads
  silently as `None`.

## Releases in this archive

### DECam — derived from Balrog synthetic source injection

Both are truth-anchored against their own injections: the supplied depth map
sets the spatial structure, and its absolute scale is shifted so the median
equals the magnitude at which the truth scatter of (obs − true) reaches
σ = 2.5/ln10/5 = 0.2171. A map delivered on any S/N convention is therefore
corrected automatically rather than by assuming a +0.753 mag shift.

| | `des/yr6` | `delve/dr3_gold` |
|---|---|---|
| bands | griz | griz |
| injections | DES Y6 Balrog (Anbajagane et al. 2025) | DELVE `BalrogOfTheStars` V4 |
| star classifier | `EXT_XGB` via surrogate + deconvolution | `bdf_extended_class_dr3gold` |
| reference band | g | g |

In the DES footprint DELVE DR3 Gold *is* DES Y6, so the two are expected to
agree in `delta_mag` space; that agreement is the primary cross-check and is
reported in `docs/source/product_verification.md`.

Detail and derivation: `docs/source/balrog_selection_functions.md`,
`docs/source/surveys/DES.md`, `docs/source/surveys/DELVE.md`.

### LSST

- `lsst_dc2` — truth-anchored g/r depth + two-curve photo-error from the DC2
  object and truth skims.
- `lsst_yr1` … `lsst_yr5` — RubinSim per-year depth maps carrying the `lsst_dc2`
  selection-function tables.
- `lsst_dp2` — DP2's own measured 5σ deep-coadd g/r depth (nside 512, not
  truth-anchored) carrying the same `lsst_dc2` tables. Real commissioning data
  with a simulation-derived selection function; see `docs/source/surveys/LSST.md`
  for what that assumes.

### Roman

- `roman_dc2` — reference HLIS depth mock, F106/F129/F158.
- `roman_hlwas_wide`, `_medium`, `_all` — HLWAS tiers.

## Shared files

- `others/ebv_sfd98_lowres_nside_512_ring_equatorial.fits` — SFD98 extinction
  map used by every release.
- `others/lsst_*` — the generic (non-DC2) LSST completeness and error model.

## Verification

`docs/source/product_verification.md` records the checks these products pass:
the column contract, depth-map/anchor consistency, physical sanity of every
curve, and the DES↔DELVE cross-check, along with a per-file sha256 manifest.
