# Product verification — `des/yr6` and `delve/dr3_gold`

Verification record for the two DECam Balrog releases as shipped in
`data.zip`. Everything below is reproducible with

```bash
python scripts/verify_products.py \
    --release des_yr6=data/surveys/des_yr6 \
    --release delve_dr3_gold=data/surveys/delve_dr3_gold \
    --figdir figs/verification --manifest artifacts/product_manifest.json
```

**Result: 72 / 72 checks pass. Test suite: 331 passed, 0 failed.**

Derivation and methodology are in {doc}`balrog_selection_functions`; per-release
summaries in {doc}`surveys/DES` and {doc}`surveys/DELVE`.

## What was verified

| Group | Checks | What it establishes |
|---|---|---|
| MANIFEST | presence + sha256 of every product | the archive is complete and identified |
| CONTRACT | column names, ranges, clamp, curve orderings | streamobs will read what it expects |
| DEPTH | each map loads; nside, area, median | the depth maps are what they claim |
| PHYSICS | anchor, inflation factor, plateau, 50% crossing | the curves are physically sensible |
| CROSS | DES vs DELVE in `delta_mag` space | the two releases are mutually consistent |

## Headline numbers

| | `des/yr6` | `delve/dr3_gold` |
|---|---|---|
| injections | 145,724,947 | 62,922,015 |
| true stars binned | 4,488,693 | 13,141,646 |
| detected / classified | 86.7% / 83.6% of detected | 63.6% / 85.4% of detected |
| bands | griz | griz |
| depth map nside | 512 | 512 |
| footprint | 5,216 deg² | 16,354 deg² |
| truth-anchored m5 (g) | 25.025 | 24.373 |
| error-inflation factor | 1.464 | 1.502 |
| bright-end detection plateau | 0.999 | 0.901 |
| combined eff crosses 50% at `delta_mag` | +0.210 | −0.144 |
| star classifier | `EXT_XGB` (surrogate + deconvolution) | `bdf_extended_class_dr3gold` |
| `sys_error` | 0.005 | 0.005 |
| `delta_saturation` | −8.4 | −5.0 |

Truth-anchored depths and the shift applied to each input map's median:

| band | DES m5 | shift | DELVE m5 | shift |
|---|---|---|---|---|
| g | 25.025 | −0.297 | 24.373 | +0.196 |
| r | 24.850 | −0.280 | 23.962 | +0.287 |
| i | 24.402 | −0.180 | 23.389 | +0.185 |
| z | 23.754 | −0.124 | 22.951 | +0.391 |

Within each survey all four shifts share a sign. That coherence is the check
that validates the anchor — a band-dependent sign flip would mean the anchor was
tracking the selection rather than the depth, which is exactly the failure mode
that ruled out the `detected` anchor sample for DES (there, g and r moved
1.03 mag in opposite directions). Both releases use `anchor_sample: nosnr`.

## The cross-check

In the DES footprint DELVE DR3 Gold *is* DES Y6, so the two must agree. They are
compared in `delta_mag` space rather than on the sky, because the sky overlap is
only 369 deg² (7.1% of DES) and is edge-dominated. Both are classified by the
same `bdf_extended_class_dr3gold` interpolation nodes, which is what makes the
comparison meaningful.

Over −4 < `delta_mag` < 0:

| quantity | median abs. difference | max |
|---|---|---|
| combined stellar efficiency | 0.063 | 0.292 |
| photo-error (sample curve) | 0.020 dex | — |

```{image} _static/verification/des_vs_delve_delta_mag.png
:alt: DES Y6 and DELVE DR3 Gold combined stellar efficiency and photometric error against delta_mag
:width: 100%
```

Regenerate with `--figdir`; the script writes this figure on every run.

The residual efficiency offset is understood: DES's bright-end plateau sits at
0.999 against DELVE's 0.901 because DELVE applies per-object quality flags
(`meas_flags`, `meas_bdf_flags`) in the efficiency numerator, the Roman/LSST
convention. DES's sparse 20″ injection grid leaves almost nothing flagged
(0.056%), so its plateau reaches unity. This is a difference in what the two
numerators count, not a disagreement about the surveys.

## Findings

Two defects were found and fixed while assembling this release.

### DES depth maps were mislabelled `nside1024`

The four `des_yr6_maglim_*` maps contained nside-512 data under an
`_nside1024` filename. The input `des_y6_5_sig_*.hsp` maps are nside 512 and
`MaglimMap.to_healpix` only ever *degrades*, so `--maglim-nside 1024` was a
no-op while the output filename took the requested value regardless.

Nothing was functionally wrong — streamobs reads nside from the map, not the
name — but the files would have gone to Zenodo misdescribed. They are renamed
to `_nside512`, and `config/surveys/des_yr6.yaml`, `surveys/DES.md` and
`scripts/des/build_des_survey_doc_figs.py` updated. The verifier now asserts
that a map's nside matches its filename.

The DES depth resolution is therefore genuinely nside 512 (≈6.9′ pixels), set by
its HealSparse inputs. DELVE's maps are degraded to the same nside 512 from
nside-16384 inputs, so both DECam releases now share a depth grid.

### The DELVE photo-error curve was quantisation-limited

The truth-scatter histogram uses 0.005 mag bins, so a binned σ can only take
multiples of 0.0025. Brightward of `delta_mag` ≈ −3.26 the DELVE curve was
pinned to that grid — stepping 0.005, 0.0075, 0.010, 0.0125 and flat across
several bins at a time — i.e. reporting the bin width rather than the scatter.
The original afterburner cut at −5.2 removed only a non-monotonic wiggle and
left this region in place, so the curve floored at exactly 0.00500 mag.

The cut moved to −3.25, the first bin whose σ reaches 0.020 — the same floor the
cleaned DES curve has. This drops 23 of 64 bins from each DELVE curve and makes
`sys_error: 0.005` safe on the same footing as DES (3.1% in quadrature). Under
the old cut the same term would have been a 41% inflation of a number that was
never a measurement.

This supersedes the `sys_error` open item in the DES handoff, which proposed
setting it to zero. Zero would have been the right response to a *real* 0.005
floor; the floor was an artefact, so the correct fix was upstream in the cut.

## Known limitations

Carried forward into the release, not fixed here.

1. **DELVE's efficiency table starts at `delta_mag` −5.0** (g = 19.375), 3.4 mag
   shallower on the bright side than DES's −8.4 (g = 16.625), because brighter
   bins fall below the 20-star minimum. Stars brighter than g ≈ 19.4 are
   flat-extrapolated at `detection_eff` = 0.90.
2. **Galaxy misclassification is noise-dominated brightward of `delta_mag` ≈ −4**
   in both releases, and has no external validation in either. SPLASH-SXDF
   cannot supply one: its `STAR_FLAG` is pure but incomplete, which validates
   completeness while making contamination unmeasurable.
3. **The photo-error curves invert faintward of the depth.** For DES the sample
   curve drops below the catalog curve in 5 bins over `delta_mag` +0.24 to
   +0.72. This is the detected-population effect described in the methodology
   doc — beyond the 50% crossing only objects that scattered bright are
   recovered, compressing the measured scatter while the reported error keeps
   rising. The verifier asserts the ordering only for `delta_mag` ≤ 0 and
   reports the inversion range.
4. **`EXT_XGB` is not computable on Balrog.** DES ships a trained surrogate plus
   a per-magnitude deconvolution; DELVE ships the exactly-reproducible
   `bdf_extended_class_dr3gold` selection instead, which is a *different*
   selection from an `EXT_XGB` cut on the real catalogue. No DELVE surrogate is
   shipped.
5. **DELVE has no figures yet.** DES has six; the generator
   `scripts/des/build_des_survey_doc_figs.py` is the template.

## Manifest

`data.zip` — 49,561,196 bytes, 92 files
sha256 `08c59cab1da451e673ea7599deececc9380ee240865b9c63ae17098b02ab2286`

Per-file sizes and sha256 for both releases are in
`artifacts/product_manifest.json`. The shipped products are:

| `des/yr6` | bytes | `delve/dr3_gold` | bytes |
|---|---|---|---|
| `des_yr6_maglim_g_nside512.fits.gz` | 1,240,134 | `delve_dr3_gold_maglim_g_nside512.fits.gz` | 4,000,839 |
| `des_yr6_maglim_r_nside512.fits.gz` | 1,235,517 | `delve_dr3_gold_maglim_r_nside512.fits.gz` | 3,941,553 |
| `des_yr6_maglim_i_nside512.fits.gz` | 1,235,153 | `delve_dr3_gold_maglim_i_nside512.fits.gz` | 3,912,308 |
| `des_yr6_maglim_z_nside512.fits.gz` | 1,239,334 | `delve_dr3_gold_maglim_z_nside512.fits.gz` | 3,981,615 |
| `des_yr6_stellar_efficiency_cutg.csv` | 2,412 | `delve_dr3_gold_stellar_efficiency_cutg.csv` | 1,524 |
| `des_yr6_photoerror_g.csv` | 1,525 | `delve_dr3_gold_photoerror_g.csv` | 828 |
| `des_yr6_photoerror_g_catalog.csv` | 1,525 | `delve_dr3_gold_photoerror_g_catalog.csv` | 828 |
| `des_yr6_galaxy_misclass_cutg.csv` | 1,647 | `delve_dr3_gold_galaxy_misclass_cutg.csv` | 1,392 |
| `des_yr6_audit.json` | 868 | `delve_dr3_gold_audit.json` | 720 |

The `*_raw.csv` photo-error provenance and the `des_y6_5_sig_*.hsp` derivation
inputs are deliberately excluded from the archive; they are build inputs, not
runtime products.

## Provenance

| | |
|---|---|
| reducer | `scripts/des/balrog_selection_function.py` |
| branch | `des_and_delve_files` |
| DES catalogs | `fiducial_injected_sof.hdf5`, `fiducial_matched_measured_sof.hdf5` |
| DELVE catalog | `BalrogOfTheStars_Catalog_V4.hdf5` (62,922,015 rows) |
| DELVE depth inputs | `delve_dr32_{g,r,i,z}_maglim_wmean.hsp` + `delve_dr311+dr312_{g,r,i,z}_maglim_Nov28th.hsp` |
| DES depth inputs | `des_y6_5_sig_maglim_band_{g,r,i,z}_nside_512.hsp` |
| environment | `streamobs` conda env, Python 3.14.5, healpy + healsparse |

Both audit JSONs ship inside the archive and record row counts, the anchor per
band, the error-inflation factor and the classifier used.

## Uploading

1. Upload `archive/data.zip` to Zenodo as a new version of the record.
2. Update `BASE_DATA_URL` in `bin/download_data.py` to the new record id
   (currently `18298544`, which still serves the *old* DES products).
3. `ARCHIVE_SIZE_MB` in the same file is already updated to 48.
