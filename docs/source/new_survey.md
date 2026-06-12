# Adding a New Photometric Survey

To add a new survey:

## 1. Create the survey data directory

Create a new directory under `data/surveys/`:

```bash
mkdir -p data/surveys/new_survey_release
```

## 2. Add survey data files

Add magnitude limit maps for each band:

* Supported formats: HEALPix (`.hsp`) or FITS (`.fits`)
* Example:

  ```
  des_y6_g_band_nside_128.hsp
  ```
* Optional: survey-specific completeness and photometric error files. If not provided, the default files in `data/others/` will be used.

## 3. Create the survey configuration

Add a configuration file in `config/surveys/`:

```yaml
name: new_survey
release: its_release
survey_files: 
    # Path to data files (leave empty for default location)
    file_path: ''

    # Band-specific magnitude limit maps
    maglim_map_g: new_survey_maglim_g_band.hsp
    maglim_map_r: new_survey_maglim_r_band.hsp

    # Band-independent maps. Keep by defaults files for completeness, ebv map, and photometric errors
    ebv_map: ebv_sfd98_fullres_nside_4096_ring_equatorial.fits
    completeness: stellar_efficiency_cutr.csv
    completeness_band: r 
    log_photo_error: photoerror_r.csv
```

### Completeness file

Required columns:

* `delta_mag`: magnitude limit minus true magnitude (before extinction and photometric errors)
* `detection_eff`: detection efficiency
* `classification_eff`: classification efficiency
* `classification_detection_eff`: combined detection and classification efficiency

### Photometric error file

Required columns:

* `delta_mag`
* `log_mag_err`: base-10 logarithm of the magnitude error

### Completeness band

`completeness_band` specifies the band used to estimate completeness and photometric uncertainties.

## 4. Define survey properties
In the same configuration file, you must add other survey properties.

```yaml
survey_properties: 
  bands: ['u', 'g', 'r', 'i', 'z', 'y'] # Bands supported by your survey

  # Extinction coefficients per band (A_band / E(B-V)) in format coeff_extinc_bandname
  coeff_extinc_g: 3.6605664439892625
  coeff_extinc_r: 2.70136780871597
  coeff_extinc_i: 2.0536599130965882
  coeff_extinc_z: 1.5900964472616756
  coeff_extinc_y: 1.3077049588254708

  # Systematic photometric errors (mag), or could be specified as sys_error_g, sys_error_r, etc
  sys_error : 0.005 # Error common to all bands
  # sys_error_i: 0.01  # Example to Override for i-band

  # Saturation limits per band, or could be specified as saturation_g, saturation_r, etc
  saturation: 16.0
  delta_saturation: -10.4
```

where:

* `bands` lists the photometric bands available in the survey.
* `coeff_extinc_<band>` gives the extinction coefficient $A_{\rm band}/E(B-V)$.
* `sys_error` is the systematic photometric uncertainty (mag). It can also be specified per band (e.g. `sys_error_g`).
* `saturation` is the magnitude below which observations are considered saturated.
* `delta_saturation` defines the saturation limit in the completeness and photometric error tables. If saturation is not reached in the input files (ie. the completeness is not very small at small delta), use the smallest available `delta_mag` in your file.

## 5. Add survey tests

Register the survey in `SURVEY_REGISTRY` in:

```text
tests/test_surveys.py
```

## 6. Run the survey tests

```bash
pytest -k "surveys"
```

All tests must pass before merging.

## 7. Update the data archive

Update the data archive and upload the new release to Zenodo.

See: [data update page](update_data.md)

## 8. Update the documentation

Document the survey in:

* [data documentation](data.md)
* [streamobs presentation](index.rst)
