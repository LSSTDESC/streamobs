# StreamObs Data Files

This directory contains large data files required for stream simulations. These files are **not** tracked in the git repository due to their size. They must be downloaded separately from Zenodo.

## Downloading Data

### Quick Start

After cloning the repository, download all required data files:

```bash
python bin/download_data.py
```

The script will:
1. Download a compressed archive from Zenodo
2. Extract it to the `data/` directory
3. Clean up temporary and system files
4. Verify the installation

### Available Commands

#### Check Current Data

View what data is currently installed:

```bash
python bin/download_data.py --list
```

Output shows:
- Subdirectories in the data folder
- Number of files per subdirectory
- Total size of installed data

#### Force Re-download

If data is corrupted or you need to update:

```bash
python bin/download_data.py --force
```

This will:
- Re-download the archive even if data exists
- Overwrite existing files
- Clean up unwanted files

#### Custom Data Location

Specify a different data directory:

```bash
python bin/download_data.py --data-dir /path/to/custom/data
```

#### Keep Archive

Save the downloaded zip file after extraction:

```bash
python bin/download_data.py --keep-archive
```

The archive will be saved as `data.zip` in the repository root.

#### Custom Data URL

Use a different data source:

```bash
python bin/download_data.py --url https://custom-server.edu/data.zip
```

### Troubleshooting Data Download

#### Problem: Download fails with "404 Not Found"

**Solution**: The data URL may have changed. Check the latest URL at:
- Zenodo record: the one `BASE_DATA_URL` in `bin/download_data.py` points at
  (that file is the single source of truth for which record is current)
- Or update `BASE_DATA_URL` in `bin/download_data.py`

#### Problem: Extraction fails

**Solution**: 
1. Check disk space: The extracted data requires ~ 800 MB
2. Check write permissions in the installation directory
3. Try re-downloading with `--force`

#### Problem: Data directory is empty after download

**Solution**:
1. Run `python bin/download_data.py --list` to check status
2. Verify the archive was extracted correctly
3. Check for error messages during extraction

#### Problem: Missing specific survey data

**Solution**:
1. Verify which surveys are included: `python bin/download_data.py --list`
2. If a survey is missing, check if it's in the Zenodo archive
3. You may need to download additional survey-specific data separately

## Data Storage and DOI

The data files are hosted on [Zenodo](https://zenodo.org) with a persistent DOI for citation and long-term access.

**DOI**: 10.5281/zenodo.17550956

**URL**: whichever record `BASE_DATA_URL` in `bin/download_data.py` names. Do not
hardcode a record id here — it has drifted from the code before. As of this
writing the code points at `18298544`, which still serves the *previous* product
set; the archive described in {doc}`product_verification` has not yet been
uploaded, and `BASE_DATA_URL` must be bumped when it is.
**Version**: 1.0  
**Last Updated**: see the record itself; the shipped product set is the one
manifested in {doc}`product_verification`.


## Data Organization

### Required Data Files

The data directory is organized into three main categories:

#### 1. **Survey-Specific Data** (`surveys/`)

Each survey subdirectory holds that release's magnitude-limit (maglim) maps plus
its selection-function tables:
- **Purpose**: define observational depth and footprint per band, and the
  detection / classification / photometric-error model keyed to `delta_mag`
- **Format**: HEALPix maps as either HealSparse `.hsp` or gzipped FITS
  `.fits.gz`, depending on the release; tables as CSV
- **Content**: 5σ point-source magnitude limits per band, and the product set
  described in {doc}`selection_function_methodology`
- **Usage**: determines which stars would be observable, and with what
  completeness and photometric error

Current releases, with the resolution and format of their maglim maps:

| directory | release | bands | maglim map |
|---|---|---|---|
| `des_yr6/` | DES Y6 Gold | griz | nside 128, `.fits.gz` |
| `delve_dr3_gold/` | DELVE DR3 Gold | griz | nside 128, `.fits.gz` |
| `lsst_dc2/` | LSST DC2 | g, r | nside 128, `.fits.gz` |
| `lsst_yr1/` … `lsst_yr5/` | LSST baseline v5.0.0, years 1–5 | g, r | nside 128, `.hsp` |
| `lsst_dp2/` | LSST DP2 | g, r | nside 128, `.hsp` |
| `roman_dc2/` | Roman DC2 | F106, F129, F158 | nside 128, `.fits.gz` |
| `roman_hlwas_wide/`, `_medium/`, `_all/` | Roman HLWAS tiers | F158 (F106 for `_all`) | nside 128, `.fits.gz` |

The DECam releases (`des_yr6`, `delve_dr3_gold`) are derived from Balrog
synthetic-source injections — see {doc}`surveys/DES`, {doc}`surveys/DELVE` and
{doc}`balrog_selection_functions`. The LSST and Roman releases are described in
{doc}`surveys/LSST` and {doc}`surveys/Roman`.

Additional surveys can be added by placing their products in a new
subdirectory — see {doc}`new_survey`.

#### 2. **Auxiliary Data** (`others/`)

Common data files required for all simulations:

- **Dust Extinction Map** (`ebv_sfd98_fullres_nside_4096_ring_equatorial.fits`):
  - E(B-V) values from Schlegel, Finkbeiner & Davis (1998)
  - Full-resolution HEALPix map (nside=4096)
  - Used to apply Galactic extinction corrections to stellar magnitudes

- **Survey Completeness** (`stellar_efficiency_cutr.csv`):
  - Detection and classification efficiencies as a function of difference between apparent magnitude and magnitude limit
  - Accounts for photometric pipeline completeness
  - Used to model realistic detection probabilities

- **Photometric Errors** (`photoerror_r.csv`):
  - Photometric uncertainties as a function of difference between apparent magnitude and magnitude limit
  - Used to add realistic observational noise to simulated photometry

#### 3. **Stream Models** (root directory)

Reference data for specific stream models:

- `erkal_2016_pal_5_input.csv` 
- `patrick_2022_splines.csv`

These are small reference files (<100 KB) and are tracked in git.

## For Developers

Informations to modify the data base can be found in [Update data page](update_data.md), which can be usefull to add [new survey](new_survey.md).