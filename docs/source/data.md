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
- Zenodo record: https://zenodo.org/records/17939098
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
**URL**: https://zenodo.org/records/17939098
**Version**: 1.0  
**Last Updated**: November 2025


## Data Organization

### Required Data Files

The data directory is organized into three main categories:

#### 1. **Survey-Specific Data** (`surveys/`)

Each survey subdirectory contains magnitude limit maps (maglim maps) in HEALPix format:
- **Purpose**: Define observational depth and survey footprint for different photometric bands
- **Format**: HEALPix maps (`.hsp` files) with nside=128
- **Content**: 5σ magnitude limits for point sources in each band
- **Usage**: Used to determine which stars would be observable in a given survey

Current surveys:
- `lsst_yr1/` - LSST baseline v5.0.0, Year 1 observations (g, r bands)
- `lsst_yr2/` - LSST baseline v5.0.0, Year 2 observations (g, r bands)
- `lsst_yr3/` - LSST baseline v5.0.0, Year 3 observations (g, r bands)
- `lsst_yr4/` - LSST baseline v5.0.0, Year 4 observations (g, r bands)
- `lsst_yr5/` - LSST baseline v5.0.0, Year 5 observations (g, r bands)
- `des_yr6/` - DES Y6 Gold

Additional surveys can be added by placing maglim maps in new subdirectories.

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