# Installation & Dependencies

This guide provides complete instructions for installing `streamobs` and its dependencies.

## Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/LSSTDESC/streamobs.git
cd streamobs

# 2a. Set environment variables
export PYTHONPATH=${PWD}:${PYTHONPATH}
export PATH=${PWD}/bin:${PATH}

# 2a. alternatively you can pip install by running the following int the base directory
pip install -e . 

# 3. Install ugali FROM GITHUB (not from PyPI -- see the warning below)
pip install git+https://github.com/DarkEnergySurvey/ugali.git

# 4. Download required data files
python bin/download_data.py
```

## Data Download

`streamobs` needs external data files (maglim maps, dust map, completeness, photometric errors). Use the downloader and refer to the Data page for details.

```bash
# Download required data (default location: data/)
python bin/download_data.py

# Useful options
python bin/download_data.py --list          # Show what's installed
python bin/download_data.py --force         # Re-download/overwrite
python bin/download_data.py --data-dir DIR  # Custom install location
```

For troubleshooting and data structure, see [StreamObs Data Files](data.md).

## Dependences

Required Python packages:

- ugali (**install from GitHub — see below**)
- numpy
- scipy
- pandas
- matplotlib
- astropy
- gala
- healpy
- healsparse

### Installing `ugali`

:::{warning}
**Install `ugali` from GitHub, not from PyPI.**

```bash
pip install git+https://github.com/DarkEnergySurvey/ugali.git
```

`pip install ugali` gives you release 1.8.0, which is many commits behind and
predates support for the current CMD 3.8 isochrone file format. With that
version the isochrone files are read against the wrong columns — `mass_init`
picks up the `logAge` column and becomes constant — so the mass PDF is all
zeros and sampling fails with:

```
RuntimeWarning: invalid value encountered in divide
ValueError: Probabilities contain NaN
```

The fix for this landed on `main` after the 1.8.0 tag, so a git install is
required until a newer release is published.
:::

## Optional Dependencies

- skyproj — sky-projection plots (`plotting.plot_stream_in_mask`); installable
  via the `plotting` extra: `pip install -e ".[plotting]"`


## Installing with pip/conda

```{note}
Package installation via pip/conda is planned for future releases.
```
