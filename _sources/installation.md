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

# 3. Download required data files
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

- ugali
- numpy
- scipy
- pandas
- matplotlib
- astropy
- gala
- healpy
- healsparse

## Optional Dependencies

- skyproj


## Installing with pip/conda

```{note}
Package installation via pip/conda is planned for future releases.
```
