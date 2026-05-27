# StreamObs: Stellar Stream Simulation & Observation Toolkit

[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://lsstdesc.github.io/stream_sim/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**StreamObs* is a Python package for stellar stream data generation and observation simulation. It provides tools to:

- **Generate mock stellar stream data** from parametric models. The goal is to provide a modular codebase that can be easily configured and run. The user should be able to generate a variety of stream morphologies by changing the configuration files. 
- **Assign photometric properties** to stars using stellar isochrones  
- **Simulate realistic survey observations** with observational effects

StreamObs bridges theoretical/dynamical stream models and realistic mock observations as they would appear in astronomical surveys like LSST.

> **Note**: This package is not intended to generate dynamic models or N-body simulations, but rather to study statistical realizations of parametrized stream geometries.

## Documentation

**Full documentation**: https://lsstdesc.github.io/stream_sim/

- [About](https://lsstdesc.github.io/stream_sim/about.html) - Overview and use cases
- [Quickstart](https://lsstdesc.github.io/stream_sim/quickstart.html) - Get started in minutes
- [Installation](https://lsstdesc.github.io/stream_sim/installation.html) - Detailed installation guide
- [API Reference](https://lsstdesc.github.io/stream_sim/modules.html) - Complete API documentation

## Installation

Installation consists of three steps:

```bash
# 1. Clone the repository
git clone https://github.com/LSSTDESC/stream_sim.git
cd stream_sim

# 2. Set environment variables
export PYTHONPATH=${PWD}:${PYTHONPATH}
export PATH=${PWD}/bin:${PATH}

# 3. Download required data files
python bin/download_data.py
```

The code requires common scientific Python packages:
- numpy, scipy, pandas, matplotlib
- astropy, healpy, gala
- ugali (for stellar isochrones)

For full setup (dependencies, optional tools, troubleshooting), see the
[Installation Guide](https://lsstdesc.github.io/stream_sim/installation.html).  
If data download fails or you want details on data structure, see the
[Data page](https://lsstdesc.github.io/stream_sim/data.html).

Eventually stream_sim will be installable through common package managers (i.e., pip and/or conda).

## Use Cases

StreamObs is particularly useful for:

- **Dynamical simulations**: Convert N-body simulation outputs into observable quantities
- **Algorithm development**: Generate test data for stream detection algorithms
- **Survey planning**: Predict stream detectability in upcoming surveys
- **Pipeline validation**: Test analysis workflows with known ground truth

See the [Quickstart Guide](https://lsstdesc.github.io/stream_sim/quickstart.html) for some examples.

## Citation

See the [Citation page](https://lsstdesc.github.io/stream_sim/citation.html) for more information.

