# About StreamObs

## Overview

**StreamObs** is a Python package designed for stellar stream data generation and observation simulation. It provides a complete pipeline to transform theoretical or dynamical stream models into realistic mock observations as they would appear in astronomical surveys.

```{note}
StreamObs is **not** intended for generating dynamical models or N-body simulations. Instead, it provides tools to study statistical realizations of parametrized stream geometries and convert simulation outputs into mock observations.
```

## What StreamObs Does

StreamObs enables to:

1. **Generate mock stellar stream data** from parametric models
   - Create positions of stars along streams using various parametric descriptions
   - Model stream morphologies including linear, sinusoidal, and spline-based geometries
   
2. **Assign photometric properties** to stream stars
   - Attribute magnitudes in multiple photometric bands based on stellar isochrones
   
3. **Simulate survey observations** of streams
   - Convert idealized mock data into observable quantities for specific surveys (e.g., LSST)
    - Apply survey selection functions and detection limits
    - Generate realistic photometric errors and observational uncertainties

## Use Cases

StreamObs is particularly useful for people who:

- Run **dynamical simulations** of stellar streams and need to convert results into observable quantities
- Develop **stream detection algorithms** and require realistic test data with known properties
- Need to generate **mock catalogs** for testing analysis pipelines
- ...


## Development Team

StreamObs is developed and maintained by members of the LSST Dark Energy Science Collaboration (DESC).

**Core Contributors:**
- Matthieu Pélissier (Université Grenoble Alpes, CNRS/IN2P3, LPSC, 53 avenue des Martyrs, F-38026 Grenoble, France)
- Peter Ferguson (DiRAC Institute and the Department of Astronomy, University of Washington, Seattle, WA, USA)
- Alex Drlica-Wagner (Department of Astronomy and Astrophysics, University of Chicago, Chicago, IL 60637, USA)


**Contact:**
- GitHub: [LSSTDESC/streamobs](https://github.com/LSSTDESC/streamobs)
- Issues: [Report bugs or request features](https://github.com/LSSTDESC/streamobs/issues)

## License

[MIT Licence](https://opensource.org/license/MIT)


## Acknowledgments

If you use StreamObs in your research, please see the [Citation](citation.md) page for how to properly cite this software.

**Funding and Support:**
- This project is supported by the LSST Dark Energy Science Collaboration

**Related Projects:**
- [ugali](https://github.com/DarkEnergySurvey/ugali)
