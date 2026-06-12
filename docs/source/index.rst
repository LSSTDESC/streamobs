.. streamobs documentation master file

Welcome to StreamObs's Documentation!
======================================

.. warning::
   This package is under active development. The API and features may change in future releases.

**StreamObs** is a Python package for stellar stream data generation and observation simulation. 
It provides a complete pipeline to transform theoretical or dynamical stream models into realistic 
mock observations as they would appear in astronomical surveys.
It now supports the following surveys (and corresponding data releases):
- **DES**: Y6
- **LSST**: simulated properties for Y1, Y2, Y3, Y4, and Y5


What StreamObs Does
-------------------

StreamObs enables you to:

1. **Generate mock stellar stream data** from parametric models
   
   - Create positions of stars along streams (phi1, phi2 or RA, Dec)
   - Model various stream morphologies (linear, sinusoidal, spline-based)
   - Assign distances and velocities to stream stars

2. **Assign photometric properties** using stellar isochrones
   
   - Generate magnitudes in multiple photometric bands (eg. g and r)

3. **Simulate realistic survey observations**
   
   - Apply survey footprint and selection functions
   - Include Galactic extinction effects
   - Add photometric errors and detection completeness
   - Simulate observational uncertainties

Use Cases
---------

StreamObs is particularly useful for people who:

- Run **dynamical simulations** and need to convert results into observable quantities
- Develop **stream detection algorithms** and require realistic test data
- Plan **survey strategies** for detecting stellar streams
- Need to generate **mock catalogs** for pipeline validation

.. note::
   StreamObs focuses on parametric stream models and observation simulation. 
   It is **not** designed for N-body simulations or dynamical modeling.


Getting Started
---------------

New to StreamObs? Start here:

1. **Installation**: :doc:`installation` - Set up StreamObs on your system
2. **About**: :doc:`about` - Understand StreamObs's purpose and design
3. **Quickstart**: :doc:`quickstart` - Learn the basics with simple examples
4. **API Reference**: :doc:`modules` - Explore the complete API documentation

.. tip::
   Having issues with data download? See the :doc:`data` page for detailed troubleshooting.

For developers
---------------

If you're interested in contributing to StreamObs or using it as a base for your own development, check out these pages:

1. **Modify Streamobs**: :doc:`modifying_streamobs` - Guidlines for modifying the core codebase
2. **Adding a new survey**: :doc:`new_survey` - Learn how to add support for a new photometric survey
3. **Update data base**: :doc:`update_data` - Learn how to update the data base

Documentation Contents
----------------------

.. toctree::
   :maxdepth: 2
   :caption: User Guide
   
   about
   installation
   quickstart
   citation

.. toctree::
   :maxdepth: 2
   :caption: For developers
   
   modifying_streamobs
   new_survey
   update_data

.. toctree::
   :maxdepth: 2
   :caption: Reference

   data
   modules

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
