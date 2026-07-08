#!/usr/bin/env python
"""
Utils for streamobs
"""

import pandas as pd
import yaml


def load_catalog(catalog):
    """Load a catalog as a :class:`pandas.DataFrame`.

    Parameters
    ----------
    catalog : pd.DataFrame or str
        A DataFrame (returned as-is) or a path to a parquet or CSV file.

    Returns
    -------
    pd.DataFrame

    Raises
    ------
    ValueError
        If ``catalog`` is neither a DataFrame nor a recognised file path.
    """
    if isinstance(catalog, pd.DataFrame):
        return catalog
    if isinstance(catalog, str):
        if catalog.endswith(".parquet"):
            return pd.read_parquet(catalog)
        if catalog.endswith(".csv"):
            return pd.read_csv(catalog)
        raise ValueError(
            f"Unrecognised file extension for catalog path '{catalog}'. "
            "Expected .parquet or .csv."
        )
    raise ValueError(
        f"catalog must be a DataFrame or a path string, got {type(catalog)}."
    )


def canonical_survey_bands(surveys, bands):
    """Sort ``(survey, band)`` pairs to a deterministic canonical order.

    Sorting key is ``(survey.name, band)`` so any permutation of the same
    inputs produces identical output.  Use this whenever building file paths,
    column names, or ordered survey/band sequences that must be commutative.

    Parameters
    ----------
    surveys : list of Survey
        One Survey per band (may repeat the same instance).
    bands : sequence of str
        Band name for each survey, parallel to *surveys*.

    Returns
    -------
    surveys_sorted : list of Survey
    bands_sorted : tuple of str
    dir_name : str
        Canonical directory name for storage paths, e.g. ``'lsst'`` or
        ``'lsst_roman'`` (survey names deduplicated, underscore-joined).
    bands_str : str
        Canonical band string for file names, e.g. ``'gr'`` or ``'gF158'``.

    Examples
    --------
    >>> canonical_survey_bands([roman, lsst], ['F158', 'g'])
    ([lsst, roman], ('g', 'F158'), 'lsst_roman', 'gF158')
    """
    raw = sorted(zip(surveys, list(bands)), key=lambda t: (t[0].name, t[1]))
    surveys_sorted = [t[0] for t in raw]
    bands_sorted = tuple(t[1] for t in raw)
    dir_name = "_".join(dict.fromkeys(s.name for s in surveys_sorted))
    bands_str = "".join(bands_sorted)
    return surveys_sorted, bands_sorted, dir_name, bands_str


def parse_config(config):
    """Parse a yaml formatted file or string into a dict.

    Parameters
    ----------
    config: yaml formatted string or file path

    Returns
    -------
    dict
    """
    try:
        # If `config` is a file
        return yaml.safe_load(open(config, "r"))
    except (OSError, FileNotFoundError):
        # Otherwise assume it is a string
        return yaml.safe_load(config)
