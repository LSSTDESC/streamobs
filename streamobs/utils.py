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
