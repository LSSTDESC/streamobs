#!/usr/bin/env python
"""
Utils for streamobs
"""
import yaml


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
