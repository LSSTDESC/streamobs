"""
Full-injection background pipeline using known catalogs.
"""

import copy

import numpy as np

from ..observed import StreamInjector
from ..surveys import Survey
from ..utils import load_catalog


class BackgroundCatalogInjector:
    """
    Thin wrapper around :class:`~streamobs.observed.StreamInjector` for
    background catalog injection.

    Provides :meth:`inject_stars` and :meth:`inject_galaxies` convenience
    methods that each call :meth:`~streamobs.observed.StreamInjector.inject`
    with the correct ``source_type``. Survey preparation — deep-copying,
    zeroing dust, setting uniform magnitude limits — is handled by the single
    :meth:`_prepare_survey` method and is shared between both injection modes.

    Parameters
    ----------
    survey : Survey
        Survey instance to inject into.
    **kwargs
        Forwarded to :meth:`_prepare_survey` at construction time (e.g.
        ``no_dust``, ``uniform_maglim``).

    Examples
    --------
    Inject a stars catalog with a uniform magnitude limit:

    >>> injector = BackgroundCatalogInjector(survey, no_dust=True,
    ...                                      uniform_maglim={'g': 26.0, 'r': 26.5})
    >>> obs_stars = injector.inject_stars(catalog_df, bands=['g', 'r'])
    """

    def __init__(self, survey: Survey, **kwargs):
        self._survey = survey
        self._kwargs = kwargs

    def inject_stars(self, catalog, bands=None, **kwargs) -> "pd.DataFrame":
        """
        Inject a catalog of stars through the full survey pipeline.

        Parameters
        ----------
        catalog : pd.DataFrame or str
            True stellar catalog. Accepts a DataFrame or a path to parquet/CSV
            (loaded via :func:`~streamobs.utils.load_catalog`). Must contain
            sky coordinates (``ra``/``dec`` or ``phi1``/``phi2``) and true
            magnitude columns.
        bands : list of str, optional
            Photometric bands to inject. Defaults to ``['r', 'g']``.
        **kwargs
            Forwarded to :meth:`~streamobs.observed.StreamInjector.inject`.

        Returns
        -------
        pd.DataFrame
            Observed catalog with survey-namespaced magnitude and flag columns.
        """
        ...

    def inject_galaxies(self, catalog, bands=None, **kwargs) -> "pd.DataFrame":
        """
        Inject a catalog of galaxies through the full survey pipeline.

        Uses :meth:`~streamobs.surveys.Survey.get_gal_misclassification` for
        the detection flag (probability a galaxy passes stellar selection).
        Requires :attr:`~streamobs.surveys.Survey.gal_misclassification` to be
        loaded on the survey.

        Parameters
        ----------
        catalog : pd.DataFrame or str
            True galaxy catalog. Accepts a DataFrame or a path to parquet/CSV.
            Must contain sky coordinates and true magnitude columns.
        bands : list of str, optional
            Photometric bands to inject. Defaults to ``['r', 'g']``.
        **kwargs
            Forwarded to :meth:`~streamobs.observed.StreamInjector.inject`.

        Returns
        -------
        pd.DataFrame
            Observed catalog with survey-namespaced magnitude and flag columns.
        """
        ...
