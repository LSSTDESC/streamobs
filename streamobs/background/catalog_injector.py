"""
Full-injection background pipeline using known catalogs.
"""

from ..observed import StreamInjector
from ..surveys import Survey
from ..utils import load_catalog


class BackgroundCatalogInjector:
    """
    Thin wrapper around :class:`~streamobs.observed.StreamInjector` for
    background catalog injection.

    Exposes :meth:`inject_stars` and :meth:`inject_galaxies` so the caller
    does not need to pass ``source_type`` explicitly.  All injection logic
    lives in :class:`~streamobs.observed.StreamInjector`.

    Parameters
    ----------
    survey : Survey
        Survey instance to inject into.
    **kwargs
        Forwarded to :class:`~streamobs.observed.StreamInjector`.

    Examples
    --------
    >>> injector = BackgroundCatalogInjector(survey)
    >>> obs_stars = injector.inject_stars(catalog_df, bands=['g', 'r'])
    >>> obs_gals  = injector.inject_galaxies(catalog_df, bands=['g', 'r'])
    """

    def __init__(self, survey: Survey, **kwargs):
        self._survey = survey
        self.streaminjector = StreamInjector(self._survey, **kwargs)

    def inject_stars(self, catalog, bands=None, **kwargs) -> "pd.DataFrame":
        """
        Inject a catalog of stars through the full survey pipeline.

        Parameters
        ----------
        catalog : pd.DataFrame or str
            True stellar catalog. Accepts a DataFrame or a path to parquet/CSV.
        bands : list of str, optional
            Photometric bands to inject. Defaults to ``['r', 'g']``.
        **kwargs
            Forwarded to :meth:`~streamobs.observed.StreamInjector.inject`.

        Returns
        -------
        pd.DataFrame
            Observed catalog with survey-namespaced magnitude and flag columns.
        """
        catalog = load_catalog(catalog)
        return self.streaminjector.inject(
            catalog, bands=bands, source_type="stars", **kwargs
        )

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
        bands : list of str, optional
            Photometric bands to inject. Defaults to ``['r', 'g']``.
        **kwargs
            Forwarded to :meth:`~streamobs.observed.StreamInjector.inject`.

        Returns
        -------
        pd.DataFrame
            Observed catalog with survey-namespaced magnitude and flag columns.
        """
        catalog = load_catalog(catalog)
        return self.streaminjector.inject(
            catalog, bands=bands, source_type="galaxies", **kwargs
        )
