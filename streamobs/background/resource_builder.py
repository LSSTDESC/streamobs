"""
Builder for precomputed background color–magnitude diagram (CMD) resources.
"""

import numpy as np

from ..surveys import Survey, SurveyFactory
from ..utils import load_catalog
from .catalog_injector import BackgroundCatalogInjector
from .storage import BackgroundStorage


class BackgroundResourceBuilder:
    """
    Build precomputed CMD histogram grids for fast background generation.

    Drives :class:`BackgroundCatalogInjector` on uniform surveys (no dust,
    constant magnitude limits) at a grid of ``(maglim_ref, delta)`` pairs,
    where ``delta = maglim_g - maglim_ref``.  The resulting 2-D color–magnitude
    histograms are stored via :class:`BackgroundStorage` and later consumed by
    :class:`~streamobs.background.generator.LightBackgroundGenerator`.

    Parameters
    ----------
    survey_name : str, optional
        Survey identifier (e.g. ``'lsst'``).
    release : str, optional
        Survey release string (e.g. ``'yr5'``).
    **kwargs
        Forwarded to :meth:`~streamobs.surveys.SurveyFactory.create_survey`.

    Examples
    --------
    >>> builder = BackgroundResourceBuilder('lsst', release='yr5')
    >>> builder.build(
    ...     catalog_stars=df_stars,
    ...     catalog_galaxies=df_gals,
    ...     bands=('g', 'r'),
    ...     maglim_ref_values=[25.5, 26.0, 26.5],
    ...     delta_range=(-1.0, 1.0),
    ...     delta_step=0.1,
    ... )
    >>> storage = BackgroundStorage(survey_name='lsst', release='yr5')
    >>> builder.save(storage)
    """

    def __init__(self, survey_name="lsst", release=None, **kwargs):
        self.survey_name = survey_name
        self.release = release
        self._kwargs = kwargs
        # Nested dict: {source_type: {(maglim_ref, delta): config_dict}}
        self.resources: dict = {}

    def build(
        self,
        catalog_stars=None,
        catalog_galaxies=None,
        bands=("g", "r"),
        maglim_ref_values=None,
        delta_range=(-1.0, 1.0),
        delta_step=0.1,
        n_bins_color=50,
        n_bins_mag=50,
        source_type="both",
        **kwargs,
    ):
        """
        Build CMD histograms for all ``(maglim_ref, delta)`` configurations.

        For each combination of source type and magnitude limit pair, injects
        the catalog into a uniform (no-dust, constant-maglim) survey and
        computes a 2-D histogram of ``(color, mag_ref_band)``.

        Parameters
        ----------
        catalog_stars : pd.DataFrame or str, optional
            True stellar catalog.  Required when ``source_type`` is
            ``'stars'`` or ``'both'``.
        catalog_galaxies : pd.DataFrame or str, optional
            True galaxy catalog.  Required when ``source_type`` is
            ``'galaxies'`` or ``'both'``.
        bands : tuple of str, optional
            Two band names ``(band1, band2)`` where color = band1 - band2 and
            ``band2`` is the reference magnitude axis. Default ``('g', 'r')``.
        maglim_ref_values : list of float, optional
            Reference magnitude limits to sweep (for the second band in
            ``bands``).
        delta_range : tuple of float, optional
            ``(min, max)`` of ``maglim_band1 - maglim_ref``.
            Default ``(-1.0, 1.0)``.
        delta_step : float, optional
            Step size between delta values. Default ``0.1``.
        n_bins_color : int, optional
            Number of color histogram bins. Default ``50``.
        n_bins_mag : int, optional
            Number of magnitude histogram bins. Default ``50``.
        source_type : str, optional
            ``'stars'``, ``'galaxies'``, or ``'both'``. Default ``'both'``.
        **kwargs
            Forwarded to
            :meth:`BackgroundCatalogInjector.inject_stars` /
            :meth:`BackgroundCatalogInjector.inject_galaxies`.
        """
        ...

    def _build_one_config(
        self,
        catalog,
        source_type,
        bands,
        maglim_ref,
        delta,
        n_bins_color,
        n_bins_mag,
        **kwargs,
    ) -> dict:
        """
        Build the CMD histogram for a single ``(maglim_ref, delta)`` pair.

        Creates a uniform survey, injects the catalog, and histograms the
        result.

        Parameters
        ----------
        catalog : pd.DataFrame
            True background catalog (stars or galaxies).
        source_type : str
            ``'stars'`` or ``'galaxies'``.
        bands : tuple of str
            ``(band1, band2)`` — color = band1_obs - band2_obs.
        maglim_ref : float
            Magnitude limit for the reference band (``band2``).
        delta : float
            ``maglim_band1 - maglim_ref``.
        n_bins_color : int
            Number of color bins.
        n_bins_mag : int
            Number of magnitude bins.
        **kwargs
            Forwarded to the injector.

        Returns
        -------
        dict
            ``{'cmd_hist': np.ndarray, 'color_edges': np.ndarray,
            'mag_edges': np.ndarray, 'n_ref': int,
            'area_ref_deg2': float}``
        """
        ...

    def _prepare_survey(
        self,
        survey: Survey,
        no_dust: bool = False,
        uniform_maglim: dict = None,
        **kwargs,
    ) -> Survey:
        """
        Return a modified copy of the survey for background injection.

        Creates a deep copy so that modifications do not affect the original.
        Handles both the full-injection case (no-dust only) and the resource-
        building case (uniform magnitude limits + no dust).

        Parameters
        ----------
        survey : Survey
            Survey to copy and modify.
        no_dust : bool, optional
            If True, zero the EBV map. Default is False.
        uniform_maglim : dict, optional
            Mapping ``{band: value}`` to replace per-pixel magnitude limit maps
            with constant arrays. E.g. ``{'g': 26.0, 'r': 26.5}``.
            When provided, ``no_dust`` is implicitly set to True since uniform
            surveys are always dust-free by construction.
        **kwargs
            Reserved for future use.

        Returns
        -------
        Survey
            Modified copy of the survey.
        """
        ...

    def _compute_cmd_histogram(
        self,
        observed_df,
        bands: tuple,
        n_bins_color: int,
        n_bins_mag: int,
    ) -> dict:
        """
        Build a 2-D color–magnitude histogram from an observed catalog.

        Color is defined as ``mag_band1_obs - mag_band2_obs`` and the
        magnitude axis is ``mag_band2_obs``.

        Parameters
        ----------
        observed_df : pd.DataFrame
            Output of the catalog injector (observed magnitudes).
        bands : tuple of str
            ``(band1, band2)``.
        n_bins_color : int
            Number of color bins.
        n_bins_mag : int
            Number of magnitude bins.

        Returns
        -------
        dict
            ``{'cmd_hist': np.ndarray (n_bins_color × n_bins_mag),
            'color_edges': np.ndarray, 'mag_edges': np.ndarray}``
        """
        ...

    def save(self, storage: BackgroundStorage, source_type="both", **kwargs):
        """
        Persist the built resources to disk via ``storage``.

        Parameters
        ----------
        storage : BackgroundStorage
            Storage backend to use.
        source_type : str, optional
            Which source types to save: ``'stars'``, ``'galaxies'``, or
            ``'both'``. Default ``'both'``.
        **kwargs
            Forwarded to :meth:`BackgroundStorage.save_data`.
        """
        ...

    @classmethod
    def load(
        cls,
        storage: BackgroundStorage,
        source_type="both",
        **kwargs,
    ) -> "BackgroundResourceBuilder":
        """
        Load precomputed resources from disk into a new builder instance.

        Parameters
        ----------
        storage : BackgroundStorage
            Storage backend to read from.
        source_type : str, optional
            Which source types to load: ``'stars'``, ``'galaxies'``, or
            ``'both'``. Default ``'both'``.

        Returns
        -------
        BackgroundResourceBuilder
            Populated instance with :attr:`resources` filled.
        """
        ...
