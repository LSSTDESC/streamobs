"""
Top-level background generation wrapper.
"""

import os

from ..surveys import Survey
from .catalog_injector import BackgroundCatalogInjector
from .generator import LightBackgroundGenerator
from .storage import BackgroundStorage


class Background:
    """
    High-level wrapper for background generation.

    Dispatches to either the full injection method
    (:class:`BackgroundCatalogInjector`) or the fast light method
    (:class:`LightBackgroundGenerator`) depending on ``method``.  Stars and
    galaxies are treated independently so users can build custom combinations.

    When ``method='light'`` and ``storage=None``, the class falls back to the
    bundled package-level resource files (if present).

    Parameters
    ----------
    survey : Survey
        Survey instance defining the observation conditions.
    source_type : str, optional
        Which background components to generate: ``'stars'``, ``'galaxies'``,
        or ``'both'``. Default ``'both'``.
    method : str, optional
        Generation method: ``'light'`` (fast, uses precomputed CMD grids) or
        ``'full'`` (runs the complete injection pipeline). Default ``'light'``.
    storage : BackgroundStorage, optional
        Resource storage for the light method.  When ``None``, package-bundled
        resources are used (via :meth:`_default_storage`).
    catalog_stars : pd.DataFrame, optional
        True stellar catalog. Required when ``method='full'`` and
        ``source_type`` includes stars.
    catalog_galaxies : pd.DataFrame, optional
        True galaxy catalog. Required when ``method='full'`` and
        ``source_type`` includes galaxies.
    **kwargs
        Forwarded to the underlying injector or generator.

    Examples
    --------
    Light method with bundled defaults::

        bg = Background(survey, source_type='both', method='light')
        catalog = bg.generate(phi1_limits=(-20, 20), phi2_limits=(-2, 2), gc_frame=frame)

    Full injection with user-supplied catalogs::

        bg = Background(
            survey,
            method='full',
            source_type='stars',
            catalog_stars=df_stars,
        )
        catalog = bg.generate(phi1_limits=(-20, 20), phi2_limits=(-2, 2))
    """

    def __init__(
        self,
        survey: Survey,
        source_type: str = "both",
        method: str = "light",
        storage: BackgroundStorage = None,
        catalog_stars=None,
        catalog_galaxies=None,
        **kwargs,
    ):
        if source_type not in ("stars", "galaxies", "both"):
            raise ValueError(
                f"source_type must be 'stars', 'galaxies', or 'both', got '{source_type}'."
            )
        if method not in ("light", "full"):
            raise ValueError(f"method must be 'light' or 'full', got '{method}'.")

        self.survey = survey
        self.source_type = source_type
        self.method = method
        self.catalog_stars = catalog_stars
        self.catalog_galaxies = catalog_galaxies
        self._kwargs = kwargs

        if storage is None and method == "light":
            storage = self._default_storage(survey)
        self.storage = storage

    def generate(
        self,
        phi1_limits,
        phi2_limits,
        gc_frame=None,
        **kwargs,
    ) -> "pd.DataFrame":
        """
        Generate a background catalog for the given stream sky region.

        Parameters
        ----------
        phi1_limits : tuple of float
            ``(phi1_min, phi1_max)`` in degrees.
        phi2_limits : tuple of float
            ``(phi2_min, phi2_max)`` in degrees.
        gc_frame : gala.coordinates.GreatCircleICRSFrame, optional
            Great-circle frame. Required for the light method and when the
            full method needs coordinate conversion.
        **kwargs
            Forwarded to the underlying generator or injector.

        Returns
        -------
        pd.DataFrame
            Background catalog. Column names follow the survey namespace
            convention for the full method; for the light method columns are
            ``phi1``, ``phi2``, ``color``, ``mag``, ``source_type``.
        """
        if self.method == "full":
            return self._generate_full(phi1_limits, phi2_limits, gc_frame, **kwargs)
        return self._generate_light(phi1_limits, phi2_limits, gc_frame, **kwargs)

    def _generate_full(
        self,
        phi1_limits,
        phi2_limits,
        gc_frame,
        **kwargs,
    ) -> "pd.DataFrame":
        """
        Generate background via full catalog injection.

        For each active source type (controlled by :attr:`source_type`), uses
        :class:`BackgroundCatalogInjector` to inject the corresponding true
        catalog into the survey. The results are concatenated.

        Parameters
        ----------
        phi1_limits, phi2_limits : tuple of float
        gc_frame : gala.coordinates.GreatCircleICRSFrame or None
        **kwargs

        Returns
        -------
        pd.DataFrame
        """
        ...

    def _generate_light(
        self,
        phi1_limits,
        phi2_limits,
        gc_frame,
        **kwargs,
    ) -> "pd.DataFrame":
        """
        Generate background via the fast light method.

        Delegates to :class:`LightBackgroundGenerator` using
        :attr:`storage`.

        Parameters
        ----------
        phi1_limits, phi2_limits : tuple of float
        gc_frame : gala.coordinates.GreatCircleICRSFrame
        **kwargs

        Returns
        -------
        pd.DataFrame
        """
        ...

    @staticmethod
    def _default_storage(survey: Survey) -> BackgroundStorage:
        """
        Return a :class:`BackgroundStorage` pointing to bundled package resources.

        Parameters
        ----------
        survey : Survey

        Returns
        -------
        BackgroundStorage
        """
        return BackgroundStorage(
            survey_name=survey.name,
            release=survey.release,
        )
