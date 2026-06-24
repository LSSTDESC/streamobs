"""
Top-level background generation wrapper.
"""

import os

import pandas as pd

from ..surveys import Survey
from .catalog_injector import BackgroundCatalogInjector
from .generator import LightBackgroundGenerator
from .storage import BackgroundStorage


class Background:
    """
    High-level wrapper for background generation.

    Dispatches to :class:`BackgroundCatalogInjector` (``method='injection'``)
    or :class:`LightBackgroundGenerator` (``method='light'``) depending on the
    chosen method.  Stars and galaxies are handled independently.

    Parameters
    ----------
    survey : Survey
        Survey instance defining the observation conditions.
    source_type : {'stars', 'galaxies', 'both'}, optional
        Which components to generate. Default ``'both'``.
    method : {'light', 'injection'}, optional
        ``'light'`` uses precomputed CMD grids (fast);
        ``'injection'`` runs the full injection pipeline. Default ``'light'``.
    storage : BackgroundStorage, optional
        Resource storage for the light method.  ``None`` falls back to
        bundled package resources via :meth:`_default_storage`.
    bands : tuple of str, optional
        Band pair forwarded to the light generator. Default ``('g', 'r')``.
    catalog_stars : pd.DataFrame, optional
        Required when ``method='injection'`` and source_type includes stars.
    catalog_galaxies : pd.DataFrame, optional
        Required when ``method='injection'`` and source_type includes galaxies.
    **kwargs
        Forwarded to the underlying injector or generator.
    """

    def __init__(
        self,
        survey: Survey,
        source_type: str = "both",
        method: str = "light",
        storage: BackgroundStorage = None,
        bands: tuple = ("g", "r"),
        catalog_stars=None,
        catalog_galaxies=None,
        **kwargs,
    ):
        if source_type not in ("stars", "galaxies", "both"):
            raise ValueError(
                f"source_type must be 'stars', 'galaxies', or 'both', got '{source_type}'."
            )
        if method not in ("light", "injection"):
            raise ValueError(f"method must be 'light' or 'injection', got '{method}'.")

        self.survey = survey
        self.source_type = source_type
        self.method = method
        self.bands = bands
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
    ):
        """
        Generate a background catalog for the given sky region.

        Parameters
        ----------
        phi1_limits : tuple of float
            ``(phi1_min, phi1_max)`` in degrees.
        phi2_limits : tuple of float
            ``(phi2_min, phi2_max)`` in degrees.
        gc_frame : gala.coordinates.GreatCircleICRSFrame, optional
            Great-circle frame. Required for ``method='light'``.
        **kwargs
            Forwarded to the underlying generator or injector.

        Returns
        -------
        pd.DataFrame
            For ``method='injection'``: survey-namespaced magnitude and flag
            columns plus ``source_type``.
        tuple of (pd.DataFrame, dict)
            For ``method='light'``: catalog with ``phi1``, ``phi2``,
            ``mag_<band>``, ``source_type`` columns, plus a metadata dict.
        """
        if self.method == "injection":
            return self._generate_injection(phi1_limits, phi2_limits, gc_frame, **kwargs)
        return self._generate_light(phi1_limits, phi2_limits, gc_frame, **kwargs)

    def _generate_injection(
        self,
        phi1_limits,
        phi2_limits,
        gc_frame,
        **kwargs,
    ) -> pd.DataFrame:
        """Inject catalogs for each active source type and concatenate."""
        inj = BackgroundCatalogInjector(self.survey)
        call_kwargs = {**self._kwargs, **kwargs}
        parts = []
        if self.source_type in ("stars", "both"):
            if self.catalog_stars is None:
                raise ValueError(
                    "catalog_stars is required when method='injection' and "
                    f"source_type='{self.source_type}'."
                )
            df = inj.inject_stars(self.catalog_stars, bands=list(self.bands), **call_kwargs)
            df = df.copy()
            df["source_type"] = "stars"
            parts.append(df)
        if self.source_type in ("galaxies", "both"):
            if self.catalog_galaxies is None:
                raise ValueError(
                    "catalog_galaxies is required when method='injection' and "
                    f"source_type='{self.source_type}'."
                )
            df = inj.inject_galaxies(self.catalog_galaxies, bands=list(self.bands), **call_kwargs)
            df = df.copy()
            df["source_type"] = "galaxies"
            parts.append(df)
        del df
        return pd.concat(parts, ignore_index=True)

    def _generate_light(
        self,
        phi1_limits,
        phi2_limits,
        gc_frame,
        **kwargs,
    ):
        """Delegate to LightBackgroundGenerator."""
        gen = LightBackgroundGenerator(self.storage, self.survey, bands=self.bands)
        return gen.generate(
            phi1_limits=phi1_limits,
            phi2_limits=phi2_limits,
            gc_frame=gc_frame,
            source_type=self.source_type,
            **{**self._kwargs, **kwargs},
        )

    @staticmethod
    def _default_storage(survey: Survey) -> BackgroundStorage:
        """Return a BackgroundStorage pointing to bundled package resources."""
        return BackgroundStorage(survey_name=survey.name)
