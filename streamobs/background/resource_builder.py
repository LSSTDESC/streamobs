"""
Builder for precomputed background color–magnitude diagram (CMD) resources.
"""

import copy

import numpy as np
import pandas as pd

from ..surveys import Survey
from ..utils import load_catalog
from .catalog_injector import BackgroundCatalogInjector
from .storage import BackgroundStorage


class BackgroundResourceBuilder:
    """
    Build precomputed CMD histogram grids for fast background generation.

    Drives :class:`BackgroundCatalogInjector` on uniform surveys (no dust,
    constant magnitude limits) at a 2D grid of ``(maglim_r, maglim_g)`` pairs
    filtered by ``|maglim_r - maglim_g| < max_delta``.  The resulting 2-D
    color–magnitude histograms are stored via :class:`BackgroundStorage` and
    later consumed by
    :class:`~streamobs.background.generator.LightBackgroundGenerator`.

    Parameters
    ----------
    survey_name : str, optional
        Survey identifier (e.g. ``'lsst'``).
    **kwargs
        Forwarded to :meth:`~streamobs.surveys.Survey.load` when ``build()``
        loads the survey (e.g. ``release='yr4'``).

    Examples
    --------
    >>> builder = BackgroundResourceBuilder('lsst', release='yr4')
    >>> builder.build(
    ...     catalog_stars=df_stars,
    ...     bands=('g', 'r'),
    ...     maglim_min=23.5,
    ...     maglim_max=27.0,
    ...     maglim_step=0.5,
    ...     max_delta=1.0,
    ...     area_ref_deg2=100.0,
    ... )
    >>> storage = BackgroundStorage(survey_name='lsst')
    >>> builder.save(storage)
    """

    def __init__(self, survey_name="lsst", **kwargs):
        self.survey_name = survey_name
        self._kwargs = kwargs
        # {source_type: {(maglim_r, maglim_g): result_dict}}
        self.resources: dict = {}
        self.bands = None

    def build(
        self,
        catalog_stars=None,
        catalog_galaxies=None,
        bands=("g", "r"),
        maglim_min=23.5,
        maglim_max=27.0,
        maglim_step=0.5,
        max_delta=1.0,
        n_bins_color=50,
        n_bins_mag=50,
        color_range=(-2, 3),
        mag_range=(14, 30),
        area_ref_deg2=None,
        source_type="both",
        **kwargs,
    ):
        """
        Build CMD histograms for all ``(maglim_r, maglim_g)`` grid pairs.

        Creates a 2D meshgrid of magnitude limits over ``[maglim_min,
        maglim_max]`` for both bands, retains only pairs satisfying
        ``|maglim_r - maglim_g| < max_delta``, then injects the catalog into
        a dust-free uniform survey for each retained pair.

        Parameters
        ----------
        catalog_stars : pd.DataFrame or str, optional
            True stellar catalog. Required when ``source_type`` is
            ``'stars'`` or ``'both'``.
        catalog_galaxies : pd.DataFrame or str, optional
            True galaxy catalog. Required when ``source_type`` is
            ``'galaxies'`` or ``'both'``.
        bands : tuple of str, optional
            ``(band_color, band_ref)`` where color = band_color - band_ref and
            ``band_ref`` is the magnitude axis. Default ``('g', 'r')``.
        maglim_min : float, optional
            Minimum magnitude limit for both bands. Default ``23.5``.
        maglim_max : float, optional
            Maximum magnitude limit for both bands. Default ``27.0``.
        maglim_step : float, optional
            Grid step for magnitude limits. Default ``0.5``.
        max_delta : float, optional
            Maximum allowed ``|maglim_r - maglim_g|``. Default ``1.0``.
        n_bins_color : int, optional
            Number of color histogram bins. Default ``50``.
        n_bins_mag : int, optional
            Number of magnitude histogram bins. Default ``50``.
        color_range : tuple of float, optional
            ``(min, max)`` for the color axis. Identical across all grid
            points so the generator can interpolate. Default ``(-2, 3)``.
        mag_range : tuple of float, optional
            ``(min, max)`` for the magnitude axis. Default ``(14, 30)``.
        area_ref_deg2 : float
            Sky area of the input catalog in deg². Stored with every grid
            point so the generator can scale counts to arbitrary pixel sizes.
        source_type : str, optional
            ``'stars'``, ``'galaxies'``, or ``'both'``. Default ``'both'``.
        **kwargs
            Forwarded to
            :meth:`BackgroundCatalogInjector.inject_stars` /
            :meth:`BackgroundCatalogInjector.inject_galaxies`.

        Returns
        -------
        self
            Returns the builder for chaining.
        """
        if area_ref_deg2 is None:
            raise ValueError("area_ref_deg2 must be provided (sky area of the input catalog in deg²).")

        survey = Survey.load(self.survey_name, **self._kwargs)
        self.bands = list(bands)

        # Resolve active (source_type, catalog) pairs
        active = []
        if source_type in ("stars", "both"):
            if catalog_stars is None:
                raise ValueError("catalog_stars is required for source_type='stars' or 'both'.")
            active.append(("stars", load_catalog(catalog_stars)))
        if source_type in ("galaxies", "both"):
            if catalog_galaxies is None:
                raise ValueError("catalog_galaxies is required for source_type='galaxies' or 'both'.")
            active.append(("galaxies", load_catalog(catalog_galaxies)))

        # Build 2D meshgrid and filter by delta threshold
        maglim_1d = np.arange(maglim_min, maglim_max + maglim_step / 2, maglim_step)
        mg_r, mg_g = np.meshgrid(maglim_1d, maglim_1d)
        mask = np.abs(mg_r - mg_g) < max_delta
        pairs = list(zip(mg_r[mask].ravel(), mg_g[mask].ravel()))

        for st, catalog in active:
            self.resources.setdefault(st, {})
            for maglim_r, maglim_g in pairs:
                result = self._build_one_config(
                    catalog,
                    survey,
                    st,
                    bands,
                    maglim_r,
                    maglim_g,
                    n_bins_color,
                    n_bins_mag,
                    color_range,
                    mag_range,
                    area_ref_deg2,
                    **kwargs,
                )
                self.resources[st][(maglim_r, maglim_g)] = result

        return self

    def _build_one_config(
        self,
        catalog,
        survey,
        source_type,
        bands,
        maglim_r,
        maglim_g,
        n_bins_color,
        n_bins_mag,
        color_range,
        mag_range,
        area_ref_deg2,
        **kwargs,
    ) -> dict:
        """
        Build the CMD histogram for a single ``(maglim_r, maglim_g)`` pair.

        Parameters
        ----------
        catalog : pd.DataFrame
            True background catalog (stars or galaxies).
        survey : Survey
            Unmodified loaded survey — ``_prepare_survey`` makes the copy.
        source_type : str
            ``'stars'`` or ``'galaxies'``.
        bands : tuple of str
            ``(band_color, band_ref)``.
        maglim_r : float
            Magnitude limit for the reference band (``bands[1]``).
        maglim_g : float
            Magnitude limit for the colour band (``bands[0]``).
        n_bins_color : int
            Number of color bins.
        n_bins_mag : int
            Number of magnitude bins.
        color_range : tuple of float
            Fixed color axis range.
        mag_range : tuple of float
            Fixed magnitude axis range.
        area_ref_deg2 : float
            Sky area of the catalog in deg².
        **kwargs
            Forwarded to the injector.

        Returns
        -------
        dict
            ``{'cmd_hist': np.ndarray, 'color_edges': np.ndarray,
            'mag_edges': np.ndarray, 'n_ref': int, 'area_ref_deg2': float}``
        """
        prepared = self._prepare_survey(
            survey, uniform_maglim={bands[0]: maglim_g, bands[1]: maglim_r}
        )
        inj = BackgroundCatalogInjector(prepared)
        if source_type == "stars":
            obs = inj.inject_stars(catalog, bands=list(bands), **kwargs)
        else:
            obs = inj.inject_galaxies(catalog, bands=list(bands), **kwargs)

        hist = self._compute_cmd_histogram(
            obs, survey.namespace, bands, n_bins_color, n_bins_mag, color_range, mag_range
        )
        return {**hist, "n_ref": len(catalog), "area_ref_deg2": area_ref_deg2}

    def _prepare_survey(
        self,
        survey: Survey,
        no_dust: bool = False,
        uniform_maglim: dict = None,
        **kwargs,
    ) -> Survey:
        """
        Return a modified deep copy of the survey for background injection.

        Parameters
        ----------
        survey : Survey
            Survey to copy and modify.
        no_dust : bool, optional
            If True, zero the EBV map. Default ``False``.
        uniform_maglim : dict, optional
            ``{band: value}`` mapping to replace per-pixel magnitude limit maps
            with constant arrays, e.g. ``{'g': 26.0, 'r': 26.5}``.
            When provided, ``no_dust`` is implicitly set to True.
        **kwargs
            Reserved for future use.

        Returns
        -------
        Survey
            Modified copy of the survey.
        """
        s = copy.deepcopy(survey)
        if uniform_maglim is not None:
            no_dust = True
            for band, value in uniform_maglim.items():
                if band in s.maglim_maps:
                    s.maglim_maps[band] = np.full_like(s.maglim_maps[band], value)
        if no_dust and s.ebv_map is not None:
            s.ebv_map = np.zeros_like(s.ebv_map)
        return s

    def _compute_cmd_histogram(
        self,
        observed_df,
        namespace: str,
        bands: tuple,
        n_bins_color: int,
        n_bins_mag: int,
        color_range: tuple,
        mag_range: tuple,
    ) -> dict:
        """
        Build a 2-D color–magnitude histogram from an observed catalog.

        Color is ``mag_bands[0]_obs - mag_bands[1]_obs``; the magnitude axis
        is ``mag_bands[1]_obs``.  Only detected sources (flag == 1) are
        histogrammed.

        Parameters
        ----------
        observed_df : pd.DataFrame
            Output of the catalog injector.
        namespace : str
            Survey column namespace, e.g. ``'lsst_yr4'``. Needed to find the
            columns of interest.
        bands : tuple of str
            ``(band_color, band_ref)``.
        n_bins_color : int
            Number of color bins.
        n_bins_mag : int
            Number of magnitude bins.
        color_range : tuple of float
            Fixed ``(min, max)`` for the color axis.
        mag_range : tuple of float
            Fixed ``(min, max)`` for the magnitude axis.

        Returns
        -------
        dict
            ``{'cmd_hist': np.ndarray (n_bins_color × n_bins_mag),
            'color_edges': np.ndarray, 'mag_edges': np.ndarray}``
        """
        flag_col = f"{namespace}_flag_observed"
        col0 = f"{namespace}_{bands[0]}_obs"
        col1 = f"{namespace}_{bands[1]}_obs"

        detected = observed_df[observed_df[flag_col] == 1]

        mag1 = pd.to_numeric(detected[col1], errors='coerce')
        mag0 = pd.to_numeric(detected[col0], errors='coerce')
        color = (mag0 - mag1).values
        mag = mag1.values

        H, xedges, yedges = np.histogram2d(
            color,
            mag,
            bins=[n_bins_color, n_bins_mag],
            range=[color_range, mag_range],
        )
        return {"cmd_hist": H, "color_edges": xedges, "mag_edges": yedges}

    def save(self, storage: BackgroundStorage, source_type="both", **kwargs):
        """
        Persist the built resources to disk via ``storage``.

        Parameters
        ----------
        storage : BackgroundStorage
            Storage backend to use.
        source_type : str, optional
            ``'stars'``, ``'galaxies'``, or ``'both'``. Default ``'both'``.
        **kwargs
            Forwarded to :meth:`BackgroundStorage.save_data`.
        """
        if self.bands is None:
            raise RuntimeError("No resources to save — call build() first.")
        types = ["stars", "galaxies"] if source_type == "both" else [source_type]
        for st in types:
            if st not in self.resources:
                continue
            storage.save_data(self.resources[st], st, self.bands, **kwargs)

    @classmethod
    def load(
        cls,
        storage: BackgroundStorage,
        source_type="both",
        bands=("g", "r"),
        **kwargs,
    ) -> "BackgroundResourceBuilder":
        """
        Load precomputed resources from disk into a new builder instance.

        Parameters
        ----------
        storage : BackgroundStorage
            Storage backend to read from.
        source_type : str, optional
            ``'stars'``, ``'galaxies'``, or ``'both'``. Default ``'both'``.
        bands : tuple of str, optional
            Band names used when the resources were built. Default ``('g', 'r')``.

        Returns
        -------
        BackgroundResourceBuilder
            Populated instance with :attr:`resources` filled.
        """
        instance = cls(survey_name=storage.survey_name)
        instance.bands = list(bands)
        types = ["stars", "galaxies"] if source_type == "both" else [source_type]
        for st in types:
            data = storage.load_all(st, bands, **kwargs)
            if data:
                instance.resources[st] = data
        return instance
