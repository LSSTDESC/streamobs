"""
Builder for precomputed background color–magnitude diagram (CMD) resources.
"""

import copy

import numpy as np
import pandas as pd

from ..columns import flag_col, obs_col
from ..surveys import Survey, SurveyFactory
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
    ...     maglim_min=24.0,
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
        # Nested dict: {source_type: {(maglim_r, maglim_g): config_dict}}
        self.resources: dict = {}
        self.bands = ("g", "r")

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
        Build CMD histograms for all ``(maglim_r, maglim_g)`` grid configurations.

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
            Two band names ``(band_g, band_r)`` where color = band_g − band_r and
            band_r is the reference magnitude axis. Default ``('g', 'r')``.
        maglim_min : float, optional
            Minimum magnitude limit for the grid. Default ``23.5``.
        maglim_max : float, optional
            Maximum magnitude limit for the grid. Default ``27.0``.
        maglim_step : float, optional
            Step size between grid values. Default ``0.5``.
        max_delta : float, optional
            Keep only pairs with ``|maglim_r − maglim_g| < max_delta``.
            Default ``1.0``.
        n_bins_color : int, optional
            Number of color histogram bins. Default ``50``.
        n_bins_mag : int, optional
            Number of magnitude histogram bins. Default ``50``.
        color_range : tuple of float, optional
            ``(min, max)`` of the color axis. Default ``(-2, 3)``.
        mag_range : tuple of float, optional
            ``(min, max)`` of the magnitude axis. Default ``(14, 30)``.
        area_ref_deg2 : float, optional
            Sky area of the reference simulation in deg².
        source_type : str, optional
            ``'stars'``, ``'galaxies'``, or ``'both'``. Default ``'both'``.
        **kwargs
            Forwarded to the catalog injector.

        Returns
        -------
        self
        """
        survey = SurveyFactory.create_survey(
            self.survey_name, release=self.release, **self._kwargs
        )
        self.bands = bands

        # Build 2-D meshgrid of (maglim_r, maglim_g) pairs
        maglim_1d = np.arange(maglim_min, maglim_max + maglim_step / 2, maglim_step)
        mg_r, mg_g = np.meshgrid(maglim_1d, maglim_1d)
        mask_delta = np.abs(mg_r - mg_g) < max_delta
        pairs = list(zip(mg_r[mask_delta].ravel(), mg_g[mask_delta].ravel()))

        active = ["stars", "galaxies"] if source_type == "both" else [source_type]

        for st in active:
            cat = load_catalog(catalog_stars if st == "stars" else catalog_galaxies)
            self.resources.setdefault(st, {})
            for maglim_r, maglim_g in pairs:
                mr_key = round(float(maglim_r), 4)
                mg_key = round(float(maglim_g), 4)
                result = self._build_one_config(
                    catalog=cat,
                    survey=survey,
                    source_type=st,
                    bands=bands,
                    maglim_r=mr_key,
                    maglim_g=mg_key,
                    n_bins_color=n_bins_color,
                    n_bins_mag=n_bins_mag,
                    color_range=color_range,
                    mag_range=mag_range,
                    area_ref_deg2=area_ref_deg2,
                    **kwargs,
                )
                self.resources[st][(mr_key, mg_key)] = result

        return self

    def _build_one_config(
        self,
        catalog,
        survey: Survey,
        source_type: str,
        bands: tuple,
        maglim_r: float,
        maglim_g: float,
        n_bins_color: int,
        n_bins_mag: int,
        color_range=(-2, 3),
        mag_range=(14, 30),
        area_ref_deg2=None,
        **kwargs,
    ) -> dict:
        """
        Build the CMD histogram for a single ``(maglim_r, maglim_g)`` pair.

        Creates a uniform survey (no dust, constant maglim), injects the catalog,
        and histograms the detected sources.

        Parameters
        ----------
        catalog : pd.DataFrame
            True background catalog (stars or galaxies).
        survey : Survey
            Loaded survey instance to prepare.
        source_type : str
            ``'stars'`` or ``'galaxies'``.
        bands : tuple of str
            ``(band_g, band_r)`` — color = band_g_obs − band_r_obs.
        maglim_r : float
            Magnitude limit for band_r.
        maglim_g : float
            Magnitude limit for band_g.
        n_bins_color : int
            Number of color bins.
        n_bins_mag : int
            Number of magnitude bins.
        color_range : tuple of float
            Color axis range.
        mag_range : tuple of float
            Magnitude axis range.
        area_ref_deg2 : float, optional
            Reference area in deg² for count normalisation.
        **kwargs
            Forwarded to the injector.

        Returns
        -------
        dict
            ``{'cmd_hist': np.ndarray, 'color_edges': np.ndarray,
            'mag_edges': np.ndarray, 'n_ref': int, 'area_ref_deg2': float}``
        """
        prepared = self._prepare_survey(
            survey,
            uniform_maglim={bands[0]: float(maglim_g), bands[1]: float(maglim_r)},
        )
        inj = BackgroundCatalogInjector(prepared)
        if source_type == "stars":
            observed = inj.inject_stars(catalog, bands=list(bands), **kwargs)
        else:
            observed = inj.inject_galaxies(catalog, bands=list(bands), **kwargs)

        namespace = prepared.namespace
        hist = self._compute_cmd_histogram(
            observed,
            namespace,
            bands,
            n_bins_color,
            n_bins_mag,
            color_range=color_range,
            mag_range=mag_range,
        )
        return {
            **hist,
            "n_ref": len(catalog),
            "area_ref_deg2": float(area_ref_deg2) if area_ref_deg2 is not None else 1.0,
        }

    def _prepare_survey(
        self,
        survey: Survey,
        no_dust: bool = False,
        uniform_maglim: dict = None,
        **kwargs,
    ) -> Survey:
        """
        Return a deep-copied, modified survey for background injection.

        Parameters
        ----------
        survey : Survey
            Survey to copy and modify.
        no_dust : bool, optional
            If True, zero the EBV map.  Default ``False``.
        uniform_maglim : dict, optional
            ``{band: value}`` mapping to replace per-pixel maglim maps with
            constant arrays.  When provided, ``no_dust`` is implied.
        **kwargs
            Reserved for future use.

        Returns
        -------
        Survey
            Modified deep copy.
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
        color_range=(-2, 3),
        mag_range=(14, 30),
    ) -> dict:
        """
        Build a 2-D color–magnitude histogram from an observed catalog.

        Color = ``mag_band_g_obs − mag_band_r_obs``; magnitude = ``mag_band_r_obs``.
        Only objects with ``flag_observed == 1`` are counted.

        Parameters
        ----------
        observed_df : pd.DataFrame
            Output of the catalog injector (observed magnitudes + flags).
        namespace : str
            Survey namespace prefix (e.g. ``'lsst_yr4'``).
        bands : tuple of str
            ``(band_g, band_r)``.
        n_bins_color, n_bins_mag : int
            Number of histogram bins on each axis.
        color_range, mag_range : tuple of float
            Axis ranges.

        Returns
        -------
        dict
            ``{'cmd_hist': np.ndarray (n_bins_color × n_bins_mag),
            'color_edges': np.ndarray, 'mag_edges': np.ndarray}``
        """
        color_edges = np.linspace(color_range[0], color_range[1], n_bins_color + 1)
        mag_edges = np.linspace(mag_range[0], mag_range[1], n_bins_mag + 1)

        mask = observed_df[flag_col(namespace)] == 1
        detected = observed_df[mask]

        if len(detected) == 0:
            return {
                "cmd_hist": np.zeros((n_bins_color, n_bins_mag)),
                "color_edges": color_edges,
                "mag_edges": mag_edges,
            }

        color = (
            detected[obs_col(bands[0], namespace)].astype(float)
            - detected[obs_col(bands[1], namespace)].astype(float)
        )
        mag = detected[obs_col(bands[1], namespace)].astype(float)

        H, xe, ye = np.histogram2d(
            color, mag,
            bins=[n_bins_color, n_bins_mag],
            range=[color_range, mag_range],
        )
        return {"cmd_hist": H, "color_edges": xe, "mag_edges": ye}

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
