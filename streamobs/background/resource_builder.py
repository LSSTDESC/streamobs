"""
Builder for precomputed background color–magnitude diagram (CMD) resources.
"""

import copy

import numpy as np

from ..columns import flag_col, obs_col
from ..surveys import Survey, SurveyFactory
from ..utils import load_catalog
from .catalog_injector import BackgroundCatalogInjector
from .storage import BackgroundStorage


class BackgroundResourceBuilder:
    """
    Build precomputed CMD histogram grids for fast background generation.

    Drives :class:`BackgroundCatalogInjector` on uniform surveys (no dust,
    constant magnitude limits) at a 2-D meshgrid of ``(maglim_b2, maglim_b1)``
    pairs (``b1 = bands[0]``, ``b2 = bands[1]``).  The resulting 2-D
    color–magnitude histograms are stored via
    :class:`BackgroundStorage` and later consumed by
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
    ...     maglim_min=24.0,
    ...     maglim_max=27.0,
    ...     maglim_step=0.5,
    ...     max_delta=1.0,
    ...     area_ref_deg2=100.0,
    ... )
    >>> storage = BackgroundStorage(survey_name='lsst', release='yr5')
    >>> builder.save(storage)
    """

    def __init__(self, survey_name="lsst", release=None, **kwargs):
        self.survey_name = survey_name
        self.release = release
        self._kwargs = kwargs
        # Nested dict: {source_type: {(maglim_b2, maglim_b1): config_dict}}
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
        Build CMD histograms for all ``(maglim_b2, maglim_b1)`` grid configurations.

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
            Two band names ``(bands[0], bands[1])`` where color = bands[0] − bands[1]
            and bands[1] is the reference magnitude axis. Default ``('g', 'r')``.
        maglim_min : float, optional
            Minimum magnitude limit for the grid. Default ``23.5``.
        maglim_max : float, optional
            Maximum magnitude limit for the grid. Default ``27.0``.
        maglim_step : float, optional
            Step size between grid values. Default ``0.5``.
        max_delta : float, optional
            Keep only pairs with ``|maglim_b2 − maglim_b1| < max_delta``.
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

        # Build 2-D meshgrid of (maglim_b2, maglim_b1) pairs
        maglim_1d = np.arange(maglim_min, maglim_max + maglim_step / 2, maglim_step)
        mg_b2, mg_b1 = np.meshgrid(maglim_1d, maglim_1d)
        mask_delta = np.abs(mg_b2 - mg_b1) < max_delta
        pairs = list(zip(mg_b2[mask_delta].ravel(), mg_b1[mask_delta].ravel()))

        active = ["stars", "galaxies"] if source_type == "both" else [source_type]

        for st in active:
            cat = load_catalog(catalog_stars if st == "stars" else catalog_galaxies)
            self.resources.setdefault(st, {})
            for maglim_b2, maglim_b1 in pairs:
                mb2_key = round(float(maglim_b2), 4)
                mb1_key = round(float(maglim_b1), 4)
                result = self._build_one_config(
                    catalog=cat,
                    survey=survey,
                    source_type=st,
                    bands=bands,
                    maglim_b2=mb2_key,
                    maglim_b1=mb1_key,
                    n_bins_color=n_bins_color,
                    n_bins_mag=n_bins_mag,
                    color_range=color_range,
                    mag_range=mag_range,
                    area_ref_deg2=area_ref_deg2,
                    **kwargs,
                )
                self.resources[st][(mb2_key, mb1_key)] = result

        return self

    def _build_one_config(
        self,
        catalog,
        survey: Survey,
        source_type: str,
        bands: tuple,
        maglim_b2: float,
        maglim_b1: float,
        n_bins_color: int,
        n_bins_mag: int,
        color_range=(-2, 3),
        mag_range=(14, 30),
        area_ref_deg2=None,
        **kwargs,
    ) -> dict:
        """
        Build the CMD histogram for a single ``(maglim_b2, maglim_b1)`` pair.

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
            ``(bands[0], bands[1])`` — color = bands[0]_obs − bands[1]_obs.
        maglim_b2 : float
            Magnitude limit for ``bands[1]`` (reference band).
        maglim_b1 : float
            Magnitude limit for ``bands[0]`` (color band).
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
            uniform_maglim={bands[0]: float(maglim_b1), bands[1]: float(maglim_b2)},
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
        no_dust: bool = True,
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
            If True, zero the EBV map.  Default ``True``.
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
        if uniform_maglim is not None or no_dust:
            if s.ebv_map is not None:
                s.ebv_map = np.zeros_like(s.ebv_map)
        if uniform_maglim is not None:
            for band, val in uniform_maglim.items():
                if band in s.maglim_maps and s.maglim_maps[band] is not None:
                    s.maglim_maps[band] = np.full_like(
                        s.maglim_maps[band], float(val), dtype=float
                    )
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

        color = detected[obs_col(bands[0], namespace)].astype(float) - detected[
            obs_col(bands[1], namespace)
        ].astype(float)
        mag = detected[obs_col(bands[1], namespace)].astype(float)

        H, xe, ye = np.histogram2d(
            color,
            mag,
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
            Which source types to save. Default ``'both'``.
        **kwargs
            Forwarded to :meth:`BackgroundStorage.save_data`.
        """
        active = ["stars", "galaxies"] if source_type == "both" else [source_type]
        for st in active:
            if st in self.resources:
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
            Which source types to load. Default ``'both'``.
        bands : tuple of str, optional
            Band names used when the resources were saved. Default ``('g', 'r')``.

        Returns
        -------
        BackgroundResourceBuilder
            Populated instance with :attr:`resources` filled.
        """
        instance = cls()
        instance.bands = bands
        active = ["stars", "galaxies"] if source_type == "both" else [source_type]
        for st in active:
            if storage.exists(st, bands):
                instance.resources[st] = storage.load_all(st, bands)
        return instance
