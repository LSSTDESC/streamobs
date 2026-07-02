"""
Builder for precomputed background color–magnitude diagram (CMD) resources.
"""

import copy
import gc
import warnings

import numpy as np

from ..columns import flag_col, obs_col, true_col
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
        maglim_max=27.5,
        maglim_step=0.5,
        max_delta=1.0,
        n_bins_color=125,
        n_bins_mag=125,
        color_range=(-0.5, 2.0),
        mag_range=(16.0, 28.0),
        area_ref_deg2=None,
        source_type="both",
        verbose=True,
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
            Number of color histogram bins. Default ``125``.
        n_bins_mag : int, optional
            Number of magnitude histogram bins. Default ``125``.
        color_range : tuple of float, optional
            ``(min, max)`` of the color axis. Default ``(-0.5, 2.0)``.
        mag_range : tuple of float, optional
            ``(min, max)`` of the magnitude axis. Default ``(16.0, 28.0)``.
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
            if verbose:
                print(f"==== Processing {st}... ====")

            cat = load_catalog(catalog_stars if st == "stars" else catalog_galaxies)
            self.resources.setdefault(st, {})
            for maglim_b2, maglim_b1 in pairs:
                if verbose:
                    print(
                        f"  Building CMD histogram for {st} with "
                        f"maglim_b2={maglim_b2:.2f}, maglim_b1={maglim_b1:.2f}..."
                    )

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

            # Release the large input catalog before loading the next source type.
            del cat
            gc.collect()

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
        color_range=(-0.5, 2.0),
        mag_range=(16.0, 28.0),
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

        catalog = self._prepare_catalog(
            catalog,
            bands,
            area_ref_deg2=area_ref_deg2,
            survey=prepared,
            uniform_maglim={bands[0]: float(maglim_b1), bands[1]: float(maglim_b2)},
        )
        n_ref = len(catalog)  # capture before deletion

        inj = BackgroundCatalogInjector(prepared)
        if source_type == "stars":
            observed = inj.inject_stars(catalog, bands=list(bands), **kwargs)
        else:
            observed = inj.inject_galaxies(catalog, bands=list(bands), **kwargs)

        # Free the catalog copy and injector (holds a ref to prepared) as soon as
        # injection is done — they are not needed for the histogram step.
        del catalog
        del inj

        namespace = prepared.namespace
        del prepared  # survey deep-copy no longer needed

        hist = self._compute_cmd_histogram(
            observed,
            namespace,
            bands,
            n_bins_color,
            n_bins_mag,
            color_range=color_range,
            mag_range=mag_range,
        )
        del observed  # large observed DataFrame freed after histogram is built

        return {
            **hist,
            "n_ref": n_ref,
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

    def _prepare_catalog(
        self,
        catalog,
        bands: tuple,
        area_ref_deg2=None,
        uniform_maglim: dict = None,
        survey: Survey = None,
    ) -> "pd.DataFrame":
        """
        Prepare the catalog for injection.
        Samples positions if they are missing, and checks that the catalog covers
        approximately the requested reference area.
        Also checks that the required true magnitude columns are present.
        """
        # Validate required magnitude columns before any allocation.
        true_band1 = true_col(bands[0], survey.namespace)
        true_band2 = true_col(bands[1], survey.namespace)
        if true_band1 not in catalog.columns or true_band2 not in catalog.columns:
            raise ValueError(
                f"True background catalog must contain true magnitudes for bands {bands} "
                f"as columns '{true_band1}' and '{true_band2}'. "
                f"Available columns: {list(catalog.columns)}"
            )

        # Precompute patch bounds with spherical projection correction.
        # A naive square [0, side]×[0, side] in (ra, dec) does not have area
        # equal to side² because equal RA intervals span less sky at higher
        # declinations.  The exact solid angle of a patch [0, ra_ext]×[0, side]
        # is  ra_ext × sin(side_rad) × (180/π)  deg², so we invert this to
        # obtain ra_ext given the target area.
        if area_ref_deg2 is not None:
            side_deg = np.sqrt(area_ref_deg2)
            sin_max = np.sin(np.radians(side_deg))
            ra_extent_deg = area_ref_deg2 / (sin_max * (180.0 / np.pi))

        # Check if positions are present; if not, assign them.
        needs_positions = "ra" not in catalog.columns or "dec" not in catalog.columns

        # Only copy the columns the injector needs (true magnitudes + positions).
        # Discarding unrelated columns here avoids doubling memory for large catalogs.
        pos_cols = [] if needs_positions else ["ra", "dec"]
        cat = catalog[[true_band1, true_band2] + pos_cols].copy()

        if needs_positions:
            if uniform_maglim is not None:
                if area_ref_deg2 is None:
                    raise ValueError(
                        "area_ref_deg2 must be provided when catalog has no positions."
                    )
                # Sample dec uniformly in solid angle (i.e. uniform in sin(dec))
                # so the distribution is isotropic on the sphere, not compressed
                # toward dec=0.  RA is sampled over the projection-corrected
                # extent so the enclosed solid angle equals area_ref_deg2.
                cat["dec"] = np.degrees(
                    np.arcsin(np.random.uniform(0.0, sin_max, size=len(cat)))
                )
                cat["ra"] = np.random.uniform(0.0, ra_extent_deg, size=len(cat))
            else:
                raise ValueError(
                    "Positions are required when the survey is not uniform."
                )

        if area_ref_deg2 is not None:
            area_estimated, pixel_area_deg2, nside_used = self._estimate_area_deg2(
                cat["ra"].to_numpy(), cat["dec"].to_numpy(), area_ref_deg2
            )
            rel_error = abs(area_estimated - area_ref_deg2) / area_ref_deg2
            tolerance = 0.05
            if rel_error > tolerance:
                msg = (
                    f"The sky area estimated from catalog positions "
                    f"({area_estimated:.2f} deg²) differs from "
                    f"area_ref_deg2={area_ref_deg2:.2f} deg² by "
                    f"{rel_error * 100:.1f} % (> {tolerance * 100:.1f} % tolerance; "
                    f"HEALPix nside={nside_used}, "
                    f"pixel area={pixel_area_deg2:.4f} deg²). "
                )
                if needs_positions:
                    warnings.warn(
                        msg
                        + "Positions were sampled within the correct patch by construction; "
                        "the discrepancy is a HEALPix boundary effect (too few objects "
                        "or too small a footprint for the pixel size).",
                        UserWarning,
                        stacklevel=3,
                    )
                else:
                    raise ValueError(
                        msg
                        + "Ensure the catalog covers approximately the requested reference area."
                    )

        return cat

    def _estimate_area_deg2(self, ra, dec, area_ref_deg2):
        """Estimate the sky area covered by (ra, dec) positions using HEALPix.

        Counts the number of unique HEALPix pixels that contain at least one
        catalog object and multiplies by the pixel solid angle.  This gives a
        correct area estimate for any footprint shape (rectangular, circular,
        irregular), unlike a bounding-box approach which would systematically
        overestimate a disc by 4/π ≈ 27 %.

        The nside is set so the expected number of objects per pixel
        (assuming the catalog uniformly fills ``area_ref_deg2``) is ~10.
        This keeps most pixels occupied (low undercounting) while keeping
        pixels small enough to bound overestimation from boundary pixels.
        At least 10 pixels are always targeted.  The result is a tuple
        ``(area_deg2, pixel_area_deg2, nside)``.

        Parameters
        ----------
        ra, dec : array-like
            Positions in degrees.
        area_ref_deg2 : float
            Expected area used to set the HEALPix resolution.

        Returns
        -------
        area_estimated : float
            Area covered by the catalog in deg².
        pixel_area_deg2 : float
            Solid angle of one HEALPix pixel at the chosen nside, in deg².
        nside : int
            HEALPix nside used for the estimate.
        """
        import healpy as hp

        # Choose nside so the footprint contains on average 10 objects per pixel.
        # This balances coverage (pixels are unlikely to be empty) against
        # boundary overestimation (fewer, larger pixels → thicker boundary ring).
        # At least 10 pixels are always targeted so tiny footprints are resolved.
        # pixel_area_deg2 ≈ 41 253 / (12 × nside²) → solve for nside.
        n = max(len(ra), 1)
        pixel_area_target = area_ref_deg2 / max(n / 10.0, 10.0)
        nside_raw = np.sqrt(41253.0 / (12.0 * max(pixel_area_target, 1e-4)))
        nside = int(2 ** np.round(np.log2(max(nside_raw, 4))))
        nside = max(4, min(nside, 4096))

        ra = np.asarray(ra)
        dec = np.asarray(dec)
        finite = np.isfinite(ra) & np.isfinite(dec)
        theta = np.radians(90.0 - dec[finite])
        phi = np.radians(ra[finite])

        pixels = hp.ang2pix(nside, theta, phi)
        n_unique = len(np.unique(pixels))
        pixel_area = hp.nside2pixarea(nside, degrees=True)

        return n_unique * pixel_area, pixel_area, nside

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
        if not mask.any():
            return {
                "cmd_hist": np.zeros((n_bins_color, n_bins_mag)),
                "color_edges": color_edges,
                "mag_edges": mag_edges,
            }

        # Extract only the two needed columns filtered to detected rows — avoids
        # materialising a full-width filtered copy of the observed DataFrame.
        # The mask is applied BEFORE converting to float so that non-numeric
        # sentinel values in undetected rows ('BAD_MAG') are never encountered.
        col_b1 = obs_col(bands[0], namespace)
        col_b2 = obs_col(bands[1], namespace)
        mag = observed_df.loc[mask, col_b2].to_numpy(dtype=float)
        color = observed_df.loc[mask, col_b1].to_numpy(dtype=float) - mag

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
