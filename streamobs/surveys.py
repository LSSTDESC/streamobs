import copy
import os
from dataclasses import dataclass
from typing import Callable, Optional

import healpy as hp
import healsparse as hsp
import numpy as np
import scipy.interpolate


@dataclass
class Survey:
    """
    Container for survey properties and data.

    This class stores all survey-specific data including HEALPix maps,
    efficiency functions, and photometric error models. Data is loaded
    once and reused to avoid redundant file I/O.

    Band-specific properties (magnitude limits, extinction, errors, etc.)
    are stored in dictionaries keyed by band name (e.g., 'g', 'r', 'i').

    Attributes
    ----------
    name : str
        Survey identifier (e.g., 'lsst').
    release : str, optional
        Survey data release version (e.g., 'yr1', 'yr5').
    bands : list of str, optional
        List of photometric bands available (e.g., ['g', 'r', 'i']).
        Default is an empty list.
    maglim_maps : dict, optional
        HEALPix maps of magnitude limits, keyed by band.
        Example: {'g': array(...), 'r': array(...)}.
        Default is an empty dictionary.
    coeff_extinc : dict, optional
        Extinction coefficients A_band/E(B-V), keyed by band.
        Example: {'g': 3.303, 'r': 2.285}.
        Default uses LSST values for standard bands.
    saturation : dict, optional
        Bright magnitude limits where detector saturates, keyed by band.
        Example: {'g': 16.0, 'r': 16.0}.
        Default is 16.0 mag for all bands.
    sys_error : dict, optional
        Systematic photometric errors (mag), keyed by band.
        Example: {'g': 0.01, 'r': 0.01}.
        Default is 0.0 for all bands.
    ebv_map : np.ndarray, optional
        HEALPix map of E(B-V) extinction values (band-independent).
    coverage : np.ndarray, optional
        HEALPix map of survey coverage (1=observed, 0=not observed).
    completeness : callable, optional
        Efficiency function f(magnitude) -> efficiency [0, 1].
        Same function used for all bands, obtained from r band.
    delta_saturation : float, optional
        Magnitude difference for saturation threshold in the initial functions.
    completeness_band : str, optional
        Band used to derive completeness function (e.g., 'r').
    log_photo_error : callable, optional
        Photometric error model f(delta_mag) -> log10(mag_error).
        Same function used for all bands, obtained from r band.

    Examples
    --------
    Load a survey:

    >>> survey = Survey.load('lsst', release='dc2')
    >>> print(survey.bands)
    ['g', 'r']

    Access magnitude limit for g-band at pixel 10000:

    >>> maglim_g = survey.get_maglim('g', pixel=10000)
    >>> print(maglim_g)

    Calculate extinction in r-band at pixel 10000:

    >>> extinction_r = survey.get_extinction('r', pixel=10000)
    >>> print(extinction_r)

    Get detection efficiency for magnitude 24.0 in g-band:

    >>> maglim_r = survey.get_maglim('r', pixel=10000)
    >>> efficiency_g = survey.get_completeness('r', 24.0, maglim_r)
    >>> print(efficiency_g)

    Calculate photometric error for magnitude 24.0 in r-band:

    >>> maglim_r = survey.get_maglim('r', pixel=10000)
    >>> photo_error_r = survey.get_photo_error('r', 24.0, maglim_r)
    >>> print(photo_error_r)
    """

    # Survey identification
    name: str
    release: Optional[str] = None

    # Available bands
    bands: list = None

    # Band-specific maps and properties (dictionaries keyed by band)
    maglim_maps: dict = None
    coeff_extinc: dict = None
    saturation: dict = None
    sys_error: dict = None

    # Band-independent functions (same for all bands)
    completeness: Optional[Callable] = None
    completeness_band: Optional[str] = None
    delta_saturation: Optional[float] = None
    log_photo_error: Optional[Callable] = None

    # Band-independent maps
    ebv_map: Optional[np.ndarray] = None
    coverage: Optional[np.ndarray] = None

    def __post_init__(self):
        """Initialize default dictionaries if not provided."""
        if self.bands is None:
            self.bands = []
        if self.maglim_maps is None:
            self.maglim_maps = {}
        if self.coeff_extinc is None:
            self.coeff_extinc = {}
        if self.saturation is None:
            self.saturation = {}
        if self.sys_error is None:
            self.sys_error = {}

    @classmethod
    def load(
        cls,
        survey: str,
        release: Optional[str] = None,
        config_file: Optional[dict] = None,
        **kwargs,
    ) -> "Survey":
        """
        Load survey data and return Survey instance.

        This is a convenience method that calls SurveyFactory.create_survey().

        Parameters
        ----------
        survey : str
            Survey name (e.g., 'lsst').
        release : str, optional
            Survey release/version (e.g., 'yr1', 'yr5').
        config_file : dict, optional
            Custom configuration dictionary.
        **kwargs
            Additional keyword arguments passed to ``SurveyFactory.create_survey()``
            to override config values (including verbose: bool = True)
        Returns
        -------
        Survey
            Loaded Survey instance.

        Examples
        --------
        Load LSST year 5 survey:

        >>> lsst_yr5 = Survey.load('lsst', release='yr5')

        Load with uniform survey variant:

        >>> lsst_yr4 = Survey.load('lsst', release='yr4', uniform_survey=True)
        >>> lsst_yr4_uniform = SurveyFactory._cached_surveys['lsst_yr4_uniform']
        """
        return SurveyFactory.create_survey(survey, release, config_file, **kwargs)

    def get_maglim(self, band: str, pixel: int = None) -> float:
        """
        Get magnitude limit for a specific band.

        Parameters
        ----------
        band : str
            Band identifier (e.g., 'g', 'r').
        pixel : int, optional
            HEALPix pixel index. If None, returns entire map.

        Returns
        -------
        float or np.ndarray
            Magnitude limit(s). If pixel is None, returns full HEALPix map array.

        Raises
        ------
        ValueError
            If the specified band is not available in the survey.
        """
        if band not in self.maglim_maps:
            raise ValueError(f"Band '{band}' not available. Available: {self.bands}")

        if pixel is None:
            return self.maglim_maps[band]
        return self.maglim_maps[band][pixel]

    def get_extinction(self, band: str, pixel: int = None) -> float:
        """
        Get extinction (A_band) for a specific band.

        Parameters
        ----------
        band : str
            Band identifier (e.g., 'g', 'r').
        pixel : int, optional
            HEALPix pixel index. If None, returns extinction map for entire sky.

        Returns
        -------
        float or np.ndarray
            Extinction in magnitudes. If pixel is None, returns full HEALPix map array.

        Raises
        ------
        ValueError
            If the specified band is not available or E(B-V) map not loaded.
        """
        if band not in self.coeff_extinc:
            raise ValueError(f"Band '{band}' not available. Available: {self.bands}")

        if self.ebv_map is None:
            raise ValueError("E(B-V) map not loaded")

        extinction = self.coeff_extinc[band] * self.ebv_map

        if pixel is None:
            return extinction
        return extinction[pixel]

    def get_photo_error(
        self, band: str, magnitude: float, maglim: float, **kwargs
    ) -> float:
        """
        Get photometric error estimate.

        Parameters
        ----------
        band : str
            Band identifier (e.g., 'g', 'r'). Used to get band-specific sys_error.
        magnitude : float or np.ndarray
            True apparent magnitude(s).
        maglim : float or np.ndarray
            Magnitude limit(s) at the position(s).
        **kwargs
            Additional keyword arguments:

            delta_saturation : float, optional
                Magnitude difference for saturation threshold in the initial error function.
                Default is -10.4.

        Returns
        -------
        float or np.ndarray
            Total photometric error (including systematic component) in magnitudes.

        Raises
        ------
        ValueError
            If photo error model is not loaded.
        """
        if self.log_photo_error is None:
            raise ValueError("Photo error model not loaded")

        delta_saturation = kwargs.get("delta_saturation", self.delta_saturation)

        # Calculate delta_mag
        delta_mag = magnitude - maglim

        # Get statistical error (same function for all bands)
        mag_err_stat = 10 ** (
            np.where(
                ((delta_mag) <= delta_saturation)
                & (magnitude >= self.saturation[band]),
                self.log_photo_error(delta_saturation),
                self.log_photo_error(delta_mag),
            )
        )
        mag_err_stat = np.where(
            magnitude < self.saturation[band],
            10 ** self.log_photo_error(delta_saturation - 1),
            mag_err_stat,
        )  # saturation at the bright end

        # Add band-specific systematic error in quadrature
        sys_err = self.sys_error.get(band, 0.005)
        total_error = np.sqrt(mag_err_stat**2 + sys_err**2)

        return total_error

    def get_completeness(
        self, band: str, magnitude: float, maglim: float, **kwargs
    ) -> float:
        """
        Get detection completeness (combined detection and classification efficiency).

        This is a convenience wrapper around :meth:`get_efficiency` with
        ``type="completeness"``.

        Parameters
        ----------
        band : str
            Band identifier (e.g., 'g', 'r').
        magnitude : float or np.ndarray
            True apparent magnitude(s).
        maglim : float or np.ndarray
            Magnitude limit(s) at the position(s).
        **kwargs
            Additional keyword arguments passed to :meth:`get_efficiency`.

        Returns
        -------
        float or np.ndarray
            Detection efficiency in range [0, 1].

        See Also
        --------
        get_efficiency : Full documentation of efficiency calculation.
        get_detection_efficiency : Detection-only efficiency.
        get_classification_efficiency : Classification-only efficiency.
        """
        return self.get_efficiency(
            band, magnitude, maglim, type="completeness", **kwargs
        )

    def get_detection_efficiency(
        self, band: str, magnitude: float, maglim: float, **kwargs
    ) -> float:
        """
        Get detection-only efficiency (ignoring classification).

        This is a convenience wrapper around :meth:`get_efficiency` with
        ``type="detection_efficiency"``.

        Parameters
        ----------
        band : str
            Band identifier (e.g., 'g', 'r').
        magnitude : float or np.ndarray
            True apparent magnitude(s).
        maglim : float or np.ndarray
            Magnitude limit(s) at the position(s).
        **kwargs
            Additional keyword arguments passed to :meth:`get_efficiency`.

        Returns
        -------
        float or np.ndarray
            Detection efficiency in range [0, 1].

        See Also
        --------
        get_efficiency : Full documentation of efficiency calculation.
        get_completeness : Combined detection and classification efficiency.
        get_classification_efficiency : Classification-only efficiency.
        """
        return self.get_efficiency(
            band, magnitude, maglim, type="detection_efficiency", **kwargs
        )

    def get_classification_efficiency(
        self, band: str, magnitude: float, maglim: float, **kwargs
    ) -> float:
        """
        Get classification-only efficiency (star/galaxy separation).

        This is a convenience wrapper around :meth:`get_efficiency` with
        ``type="classification_efficiency"``.

        Parameters
        ----------
        band : str
            Band identifier (e.g., 'g', 'r').
        magnitude : float or np.ndarray
            True apparent magnitude(s).
        maglim : float or np.ndarray
            Magnitude limit(s) at the position(s).
        **kwargs
            Additional keyword arguments passed to :meth:`get_efficiency`.

        Returns
        -------
        float or np.ndarray
            Classification efficiency in range [0, 1].

        See Also
        --------
        get_efficiency : Full documentation of efficiency calculation.
        get_completeness : Combined detection and classification efficiency.
        get_detection_efficiency : Detection-only efficiency.
        """
        return self.get_efficiency(
            band, magnitude, maglim, type="classification_efficiency", **kwargs
        )

    def get_efficiency(
        self, band: str, magnitude: float, maglim: float, type="completeness", **kwargs
    ) -> float:
        """
        Get efficiency function values.

        Parameters
        ----------
        band : str
            Band identifier (e.g., 'g', 'r'). Currently unused since same function
            is used for all bands.
        magnitude : float or np.ndarray
            True apparent magnitude(s).
        maglim : float or np.ndarray
            Magnitude limit(s) at the position(s).

        **kwargs
            Additional keyword arguments:
            type : str, optional
                Type of efficiency function to use. Options are "completeness", "detection_efficiency", or "classification_efficiency".
                Default is "completeness".
            delta_saturation : float, optional
                Magnitude difference for saturation threshold in the initial completeness function.
                Default is -10.4.

        Returns
        -------
        float or np.ndarray
            Detection efficiency in range [0, 1]. Value of 0 means no detection,
            1 means certain detection.

        Raises
        ------
        ValueError
            If the efficiency function for the requested ``type`` is not loaded,
            or if ``type`` is not one of the recognized values.

        Notes
        -----
        The efficiency is set to 0 in the following cases:

        - Source is brighter than saturation limit (``magnitude < saturation``)
        - Survey does not cover the position (``maglim < saturation`` or ``maglim`` is NaN)
        - Source is too faint for reliable detection

        Examples
        --------
        Get completeness for a single source:

        >>> maglim = survey.get_maglim('r', pixel=10000)
        >>> eff = survey.get_efficiency('r', 24.0, maglim)

        Get detection efficiency for an array of magnitudes:

        >>> mags = np.linspace(18, 28, 100)
        >>> eff = survey.get_efficiency('r', mags, maglim, type="detection_efficiency")

        See Also
        --------
        get_completeness : Convenience wrapper for combined efficiency.
        get_detection_efficiency : Convenience wrapper for detection-only.
        get_classification_efficiency : Convenience wrapper for classification-only.
        """
        if type == "completeness":
            func = self.completeness
        elif type == "detection_efficiency":
            func = self.efficiency_detection
        elif type == "classification_efficiency":
            func = self.efficiency_classification
        else:
            raise ValueError(f"Efficiency type '{type}' not recognized.")
        if func is None:
            raise ValueError(f"Efficiency function for type '{type}' not loaded")

        delta_saturation = kwargs.get("delta_saturation", self.delta_saturation)

        # Calculate delta_mag
        delta_mag = magnitude - maglim

        # Apply saturation condition: 1 padding for objects fainter than saturation but equivalent mag brighter than saturation0
        compl = np.where(
            (magnitude > self.saturation[band]) & (delta_mag <= delta_saturation),
            1.0,
            func(delta_mag),
        )  # 1 padded

        compl = np.where(
            magnitude < self.saturation[band], 0.0, compl
        )  # saturation at the bright end

        compl = np.where(
            (maglim < self.saturation[band]) | np.isnan(maglim), 0.0, compl
        )  # not observed if the area is not covered

        return compl


class SurveyFactory:
    """
    Factory for creating and caching Survey instances.

    This class handles:

    - Loading survey configurations from YAML/JSON files
    - Caching loaded surveys to avoid redundant file I/O
    - Validating survey configurations
    - Populating Survey objects with all necessary data

    Attributes
    ----------
    _cached_surveys : dict
        Cache of previously loaded Survey objects, keyed by "survey_release".
    """

    _cached_surveys = {}

    @classmethod
    def create_survey(
        cls,
        survey: str,
        release: Optional[str] = None,
        config_file: Optional[dict] = None,
        **kwargs,
    ) -> Survey:
        """
        Create or retrieve a cached Survey instance.

        Parameters
        ----------
        survey : str
            Survey name (e.g., 'lsst').
        release : str, optional
            Survey release/data version (e.g., 'yr1', 'yr5').
            If None, loads the base survey configuration.
        config_file : dict, optional
            Custom configuration dictionary to use instead of loading from file.
            If provided, bypasses the standard config file loading.
        **kwargs
            uniform_survey : bool, optional
                If True, also creates and caches a uniform version of the survey
                with constant magnitude limits. The uniform survey can be accessed
                via ``SurveyFactory._cached_surveys['{survey}_{release}_uniform']``.
                Default is False.
            verbose : bool, optional
                Whether to print progress messages. Default is True.
            Additional keyword arguments override config values

        Returns
        -------
        Survey
            Fully loaded Survey instance (from cache if available).

        Notes
        -----
        - First call: Loads all data files and caches the result.
        - Subsequent calls: Returns cached instance (much faster).
        - Use ``clear_cache()`` to force reloading from files.
        - When ``uniform_survey=True``, the original survey is returned, but a
          uniform variant is also cached for later retrieval.

        See Also
        --------
        _save_uniform_survey : Creates the uniform survey variant.
        clear_cache : Clears cached surveys.
        list_cached_surveys : Lists all cached survey keys.
        """
        verbose = kwargs.get("verbose", True)
        uniform_survey = kwargs.pop("uniform_survey", False)

        # Generate unique cache key
        cache_key = f"{survey}_{release}" if release else survey

        # Return cached survey if available
        if cache_key in cls._cached_surveys:
            if verbose:
                print(f"✓ Using cached survey data for '{cache_key}'")
        else:
            if verbose:
                print(f"Loading survey data for '{cache_key}'...")

            # Load or use provided configuration
            if config_file:
                config = copy.deepcopy(config_file)
            else:
                config = cls._load_config(survey, release, **kwargs)

            # Apply any custom overrides
            config.update(**kwargs)

            # Create empty survey object and populate with data
            survey_obj = Survey(name=survey, release=release)
            cls._load_survey_data(survey_obj, config, **kwargs)

            # Cache for future use
            cls._cached_surveys[cache_key] = survey_obj
            if verbose:
                print(f"✓ Survey '{cache_key}' loaded and cached successfully")

        if uniform_survey:
            cls._save_uniform_survey(cls._cached_surveys[cache_key], **kwargs)

        return cls._cached_surveys[cache_key]

    @classmethod
    def clear_cache(
        cls, survey: Optional[str] = None, release: Optional[str] = None, **kwargs
    ):
        """
        Clear cached survey data to free memory or force reload.

        Parameters
        ----------
        survey : str, optional
            Survey name to clear. If None, clears all cached surveys.
        release : str, optional
            Survey release to clear. Only used if survey is specified.
        **kwargs
            verbose : bool, optional
                Whether to print progress messages. Default is True.

        Examples
        --------
        Clear all cached surveys:

        >>> SurveyFactory.clear_cache()

        Clear all LSST surveys:

        >>> SurveyFactory.clear_cache('lsst')

        Clear only LSST year 1:

        >>> SurveyFactory.clear_cache('lsst', 'yr1')
        """
        verbose = kwargs.get("verbose", True)
        if survey is None:
            # Clear entire cache
            num_cleared = len(cls._cached_surveys)
            cls._cached_surveys.clear()
            if verbose:
                print(f"✓ Cleared {num_cleared} cached survey(s)")
        else:
            # Clear specific survey(s)
            if release is None:
                keys_to_clear = [
                    key for key in cls._cached_surveys if key.startswith(f"{survey}_")
                ] + ([survey] if survey in cls._cached_surveys else [])
                for key in keys_to_clear:
                    del cls._cached_surveys[key]
                    if verbose:
                        print(f"✓ Cleared cached survey '{key}'")
            else:
                cache_key = f"{survey}_{release}"
                if cache_key in cls._cached_surveys:
                    del cls._cached_surveys[cache_key]
                    if verbose:
                        print(f"✓ Cleared cached survey '{cache_key}'")
                else:
                    if verbose:
                        print(f"⚠ Survey '{cache_key}' not found in cache")

    @classmethod
    def list_cached_surveys(cls) -> list:
        """
        Get list of currently cached survey names.

        Returns
        -------
        list of str
            Names of cached surveys in format "survey_release" or "survey".
        """
        return list(cls._cached_surveys.keys())

    @classmethod
    def _save_uniform_survey(cls, survey: Survey, **kwargs):
        """
        Create and cache a uniform version of the survey with constant magnitude limits.

        This method creates a copy of the survey where the magnitude limit maps
        are replaced with a single uniform value (the median of valid pixels).
        This is useful for studying selection effects independent of survey depth
        variations.

        Parameters
        ----------
        survey : Survey
            Survey object to convert to uniform.
        **kwargs
            Additional keyword arguments:

            verbose : bool, optional
                Whether to print progress messages. Default is True.
            uniform_survey_mask_type : list of str, optional
                Mask types to apply when computing median. Options are 'nan'
                (exclude NaN pixels) and 'dust' (exclude high-extinction regions).
                Default is ['nan', 'dust'].

        See Also
        --------
        _get_median_value : Computes the median magnitude limit.
        """
        verbose = kwargs.get("verbose", True)
        uniform_key = (
            f"{survey.name}_{survey.release}_uniform"
            if survey.release
            else f"{survey.name}_uniform"
        )

        if uniform_key in cls._cached_surveys:
            if verbose:
                print(f"✓ Uniform survey '{uniform_key}' already cached")
            return

        if verbose:
            print(f"Creating uniform survey '{uniform_key}'...")

        # Make a deep copy to avoid modifying the original
        survey_uniform = copy.deepcopy(survey)

        # Set uniform magnitude limits for each band
        for band, maglim_map in survey_uniform.maglim_maps.items():
            uniform_limit = cls._get_median_value(
                maglim_map, survey=survey_uniform, **kwargs
            )
            survey_uniform.maglim_maps[band] = np.where(
                np.isnan(maglim_map), np.nan, uniform_limit
            )  # Set uniform value, keep NaNs - not observed areas
            if verbose:
                print(f"  {band}-band: uniform maglim = {uniform_limit:.2f} mag")

        cls._cached_surveys[uniform_key] = survey_uniform
        if verbose:
            print(f"✓ Uniform survey cached as '{uniform_key}'")

    @classmethod
    def _get_median_value(cls, maglim_map: np.ndarray, **kwargs) -> float:
        """
        Calculate median magnitude limit, optionally masking invalid regions.

        Parameters
        ----------
        maglim_map : np.ndarray
            HEALPix map of magnitude limits.
        **kwargs
            Additional keyword arguments:

            survey : Survey, optional
                Survey object (required if 'dust' is in mask_type).
            verbose : bool, optional
                Whether to print progress messages. Default is True.
            uniform_survey_mask_type : list of str or str, optional
                Types of masking to apply before computing the median:

                - 'nan' : Exclude NaN pixels (always recommended).
                - 'dust' : Exclude high-extinction regions (E(B-V) >= 0.2).

                Default is ['nan', 'dust'].

        Returns
        -------
        float
            Median magnitude limit of valid pixels.

        Raises
        ------
        ValueError
            If 'dust' mask is requested but survey object or E(B-V) map is unavailable.
        """
        verbose = kwargs.get("verbose", True)
        mask_type = kwargs.get("uniform_survey_mask_type", ["nan", "dust"])

        if isinstance(mask_type, str):
            mask_type = [mask_type]

        # Build valid pixel mask
        valid_pixels = np.ones_like(maglim_map, dtype=bool)

        for mtype in mask_type:
            if mtype == "nan":
                valid_pixels &= ~np.isnan(maglim_map)

            elif mtype == "dust":
                survey = kwargs.get("survey")
                if survey is None:
                    raise ValueError("Survey object required for 'dust' mask type.")
                if not hasattr(survey, "ebv_map") or survey.ebv_map is None:
                    raise ValueError("E(B-V) map not available for 'dust' mask type.")

                # Resample E(B-V) map if needed
                nside_maglim = hp.get_nside(maglim_map)
                nside_ebv = hp.get_nside(survey.ebv_map)
                if nside_maglim != nside_ebv:
                    ebv_map = hp.ud_grade(survey.ebv_map, nside_out=nside_maglim)
                else:
                    ebv_map = survey.ebv_map

                valid_pixels &= ebv_map < 0.2

        if verbose:
            n_valid = np.sum(valid_pixels)
            n_total = len(maglim_map)
            print(
                f"  Computing median from {n_valid}/{n_total} pixels (mask: {mask_type})"
            )

        return np.nanmedian(maglim_map[valid_pixels])

    @classmethod
    def _load_config(cls, survey: str, release: Optional[str] = None, **kwargs) -> dict:
        """
        Load survey configuration from YAML or JSON file.

        Parameters
        ----------
        survey : str
            Survey name.
        release : str, optional
            Survey release identifier.
        **kwargs
            verbose : bool, optional
                Whether to print progress messages. Default is True.

        Returns
        -------
        dict
            Configuration dictionary containing survey_files and survey_properties.

        Raises
        ------
        FileNotFoundError
            If no configuration file is found for the survey.
        ValueError
            If the loaded config doesn't match the requested survey/release.
        """
        verbose = kwargs.get("verbose", True)
        # Construct path to config directory
        current_dir = os.path.dirname(os.path.abspath(__file__))
        config_dir = os.path.join(current_dir, "..", "config/surveys/")

        # Build config filename
        config_name = f"{survey}_{release}" if release else survey

        # Try YAML first, then JSON
        yaml_path = os.path.join(config_dir, f"{config_name}.yaml")
        json_path = os.path.join(config_dir, f"{config_name}.json")

        config_path = None
        if os.path.exists(yaml_path):
            config_path = yaml_path
        elif os.path.exists(json_path):
            config_path = json_path
        else:
            raise FileNotFoundError(
                f"Configuration file not found for '{config_name}'\n"
                f"Searched in: {config_dir}\n"
                f"Expected: {config_name}.yaml or {config_name}.json"
            )

        # Load configuration file
        if verbose:
            print(f"  Loading config from: {os.path.basename(config_path)}")
        if config_path.endswith(".yaml"):
            import yaml

            with open(config_path, "r") as file:
                config_data = yaml.safe_load(file)
        else:
            import json

            with open(config_path, "r") as file:
                config_data = json.load(file)

        # Validate configuration matches request
        config_survey = config_data.get("name")
        config_release = config_data.get("release", None)

        if survey != config_survey:
            raise ValueError(
                f"Survey name mismatch:\n"
                f"  Requested: '{survey}'\n"
                f"  In config: '{config_survey}'"
            )

        if release != config_release:
            raise ValueError(
                f"Release mismatch for survey '{survey}':\n"
                f"  Requested: '{release}'\n"
                f"  In config: '{config_release}'"
            )

        return config_data

    @classmethod
    def _load_survey_data(cls, survey: Survey, config: dict, **kwargs):
        """
        Load all survey data files into the Survey object.

        This method performs the following operations:

        1. Determines available bands from config
        2. Determines the data directory from config
        3. Loads HEALPix maps for each band (magnitude limits)
        4. Loads band-independent maps (extinction, coverage)
        5. Loads efficiency and photometric error functions per band
        6. Sets survey properties per band (extinction, saturation, sys_error)

        Parameters
        ----------
        survey : Survey
            Empty Survey object to populate.
        config : dict
            Configuration dictionary with 'survey_files' and 'survey_properties' keys.
        **kwargs
            verbose : bool, optional
                Whether to print progress messages. Default is True.
        """
        verbose = kwargs.get("verbose", True)
        if verbose:
            print("\n" + "=" * 70)
            print("LOADING SURVEY DATA FILES")
            print("=" * 70)

        # Get file configuration and determine data path
        survey_config = config.get("survey_files", {})
        data_path_survey = survey_config.get("file_path", "")

        # Find data directory where might be stored additional files
        current_dir = os.path.dirname(os.path.abspath(__file__))
        survey_name = (
            f"{survey.name}_{survey.release}" if survey.release else survey.name
        )
        data_path_others = os.path.join(current_dir, "..", "data", "others")
        # Use default path if not specified
        if not data_path_survey:
            data_path_survey = os.path.join(
                current_dir, "..", "data", "surveys", survey_name
            )

        if verbose:
            print(f"Survey data directory: {data_path_survey}\n")
            print(f"Fallback directory for shared data files: {data_path_others}\n")

        # Determine available bands from config
        # Look for maglim_map_* entries to discover bands
        available_bands = []
        for key in survey_config.keys():
            if key.startswith("maglim_map_"):
                band = key.replace("maglim_map_", "")
                if band and band not in available_bands:
                    available_bands.append(band)

        # Also check survey_properties for bands
        props = config.get("survey_properties", {})
        bands_from_props = props.get("bands", [])
        for band in bands_from_props:
            if band not in available_bands:
                available_bands.append(band)

        # Default to ['g', 'r'] if no bands specified
        if not available_bands:
            available_bands = ["g", "r"]

        survey.bands = sorted(available_bands)
        if verbose:
            print(f"Available bands: {', '.join(survey.bands)}\n")

        # Load survey properties per band
        if verbose:
            print("\nLoading survey properties...")

        # Extinction coefficients per band
        default_extinc = {"g": 3.303, "r": 2.285, "i": 1.682, "z": 1.322, "y": 1.087}
        for band in survey.bands:
            coeff_key = f"coeff_extinc_{band}"
            survey.coeff_extinc[band] = props.get(
                coeff_key, default_extinc.get(band, None)
            )

        # Saturation per band
        default_saturation = props.get("saturation", 16.0)
        for band in survey.bands:
            sat_key = f"saturation_{band}"
            survey.saturation[band] = props.get(sat_key, default_saturation)

        # Systematic error per band
        default_sys_error = props.get("sys_error", 0.005)
        for band in survey.bands:
            sys_key = f"sys_error_{band}"
            survey.sys_error[band] = props.get(sys_key, default_sys_error)

        # Load band-specific magnitude limit maps
        if verbose:
            print("Loading magnitude limit maps...")
        for band in survey.bands:
            attr_name = f"maglim_map_{band}"
            filename = survey_config.get(attr_name)
            extension = os.path.splitext(filename)[1] if filename else ""

            if filename is None:
                if verbose:
                    print(
                        f"  ⚠ Warning: '{attr_name}' not specified in config (skipping {band}-band)"
                    )
                continue

            full_path = cls._find_file(filename, data_path_survey, data_path_others)

            if full_path is None:
                if verbose:
                    print(
                        f"  ✗ {band}-band magnitude limit: File not found - {filename}"
                    )
                continue

            try:
                if extension.lower() == ".hsp":
                    survey.maglim_maps[band] = cls.open_map_sparse(
                        full_path, map_min=survey.saturation[band]
                    )
                else:
                    survey.maglim_maps[band] = hp.read_map(full_path, verbose=False)
                if verbose:
                    print(f"  ✓ Success for {band}-band magnitude limit")
            except Exception as e:
                if verbose:
                    print(f"    ✗ Failed to load {band}-band maglim map: {e}")
                survey.maglim_maps[band] = None

        # Load completeness function (same for all bands)
        if verbose:
          print("\nLoading completeness/efficiency function...")
        survey.completeness_band = survey_config.get("completeness_band", "r") # default to r band
        survey.delta_saturation = props.get("delta_saturation", -10.4)
        if "completeness" in survey_config:
            # Use default saturation if not band-specific
            cls._load_file(
                survey,
                survey_config,
                "completeness",
                "Completeness/efficiency function",
                lambda f: cls.set_completeness(
                    f, delta_saturation=survey.delta_saturation
                ),
                data_path_survey,
                data_path_others,
                **kwargs,
            )

            # Load detection and classification efficiencies separately, which can be usefull to study selection effects
            # eg. perfect star-galaxy separation vs imperfect one
            # Try to load detection efficiency
            try:
                cls._load_file(
                    survey,
                    survey_config,
                    "efficiency_detection",
                    "Detection efficiency function",
                    lambda f: cls.set_completeness(
                        f,
                        delta_saturation=survey.delta_saturation,
                        selection="detected",
                    ),
                    data_path_survey,
                    data_path_others,
                    filename=survey_config.get(
                        "completeness"
                    ),  # use same file as completeness
                    **kwargs,
                )
            except:
                if verbose:
                    print("No detection efficiency file found, skipping.")

            # Try to load classification efficiency
            try:
                cls._load_file(
                    survey,
                    survey_config,
                    "efficiency_classification",
                    "Classification efficiency function",
                    lambda f: cls.set_completeness(
                        f,
                        delta_saturation=survey.delta_saturation,
                        selection="classified",
                    ),
                    data_path_survey,
                    data_path_others,
                    filename=survey_config.get(
                        "completeness"
                    ),  # use same file as completeness
                    **kwargs,
                )
            except:
                if verbose:
                    print("No classification efficiency file found, skipping.")

        # Load photometric error model (same for all bands)
        if verbose:
            print("\nLoading photometric error model...")
        if "log_photo_error" in survey_config:
            # Use default saturation if not band-specific
            cls._load_file(
                survey,
                survey_config,
                "log_photo_error",
                "Photometric error model",
                lambda f: cls.set_photo_error(
                    f, delta_saturation=survey.delta_saturation
                ),
                data_path_survey,
                data_path_others,
                **kwargs,
            )

        # Load band-independent maps
        if verbose:
            print("\nLoading band-independent maps...")
        band_independent = {
            "ebv_map": ("E(B-V) extinction map", lambda f: hp.read_map(f)),
            "coverage": ("Survey coverage map", lambda f: hp.read_map(f)),
        }
        for attr_name, (description, loader_func) in band_independent.items():
            cls._load_file(
                survey,
                survey_config,
                attr_name,
                description,
                loader_func,
                data_path_survey,
                data_path_others,
                **kwargs,
            )

        if survey.coverage is None:
            if verbose:
                print("\nBuilding coverage map from magnitude limit maps...")
            cls._build_coverage_map(survey, **kwargs)

        # Print summary
        if verbose:
            print("\nSurvey properties summary:")
            for band in survey.bands:
                print(f"  {band}-band:")
                print(f"    Extinction coefficient: {survey.coeff_extinc[band]:.3f}")
                print(f"    Saturation limit: {survey.saturation[band]:.1f} mag")
                print(f"    Systematic error: {survey.sys_error[band]:.4f} mag")

            print("\n" + "=" * 70)
            print("SURVEY DATA LOADED SUCCESSFULLY")
            print("=" * 70 + "\n")

    @classmethod
    def _find_file(
        cls, filename: str, data_path_survey: str, data_path_others: str
    ) -> Optional[str]:
        """
        Find a file by checking data_path_survey first, then data_path_others.

        Parameters
        ----------
        filename : str
            Name of the file to find.
        data_path_survey : str
            Primary data directory (survey-specific).
        data_path_others : str
            Fallback data directory (shared files).

        Returns
        -------
        str or None
            Full path to the file if found, None otherwise.
        """
        # Check survey-specific directory first
        full_path = os.path.join(data_path_survey, filename)
        if os.path.exists(full_path):
            return full_path

        # Check shared directory as fallback
        full_path = os.path.join(data_path_others, filename)
        if os.path.exists(full_path):
            return full_path

        return None

    @classmethod
    def _load_file(
        cls,
        survey: Survey,
        config: dict,
        attr_name: str,
        description: str,
        loader_func: Callable,
        data_path_survey: str,
        data_path_others: str = None,
        filename: str = None,
        **kwargs,
    ):
        """
        Load a single data file and attach it to the survey object.

        Parameters
        ----------
        survey : Survey
            Survey object to populate.
        config : dict
            Configuration dictionary containing filenames.
        attr_name : str
            Attribute name to set on the survey object (e.g., 'ebv_map').
        description : str
            Human-readable description for logging.
        loader_func : callable
            Function to load the file (e.g., hp.read_map).
        data_path_survey : str
            Primary directory containing the data files.
        data_path_others : str, optional
            Fallback directory for shared data files.
        filename : str, optional
            Override filename from config.
        **kwargs
            verbose : bool, optional
                Whether to print progress messages. Default is True.
        """
        verbose = kwargs.get("verbose", True)
        # Get filename from config
        if filename is None:
            filename = config.get(attr_name)

        if filename is None:
            if verbose:
                print(f"  ⚠ Warning: '{attr_name}' not specified in config (skipping)")
            return

        # Find file in data paths
        full_path = cls._find_file(
            filename, data_path_survey, data_path_others or data_path_survey
        )

        if full_path is None:
            if verbose:
                print(f"  ✗ {description}: File not found - {filename}")
            return

        # Load the file
        try:
            if verbose:
                print(f"  Loading {description}...")
                print(f"    File: {filename}")
            data = loader_func(full_path)
            setattr(survey, attr_name, data)
            if verbose:
                print(f"    ✓ Success")
        except Exception as e:
            if verbose:
                print(f"    ✗ Failed to load {attr_name}: {e}")
            setattr(survey, attr_name, None)

    @staticmethod
    def open_map_sparse(file_path, **kwargs):
        """
        Open a sparse HEALPix map and return it as a standard HEALPix map.

        Parameters
        ----------
        file_path : str
            Full path to the sparse HEALPix file (.hsp).
        **kwargs
            Additional keyword arguments (currently unused).

        Returns
        -------
        np.ndarray
            Standard HEALPix map with NaN for unseen pixels.
        """
        hsp_map = hsp.HealSparseMap.read(file_path)
        nside_sparse = hsp_map.nside_sparse
        nest = False  # by default set to True in hsp, but to False in healpy
        map_hpx = hsp_map.generate_healpix_map(nside=nside_sparse, nest=nest)
        map_hpx = np.where(map_hpx == hp.UNSEEN, np.nan, map_hpx)

        map_min, map_max = kwargs.get("map_min", None), kwargs.get("map_max", None)
        if map_min is not None:
            map_hpx = np.where(map_hpx < map_min, np.nan, map_hpx)
        if map_max is not None:
            map_hpx = np.where(map_hpx > map_max, np.nan, map_hpx)

        return map_hpx

    @staticmethod
    def set_completeness(filename, delta_saturation=-10.4, selection="both"):
        """
        Load and interpolate completeness/efficiency function from file.

        Parameters
        ----------
        filename : str
            Path to CSV file with columns:

            - 'delta_mag' : Magnitude difference from limit (mag - maglim)
            - 'eff_star' : Detection efficiency only
            - 'classifiction_eff' : Classification efficiency only
            - 'classification_detection_eff' : Combined efficiency

        delta_saturation : float, optional
            Magnitude difference threshold for saturation. Default is -10.4.
        selection : str, optional
            Which efficiency to use:

            - 'detected' : Detection efficiency only (column 'eff_star')
            - 'classified' : Classification efficiency only (column 'classifiction_eff')
            - 'both' : Combined detection and classification (column 'classification_detection_eff')

            Default is 'both'.

        Returns
        -------
        callable
            Interpolation function f(delta_mag) -> efficiency [0, 1].

        Raises
        ------
        ValueError
            If selection is not one of 'detected', 'classified', or 'both'.

        Notes
        -----
        - Bright stars (delta_mag <= delta_saturation): Efficiency forced to 0.
        - Faint stars (beyond data): Returns efficiency = 0.0.
        - The saturation parameter is automatically passed from the survey object.
        """
        # Load photometric error data
        data = np.genfromtxt(filename, delimiter=",", names=True)
        delta_mags = data["delta_mag"]

        # Select efficiency column based on user choice
        if selection == "detected":
            efficiencies = data["detection_eff"]
        elif selection == "classified":
            efficiencies = data["classifiction_eff"]
        elif selection == "both":
            efficiencies = data["classification_detection_eff"]
        else:
            raise ValueError(
                f"Invalid selection '{selection}'. Must be 'detected', 'classified', or 'both'"
            )

        # Extend efficiency to bright end (force to zero at saturation)
        if delta_mags.min() > delta_saturation:
            delta_mags = np.insert(delta_mags, 0, delta_saturation)
            efficiencies = np.insert(efficiencies, 0, 0.0)
        elif delta_mags.min() == delta_saturation:
            # Ensure efficiency is zero at magnitude very near saturation
            delta_mags = np.insert(delta_mags, 0, delta_saturation-1e-5)
            efficiencies = np.insert(efficiencies, 0, 0.0)
        else:
            # Ensure efficiency is zero at saturation
            efficiencies[delta_mags <= delta_saturation] = 0.0

        # Extend efficiency to faint end (force to zero)
        delta_mags = np.append(delta_mags, delta_mags[-1] + 1)
        efficiencies = np.append(efficiencies, 0.0)

        # Create interpolation function
        interpolator = scipy.interpolate.interp1d(
            delta_mags,
            efficiencies,
            bounds_error=False,
            fill_value=0.0,  # Return 0 for very faint/bright stars
        )

        return interpolator

    @staticmethod
    def set_photo_error(filename, delta_saturation=-10.4):
        """
        Load photometric error model from file.

        This creates an interpolation function that estimates photometric uncertainty
        as a function of the magnitude difference from the survey limit (delta_mag).

        Parameters
        ----------
        filename : str
            Path to CSV file with columns:

            - 'delta_mag' : Magnitude difference from limit (mag - maglim)
            - 'log_mag_err' : Logarithm (base 10) of magnitude error

        delta_saturation : float, optional
            Bright magnitude difference threshold. Used to determine the bright-end
            extension point. Default is -10.4.

        Returns
        -------
        callable
            Interpolation function f(delta_mag) -> log10(magnitude_error).

        Notes
        -----
        - Bright stars (delta_mag < delta_saturation): Extended with constant error.
        - Faint stars (beyond data): Returns log10(error) = 1.0 (error = 10 mag).
        - The saturation parameter is automatically passed from the survey object.
        """
        # Load photometric error data
        data = np.genfromtxt(filename, delimiter=",", names=True)
        delta_mags = data["delta_mag"]
        log_errors = data["log_mag_err"]

        # Extend to bright end (keep constant for very bright stars)
        if delta_mags.min() > delta_saturation:
            delta_mags = np.insert(delta_mags, 0, delta_saturation)
            log_errors = np.insert(log_errors, 0, log_errors[0])

        # Create interpolation function
        interpolator = scipy.interpolate.interp1d(
            delta_mags,
            log_errors,
            bounds_error=False,
            fill_value=1.0,  # Return log10(error)=1 (10 mag error) for very faint stars
        )

        return interpolator

    @classmethod
    def _build_coverage_map(cls, survey: Survey, **kwargs):
        """
        Build coverage map from magnitude limit maps.

        Creates a combined coverage map by checking all available magnitude limit
        maps and determining which pixels are covered in all bands (logical AND).

        Parameters
        ----------
        survey : Survey
            Survey object with loaded magnitude limit maps.
        **kwargs
            verbose : bool, optional
                Whether to print progress messages. Default is True.

        Returns
        -------
        np.ndarray or None
            HEALPix coverage map (1=observed, 0=not observed), or None if no
            magnitude limit maps are available.
        """
        verbose = kwargs.get("verbose", True)
        nside = None
        coverage_map = None

        # Find the minimum nside among all maglim maps
        for band, maglim_map in survey.maglim_maps.items():
            if maglim_map is None:
                continue

            nside_candidate = hp.get_nside(maglim_map)
            if nside is None or nside_candidate < nside:
                nside = nside_candidate

        if nside is None:
            if verbose:
                print(
                    "  ⚠ Warning: No magnitude limit maps found, cannot build coverage"
                )
            return None

        # Initialize coverage map as all True
        npix = hp.nside2npix(nside)
        coverage_map = np.ones(npix, dtype=bool)

        # Accumulate coverage from all bands (logical AND)
        for band, maglim_map in survey.maglim_maps.items():
            if maglim_map is None:
                continue

            # Convert to common nside if needed
            nside_band = hp.get_nside(maglim_map)
            if nside_band != nside:
                maglim_map = hp.ud_grade(maglim_map, nside_out=nside)

            # Update coverage: pixel is covered if maglim is finite and above saturation
            band_coverage = np.isfinite(maglim_map) & (
                maglim_map > survey.saturation[band]
            )
            coverage_map = coverage_map & band_coverage

        # Convert boolean to float (1.0 = covered, 0.0 = not covered)
        survey.coverage = coverage_map.astype(float)
        if verbose:
            print(
                f"  ✓ Built coverage map (nside={nside}, {np.sum(coverage_map)} pixels covered)"
            )

        return survey.coverage
