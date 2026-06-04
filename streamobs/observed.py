#!/usr/bin/env python

import copy
import os

import astropy.coordinates as coord
import astropy.units as u
import gala.coordinates as gc
import healpy as hp
import numpy as np
import pandas as pd

from .model import StreamModel
from .plotting import plot_stream_in_mask
from .surveys import Survey


class StreamInjector:
    """
    Inject observational effects into stream data for a given survey.

    This class handles the core injection logic while keeping survey data separate.
    All survey data is loaded once and cached, making multiple injections efficient.

    Attributes
    ----------
    survey : Survey
        Survey instance containing all survey-specific data and functions.
    mask_cache : dict (class attribute)
        Cache of previously created HEALPix masks to avoid recomputation.
    _last_gc_frame : GreatCircleICRSFrame or None
        The most recently used great circle frame. This allows reusing the same
        sky location across multiple inject() calls when gc_frame='last'.

    Examples
    --------
    Initialize with a survey:

    >>> injector = StreamInjector('lsst', release='dc2')

    Or with a pre-loaded Survey object:

    >>> survey = Survey.load('lsst', release='dc2')
    >>> injector = StreamInjector(survey)
    """

    mask_cache = {}

    def __init__(self, survey, **kwargs):
        """
        Initialize with survey configuration.

        Parameters
        ----------
        survey : str or Survey
            Either a survey name string (e.g., 'lsst') or a pre-loaded Survey instance.
        **kwargs
            Additional keyword arguments passed to Survey.load() if survey is a string.
            Common options include:

            release : str, optional
                Survey release version (e.g., 'dc2', 'yr1', 'yr10').

        Raises
        ------
        ValueError
            If survey is neither a string nor a Survey instance.
        """

        if isinstance(survey, str):
            self.survey = Survey.load(survey=survey, **kwargs)
        elif isinstance(survey, Survey):
            self.survey = survey
        else:
            raise ValueError("survey must be a string or Survey instance.")

        # Instance attribute to store the last used gc_frame
        self._last_gc_frame = None

    def inject(self, data, bands=["r", "g"], **kwargs):
        """
        Add observed quantities from the survey to the given data.

        This method applies observational effects including photometric errors,
        magnitude measurements, and detection flags based on survey properties.

        Parameters
        ----------
        data : str or pd.DataFrame
            Input data as DataFrame or path to the file (CSV or Excel).
            Must contain either (ra, dec) or (phi1, phi2) coordinates, and
            magnitude columns for the specified bands.
        bands : list of str, optional
            List of photometric bands to process. Default is ['r', 'g'].
        **kwargs
            Additional keyword arguments:

            seed : int, optional
                Random seed for reproducibility.
            nside : int, optional
                HEALPix nside parameter. Default is 4096.
            detection_mag_cut : list of str, optional
                Bands to apply SNR detection cut. Default is ['g'].
            save : bool, optional
                Whether to save the output data. Default is False.
            folder : str or Path, optional
                Output folder path if save=True.
            dust_correction : bool, optional
                Whether to apply dust correction to observed magnitudes. Default is True.
            perfect_galstarsep : bool, optional
                If True, also computes a flag assuming perfect star/galaxy separation
                (detection efficiency only, no classification losses). Default is False.
            verbose : bool, optional
                Whether to print progress information. Default is True.

        Returns
        -------
        pd.DataFrame
            DataFrame with the following added columns:

            - mag_<band>_obs : Observed magnitudes for each band
            - magerr_<band> : Photometric errors for each band
            - flag_observed : Boolean flag (True=detected, False=not detected). Includes both detection and classification efficiencies.
            - flag_perfect_galstarsep : Boolean flag assuming perfect star/galaxy separation (detection efficiency only, no classification losses; only if requested with perfect_galstarsep=True)
            - ra, dec : Sky coordinates (if not already present)

        Raises
        ------
        ValueError
            If required columns are missing or bands are not supported.
        """
        verbose = kwargs.get("verbose", True)
        perfect_galstarsep = kwargs.pop("perfect_galstarsep", False)

        # Load data
        data = self._load_data(data)

        # Set the seed for reproducibility
        seed = kwargs.pop("seed", None)
        rng = np.random.default_rng(seed)

        # Complete missing columns (ra/dec, magnitudes)
        data = self.complete_data(data, rng=rng, seed=seed, bands=bands, **kwargs)

        # Get HEALPix pixel indices
        nside = kwargs.pop("nside", 4096)
        pix = hp.ang2pix(nside, data["ra"], data["dec"], lonlat=True)

        # Initialize detection flags (will be updated per band)
        flag_completeness_band = None
        flag_detection_only_band = None

        # Process each band
        for band in bands:
            if band not in ["r", "g"]:
                raise ValueError("Currently only 'r' and 'g' bands are supported.")

            # Get extinction for this band
            nside_ebv = hp.get_nside(self.survey.ebv_map)
            if nside_ebv != nside:
                pix_ebv = hp.ang2pix(nside_ebv, data["ra"], data["dec"], lonlat=True)
            else:
                pix_ebv = pix
            extinction_band = self.survey.get_extinction(band, pixel=pix_ebv)

            # Calculate true apparent magnitudes (including extinction)
            apparent_mag_true = data["mag_" + band] + extinction_band

            # Get magnitude limits
            nside_maglim = hp.get_nside(self.survey.maglim_maps[band])
            if nside_maglim != nside:
                pix_maglim = hp.ang2pix(
                    nside_maglim, data["ra"], data["dec"], lonlat=True
                )
            else:
                pix_maglim = pix

            # Calculate photometric errors
            mag_err = self.survey.get_photo_error(
                band,
                apparent_mag_true,
                self.survey.get_maglim(band, pixel=pix_maglim),
            )

            # Sample measured magnitudes
            mag_obs = self.sample_measured_magnitudes(
                apparent_mag_true,
                mag_err,
                rng=rng,
                seed=seed,
                **kwargs,
            )

            dust_correction = kwargs.get("dust_correction", True)
            if dust_correction:
                if verbose:
                    print(
                        f"Applying dust correction for {band}-band on observed magnitudes."
                    )
                # Correct observed magnitudes for extinction (only for valid detections)
                valid_mask = mag_obs != "BAD_MAG"
                # Create a float array for corrected magnitudes
                mag_obs_corrected = np.empty(len(mag_obs), dtype=object)
                mag_obs_corrected[~valid_mask] = "BAD_MAG"
                mag_obs_corrected[valid_mask] = (
                    mag_obs[valid_mask].astype(float) - extinction_band[valid_mask]
                )
                mag_obs = mag_obs_corrected

            # Add new columns
            new_columns = pd.DataFrame(
                {
                    "mag_" + band + "_obs": mag_obs,
                    "magerr_" + band: mag_err,
                }
            )

            # Reset indices and concatenate
            data = data.reset_index(drop=True)
            new_columns = new_columns.reset_index(drop=True)
            data = pd.concat([data, new_columns], axis=1)

            # Compute detection flag for completeness-band (reference band)
            if band == self.survey.completeness_band:
                flag_completeness_band = self.detect_flag(
                    pix_maglim,
                    mag=apparent_mag_true,
                    band=band,
                    rng=rng,
                    seed=seed,
                    perfect_galstarsep=False,
                    **kwargs,
                )
                if perfect_galstarsep:
                    flag_detection_only_band = self.detect_flag(
                        pix_maglim,
                        mag=apparent_mag_true,
                        band=band,
                        rng=rng,
                        seed=seed,
                        perfect_galstarsep=True,
                        **kwargs,
                    )

        # Apply detection threshold
        if flag_completeness_band is None:
            if self.survey.completeness_band in bands:
                raise ValueError(
                    f"flag_completeness_{self.survey.completeness_band} must be computed for detection in {self.survey.completeness_band} band."
                )
            else:
                raise ValueError(f"Detection flag requires '{self.survey.completeness_band}' band to be in bands.")

        # Build combined detection flags
        # Start with flux validity check (not BAD_MAG)
        flag_valid_flux = data["mag_r_obs"] != "BAD_MAG"
        if "g" in bands:
            flag_valid_flux &= data["mag_g_obs"] != "BAD_MAG"

        # Combine with completeness
        flag_observed = flag_valid_flux & flag_completeness_band if flag_completeness_band is not None else flag_valid_flux 
        if perfect_galstarsep:
            flag_perfect = flag_valid_flux & flag_detection_only_band if flag_detection_only_band is not None else flag_valid_flux

        # Apply SNR cuts
        detection_mag_cut = kwargs.get("detection_mag_cut", ["g"])
        SNR_min = 5.0

        for band in detection_mag_cut:
            if band not in bands:
                if verbose:
                    print(
                        f"Warning: SNR cut requested for {band}-band but it's not in bands list. Skipping."
                    )
                continue
            if verbose:
                print(f"Applying detection cut on {band}-band with SNR >= {SNR_min}")
            SNR = 1.0 / data["magerr_" + band]
            flag_observed &= SNR >= SNR_min
            if perfect_galstarsep:
                flag_perfect &= SNR >= SNR_min

        # Force SNR cut on r-band for perfect gal/star separation (because it's not included in the detection efficiency functions, while it is for the completeness functions)
        if perfect_galstarsep:
            SNR_r = 1.0 / data["magerr_r"]
            flag_perfect &= SNR_r >= SNR_min

        # Store flags in DataFrame
        data["flag_observed"] = flag_observed
        if perfect_galstarsep:
            data["flag_perfect_galstarsep"] = flag_perfect

        # Save if requested
        if kwargs.get("save"):
            self._save_injected_data(data, kwargs.get("folder", None))

        # Return data (do NOT store as instance attribute to avoid conflicts between runs)
        return data

    def _load_data(self, data):
        """
        Load data from file or return the provided DataFrame.

        Parameters
        ----------
        data : str or pd.DataFrame
            Path to the file or pandas DataFrame.

        Returns
        -------
        pd.DataFrame
            Loaded DataFrame.

        Raises
        ------
        ValueError
            If file format is unsupported or data type is invalid.
        """
        if isinstance(data, pd.DataFrame):
            return data
        elif isinstance(data, str):
            extension = data.split(".")[-1]
            if extension == "csv":
                return pd.read_csv(data)
            elif extension in ["xls", "xlsx"]:
                return pd.read_excel(data)
            else:
                raise ValueError(f"Unsupported file format: {extension}")
        else:
            raise ValueError(f"Unsupported file format")

    def complete_data(self, data, bands=["r", "g"], **kwargs):
        """
        Validate and complete the columns needed for injection.

        This function ensures sky coordinates are present and, if requested,
        fills missing photometric magnitude columns using a stream model.
        Specifically:

        - If (ra, dec) are missing but (phi1, phi2) exist, it converts to
          sky coordinates via ``phi_to_radec``.
        - It checks for ``mag_<band>`` for each requested band; any missing
          bands are generated via ``StreamModel.complete_catalog`` when a
          ``stream_config`` is provided.
        - Existing columns are preserved and are not overwritten.

        Parameters
        ----------
        data : pandas.DataFrame
            Input catalog. Must contain either (ra, dec) or (phi1, phi2).
            The DataFrame is copied; the input object is not modified in place.
        bands : list[str], optional
            Photometric bands to ensure (as ``mag_<band>``). Default is
            ``['r', 'g']``.
        **kwargs
            Additional options:

            rng : numpy.random.Generator, optional
                Random number generator instance.
            seed : int, optional
                Random seed used to initialize an RNG if ``rng`` is not given.
            gc_frame : gala.coordinates.GreatCircleICRSFrame or 'last', optional
                Great circle frame to use. If 'last', uses the frame from the
                previous call. If None, generates a new random frame.
            stream_config : dict, optional
                Configuration passed to ``StreamModel``. Required only when at
                least one requested magnitude band is missing. See
                ``StreamModel`` for the expected schema (e.g., track, distance
                modulus, isochrone).

            Keyword arguments forwarded to ``phi_to_radec`` (only used when
            converting from (phi1, phi2)):

            mask_type : list[str] or str, optional
                Mask(s) used when searching a frame (e.g., 'footprint',
                'maglim_r', 'ebv').
            percentile_threshold : float, optional
                Fraction of points that must lie inside the mask when selecting
                a frame. Default defined in ``phi_to_radec``.
            max_iter : int, optional
                Maximum trials when searching a frame. Default defined in
                ``phi_to_radec``.

        Returns
        -------
        pandas.DataFrame
            A copy of the input with the required columns present:
            ``ra``, ``dec``, and one ``mag_<band>`` column per requested band.

        Raises
        ------
        ValueError
            If neither (ra, dec) nor (phi1, phi2) are present; or if one or
            more magnitude bands are missing and ``stream_config`` is not
            provided; or if after completion a required column is still absent.

        Notes
        -----
        - Only missing magnitude columns are synthesized; existing values are
          left untouched.
        - Magnitudes are produced by ``StreamModel.complete_catalog`` and rely
          on the model's configuration (e.g., isochrone and distance modulus).
        - When a new gc_frame is created or used, it is stored in ``self._last_gc_frame``
          for potential reuse via ``gc_frame='last'``.

        Examples
        --------
        Convert from (phi1, phi2) and generate both bands:

        >>> df = pd.DataFrame({'phi1': [-5, 0, 5], 'phi2': [0, 0, 0]})
        >>> cfg = {...}  # stream model config
        >>> out = injector.complete_data(df, bands=['r', 'g'], stream_config=cfg)

        Use already existing ra/dec and fill a missing 'g' band from a model:

        >>> df = pd.DataFrame({'phi1': [-5, 0, 5], 'phi2': [0, 0, 0], 'ra': [10.1, 12.5, 24.7], 'dec': [5.2, 6.2, 7.2],})
        >>> out = injector.complete_data(df, bands=['r', 'g'], stream_config=cfg, seed=42)

        Reuse the same gc_frame for multiple streams:

        >>> data1 = injector.complete_data(df1, seed=42)  # Creates new frame
        >>> data2 = injector.complete_data(df2, gc_frame='last')  # Reuses frame from data1
        """

        required_columns = ["ra", "dec"] + [f"mag_{band}" for band in bands]

        rng = kwargs.pop("rng", None)
        seed = kwargs.pop("seed", None)

        # Make explicit copy to avoid SettingWithCopyWarning
        data = data.copy()

        if not ("ra" in data.columns and "dec" in data.columns):
            if "phi1" not in data.columns or "phi2" not in data.columns:
                raise ValueError(
                    "Input data must contain either (ra, dec) or (phi1, phi2) columns."
                )

            # Handle 'last' keyword for gc_frame
            gc_frame_param = kwargs.get("gc_frame", None)
            if gc_frame_param == "last":
                kwargs["gc_frame"] = self._last_gc_frame

            # Convert coordinates (Phi1, Phi2) into (ra,dec)
            stream_coord = self.phi_to_radec(
                data["phi1"],
                data["phi2"],
                seed=seed,
                rng=rng,
                **kwargs,
            )
            data.loc[:, "ra"] = stream_coord.icrs.ra.deg
            data.loc[:, "dec"] = stream_coord.icrs.dec.deg

        # Sample missing magnitudes if needed
        mag_bands_missing = [
            f"mag_{band}" for band in bands if f"mag_{band}" not in data.columns
        ]
        if mag_bands_missing:
            stream_config = kwargs.get("stream_config", None)
            if stream_config is None:
                raise ValueError(
                    "`stream_config` must be provided to sample missing magnitudes."
                )
            stream_model = StreamModel(stream_config)
            data = stream_model.complete_catalog(
                data,
                inplace=True,
                columns_to_add=mag_bands_missing,
            )

        for col in required_columns:
            if col not in data.columns:
                raise ValueError(f"Input data must contain '{col}' column.")

        return data

    def phi_to_radec(
        self,
        phi1,
        phi2,
        gc_frame=None,
        seed=None,
        rng=None,
        mask_type=["footprint"],
        **kwargs,
    ):
        """
        Transform stream coordinates (phi1, phi2) to sky coordinates (RA, Dec).

        This method converts stream coordinates to celestial coordinates using a great circle
        frame. If no frame is provided, it automatically finds one randomly chosen such that a given percentile
        of the points lie within the mask defined with mask_type.

        The frame used (whether provided or generated) is stored in ``self._last_gc_frame``
        for potential reuse via ``gc_frame='last'`` in subsequent calls.

        Parameters
        ----------
        phi1, phi2 : array-like
            Stream coordinates in degrees.
        gc_frame : gala.coordinates.GreatCircleICRSFrame or 'last', optional
            Great circle coordinate frame. If None, will be automatically determined.
            If 'last', uses the frame from the previous call (stored in self._last_gc_frame).
        seed : int, optional
            Random seed for reproducible frame selection.
        rng : numpy.random.Generator, optional
            Random number generator instance.
        mask_type : list of str, optional
            Types of masks to use for footprint validation.
            Options: ["footprint", "maglim_g", "maglim_r", "ebv"].
            Default is ["footprint"].
        **kwargs
            Additional keyword arguments passed to _find_gc_frame():

            percentile_threshold : float, optional
                Minimum fraction of points that must be in mask. Default is 0.99.
            max_iter : int, optional
                Maximum number of random trials. Default is 1000.

        Returns
        -------
        astropy.coordinates.SkyCoord
            Sky coordinates in ICRS frame.

        Raises
        ------
        ValueError
            If phi1 and phi2 have different lengths or contain invalid values.
        RuntimeError
            If no suitable great circle frame could be found.

        Examples
        --------
        Convert stream coordinates to sky coordinates:

        >>> phi1 = np.linspace(-10, 10, 1000)
        >>> phi2 = np.zeros_like(phi1)
        >>> coords = injector.phi_to_radec(phi1, phi2, seed=42)

        Reuse the frame from a previous call:

        >>> coords2 = injector.phi_to_radec(phi1_2, phi2_2, gc_frame='last')
        """
        # Input validation
        phi1_arr = np.asarray(phi1, dtype=float)
        phi2_arr = np.asarray(phi2, dtype=float)

        if phi1_arr.size == 0 or phi2_arr.size == 0:
            raise ValueError("phi1 and phi2 cannot be empty arrays")

        # Handle 'last' keyword for gc_frame
        if gc_frame == "last":
            gc_frame = self._last_gc_frame

        # Find or use provided great circle frame
        if gc_frame is None:
            gc_frame = self._find_gc_frame(
                rng=rng,
                seed=seed,
                mask_type=mask_type,
                phi1=phi1_arr,
                phi2=phi2_arr,
                **kwargs,
            )

        # Store the frame for potential reuse
        self._last_gc_frame = gc_frame

        # Transform to sky coordinates
        phi1_deg = phi1_arr * u.deg
        phi2_deg = phi2_arr * u.deg
        stream_coord = coord.SkyCoord(phi1=phi1_deg, phi2=phi2_deg, frame=gc_frame)

        return stream_coord

    def _find_gc_frame(
        self,
        phi1=None,
        phi2=None,
        mask=None,
        mask_type=["footprint"],
        percentile_threshold=0.99,
        max_iter=1000,
        rng=None,
        seed=None,
        verbose=True,
        **kwargs,
    ):
        """
        Find a great circle frame such that a chosen fraction of points lie within the chosen mask.

        This method iteratively tries random great circle orientations until it finds
        one where at least `percentile_threshold` of the stream points fall within
        the survey mask.

        Parameters
        ----------
        phi1, phi2 : array-like, optional
            Stream coordinates to validate against the mask.
        mask : np.ndarray, optional
            Pre-computed HEALPix mask. If None, will be created from mask_type.
        mask_type : list of str, optional
            Types of masks to combine for footprint validation.
            Default is ["footprint"].
        percentile_threshold : float, optional
            Minimum fraction of points that must be within the mask.
            Default is 0.99.
        max_iter : int, optional
            Maximum number of random trials. Default is 1000.
        rng : numpy.random.Generator, optional
            Random number generator instance.
        seed : int, optional
            Random seed if rng is not provided.
        verbose : bool, optional
            Whether to print progress information. Default is True.
        **kwargs
            Additional keyword arguments (currently unused).

        Returns
        -------
        gala.coordinates.GreatCircleICRSFrame or None
            Great circle frame, or None if no suitable frame found after max_iter attempts.
        """
        if rng is None:
            rng = np.random.default_rng(seed)

        # Create the mask if not provided
        if mask is None:
            healpix_mask = self._create_mask(mask_type, verbose=verbose)
        else:
            healpix_mask = mask
            if verbose:
                print("Using provided HEALPix mask for footprint checking.")

        # Do NOT store mask as instance attribute to avoid conflicts between runs
        # (each run may need a different mask)

        # If no mask is available, return a random great circle frame
        if healpix_mask is None:
            if verbose:
                print("No mask available, returning a random great circle frame.")
            end1 = self._random_uniform_skycoord(rng)
            end2 = self._random_uniform_skycoord(rng)
            gc_frame = gc.GreatCircleICRSFrame.from_endpoints(end1, end2)
            return gc_frame

        if phi1 is None or phi2 is None:
            raise ValueError("phi1 and phi2 must be provided if no mask is given.")

        # Iteratively try random great circle frames
        trials = 0
        while trials < max_iter:
            trials += 1

            # Generate random endpoints for the great circle
            end1 = self._random_uniform_skycoord(rng)
            end2 = self._random_uniform_skycoord(rng)
            gc_frame = gc.GreatCircleICRSFrame.from_endpoints(end1, end2)

            # Transform stream points to ICRS and check mask coverage
            pts_gc = coord.SkyCoord(
                phi1=phi1 * u.deg, phi2=phi2 * u.deg, frame=gc_frame
            )
            pts_icrs = pts_gc.transform_to("icrs")
            ra_all = pts_icrs.ra.deg
            dec_all = pts_icrs.dec.deg
            frac = self._fraction_inside_mask(ra_all, dec_all, healpix_mask)
            if frac >= percentile_threshold:
                if verbose:
                    print(
                        f"Found suitable great circle frame after {trials} trials with {frac*100:.2f}% points inside the mask."
                    )
                return gc_frame
        if verbose:
            print(
                f"Could not find a suitable great circle frame after {max_iter} trials."
            )
        return None

    def _random_uniform_skycoord(self, rng):
        """
        Generate a random point uniformly distributed on the sky.

        Parameters
        ----------
        rng : numpy.random.Generator
            Random number generator instance.

        Returns
        -------
        astropy.coordinates.SkyCoord
            Random sky coordinate in ICRS frame.
        """
        z = rng.uniform(-1.0, 1.0)
        dec = np.degrees(np.arcsin(z))
        ra = rng.uniform(0.0, 360.0)
        return coord.SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")

    def _fraction_inside_mask(self, ra_deg, dec_deg, healpix_mask):
        """
        Calculate the fraction of points that fall within valid mask regions.

        Parameters
        ----------
        ra_deg, dec_deg : array-like
            Coordinates in degrees.
        healpix_mask : np.ndarray
            Boolean HEALPix mask array (1=valid, 0=invalid).

        Returns
        -------
        float
            Fraction of points inside mask (range 0.0 to 1.0).
        """
        nside = hp.get_nside(healpix_mask)
        pix_indices = hp.ang2pix(nside, ra_deg, dec_deg, lonlat=True)
        return np.count_nonzero(healpix_mask[pix_indices]) / len(pix_indices)

    def _create_mask(self, mask_type, verbose=True, ebv_threshold=0.2):
        """
        Create a combined boolean mask from specified mask types.

        This method uses a class-level cache to avoid recomputing masks. The cache key
        includes the survey name, mask types, and ebv_threshold to ensure correct cache hits.

        Parameters
        ----------
        mask_type : str, list of str, or None
            Type(s) of masks to combine. Options: ["footprint", "coverage",
            "maglim_<band>", "ebv"]. If None, returns None.
        verbose : bool, optional
            Whether to print status messages. Default is True.
        ebv_threshold : float, optional
            E(B-V) threshold for extinction mask (only used if 'ebv' in mask_type).
            Pixels with E(B-V) > ebv_threshold are masked out. Default is 0.2.

        Returns
        -------
        np.ndarray or None
            Combined boolean mask array (1=valid, 0=invalid), or None if mask_type is None.

        Raises
        ------
        ValueError
            If mask_type is invalid or required maps are missing.
        """
        # Normalize mask_type to list
        if mask_type is None:
            if verbose:
                print("⚠ No mask_type provided to build mask.")
            return None

        if isinstance(mask_type, str):
            mask_type = [mask_type]
        elif not isinstance(mask_type, list):
            raise ValueError("mask_type must be a string, list of strings, or None.")

        # Sort mask_type for consistent cache keys
        mask_type = sorted(mask_type)

        # Create cache key that includes survey name and ebv_threshold if relevant
        survey_name = getattr(self.survey, "name", "unknown")
        uses_ebv = "ebv" in mask_type
        cache_key = (survey_name, tuple(mask_type), ebv_threshold if uses_ebv else None)

        # Check cache first
        if cache_key in self.mask_cache:
            if verbose:
                print(f"✓ Using cached mask for {mask_type}")
            return self.mask_cache[cache_key]

        if verbose:
            print(f"Building new mask for {mask_type}...")

        # Find the minimum nside among the needed maps and collect maps
        nside_target = []
        maps = {}

        for m in mask_type:
            if "maglim" in m:
                band = m.split("_")[-1]
                if band not in self.survey.maglim_maps:
                    raise ValueError(
                        f"Band '{band}' not found in survey magnitude limit maps. Available: {list(self.survey.maglim_maps.keys())}"
                    )
                nside = hp.get_nside(self.survey.maglim_maps[band])
                maps[m] = self.survey.maglim_maps[band]
                nside_target.append(nside)

            elif m in ["coverage", "footprint"]:
                if self.survey.coverage is None:
                    raise ValueError("Survey coverage map is not available.")
                nside = hp.get_nside(self.survey.coverage)
                maps[m] = self.survey.coverage
                nside_target.append(nside)

            elif m == "ebv":
                if self.survey.ebv_map is None:
                    raise ValueError("Survey E(B-V) extinction map is not available.")
                nside = hp.get_nside(self.survey.ebv_map)
                maps[m] = self.survey.ebv_map
                nside_target.append(nside)
            else:
                raise ValueError(
                    f"Unknown mask type: '{m}'. Valid options: 'footprint', 'coverage', 'maglim_<band>', 'ebv'"
                )

        if not nside_target:
            raise ValueError(f"No valid maps found for mask_type: {mask_type}")

        nside_min = min(nside_target)

        # Upgrade/downgrade all maps to the same nside
        for m in maps:
            nside = hp.get_nside(maps[m])
            if nside != nside_min:
                if verbose:
                    print(f"  Resampling {m} from nside={nside} to nside={nside_min}")
                maps[m] = hp.ud_grade(maps[m], nside_min)

        # Initialize combined mask (start with all True)
        npix = hp.nside2npix(nside_min)
        mask_map = np.ones(npix, dtype=bool)

        # Combine the masks with appropriate thresholds
        for m in mask_type:
            if "maglim" in m:
                band = m.split("_")[-1]
                # Valid regions are where magnitude limit is above saturation
                if band in self.survey.saturation:
                    mask_map &= maps[m] > self.survey.saturation[band]
                    if verbose:
                        print(
                            f"  Applied saturation cut for {band} band (> {self.survey.saturation[band]} mag)"
                        )
                else:
                    # If no saturation defined, just check for positive values
                    mask_map &= maps[m] > 0
                    if verbose:
                        print(
                            f"  Applied positivity cut for {band} band (no saturation defined)"
                        )
            elif m in ["coverage", "footprint"]:
                # Valid regions have coverage > 0.5
                mask_map &= maps[m] > 0.5
            elif m == "ebv":
                # Valid regions have low extinction
                mask_map &= maps[m] < ebv_threshold

        # Store in cache
        self.mask_cache[cache_key] = mask_map

        if verbose:
            total_pixels = len(mask_map)
            valid_pixels = np.sum(mask_map)
            coverage_fraction = valid_pixels / total_pixels
            print(f"✓ Mask created: valid pixels fraction = {coverage_fraction:.1f}")
            print(f"  Cached with key: {cache_key}")

        return mask_map

    def sample_measured_magnitudes(self, mag_true, mag_err, **kwargs):
        """
        Sample measured magnitudes from true apparent magnitudes and errors.

        This method adds photometric noise to true magnitudes by sampling
        from a Gaussian distribution in flux space.

        Parameters
        ----------
        mag_true : float or np.ndarray
            True apparent magnitude(s).
        mag_err : float or np.ndarray
            Magnitude error(s).
        **kwargs
            Additional keyword arguments:

            rng : numpy.random.Generator, optional
                Random number generator instance.
            seed : int, optional
                Random seed if rng is not provided.

        Returns
        -------
        np.ndarray or str
            Measured magnitude(s). Returns "BAD_MAG" for objects with negative flux.
        """
        rng = kwargs.pop("rng", None)
        if rng is None:
            seed = kwargs.pop("seed", None)
            rng = np.random.default_rng(seed)

        # Sample the fluxes their errors
        flux_obs = StreamInjector.magToFlux(mag_true) + rng.normal(
            scale=self.getFluxError(mag_true, mag_err)
        )

        # If the flux is negative, set the magnitude to "BAD_MAG" (not detected). Otherwise, convert the flux back to magnitude
        mag_obs = np.where(
            flux_obs > 0.0, StreamInjector.fluxToMag(flux_obs), "BAD_MAG"
        )

        return mag_obs

    def detect_flag(self, pix, mag=None, band="r", **kwargs):
        """
        Apply the survey selection to determine detection flags for stars.

        This method uses the survey completeness function and random sampling
        to determine which stars would be detected by the survey.

        Parameters
        ----------
        pix : int or np.ndarray
            HEALPix pixel index/indices.
        mag : float or np.ndarray, optional
            Magnitude(s). Default is None.
        band : str, optional
            Band to consider for detection. Default is 'r'.
        **kwargs
            Additional keyword arguments:

            rng : numpy.random.Generator, optional
                Random number generator instance.
            seed : int, optional
                Random seed if rng is not provided.
            perfect_galstarsep : bool, optional
                If True, assumes perfect star/galaxy separation. Default is False.

        Returns
        -------
        np.ndarray
            Boolean array: True for detected stars, False otherwise.

        Raises
        ------
        ValueError
            If magnitude values are not provided.
        """

        rng = kwargs.pop("rng", None)
        if rng is None:
            seed = kwargs.pop("seed", None)
            rng = np.random.default_rng(seed)

        # Select the appropriate magnitude and map depending on the band
        maglim = self.survey.get_maglim(band, pixel=pix)

        perfect_galstarsep = kwargs.get("perfect_galstarsep", False)
        if perfect_galstarsep:
            compl = self.survey.get_detection_efficiency(band, mag, maglim)
        else:
            compl = self.survey.get_completeness(band, mag, maglim)

        # Set the threshold using completeness
        threshold = rng.uniform(size=len(mag)) <= compl

        return threshold

    def _save_injected_data(self, data, folder):
        """
        Save the injected data to a CSV file.

        Parameters
        ----------
        data : pd.DataFrame
            Data to save.
        folder : str or Path, optional
            Path to the folder where the file will be saved. If None, uses default
            location in package's data/outputs directory.
        """

        if folder is None:
            current_dir = os.path.dirname(os.path.abspath(__file__))
            folder = os.path.join(current_dir, "..", "data/outputs/")
        if not os.path.exists(folder):
            os.makedirs(folder)

        file_name = folder + "data_injected.csv"
        print(f"Saving injected data to {file_name}")
        data.to_csv(file_name, index=False)

    @staticmethod
    def magToFlux(mag):
        """
        Convert from AB magnitude to flux.

        Parameters
        ----------
        mag : float or np.ndarray
            AB magnitude(s).

        Returns
        -------
        float or np.ndarray
            Flux in Janskys (Jy).
        """
        return 3631.0 * 10 ** (-0.4 * mag)

    @staticmethod
    def fluxToMag(flux):
        """
        Convert from flux to AB magnitude.

        Parameters
        ----------
        flux : float or np.ndarray
            Flux in Janskys (Jy).

        Returns
        -------
        float or np.ndarray
            AB magnitude(s).
        """
        return -2.5 * np.log10(flux / 3631.0)

    @staticmethod
    def getFluxError(mag, mag_error):
        """
        Convert magnitude error to flux error.

        Parameters
        ----------
        mag : float or np.ndarray
            Magnitude(s).
        mag_error : float or np.ndarray
            Magnitude error(s).

        Returns
        -------
        float or np.ndarray
            Flux error in Janskys (Jy).
        """
        return StreamInjector.magToFlux(mag) * mag_error / 1.0857362

    @classmethod
    def clear_mask_cache(cls):
        """
        Clear the mask cache.

        This can be useful if you want to free memory or force masks to be recomputed.

        Examples
        --------
        >>> StreamInjector.clear_mask_cache()
        """
        cls.mask_cache.clear()
        print("✓ Mask cache cleared")

    @classmethod
    def list_cached_masks(cls):
        """
        List all cached masks.

        Returns
        -------
        list of tuples
            List of cache keys (survey_name, mask_types, ebv_threshold)

        Examples
        --------
        >>> StreamInjector.list_cached_masks()
        [('LSST', ('footprint', 'maglim_r'), None),
         ('LSST', ('ebv', 'footprint'), 0.2)]
        """
        if not cls.mask_cache:
            print("No masks cached")
            return []

        print(f"Cached masks ({len(cls.mask_cache)}):")
        for key in cls.mask_cache.keys():
            survey_name, mask_types, ebv_thresh = key
            ebv_str = f", ebv_threshold={ebv_thresh}" if ebv_thresh is not None else ""
            print(f"  - {survey_name}: {list(mask_types)}{ebv_str}")
        return list(cls.mask_cache.keys())

    def plot_stream_in_mask(self, data, mask_type, ebv_threshold=0.2, **kwargs):
        """
        Plot the stream over the footprint mask.

        Creates a visualization showing the stream's position relative to the
        survey footprint or other masks.

        Parameters
        ----------
        data : pd.DataFrame
            Data containing 'ra' and 'dec' columns.
        mask_type : str or list of str
            Type(s) of masks to plot. Options: ["footprint", "coverage",
            "maglim_<band>", "ebv"].
        ebv_threshold : float, optional
            E(B-V) threshold (only used if 'ebv' in mask_type). Default is 0.2.
        **kwargs
            Additional arguments passed to plotting function:

            output_folder : str, optional
                Path to save the figure.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure object.
        ax : matplotlib.axes.Axes
            The axes object.

        Raises
        ------
        ValueError
            If mask cannot be created from mask_type parameter.

        Examples
        --------
        Plot stream in footprint:

        >>> fig, ax = injector.plot_stream_in_mask(data, ['footprint', 'maglim_r'])

        Plot with custom E(B-V) threshold:

        >>> fig, ax = injector.plot_stream_in_mask(
        ...     data, ['footprint', 'ebv'], ebv_threshold=0.15
        ... )
        """
        # Get or create the mask
        mask = self._create_mask(mask_type, verbose=False, ebv_threshold=ebv_threshold)

        if mask is None:
            raise ValueError("Could not create mask. Check mask_type parameter.")

        # Call the plotting function
        fig, ax = plot_stream_in_mask(
            data["ra"], data["dec"], mask, output_folder=kwargs.get("output_folder")
        )
        return fig, ax
