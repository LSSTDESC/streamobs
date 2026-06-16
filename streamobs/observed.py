#!/usr/bin/env python

import copy
import os

import astropy.coordinates as coord
import astropy.units as u
import gala.coordinates as gc
import healpy as hp
import numpy as np
import pandas as pd

from .columns import err_col, flag_col, obs_col, perfect_flag_col, true_col
from .model import StreamModel
from .plotting import plot_stream_in_mask
from .surveys import Survey


class StreamInjector:
    """
    Inject observational effects into stream data for one or more surveys.

    A single injector handles both the single- and multi-survey cases: pass one
    survey or several. The same shared sky placement and a single shared draw of
    true magnitudes (masses sampled once via the isochrone and interpolated into
    every survey's bands) guarantee the **same physical star** gets consistent
    magnitudes across surveys. Each survey contributes its own
    ``<survey>_<band>_obs`` / ``<survey>_<band>_err`` / ``<survey>_flag_observed``
    columns, computed with that survey's own maglim maps and completeness
    functions. Output columns are **always** survey-namespaced, even for a single
    survey.

    All survey data is loaded once and cached, making multiple injections
    efficient.

    Attributes
    ----------
    surveys : dict
        ``{namespace: Survey}`` for every survey this injector serves. The
        namespace is the column prefix (``lsst_r_obs``, ``roman_F158_obs``, ...).
    primary : str
        Namespace of the survey whose footprint drives the shared sky placement
        and whose ``_save_injected_data`` is used.
    survey : Survey
        The primary :class:`~streamobs.surveys.Survey` (convenience accessor used
        by the mask / coordinate helpers; equal to ``surveys[primary]``).
    mask_cache : dict (class attribute)
        Cache of previously created HEALPix masks to avoid recomputation.
    _last_gc_frame : GreatCircleICRSFrame or None
        The most recently used great circle frame. This allows reusing the same
        sky location across multiple inject() calls when gc_frame='last'.

    Examples
    --------
    Single survey (columns namespaced by the survey's name):

    >>> injector = StreamInjector('lsst', release='dc2')
    >>> out = injector.inject(df, bands=['r', 'g'])  # -> lsst_r_obs, ...

    Several surveys at once:

    >>> injector = StreamInjector({'lsst': 'lsst', 'roman': 'roman'})
    >>> out = injector.inject(
    ...     df, survey_bands={'lsst': ['r', 'g'], 'roman': ['F106', 'F158']},
    ...     stream_config=scene['stream'], seed=42,
    ... )
    """

    mask_cache = {}

    def __init__(self, survey, primary=None, **kwargs):
        """
        Initialize with one or more survey configurations.

        Parameters
        ----------
        survey : str, Survey, dict, or list
            One survey or several. Accepted forms:

            - a survey-name string (e.g. ``'lsst'``) or a pre-loaded
              :class:`~streamobs.surveys.Survey` — a single survey, namespaced
              by its own name;
            - a ``{namespace: spec}`` dict, where ``spec`` is a name string or a
              ``Survey`` and the key is the column namespace;
            - a list/tuple of specs, each namespaced by its survey's name.
        primary : str, optional
            Namespace of the survey that drives the shared sky placement.
            Defaults to the first survey.
        **kwargs
            Forwarded to :meth:`Survey.load` for any ``spec`` given as a name
            string (e.g. ``release``).

        Raises
        ------
        ValueError
            If ``surveys`` is empty or of an unsupported type, or if ``primary``
            is not one of the provided surveys.
        """
        self.surveys = self._normalize_surveys(survey, **kwargs)
        self.survey_names = list(self.surveys)
        if not self.survey_names:
            raise ValueError("At least one survey is required.")
        self.primary = primary if primary is not None else self.survey_names[0]
        if self.primary not in self.surveys:
            raise ValueError(
                f"primary='{self.primary}' is not one of {self.survey_names}."
            )

        # Instance attribute to store the last used gc_frame
        self._last_gc_frame = None

    @property
    def survey(self):
        """The primary :class:`~streamobs.surveys.Survey`.

        Mask, coordinate and footprint helpers operate on this survey.
        """
        return self.surveys[self.primary]

    @classmethod
    def _normalize_surveys(cls, surveys, **kwargs):
        """Coerce the ``surveys`` argument into a ``{namespace: Survey}`` dict."""
        if isinstance(surveys, (str, Survey)):
            survey = cls._load_survey(surveys, **kwargs)
            return {survey.name: survey}
        if isinstance(surveys, (list, tuple)):
            out = {}
            for spec in surveys:
                survey = cls._load_survey(spec, **kwargs)
                out[survey.name] = survey
            return out
        if isinstance(surveys, dict):
            return {
                name: cls._load_survey(spec, **kwargs) for name, spec in surveys.items()
            }
        raise ValueError(
            "surveys must be a survey name, a Survey, a {name: spec} dict, "
            "or a list of specs."
        )

    @staticmethod
    def _load_survey(spec, **kwargs):
        """Resolve a single survey spec (name string or Survey) to a Survey."""
        if isinstance(spec, Survey):
            return spec
        if isinstance(spec, str):
            return Survey.load(survey=spec, **kwargs)
        raise ValueError("Each survey spec must be a string or Survey instance.")

    def inject(self, data, survey_bands=None, bands=None, stream_config=None, **kwargs):
        """
        Add observed quantities from every survey into a single catalog.

        Applies observational effects (photometric errors, measured magnitudes,
        detection flags) for each survey this injector serves. A single shared
        sky placement and a single shared true-magnitude fill (masses sampled
        once and interpolated into every survey's bands) ensure the same physical
        star is consistent across surveys. Output columns are always
        survey-namespaced (``<survey>_<band>_obs`` etc.).

        Parameters
        ----------
        data : str or pd.DataFrame
            Input data as DataFrame or path to the file (CSV or Excel). May
            contain only stream coordinates (``phi1``/``phi2`` or ``ra``/``dec``);
            anything missing is sampled from ``stream_config``. An all-empty frame
            of length ``N`` is accepted (geometry and magnitudes are then sampled
            for ``N`` rows).
        survey_bands : dict, optional
            ``{survey_name: [bands]}`` — bands to inject for each survey. Keys
            must match the surveys this injector was built with. For a single
            survey you may instead pass ``bands`` (below).
        bands : list of str, optional
            Convenience shorthand for the single-survey case: the bands to inject
            for the (only) survey. Ignored when ``survey_bands`` is given. If
            neither is provided and there is exactly one survey, defaults to
            ``['r', 'g']``.
        stream_config : dict, optional
            The ``stream`` section consumed by
            :class:`~streamobs.model.StreamModel`. Required when any coordinate
            or true-magnitude column is missing. Its isochrone produces the
            shared ``<survey>_<band>_true`` columns.
        **kwargs
            Additional keyword arguments:

            seed : int, optional
                Random seed for reproducibility.
            dist : float or array-like, optional
                Distance modulus used directly (scalar broadcast or per-row
                vector) instead of sampling from the config's ``distance_modulus``
                model — lets magnitudes be filled without ``phi1``.
            nside : int, optional
                HEALPix nside parameter. Default is 4096.
            detection_mag_cut : list of str, optional
                Non-reference bands to apply the explicit SNR>=5 cut to. The
                reference band (``survey.completeness_band``) is always cut via
                the selection functions, so the default here is every injected
                band *except* the reference band. Net effect: every injected band
                must have SNR >= 5, with the reference band counted once.
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
            DataFrame with shared ``ra``/``dec`` and, per survey:

            - <survey>_<band>_true : True (noiseless) apparent magnitudes
            - <survey>_<band>_obs : Observed (noisy) magnitudes
            - <survey>_<band>_err : Reported photometric errors
            - <survey>_flag_observed : Boolean detection flag (detection and
              classification efficiencies)
            - <survey>_flag_perfect_galstarsep : Boolean flag assuming perfect
              star/galaxy separation (only if perfect_galstarsep=True)

        Raises
        ------
        ValueError
            If required columns are missing, or if ``survey_bands`` references an
            unknown survey, or if neither ``survey_bands`` nor ``bands`` can be
            resolved.
        """
        survey_bands = self._resolve_survey_bands(survey_bands, bands)

        # Load data
        data = self._load_data(data)

        # Set the seed for reproducibility
        seed = kwargs.pop("seed", None)
        rng = np.random.default_rng(seed)

        # Shared sky placement + shared true-magnitude fill (masses sampled once
        # across all surveys).
        data = self._complete_shared(
            data,
            survey_bands,
            stream_config=stream_config,
            rng=rng,
            seed=seed,
            **kwargs,
        )

        # Per-survey observational injection. Independent child RNGs make the
        # result independent of survey ordering and reproducible from `seed`.
        children = rng.spawn(len(self.survey_names))
        for child_rng, name in zip(children, self.survey_names):
            if name not in survey_bands:
                continue
            data = self._inject_one_survey(
                data,
                list(survey_bands[name]),
                survey=self.surveys[name],
                survey_namespace=name,
                rng=child_rng,
                **kwargs,
            )

        # Save if requested
        if kwargs.get("save"):
            self._save_injected_data(data, kwargs.get("folder", None))

        # Return data (do NOT store as instance attribute to avoid conflicts between runs)
        return data

    def _resolve_survey_bands(self, survey_bands, bands):
        """Normalize the ``survey_bands`` / ``bands`` arguments to a dict.

        ``survey_bands`` wins if given (validated against the known surveys).
        Otherwise ``bands`` is accepted only when there is a single survey; if
        neither is given and there is one survey, defaults to ``['r', 'g']``.
        """
        if survey_bands is not None:
            unknown = set(survey_bands) - set(self.surveys)
            if unknown:
                raise ValueError(
                    f"survey_bands references unknown surveys {sorted(unknown)}; "
                    f"available: {self.survey_names}."
                )
            return {name: list(b) for name, b in survey_bands.items()}

        if len(self.survey_names) != 1:
            raise ValueError(
                "Pass survey_bands={name: [bands]} when the injector serves "
                f"multiple surveys ({self.survey_names})."
            )
        if bands is None:
            bands = ["r", "g"]
        return {self.primary: list(bands)}

    def _inject_one_survey(
        self, data, bands, survey, survey_namespace, rng=None, seed=None, **kwargs
    ):
        """Add one survey's observed magnitudes, errors and detection flags.

        This holds the per-band observational logic. It assumes ``data`` already
        carries ``ra``/``dec`` and the true-magnitude columns
        (``true_col(band, survey_namespace)``) for the requested bands; it does
        **not** sample positions or true magnitudes.

        Parameters
        ----------
        data : pandas.DataFrame
            Catalog with ``ra``/``dec`` and ``true_col(band, survey_namespace)``
            for every requested band.
        bands : list of str
            Bands to process for this survey.
        survey : Survey
            The survey supplying maglim maps, completeness and error curves.
        survey_namespace : str
            Column-naming namespace ⇒ ``<survey>_<band>_obs`` /
            ``<survey>_<band>_err`` / ``<survey>_flag_observed``.
        rng : numpy.random.Generator, optional
            Random generator for the noise draw and detection sampling.
        seed : int, optional
            Seed used to build an RNG when ``rng`` is None.
        **kwargs
            ``nside``, ``detection_mag_cut``, ``dust_correction``,
            ``perfect_galstarsep``, ``verbose`` (see :meth:`inject`).

        Returns
        -------
        pandas.DataFrame
            ``data`` with this survey's observed columns and detection flag(s).
        """
        if rng is None:
            rng = np.random.default_rng(seed)

        verbose = kwargs.get("verbose", True)
        perfect_galstarsep = kwargs.pop("perfect_galstarsep", False)

        # Get HEALPix pixel indices
        nside = kwargs.pop("nside", 4096)
        pix = hp.ang2pix(nside, data["ra"], data["dec"], lonlat=True)

        # Initialize detection flags (will be updated per band)
        flag_completeness_band = None
        flag_detection_only_band = None

        # Process each band
        for band in bands:
            # Get extinction for this band
            nside_ebv = hp.get_nside(survey.ebv_map)
            if nside_ebv != nside:
                pix_ebv = hp.ang2pix(nside_ebv, data["ra"], data["dec"], lonlat=True)
            else:
                pix_ebv = pix
            extinction_band = survey.get_extinction(band, pixel=pix_ebv)

            # Calculate true apparent magnitudes (including extinction)
            apparent_mag_true = data[true_col(band, survey_namespace)] + extinction_band

            # Get magnitude limits
            nside_maglim = hp.get_nside(survey.maglim_maps[band])
            if nside_maglim != nside:
                pix_maglim = hp.ang2pix(
                    nside_maglim, data["ra"], data["dec"], lonlat=True
                )
            else:
                pix_maglim = pix

            # Calculate photometric errors. The *sample* error (true scatter)
            # drives the noise draw; the *catalog* error (reported) is written
            # as magerr and used for the S/N cut. When no sample curve is loaded,
            # the sample error falls back to the catalog error, so the two are
            # identical (outputs unchanged from the single-curve behaviour).
            maglim_vals = survey.get_maglim(band, pixel=pix_maglim)
            mag_err_sample = survey.get_photo_error(
                band, apparent_mag_true, maglim_vals, kind="sample"
            )
            mag_err = survey.get_photo_error(
                band, apparent_mag_true, maglim_vals, kind="catalog"
            )

            # Sample measured magnitudes using the sample (true-scatter) error
            mag_obs = self.sample_measured_magnitudes(
                apparent_mag_true,
                mag_err_sample,
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
                    obs_col(band, survey_namespace): mag_obs,
                    err_col(band, survey_namespace): mag_err,
                }
            )

            # Reset indices and concatenate
            data = data.reset_index(drop=True)
            new_columns = new_columns.reset_index(drop=True)
            data = pd.concat([data, new_columns], axis=1)

            # Compute detection flag for completeness-band (reference band)
            if band == survey.completeness_band:
                flag_completeness_band = self.detect_flag(
                    pix_maglim,
                    survey=survey,
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
                        survey=survey,
                        mag=apparent_mag_true,
                        band=band,
                        rng=rng,
                        seed=seed,
                        perfect_galstarsep=True,
                        **kwargs,
                    )

        # Apply detection threshold
        if flag_completeness_band is None:
            if survey.completeness_band in bands:
                raise ValueError(
                    f"flag_completeness_{survey.completeness_band} must be computed for detection in {survey.completeness_band} band."
                )
            else:
                raise ValueError(
                    f"Detection flag requires '{survey.completeness_band}' band to be in bands."
                )

        # Build combined detection flags
        # Start with flux validity check (not BAD_MAG) across all injected bands
        flag_valid_flux = None
        for band in bands:
            band_valid = data[obs_col(band, survey_namespace)] != "BAD_MAG"
            flag_valid_flux = (
                band_valid if flag_valid_flux is None else flag_valid_flux & band_valid
            )

        # Combine with completeness
        flag_observed = (
            flag_valid_flux & flag_completeness_band
            if flag_completeness_band is not None
            else flag_valid_flux
        )
        if perfect_galstarsep:
            flag_perfect = (
                flag_valid_flux & flag_detection_only_band
                if flag_detection_only_band is not None
                else flag_valid_flux
            )

        # Apply SNR cuts.
        #
        # The SNR>=SNR_min cut on the *reference* band (``survey.completeness_band``)
        # is owned by the survey's selection functions — they are estimated with
        # that cut already applied. We therefore apply it to the reference band
        # exactly **once** here, to both flags. (``get_completeness`` bakes it in,
        # so for ``flag_observed`` this is idempotent; ``get_detection_efficiency``
        # does not, so for ``flag_perfect`` this line is what supplies it.) See
        # the "S/N cut ownership" note in docs/source/roman_multisurvey_plan.md
        # for the path to folding it into the efficiency curve itself (option a).
        #
        # Every *other* injected band gets the cut applied explicitly below.
        SNR_min = 5.0
        ref_band = survey.completeness_band

        if ref_band in bands:
            SNR_ref = 1.0 / data[err_col(ref_band, survey_namespace)]
            flag_observed &= SNR_ref >= SNR_min
            if perfect_galstarsep:
                flag_perfect &= SNR_ref >= SNR_min

        # Non-reference bands. Default: every injected band except the reference
        # band (whose cut is handled above via the selection functions).
        detection_mag_cut = kwargs.get(
            "detection_mag_cut", [b for b in bands if b != ref_band]
        )
        for band in detection_mag_cut:
            if band == ref_band:
                # The reference band's SNR cut is already applied above.
                continue
            if band not in bands:
                if verbose:
                    print(
                        f"Warning: SNR cut requested for {band}-band but it's not in bands list. Skipping."
                    )
                continue
            if verbose:
                print(f"Applying detection cut on {band}-band with SNR >= {SNR_min}")
            SNR = 1.0 / data[err_col(band, survey_namespace)]
            flag_observed &= SNR >= SNR_min
            if perfect_galstarsep:
                flag_perfect &= SNR >= SNR_min

        # Store flags in DataFrame
        data[flag_col(survey_namespace)] = flag_observed
        if perfect_galstarsep:
            data[perfect_flag_col(survey_namespace)] = flag_perfect

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

    def complete_data(
        self,
        data,
        survey_bands=None,
        bands=None,
        stream_config=None,
        dist=None,
        **kwargs,
    ):
        """Complete the columns the injector needs, filling the rest from the config.

        Public helper: give it a (possibly partial) catalog and it returns one
        with everything the injector requires present — sky coordinates
        (``ra``/``dec``, converting from ``phi1``/``phi2`` if needed) and the
        per-survey true-magnitude columns ``<survey>_<band>_true``. Anything
        already present is preserved; only missing columns are sampled, using
        ``stream_config`` (a :class:`~streamobs.model.StreamModel` config). The
        stellar masses are drawn **once** and interpolated into every survey's
        bands, so the same physical star is consistent across surveys.

        This is the same completion :meth:`inject` runs internally, exposed so
        you can build/inspect a completed catalog without injecting noise.

        Parameters
        ----------
        data : str or pandas.DataFrame
            Input catalog (or path). May contain only stream coordinates
            (``phi1``/``phi2`` or ``ra``/``dec``), an all-empty frame of length
            ``N``, or any subset of the target columns.
        survey_bands : dict, optional
            ``{survey_name: [bands]}`` — the true-magnitude columns to ensure per
            survey. For a single survey you may instead pass ``bands``; if neither
            is given and there is exactly one survey, defaults to ``['r', 'g']``.
        bands : list of str, optional
            Single-survey shorthand for ``survey_bands={primary: bands}``.
        stream_config : dict, optional
            :class:`~streamobs.model.StreamModel` config used to sample any
            missing geometry / true magnitudes. Required only when something is
            missing.
        dist : float or array-like or None, optional
            Distance modulus to use directly (scalar broadcast or per-row
            vector) instead of sampling from the config's ``distance_modulus``
            model. Lets magnitudes be filled without ``phi1`` / a distance model.
        **kwargs
            ``rng`` / ``seed`` for reproducibility, plus ``gc_frame``,
            ``mask_type``, ``percentile_threshold``, ``max_iter`` forwarded to
            :meth:`phi_to_radec` when converting ``phi1``/``phi2`` → ``ra``/``dec``.

        Returns
        -------
        pandas.DataFrame
            A copy of the input with ``ra``/``dec`` and the requested
            ``<survey>_<band>_true`` columns present.

        Raises
        ------
        ValueError
            If neither (ra, dec) nor (phi1, phi2) are present, or if columns are
            missing and ``stream_config`` is not provided.

        Examples
        --------
        >>> df = pd.DataFrame({'phi1': [-5, 0, 5], 'phi2': [0, 0, 0]})
        >>> out = injector.complete_data(df, bands=['r', 'g'], stream_config=cfg)
        >>> out = injector.complete_data(
        ...     df, survey_bands={'lsst': ['r', 'g'], 'roman': ['F106', 'F158']},
        ...     stream_config=cfg, seed=42,
        ... )
        """
        survey_bands = self._resolve_survey_bands(survey_bands, bands)
        data = self._load_data(data).copy()

        rng = kwargs.pop("rng", None)
        seed = kwargs.pop("seed", None)
        if rng is None:
            rng = np.random.default_rng(seed)

        return self._complete_shared(
            data,
            survey_bands,
            stream_config=stream_config,
            rng=rng,
            seed=seed,
            dist=dist,
            **kwargs,
        )

    def _ensure_radec(self, data, rng=None, seed=None, **kwargs):
        """Ensure ``ra``/``dec`` are present, converting from (phi1, phi2) if needed.

        Existing ``ra``/``dec`` are left untouched. When converting, the great
        circle frame is found (or reused via ``gc_frame='last'``) through
        :meth:`phi_to_radec`, using the primary survey's footprint.
        """
        if "ra" in data.columns and "dec" in data.columns:
            return data

        if "phi1" not in data.columns or "phi2" not in data.columns:
            raise ValueError(
                "Input data must contain either (ra, dec) or (phi1, phi2) columns."
            )

        # Handle 'last' keyword for gc_frame
        if kwargs.get("gc_frame", None) == "last":
            kwargs["gc_frame"] = self._last_gc_frame

        # Convert coordinates (phi1, phi2) into (ra, dec)
        stream_coord = self.phi_to_radec(
            data["phi1"], data["phi2"], seed=seed, rng=rng, **kwargs
        )
        data.loc[:, "ra"] = stream_coord.icrs.ra.deg
        data.loc[:, "dec"] = stream_coord.icrs.dec.deg
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

    def detect_flag(self, pix, mag=None, band="r", survey=None, **kwargs):
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
        survey : Survey, optional
            Survey whose completeness/detection-efficiency curves to use.
            Defaults to the primary survey.
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

        if survey is None:
            survey = self.survey

        # Select the appropriate magnitude and map depending on the band
        maglim = survey.get_maglim(band, pixel=pix)

        perfect_galstarsep = kwargs.get("perfect_galstarsep", False)
        if perfect_galstarsep:
            compl = survey.get_detection_efficiency(band, mag, maglim)
        else:
            compl = survey.get_completeness(band, mag, maglim)

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

    def _complete_shared(
        self,
        data,
        survey_bands,
        stream_config=None,
        rng=None,
        seed=None,
        dist=None,
        **kwargs,
    ):
        """Fill shared geometry, ``ra``/``dec`` and per-survey true magnitudes.

        Positions and masses are drawn *once* so all surveys describe the same
        physical stars (the isochrone produces every survey's
        ``<survey>_<band>_true`` column from one shared mass draw). Existing
        columns are preserved (only missing values are filled). ``ra``/``dec``
        are placed using the primary survey's footprint. ``dist`` (a float or
        per-row vector) overrides the model's distance sampling when given.
        """
        verbose = kwargs.get("verbose", True)

        true_cols = []
        for name, bands in survey_bands.items():
            true_cols += [true_col(b, name) for b in bands]

        have_radec = "ra" in data.columns and "dec" in data.columns
        have_phi = "phi1" in data.columns and "phi2" in data.columns
        missing_true = [c for c in true_cols if c not in data.columns]

        # Sample stream geometry and/or the shared true magnitudes from the model.
        need_phi = not have_radec and not have_phi
        if need_phi or missing_true:
            if stream_config is None:
                raise ValueError(
                    "stream_config is required to sample stream geometry/magnitudes."
                )
            stream_model = StreamModel(stream_config)
            cols_to_add = []
            if need_phi:
                cols_to_add += ["phi1", "phi2"]
            # `dist` is needed before magnitudes; the model fills it (from the
            # distance_modulus model or the supplied `dist`) if absent.
            cols_to_add += ["dist"] + missing_true
            data = stream_model.complete_catalog(
                data,
                columns_to_add=cols_to_add,
                inplace=True,
                verbose=verbose,
                dist=dist,
            )

        # Convert (phi1, phi2) -> (ra, dec) using the primary survey footprint.
        if not have_radec:
            data = self._ensure_radec(data, rng=rng, seed=seed, **kwargs)

        return data
