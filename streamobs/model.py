#!/usr/bin/env python
"""
Models for simulating streams.
"""

import copy
import warnings

import numpy as np
import pandas as pd

from streamobs.columns import true_col
from streamobs.functions import function_factory
from streamobs.samplers import sampler_factory

# Roman ugali isochrones are delivered in Vega magnitudes, but our catalogs are
# AB. These are the per-band AB - Vega offsets (AB = Vega + diff) for the Roman
# WFI filters — the mode of the by-chip zeropoints from the Roman technical
# information (Roman_zeropoints_20240301.ecsv), as used by the
# rubin_roman_object_classification prototype. The conversion is applied
# unconditionally to any Roman band (a no-op for non-Roman bands).
#
# TODO: this Vega->AB conversion really belongs in ugali, so isochrones are
# returned natively in AB. Move it upstream and delete this table once that lands.
ROMAN_VEGA_TO_AB = {
    "F062": 0.153,
    "F087": 0.481,
    "F106": 0.660,
    "F129": 1.051,
    "F146": 1.164,
    "F158": 1.315,
    "F184": 1.556,
    "F213": 1.837,
}


class ConfigurableModel(object):
    """Baseclass for models built from configs."""

    def __init__(self, config, **kwargs):
        """Initialize with configuration.

        Parameters
        ----------
        config : dict or None
            Configuration used by subclasses to build internal components.
        **kwargs
            Optional overrides merged into ``config`` before building.
        """
        self._config = copy.deepcopy(config)
        if self._config is not None:
            self._config.update(**kwargs)
            self._create_model()

    def _create_model(self):
        pass

    def sample(self, size):
        pass


class StreamModel(ConfigurableModel):
    """High-level object for the various components of the stream model."""

    def __init__(self, config, **kwargs):
        """Create the stream from the config object.

        Parameters
        ----------
        config : dict
            Configuration sections: ``density``, ``track``, ``distance_modulus``,
            ``isochrone``, and optionally ``velocity``.
        **kwargs
            Optional overrides merged into ``config`` before building.
        """
        super().__init__(config, **kwargs)

    def _create_model(self):
        """Instantiate sub-models from configuration sections."""
        self.density = self._create_density()
        self.track = self._create_track()
        self.distance_modulus = self._create_distance_modulus()
        self.isochrone = self._create_isochrone()
        self.velocity = self._create_velocity()

    def _create_density(self):
        """Build density sampler from ``config['density']``."""
        config = self._config.get("density")
        return DensityModel(config)

    def _create_linear_density(self):
        """to be used with the cubic spline methods"""
        config = self._config.get("linear_density")
        return DensityModel(config)

    def _create_track(self):
        """Build track model (center and spread functions) from ``config['track']``."""
        config = self._config.get("track")
        return TrackModel(config)

    def _create_distance_modulus(self):
        """Build distance-modulus track from ``config['distance_modulus']`` if present."""
        config = self._config.get("distance_modulus")
        if config:
            return TrackModel(config)
        else:
            return None

    def _create_isochrone(self):
        """Build isochrone model from ``config['isochrone']`` if present."""
        config = self._config.get("isochrone")
        if config:
            iso = IsochroneModel(config)
            iso.create_isochrone(config)
            return iso
        else:
            return None

    def _create_velocity(self):
        """Build velocity model from ``config['velocity']`` if present."""
        config = self._config.get("velocity")
        if config:
            return VelocityModel(config)
        else:
            return None

    def sample(self, size):
        """Sample stream stars and derived quantities.

        Parameters
        ----------
        size : int
            Number of stars to generate.

        Returns
        -------
        pandas.DataFrame
            Columns include: ``phi1``, ``phi2``, ``dist``, ``mu1``, ``mu2``,
            ``rv``, and the isochrone magnitude columns ``<survey>_<band>_true``
            (per survey/band). Some may be None if the sub-model is absent.
        """

        # Sample phi1 and phi2
        phi1 = self.density.sample(size)
        phi2 = self.track.sample(phi1)

        # Sample distances
        if self.distance_modulus:
            dist = self.distance_modulus.sample(phi1)
        else:
            dist = None

        # Sample kinematics
        if self.velocity:
            mu1, mu2, rv = self.velocity.sample(phi1)
        else:
            mu1, mu2, rv = None, None, None

        # Create the DataFrame of stream stars
        df = pd.DataFrame(
            {
                "phi1": phi1,
                "phi2": phi2,
                "dist": dist,
                "mu1": mu1,
                "mu2": mu2,
                "rv": rv,
            }
        )

        # Sample magnitudes from isochrone (band-/survey-general columns)
        if self.isochrone:
            mag_data = self._sample_iso_mags(size, dist)
        else:
            mag_data = {}  # no isochrone -> no magnitude columns
        for col, vals in mag_data.items():
            df[col] = vals

        return df

    def _iso_mag_columns(self):
        """Names of the magnitude columns produced by the isochrone model.

        Always ``[<survey>_<band>_true, ...]`` for every survey/band the
        isochrone carries (a single-survey isochrone simply has one survey); no
        isochrone ⇒ ``[]``. ``IsochroneModel`` tracks ``surveys`` /
        ``survey_bands`` in both configuration forms, so the naming is uniform.
        """
        iso = self.isochrone
        if iso is None:
            return []
        cols = []
        for name in iso.surveys:
            band_1, band_2 = iso.survey_bands[name]
            cols += [true_col(band_1, name), true_col(band_2, name)]
        return cols

    def _sample_iso_mags(self, n, dist, masses=None):
        """Sample isochrone magnitudes as a ``{column: values}`` dict.

        Returns each survey's ``<namespace>_<band>_true`` columns plus the shared
        ``mass`` column (the initial masses used for every band). When ``masses``
        is given it is used directly instead of an IMF draw, so the sampled
        magnitudes reproduce those exact stars.
        """
        mags, masses = self.isochrone.sample_multisurvey(n, dist, masses=masses)
        cols = {true_col(band, name): vals for (name, band), vals in mags.items()}
        cols["mass"] = masses
        return cols

    def complete_catalog(
        self,
        catalog,
        columns_to_add=None,
        size=None,
        inplace=False,
        save_path=None,
        verbose=True,
        dist=None,
    ):
        """Complete only the requested columns in a catalog.

        This method takes an input catalog (or a desired size when no catalog
        is provided) and fills in only the requested stream-model columns
        while preserving pre-existing non-null values. Columns are generated
        using the configured sub-models (density, track, distance modulus,
        isochrone, velocity) and only if those capabilities are available.

        Pre-existing values are never overwritten: for every column only the
        missing (absent or NaN) rows are filled. In particular, supplying some
        of an isochrone's bands and requesting the others fills only the missing
        bands and leaves the provided ones untouched.

        Parameters
        ----------
        catalog : pandas.DataFrame or str or dict or None
            Input catalog. If a string, it is interpreted as a CSV filepath
            to read. If a dict, it will be converted to a DataFrame.
            If None, ``size`` must be provided to create an empty frame of
            that length.
        columns_to_add : sequence of str or None, optional
            The columns to ensure in the output. Valid entries are
            {'phi1','phi2','dist','mu1','mu2','rv'} plus the isochrone magnitude
            columns (``<survey>_<band>_true``). If None, all valid columns
            supported by the configured model are considered.
        size : int or None, optional
            Required when ``catalog`` is None or an empty table; ignored
            otherwise.
        inplace : bool, default False
            If True and a DataFrame or CSV path is provided, modify that
            object in place (for CSV, overwrite the input file).
        save_path : str or None, optional
            If provided, write the completed catalog to this CSV path.
        verbose : bool, default True
            If True, print progress/status messages.
        dist : float or array-like or None, optional
            Distance modulus to use directly instead of sampling one from the
            ``distance_modulus`` sub-model. A scalar is broadcast to every row
            that needs a ``dist`` value; an array must have one entry per row.
            When given, ``phi1`` and a ``distance_modulus`` model are **not**
            required to fill magnitudes. Only missing ``dist`` rows are set.

        Returns
        -------
        pandas.DataFrame
            The completed catalog. If ``inplace`` is True and a DataFrame was
            provided, the same object is returned after modification.

        Raises
        ------
        ValueError
            If ``size`` is required but not provided, or when dependencies are
            missing (e.g., requesting 'phi2' without available 'phi1').

        Notes
        -----
        - Dependencies: 'phi2' and 'dist' require 'phi1'. Magnitudes require
          'dist' and an isochrone model. Velocities require 'phi1' and a
          velocity model.
        - Existing non-null values are preserved: only the missing rows are
          filled for ``phi1``/``phi2``/``dist``, the magnitude columns, and the
          shared ``mass`` column (supplying some bands and requesting others
          fills only the missing ones, colour-consistently). Velocities are the
          exception — ``mu1``/``mu2``/``rv`` are recomputed for the whole columns
          to keep kinematic coherence across rows.
        - When ``catalog`` is a CSV path and ``inplace`` is True, the original
          file is overwritten.
        """
        # Supported outputs and capability filtering
        # Columns this method can fill using the configured model
        # Magnitude columns are survey-namespaced (<survey>_<band>_true).
        mag_cols = self._iso_mag_columns()
        # The isochrone also produces the shared initial-mass column.
        mass_cols = ("mass",) if self.isochrone is not None else ()
        all_cols = (
            ("phi1", "phi2", "dist")
            + tuple(mag_cols)
            + mass_cols
            + ("mu1", "mu2", "rv")
        )
        target_cols = (
            list(all_cols)
            if columns_to_add is None
            else [c for c in columns_to_add if c in all_cols]
        )
        unknown = (
            []
            if columns_to_add is None
            else sorted(set(columns_to_add) - set(all_cols))
        )
        if unknown:
            warnings.warn(f"Ignoring unknown columns: {unknown}")

        if self.velocity is None:
            removed = [c for c in target_cols if c in ("mu1", "mu2", "rv")]
            target_cols = [c for c in target_cols if c not in ("mu1", "mu2", "rv")]
            if removed:
                self._info(verbose, "Velocity model not defined; skipping velocities.")
        if self.distance_modulus is None and dist is None:
            removed = [c for c in target_cols if c == "dist"]
            target_cols = [c for c in target_cols if c != "dist"]
            if removed:
                self._info(
                    verbose, "Distance modulus model not defined; skipping distances."
                )

        # Load/normalize input catalog
        df, src_path = self._open_catalog(catalog, size=size, inplace=inplace)
        N = len(df)

        # phi1
        if "phi1" in target_cols:
            idx = self._missing_idx(df, "phi1")
            if len(idx) > 0:
                if self.density is None:
                    raise ValueError("Density model is required to sample phi1")
                df.loc[idx, "phi1"] = self.density.sample(len(idx))
                self._info(verbose, f"Filled {len(idx)} phi1 values.")

        # phi2 (needs phi1)
        if "phi2" in target_cols:
            if "phi1" not in df.columns or df["phi1"].isna().any():
                raise ValueError(
                    "phi1 required to sample phi2; include 'phi1' in columns_to_add or provide it in catalog"
                )
            idx = self._missing_idx(df, "phi2")
            if len(idx) > 0:
                if self.track is None:
                    raise ValueError("Track model is required to sample phi2")
                df.loc[idx, "phi2"] = self.track.sample(df.loc[idx, "phi1"].to_numpy())
                self._info(verbose, f"Filled {len(idx)} phi2 values.")

        # dist (needs phi1)
        if "dist" in target_cols or (
            any(c in target_cols for c in mag_cols)
            and any(c not in df.columns for c in mag_cols)
        ):
            idx = self._missing_idx(df, "dist")
            if len(idx) > 0:
                if dist is not None:
                    # Use the distance supplied directly (scalar broadcast or
                    # per-row vector); no phi1 / distance_modulus model needed.
                    dist_arr = np.asarray(dist, dtype=float)
                    if dist_arr.ndim == 0:
                        df.loc[idx, "dist"] = float(dist_arr)
                    else:
                        if dist_arr.shape[0] != N:
                            raise ValueError(
                                f"dist vector has length {dist_arr.shape[0]} but "
                                f"the catalog has {N} rows."
                            )
                        pos = df.index.get_indexer(idx)
                        df.loc[idx, "dist"] = dist_arr[pos]
                    self._info(verbose, f"Set {len(idx)} dist values from `dist`.")
                else:
                    if self.distance_modulus is None:
                        raise ValueError(
                            "No distance_modulus model is configured; pass `dist` "
                            "(a float or per-row vector) to set distances."
                        )
                    if "phi1" not in df.columns or df["phi1"].isna().any():
                        raise ValueError(
                            "phi1 required to sample dist; include 'phi1' in columns_to_add or provide it in catalog"
                        )

                    df.loc[idx, "dist"] = self.distance_modulus.sample(
                        df.loc[idx, "phi1"].to_numpy()
                    )
                    self._info(verbose, f"Filled {len(idx)} dist values.")

        # magnitudes + shared initial mass (need dist and isochrone)
        requested_mags = [c for c in mag_cols if c in target_cols]
        want_mass = "mass" in target_cols and self.isochrone is not None
        fill_targets = requested_mags + (["mass"] if want_mass else [])
        if fill_targets:
            # Only touch rows that are missing a requested column; existing values
            # are preserved (never overwritten).
            missing = {c: self._missing_idx(df, c) for c in fill_targets}
            to_fill = {c: idx for c, idx in missing.items() if len(idx) > 0}
            if not to_fill:
                self._info(
                    verbose,
                    f"{fill_targets} already present; no sampling performed.",
                )
            else:
                # Verify distance availability
                if "dist" not in df.columns:
                    raise ValueError(
                        "dist is required to sample apparent magnitudes; include 'dist' in `columns_to_add`, provide it in the catalog, or pass `dist=`."
                    )
                dist_vals = df["dist"].to_numpy()
                # Reuse a fully-present input `mass` column as the initial masses
                # so the sampled magnitudes reproduce the user's simulation stars;
                # otherwise draw fresh masses from the IMF.
                masses_in = None
                if "mass" in df.columns and df["mass"].notna().all():
                    masses_in = df["mass"].to_numpy()
                # One shared mass draw -> all bands; the newly filled cells are
                # mutually colour-consistent. Assign only the missing rows so any
                # bands/values already present are left untouched.
                mags = self._sample_iso_mags(N, dist_vals, masses=masses_in)
                for col, idx in to_fill.items():
                    pos = df.index.get_indexer(idx)
                    if col not in df.columns:
                        df[col] = np.nan
                    df.loc[idx, col] = np.asarray(mags[col])[pos]
                self._info(
                    verbose,
                    f"Filled {sorted(to_fill)} (missing rows only).",
                )

        # velocities (need phi1 and velocity model)
        if any(c in target_cols for c in ("mu1", "mu2", "rv")):
            if any(c in df.columns for c in ("mu1", "mu2", "rv")):
                self._info(
                    verbose, "Velocity components already exist; no sampling performed."
                )
            else:
                if "phi1" not in df.columns or df["phi1"].isna().any():
                    raise ValueError("phi1 required to sample velocities")
                mu1, mu2, rv = self.velocity.sample(df["phi1"].to_numpy())
                if "rv" in df.columns or "mu1" in df.columns or "mu2" in df.columns:
                    self._info(
                        verbose,
                        "Overwriting existing velocity components to keep consistency.",
                    )
                df["mu1"] = mu1
                df["mu2"] = mu2
                df["rv"] = rv
                self._info(verbose, f"Filled velocities for {N} rows.")

        if save_path is not None:
            df.to_csv(save_path, index=False)
            self._info(verbose, f"Saved completed catalog to {save_path}.")
        elif isinstance(catalog, str) and inplace:
            df.to_csv(src_path, index=False)
            self._info(verbose, f"Overwrote original catalog at {src_path}.")

        return df

    def _missing_idx(self, df: pd.DataFrame, col: str):
        """Return indices of rows needing a given column.

        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame to inspect.
        col : str
            Column name to check.

        Returns
        -------
        pandas.Index
            Index of rows where ``col`` is missing or NaN. If ``col`` is not
            present in ``df``, returns all row indices.
        """
        if col not in df.columns:
            return df.index
        return df.index[df[col].isna()]

    def _info(self, verbose: bool, msg: str):
        """Conditional verbose print.

        Parameters
        ----------
        verbose : bool
            When True, print ``msg``; otherwise, do nothing.
        msg : str
            Message to print.
        """
        if verbose:
            print(msg)

    def _open_catalog(self, catalog, size=None, inplace=False):
        """Load or normalize the input catalog to a DataFrame.

        Parameters
        ----------
        catalog : pandas.DataFrame or str or dict or None
            Input catalog. If str, treated as a CSV file path. If dict,
            converted to a DataFrame. If None, ``size`` must be provided.
        size : int or None, optional
            Required when ``catalog`` is None or an empty table; ignored
            otherwise.
        inplace : bool, default False
            If True and ``catalog`` is a DataFrame, return it as-is; otherwise
            return a copy to avoid side effects.

        Returns
        -------
        df : pandas.DataFrame
            The loaded or constructed DataFrame.
        src_path : str or None
            The source CSV path if ``catalog`` was a string; otherwise None.

        Raises
        ------
        ValueError
            If ``size`` is required but not provided.
        TypeError
            If ``catalog`` is not one of the supported types.
        """
        src_path = None
        if catalog is None:
            if size is None:
                raise ValueError("size must be provided when catalog is None")
            df = pd.DataFrame(index=np.arange(int(size)))
        elif isinstance(catalog, str):
            src_path = catalog
            df = pd.read_csv(catalog)
        elif isinstance(catalog, pd.DataFrame):
            df = catalog if inplace else catalog.copy()
        elif isinstance(catalog, dict):
            df = pd.DataFrame(catalog)
        else:
            raise TypeError("catalog must be None, path, DataFrame, or dict")

        if len(df) == 0:
            if size is None:
                raise ValueError("Empty catalog; provide size")
            df = pd.DataFrame(index=np.arange(int(size)))

        df = self._standardize_columns_name(df)

        return df, src_path

    def _standardize_columns_name(self, catalog):
        """Standardize column names in the catalog DataFrame.

        Parameters
        ----------
        catalog : pandas.DataFrame
            Input catalog.

        Returns
        -------
        pandas.DataFrame
            Catalog with standardized column names.
        """
        # Mapping of possible column name variants to standard names
        col_mapping = {
            "dist": ["dist", "distance", "distance_modulus"],
            "g_true": ["g_true", "mag_g", "g_mag", "g", "gmag", "magnitude_g"],
            "r_true": ["r_true", "mag_r", "r_mag", "r", "rmag", "magnitude_r"],
            "phi1": ["phi1", "phi_1", "Phi1", "Phi_1"],
            "phi2": ["phi2", "phi_2", "Phi2", "Phi_2"],
            "mu1": ["mu1", "mu_1"],
            "mu2": ["mu2", "mu_2"],
            "rv": ["rv", "radial_velocity", "v_radial"],
        }

        # Create reverse mapping for renaming
        reverse_mapping = {}
        for standard_name, variants in col_mapping.items():
            for var in variants:
                reverse_mapping[var.lower()] = standard_name

        catalog = catalog.rename(columns=reverse_mapping)

        return catalog


class DensityModel(ConfigurableModel):
    """Density along the stream; samples ``phi1`` positions."""

    def _create_model(self):
        """Instantiate the density sampler from configuration."""
        kwargs = copy.deepcopy(self._config)
        type_ = kwargs.pop("type").lower()
        self.density = sampler_factory(type_, **kwargs)

    def sample(self, size):
        """Draw ``phi1`` samples.

        Parameters
        ----------
        size : int
            Number of samples.

        Returns
        -------
        numpy.ndarray
            Sampled ``phi1`` values.
        """
        return self.density.sample(size)


class TrackModel(ConfigurableModel):
    """Transverse track model; samples ``phi2`` given ``phi1``."""

    def _create_model(self):
        """Build center/spread functions from configuration."""
        kwargs = copy.deepcopy(self._config["center"])
        type_ = kwargs.pop("type").lower()
        self.center = function_factory(type_, **kwargs)

        kwargs = copy.deepcopy(self._config["spread"])
        type_ = kwargs.pop("type").lower()
        self.spread = function_factory(type_, **kwargs)

    def _create_sampler(self, x):
        """Create the sampler (Gaussian or Uniform) at positions ``x``."""
        type_ = self._config.get("sampler", "Gaussian").lower()
        if type_ == "gaussian":
            mu = self.center(x)
            sigma = self.spread(x)
            kwargs = dict(mu=mu, sigma=sigma)
        elif type_ == "uniform":
            xmin = self.center(x) - self.spread(x)
            xmax = xmin + 2 * self.spread(x)
            kwargs = dict(xmin=xmin, xmax=xmax)
        else:
            raise Exception(f"Unrecognized sampler: {type_}")

        self._sampler = sampler_factory(type_, **kwargs)

    def sample(self, x):
        """Sample ``phi2`` at given ``phi1`` positions ``x``.

        Parameters
        ----------
        x : array-like
            ``phi1`` positions where to sample ``phi2``.

        Returns
        -------
        numpy.ndarray
            Sampled ``phi2`` values.
        """
        size = len(x)
        self._create_sampler(x)
        return self._sampler.sample(size)


class DistanceModel(ConfigurableModel):
    pass


class IsochroneModel(ConfigurableModel):
    """Isochrone wrapper using ``ugali`` for CMD sampling.

    Two configuration forms are supported:

    - **Single-survey** (legacy): the isochrone section carries the ``ugali``
      factory keys directly (``name``, ``survey``, ``age``, ``z``, ``band_1``,
      ``band_2``, ...). :meth:`sample` returns ``(mag_band_1, mag_band_2)`` and
      reproduces the previous behaviour exactly.
    - **Multi-survey**: a ``surveys`` mapping
      ``{survey_name: {survey, band_1, band_2}}`` plus shared keys
      (``name``, ``age``, ``z``, ...) at the top level. One ``ugali`` isochrone
      is built per survey from the *same* stellar population, so a single shared
      draw of initial masses (:meth:`sample_masses`) is interpolated into every
      survey's bands — giving the same physical star consistent magnitudes
      across surveys. :meth:`sample_multisurvey` returns
      ``{(survey, band): apparent_mag}``.

    Roman bands are always converted from Vega to AB (see
    :data:`ROMAN_VEGA_TO_AB`); other bands pass through unchanged.
    """

    # Defaults for the shared isochrone mass grid (see ugali Isochrone.sample).
    # Masses are drawn from this discretized grid, so ``_MASS_STEPS`` bounds the
    # number of *distinct* masses a stream can contain. A convergence check
    # (g/r magnitude percentiles + distinct-mass count for a 5000-star stream)
    # showed 1000 steps yields only ~220 distinct masses (granular CMD) and
    # ~0.03 mag median scatter, while 4000 steps roughly triples the distinct
    # masses (~600) and tightens convergence to <0.015 mag — for a negligible
    # one-time grid-sampling cost. Override per call via ``mass_steps=``.
    _MASS_MIN = 0.1
    _MASS_STEPS = 4000

    def __init__(self, config, **kwargs):
        super().__init__(config, **kwargs)

    def create_isochrone(self, config):
        """Construct the underlying ``ugali`` isochrone(s) from configuration.

        Both configuration forms are normalized to a single ``{namespace:
        factory_cfg}`` mapping and built through the same loop — a legacy flat
        config simply becomes a one-entry mapping — so there is no separate
        single- vs. multi-survey code path.

        Parameters
        ----------
        config : dict
            Isochrone factory configuration. A ``surveys`` key selects the
            multi-survey form; otherwise the single-survey (legacy) form is used.
        """
        survey_configs, shared = self._normalize_iso_config(config)
        self._build_isochrones(survey_configs, shared)

    def _normalize_iso_config(self, config):
        """Coerce either config form into ``({namespace: factory_cfg}, shared)``.

        Multi-survey: the ``surveys`` mapping is returned verbatim (its keys are
        the column namespaces) with the top-level keys as shared stellar params.
        Single-survey (legacy flat): a one-entry mapping keyed by
        ``{survey}_{release}`` (or just ``{survey}``), matching
        :attr:`streamobs.surveys.Survey.namespace`, with no shared params.
        """
        if "surveys" in config:
            self.multi_survey = True
            shared = {k: v for k, v in config.items() if k != "surveys"}
            return dict(config["surveys"]), shared

        self.multi_survey = False
        if "distance_modulus" in config:
            warnings.warn(
                'Please use the "distance_modulus" section of the configuration '
                "file, instead of the isochrone section, to define a distance modulus."
            )
        survey = config.get("survey")
        release = config.get("release")
        namespace = f"{survey}_{release}" if release else survey
        return {namespace: dict(config)}, {}

    def _build_iso(self, factory_config):
        """Build one ``ugali`` isochrone with its distance modulus reset to 0.

        ``release`` is a column-namespacing concept (it distinguishes survey
        versions in the output column names), not a ``ugali`` factory argument,
        so it is stripped before the isochrone is constructed.
        """
        import ugali.isochrone

        factory_config = {k: v for k, v in factory_config.items() if k != "release"}
        iso = ugali.isochrone.factory(**factory_config)
        iso.params["distance_modulus"].set_bounds([0, 50])
        iso.distance_modulus = 0
        return iso

    def _build_isochrones(self, survey_configs, shared):
        """Build one isochrone per namespace, sharing the top-level stellar params.

        Drives both configuration forms (a legacy flat config is just a
        one-entry ``survey_configs``). The first entry is the primary isochrone
        that drives the shared mass draw and the legacy :meth:`sample`.
        """
        self.isos = {}
        self.survey_bands = {}
        self.surveys = []
        for name, scfg in survey_configs.items():
            factory_config = {**shared, **scfg}
            self.isos[name] = self._build_iso(factory_config)
            self.survey_bands[name] = (scfg.get("band_1"), scfg.get("band_2"))
            self.surveys.append(name)
        # Primary isochrone drives the shared mass draw and the legacy sample().
        self.survey_name = self.surveys[0]
        self.iso = self.isos[self.survey_name]
        self.band_1, self.band_2 = self.survey_bands[self.survey_name]

    def sample_masses(self, nstars, rng=None, mass_min=None, mass_steps=None):
        """Draw ``nstars`` initial stellar masses from the shared isochrone IMF.

        The masses are drawn once from the primary isochrone's mass PDF and are
        meant to be interpolated into each survey's bands, so the same physical
        star gets consistent magnitudes across surveys.

        Parameters
        ----------
        nstars : int
            Number of stars to draw (returns exactly this many).
        rng : numpy.random.Generator, optional
            Random number generator (a default one is created if omitted).
        mass_min, mass_steps : float, int, optional
            Passed to the ``ugali`` isochrone sampler.

        Returns
        -------
        numpy.ndarray
            Initial masses, shape ``(nstars,)``.
        """
        rng = np.random.default_rng() if rng is None else rng
        mass_min = self._MASS_MIN if mass_min is None else mass_min
        mass_steps = self._MASS_STEPS if mass_steps is None else mass_steps
        grid = self.iso.sample(mass_min=mass_min, mass_steps=mass_steps)
        mass_init, mass_pdf = grid[0], grid[1]
        pdf = mass_pdf / mass_pdf.sum()
        return rng.choice(mass_init, size=int(nstars), p=pdf)

    def _absolute_mags(self, iso, masses, mass_min=None, mass_steps=None):
        """Interpolate a survey isochrone's absolute mags at given init masses."""
        mass_min = self._MASS_MIN if mass_min is None else mass_min
        mass_steps = self._MASS_STEPS if mass_steps is None else mass_steps
        grid = iso.sample(mass_min=mass_min, mass_steps=mass_steps)
        mass_init, mag_1, mag_2 = grid[0], grid[3], grid[4]
        order = np.argsort(mass_init)
        mass_init = mass_init[order]
        return (
            np.interp(masses, mass_init, mag_1[order]),
            np.interp(masses, mass_init, mag_2[order]),
        )

    def _to_ab(self, band, mag):
        """Convert Roman bands from Vega to AB (``AB = Vega + offset``).

        Applied unconditionally: Roman bands use the offset in
        :data:`ROMAN_VEGA_TO_AB`, every other band passes through unchanged.
        ``ugali`` delivers Roman isochrones in Vega while our catalogs are AB.

        TODO: this really belongs in ugali (return AB natively); remove once it
        does.
        """
        return mag + ROMAN_VEGA_TO_AB.get(band, 0.0)

    @staticmethod
    def _add_distance_modulus(abs_mag, distance_modulus):
        """Add a scalar or per-star distance modulus to absolute magnitudes."""
        if distance_modulus is None:
            return abs_mag
        return abs_mag + np.asarray(distance_modulus, dtype=float)

    def sample_multisurvey(
        self, nstars, distance_modulus, rng=None, masses=None, **kwargs
    ):
        """Sample apparent magnitudes for every ``(survey, band)``.

        A single shared set of initial masses is interpolated into each survey's
        bands, so the same physical star is consistent across surveys. The masses
        are drawn from the shared IMF (:meth:`sample_masses`) unless supplied via
        ``masses``.

        Parameters
        ----------
        nstars : int
            Number of stars (and the required length of ``masses`` if given).
        distance_modulus : float or array-like
            Distance modulus per star (broadcast if scalar).
        rng : numpy.random.Generator, optional
            Used only when ``masses`` is None.
        masses : array-like, optional
            Initial stellar masses to use directly — e.g. an external
            simulation's masses — instead of drawing from the IMF. Must have
            length ``nstars``.

        Returns
        -------
        dict, numpy.ndarray
            ``{(survey_name, band): apparent_magnitude_array}`` and the initial
            masses used (shape ``(nstars,)``).
        """
        if masses is None:
            masses = self.sample_masses(
                nstars,
                rng=rng,
                mass_min=kwargs.get("mass_min"),
                mass_steps=kwargs.get("mass_steps"),
            )
        else:
            masses = np.asarray(masses, dtype=float)
            if len(masses) != int(nstars):
                raise ValueError(
                    f"masses has length {len(masses)} but nstars={int(nstars)}."
                )
        out = {}
        for name in self.surveys:
            band_1, band_2 = self.survey_bands[name]
            abs_1, abs_2 = self._absolute_mags(self.isos[name], masses)
            abs_1 = self._to_ab(band_1, abs_1)
            abs_2 = self._to_ab(band_2, abs_2)
            out[(name, band_1)] = self._add_distance_modulus(abs_1, distance_modulus)
            out[(name, band_2)] = self._add_distance_modulus(abs_2, distance_modulus)
        return out, masses

    def sample(self, nstars, distance_modulus, rng=None, masses=None, **kwargs):
        """Simulate magnitudes in the two bands of the (primary) isochrone.

        Draws *exactly* ``nstars`` stars: a fixed set of initial masses is drawn
        once from the isochrone IMF (:meth:`sample_masses`) and interpolated into
        the two bands. This differs from the historical behaviour, where
        ``nstars`` was converted to a total stellar mass and ``ugali``'s
        ``simulate`` returned a random-length IMF realization (count generally
        ``!= nstars``). The fixed-count semantics are required so the *same
        physical star* can be shared across surveys (see
        :meth:`sample_multisurvey`).

        Parameters
        ----------
        nstars : int
            Number of stars to simulate (returns exactly this many).
        distance_modulus : float or array-like
            Distance modulus per star (broadcast if scalar).

        Returns
        -------
        tuple of numpy.ndarray
            ``(mag_band_1, mag_band_2)`` arrays. For a multi-survey isochrone
            this returns the primary survey's two bands; use
            :meth:`sample_multisurvey` to get every survey's bands.
        """
        mags, _ = self.sample_multisurvey(
            nstars, distance_modulus, rng=rng, masses=masses, **kwargs
        )
        return (
            mags[(self.survey_name, self.band_1)],
            mags[(self.survey_name, self.band_2)],
        )

    def _dist_to_modulus(self, distance):
        """
        Convert physical distances in pc into distance modulus
        """
        if distance is None:
            return 0
        elif np.all(distance == 0):
            warnings.warn(
                "Distances are equal to 0, distance modulus has been set to 0."
            )
            return 0
        else:
            return 5 * np.log10(distance) - 5


class VelocityModel(ConfigurableModel):
    """Placeholder for velocity model."""

    def sample(self, phi1):
        """Placeholder"""
        warnings.warn("VelocityModel not implemented!")

        if np.isscalar(phi1):
            mu1, mu2, rv = np.nan, np.nan, np.nan
        else:
            mu1, mu2, rv = np.nan * np.ones_like([phi1, phi1, phi1])

        return mu1, mu2, rv


class BackgroundModel(StreamModel):
    """Background model."""

    pass


class SplineStreamModel(StreamModel):
    """Spline-based stream model with linear-density component."""

    def __init__(self, config, **kwargs):
        """Create spline stream from configuration.

        Parameters
        ----------
        config : dict
            Must include ``linear_density`` and ``track`` sections. Optional
            ``stream_name`` is propagated to sub-sections.
        **kwargs
            Optional overrides merged into ``config`` before building.
        """

        stream_name = None
        if config["stream_name"]:
            stream_name = config["stream_name"]

        config["linear_density"]["stream_name"] = stream_name
        config["track"]["center"]["stream_name"] = stream_name
        config["track"]["spread"]["stream_name"] = stream_name
        super().__init__(config, **kwargs)

    def _create_model(self):
        """Instantiate spline-specific components and common sub-models."""
        self.density = self._create_linear_density()
        self.track = self._create_track()
        self.distance_modulus = self._create_distance_modulus()
        self.isochrone = self._create_isochrone()
        self.velocity = self._create_velocity()
