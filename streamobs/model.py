#!/usr/bin/env python
"""
Models for simulating streams.
"""

import copy
import warnings

import numpy as np
import pandas as pd

from streamobs.functions import function_factory
from streamobs.samplers import sampler_factory


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
            ``rv``, ``mag_g``, ``mag_r`` (some may be None if the sub-model is absent).
        """

        # Sample phi1 and phi2
        phi1 = self.density.sample(size)
        phi2 = self.track.sample(phi1)

        # Sample distances
        if self.distance_modulus:
            dist = self.distance_modulus.sample(phi1)
        else:
            dist = None

        # Sample magnitudes from isochrone
        if self.isochrone:

            mag_g, mag_r = self.isochrone.sample(size, dist)
        else:
            mag_g, mag_r = None, None

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
                "mag_g": mag_g,
                "mag_r": mag_r,
            }
        )
        return df

    def complete_catalog(
        self,
        catalog,
        columns_to_add=None,
        size=None,
        inplace=False,
        save_path=None,
        verbose=True,
    ):
        """Complete only the requested columns in a catalog.

        This method takes an input catalog (or a desired size when no catalog
        is provided) and fills in only the requested stream-model columns
        while preserving pre-existing non-null values. Columns are generated
        using the configured sub-models (density, track, distance modulus,
        isochrone, velocity) and only if those capabilities are available.

        Parameters
        ----------
        catalog : pandas.DataFrame or str or dict or None
            Input catalog. If a string, it is interpreted as a CSV filepath
            to read. If a dict, it will be converted to a DataFrame.
            If None, ``size`` must be provided to create an empty frame of
            that length.
        columns_to_add : sequence of str or None, optional
            The columns to ensure in the output. Valid entries are
            {'phi1','phi2','dist','mag_g','mag_r','mu1','mu2','rv'}.
            If None, all valid columns supported by the configured model are
            considered.
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
        - Existing non-null values are preserved; only missing rows are filled,
          except for magnitudes and velocities where the method intentionally
          overwrites the whole columns to keep internal consistency (e.g.,
          colors and kinematic coherence across rows).
        - When ``catalog`` is a CSV path and ``inplace`` is True, the original
          file is overwritten.
        """
        # Supported outputs and capability filtering
        # Columns this method can fill using the configured model
        all_cols = ("phi1", "phi2", "dist", "mag_g", "mag_r", "mu1", "mu2", "rv")
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

        if self.isochrone is None:
            removed = [c for c in target_cols if c in ("mag_g", "mag_r")]
            target_cols = [c for c in target_cols if c not in ("mag_g", "mag_r")]
            if removed:
                self._info(verbose, "Isochrone model not defined; skipping magnitudes.")
        if self.velocity is None:
            removed = [c for c in target_cols if c in ("mu1", "mu2", "rv")]
            target_cols = [c for c in target_cols if c not in ("mu1", "mu2", "rv")]
            if removed:
                self._info(verbose, "Velocity model not defined; skipping velocities.")
        if self.distance_modulus is None:
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
            any(c in target_cols for c in ("mag_g", "mag_r"))
            and any(c not in df.columns for c in ("mag_g", "mag_r"))
        ):
            idx = self._missing_idx(df, "dist")
            if len(idx) > 0:
                if "phi1" not in df.columns or df["phi1"].isna().any():
                    raise ValueError(
                        "phi1 required to sample dist; include 'phi1' in columns_to_add or provide it in catalog"
                    )

                df.loc[idx, "dist"] = self.distance_modulus.sample(
                    df.loc[idx, "phi1"].to_numpy()
                )
                self._info(verbose, f"Filled {len(idx)} dist values.")

        # magnitudes (need dist and isochrone)
        if any(c in target_cols for c in ("mag_g", "mag_r")):
            if "mag_g" in df.columns and "mag_r" in df.columns:
                self._info(
                    verbose, "'mag_g' and 'mag_r' already exist; no sampling performed."
                )
            else:
                # Verify distance modulus availability
                if "dist" not in df.columns:
                    raise ValueError(
                        "dist is required to sample apparent magnitudes; include 'dist' in `columns_to_add` or provide it in catalog"
                    )
                dist_vals = df["dist"].to_numpy()
                mag_g, mag_r = self.isochrone.sample(N, dist_vals)
                # Add both magnitudes to keep colors consistent across rows
                if "mag_g" in df.columns or "mag_r" in df.columns:
                    self._info(
                        verbose,
                        "Overwriting existing mag_g and/or mag_r to keep colors consistent.",
                    )
                df["mag_g"] = mag_g
                df["mag_r"] = mag_r
                self._info(verbose, f"Filled magnitudes for {N} rows.")

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
            "mag_g": ["mag_g", "g_mag", "g", "gmag", "magnitude_g"],
            "mag_r": ["mag_r", "r_mag", "r", "rmag", "magnitude_r"],
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
    """Isochrone wrapper using ``ugali`` for CMD sampling."""

    def __init__(self, config, **kwargs):
        super().__init__(config, **kwargs)

    def create_isochrone(self, config):
        """Construct the underlying ``ugali`` isochrone from configuration.

        Parameters
        ----------
        config : dict
            Isochrone factory configuration.
        """
        import ugali.isochrone

        self.iso = ugali.isochrone.factory(**config)
        self.iso.params["distance_modulus"].set_bounds([0, 50])
        if "distance_modulus" in config:
            warnings.warn(
                'Please use the "distance_modulus" section of the configuration file, instead of the isochrone section, to define a distance module.'
            )
        self.iso.distance_modulus = 0

    def sample(self, nstars, distance_modulus, **kwargs):
        """Simulate magnitudes in g and r bands.

        Parameters
        ----------
        nstars : int
            Number of stars to simulate.
        distance_modulus : float or array-like
            Distance modulus per star (broadcast if scalar).

        Returns
        -------
        tuple of numpy.ndarray
            ``(mag_g, mag_r)`` arrays.
        """
        stellar_mass = nstars * self.iso.stellar_mass()
        if np.isscalar(distance_modulus):
            mag_g, mag_r = self.iso.simulate(
                stellar_mass, distance_modulus=self.iso.distance_modulus
            )
            mag_g, mag_r = [
                mag + np.ones_like(mag) * distance_modulus for mag in (mag_g, mag_r)
            ]
        else:
            mag_g, mag_r = self.iso.simulate(
                stellar_mass, distance_modulus=self.iso.distance_modulus
            )
            mag_g, mag_r = [mag + distance_modulus for mag in (mag_g, mag_r)]

        return mag_g, mag_r

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
        self.distance_modulus = self._create_distance()
        self.isochrone = self._create_isochrone()
        self.velocity = self._create_velocity()
