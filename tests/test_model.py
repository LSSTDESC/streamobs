"""
tests/test_model.py
===================
Tests for ``streamobs.model.StreamModel`` and its sub-models.

Focus: verify that every quantity the model can produce is actually sampled
and has the right shape / dtype, without requiring the optional ugali /
velocity backends.
"""

import numpy as np
import pandas as pd
import pytest

from streamobs.model import (
    DensityModel,
    StreamModel,
    TrackModel,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

N = 200  # number of stars used in every sampling test


# ---------------------------------------------------------------------------
# DensityModel
# ---------------------------------------------------------------------------

class TestDensityModel:
    """Unit tests for DensityModel (phi1 sampler)."""

    def test_uniform_sample_shape(self):
        cfg = {"type": "uniform", "xmin": -10.0, "xmax": 10.0}
        model = DensityModel(cfg)
        samples = model.sample(N)
        assert samples.shape == (N,), "Sample array must have length N"

    def test_uniform_sample_range(self):
        cfg = {"type": "uniform", "xmin": -5.0, "xmax": 5.0}
        model = DensityModel(cfg)
        samples = model.sample(1_000)
        assert np.all(samples >= -5.0), "All samples must be >= xmin"
        assert np.all(samples <= 5.0),  "All samples must be <= xmax"

    def test_gaussian_sample_shape(self):
        cfg = {"type": "gaussian", "mu": 0.0, "sigma": 1.0}
        model = DensityModel(cfg)
        samples = model.sample(N)
        assert samples.shape == (N,)

    def test_gaussian_sample_dtype(self):
        cfg = {"type": "gaussian", "mu": 0.0, "sigma": 1.0}
        model = DensityModel(cfg)
        samples = model.sample(N)
        assert np.issubdtype(samples.dtype, np.floating)


# ---------------------------------------------------------------------------
# TrackModel
# ---------------------------------------------------------------------------

class TestTrackModel:
    """Unit tests for TrackModel (phi2 sampler given phi1)."""

    @pytest.fixture
    def flat_track_config(self):
        return {
            "center": {"type": "polynomial", "coeffs": [0.0]},
            "spread": {"type": "polynomial", "coeffs": [0.1]},
            "sampler": "gaussian",
        }

    def test_sample_shape_matches_input(self, flat_track_config):
        model = TrackModel(flat_track_config)
        phi1 = np.linspace(-10, 10, N)
        phi2 = model.sample(phi1)
        assert phi2.shape == (N,)

    def test_flat_track_mean_near_zero(self, flat_track_config):
        """With a flat center at 0, the sample mean should be close to 0."""
        model = TrackModel(flat_track_config)
        phi1 = np.linspace(-10, 10, 5_000)
        phi2 = model.sample(phi1)
        assert abs(phi2.mean()) < 0.05, "Mean phi2 should be near 0 for a flat track"

    def test_uniform_sampler(self):
        cfg = {
            "center": {"type": "polynomial", "coeffs": [1.0]},
            "spread": {"type": "polynomial", "coeffs": [0.5]},
            "sampler": "uniform",
        }
        model = TrackModel(cfg)
        phi1 = np.linspace(-5, 5, N)
        phi2 = model.sample(phi1)
        # With center=1 and spread=0.5, phi2 must lie in [0.5, 1.5]
        assert np.all(phi2 >= 0.5 - 1e-9)
        assert np.all(phi2 <= 1.5 + 1e-9)

    def test_unknown_sampler_raises(self):
        cfg = {
            "center": {"type": "polynomial", "coeffs": [0.0]},
            "spread": {"type": "polynomial", "coeffs": [0.1]},
            "sampler": "nonexistent_sampler",
        }
        model = TrackModel(cfg)
        with pytest.raises(Exception):
            model.sample(np.array([0.0]))


# ---------------------------------------------------------------------------
# StreamModel — geometry only (no isochrone / velocity)
# ---------------------------------------------------------------------------

class TestStreamModelGeometry:
    """StreamModel sampling tests that do not require ugali or velocity."""

    def test_sample_returns_dataframe(self, minimal_stream_config):
        model = StreamModel(minimal_stream_config)
        df = model.sample(N)
        assert isinstance(df, pd.DataFrame)

    def test_sample_has_required_columns(self, minimal_stream_config):
        model = StreamModel(minimal_stream_config)
        df = model.sample(N)
        for col in ("phi1", "phi2"):
            assert col in df.columns, f"Column '{col}' missing from output"

    def test_sample_size(self, minimal_stream_config):
        model = StreamModel(minimal_stream_config)
        df = model.sample(N)
        assert len(df) == N

    def test_phi1_within_bounds(self, minimal_stream_config):
        """phi1 should respect the density's uniform bounds."""
        model = StreamModel(minimal_stream_config)
        df = model.sample(2_000)
        assert df["phi1"].between(-10.0, 10.0).all()

    def test_phi2_is_numeric(self, minimal_stream_config):
        model = StreamModel(minimal_stream_config)
        df = model.sample(N)
        assert pd.api.types.is_float_dtype(df["phi2"])

    def test_no_distance_column_is_none(self, minimal_stream_config):
        """Without a distance_modulus section, 'dist' column should be all NaN/None."""
        model = StreamModel(minimal_stream_config)
        df = model.sample(N)
        assert df["dist"].isna().all(), "'dist' should be NaN when no distance_modulus configured"

    def test_no_magnitudes_when_no_isochrone(self, minimal_stream_config):
        model = StreamModel(minimal_stream_config)
        df = model.sample(N)
        assert df["mag_g"].isna().all()
        assert df["mag_r"].isna().all()

    def test_no_kinematics_when_no_velocity(self, minimal_stream_config):
        model = StreamModel(minimal_stream_config)
        df = model.sample(N)
        for col in ("mu1", "mu2", "rv"):
            assert df[col].isna().all(), f"'{col}' should be NaN when no velocity model"

    @pytest.mark.parametrize("size", [1, 10, 100, 1_000])
    def test_various_sample_sizes(self, minimal_stream_config, size):
        model = StreamModel(minimal_stream_config)
        df = model.sample(size)
        assert len(df) == size


# ---------------------------------------------------------------------------
# StreamModel — with distance modulus
# ---------------------------------------------------------------------------

class TestStreamModelWithDistance:

    def test_dist_column_populated(self, stream_config_with_distance):
        model = StreamModel(stream_config_with_distance)
        df = model.sample(N)
        assert "dist" in df.columns
        assert not df["dist"].isna().all(), "'dist' should contain values"

    def test_dist_values_positive(self, stream_config_with_distance):
        """Distance modulus values should be positive (in mag)."""
        model = StreamModel(stream_config_with_distance)
        df = model.sample(N)
        assert (df["dist"] > 0).all()


# ---------------------------------------------------------------------------
# StreamModel — complete_catalog
# ---------------------------------------------------------------------------

class TestStreamModelCompleteCatalog:

    def test_complete_catalog_from_size(self, minimal_stream_config):
        """complete_catalog should build a catalog of the requested size."""
        model = StreamModel(minimal_stream_config)
        df = model.complete_catalog(catalog=None, size=N, columns_to_add=["phi1", "phi2"], verbose=False)
        assert len(df) == N

    def test_complete_catalog_preserves_existing(self, minimal_stream_config):
        """Pre-existing non-null phi1 values must be kept intact."""
        model = StreamModel(minimal_stream_config)
        existing_phi1 = np.linspace(-5, 5, N)
        catalog = pd.DataFrame({"phi1": existing_phi1})
        df = model.complete_catalog(
            catalog=catalog,
            columns_to_add=["phi1", "phi2"],
            verbose=False,
        )
        np.testing.assert_array_equal(df["phi1"].values, existing_phi1)

    def test_complete_catalog_warns_unknown_columns(self, minimal_stream_config):
        model = StreamModel(minimal_stream_config)
        with pytest.warns(UserWarning):
            model.complete_catalog(
                catalog=None,
                size=N,
                columns_to_add=["phi1", "totally_fake_column"],
                verbose=False,
            )

    def test_complete_catalog_raises_without_size(self, minimal_stream_config):
        model = StreamModel(minimal_stream_config)
        with pytest.raises(ValueError, match="size"):
            model.complete_catalog(catalog=None, size=None, verbose=False)