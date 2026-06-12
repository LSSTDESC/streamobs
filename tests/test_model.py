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

@pytest.mark.model
class TestDensityModel:
    """Unit tests for DensityModel (phi1 sampler)."""

    def _verify_sampling(self, model):
        """Helper to verify that sampling produces a 1D array of the right length."""
        samples = model.sample(N)
        assert samples.shape == (N,), "Sample array must have length N"
        assert np.issubdtype(samples.dtype, np.floating), "Samples must be floats"
        assert np.all(np.isfinite(samples)), "Samples must be finite numbers"
        assert np.all(~np.isnan(samples)), "Samples must not contain NaNs"
        return samples

    def test_uniform_sample(self):
        cfg = {"type": "uniform", "xmin": -10.0, "xmax": 10.0}
        model = DensityModel(cfg)
        samples = self._verify_sampling(model)
        assert np.all(samples>cfg["xmin"]) and np.all(samples<cfg["xmax"]), "Uniform samples must be within [xmin, xmax]"


    def test_gaussian_sample_shape(self):
        cfg = {"type": "gaussian", "mu": 0.0, "sigma": 1.0}
        model = DensityModel(cfg)
        samples = self._verify_sampling(model)




# ---------------------------------------------------------------------------
# TrackModel
# ---------------------------------------------------------------------------
@pytest.mark.model
class TestTrackModel:
    """Unit tests for TrackModel (phi2 sampler given phi1)."""

    def _verify_sampling(self, model, phi1):
        """Helper to verify that sampling produces a 1D array of the right length."""
        samples = model.sample(phi1)
        assert samples.shape == (N,), "Sample array must have length N"
        assert np.issubdtype(samples.dtype, np.floating), "Samples must be floats"
        assert np.all(np.isfinite(samples)), "Samples must be finite numbers"
        assert np.all(~np.isnan(samples)), "Samples must not contain NaNs"
        return samples
    
    def test_constant_track(self):
        cfg = {
            "center": {"type": "constant", "value": 0.0},
            "spread": {"type": "constant", "value": 0.2},
            "sampler": "gaussian"
        }
        model = TrackModel(cfg)
        phi1 = np.linspace(-10, 10, N)
        samples = self._verify_sampling(model, phi1)
        assert np.all(np.abs(samples) < 10*cfg["spread"]['value']), "Samples should be within a few sigma of the center"

    def test_sinusoidal_track(self):
        cfg = {
            "center": {"type": "sinusoid", "amplitude": 0.5, "period": 2.},
            "spread": {"type": "constant", "value": 0.2},
            "sampler": "gaussian"
        }
        model = TrackModel(cfg)
        phi1 = np.linspace(-10, 10, N)
        samples = self._verify_sampling(model, phi1)
        expected_center = cfg["center"]["amplitude"] * np.sin(phi1*2*np.pi/cfg["center"]["period"])
        assert np.all(np.abs(samples - expected_center) < 10*cfg["spread"]['value']), "Samples should be within a few sigma of the sinusoidal center"


# ---------------------------------------------------------------------------
# StreamModel — full
# ---------------------------------------------------------------------------

@pytest.mark.model
class TestStreamModelFull:
    """Tests for StreamModel when having a complete config"""
    def test_full_model(self, stream_config_with_distance):
        """Test that StreamModel can be instantiated and sampled with a full config."""
        model = StreamModel(stream_config_with_distance)
        samples = model.sample(N)
        assert isinstance(samples, pd.DataFrame), "Samples should be returned as a DataFrame"
        expected_columns = {"phi1", "phi2", "dist", "mag_g", "mag_r"} # Not adding mu1, mu2, rv since not implemented yet
        assert expected_columns.issubset(samples.columns), f"Sample DataFrame should contain columns {expected_columns}"
        assert len(samples) == N, f"Sample DataFrame should have {N} rows"
        for col in expected_columns:
            assert np.issubdtype(samples[col].dtype, np.floating), f"Column {col} should contain floats"
            assert np.all(np.isfinite(samples[col])), f"Column {col} should contain finite numbers"
            assert np.all(~np.isnan(samples[col])), f"Column {col} should not contain NaNs"