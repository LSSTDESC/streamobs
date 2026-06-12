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

    pass  # No tests yet — just a placeholder to mark the section

# ---------------------------------------------------------------------------
# Distance modulus model
# ---------------------------------------------------------------------------

@pytest.mark.model
class TestDistanceModulus:
    """Tests for Distance modulus"""
    pass  # No tests yet — just a placeholder to mark the section


# ---------------------------------------------------------------------------
# StreamModel — full
# ---------------------------------------------------------------------------

@pytest.mark.model
class TestStreamModelFull:
    """Tests for StreamModel when having a complete config"""
    pass  # No tests yet — just a placeholder to mark the section
