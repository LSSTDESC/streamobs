"""
tests/test_model.py
===================
Tests for ``streamobs.model.StreamModel`` and its sub-models.

Focus: verify that every quantity the model can produce is actually sampled
and has the right shape / dtype, without requiring the optional ugali.
"""

import numpy as np
import pandas as pd
import pytest

from streamobs.model import DensityModel, StreamModel, TrackModel

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
        assert np.all(samples > cfg["xmin"]) and np.all(
            samples < cfg["xmax"]
        ), "Uniform samples must be within [xmin, xmax]"

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
            "sampler": "gaussian",
        }
        model = TrackModel(cfg)
        phi1 = np.linspace(-10, 10, N)
        samples = self._verify_sampling(model, phi1)
        assert np.all(
            np.abs(samples) < 10 * cfg["spread"]["value"]
        ), "Samples should be within a few sigma of the center"

    def test_sinusoidal_track(self):
        cfg = {
            "center": {"type": "sinusoid", "amplitude": 0.5, "period": 2.0},
            "spread": {"type": "constant", "value": 0.2},
            "sampler": "gaussian",
        }
        model = TrackModel(cfg)
        phi1 = np.linspace(-10, 10, N)
        samples = self._verify_sampling(model, phi1)
        expected_center = cfg["center"]["amplitude"] * np.sin(
            phi1 * 2 * np.pi / cfg["center"]["period"]
        )
        assert np.all(
            np.abs(samples - expected_center) < 10 * cfg["spread"]["value"]
        ), "Samples should be within a few sigma of the sinusoidal center"


# ---------------------------------------------------------------------------
# StreamModel — full
# ---------------------------------------------------------------------------


@pytest.mark.model
class TestStreamModelFull:
    """Tests for StreamModel when having a complete config"""

    def _verify_catalogue_content(self, catalog, expected_columns):
        """Helper to verify that a completed catalog contains the expected columns with valid data."""
        assert expected_columns.issubset(
            catalog.columns
        ), f"Catalog should contain columns {expected_columns}"
        for col in expected_columns:
            assert np.issubdtype(
                catalog[col].dtype, np.floating
            ), f"Column {col} should contain floats"
            assert np.all(
                np.isfinite(catalog[col])
            ), f"Column {col} should contain finite numbers"
            assert np.all(
                ~np.isnan(catalog[col])
            ), f"Column {col} should not contain NaNs"

    def test_full_model(self, stream_config_with_distance):
        """Test that StreamModel can be instantiated and sampled with a full config."""
        model = StreamModel(stream_config_with_distance)
        samples = model.sample(N)
        assert isinstance(
            samples, pd.DataFrame
        ), "Samples should be returned as a DataFrame"
        expected_columns = {
            "phi1",
            "phi2",
            "dist",
            "lsst_g_true",
            "lsst_r_true",
        }  # Not adding mu1, mu2, rv since not implemented yet
        self._verify_catalogue_content(samples, expected_columns)

    def test_complete_catalog(self, sample_catalog_phi, stream_config_with_distance):
        """Test that complete_catalog produces a catalog with the expected columns and valid data."""
        model = StreamModel(stream_config_with_distance)
        completed_catalog = model.complete_catalog(catalog=sample_catalog_phi)
        expected_columns = {
            "phi1",
            "phi2",
            "dist",
            "lsst_g_true",
            "lsst_r_true",
        }  # Not adding mu1, mu2, rv since not implemented yet
        self._verify_catalogue_content(completed_catalog, expected_columns)
        assert len(completed_catalog) == len(
            sample_catalog_phi
        ), f"Completed catalog should have the same number of rows as the input catalog"

        # Verify that I can add a targeted column (e.g. dist) to the input catalog and complete the rest
        partial_catalog = completed_catalog.drop(
            columns=["lsst_r_true", "dist", "lsst_g_true"]
        ).reset_index(drop=True)
        completed_catalog = model.complete_catalog(
            catalog=partial_catalog,
            columns_to_add=["dist"],
        )
        expected_columns = {
            "phi1",
            "phi2",
            "dist",
        }  # Not adding mu1, mu2, rv since not implemented yet
        self._verify_catalogue_content(completed_catalog, expected_columns)
        assert len(completed_catalog) == len(
            partial_catalog
        ), f"Completed catalog should have the same number of rows as the input catalog"
        assert (
            "lsst_r_true" not in completed_catalog.columns
        ), "Column 'lsst_r_true' should not be added when not requested"


# ---------------------------------------------------------------------------
# StreamModel.complete_catalog — permutations of input columns / dist source
# ---------------------------------------------------------------------------


@pytest.mark.model
class TestCompleteCatalogPermutations:
    """Exercise complete_catalog across the ways a catalog can be partially filled.

    Covers: which columns are supplied (empty frame, phi-only, ra/dec+dist),
    where the distance comes from (distance_modulus model vs. a directly supplied
    scalar/vector ``dist``), and the preserve-existing-values contract (whole
    columns and individual NaN rows are never overwritten).
    """

    MAGS = {"lsst_g_true", "lsst_r_true"}

    def test_empty_frame_fills_all_model_columns(self, stream_config_with_distance):
        """size=N with no catalog -> geometry, dist, and both bands are filled."""
        model = StreamModel(stream_config_with_distance)
        out = model.complete_catalog(catalog=None, size=12, verbose=False)
        assert len(out) == 12
        assert ({"phi1", "phi2", "dist"} | self.MAGS).issubset(out.columns)
        assert out[list(self.MAGS)].notna().all().all()

    def test_phi_only_fills_dist_then_mags(
        self, sample_catalog_phi, stream_config_with_distance
    ):
        """phi1/phi2 present -> dist sampled from the model, then magnitudes."""
        model = StreamModel(stream_config_with_distance)
        out = model.complete_catalog(catalog=sample_catalog_phi.copy(), verbose=False)
        assert ({"dist"} | self.MAGS).issubset(out.columns)
        assert out[list(self.MAGS)].notna().all().all()

    def test_radec_plus_dist_fills_mags_without_phi(self, stream_config_with_distance):
        """With dist already present, magnitudes fill even when phi1 is absent."""
        model = StreamModel(stream_config_with_distance)
        df = pd.DataFrame(
            {"ra": [10.0, 11.0, 12.0], "dec": [-1.0, 0.0, 1.0], "dist": [16.0] * 3}
        )
        out = model.complete_catalog(
            catalog=df, columns_to_add=list(self.MAGS), verbose=False
        )
        assert self.MAGS.issubset(out.columns)
        assert out[list(self.MAGS)].notna().all().all()
        assert "phi1" not in out.columns

    def test_existing_band_preserved_when_filling_other(
        self, stream_config_with_distance
    ):
        """Providing one band and requesting both leaves the provided one intact."""
        model = StreamModel(stream_config_with_distance)
        g = np.array([20.0, 21.0, 22.0])
        df = pd.DataFrame({"dist": [16.0] * 3, "lsst_g_true": g.copy()})
        out = model.complete_catalog(
            catalog=df, columns_to_add=list(self.MAGS), verbose=False
        )
        assert np.allclose(out["lsst_g_true"].to_numpy(), g), "present band overwritten"
        assert out["lsst_r_true"].notna().all(), "missing band not filled"

    def test_present_columns_skip_sampling(self, stream_config_with_distance):
        """Both bands present -> values are returned untouched."""
        model = StreamModel(stream_config_with_distance)
        g = np.array([20.0, 21.0])
        r = np.array([19.0, 19.5])
        df = pd.DataFrame(
            {"dist": [16.0, 16.0], "lsst_g_true": g.copy(), "lsst_r_true": r.copy()}
        )
        out = model.complete_catalog(
            catalog=df, columns_to_add=list(self.MAGS), verbose=False
        )
        assert np.allclose(out["lsst_g_true"].to_numpy(), g)
        assert np.allclose(out["lsst_r_true"].to_numpy(), r)

    def test_partial_rows_only_missing_filled(self, stream_config_with_distance):
        """A band with some NaN rows keeps its finite rows; only NaNs are filled."""
        model = StreamModel(stream_config_with_distance)
        g = np.array([20.0, np.nan, 22.0])
        df = pd.DataFrame({"dist": [16.0] * 3, "lsst_g_true": g.copy()})
        out = model.complete_catalog(
            catalog=df, columns_to_add=["lsst_g_true"], verbose=False
        )
        filled = out["lsst_g_true"].to_numpy()
        assert filled[0] == 20.0 and filled[2] == 22.0, "finite rows overwritten"
        assert np.isfinite(filled[1]), "NaN row not filled"

    def test_dist_scalar_broadcast_without_distance_model(self, minimal_stream_config):
        """A scalar `dist` lets mags fill with no distance_modulus model / no phi1."""
        model = StreamModel(minimal_stream_config)  # no distance_modulus section
        df = pd.DataFrame({"ra": [10.0, 11.0], "dec": [0.0, 1.0]})
        out = model.complete_catalog(
            catalog=df,
            columns_to_add=["dist"] + list(self.MAGS),
            dist=16.5,
            verbose=False,
        )
        assert np.allclose(out["dist"].to_numpy(), 16.5)
        assert out[list(self.MAGS)].notna().all().all()

    def test_dist_vector_assigned_per_row(self, minimal_stream_config):
        """A per-row `dist` vector is assigned row-wise."""
        model = StreamModel(minimal_stream_config)
        df = pd.DataFrame({"ra": [10.0, 11.0, 12.0], "dec": [0.0, 1.0, 2.0]})
        dvec = np.array([15.0, 16.0, 17.0])
        out = model.complete_catalog(
            catalog=df,
            columns_to_add=["dist"] + list(self.MAGS),
            dist=dvec,
            verbose=False,
        )
        assert np.allclose(out["dist"].to_numpy(), dvec)
        assert out[list(self.MAGS)].notna().all().all()

    def test_dist_vector_wrong_length_raises(self, minimal_stream_config):
        model = StreamModel(minimal_stream_config)
        df = pd.DataFrame({"ra": [10.0, 11.0, 12.0], "dec": [0.0, 1.0, 2.0]})
        with pytest.raises(ValueError):
            model.complete_catalog(
                catalog=df,
                columns_to_add=["dist"],
                dist=np.array([1.0, 2.0]),
                verbose=False,
            )

    def test_dist_overrides_distance_model(self, stream_config_with_distance):
        """When given, `dist` is used instead of the configured distance model."""
        model = StreamModel(stream_config_with_distance)  # model would give 16.8
        df = pd.DataFrame({"phi1": [0.0, 1.0], "phi2": [0.0, 0.0]})
        out = model.complete_catalog(
            catalog=df, columns_to_add=["dist"], dist=20.0, verbose=False
        )
        assert np.allclose(out["dist"].to_numpy(), 20.0)

    def test_mags_without_dist_or_phi_raises(self, minimal_stream_config):
        """No distance_modulus model, no dist, no phi1 -> cannot fill magnitudes."""
        model = StreamModel(minimal_stream_config)
        df = pd.DataFrame({"ra": [10.0, 11.0], "dec": [0.0, 1.0]})
        with pytest.raises(ValueError):
            model.complete_catalog(
                catalog=df, columns_to_add=list(self.MAGS), verbose=False
            )
