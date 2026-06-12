"""
tests/test_observed.py
=======================
Tests for ``streamobs.observed.StreamInjector``.

All tests use the ``mock_survey`` / ``mock_injector`` fixtures from
``conftest.py``, so no real survey files are required.

The injection pipeline is tested end-to-end by calling
``StreamInjector.inject()`` with controlled input DataFrames and verifying
the output columns, dtypes, and flag semantics.

Structure
---------
- ``TestStreamInjectorInit``          — construction from Survey object / string
- ``TestStaticHelpers``               — magToFlux / fluxToMag / getFluxError
- ``TestSampleMeasuredMagnitudes``    — photometric noise sampling
- ``TestDetectFlag``                  — detection threshold logic
- ``TestInjectOutputColumns``         — inject() column existence and types
- ``TestInjectFlagBehaviour``         — flag_observed semantics
- ``TestInjectReproducibility``       — seed / RNG determinism
- ``TestInjectDifferentDatasets``     — varying input shapes and column sets
"""

import numpy as np
import pandas as pd
import pytest

from streamobs.observed import StreamInjector


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

BANDS = ["r", "g"]


def _bright_catalog(n=40, ra_center=45.0, dec_center=-10.0):
    """Catalog with stars well above the mock survey's 25-mag limit → all detected."""
    rng = np.random.default_rng(42)
    return pd.DataFrame({
        "ra":    np.full(n, ra_center) + rng.uniform(-1, 1, n),
        "dec":   np.full(n, dec_center) + rng.uniform(-1, 1, n),
        "mag_g": np.full(n, 20.0),   # 5 mag above limit
        "mag_r": np.full(n, 20.0),
    })


def _faint_catalog(n=40, ra_center=45.0, dec_center=-10.0):
    """Catalog with stars well below the 25-mag limit → mostly undetected."""
    rng = np.random.default_rng(99)
    return pd.DataFrame({
        "ra":    np.full(n, ra_center) + rng.uniform(-1, 1, n),
        "dec":   np.full(n, dec_center) + rng.uniform(-1, 1, n),
        "mag_g": np.full(n, 28.0),   # 3 mag below limit
        "mag_r": np.full(n, 28.0),
    })


# ---------------------------------------------------------------------------
# Init
# ---------------------------------------------------------------------------

class TestStreamInjectorInit:

    def test_init_from_survey_object(self, mock_survey):
        inj = StreamInjector(mock_survey)
        assert inj.survey is mock_survey

    def test_init_invalid_type_raises(self):
        with pytest.raises(ValueError, match="survey must be"):
            StreamInjector(12345)

    def test_last_gc_frame_initially_none(self, mock_injector):
        assert mock_injector._last_gc_frame is None


# ---------------------------------------------------------------------------
# Static helpers
# ---------------------------------------------------------------------------

class TestStaticHelpers:

    def test_mag_to_flux_scalar(self):
        flux = StreamInjector.magToFlux(0.0)
        assert flux == pytest.approx(3631.0)

    def test_flux_to_mag_inverse(self):
        mag = 22.5
        assert StreamInjector.fluxToMag(StreamInjector.magToFlux(mag)) == pytest.approx(mag)

    def test_get_flux_error_positive(self):
        err = StreamInjector.getFluxError(22.0, 0.01)
        assert float(err) > 0

    def test_mag_flux_roundtrip_array(self):
        mags = np.array([18.0, 20.0, 22.0, 24.0])
        recovered = StreamInjector.fluxToMag(StreamInjector.magToFlux(mags))
        np.testing.assert_allclose(recovered, mags, rtol=1e-6)


# ---------------------------------------------------------------------------
# sample_measured_magnitudes
# ---------------------------------------------------------------------------

class TestSampleMeasuredMagnitudes:

    def test_output_shape(self, mock_injector):
        n = 100
        mag_true = np.full(n, 22.0)
        mag_err  = np.full(n, 0.05)
        mag_obs  = mock_injector.sample_measured_magnitudes(mag_true, mag_err, seed=0)
        assert len(mag_obs) == n

    def test_bright_stars_not_bad_mag(self, mock_injector):
        """Very bright stars should nearly always produce valid (non-BAD_MAG) fluxes."""
        n = 200
        mag_true = np.full(n, 16.5)   # just above saturation, large flux
        mag_err  = np.full(n, 0.001)
        mag_obs  = mock_injector.sample_measured_magnitudes(mag_true, mag_err, seed=7)
        bad = np.sum(mag_obs == "BAD_MAG")
        assert bad == 0, f"Bright stars produced {bad} BAD_MAG values"

    def test_faint_stars_may_be_bad_mag(self, mock_injector):
        """Very faint stars (large errors) may occasionally return BAD_MAG."""
        n = 500
        mag_true = np.full(n, 30.0)
        mag_err  = np.full(n, 5.0)    # enormous error → many negative fluxes
        mag_obs  = mock_injector.sample_measured_magnitudes(mag_true, mag_err, seed=3)
        bad = np.sum(mag_obs == "BAD_MAG")
        assert bad > 0, "Expected some BAD_MAG for stars with huge errors"

    def test_reproducible_with_seed(self, mock_injector):
        mag_true = np.linspace(20.0, 24.0, 50)
        mag_err  = np.full(50, 0.05)
        obs1 = mock_injector.sample_measured_magnitudes(mag_true, mag_err, seed=123)
        obs2 = mock_injector.sample_measured_magnitudes(mag_true, mag_err, seed=123)
        np.testing.assert_array_equal(obs1, obs2)

    def test_different_seeds_differ(self, mock_injector):
        mag_true = np.linspace(20.0, 24.0, 50)
        mag_err  = np.full(50, 0.05)
        obs1 = mock_injector.sample_measured_magnitudes(mag_true, mag_err, seed=1)
        obs2 = mock_injector.sample_measured_magnitudes(mag_true, mag_err, seed=2)
        assert not np.array_equal(obs1, obs2)


# ---------------------------------------------------------------------------
# detect_flag
# ---------------------------------------------------------------------------

class TestDetectFlag:

    def _pixel_array(self, n, nside=64):
        """Return an array of valid pixel indices for the mock survey's nside."""
        import healpy as hp
        npix = hp.nside2npix(nside)
        return np.random.default_rng(0).integers(0, npix, n)

    def test_output_is_boolean_array(self, mock_injector):
        pix  = self._pixel_array(100)
        mags = np.full(100, 22.0)
        flag = mock_injector.detect_flag(pix, mag=mags, band="r", seed=0)
        assert flag.dtype == bool

    def test_output_shape(self, mock_injector):
        n    = 80
        pix  = self._pixel_array(n)
        mags = np.full(n, 22.0)
        flag = mock_injector.detect_flag(pix, mag=mags, band="r", seed=0)
        assert flag.shape == (n,)

    def test_bright_stars_mostly_detected(self, mock_injector):
        """Stars well within the limit (bright end) should almost all be flagged True."""
        n    = 500
        pix  = self._pixel_array(n)
        mags = np.full(n, 20.0)   # 5 mag above 25-mag limit
        flag = mock_injector.detect_flag(pix, mag=mags, band="r", seed=42)
        assert flag.mean() > 0.95, "Expected >95% detection for bright stars"

    def test_faint_stars_mostly_undetected(self, mock_injector):
        """Stars far below the limit should nearly never be flagged True."""
        n    = 500
        pix  = self._pixel_array(n)
        mags = np.full(n, 28.0)   # 3 mag below limit
        flag = mock_injector.detect_flag(pix, mag=mags, band="r", seed=99)
        assert flag.mean() < 0.05, "Expected <5% detection for very faint stars"


# ---------------------------------------------------------------------------
# inject — output columns
# ---------------------------------------------------------------------------

class TestInjectOutputColumns:

    def test_inject_adds_mag_obs_columns(self, mock_injector):
        data = _bright_catalog()
        result = mock_injector.inject(data, bands=BANDS, seed=0, verbose=False)
        assert "mag_g_obs" in result.columns
        assert "mag_r_obs" in result.columns

    def test_inject_adds_magerr_columns(self, mock_injector):
        data = _bright_catalog()
        result = mock_injector.inject(data, bands=BANDS, seed=0, verbose=False)
        assert "magerr_g" in result.columns
        assert "magerr_r" in result.columns

    def test_inject_adds_flag_observed(self, mock_injector):
        data = _bright_catalog()
        result = mock_injector.inject(data, bands=BANDS, seed=0, verbose=False)
        assert "flag_observed" in result.columns

    def test_inject_no_flag_perfect_by_default(self, mock_injector):
        data = _bright_catalog()
        result = mock_injector.inject(data, bands=BANDS, seed=0, verbose=False)
        assert "flag_perfect_galstarsep" not in result.columns

    def test_inject_flag_perfect_when_requested(self, mock_injector):
        data = _bright_catalog()
        result = mock_injector.inject(
            data, bands=BANDS, seed=0, verbose=False, perfect_galstarsep=True
        )
        assert "flag_perfect_galstarsep" in result.columns

    def test_inject_preserves_row_count(self, mock_injector):
        data = _bright_catalog(n=60)
        result = mock_injector.inject(data, bands=BANDS, seed=0, verbose=False)
        assert len(result) == 60

    def test_inject_preserves_original_columns(self, mock_injector):
        data = _bright_catalog()
        result = mock_injector.inject(data, bands=BANDS, seed=0, verbose=False)
        for col in ("ra", "dec", "mag_g", "mag_r"):
            assert col in result.columns


# ---------------------------------------------------------------------------
# inject — flag behaviour
# ---------------------------------------------------------------------------

class TestInjectFlagBehaviour:

    def test_bright_catalog_high_detection_rate(self, mock_injector):
        data   = _bright_catalog(n=200)
        result = mock_injector.inject(data, bands=BANDS, seed=42, verbose=False)
        rate   = result["flag_observed"].mean()
        assert rate > 0.5, f"Expected >50% detection for bright catalog, got {rate:.2%}"

    def test_faint_catalog_low_detection_rate(self, mock_injector):
        data   = _faint_catalog(n=200)
        result = mock_injector.inject(data, bands=BANDS, seed=42, verbose=False)
        rate   = result["flag_observed"].mean()
        assert rate < 0.5, f"Expected <50% detection for faint catalog, got {rate:.2%}"

    def test_flag_observed_is_boolean(self, mock_injector):
        data   = _bright_catalog()
        result = mock_injector.inject(data, bands=BANDS, seed=0, verbose=False)
        assert result["flag_observed"].dtype == bool

    def test_magerr_positive(self, mock_injector):
        data   = _bright_catalog()
        result = mock_injector.inject(data, bands=BANDS, seed=0, verbose=False)
        assert (result["magerr_r"] > 0).all()
        assert (result["magerr_g"] > 0).all()


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------

class TestInjectReproducibility:

    def test_same_seed_same_flags(self, mock_injector):
        data = _bright_catalog(n=80)
        r1   = mock_injector.inject(data.copy(), bands=BANDS, seed=7, verbose=False)
        r2   = mock_injector.inject(data.copy(), bands=BANDS, seed=7, verbose=False)
        np.testing.assert_array_equal(
            r1["flag_observed"].values,
            r2["flag_observed"].values,
        )

    def test_different_seeds_give_different_flags(self, mock_injector):
        data = _bright_catalog(n=200)
        r1   = mock_injector.inject(data.copy(), bands=BANDS, seed=1, verbose=False)
        r2   = mock_injector.inject(data.copy(), bands=BANDS, seed=2, verbose=False)
        # Flags may occasionally agree but should differ for a large catalog
        assert not np.array_equal(
            r1["flag_observed"].values, r2["flag_observed"].values
        ), "Different seeds should produce different flags"


# ---------------------------------------------------------------------------
# Varying datasets
# ---------------------------------------------------------------------------

class TestInjectDifferentDatasets:

    def test_inject_single_band_r(self, mock_injector, sample_catalog_radec):
        result = mock_injector.inject(
            sample_catalog_radec, bands=["r"], seed=0, verbose=False
        )
        assert "mag_r_obs" in result.columns
        assert "magerr_r"  in result.columns

    def test_inject_single_band_g(self, mock_injector, sample_catalog_radec):
        result = mock_injector.inject(
            sample_catalog_radec, bands=["g"], seed=0, verbose=False,
            detection_mag_cut=[]  # skip SNR cut on g when r is absent
        )
        assert "mag_g_obs" in result.columns

    def test_inject_n1_catalog(self, mock_injector):
        """Injection must work for a catalog with a single star."""
        data = pd.DataFrame({
            "ra":    [45.0],
            "dec":   [-10.0],
            "mag_g": [22.0],
            "mag_r": [21.5],
        })
        result = mock_injector.inject(data, bands=BANDS, seed=0, verbose=False)
        assert len(result) == 1

    def test_inject_large_catalog(self, mock_injector):
        """Injection must scale to larger catalogs without error."""
        rng = np.random.default_rng(55)
        n   = 2_000
        data = pd.DataFrame({
            "ra":    rng.uniform(30, 60, n),
            "dec":   rng.uniform(-20, 0, n),
            "mag_g": rng.uniform(21, 24, n),
            "mag_r": rng.uniform(20, 23, n),
        })
        result = mock_injector.inject(data, bands=BANDS, seed=0, verbose=False)
        assert len(result) == n
        assert "flag_observed" in result.columns

    def test_inject_with_extra_columns_preserved(self, mock_injector):
        """Unrelated columns in the catalog should pass through unchanged."""
        data = _bright_catalog(n=30)
        data["stream_id"] = np.arange(30)
        result = mock_injector.inject(data, bands=BANDS, seed=0, verbose=False)
        assert "stream_id" in result.columns
        np.testing.assert_array_equal(result["stream_id"].values, np.arange(30))

    def test_inject_unsupported_band_raises(self, mock_injector):
        data = _bright_catalog()
        data["mag_z"] = 22.0
        with pytest.raises(ValueError):
            mock_injector.inject(data, bands=["z"], seed=0, verbose=False)


# ---------------------------------------------------------------------------
# Mask cache
# ---------------------------------------------------------------------------

class TestMaskCache:

    def test_clear_mask_cache(self):
        StreamInjector.clear_mask_cache()
        assert StreamInjector.mask_cache == {}

    def test_list_cached_masks_empty(self):
        StreamInjector.clear_mask_cache()
        result = StreamInjector.list_cached_masks()
        assert result == []