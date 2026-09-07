"""
tests/test_match_filter.py
==========================
Tests for ``streamobs.match_filter``.

Integration test: inject a mock stream with known isochrone parameters, then
verify that a match filter built with the same parameters selects the majority
of detected stream stars.
"""

import numpy as np
import pytest

from streamobs.match_filter import build_match_filter, is_in_match_filter

# ---------------------------------------------------------------------------
# Smoke tests — return types and shapes
# ---------------------------------------------------------------------------


@pytest.mark.match_filter
class TestMatchFilterAPI:
    def test_build_match_filter_return_shape(self):
        polygon = build_match_filter(
            distance_modulus=16.8, age=12.0, metallicity=0.0006
        )
        assert isinstance(polygon, np.ndarray), "polygon_vertices must be an ndarray"
        assert polygon.ndim == 2 and polygon.shape[1] == 2, "shape must be (N, 2)"
        assert polygon.shape[0] >= 3, "polygon must have at least 3 vertices"
        assert np.all(np.isfinite(polygon)), "polygon vertices must be finite"

    def test_is_in_match_filter_return_type(self):
        polygon = build_match_filter(
            distance_modulus=16.8, age=12.0, metallicity=0.0006
        )
        mag_g = np.array([20.0, 22.0, 24.0])
        mag_r = np.array([19.5, 21.5, 23.5])
        mask = is_in_match_filter(mag_g, mag_r, polygon_vertices=polygon)
        assert isinstance(mask, np.ndarray), "output must be an ndarray"
        assert mask.dtype == bool, "output must be boolean"
        assert len(mask) == len(mag_g), "output length must match input length"

    def test_is_in_match_filter_handles_nan(self):
        polygon = build_match_filter(
            distance_modulus=16.8, age=12.0, metallicity=0.0006
        )
        mag_g = np.array([20.0, float("nan"), 24.0])
        mag_r = np.array([19.5, 21.5, float("nan")])
        mask = is_in_match_filter(mag_g, mag_r, polygon_vertices=polygon)
        assert not mask[1], "NaN g should not be selected"
        assert not mask[2], "NaN r should not be selected"


# ---------------------------------------------------------------------------
# Integration test — match filter selects most stream stars
# ---------------------------------------------------------------------------


@pytest.mark.match_filter
class TestMatchFilterIntegration:
    def test_match_filter_selects_stream_stars(
        self, mock_injector, stream_catalog, stream_config_with_distance
    ):
        """
        A match filter built with the same isochrone as the stream should
        select at least 80% of detected stream stars.
        """
        injected = mock_injector.inject(
            stream_catalog, seed=42, perfect_galstarsep=True
        )
        detected = injected[injected["lsst_yr4_flag_observed"]].copy()

        if len(detected) < 30:
            pytest.skip(
                f"Too few detected stars ({len(detected)}) to test match filter."
            )

        iso = stream_config_with_distance["isochrone"]
        dm = stream_config_with_distance["distance_modulus"]["center"]["value"]

        # Geometry/isochrone check with the constant red cap: the polygon itself
        # (blue edge, magnitude range, isochrone match) must select most
        # detected stream stars.
        polygon = build_match_filter(
            distance_modulus=dm,
            age=iso["age"],
            metallicity=iso["z"],
            red_cap_mode="constant",
        )
        mask = is_in_match_filter(
            detected["lsst_yr4_g_obs"],
            detected["lsst_yr4_r_obs"],
            polygon_vertices=polygon,
        )
        fraction = mask.sum() / len(detected)
        assert fraction > 0.80, (
            f"Match filter selected only {fraction*100:.1f}% of detected stream "
            f"stars (expected > 80%). This may indicate the filter is too narrow "
            f"or the isochrone parameters do not match."
        )

        # The default error-shrunk red cap intentionally trades faint-red
        # completeness for background suppression (~7% better stream S/N on the
        # DC2 background): the cap tightens as errors grow, cutting faint red
        # MS stars near/past the maglim. Assert the documented behavior: a
        # strict subset of the constant-cap selection, but still the majority
        # of detected stream stars.
        polygon_default = build_match_filter(
            distance_modulus=dm, age=iso["age"], metallicity=iso["z"]
        )
        mask_default = is_in_match_filter(
            detected["lsst_yr4_g_obs"],
            detected["lsst_yr4_r_obs"],
            polygon_vertices=polygon_default,
        )
        assert not np.any(mask_default & ~mask), (
            "error_shrunk selection should be a subset of the constant-cap " "selection"
        )
        fraction_default = mask_default.sum() / len(detected)
        assert fraction_default > 0.60, (
            f"Default (error_shrunk) filter selected only "
            f"{fraction_default*100:.1f}% of detected stream stars (expected "
            f"> 60%; it trades faint-red completeness for S/N)."
        )

    def test_match_filter_rejects_offset_stars(self, stream_config_with_distance):
        """
        Stars shifted 2 mag blueward of the isochrone should almost all fall
        outside the match filter, showing the filter is not trivially permissive.
        """
        iso = stream_config_with_distance["isochrone"]
        dm = stream_config_with_distance["distance_modulus"]["center"]["value"]
        polygon = build_match_filter(
            distance_modulus=dm, age=iso["age"], metallicity=iso["z"]
        )

        rng = np.random.default_rng(0)
        n = 500
        # Sample magnitudes spanning the filter's magnitude range
        mag_g = rng.uniform(dm + 1.0, dm + 9.0, n)
        # Shift 2 mag blueward of any reasonable isochrone color
        mag_r = mag_g + 2.5

        mask = is_in_match_filter(mag_g, mag_r, polygon_vertices=polygon)
        fraction = mask.sum() / n
        assert fraction < 0.05, (
            f"Match filter selected {fraction*100:.1f}% of stars shifted 2 mag "
            f"blueward (expected < 5%)."
        )
