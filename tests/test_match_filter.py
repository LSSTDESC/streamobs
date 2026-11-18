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
from matplotlib.path import Path

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


# ---------------------------------------------------------------------------
# Colour across surveys (DES g - Euclid VIS)
# ---------------------------------------------------------------------------

_DES_EUCLID = dict(band_1="g", survey_1="des", band_2="VIS", survey_2="euclid")


@pytest.mark.match_filter
class TestCrossSurveyMatchFilter:
    def test_same_survey_spelled_out_is_unchanged(self):
        """survey_1 == survey_2 is exactly the single-survey filter."""
        kw = dict(distance_modulus=16.8, age=12.0, metallicity=0.0006)
        assert np.array_equal(
            build_match_filter(**kw, survey="des", band_1="g", band_2="r"),
            build_match_filter(
                **kw, band_1="g", survey_1="des", band_2="r", survey_2="des"
            ),
        )

    def test_selects_the_cross_survey_locus(self, stream_config_with_distance):
        """Noiseless DES g / Euclid VIS stream stars lie on the DES g - Euclid
        VIS locus: a narrow (+/- 0.02 mag) polygon around it selects them, one
        around DES g - r does not."""
        from streamobs.model import StreamModel

        dm = stream_config_with_distance["distance_modulus"]["center"]["value"]
        isochrone = {
            "name": "Marigo2017",
            "age": 12.0,
            "z": 0.0006,
            "surveys": {"des": {"bands": ["g"]}, "euclid": {"bands": ["VIS"]}},
        }
        stars = StreamModel(
            {**stream_config_with_distance, "isochrone": isochrone}
        ).sample(5000, seed=1)
        # Main-sequence stars, below the RGB clip and within the polygon range
        g, vis = stars["des_g_true"], stars["euclid_VIS_true"]
        main_sequence = (g > dm + 3.5) & (g < dm + 8.5)
        g, vis = g[main_sequence], vis[main_sequence]
        assert len(g) > 100

        # The default polygon is wider than the g-VIS / g-r difference; a
        # narrow one checks that band_2 is paired with band_1 at the same mass.
        narrow = dict(
            age=12.0,
            metallicity=0.0006,
            color_spread=[0.02, 0.02],
            error_multiplier=[0.0, 0.0],
            distance_modulus_spread=0.0,
        )
        polygon = build_match_filter(dm, **narrow, **_DES_EUCLID)
        inside = is_in_match_filter(g, vis, polygon_vertices=polygon)
        assert inside.mean() > 0.95, f"only {inside.mean():.1%} selected"
        g_r = build_match_filter(dm, **narrow, survey="des", band_1="g", band_2="r")
        wrong = is_in_match_filter(g, vis, polygon_vertices=g_r)
        assert wrong.mean() < 0.5, f"g-r polygon selected {wrong.mean():.1%}"


# ---------------------------------------------------------------------------
# Isochrone locus: only the points before ugali's ``index``
# ---------------------------------------------------------------------------


class _FakeIsochrone:
    """An ugali-like isochrone: colour, magnitude, and the index of its late stages."""

    color = np.array([0.1, 0.2, 0.3, -0.4, -0.5])
    mag = np.array([4.0, 5.0, 6.0, 10.0, 12.0])
    index = 3  # the last two points: late stages (e.g. a white dwarf sequence)


@pytest.mark.match_filter
class TestIsochroneLocus:
    def test_same_survey_stops_at_index(self, monkeypatch):
        """One survey: the points before ``index`` only, as with two surveys."""
        from streamobs import match_filter

        monkeypatch.setattr(match_filter, "isochrone_factory", lambda *a, **k: _FakeIsochrone())
        color, mag = match_filter._isochrone_locus("Dotter2016", 13.0, 1e-4, "VIS", "euclid", "H", "euclid")
        np.testing.assert_array_equal(color, _FakeIsochrone.color[:3])
        np.testing.assert_array_equal(mag, _FakeIsochrone.mag[:3])

    def test_dotter_white_dwarfs_left_out(self):
        """A Dotter2016 Euclid isochrone with a white dwarf sequence (13 Gyr,
        Z = 1e-4): its match filter is as narrow as without it, and its locus
        has none of the late stages."""
        from ugali.analysis.isochrone import factory as isochrone_factory

        from streamobs.match_filter import _isochrone_locus

        try:
            isochrone = isochrone_factory("Dotter2016", age=13.0, z=1e-4, survey="euclid", band_1="VIS", band_2="H")
        except Exception as error:  # the Dotter2016 Euclid grid is not installed
            pytest.skip(f"no Dotter2016 Euclid isochrone: {error}")
        if (isochrone.stage[isochrone.index:] < 6).all():
            pytest.skip("this Dotter2016 file has no white dwarf sequence")
        color, mag = _isochrone_locus("Dotter2016", 13.0, 1e-4, "VIS", "euclid", "H", "euclid")
        assert len(mag) == isochrone.index

        polygon = build_match_filter(
            distance_modulus=16.5, age=13.0, metallicity=1e-4, isochrone_model="Dotter2016",
            band_1="VIS", band_2="H", survey_1="euclid", survey_2="euclid",
        )
        # Main sequence: one narrow interval of colour at each magnitude
        main_sequence_colors = np.linspace(-1.0, 1.5, 2501)
        for mag_vis in np.arange(21.0, 25.01, 0.25):
            inside = Path(polygon).contains_points(
                np.c_[main_sequence_colors, np.full_like(main_sequence_colors, mag_vis)]
            )
            width = inside.sum() * (main_sequence_colors[1] - main_sequence_colors[0])
            assert 0.05 < width < 0.5, f"width {width:.2f} mag at VIS = {mag_vis}"
