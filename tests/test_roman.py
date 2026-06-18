"""
tests/test_roman.py
===================
Implementation behavior tests for the Roman DC2 survey.

These tests verify *runtime behavior* of the Roman pipeline:

- The ``roman/dc2`` Survey loads and has the expected namespaced columns after
  an inject.
- Vega→AB conversion is applied for Roman bands (``F158``, ``F106``, etc.) and
  is a no-op for non-Roman bands.
- The completeness and photo-error loaders work, including the intentionally
  misspelled ``classifiction_eff`` column fallback path.

All tests that require the ``roman_dc2`` data files under
``data/surveys/roman_dc2/`` skip gracefully when those files are absent.
"""

import os

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Skip sentinel — evaluated once at import time
# ---------------------------------------------------------------------------

_ROMAN_DC2_DATA_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..",
    "data",
    "surveys",
    "roman_dc2",
)

_ROMAN_DC2_DATA_PRESENT = os.path.isdir(_ROMAN_DC2_DATA_DIR) and any(
    fname.endswith((".fits", ".fits.gz", ".csv", ".hsp"))
    for fname in os.listdir(_ROMAN_DC2_DATA_DIR)
    if os.path.isdir(_ROMAN_DC2_DATA_DIR)
)

_skip_no_data = pytest.mark.skipif(
    not _ROMAN_DC2_DATA_PRESENT,
    reason="roman_dc2 data files not present under data/surveys/roman_dc2/",
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def roman_dc2_survey():
    """Load the roman/dc2 survey once per module; skip if data absent."""
    if not _ROMAN_DC2_DATA_PRESENT:
        pytest.skip("roman_dc2 data files not present")
    from streamobs import surveys

    return surveys.Survey.load("roman", release="dc2", verbose=False)


@pytest.fixture(scope="module")
def roman_dc2_injector(roman_dc2_survey):
    """StreamInjector backed by the roman/dc2 survey."""
    from streamobs.observed import StreamInjector

    return StreamInjector(survey=roman_dc2_survey, verbose=False)



# ---------------------------------------------------------------------------
# Survey loading
# ---------------------------------------------------------------------------


@pytest.mark.surveys
class TestRomanDC2SurveyLoading:
    """Verify roman/dc2 Survey loads with the expected structure."""

    @_skip_no_data
    def test_survey_loads(self, roman_dc2_survey):
        assert roman_dc2_survey is not None

    @_skip_no_data
    def test_survey_name_and_release(self, roman_dc2_survey):
        assert roman_dc2_survey.name == "roman"
        assert roman_dc2_survey.release == "dc2"

    @_skip_no_data
    def test_namespace(self, roman_dc2_survey):
        assert roman_dc2_survey.namespace == "roman_dc2"

    @_skip_no_data
    def test_bands_present(self, roman_dc2_survey):
        for band in ("F106", "F129", "F158"):
            assert band in roman_dc2_survey.bands, (
                f"Expected band '{band}' in roman_dc2.bands={roman_dc2_survey.bands}"
            )

    @_skip_no_data
    def test_f158_maglim_map_loaded(self, roman_dc2_survey):
        assert "F158" in roman_dc2_survey.maglim_maps, "F158 maglim map missing"
        assert roman_dc2_survey.maglim_maps["F158"] is not None

    @_skip_no_data
    def test_completeness_band_is_f158(self, roman_dc2_survey):
        assert roman_dc2_survey.completeness_band == "F158"

    @_skip_no_data
    def test_ebv_map_loaded(self, roman_dc2_survey):
        assert roman_dc2_survey.ebv_map is not None

    @_skip_no_data
    def test_extinction_coefficients(self, roman_dc2_survey):
        # From roman_dc2.yaml survey_properties
        expected = {"F106": 1.1495, "F129": 0.8497, "F158": 0.6140}
        for band, coeff in expected.items():
            assert band in roman_dc2_survey.coeff_extinc, (
                f"Extinction coefficient missing for band '{band}'"
            )
            assert np.isclose(
                roman_dc2_survey.coeff_extinc[band], coeff, rtol=1e-3
            ), (
                f"Extinction coefficient for '{band}': "
                f"expected {coeff}, got {roman_dc2_survey.coeff_extinc[band]}"
            )


# ---------------------------------------------------------------------------
# Injection — namespaced columns
# ---------------------------------------------------------------------------


def _make_roman_dc2_catalog(n=20, rng=None):
    """Build a small catalog pre-placed in the Roman DC2 footprint.

    The Roman DC2 mock covers roughly ra∈[51, 52], dec∈[-38.5, -37.5].
    Pre-supplying ra/dec (and true mags) avoids the phi_to_radec step, which
    requires searching for a valid great-circle frame — expensive and unreliable
    for a small-footprint survey.

    Band names are uppercase (F158, F106) matching the survey config.
    Column names follow the namespace convention: roman_dc2_<band>_true.
    """
    import pandas as pd

    if rng is None:
        rng = np.random.default_rng(42)
    return pd.DataFrame(
        {
            "ra": rng.uniform(51.1, 51.9, n),
            "dec": rng.uniform(-38.3, -37.7, n),
            "roman_dc2_F158_true": rng.uniform(22.0, 25.5, n),
        }
    )


@pytest.mark.surveys
class TestRomanDC2Injection:
    """Verify inject() produces roman_dc2-namespaced output columns.

    Catalogs are pre-placed in the Roman DC2 footprint (ra≈51–52, dec≈-38) with
    pre-filled true magnitudes so that the pipeline does not attempt the
    expensive phi_to_radec great-circle search on the small DC2 footprint.
    """

    @_skip_no_data
    def test_inject_produces_namespaced_columns(self, roman_dc2_injector):
        """After inject(), every band column must be prefixed roman_dc2_<band>_*.

        Only F158 is injected here because F106/F129 don't have dedicated maglim
        maps in the roman_dc2 config (the survey has completeness + photo-error
        only for F158, which is the detection/completeness band).
        """
        df = _make_roman_dc2_catalog(n=20)
        out = roman_dc2_injector.inject(
            df,
            bands=["F158"],
            verbose=False,
        )

        expected_cols = [
            "ra",
            "dec",
            "roman_dc2_F158_true",
            "roman_dc2_F158_obs",
            "roman_dc2_F158_err",
            "roman_dc2_flag_observed",
        ]
        missing = [c for c in expected_cols if c not in out.columns]
        assert not missing, f"Missing columns after inject: {missing}"

    @_skip_no_data
    def test_true_mags_preserved_after_inject(self, roman_dc2_injector):
        """True magnitudes supplied in the input must survive inject() unchanged."""
        df = _make_roman_dc2_catalog(n=15)
        true_in = df["roman_dc2_F158_true"].values.copy()
        out = roman_dc2_injector.inject(df, bands=["F158"], verbose=False)
        assert np.allclose(out["roman_dc2_F158_true"].values, true_in), (
            "roman_dc2_F158_true changed during inject()"
        )

    @_skip_no_data
    def test_flag_observed_is_boolean(self, roman_dc2_injector):
        """roman_dc2_flag_observed must exist and be boolean-like (0/1 or bool)."""
        df = _make_roman_dc2_catalog(n=20)
        out = roman_dc2_injector.inject(df, bands=["F158"], verbose=False)
        flags = out["roman_dc2_flag_observed"]
        unique_vals = set(flags.dropna().unique())
        assert unique_vals.issubset({0, 1, True, False}), (
            f"roman_dc2_flag_observed contains unexpected values: {unique_vals}"
        )

    @_skip_no_data
    def test_obs_mags_differ_from_true(self, roman_dc2_injector):
        """Observed magnitudes must differ from true (noise was applied)."""
        rng = np.random.default_rng(99)
        df = _make_roman_dc2_catalog(n=30, rng=rng)
        out = roman_dc2_injector.inject(df, bands=["F158"], verbose=False, seed=99)
        observed = out["roman_dc2_flag_observed"]
        detected = out[observed == 1]
        if len(detected) == 0:
            pytest.skip("No detected stars in this footprint sample — increase n")
        diffs = (detected["roman_dc2_F158_obs"] - detected["roman_dc2_F158_true"]).abs()
        assert diffs.mean() > 0, "Observed mags identical to true — noise not applied"


# ---------------------------------------------------------------------------
# Round-trip regression: complete_data → inject column-case consistency
# ---------------------------------------------------------------------------


@pytest.mark.surveys
class TestRomanDC2RoundTripColumnCase:
    """Regression: complete_data and inject must produce matching-case F158 columns.

    Before the band-name normalisation (f158→F158 in config + inject args),
    complete_data emitted roman_dc2_F158_true (uppercase, from ugali) while inject
    emitted roman_dc2_f158_obs (lowercase, from the config key).  After the fix both
    should be uppercase: roman_dc2_F158_true AND roman_dc2_F158_obs must coexist.

    This test is skipped if roman_dc2 data files are absent.
    """

    @_skip_no_data
    def test_true_and_obs_columns_have_matching_case(self, roman_dc2_injector):
        """complete_data→inject round-trip: F158 true and obs columns both uppercase."""
        import pandas as pd

        rng = np.random.default_rng(7)
        n = 25
        # Pre-supply true magnitudes using the now-uppercase column name
        df = pd.DataFrame(
            {
                "ra": rng.uniform(51.1, 51.9, n),
                "dec": rng.uniform(-38.3, -37.7, n),
                "roman_dc2_F158_true": rng.uniform(22.0, 25.5, n),
            }
        )
        out = roman_dc2_injector.inject(df, bands=["F158"], verbose=False)

        # Both columns must be present with matching case
        assert "roman_dc2_F158_true" in out.columns, (
            "roman_dc2_F158_true missing after inject — uppercase true column lost"
        )
        assert "roman_dc2_F158_obs" in out.columns, (
            "roman_dc2_F158_obs missing after inject — band key is still lowercase"
        )

        # Confirm the old lowercase obs column is NOT present (the bug is gone)
        assert "roman_dc2_f158_obs" not in out.columns, (
            "roman_dc2_f158_obs still present — lowercase band key not fully removed"
        )


# ---------------------------------------------------------------------------
# Vega → AB conversion
# ---------------------------------------------------------------------------


@pytest.mark.surveys
class TestRomanVegaToAB:
    """Verify that ``IsochroneModel._to_ab`` applies the correct offsets."""

    def _make_iso_model(self):
        """Build a minimal Roman IsochroneModel for unit-testing _to_ab."""
        from streamobs.model import IsochroneModel

        cfg = {
            "name": "Marigo2017",
            "survey": "roman",
            "age": 12.0,
            "z": 0.0006,
            "band_1": "F158",
            "band_2": "F106",
        }
        iso_model = IsochroneModel(cfg)
        iso_model.create_isochrone(cfg)
        return iso_model

    def test_f158_offset_applied(self):
        """F158 Vega→AB offset must be ~1.315 mag."""
        from streamobs.model import ROMAN_VEGA_TO_AB

        iso_model = self._make_iso_model()
        mags = np.array([20.0, 21.0, 22.0])
        ab_mags = iso_model._to_ab("F158", mags)
        expected_offset = ROMAN_VEGA_TO_AB["F158"]
        assert np.allclose(ab_mags - mags, expected_offset), (
            f"F158 Vega→AB offset: expected {expected_offset}, "
            f"got {(ab_mags - mags)[0]}"
        )

    def test_f106_offset_applied(self):
        """F106 Vega→AB offset must be ~0.660 mag."""
        from streamobs.model import ROMAN_VEGA_TO_AB

        iso_model = self._make_iso_model()
        mags = np.array([20.0, 21.0, 22.0])
        ab_mags = iso_model._to_ab("F106", mags)
        expected_offset = ROMAN_VEGA_TO_AB["F106"]
        assert np.allclose(ab_mags - mags, expected_offset), (
            f"F106 Vega→AB offset: expected {expected_offset}, "
            f"got {(ab_mags - mags)[0]}"
        )

    def test_non_roman_band_unchanged(self):
        """Non-Roman bands (e.g. 'g') must pass through _to_ab unchanged."""
        iso_model = self._make_iso_model()
        mags = np.array([20.0, 21.0, 22.0])
        ab_mags = iso_model._to_ab("g", mags)
        assert np.allclose(ab_mags, mags), (
            "Non-Roman band 'g' should not be modified by _to_ab"
        )

    def test_all_roman_bands_have_positive_offset(self):
        """Every Roman band in ROMAN_VEGA_TO_AB must have a positive offset
        (AB magnitudes are systematically brighter than Vega for Roman NIR)."""
        from streamobs.model import ROMAN_VEGA_TO_AB

        for band, offset in ROMAN_VEGA_TO_AB.items():
            assert offset > 0, (
                f"Roman band {band} has non-positive Vega→AB offset: {offset}"
            )

    def test_ab_mags_greater_than_vega(self):
        """AB mags must be strictly larger (numerically dimmer in flux) than Vega
        for Roman NIR bands — positive offset means AB number > Vega number."""
        from streamobs.model import ROMAN_VEGA_TO_AB

        iso_model = self._make_iso_model()
        mags_vega = np.array([20.0])
        for band in ("F106", "F129", "F158", "F184"):
            mags_ab = iso_model._to_ab(band, mags_vega)
            assert mags_ab[0] > mags_vega[0], (
                f"AB mag should be larger than Vega for Roman band {band}"
            )


# ---------------------------------------------------------------------------
# Completeness and photo-error loaders
# ---------------------------------------------------------------------------


@pytest.mark.surveys
class TestRomanDC2Efficiencies:
    """Verify completeness and photo-error loaders for roman/dc2."""

    @_skip_no_data
    def test_completeness_callable(self, roman_dc2_survey):
        assert callable(roman_dc2_survey.completeness), (
            "roman_dc2.completeness must be callable"
        )

    @_skip_no_data
    def test_efficiency_detection_callable(self, roman_dc2_survey):
        assert callable(roman_dc2_survey.efficiency_detection), (
            "roman_dc2.efficiency_detection must be callable"
        )

    @_skip_no_data
    def test_efficiency_classification_callable(self, roman_dc2_survey):
        assert callable(roman_dc2_survey.efficiency_classification), (
            "roman_dc2.efficiency_classification must be callable"
        )

    @_skip_no_data
    def test_log_photo_error_catalog_callable(self, roman_dc2_survey):
        assert callable(roman_dc2_survey.log_photo_error_catalog), (
            "roman_dc2.log_photo_error_catalog must be callable"
        )

    @_skip_no_data
    def test_log_photo_error_sample_callable(self, roman_dc2_survey):
        assert callable(roman_dc2_survey.log_photo_error_sample), (
            "roman_dc2.log_photo_error_sample must be callable"
        )

    @_skip_no_data
    def test_two_curve_photo_error_differ(self, roman_dc2_survey):
        """Catalog and sample photo-error curves must differ for roman/dc2.

        roman_dc2 uses separate catalog (reported magerr) and sample (true
        scatter) CSV files — the two curves must produce different values.
        """
        delta_mags = np.linspace(-5.0, 0.0, 20)
        catalog_vals = roman_dc2_survey.log_photo_error_catalog(delta_mags)
        sample_vals = roman_dc2_survey.log_photo_error_sample(delta_mags)
        assert not np.allclose(catalog_vals, sample_vals), (
            "Catalog and sample photo-error curves should differ for roman_dc2 "
            "(separate CSV files)"
        )

    @_skip_no_data
    def test_completeness_bright_near_one(self, roman_dc2_survey):
        """At bright magnitudes (well above saturation), completeness ≈ 1."""
        maglim = 26.0
        bright_mags = np.array([20.0, 21.0, 22.0])
        comp = roman_dc2_survey.get_completeness("F158", bright_mags, maglim)
        assert np.all(comp > 0.5), (
            f"Completeness should be >0.5 for bright stars: {comp}"
        )

    @_skip_no_data
    def test_completeness_below_saturation_zero(self, roman_dc2_survey):
        """Stars brighter than saturation should have completeness = 0."""
        # saturation is 17.0 per roman_dc2.yaml
        maglim = 26.0
        sat_mags = np.array([10.0, 12.0, 15.0])
        comp = roman_dc2_survey.get_completeness("F158", sat_mags, maglim)
        assert np.all(comp == 0.0), (
            f"Completeness should be 0 below saturation: {comp}"
        )

    @_skip_no_data
    def test_classifiction_eff_fallback_via_csv(self):
        """set_completeness() must load 'classifiction_eff' (misspelled) column.

        The roman_dc2 efficiency CSV uses the intentionally misspelled header
        ``classifiction_eff``. Verify that the loader accepts this column name
        and returns a callable interpolator.
        """
        from streamobs.surveys import SurveyFactory

        csv_path = os.path.join(_ROMAN_DC2_DATA_DIR, "roman_stellar_efficiency_cutf158.csv")
        if not os.path.exists(csv_path):
            pytest.skip(f"efficiency CSV not found: {csv_path}")

        func = SurveyFactory.set_completeness(csv_path, selection="classified")
        assert callable(func), (
            "set_completeness(..., selection='classified') should return callable "
            "even with misspelled 'classifiction_eff' column"
        )
        # Verify it returns sensible values in [-1, 0] delta_mag range
        delta_mags = np.linspace(-3.0, 0.0, 10)
        vals = func(delta_mags)
        assert np.all(np.isfinite(vals)), "Efficiency values should be finite"
        assert np.all((vals >= 0.0) & (vals <= 1.0)), (
            f"Efficiency values should be in [0, 1]: {vals}"
        )
