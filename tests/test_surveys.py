"""
tests/test_surveys.py
=====================
Tests for ``streamobs.surveys.Survey``.

To add a new survey, just append an entry to ``SURVEY_REGISTRY`` at the top
of this file. Every test class below will automatically run against it.

Registry format
---------------
Each entry is a dict with:
    survey  : str           — name passed to Survey.load()
    release : str or None   — release passed to Survey.load()
    bands   : list[str]     — bands expected to be present
    expected_maglim : list[str] — bands for which maglim maps should be present

Example
-------
    {"survey": "lsst", "release": "yr5", "bands": ["g", "r"], "expected_maglim": ["g", "r"]},
"""

import pytest
from streamobs import surveys
import numpy as np


# ---------------------------------------------------------------------------
# Registry — add new surveys here
# ---------------------------------------------------------------------------

SURVEY_REGISTRY = [
    {"survey": "lsst", "release": "yr4",  "expected_bands": ["g", "r"], "expected_maglim": ['g', 'r']},
    {"survey": "des", "release": "yr6", "expected_bands": ["g", "r"], "expected_maglim": ['g', 'r']},
]

# IDs shown in pytest output, e.g. "lsst_yr4"
def _survey_id(entry):
    return f"{entry['survey']}_{entry['release'] or 'base'}"

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.mark.surveys
@pytest.fixture(scope="module", params=SURVEY_REGISTRY, ids=[_survey_id(e) for e in SURVEY_REGISTRY])
def loaded_survey(request):
    """Load each registered survey once per module and cache it."""
    entry = request.param
    survey = surveys.Survey.load(
        survey=entry["survey"],
        release=entry["release"],
        verbose=False,
    )
    # Attach the registry entry so tests can access expected values
    survey._test_entry = entry
    return survey


# ---------------------------------------------------------------------------
# Class 1 — Can the survey be created?
# ---------------------------------------------------------------------------
@pytest.mark.surveys
class TestSurveyCreation:
    """Verify that every registered survey can be loaded without error."""

    def test_survey_loads(self, loaded_survey):
        assert loaded_survey is not None

    def test_returns_survey_instance(self, loaded_survey):
        assert isinstance(loaded_survey, surveys.Survey)

    def test_name_matches_registry(self, loaded_survey):
        assert loaded_survey.name == loaded_survey._test_entry["survey"]

    def test_release_matches_registry(self, loaded_survey):
        assert loaded_survey.release == loaded_survey._test_entry["release"]


# ---------------------------------------------------------------------------
# Class 2 — Does the loaded survey have the expected properties?
# ---------------------------------------------------------------------------
@pytest.mark.surveys
class TestSurveyProperties:
    """
    Check that each loaded survey satisfies the properties we rely on
    throughout the pipeline.
    """

    def test_expected_bands_present(self, loaded_survey):
        expected = set(loaded_survey._test_entry["expected_bands"])
        assert expected.issubset(set(loaded_survey.bands)), (
            f"Missing bands: {expected - set(loaded_survey.bands)}"
        )

    def test_maglim_maps_loaded_for_each_band(self, loaded_survey):
        expected_maglim = set(loaded_survey._test_entry["expected_maglim"])
        for band in loaded_survey.bands:
            if band in expected_maglim:
                assert band in loaded_survey.maglim_maps, f"No maglim map for band '{band}'"
                assert loaded_survey.maglim_maps[band] is not None, (
                    f"maglim_maps['{band}'] is None"
                )

    def test_ebv_map_loaded(self, loaded_survey):
        assert loaded_survey.coeff_extinc is not None, "Extinction coefficients are None"
        assert loaded_survey.ebv_map is not None, "EBV map is None"

    def test_coverage_map_loaded(self, loaded_survey):
        assert loaded_survey.coverage is not None, "Coverage map is None"

    def test_errors_loaded(self, loaded_survey):
        assert loaded_survey.sys_error is not None, "Systematic error dict is None"
        assert loaded_survey.log_photo_error is not None, "Log photo error function is None"

    def test_delta_saturation_loaded(self, loaded_survey):
        assert loaded_survey.delta_saturation is not None, "delta_saturation is None"
        assert loaded_survey.saturation is not None, "Saturation dict is None"


    def test_efficiencies_loaded(self, loaded_survey):
        assert loaded_survey.completeness is not None, "Completeness function is None"
        assert loaded_survey.completeness_band is not None, "Completeness band is None"
        assert loaded_survey.completeness_band in loaded_survey.bands, (
            f"Completeness band '{loaded_survey.completeness_band}' not in survey bands {loaded_survey.bands}"
        )
        assert hasattr(loaded_survey.completeness, "__call__"), "Completeness is not callable"
        assert hasattr(loaded_survey.efficiency_classification, "__call__"), "Efficiency classification is not callable"
        assert hasattr(loaded_survey.efficiency_detection, "__call__"), "Efficiency detection is not callable"


    def test_extinction_behavior(self, loaded_survey):

        # Points within galactic plane
        la_plane = np.linspace(0, 360, 10)
        lb_plane = np.zeros(len(la_plane))

        # Clearly outside the Galactic plane
        la_high = np.linspace(0, 360, 10)
        lb_high = np.linspace(80, 90, len(la_high))

        # Convert those coordinates to healpix pixels
        from astropy.coordinates import SkyCoord
        from astropy import units as u
        import healpy as hp

        nside_ebv = hp.npix2nside(len(loaded_survey.ebv_map))
        # Convert coordinates to healpix pixels

        coord_plane = SkyCoord(
            l=la_plane * u.deg,
            b=lb_plane * u.deg,
            frame="galactic",
        )

        coord_high = SkyCoord(
            l=la_high * u.deg,
            b=lb_high * u.deg,
            frame="galactic",
        )

        pix_plane = hp.ang2pix(nside_ebv, coord_plane.l.degree, coord_plane.b.degree, lonlat=True)
        pix_high = hp.ang2pix(nside_ebv, coord_high.l.degree, coord_high.b.degree, lonlat=True)

        for band in loaded_survey.bands:
            ext_plane = loaded_survey.get_extinction(band, pix_plane)
            ext_high = loaded_survey.get_extinction(band, pix_high)
            mean_plane = np.mean(ext_plane)
            mean_high = np.mean(ext_high)
            assert mean_plane > mean_high, (
                f"Extinction should be higher in the galactic plane than at high latitude for band '{band}'"
            )

    def test_completeness_behavior(self, loaded_survey, saturation_magnitudes, bright_magnitudes, faint_magnitudes, base_maglim):
        assert loaded_survey.completeness is not None, "Completeness function is None"

        # Test that completeness is ~1 for bright stars and ~0 for faint stars
        completeness_band = loaded_survey.completeness_band
        assert completeness_band is not None, "Completeness band is None"
        assert completeness_band in loaded_survey.bands, (
            f"Completeness band '{completeness_band}' not in survey bands {loaded_survey.bands}"
        )
        assert completeness_band in loaded_survey.maglim_maps, (
            f"Completeness band '{completeness_band}' does not have a maglim map"
        )
        
        sat = loaded_survey.saturation[completeness_band]

        # work only with magnitudes above or below saturation
        sat_mag = saturation_magnitudes[saturation_magnitudes < sat]
        bright_mag = bright_magnitudes[bright_magnitudes > sat]
        faint_mag = faint_magnitudes[faint_magnitudes > sat]

        # Verify completeness behavior in each regime
        if len(sat_mag) > 0:
            comp_sat = loaded_survey.get_completeness(completeness_band, sat_mag, base_maglim)
            assert np.all(comp_sat == 0.0), (
                f"Completeness should be 0 for magnitudes below saturation in band '{completeness_band}'"
            )
        if len(bright_mag) > 0:
            comp_bright = loaded_survey.get_completeness(completeness_band, bright_mag, base_maglim)
            assert np.all(comp_bright > 0.9), (
                f"Completeness should be near 1 for magnitudes well above saturation in band '{completeness_band}'"
            )
        if len(faint_mag) > 0:
            comp_faint = loaded_survey.get_completeness(completeness_band, faint_mag, base_maglim)
            assert np.all(comp_faint < 0.1), (
                f"Completeness should be near 0 for magnitudes well below saturation in band '{completeness_band}'"
            )

    def test_log_photo_error_behavior(self, loaded_survey, saturation_magnitudes, bright_magnitudes, faint_magnitudes, base_maglim):
        assert loaded_survey.log_photo_error is not None, "Log photo error function is None"

        # Test that log photo error behaves reasonably across magnitude ranges
        for band in loaded_survey.bands:
            sat = loaded_survey.saturation[band]
            sys_error = loaded_survey.sys_error[band]
            assert sys_error > 0, f"Systematic error should be positive for band '{band}'"

            # work only with magnitudes above or below saturation
            sat_mag = saturation_magnitudes[saturation_magnitudes < sat]
            bright_mag = bright_magnitudes[bright_magnitudes > sat]
            faint_mag = faint_magnitudes[faint_magnitudes > sat]

            bright_mag_mean, faint_mag_mean  = None, None
            if len(sat_mag) > 0:
                err_sat = loaded_survey.get_photo_error(band, sat_mag, base_maglim)
                assert np.all(err_sat > 5.0), (
                    f"Photo errors should be large for magnitudes below saturation in band '{band}'"
                )
            if len(bright_mag) > 0:
                err_bright = loaded_survey.get_photo_error(band, bright_mag, base_maglim)
                assert np.all(err_bright < 3*sys_error), (
                    f"Photo errors should be close to systematic error for bright magnitudes in band '{band}'"
                )
                bright_mag_mean = np.mean(err_bright)
            if len(faint_mag) > 0:
                err_faint = loaded_survey.get_photo_error(band, faint_mag, base_maglim)
                assert np.all(err_faint > 20*sys_error), (
                    f"Photo errors should be large for faint magnitudes in band '{band}'"
                )
                faint_mag_mean = np.mean(err_faint)
            
            if faint_mag_mean is not None and bright_mag_mean is not None:
                assert faint_mag_mean > bright_mag_mean, (
                    f"Mean photo error should increase with magnitude in band '{band}'"
                )

            error_at_maglim = loaded_survey.get_photo_error(band, base_maglim, base_maglim)
            snr_at_maglim = 1 / error_at_maglim
            assert np.isclose(snr_at_maglim, 5.0, atol = 0.25), (
                f"Photo error at maglim should correspond to SNR=5 for band '{band}'"
            )

