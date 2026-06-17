"""
tests/test_observed.py
=======================
Tests for ``streamobs.observed.StreamInjector``.

The injection pipeline is tested end-to-end by calling
``StreamInjector.inject()`` with controlled input DataFrames and verifying
the output columns, dtypes, and flag semantics.
"""

import numpy as np
import pandas as pd
import pytest

from streamobs.observed import StreamInjector
from streamobs.surveys import Survey

# ---------------------------------------------------------------------------
# Injector properties
# ---------------------------------------------------------------------------


@pytest.mark.observed
class TestStreamInjectorProperties:
    """Tests for StreamInjector properties and basic behavior."""

    def test_injector_initialization(self, mock_injector):
        """Test that the injector initializes with the expected properties."""
        assert isinstance(
            mock_injector, StreamInjector
        ), "Injector must be an instance of StreamInjector"
        assert hasattr(
            mock_injector, "survey"
        ), "Injector must have a 'survey' property"
        assert isinstance(
            mock_injector.survey, Survey
        ), "Survey property must be a Survey instance"
        assert hasattr(
            mock_injector, "mask_cache"
        ), "Injector must have a 'mask_cache' property"

        # Initialize injector directly with survey name and release
        injector_direct = StreamInjector(survey="lsst", release="yr4")
        assert isinstance(
            injector_direct, StreamInjector
        ), "Directly initialized injector must be an instance of StreamInjector"
        assert injector_direct.survey.name == "lsst", "Survey name must be 'lsst'"
        assert injector_direct.survey.release == "yr4", "Survey release must be 'yr4'"

    def test_injector_multisurveys(self):
        """Test that the injector can handle multiple surveys."""
        # Create a second survey and injector
        survey1_dict = {"survey": "lsst", "release": "yr4"}
        survey2_dict = {"survey": "lsst", "release": "yr5"}
        survey1 = Survey.load(**survey1_dict)
        survey2 = Survey.load(**survey2_dict)

        def test_injector_initialization_with_multiple_surveys(injector):
            """Test that the injector initializes with multiple surveys."""
            assert isinstance(
                injector, StreamInjector
            ), "Injector must be an instance of StreamInjector"
            assert hasattr(
                injector, "primary"
            ), "Injector must have a 'primary' property"
            assert isinstance(
                injector.surveys, dict
            ), "Survey property must be a dict of Survey instances"
            assert all(
                isinstance(value, Survey) for key, value in injector.surveys.items()
            ), "All elements in survey dict must be Survey instances"
            assert (
                len(injector.surveys) == 2
            ), f"Injector must have 2 surveys, got {len(injector.surveys)} with keys: {list(injector.surveys.keys())}"
            assert isinstance(
                injector.primary, Survey
            ), "Primary survey must be a Survey instance"

        # Intialize injectors for both surveys using list of survey objects
        injector1 = StreamInjector(survey=[survey1, survey2])
        test_injector_initialization_with_multiple_surveys(injector1)

        # Initialize injectors for both surveys using list of survey dicts
        injector3 = StreamInjector(survey=[survey1_dict, survey2_dict])
        test_injector_initialization_with_multiple_surveys(injector3)

        # Initialize injectors for both surveys using dicts # but I do not think
        # this is more useful than the previous approach.
        injector3 = StreamInjector({"lsst_yr4": survey1_dict, "lsst_yr5": survey2_dict})
        test_injector_initialization_with_multiple_surveys(injector3)


# ---------------------------------------------------------------------------
# Injector behavior
# ---------------------------------------------------------------------------


@pytest.mark.observed
class TestStreamInjectorBehavior:
    """Tests for StreamInjector behavior and output structure."""

    def _verify_injected_catalog_content(
        self,
        injected_catalog,
        expected_columns=[
            "ra",
            "dec",
            "lsst_yr4_g_true",
            "lsst_yr4_r_true",
            "lsst_yr4_g_obs",
            "lsst_yr4_r_obs",
            "lsst_yr4_flag_observed",
            "lsst_yr4_flag_perfect_galstarsep",
        ],
    ):
        """Helper method to verify the content of the injected catalog."""
        # Verify expected columns are present
        assert set(expected_columns).issubset(
            injected_catalog.columns
        ), "Injected catalog must contain all expected columns"

    def test_full_injection_pipeline(self, mock_injector, stream_catalog, verbose):
        """Test the full injection pipeline with a controlled input catalog."""
        # Perform injection
        injected_catalog = mock_injector.inject(
            stream_catalog, perfect_galstarsep=True, verbose=verbose
        )

        # Minimal expected columns in the injected catalog (position, magnitude, and flags)
        self._verify_injected_catalog_content(
            injected_catalog,
        )

    def test_injection_partialinput(
        self, mock_injector, stream_catalog, stream_config_with_distance, verbose
    ):
        """Test injection with a catalog that has some missing columns."""
        data_without_mag = stream_catalog.drop(
            columns=["lsst_yr4_g_true", "lsst_yr4_r_true"]
        )
        injected_catalog = mock_injector.inject(
            data_without_mag,
            perfect_galstarsep=True,
            stream_config=stream_config_with_distance,
            verbose=verbose,
        )
        self._verify_injected_catalog_content(injected_catalog)

    def test_random_injection(self, mock_injector, stream_catalog, seed, verbose):
        """Test random sky injection"""
        mask_type = ["footprint", "ebv", "maglim_g"]
        # Inject a first time
        stream_coord_1 = mock_injector.phi_to_radec(
            stream_catalog["phi1"],
            stream_catalog["phi2"],
            seed=seed,
            gc_frame=None,
            mask_type=mask_type,
            verbose=verbose,
        )
        gc_1 = mock_injector._last_gc_frame

        # Inject a second time with the same random seed
        stream_coord_2 = mock_injector.phi_to_radec(
            stream_catalog["phi1"],
            stream_catalog["phi2"],
            seed=seed,
            gc_frame=None,
            mask_type=mask_type,
            verbose=verbose,
        )
        gc_2 = mock_injector._last_gc_frame

        # Inject a 3rd time using the existing gc_frame (should not use a random
        # seed)
        stream_coord_3 = mock_injector.phi_to_radec(
            stream_catalog["phi1"],
            stream_catalog["phi2"],
            seed=None,
            gc_frame=gc_1,
            mask_type=mask_type,
            verbose=verbose,
        )
        gc_3 = mock_injector._last_gc_frame

        def compare_coords(coord1, coord2):
            """Helper function to compare two coordinate DataFrames."""
            assert np.allclose(
                coord1.icrs.ra.deg, coord2.icrs.ra.deg
            ), "RA values should be the same"
            assert np.allclose(
                coord1.icrs.dec.deg, coord2.icrs.dec.deg
            ), "Dec values should be the same"

        def get_gc_frame_dict(gc_frame):
            origin = gc_frame.origin
            pole = gc_frame.pole
            priority = gc_frame.priority
            gc_frame_params = {
                "origin": {
                    "ra": float(origin.ra.deg),
                    "dec": float(origin.dec.deg),
                    "unit": "deg",
                },
                "pole": {
                    "ra": float(pole.ra.deg),
                    "dec": float(pole.dec.deg),
                    "unit": "deg",
                },
                "priority": str(priority),
            }
            return gc_frame_params

        def compare_gc_frames(gc1, gc2):
            """Helper function to compare two gc_frame objects."""
            assert get_gc_frame_dict(gc1) == get_gc_frame_dict(
                gc2
            ), "gc_frame parameters should be the same"

        # Verify that the first two injections produce the same coordinates (same random seed)
        compare_coords(stream_coord_1, stream_coord_2)
        compare_gc_frames(gc_1, gc_2)

        # Verify that the 3rd injection produces the same coordinates (same gc_frame, no new random seed)
        compare_coords(stream_coord_1, stream_coord_3)
        compare_gc_frames(gc_1, gc_3)

        masks_before = mock_injector.list_cached_masks()
        mock_injector.clear_mask_cache()
        masks = mock_injector.list_cached_masks()
        assert len(masks) == 0, "Mask cache should be empty after clearing"
        assert len(masks_before) > len(
            masks
        ), "Mask cache should have had entries before clearing"


# ---------------------------------------------------------------------------
# complete_data + bands/survey API
# ---------------------------------------------------------------------------


@pytest.mark.observed
class TestCompleteDataAndAPI:
    """complete_data (single + multi survey) and the unified bands/survey API."""

    def test_complete_data_single_survey(self, mock_injector, stream_config_with_distance):
        """complete_data fills ra/dec + true mags and preserves existing columns."""
        df = pd.DataFrame({"phi1": [-3.0, 0.0, 3.0], "phi2": [0.0, 0.0, 0.0]})
        out = mock_injector.complete_data(
            df, bands=["g", "r"], stream_config=stream_config_with_distance, seed=1
        )
        for col in ["ra", "dec", "lsst_yr4_g_true", "lsst_yr4_r_true"]:
            assert col in out.columns, f"missing {col}"
        assert out[["lsst_yr4_g_true", "lsst_yr4_r_true"]].notna().all().all()
        # Input is not mutated (complete_data works on a copy).
        assert "ra" not in df.columns

    def test_complete_data_multisurvey(
        self, mock_multisurvey_injector, multisurvey_stream_config
    ):
        """Multi-survey complete_data fills every namespace from ONE mass draw.

        Both namespaces are backed by the same ugali survey/bands, so the shared
        masses must yield identical true magnitudes across them.
        """
        df = pd.DataFrame({"phi1": [-3.0, 0.0, 3.0], "phi2": [0.0, 0.0, 0.0]})
        out = mock_multisurvey_injector.complete_data(
            df,
            bands={"lsst_yr4": ["g", "r"], "lsst_yr5": ["g", "r"]},
            stream_config=multisurvey_stream_config,
            seed=7,
        )
        for col in [
            "ra",
            "dec",
            "lsst_yr4_g_true",
            "lsst_yr4_r_true",
            "lsst_yr5_g_true",
            "lsst_yr5_r_true",
        ]:
            assert col in out.columns, f"missing {col}"
        # Same physical star across surveys -> identical true mags.
        assert np.allclose(out["lsst_yr4_g_true"], out["lsst_yr5_g_true"])
        assert np.allclose(out["lsst_yr4_r_true"], out["lsst_yr5_r_true"])

    def test_bands_list_rejected_for_multisurvey(self, mock_multisurvey_injector):
        """A plain list of bands is ambiguous for a multi-survey injector."""
        df = pd.DataFrame({"phi1": [0.0], "phi2": [0.0]})
        with pytest.raises(ValueError):
            mock_multisurvey_injector.complete_data(df, bands=["g", "r"])

    def test_bands_dict_unknown_survey_raises(self, mock_injector):
        """A bands dict referencing an unknown namespace is rejected."""
        df = pd.DataFrame({"phi1": [0.0], "phi2": [0.0]})
        with pytest.raises(ValueError):
            mock_injector.complete_data(df, bands={"nope": ["g"]})

    def test_detect_flag_requires_survey(self, mock_injector):
        """`survey` is now a required argument of detect_flag."""
        with pytest.raises(TypeError):
            mock_injector.detect_flag(0, mag=np.array([20.0]), band="r")
