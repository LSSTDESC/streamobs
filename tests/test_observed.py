"""
tests/test_observed.py
=======================
Tests for ``streamobs.observed.StreamInjector``.

All tests use the ``mock_survey`` / ``mock_injector`` fixtures from
``conftest.py``, so no real survey files are required.

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
        assert isinstance(mock_injector, StreamInjector), "Injector must be an instance of StreamInjector"
        assert hasattr(mock_injector, "survey"), "Injector must have a 'survey' property"
        assert isinstance(mock_injector.survey, Survey), "Survey property must be a Survey instance"
        assert hasattr(mock_injector, "mask_cache"), "Injector must have a 'mask_cache' property"

        # Initialize injector directly with survey name and release
        injector_direct = StreamInjector(survey="lsst", release="yr4")
        assert isinstance(injector_direct, StreamInjector), "Directly initialized injector must be an instance of StreamInjector"
        assert injector_direct.survey.name == "lsst", "Survey name must be 'lsst'"
        assert injector_direct.survey.release == "yr4", "Survey release must be 'yr4'"


# ---------------------------------------------------------------------------
# Injector behavior
# ---------------------------------------------------------------------------

@pytest.mark.observed
class TestStreamInjectorBehavior:
    """Tests for StreamInjector behavior and output structure."""
    def _verify_injected_catalog_content(self, injected_catalog, expected_columns=["ra", "dec","mag_g", "mag_r", "mag_g_obs", "mag_r_obs","flag_observed", "flag_perfect_galstarsep"]):
        """Helper method to verify the content of the injected catalog."""
        # Verify expected columns are present
        assert set(expected_columns).issubset(injected_catalog.columns), "Injected catalog must contain all expected columns"


    def test_full_injection_pipeline(self, mock_injector, stream_catalog):
        """Test the full injection pipeline with a controlled input catalog."""
        # Perform injection
        injected_catalog = mock_injector.inject(stream_catalog, perfect_galstarsep=True)

        # Minimal expected columns in the injected catalog (position, magnitude, and flags)
        self._verify_injected_catalog_content(injected_catalog,)

    def test_injection_partialinput(self, mock_injector,stream_catalog, stream_config_with_distance):
        """Test injection with a catalog that has some missing columns."""
        data_without_mag = stream_catalog.drop(columns=["mag_g", "mag_r"])
        injected_catalog = mock_injector.inject(data_without_mag, perfect_galstarsep=True, stream_config=stream_config_with_distance)
        self._verify_injected_catalog_content(injected_catalog)

    def test_random_injection(self, mock_injector, stream_catalog,seed):
        """Test random sky injection"""
        mask_type = ["footprint", "ebv", "maglim_g"]
        # Inject a first time
        stream_coord_1 = mock_injector.phi_to_radec(stream_catalog['phi1'], stream_catalog['phi2'], seed=seed,gc_frame=None,mask_type=mask_type)
        gc_1 = mock_injector._last_gc_frame

        # Inject a second time with the same random seed
        stream_coord_2 = mock_injector.phi_to_radec(stream_catalog['phi1'], stream_catalog['phi2'], seed=seed,gc_frame=None,mask_type=mask_type)
        gc_2 = mock_injector._last_gc_frame

        # Inject a 3rd time using the existing gc_frame (should not use a random
        # seed)
        stream_coord_3 = mock_injector.phi_to_radec(stream_catalog['phi1'], stream_catalog['phi2'], seed=None,gc_frame=gc_1,mask_type=mask_type)
        gc_3 = mock_injector._last_gc_frame

        def compare_coords(coord1, coord2):
            """Helper function to compare two coordinate DataFrames."""
            assert np.allclose(coord1.icrs.ra.deg, coord2.icrs.ra.deg), "RA values should be the same"
            assert np.allclose(coord1.icrs.dec.deg, coord2.icrs.dec.deg), "Dec values should be the same"
        
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
                "pole": {"ra": float(pole.ra.deg), "dec": float(pole.dec.deg), "unit": "deg"},
                "priority": str(priority),
            }
            return gc_frame_params

        def compare_gc_frames(gc1, gc2):
            """Helper function to compare two gc_frame objects."""
            assert get_gc_frame_dict(gc1) == get_gc_frame_dict(gc2), "gc_frame parameters should be the same"

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
        assert len(masks_before) > len(masks), "Mask cache should have had entries before clearing"






        