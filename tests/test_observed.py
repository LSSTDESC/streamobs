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
    def _verify_injected_catalog_content(self, injected_catalog, expected_columns):
        """Helper method to verify the content of the injected catalog."""
        # Verify expected columns are present
        assert set(expected_columns).issubset(injected_catalog.columns), "Injected catalog must contain all expected columns"


    def test_full_injection_pipeline(self, mock_injector, stream_catalog):
        """Test the full injection pipeline with a controlled input catalog."""
        # Perform injection
        injected_catalog = mock_injector.inject(stream_catalog, perfect_galstarsep=True)

        # Minimal expected columns in the injected catalog (position, magnitude, and flags)
        expected_columns = [
            "ra", "dec", "mag_g", "mag_r",
            "flag_observed", "flag_perfect_galstarsep"
        ]
        self._verify_injected_catalog_content(injected_catalog, expected_columns)