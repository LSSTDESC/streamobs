"""
tests/test_spline_model.py
==========================
Regression test for ``streamobs.model.SplineStreamModel``.

Guards against the bug where ``_create_model`` called an undefined
``_create_distance()`` (now ``_create_distance_modulus()``), which made every
spline-stream instantiation raise ``AttributeError``.
"""

import os

import pytest
import yaml

from streamobs.model import SplineStreamModel

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CFG = os.path.join(REPO, "config", "atlas_spline_config.yaml")


@pytest.mark.model
class TestSplineStreamModel:
    def test_instantiate_and_sample(self):
        """Constructing + sampling a spline stream must not raise (was AttributeError)."""
        if not os.path.exists(CFG):
            pytest.skip("atlas_spline_config.yaml not present")
        cfg = yaml.safe_load(open(CFG))["stream"]
        # The spline density/track read an external interpolation file; skip if
        # the (gitignored / downloaded) data is not available locally.
        data_file = cfg.get("filename", "").lstrip("./")
        if data_file and not os.path.exists(os.path.join(REPO, data_file)):
            pytest.skip(f"spline data file {data_file} not present")

        model = SplineStreamModel(cfg)  # previously raised AttributeError here
        df = model.sample(25)
        assert len(df) == 25
        for col in ("phi1", "phi2"):
            assert col in df.columns
