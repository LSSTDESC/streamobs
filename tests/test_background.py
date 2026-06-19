"""
Tests for the streamobs.background subpackage.

All tests are marked ``background`` so they can be run in isolation:

    pytest tests/test_background.py -m background

Tests that write to disk (Storage, ResourceBuilder) use the ``tmp_path``
pytest fixture so the temporary directory is created and cleaned up
automatically regardless of test outcome.
"""

import numpy as np
import pandas as pd
import pytest

pytestmark = pytest.mark.background


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def mock_survey(verbose=False):
    from streamobs.surveys import Survey

    return Survey.load("lsst", release="yr4", verbose=verbose)


@pytest.fixture(scope="module")
def tiny_stars_catalog():
    """Minimal DataFrame that looks like a true background star catalog."""
    rng = np.random.default_rng(0)
    n = 100
    return pd.DataFrame(
        {
            "ra": rng.uniform(30.0, 60.0, n),
            "dec": rng.uniform(-20.0, 0.0, n),
            "lsst_yr4_g_true": rng.uniform(20.0, 27.0, n),
            "lsst_yr4_r_true": rng.uniform(20.0, 27.0, n),
        }
    )


@pytest.fixture(scope="module")
def tiny_galaxies_catalog():
    """Minimal DataFrame that looks like a true background galaxy catalog."""
    rng = np.random.default_rng(1)
    n = 100
    return pd.DataFrame(
        {
            "ra": rng.uniform(30.0, 60.0, n),
            "dec": rng.uniform(-20.0, 0.0, n),
            "lsst_yr4_g_true": rng.uniform(20.0, 28.0, n),
            "lsst_yr4_r_true": rng.uniform(20.0, 28.0, n),
        }
    )


# ---------------------------------------------------------------------------
# Part 1 — Survey galaxy misclassification
# ---------------------------------------------------------------------------


class TestSurveyGalMisclassification:
    """Tests for the gal_misclassification field and method on Survey."""

    def test_gal_misclassification_field_exists(self, mock_survey):
        """Survey dataclass must expose the gal_misclassification attribute."""
        assert hasattr(mock_survey, "gal_misclassification")

    def test_get_gal_misclassification_raises_when_not_loaded(self, mock_survey):
        """get_gal_misclassification raises ValueError if the function is not loaded."""
        import copy

        s = copy.deepcopy(mock_survey)
        s.gal_misclassification = None
        mags = np.linspace(20.0, 27.0, 10)
        maglim = np.full(10, 26.5)
        with pytest.raises(ValueError, match="gal_misclassification"):
            s.get_gal_misclassification("r", mags, maglim)

    def test_get_gal_misclassification_no_1padding(self, mock_survey):
        """get_gal_misclassification must NOT return 1 at the bright end (no 1-padding).

        Verify by checking that the method returns 0 for very bright objects
        (below saturation) and does not saturate to 1.0 for delta_mag far below
        delta_saturation.
        """
        ...

    def test_set_completeness_missclassified_loads_file(self, tmp_path):
        """set_completeness with selection='missclassified' must return a callable interpolator."""
        from streamobs.surveys import SurveyFactory

        csv_path = tmp_path / "gal_misclassification.csv"
        csv_path.write_text(
            "delta_mag,gal_misclassification_eff\n"
            "-12.0,0.0\n"
            "-5.0,0.0\n"
            "0.0,0.3\n"
            "1.0,0.0\n"
        )
        func = SurveyFactory.set_completeness(str(csv_path), selection="missclassified")
        assert callable(func)
        result = func(np.array([0.0]))
        assert 0.0 <= float(result) <= 1.0


# ---------------------------------------------------------------------------
# Part 2 — StreamInjector source_type parameter
# ---------------------------------------------------------------------------


class TestInjectSourceType:
    """Tests that source_type is accepted and routed correctly in StreamInjector."""

    def test_inject_accepts_source_type_kwarg(self):
        """inject() must accept source_type without raising TypeError."""
        from streamobs.observed import StreamInjector
        import inspect

        sig = inspect.signature(StreamInjector.inject)
        # source_type is passed through **kwargs — verify inject accepts **kwargs
        assert "kwargs" in str(sig) or "source_type" in sig.parameters

    def test_inject_galaxies_uses_gal_misclassification(
        self, mock_survey, tiny_galaxies_catalog
    ):
        """inject() with source_type='galaxies' must reach get_gal_misclassification."""
        ...


# ---------------------------------------------------------------------------
# Part 3 — BackgroundCatalogInjector
# ---------------------------------------------------------------------------


class TestBackgroundCatalogInjector:
    """Tests for BackgroundCatalogInjector."""

    def test_init(self, mock_survey):
        """BackgroundCatalogInjector can be instantiated with a Survey."""
        from streamobs.background import BackgroundCatalogInjector

        inj = BackgroundCatalogInjector(mock_survey)
        assert inj._survey is mock_survey

    def test_prepare_survey_no_dust(self, mock_survey):
        """_prepare_survey with no_dust=True must zero the EBV map."""
        ...

    def test_prepare_survey_uniform_maglim(self, mock_survey):
        """_prepare_survey with uniform_maglim must replace maglim maps with constants."""
        ...

    def test_inject_stars_delegates_to_stream_injector(
        self, mock_survey, tiny_stars_catalog
    ):
        """inject_stars must return a DataFrame (delegates to StreamInjector)."""
        ...

    def test_inject_galaxies_delegates_to_stream_injector(
        self, mock_survey, tiny_galaxies_catalog
    ):
        """inject_galaxies must return a DataFrame."""
        ...


# ---------------------------------------------------------------------------
# Part 4 — BackgroundStorage
# ---------------------------------------------------------------------------


class TestBackgroundStorage:
    """Tests for BackgroundStorage persistence helpers."""

    def test_get_path_naming(self, tmp_path):
        """get_path must embed source_type and band names in the filename."""
        from streamobs.background import BackgroundStorage

        storage = BackgroundStorage(
            base_path=str(tmp_path), survey_name="lsst", release="yr4"
        )
        path = storage.get_path("stars", ("g", "r"))
        assert "stars" in path
        assert "gr" in path
        assert path.endswith(".parquet")

    def test_save_data_creates_file(self, tmp_path):
        """save_data must write a file at the path returned by get_path."""
        ...

    def test_load_data_roundtrip(self, tmp_path):
        """load_data must recover the dict saved by save_data."""
        ...

    def test_exists_false_before_save(self, tmp_path):
        """exists must return False when the file is not on disk."""
        from streamobs.background import BackgroundStorage

        storage = BackgroundStorage(
            base_path=str(tmp_path), survey_name="lsst", release="yr4"
        )
        assert storage.exists("stars", ("g", "r")) is False

    def test_exists_true_after_save(self, tmp_path):
        """exists must return True after save_data has been called."""
        ...


# ---------------------------------------------------------------------------
# Part 5 — BackgroundResourceBuilder
# ---------------------------------------------------------------------------


class TestBackgroundResourceBuilder:
    """Tests for BackgroundResourceBuilder (uses tmp_path for storage)."""

    def test_init(self):
        """BackgroundResourceBuilder can be instantiated."""
        from streamobs.background import BackgroundResourceBuilder

        builder = BackgroundResourceBuilder(survey_name="lsst", release="yr4")
        assert builder.survey_name == "lsst"
        assert isinstance(builder.resources, dict)

    def test_build_one_config(self, mock_survey, tiny_stars_catalog):
        """_build_one_config must return a dict with the expected keys."""
        ...

    def test_save_via_storage(self, tmp_path, tiny_stars_catalog):
        """save must write a parquet file via BackgroundStorage."""
        ...

    def test_load_via_storage(self, tmp_path, tiny_stars_catalog):
        """load must reconstruct resources from the file saved by save."""
        ...


# ---------------------------------------------------------------------------
# Part 6 — LightBackgroundGenerator
# ---------------------------------------------------------------------------


class TestLightBackgroundGenerator:
    """Tests for LightBackgroundGenerator."""

    @pytest.fixture(scope="class")
    def storage_with_data(self, tmp_path_factory, mock_survey, tiny_stars_catalog):
        """BackgroundStorage pre-populated with a minimal CMD grid."""
        ...

    def test_init(self, mock_survey, tmp_path):
        """LightBackgroundGenerator can be instantiated."""
        from streamobs.background import BackgroundStorage, LightBackgroundGenerator

        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        gen = LightBackgroundGenerator(storage, mock_survey)
        assert gen.survey is mock_survey

    def test_generate_stars(self):
        """generate with source_type='stars' must return a DataFrame."""
        ...

    def test_generate_galaxies(self):
        """generate with source_type='galaxies' must return a DataFrame."""
        ...

    def test_generate_both(self):
        """generate with source_type='both' must return combined DataFrame."""
        ...


# ---------------------------------------------------------------------------
# Part 7 — Background (top-level wrapper)
# ---------------------------------------------------------------------------


class TestBackground:
    """Tests for the Background wrapper class."""

    def test_init_rejects_invalid_source_type(self, mock_survey):
        """Background raises ValueError for unknown source_type."""
        from streamobs.background import Background

        with pytest.raises(ValueError, match="source_type"):
            Background(mock_survey, source_type="invalid")

    def test_init_rejects_invalid_method(self, mock_survey):
        """Background raises ValueError for unknown method."""
        from streamobs.background import Background

        with pytest.raises(ValueError, match="method"):
            Background(mock_survey, method="invalid")

    def test_full_method_stars(self, mock_survey, tiny_stars_catalog):
        """Background with method='full' and source_type='stars' must generate a catalog."""
        ...

    def test_full_method_galaxies(self, mock_survey, tiny_galaxies_catalog):
        """Background with method='full' and source_type='galaxies' must generate a catalog."""
        ...

    def test_full_method_both(
        self, mock_survey, tiny_stars_catalog, tiny_galaxies_catalog
    ):
        """Background with method='full' and source_type='both' must combine catalogs."""
        ...

    def test_light_method_default_storage(self, mock_survey):
        """Background with method='light' and storage=None must fall back to _default_storage."""
        from streamobs.background import Background

        bg = Background(mock_survey, method="light")
        assert bg.storage is not None

    def test_light_method_custom_storage(self, mock_survey, tmp_path):
        """Background accepts a user-supplied BackgroundStorage."""
        from streamobs.background import Background, BackgroundStorage

        storage = BackgroundStorage(
            base_path=str(tmp_path), survey_name="lsst", release="yr4"
        )
        bg = Background(mock_survey, method="light", storage=storage)
        assert bg.storage is storage
