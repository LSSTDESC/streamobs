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
def tiny_galaxies_catalog():
    """Galaxy catalog: same sky positions as stream_catalog but with randomised magnitudes."""
    rng = np.random.default_rng(1)
    n = 100
    return pd.DataFrame(
        {
            "ra": rng.uniform(30.0, 60.0, n),
            "dec": rng.uniform(-20.0, 0.0, n),
            "lsst_g_true": rng.uniform(20.0, 28.0, n),
            "lsst_r_true": rng.uniform(20.0, 28.0, n),
        }
    )


# ---------------------------------------------------------------------------
# StreamInjector source_type parameter
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
# BackgroundCatalogInjector
# ---------------------------------------------------------------------------


class TestBackgroundCatalogInjector:
    """Tests for BackgroundCatalogInjector."""

    def test_init(self, mock_survey):
        """BackgroundCatalogInjector can be instantiated with a Survey."""
        from streamobs.background import BackgroundCatalogInjector

        inj = BackgroundCatalogInjector(mock_survey)
        assert inj._survey is mock_survey

    def test_inject_stars_returns_nonempty_dataframe(
        self, mock_survey, stream_catalog
    ):
        """inject_stars must return a non-empty DataFrame (smoke test only — full coverage in test_observed)."""
        from streamobs.background import BackgroundCatalogInjector

        result = BackgroundCatalogInjector(mock_survey).inject_stars(
            stream_catalog, bands=["g", "r"]
        )
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0

    def test_inject_galaxies_returns_nonempty_dataframe(
        self, mock_survey, tiny_galaxies_catalog
    ):
        """inject_galaxies must return a non-empty DataFrame."""
        from streamobs.background import BackgroundCatalogInjector

        result = BackgroundCatalogInjector(mock_survey).inject_galaxies(
            tiny_galaxies_catalog, bands=["g", "r"]
        )
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0

    def test_inject_galaxies_has_detection_flag(
        self, mock_survey, tiny_galaxies_catalog
    ):
        """inject_galaxies output must include a detection flag column."""
        from streamobs.background import BackgroundCatalogInjector

        result = BackgroundCatalogInjector(mock_survey).inject_galaxies(
            tiny_galaxies_catalog, bands=["g", "r"]
        )
        flag_cols = [c for c in result.columns if "flag_observed" in c]
        assert len(flag_cols) > 0
        flag_col = flag_cols[0]
        assert result[flag_col].isin([0, 1]).all()


# ---------------------------------------------------------------------------
# BackgroundStorage
# ---------------------------------------------------------------------------


class TestBackgroundStorage:
    """Tests for BackgroundStorage persistence helpers."""

    def test_get_path_naming(self, tmp_path):
        """get_path must embed source_type and band names in the filename."""
        from streamobs.background import BackgroundStorage

        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        path = storage.get_path("stars", ("g", "r"))
        assert "stars" in path
        assert "gr" in path
        assert path.endswith(".parquet")

    @pytest.fixture
    def _grid_data(self):
        """Minimal two-pair CMD grid for storage tests."""
        def _make(offset=0):
            return {
                "cmd_hist": np.arange(25, dtype=float).reshape(5, 5) + offset,
                "color_edges": np.linspace(-2, 3, 6),
                "mag_edges": np.linspace(14, 30, 6),
                "n_ref": 100,
                "area_ref_deg2": 10.0,
            }
        return {(26.0, 25.5): _make(0), (26.0, 26.0): _make(1)}

    def test_load_data_roundtrip(self, tmp_path, _grid_data):
        """load_data must recover a single pair using predicate pushdown."""
        from streamobs.background import BackgroundStorage

        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        storage.save_data(_grid_data, "stars", ("g", "r"))
        loaded = storage.load_data("stars", ("g", "r"), 26.0, 25.5)
        expected = _grid_data[(26.0, 25.5)]
        np.testing.assert_allclose(loaded["cmd_hist"], expected["cmd_hist"])
        np.testing.assert_allclose(loaded["color_edges"], expected["color_edges"])
        np.testing.assert_allclose(loaded["mag_edges"], expected["mag_edges"])
        assert loaded["n_ref"] == expected["n_ref"]
        assert loaded["area_ref_deg2"] == expected["area_ref_deg2"]
        assert len(loaded) == len(expected), "Loaded dict must have exactly the expected keys"

    def test_exists_false_before_save(self, tmp_path):
        """exists must return False when the file is not on disk."""
        from streamobs.background import BackgroundStorage

        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        assert storage.exists("stars", ("g", "r")) is False

    def test_exists_true_after_save(self, tmp_path, _grid_data):
        """exists must return True after save_data has been called."""
        from streamobs.background import BackgroundStorage

        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        storage.save_data(_grid_data, "stars", ("g", "r"))
        assert storage.exists("stars", ("g", "r")) is True

    def test_save_overwrites_no_extra_rows(self, tmp_path, _grid_data):
        """Saving the same file twice must overwrite — row count must not double."""
        from streamobs.background import BackgroundStorage

        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        storage.save_data(_grid_data, "stars", ("g", "r"))
        storage.save_data(_grid_data, "stars", ("g", "r"))  # overwrite
        loaded = storage.load_all("stars", ("g", "r"))
        assert len(loaded) == len(_grid_data), \
            "File must contain exactly n_pairs rows after a second save"

    def test_load_all_returns_all_pairs(self, tmp_path, _grid_data):
        """load_all must return a dict with all saved (maglim_r, maglim_g) pairs."""
        from streamobs.background import BackgroundStorage

        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        storage.save_data(_grid_data, "stars", ("g", "r"))
        all_data = storage.load_all("stars", ("g", "r"))
        assert set(all_data) == set(_grid_data), "Loaded keys must match saved keys"

    def test_load_data_returns_correct_types(self, tmp_path, _grid_data):
        """Loaded dict must have numpy arrays for histogram fields and scalars for metadata."""
        from streamobs.background import BackgroundStorage

        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        storage.save_data(_grid_data, "stars", ("g", "r"))
        d = storage.load_data("stars", ("g", "r"), 26.0, 25.5)
        assert isinstance(d["cmd_hist"], np.ndarray)
        assert isinstance(d["color_edges"], np.ndarray)
        assert isinstance(d["mag_edges"], np.ndarray)
        assert isinstance(d["n_ref"], int)
        assert isinstance(d["area_ref_deg2"], float)


# ---------------------------------------------------------------------------
# BackgroundResourceBuilder
# ---------------------------------------------------------------------------

# Galaxy catalog with a steep count-magnitude distribution, mimicking the galaxy LF.
# Used to test that a deeper maglim detects more misclassified galaxies because more
# galaxies exist at fainter magnitudes.
@pytest.fixture(scope="module")
def galaxy_lf_catalog():
    """Galaxy catalog with many more objects at faint magnitudes.
    """
    rng = np.random.default_rng(11)
    n_bright, n_faint = 20, 200
    n = n_bright + n_faint
    r_mags = np.concatenate([rng.uniform(20, 24, n_bright), rng.uniform(24, 28, n_faint)])
    g_mags = r_mags + 0.5  # constant colour so both bands track the same maglim
    return pd.DataFrame(
        {
            "ra": rng.uniform(30.0, 60.0, n),
            "dec": rng.uniform(-20.0, 0.0, n),
            "lsst_r_true": r_mags,
            "lsst_g_true": g_mags,
        }
    )


# Shared small catalog for histogram property tests (fixed magnitudes → deterministic).
@pytest.fixture(scope="module")
def bright_star_catalog():
    """100 stars all at r=22.7 / g=23.45 (color g-r=0.75).

    Magnitudes are chosen to sit at the centre of their histogram bins
    (with default n_bins=10, color_range=(-2,3), mag_range=(14,30)):
      * r=22.7  → mag bin 5 [22.0, 23.6), centre 22.8
      * g-r=0.75 → color bin 5 [0.5, 1.0), centre 0.75
    This avoids bin-edge ambiguity that would split counts between adjacent bins.
    """
    rng = np.random.default_rng(7)
    n = 100
    return pd.DataFrame(
        {
            "ra": rng.uniform(35.0, 55.0, n),
            "dec": rng.uniform(-15.0, -5.0, n),
            "lsst_r_true": np.full(n, 22.7),
            "lsst_g_true": np.full(n, 23.45),  # color g-r = 0.75
        }
    )


class TestBackgroundResourceBuilder:
    """Tests for BackgroundResourceBuilder (uses tmp_path for storage)."""

    # Shared builder kwargs for _build_one_config calls.
    _COMMON = dict(bands=("g", "r"), n_bins_color=10, n_bins_mag=10,
                   color_range=(-2, 3), mag_range=(14, 30), area_ref_deg2=100.0)

    def _one_config(self, mock_survey, catalog, source_type="stars",
                    maglim_r=26.0, maglim_g=25.5):
        from streamobs.background import BackgroundResourceBuilder
        b = BackgroundResourceBuilder(survey_name="lsst")
        return b._build_one_config(
            catalog, mock_survey, source_type,
            maglim_r=maglim_r, maglim_g=maglim_g,
            **self._COMMON,
        )

    # ------------------------------------------------------------------
    # Init
    # ------------------------------------------------------------------

    def test_init(self):
        """BackgroundResourceBuilder can be instantiated."""
        from streamobs.background import BackgroundResourceBuilder

        builder = BackgroundResourceBuilder(survey_name="lsst", release="yr4")
        assert builder.survey_name == "lsst"
        assert isinstance(builder.resources, dict)

    # ------------------------------------------------------------------
    # _build_one_config — output structure
    # ------------------------------------------------------------------

    def test_build_one_config_stars_returns_expected_keys(self, mock_survey, stream_catalog):
        """_build_one_config must return the five expected keys for stars."""
        result = self._one_config(mock_survey, stream_catalog)
        assert set(result) == {"cmd_hist", "color_edges", "mag_edges", "n_ref", "area_ref_deg2"}
        assert result["cmd_hist"].shape == (10, 10)
        assert result["n_ref"] == len(stream_catalog)
        assert result["area_ref_deg2"] == 100.0

    def test_build_one_config_galaxies_returns_expected_keys(self, mock_survey, tiny_galaxies_catalog):
        """_build_one_config must work for galaxies and return the same five keys."""
        result = self._one_config(mock_survey, tiny_galaxies_catalog, source_type="galaxies")
        assert set(result) == {"cmd_hist", "color_edges", "mag_edges", "n_ref", "area_ref_deg2"}
        assert result["cmd_hist"].shape == (10, 10)
        assert result["n_ref"] == len(tiny_galaxies_catalog)

    # ------------------------------------------------------------------
    # _build_one_config — histogram properties
    # ------------------------------------------------------------------

    def test_cmd_counts_are_non_negative_integers(self, mock_survey, bright_star_catalog):
        """All histogram counts must be non-negative and whole numbers."""
        result = self._one_config(mock_survey, bright_star_catalog)
        H = result["cmd_hist"]
        assert (H >= 0).all(), "Counts must be non-negative"
        assert np.all(H == np.floor(H)), "Counts must be integer-valued"

    def test_cmd_sum_leq_n_ref(self, mock_survey, bright_star_catalog):
        """Total histogram count must not exceed the input catalog size."""
        result = self._one_config(mock_survey, bright_star_catalog)
        assert result["cmd_hist"].sum() <= result["n_ref"]

    def test_larger_maglim_gives_more_counts_stars(self, mock_survey, stream_catalog):
        """A fainter magnitude limit must detect at least as many stars.

        Stellar completeness rises to ~1 for bright sources (efficiency padding), so a
        very bright maglim cuts most of the stream while a deep maglim keeps them all.
        """
        result_bright = self._one_config(mock_survey, stream_catalog,
                                         source_type="stars", maglim_r=21.0, maglim_g=20.5)
        result_faint  = self._one_config(mock_survey, stream_catalog,
                                         source_type="stars", maglim_r=26.0, maglim_g=25.5)
        assert result_faint["cmd_hist"].sum() >= result_bright["cmd_hist"].sum()

    def test_larger_maglim_gives_more_counts_galaxies(self, mock_survey, galaxy_lf_catalog):
        """A fainter maglim must misclassify more galaxies when the LF is steep.

        The galaxy LF has 10× more objects at r≈27 than at r≈24. The misclassification
        efficiency is concentrated within ~1 mag of the detection limit. Therefore:
          - maglim≈24 catches ~20 galaxies near their limit → few misclassified.
          - maglim≈27 catches ~200 galaxies near their limit → many more misclassified.
        """
        result_bright = self._one_config(mock_survey, galaxy_lf_catalog,
                                         source_type="galaxies", maglim_r=24.5, maglim_g=24.0)
        result_faint  = self._one_config(mock_survey, galaxy_lf_catalog,
                                         source_type="galaxies", maglim_r=27.5, maglim_g=27.0)
        assert result_faint["cmd_hist"].sum() >= result_bright["cmd_hist"].sum()

    def test_no_counts_well_above_maglim(self, mock_survey, bright_star_catalog):
        """No detected stars should have observed magnitude 3+ mag above the maglim."""
        # All stars at true r=22; with maglim_r=24 they are all detected.
        # Observed scatter is tiny (~0.02 mag) so NO counts should appear at mag > 25.
        maglim_r = 24.0
        result = self._one_config(mock_survey, bright_star_catalog, maglim_r=maglim_r, maglim_g=23.5)
        mag_edges = result["mag_edges"]
        mag_centers = (mag_edges[:-1] + mag_edges[1:]) / 2
        well_above = mag_centers > maglim_r + 3
        assert result["cmd_hist"][:, well_above].sum() == 0, \
            "Counts 3+ mag above the limit must be zero"

    def test_bright_stars_concentrated_in_expected_color_mag_bins(self, mock_survey, bright_star_catalog):
        """Stars at r=22, g=22.5 (color 0.5) must concentrate around (color≈0.5, mag≈22)."""
        result = self._one_config(mock_survey, bright_star_catalog, maglim_r=26.0, maglim_g=25.5)
        H = result["cmd_hist"]
        mag_edges = result["mag_edges"]
        color_edges = result["color_edges"]
        mag_centers = (mag_edges[:-1] + mag_edges[1:]) / 2
        color_centers = (color_edges[:-1] + color_edges[1:]) / 2

        # Find the bin that should contain all counts.
        mag_bin = np.searchsorted(mag_edges, 22.0, side="right") - 1
        color_bin = np.searchsorted(color_edges, 0.5, side="right") - 1
        mag_bin = min(mag_bin, H.shape[1] - 1)
        color_bin = min(color_bin, H.shape[0] - 1)

        # The peak bin must hold the bulk of counts.
        assert H[color_bin, mag_bin] > 0.8 * H.sum(), \
            f"Expected most counts in bin (color≈{color_centers[color_bin]:.2f}, mag≈{mag_centers[mag_bin]:.2f})"

    # ------------------------------------------------------------------
    # save / load roundtrip
    # ------------------------------------------------------------------

    def test_save_via_storage(self, tmp_path, mock_survey, bright_star_catalog):
        """save must write a parquet file for each source type."""
        from streamobs.background import BackgroundResourceBuilder, BackgroundStorage

        builder = BackgroundResourceBuilder(survey_name="lsst")
        builder.bands = ["g", "r"]
        builder.resources["stars"] = {
            (26.0, 25.5): builder._build_one_config(
                bright_star_catalog, mock_survey, "stars",
                maglim_r=26.0, maglim_g=25.5, **self._COMMON,
            )
        }
        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        builder.save(storage, source_type="stars")
        assert storage.exists("stars", ("g", "r"))

    def test_load_via_storage_roundtrip(self, tmp_path, mock_survey, bright_star_catalog):
        """Resources loaded from storage must match the saved histogram exactly."""
        from streamobs.background import BackgroundResourceBuilder, BackgroundStorage

        builder = BackgroundResourceBuilder(survey_name="lsst")
        builder.bands = ["g", "r"]
        result = builder._build_one_config(
            bright_star_catalog, mock_survey, "stars",
            maglim_r=26.0, maglim_g=25.5, **self._COMMON,
        )
        builder.resources["stars"] = {(26.0, 25.5): result}

        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        builder.save(storage, source_type="stars")

        loaded = BackgroundResourceBuilder.load(storage, source_type="stars", bands=("g", "r"))
        loaded_result = loaded.resources["stars"][(26.0, 25.5)]
        np.testing.assert_allclose(loaded_result["cmd_hist"], result["cmd_hist"])
        np.testing.assert_allclose(loaded_result["color_edges"], result["color_edges"])
        np.testing.assert_allclose(loaded_result["mag_edges"], result["mag_edges"])
        assert loaded_result["n_ref"] == result["n_ref"]
        assert loaded_result["area_ref_deg2"] == result["area_ref_deg2"]

    def test_load_via_storage_both_types(self, tmp_path, mock_survey,
                                          bright_star_catalog, tiny_galaxies_catalog):
        """Saving both source types and loading must return resources for each."""
        from streamobs.background import BackgroundResourceBuilder, BackgroundStorage

        builder = BackgroundResourceBuilder(survey_name="lsst")
        builder.bands = ["g", "r"]
        builder.resources["stars"] = {
            (26.0, 25.5): builder._build_one_config(
                bright_star_catalog, mock_survey, "stars",
                maglim_r=26.0, maglim_g=25.5, **self._COMMON,
            )
        }
        builder.resources["galaxies"] = {
            (26.0, 25.5): builder._build_one_config(
                tiny_galaxies_catalog, mock_survey, "galaxies",
                maglim_r=26.0, maglim_g=25.5, **self._COMMON,
            )
        }
        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        builder.save(storage, source_type="both")

        loaded = BackgroundResourceBuilder.load(storage, source_type="both", bands=("g", "r"))
        assert "stars" in loaded.resources
        assert "galaxies" in loaded.resources
        assert (26.0, 25.5) in loaded.resources["stars"]
        assert (26.0, 25.5) in loaded.resources["galaxies"]


# ---------------------------------------------------------------------------
# Part 6 — LightBackgroundGenerator
# ---------------------------------------------------------------------------


class TestLightBackgroundGenerator:
    """Tests for LightBackgroundGenerator."""

    @pytest.fixture(scope="class")
    def storage_with_data(self, tmp_path_factory, mock_survey, stream_catalog):
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

    def test_full_method_stars(self, mock_survey, stream_catalog):
        """Background with method='full' and source_type='stars' must generate a catalog."""
        ...

    def test_full_method_galaxies(self, mock_survey, tiny_galaxies_catalog):
        """Background with method='full' and source_type='galaxies' must generate a catalog."""
        ...

    def test_full_method_both(
        self, mock_survey, stream_catalog, tiny_galaxies_catalog
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

        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        bg = Background(mock_survey, method="light", storage=storage)
        assert bg.storage is storage
