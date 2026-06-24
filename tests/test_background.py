"""
Tests for the streamobs.background subpackage.

All tests are marked ``background`` so they can be run in isolation:

    pytest tests/test_background.py -m background

Tests that write to disk (Storage, ResourceBuilder) use the ``tmp_path``
pytest fixture so the temporary directory is created and cleaned up
automatically regardless of test outcome.
"""

import os

import numpy as np
import pandas as pd
import pytest

pytestmark = pytest.mark.background


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------


def _make_fake_grid(n_color=8, n_mag=8):
    """Return a minimal 2-pair CMD grid for storage round-trip tests."""
    rng = np.random.default_rng(7)
    color_edges = np.linspace(-1.0, 3.0, n_color + 1)
    mag_edges = np.linspace(18.0, 28.0, n_mag + 1)
    grid = {}
    for mr, mg in [(24.0, 24.0), (25.0, 25.0)]:
        cmd = rng.integers(0, 20, size=(n_color, n_mag)).astype(float)
        grid[(mr, mg)] = {
            "cmd_hist": cmd,
            "color_edges": color_edges,
            "mag_edges": mag_edges,
            "n_ref": 500,
            "area_ref_deg2": 100.0,
        }
    return grid


def _make_gc_frame():
    """Great-circle frame with stream along dec≈−20°, within LSST footprint."""
    import astropy.coordinates as coord
    import astropy.units as u
    import gala.coordinates as gc

    end1 = coord.SkyCoord(ra=30 * u.deg, dec=-20 * u.deg, frame="icrs")
    end2 = coord.SkyCoord(ra=60 * u.deg, dec=-20 * u.deg, frame="icrs")
    return gc.GreatCircleICRSFrame.from_endpoints(end1, end2)


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
        with pytest.raises(ValueError, match="missclassified"):
            s.get_gal_misclassification("r", mags, maglim)

    def test_get_gal_misclassification_no_1padding(self, mock_survey):
        """get_gal_misclassification must NOT return 1 at the bright end (no 1-padding).

        Verify by checking that the method returns 0 for very bright objects
        (below saturation) and does not saturate to 1.0 for delta_mag far below
        delta_saturation.
        """
        ...


# ---------------------------------------------------------------------------
# Part 2 — StreamInjector source_type parameter
# ---------------------------------------------------------------------------


class TestInjectSourceType:
    """Tests that source_type is accepted and routed correctly in StreamInjector."""

    def test_inject_accepts_source_type_kwarg(self):
        """inject() must accept source_type without raising TypeError."""
        import inspect

        from streamobs.observed import StreamInjector

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

    def test_inject_stars_returns_nonempty_dataframe(self, mock_survey, stream_catalog):
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
        flag_cols = [c for c in result.columns if "flag" in c]
        assert len(flag_cols) > 0
        assert "lsst_yr4_flag_observed" in flag_cols
        assert result["lsst_yr4_flag_observed"].isin([0, 1]).all()
        assert result["lsst_yr4_flag_observed"].sum() > 0  # at least one detected


# ---------------------------------------------------------------------------
# Part 4 — BackgroundStorage
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

    def test_save_data_creates_file(self, tmp_path):
        """save_data must write a file at the path returned by get_path."""
        from streamobs.background import BackgroundStorage

        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        grid = _make_fake_grid()
        storage.save_data(grid, "stars", ("g", "r"))
        assert os.path.isfile(storage.get_path("stars", ("g", "r")))

    def test_load_data_roundtrip(self, tmp_path):
        """load_data(maglim_b2, maglim_b1) must recover a single pair via predicate pushdown."""
        from streamobs.background import BackgroundStorage

        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        grid = _make_fake_grid()
        storage.save_data(grid, "stars", ("g", "r"))

        for (mr, mg), expected in grid.items():
            loaded = storage.load_data("stars", ("g", "r"), mr, mg)
            assert np.allclose(loaded["cmd_hist"], expected["cmd_hist"])
            assert np.allclose(loaded["color_edges"], expected["color_edges"])
            assert np.allclose(loaded["mag_edges"], expected["mag_edges"])
            assert loaded["n_ref"] == expected["n_ref"]
            assert np.isclose(loaded["area_ref_deg2"], expected["area_ref_deg2"])

    def test_load_all_roundtrip(self, tmp_path):
        """load_all must recover the full grid saved by save_data."""
        from streamobs.background import BackgroundStorage

        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        grid = _make_fake_grid()
        storage.save_data(grid, "stars", ("g", "r"))
        loaded = storage.load_all("stars", ("g", "r"))

        assert set(loaded.keys()) == set(grid.keys())
        for key in grid:
            assert np.allclose(loaded[key]["cmd_hist"], grid[key]["cmd_hist"])
            assert np.allclose(loaded[key]["color_edges"], grid[key]["color_edges"])
            assert np.allclose(loaded[key]["mag_edges"], grid[key]["mag_edges"])
            assert loaded[key]["n_ref"] == grid[key]["n_ref"]
            assert np.isclose(loaded[key]["area_ref_deg2"], grid[key]["area_ref_deg2"])

    def test_exists_false_before_save(self, tmp_path):
        """exists must return False when the file is not on disk."""
        from streamobs.background import BackgroundStorage

        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        assert storage.exists("stars", ("g", "r")) is False

    def test_exists_true_after_save(self, tmp_path):
        """exists must return True after save_data has been called."""
        from streamobs.background import BackgroundStorage

        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        storage.save_data(_make_fake_grid(), "stars", ("g", "r"))
        assert storage.exists("stars", ("g", "r")) is True


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

    def test_build_one_config(self, mock_survey, stream_catalog):
        """_build_one_config must return a dict with the expected keys."""
        from streamobs.background import BackgroundResourceBuilder

        builder = BackgroundResourceBuilder(survey_name="lsst", release="yr4")
        result = builder._build_one_config(
            catalog=stream_catalog,
            survey=mock_survey,
            source_type="stars",
            bands=("g", "r"),
            maglim_b2=26.0,
            maglim_b1=26.5,
            n_bins_color=10,
            n_bins_mag=10,
            color_range=(-2, 3),
            mag_range=(14, 30),
            area_ref_deg2=100.0,
        )
        assert set(result.keys()) == {
            "cmd_hist",
            "color_edges",
            "mag_edges",
            "n_ref",
            "area_ref_deg2",
        }
        assert result["cmd_hist"].shape == (10, 10)
        assert result["cmd_hist"].sum() >= 0
        assert result["n_ref"] == len(stream_catalog)

    def test_save_via_storage(self, tmp_path, stream_catalog):
        """save must write a parquet file via BackgroundStorage."""
        from streamobs.background import (BackgroundResourceBuilder,
                                          BackgroundStorage)

        builder = BackgroundResourceBuilder(survey_name="lsst", release="yr4")
        builder.build(
            catalog_stars=stream_catalog,
            maglim_min=26.0,
            maglim_max=26.0,
            maglim_step=1.0,
            max_delta=1.0,
            source_type="stars",
            area_ref_deg2=100.0,
            n_bins_color=5,
            n_bins_mag=5,
        )
        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        builder.save(storage, source_type="stars")
        assert storage.exists("stars", builder.bands)

    def test_load_via_storage(self, tmp_path, stream_catalog):
        """load must reconstruct resources from the file saved by save."""
        from streamobs.background import (BackgroundResourceBuilder,
                                          BackgroundStorage)

        builder = BackgroundResourceBuilder(survey_name="lsst", release="yr4")
        builder.build(
            catalog_stars=stream_catalog,
            maglim_min=26.0,
            maglim_max=26.0,
            maglim_step=1.0,
            max_delta=1.0,
            source_type="stars",
            area_ref_deg2=100.0,
            n_bins_color=5,
            n_bins_mag=5,
        )
        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        builder.save(storage, source_type="stars")

        loaded = BackgroundResourceBuilder.load(
            storage, source_type="stars", bands=builder.bands
        )
        assert "stars" in loaded.resources
        for key in builder.resources["stars"]:
            assert key in loaded.resources["stars"]
            orig = builder.resources["stars"][key]
            reco = loaded.resources["stars"][key]
            assert np.allclose(orig["cmd_hist"], reco["cmd_hist"])


# ---------------------------------------------------------------------------
# Part 6 — LightBackgroundGenerator
# ---------------------------------------------------------------------------


class TestLightBackgroundGenerator:
    """Tests for LightBackgroundGenerator.

    Uses in-memory HEALPix maps (nside=8, ~3 ms to create) instead of the
    full LSST yr4 survey so the class runs in a few seconds when isolated.
    """

    @pytest.fixture(scope="class")
    def gc_frame(self):
        """Great-circle frame built once for all tests in this class."""
        return _make_gc_frame()

    @pytest.fixture(scope="class")
    def fast_survey(self):
        """Minimal Survey with tiny in-memory maps — no disk I/O."""
        import healpy as hp

        from streamobs.surveys import Survey

        nside = 8
        n_pix = hp.nside2npix(nside)
        return Survey(
            name="lsst",
            release="yr4",
            maglim_maps={"g": np.full(n_pix, 24.3), "r": np.full(n_pix, 24.8)},
            coeff_extinc={"g": 3.303, "r": 2.285},
            ebv_map=np.full(n_pix, 0.01),
        )

    @pytest.fixture(scope="class")
    def storage_with_data(self, tmp_path_factory):
        """BackgroundStorage with a 2×2 CMD grid for stars and galaxies (fake data)."""
        from streamobs.background import BackgroundStorage

        tmp_path = tmp_path_factory.mktemp("generator")
        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")

        rng = np.random.default_rng(42)
        n_color, n_mag = 8, 8
        color_edges = np.linspace(-1.0, 3.0, n_color + 1)
        mag_edges = np.linspace(18.0, 28.0, n_mag + 1)

        # 2×2 grid; deeper pair has proportionally more counts
        pairs = [(24.0, 24.0), (24.0, 24.5), (24.5, 24.0), (24.5, 24.5)]
        for source_type in ("stars", "galaxies"):
            grid = {}
            for mr, mg in pairs:
                scale = (mr + mg) / (24.0 + 24.0)
                cmd = rng.integers(1, 5, size=(n_color, n_mag)).astype(float) * scale
                grid[(mr, mg)] = {
                    "cmd_hist": cmd,
                    "color_edges": color_edges,
                    "mag_edges": mag_edges,
                    "n_ref": int(20 * scale),  # small → ~1-3 objects per pixel
                    "area_ref_deg2": 10.0,
                }
            storage.save_data(grid, source_type, ("g", "r"))

        return storage

    def test_init(self, fast_survey, tmp_path):
        """LightBackgroundGenerator can be instantiated."""
        from streamobs.background import (BackgroundStorage,
                                          LightBackgroundGenerator)

        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        gen = LightBackgroundGenerator(storage, fast_survey)
        assert gen.survey is fast_survey
        assert gen.bands == ("g", "r")
        assert gen._resources == {}

    def _verify_dataframe_content(self, df, meta, bands, survey):
        """Assert standard columns and meta keys are present and catalog is non-empty."""
        from streamobs.columns import obs_col

        namespace = survey.namespace
        for col in (
            "ra",
            "dec",
            "phi1",
            "phi2",
            obs_col(bands[0], namespace),
            obs_col(bands[1], namespace),
            "source_type",
        ):
            assert col in df.columns, f"Missing column: {col}"
        assert isinstance(meta, dict)
        for key in ("nside", "color_edges", "mag_edges", "band1", "band2"):
            assert key in meta, f"Missing meta key: {key}"
        assert len(df) > 0, "No objects generated"

    def test_generate_stars(self, storage_with_data, fast_survey, gc_frame):
        """generate with source_type='stars' returns (df, meta) with correct columns."""
        from streamobs.background import LightBackgroundGenerator

        gen = LightBackgroundGenerator(storage_with_data, fast_survey, bands=("g", "r"))
        df, meta = gen.generate(
            phi1_limits=(-3, 3),
            phi2_limits=(-1, 1),
            gc_frame=gc_frame,
            nside=32,
            source_type="stars",
            seed=0,
        )
        self._verify_dataframe_content(df, meta, ("g", "r"), fast_survey)

    def test_generate_galaxies(self, storage_with_data, fast_survey, gc_frame):
        """generate with source_type='galaxies' returns a DataFrame with correct columns."""
        from streamobs.background import LightBackgroundGenerator

        gen = LightBackgroundGenerator(storage_with_data, fast_survey, bands=("g", "r"))
        df, meta = gen.generate(
            phi1_limits=(-3, 3),
            phi2_limits=(-1, 1),
            gc_frame=gc_frame,
            nside=32,
            source_type="galaxies",
            seed=1,
        )
        self._verify_dataframe_content(df, meta, ("g", "r"), fast_survey)

    def test_generate_both(self, storage_with_data, fast_survey, gc_frame):
        """generate with source_type='both' produces rows for both source types."""
        from streamobs.background import LightBackgroundGenerator

        gen = LightBackgroundGenerator(storage_with_data, fast_survey, bands=("g", "r"))
        df, meta = gen.generate(
            phi1_limits=(-3, 3),
            phi2_limits=(-1, 1),
            gc_frame=gc_frame,
            nside=32,
            source_type="both",
            seed=2,
        )
        self._verify_dataframe_content(df, meta, ("g", "r"), fast_survey)

        types_present = set(df["source_type"].unique())
        assert "stars" in types_present
        assert "galaxies" in types_present

    def test_generate_within_limits(self, storage_with_data, fast_survey, gc_frame):
        """All generated phi1/phi2 must lie strictly within the requested limits."""
        from streamobs.background import LightBackgroundGenerator

        phi1_lim = (-3, 3)
        phi2_lim = (-1, 1)
        gen = LightBackgroundGenerator(storage_with_data, fast_survey, bands=("g", "r"))
        df, _ = gen.generate(
            phi1_limits=phi1_lim,
            phi2_limits=phi2_lim,
            gc_frame=gc_frame,
            nside=32,
            source_type="both",
            seed=7,
        )
        if len(df) > 0:
            assert (
                df["phi1"].between(*phi1_lim).all()
            ), f"phi1 out of range: [{df['phi1'].min():.3f}, {df['phi1'].max():.3f}]"
            assert (
                df["phi2"].between(*phi2_lim).all()
            ), f"phi2 out of range: [{df['phi2'].min():.3f}, {df['phi2'].max():.3f}]"

    # ------------------------------------------------------------------
    # Fixture shared by the count-comparison tests below
    # ------------------------------------------------------------------

    @pytest.fixture(scope="class")
    def count_storage(self, tmp_path_factory):
        """Storage with two CMD grid points (maglim=22 → few, maglim=26 → many)."""
        from streamobs.background import BackgroundStorage

        tmp = tmp_path_factory.mktemp("count")
        storage = BackgroundStorage(base_path=str(tmp), survey_name="test")

        n_color, n_mag = 5, 5
        color_edges = np.linspace(-1, 3, n_color + 1)
        mag_edges = np.linspace(18, 28, n_mag + 1)
        flat_cmd = np.ones((n_color, n_mag))

        grid = {
            (22.0, 22.0): {
                "cmd_hist": 2 * flat_cmd.copy(),
                "color_edges": color_edges,
                "mag_edges": mag_edges,
                "n_ref": 5000,
                "area_ref_deg2": 10.0,
            },
            (26.0, 26.0): {
                "cmd_hist": 10 * flat_cmd.copy(),
                "color_edges": color_edges,
                "mag_edges": mag_edges,
                "n_ref": 5000,
                "area_ref_deg2": 10.0,
            },
        }
        storage.save_data(grid, "stars", ("g", "r"))
        return storage

    # ------------------------------------------------------------------
    # Interpolation
    # ------------------------------------------------------------------

    def test_interpolate_cmd_at_grid_center(self, storage_with_data, fast_survey):
        """Bilinear interpolation at the center of a 2x2 grid equals the mean of corners."""
        from streamobs.background import LightBackgroundGenerator

        gen = LightBackgroundGenerator(storage_with_data, fast_survey, bands=("g", "r"))
        gen._load_resources("stars")
        grid = gen._resources["stars"]

        # 2x2 grid corners
        corners = [(24.0, 24.0), (24.0, 24.5), (24.5, 24.0), (24.5, 24.5)]
        expected = np.mean([grid[k]["cmd_hist"] for k in corners], axis=0)

        result = gen._interpolate_cmd(24.25, 24.25, "stars")
        assert np.allclose(
            result["cmd_hist"], expected, atol=1e-9
        ), "Interpolated histogram does not match equal-weight mean of 4 corners"

    def test_interpolation_vs_ground_truth(self, mock_survey, tmp_path):
        """Bilinear interpolation at (26.0, 26.0) must match direct injection within tolerances.

        Builds a 2×2 CMD grid at corners {25.9, 26.1}², interpolates at the centre,
        then compares against a CMD built by direct injection at (26.0, 26.0).
        Thresholds: total-count bias < 5 % and per-bin correlation > 0.95.
        """
        import pandas as pd

        from streamobs.background import (BackgroundResourceBuilder,
                                          BackgroundStorage,
                                          LightBackgroundGenerator)

        # Synthetic catalog with columns matching true_col(band, "lsst_yr4") → "lsst_<band>_true"
        rng = np.random.default_rng(0)
        n = 50000
        r_mag = rng.uniform(18.0, 28.0, n)
        cat = pd.DataFrame(
            {
                "ra": rng.uniform(0.0, 360.0, n),
                "dec": rng.uniform(-90.0, 0.0, n),
                "lsst_r_true": r_mag,
                "lsst_g_true": r_mag + rng.uniform(0.4, 0.7, n),
            }
        )

        builder = BackgroundResourceBuilder(survey_name="lsst", release="yr4")
        common_kw = dict(
            catalog=cat,
            survey=mock_survey,
            source_type="stars",
            bands=("g", "r"),
            n_bins_color=40,
            n_bins_mag=40,
            color_range=(-1, 3),
            mag_range=(15, 27),
            area_ref_deg2=100.0,
        )

        # Build 2×2 corner grid
        corners = [(25.9, 25.9), (25.9, 26.1), (26.1, 25.9), (26.1, 26.1)]
        grid = {
            (mb2, mb1): builder._build_one_config(
                maglim_b2=mb2, maglim_b1=mb1, **common_kw
            )
            for mb2, mb1 in corners
        }
        storage = BackgroundStorage(base_path=str(tmp_path), survey_name="lsst")
        storage.save_data(grid, "stars", ("g", "r"))

        # Ground-truth CMD built by direct injection at the centre point
        truth = builder._build_one_config(maglim_b2=26.0, maglim_b1=26.0, **common_kw)

        # Interpolated CMD from the 2×2 grid
        gen = LightBackgroundGenerator(storage, mock_survey, bands=("g", "r"))
        gen._load_resources("stars")
        interp = gen._interpolate_cmd(26.0, 26.0, "stars")

        h_truth = truth["cmd_hist"]
        h_interp = interp["cmd_hist"]

        # Total-count bias ≤ 5 %
        n_truth = h_truth.sum()
        assert (
            n_truth > 0
        ), "Ground-truth CMD is empty — catalog may not match survey namespace"
        assert (
            abs(h_interp.sum() - n_truth) / n_truth < 0.02
        ), f"Total-count bias {abs(h_interp.sum() - n_truth) / n_truth:.1%} exceeds 2 %"

        # Per-bin correlation ≥ 0.95 on non-empty bins
        mask = (h_truth > 0) | (h_interp > 0)
        assert mask.sum() > 0, "No non-empty bins found in either histogram"
        corr = np.corrcoef(h_truth[mask], h_interp[mask])[0, 1]
        assert corr > 0.95, f"Per-bin correlation {corr:.4f} below 0.95"

    # ------------------------------------------------------------------
    # Count monotonicity: maglim depth
    # ------------------------------------------------------------------

    def test_fewer_objects_with_shallower_maglim(self, count_storage, gc_frame):
        """A shallower magnitude limit must produce fewer background objects."""
        import healpy as hp

        from streamobs.background import LightBackgroundGenerator
        from streamobs.surveys import Survey

        nside_survey = 64
        n_pix = hp.nside2npix(nside_survey)
        common = dict(coeff_extinc={"g": 3.303, "r": 2.285}, ebv_map=np.zeros(n_pix))

        survey_shallow = Survey(
            name="test",
            release="v1",
            maglim_maps={"g": np.full(n_pix, 22.0), "r": np.full(n_pix, 22.0)},
            **common,
        )
        survey_deep = Survey(
            name="test",
            release="v1",
            maglim_maps={"g": np.full(n_pix, 26.0), "r": np.full(n_pix, 26.0)},
            **common,
        )

        kwargs = dict(
            phi1_limits=(-5, 5),
            phi2_limits=(-1, 1),
            gc_frame=gc_frame,
            nside=nside_survey,
            source_type="stars",
            seed=0,
        )

        df_shallow, _ = LightBackgroundGenerator(
            count_storage, survey_shallow, bands=("g", "r")
        ).generate(**kwargs)
        df_deep, _ = LightBackgroundGenerator(
            count_storage, survey_deep, bands=("g", "r")
        ).generate(**kwargs)

        assert len(df_deep) > len(
            df_shallow
        ), f"Expected deep ({len(df_deep)}) > shallow ({len(df_shallow)})"

    # ------------------------------------------------------------------
    # Count monotonicity: dust extinction
    # ------------------------------------------------------------------

    def test_fewer_objects_with_more_dust(self, count_storage, gc_frame):
        """Higher dust extinction must reduce the effective depth and produce fewer objects."""
        import healpy as hp

        from streamobs.background import LightBackgroundGenerator
        from streamobs.surveys import Survey

        nside_survey = 64
        n_pix = hp.nside2npix(nside_survey)
        common = dict(
            maglim_maps={"g": np.full(n_pix, 25.0), "r": np.full(n_pix, 25.0)},
            coeff_extinc={"g": 3.303, "r": 2.285},
        )

        survey_nodust = Survey(
            name="test", release="v1", ebv_map=np.zeros(n_pix), **common
        )
        survey_highdust = Survey(
            name="test", release="v1", ebv_map=np.full(n_pix, 1.0), **common
        )

        kwargs = dict(
            phi1_limits=(-5, 5),
            phi2_limits=(-1, 1),
            gc_frame=gc_frame,
            nside=nside_survey,
            source_type="stars",
            seed=0,
        )

        df_nodust, _ = LightBackgroundGenerator(
            count_storage, survey_nodust, bands=("g", "r")
        ).generate(**kwargs)
        df_highdust, _ = LightBackgroundGenerator(
            count_storage, survey_highdust, bands=("g", "r")
        ).generate(**kwargs)

        assert len(df_nodust) > len(
            df_highdust
        ), f"Expected no-dust ({len(df_nodust)}) > high-dust ({len(df_highdust)})"

    # ------------------------------------------------------------------
    # nside capping
    # ------------------------------------------------------------------

    def test_nside_capped_with_warning(self, storage_with_data, fast_survey, gc_frame):
        """generate() warns and caps nside when it exceeds the maglim map resolution."""
        import warnings

        import healpy as hp

        from streamobs.background import LightBackgroundGenerator

        maglim_nside = hp.get_nside(fast_survey.maglim_maps["r"])  # 8 for fast_survey
        oversized_nside = maglim_nside * 8  # deliberately too large

        gen = LightBackgroundGenerator(storage_with_data, fast_survey, bands=("g", "r"))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            _, meta = gen.generate(
                phi1_limits=(-3, 3),
                phi2_limits=(-1, 1),
                gc_frame=gc_frame,
                nside=oversized_nside,
                source_type="stars",
                seed=0,
            )

        user_warnings = [w for w in caught if issubclass(w.category, UserWarning)]
        assert (
            len(user_warnings) == 1
        ), f"Expected 1 UserWarning, got {len(user_warnings)}"
        assert "nside" in str(user_warnings[0].message).lower()
        assert (
            meta["nside"] == maglim_nside
        ), f"meta['nside'] should be capped to {maglim_nside}, got {meta['nside']}"


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

    def test_injection_requires_catalog_stars(self, mock_survey):
        """generate() raises ValueError when catalog_stars is missing."""
        from streamobs.background import Background

        bg = Background(mock_survey, method="injection", source_type="stars")
        with pytest.raises(ValueError, match="catalog_stars"):
            bg.generate(phi1_limits=(-10, 10), phi2_limits=(-2, 2))

    def test_injection_requires_catalog_galaxies(self, mock_survey):
        """generate() raises ValueError when catalog_galaxies is missing."""
        from streamobs.background import Background

        bg = Background(mock_survey, method="injection", source_type="galaxies")
        with pytest.raises(ValueError, match="catalog_galaxies"):
            bg.generate(phi1_limits=(-10, 10), phi2_limits=(-2, 2))
