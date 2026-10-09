"""
tests/test_model.py
===================
Tests for ``streamobs.model.StreamModel`` and its sub-models.

Focus: verify that every quantity the model can produce is actually sampled
and has the right shape / dtype, without requiring the optional ugali.
"""

import numpy as np
import pandas as pd
import pytest

from streamobs.model import (DensityModel, IsochroneModel, StreamModel,
                             TrackModel)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

N = 200  # number of stars used in every sampling test


# ---------------------------------------------------------------------------
# DensityModel
# ---------------------------------------------------------------------------


@pytest.mark.model
class TestDensityModel:
    """Unit tests for DensityModel (phi1 sampler)."""

    def _verify_sampling(self, model):
        """Helper to verify that sampling produces a 1D array of the right length."""
        samples = model.sample(N)
        assert samples.shape == (N,), "Sample array must have length N"
        assert np.issubdtype(samples.dtype, np.floating), "Samples must be floats"
        assert np.all(np.isfinite(samples)), "Samples must be finite numbers"
        assert np.all(~np.isnan(samples)), "Samples must not contain NaNs"
        return samples

    def test_uniform_sample(self):
        cfg = {"type": "uniform", "xmin": -10.0, "xmax": 10.0}
        model = DensityModel(cfg)
        samples = self._verify_sampling(model)
        assert np.all(samples > cfg["xmin"]) and np.all(
            samples < cfg["xmax"]
        ), "Uniform samples must be within [xmin, xmax]"

    def test_gaussian_sample_shape(self):
        cfg = {"type": "gaussian", "mu": 0.0, "sigma": 1.0}
        model = DensityModel(cfg)
        samples = self._verify_sampling(model)


# ---------------------------------------------------------------------------
# TrackModel
# ---------------------------------------------------------------------------
@pytest.mark.model
class TestTrackModel:
    """Unit tests for TrackModel (phi2 sampler given phi1)."""

    def _verify_sampling(self, model, phi1):
        """Helper to verify that sampling produces a 1D array of the right length."""
        samples = model.sample(phi1)
        assert samples.shape == (N,), "Sample array must have length N"
        assert np.issubdtype(samples.dtype, np.floating), "Samples must be floats"
        assert np.all(np.isfinite(samples)), "Samples must be finite numbers"
        assert np.all(~np.isnan(samples)), "Samples must not contain NaNs"
        return samples

    def test_constant_track(self):
        cfg = {
            "center": {"type": "constant", "value": 0.0},
            "spread": {"type": "constant", "value": 0.2},
            "sampler": "gaussian",
        }
        model = TrackModel(cfg)
        phi1 = np.linspace(-10, 10, N)
        samples = self._verify_sampling(model, phi1)
        assert np.all(
            np.abs(samples) < 10 * cfg["spread"]["value"]
        ), "Samples should be within a few sigma of the center"

    def test_sinusoidal_track(self):
        cfg = {
            "center": {"type": "sinusoid", "amplitude": 0.5, "period": 2.0},
            "spread": {"type": "constant", "value": 0.2},
            "sampler": "gaussian",
        }
        model = TrackModel(cfg)
        phi1 = np.linspace(-10, 10, N)
        samples = self._verify_sampling(model, phi1)
        expected_center = cfg["center"]["amplitude"] * np.sin(
            phi1 * 2 * np.pi / cfg["center"]["period"]
        )
        assert np.all(
            np.abs(samples - expected_center) < 10 * cfg["spread"]["value"]
        ), "Samples should be within a few sigma of the sinusoidal center"


# ---------------------------------------------------------------------------
# StreamModel — full
# ---------------------------------------------------------------------------


@pytest.mark.model
class TestStreamModelFull:
    """Tests for StreamModel when having a complete config"""

    def _verify_catalogue_content(self, catalog, expected_columns):
        """Helper to verify that a completed catalog contains the expected columns with valid data."""
        assert expected_columns.issubset(
            catalog.columns
        ), f"Catalog should contain columns {expected_columns}"
        for col in expected_columns:
            assert np.issubdtype(
                catalog[col].dtype, np.floating
            ), f"Column {col} should contain floats"
            assert np.all(
                np.isfinite(catalog[col])
            ), f"Column {col} should contain finite numbers"
            assert np.all(
                ~np.isnan(catalog[col])
            ), f"Column {col} should not contain NaNs"

    def test_full_model(self, stream_config_with_distance):
        """Test that StreamModel can be instantiated and sampled with a full config."""
        model = StreamModel(stream_config_with_distance)
        samples = model.sample(N)
        assert isinstance(
            samples, pd.DataFrame
        ), "Samples should be returned as a DataFrame"
        expected_columns = {
            "phi1",
            "phi2",
            "dist",
            "lsst_g_true",
            "lsst_r_true",
        }  # Not adding mu1, mu2, rv since not implemented yet
        self._verify_catalogue_content(samples, expected_columns)

    def test_complete_catalog(self, sample_catalog_phi, stream_config_with_distance):
        """Test that complete_catalog produces a catalog with the expected columns and valid data."""
        model = StreamModel(stream_config_with_distance)
        completed_catalog = model.complete_catalog(catalog=sample_catalog_phi)
        expected_columns = {
            "phi1",
            "phi2",
            "dist",
            "lsst_g_true",
            "lsst_r_true",
        }  # Not adding mu1, mu2, rv since not implemented yet
        self._verify_catalogue_content(completed_catalog, expected_columns)
        assert len(completed_catalog) == len(
            sample_catalog_phi
        ), f"Completed catalog should have the same number of rows as the input catalog"

        # Verify that I can add a targeted column (e.g. dist) to the input catalog and complete the rest
        partial_catalog = completed_catalog.drop(
            columns=["lsst_r_true", "dist", "lsst_g_true"]
        ).reset_index(drop=True)
        completed_catalog = model.complete_catalog(
            catalog=partial_catalog,
            columns_to_add=["dist"],
        )
        expected_columns = {
            "phi1",
            "phi2",
            "dist",
        }  # Not adding mu1, mu2, rv since not implemented yet
        self._verify_catalogue_content(completed_catalog, expected_columns)
        assert len(completed_catalog) == len(
            partial_catalog
        ), f"Completed catalog should have the same number of rows as the input catalog"
        assert (
            "lsst_r_true" not in completed_catalog.columns
        ), "Column 'lsst_r_true' should not be added when not requested"


# ---------------------------------------------------------------------------
# StreamModel.complete_catalog — permutations of input columns / dist source
# ---------------------------------------------------------------------------


@pytest.mark.model
class TestCompleteCatalogPermutations:
    """Exercise complete_catalog across the ways a catalog can be partially filled.

    Covers: which columns are supplied (empty frame, phi-only, ra/dec+dist),
    where the distance comes from (distance_modulus model vs. a directly supplied
    scalar/vector ``dist``), and the preserve-existing-values contract (whole
    columns and individual NaN rows are never overwritten).
    """

    MAGS = {"lsst_g_true", "lsst_r_true"}

    def test_empty_frame_fills_all_model_columns(self, stream_config_with_distance):
        """size=N with no catalog -> geometry, dist, and both bands are filled."""
        model = StreamModel(stream_config_with_distance)
        out = model.complete_catalog(catalog=None, size=12, verbose=False)
        assert len(out) == 12
        assert ({"phi1", "phi2", "dist"} | self.MAGS).issubset(out.columns)
        assert out[list(self.MAGS)].notna().all().all()

    def test_phi_only_fills_dist_then_mags(
        self, sample_catalog_phi, stream_config_with_distance
    ):
        """phi1/phi2 present -> dist sampled from the model, then magnitudes."""
        model = StreamModel(stream_config_with_distance)
        out = model.complete_catalog(catalog=sample_catalog_phi.copy(), verbose=False)
        assert ({"dist"} | self.MAGS).issubset(out.columns)
        assert out[list(self.MAGS)].notna().all().all()

    def test_radec_plus_dist_fills_mags_without_phi(self, stream_config_with_distance):
        """With dist already present, magnitudes fill even when phi1 is absent."""
        model = StreamModel(stream_config_with_distance)
        df = pd.DataFrame(
            {"ra": [10.0, 11.0, 12.0], "dec": [-1.0, 0.0, 1.0], "dist": [16.0] * 3}
        )
        out = model.complete_catalog(
            catalog=df, columns_to_add=list(self.MAGS), verbose=False
        )
        assert self.MAGS.issubset(out.columns)
        assert out[list(self.MAGS)].notna().all().all()
        assert "phi1" not in out.columns

    def test_existing_band_preserved_when_filling_other(
        self, stream_config_with_distance
    ):
        """Providing one band and requesting both leaves the provided one intact."""
        model = StreamModel(stream_config_with_distance)
        g = np.array([20.0, 21.0, 22.0])
        df = pd.DataFrame({"dist": [16.0] * 3, "lsst_g_true": g.copy()})
        out = model.complete_catalog(
            catalog=df, columns_to_add=list(self.MAGS), verbose=False
        )
        assert np.allclose(out["lsst_g_true"].to_numpy(), g), "present band overwritten"
        assert out["lsst_r_true"].notna().all(), "missing band not filled"

    def test_present_columns_skip_sampling(self, stream_config_with_distance):
        """Both bands present -> values are returned untouched."""
        model = StreamModel(stream_config_with_distance)
        g = np.array([20.0, 21.0])
        r = np.array([19.0, 19.5])
        df = pd.DataFrame(
            {"dist": [16.0, 16.0], "lsst_g_true": g.copy(), "lsst_r_true": r.copy()}
        )
        out = model.complete_catalog(
            catalog=df, columns_to_add=list(self.MAGS), verbose=False
        )
        assert np.allclose(out["lsst_g_true"].to_numpy(), g)
        assert np.allclose(out["lsst_r_true"].to_numpy(), r)

    def test_partial_rows_only_missing_filled(self, stream_config_with_distance):
        """A band with some NaN rows keeps its finite rows; only NaNs are filled."""
        model = StreamModel(stream_config_with_distance)
        g = np.array([20.0, np.nan, 22.0])
        df = pd.DataFrame({"dist": [16.0] * 3, "lsst_g_true": g.copy()})
        out = model.complete_catalog(
            catalog=df, columns_to_add=["lsst_g_true"], verbose=False
        )
        filled = out["lsst_g_true"].to_numpy()
        assert filled[0] == 20.0 and filled[2] == 22.0, "finite rows overwritten"
        assert np.isfinite(filled[1]), "NaN row not filled"

    def test_dist_scalar_broadcast_without_distance_model(self, minimal_stream_config):
        """A scalar `dist` lets mags fill with no distance_modulus model / no phi1."""
        model = StreamModel(minimal_stream_config)  # no distance_modulus section
        df = pd.DataFrame({"ra": [10.0, 11.0], "dec": [0.0, 1.0]})
        out = model.complete_catalog(
            catalog=df,
            columns_to_add=["dist"] + list(self.MAGS),
            dist=16.5,
            verbose=False,
        )
        assert np.allclose(out["dist"].to_numpy(), 16.5)
        assert out[list(self.MAGS)].notna().all().all()

    def test_dist_vector_assigned_per_row(self, minimal_stream_config):
        """A per-row `dist` vector is assigned row-wise."""
        model = StreamModel(minimal_stream_config)
        df = pd.DataFrame({"ra": [10.0, 11.0, 12.0], "dec": [0.0, 1.0, 2.0]})
        dvec = np.array([15.0, 16.0, 17.0])
        out = model.complete_catalog(
            catalog=df,
            columns_to_add=["dist"] + list(self.MAGS),
            dist=dvec,
            verbose=False,
        )
        assert np.allclose(out["dist"].to_numpy(), dvec)
        assert out[list(self.MAGS)].notna().all().all()

    def test_dist_vector_wrong_length_raises(self, minimal_stream_config):
        model = StreamModel(minimal_stream_config)
        df = pd.DataFrame({"ra": [10.0, 11.0, 12.0], "dec": [0.0, 1.0, 2.0]})
        with pytest.raises(ValueError):
            model.complete_catalog(
                catalog=df,
                columns_to_add=["dist"],
                dist=np.array([1.0, 2.0]),
                verbose=False,
            )

    def test_dist_overrides_distance_model(self, stream_config_with_distance):
        """When given, `dist` is used instead of the configured distance model."""
        model = StreamModel(stream_config_with_distance)  # model would give 16.8
        df = pd.DataFrame({"phi1": [0.0, 1.0], "phi2": [0.0, 0.0]})
        out = model.complete_catalog(
            catalog=df, columns_to_add=["dist"], dist=20.0, verbose=False
        )
        assert np.allclose(out["dist"].to_numpy(), 20.0)

    def test_mags_without_dist_or_phi_raises(self, minimal_stream_config):
        """No distance_modulus model, no dist, no phi1 -> cannot fill magnitudes."""
        model = StreamModel(minimal_stream_config)
        df = pd.DataFrame({"ra": [10.0, 11.0], "dec": [0.0, 1.0]})
        with pytest.raises(ValueError):
            model.complete_catalog(
                catalog=df, columns_to_add=list(self.MAGS), verbose=False
            )

    def test_complete_catalog_exposes_mass_column(self, stream_config_with_distance):
        """A completed catalog carries the shared `mass` column."""
        model = StreamModel(stream_config_with_distance)
        out = model.complete_catalog(catalog=None, size=20, verbose=False)
        assert "mass" in out.columns
        assert out["mass"].notna().all()
        assert np.issubdtype(out["mass"].dtype, np.floating)

    def test_input_mass_column_drives_magnitudes(self, stream_config_with_distance):
        """Providing a `mass` column makes the sampled mags match those masses."""
        model = StreamModel(stream_config_with_distance)
        iso = model.isochrone
        _, masses = iso.sample(15, 16.8, rng=np.random.default_rng(3))
        df = pd.DataFrame({"dist": [16.8] * 15, "mass": masses})
        columns_to_add = list(self.MAGS)
        out = model.complete_catalog(
            catalog=df, columns_to_add=columns_to_add, verbose=False
        )
        # masses is supplied directly, so rng is unused here -- omit it rather
        # than implying it affects the result.
        direct, _ = iso.sample(15, 16.8, masses=masses)
        # MAGS is a set, so list(self.MAGS) order is not guaranteed; look up
        # each column by its known name instead of relying on position.
        assert np.allclose(out["lsst_g_true"].to_numpy(), direct[("lsst_yr4", "g")])
        assert np.allclose(out["lsst_r_true"].to_numpy(), direct[("lsst_yr4", "r")])

    def test_reproductivility_complete_catalog(self, stream_config_with_distance):
        """Repeated calls to complete_catalog with the same RNG seed yield identical results."""
        model1 = StreamModel(stream_config_with_distance)
        model2 = StreamModel(stream_config_with_distance)
        # Two independently-seeded generators, not the same shared instance --
        # reusing one `rng` object across both calls would advance its state
        # between calls and never agree.
        out1 = model1.complete_catalog(
            catalog=None, size=10, verbose=False, rng=np.random.default_rng(0)
        )
        out2 = model2.complete_catalog(
            catalog=None, size=10, verbose=False, rng=np.random.default_rng(0)
        )

        for col in out1.columns:
            assert np.allclose(
                out1[col].to_numpy(), out2[col].to_numpy()
            ), f"Column {col} differs between runs"


# ---------------------------------------------------------------------------
# IsochroneModel — shared initial masses (user-supplied + `mass` column)
# ---------------------------------------------------------------------------


@pytest.fixture
def single_survey_iso_config(stream_config_with_distance):
    """Legacy flat config: single survey/release, namespace = 'lsst_yr4'."""
    return stream_config_with_distance["isochrone"]


@pytest.fixture
def multi_survey_iso_config(multisurvey_stream_config):
    """Multi-survey config: two namespaces sharing the same ugali survey/bands."""
    return multisurvey_stream_config["isochrone"]


@pytest.mark.model
class TestIsochroneModel:
    """All standalone IsochroneModel behaviour: construction, mass sampling,
    and apparent-magnitude sampling for both the single-survey (legacy) and
    multi-survey config forms.
    """

    # -- construction / config normalization ------------------------------

    def test_single_survey_sets_legacy_attrs(self, single_survey_iso_config):
        iso = IsochroneModel(single_survey_iso_config)
        assert iso.multi_survey is False
        assert iso.surveys == ["lsst_yr4"]
        assert iso.survey_name == "lsst_yr4"
        assert iso.band_1 == "g", "Legacy band_1 must be the first survey's band_1"
        assert iso.band_2 == "r", "Legacy band_2 must be the first survey's band_2"

    def test_single_survey_without_release_uses_bare_namespace(self):
        iso = IsochroneModel(
            {
                "name": "Marigo2017",
                "survey": "lsst",
                "age": 12.0,
                "z": 0.0006,
                "band_1": "g",
                "band_2": "r",
            }
        )
        assert iso.surveys == [
            "lsst"
        ], "Legacy single-survey config without release should yield bare namespace"
        assert (
            iso.survey_name == "lsst"
        ), "Legacy survey_name must be the bare namespace when no release is given"

    def test_old_ugali_raises_actionable_import_error(
        self, single_survey_iso_config, monkeypatch
    ):
        """ugali < 1.9 (no Vega->AB support) is rejected with an install hint."""
        import ugali.isochrone.model

        monkeypatch.delattr(ugali.isochrone.model.Isochrone, "vega_to_ab")
        with pytest.raises(ImportError, match="ugali >= 1.9"):
            IsochroneModel(single_survey_iso_config)

    def test_multi_survey_sets_namespace_attrs(self, multi_survey_iso_config):
        iso = IsochroneModel(multi_survey_iso_config)
        assert iso.multi_survey is True
        assert set(iso.surveys) == {"lsst_yr4", "lsst_yr5"}
        assert iso.survey_bands["lsst_yr4"] == ("g", "r")
        assert iso.survey_bands["lsst_yr5"] == ("g", "r")

    def test_multi_survey_primary_iso_is_first_entry(self, multi_survey_iso_config):
        """The first `surveys` entry drives the legacy/shared attrs (survey_name, iso, band_1/2)."""
        iso = IsochroneModel(multi_survey_iso_config)
        assert (
            iso.survey_name == iso.surveys[0]
        ), "Legacy survey_name must be the first survey's"
        assert (
            iso.iso is iso.isos[iso.survey_name]
        ), "Primary isochrone must be the first survey's"

    # -- mass sampling ------------------------------------------------------

    def test_sample_masses_returns_requested_length(
        self, single_survey_iso_config, rng
    ):
        iso = IsochroneModel(single_survey_iso_config)
        masses = iso.sample_masses(50, rng=rng)
        assert (
            len(masses) == 50
        ), "sample_masses must return the requested number of masses"

    def test_sample_masses_reproducible_with_same_rng_state(
        self, single_survey_iso_config
    ):
        iso = IsochroneModel(single_survey_iso_config)
        m1 = iso.sample_masses(30, rng=np.random.default_rng(0))
        m2 = iso.sample_masses(30, rng=np.random.default_rng(0))
        assert np.array_equal(
            m1, m2
        ), "Same RNG state should produce identical mass draws"

    def test_sample_masses_positive(self, single_survey_iso_config, rng):
        iso = IsochroneModel(single_survey_iso_config)
        masses = iso.sample_masses(100, rng=rng)
        assert np.all(masses > 0), "All sampled masses must be positive"

    # -- sample(): single-survey (legacy) form ------------------------------

    def test_sample_returns_dict_keyed_by_namespace_and_band(
        self, single_survey_iso_config, rng
    ):
        iso = IsochroneModel(single_survey_iso_config)
        mags, masses = iso.sample(20, distance_modulus=16.8, rng=rng)
        assert set(mags.keys()) == {
            ("lsst_yr4", "g"),
            ("lsst_yr4", "r"),
        }, "Returned mags dict must have keys for each namespace and band"
        assert len(masses) == 20, "Returned masses array must have the requested length"
        for arr in mags.values():
            assert arr.shape == (
                20,
            ), "Each magnitude array must have the requested length"

    def test_sample_with_explicit_masses_is_deterministic(
        self, single_survey_iso_config, rng
    ):
        """No rng needed when masses are supplied; same masses -> same mags."""
        iso = IsochroneModel(single_survey_iso_config)
        masses = iso.sample_masses(25, rng=rng)

        mags1, masses1 = iso.sample(25, distance_modulus=16.8, masses=masses)
        mags2, masses2 = iso.sample(25, distance_modulus=16.8, masses=masses)

        assert np.array_equal(masses1, masses2)
        for key in mags1:
            assert np.allclose(
                mags1[key], mags2[key]
            ), f"Magnitudes for {key} differ when using the same masses"

    def test_sample_rejects_mismatched_mass_length(self, single_survey_iso_config):
        iso = IsochroneModel(single_survey_iso_config)
        with pytest.raises(ValueError):
            iso.sample(50, distance_modulus=16.8, masses=np.ones(49))

    def test_distance_modulus_shifts_magnitudes(self, single_survey_iso_config):
        """Apparent mag at dm=17 must be exactly 1 mag fainter than at dm=16, same masses."""
        iso = IsochroneModel(single_survey_iso_config)
        masses = iso.sample_masses(20, rng=np.random.default_rng(1))

        mags_near, _ = iso.sample(20, distance_modulus=16.0, masses=masses)
        mags_far, _ = iso.sample(20, distance_modulus=17.0, masses=masses)

        for key in mags_near:
            assert np.allclose(
                mags_far[key] - mags_near[key], 1.0
            ), f"Distance modulus shift failed for {key}: expected 1.0, got {mags_far[key] - mags_near[key]}"

    def test_distance_modulus_none_returns_absolute_magnitudes(
        self, single_survey_iso_config
    ):
        """distance_modulus=None must be a no-op (apparent == absolute mag)."""
        iso = IsochroneModel(single_survey_iso_config)
        masses = iso.sample_masses(20, rng=np.random.default_rng(2))

        mags_dm0, _ = iso.sample(20, distance_modulus=0.0, masses=masses)
        mags_none, _ = iso.sample(20, distance_modulus=None, masses=masses)

        for key in mags_dm0:
            assert np.allclose(
                mags_dm0[key], mags_none[key]
            ), f"distance_modulus=None should yield same mags as dm=0 for {key}"

    # -- sample(): multi-survey form, "same physical star" invariant -------

    def test_sample_returns_all_namespace_band_keys(self, multi_survey_iso_config, rng):
        iso = IsochroneModel(multi_survey_iso_config)
        mags, masses = iso.sample(30, distance_modulus=16.8, rng=rng)
        assert set(mags.keys()) == {
            ("lsst_yr4", "g"),
            ("lsst_yr4", "r"),
            ("lsst_yr5", "g"),
            ("lsst_yr5", "r"),
        }, f"Expected keys for all namespaces and bands, got {set(mags.keys())}"

    def test_shared_mass_draw_gives_identical_mags_across_namespaces(
        self, multi_survey_iso_config, rng
    ):
        """lsst_yr4 and lsst_yr5 wrap the same ugali survey/bands, so a single
        shared mass draw must give bit-for-bit identical magnitudes in both --
        this is the core 'same physical star, consistent across surveys'
        contract the multi-survey path exists for."""
        iso = IsochroneModel(multi_survey_iso_config)
        mags, masses = iso.sample(30, distance_modulus=16.8, rng=rng)

        assert np.array_equal(
            mags[("lsst_yr4", "g")], mags[("lsst_yr5", "g")]
        ), "Shared mass draw should yield identical g mags across namespaces"
        assert np.array_equal(
            mags[("lsst_yr4", "r")], mags[("lsst_yr5", "r")]
        ), "Shared mass draw should yield identical r mags across namespaces"

    def test_explicit_masses_shared_across_namespaces(self, multi_survey_iso_config):
        iso = IsochroneModel(multi_survey_iso_config)
        masses = iso.sample_masses(15, rng=np.random.default_rng(5))
        mags, _ = iso.sample(15, distance_modulus=16.8, masses=masses)

        assert np.array_equal(
            mags[("lsst_yr4", "g")], mags[("lsst_yr5", "g")]
        ), "Explicit masses should yield identical g mags across namespaces"
        assert np.array_equal(
            mags[("lsst_yr4", "r")], mags[("lsst_yr5", "r")]
        ), "Explicit masses should yield identical r mags across namespaces"


# ---------------------------------------------------------------------------
# IsochroneModel — any number of bands, on demand, with the HB spread
# ---------------------------------------------------------------------------

# Shared stellar population of the multi-band tests.
_POP = {"name": "Marigo2017", "age": 12.0, "z": 0.0006}


def _with_isochrone(stream_config, isochrone):
    """Copy of a stream config with its isochrone section replaced."""
    return {**stream_config, "isochrone": isochrone}


def _des_euclid(des_bands=None, euclid_bands=None):
    """Multi-survey isochrone config, DES first (primary), optional bands."""
    surveys = {"des": {"survey": "des"}, "euclid": {"survey": "euclid"}}
    if des_bands is not None:
        surveys["des"]["bands"] = list(des_bands)
    if euclid_bands is not None:
        surveys["euclid"]["bands"] = list(euclid_bands)
    return {**_POP, "surveys": surveys}


@pytest.mark.model
class TestMultiBandIsochrone:
    """Bands beyond the legacy pair: config, on-demand sampling, HB spread and
    reproducibility (one IMF draw whatever the surveys and bands)."""

    # -- configuration --------------------------------------------------------

    def test_bands_list_gives_one_column_per_band(self, stream_config_with_distance):
        """`bands: [g, r, i]` -> three columns, g/r identical to `band_1/band_2`."""
        legacy = StreamModel(
            _with_isochrone(
                stream_config_with_distance,
                {**_POP, "survey": "des", "band_1": "g", "band_2": "r"},
            )
        ).sample(500, seed=4)
        listed = StreamModel(
            _with_isochrone(
                stream_config_with_distance,
                {**_POP, "survey": "des", "bands": ["g", "r", "i"]},
            )
        ).sample(500, seed=4)
        assert {"des_g_true", "des_r_true", "des_i_true"}.issubset(listed.columns)
        assert "des_i_true" not in legacy.columns
        pd.testing.assert_frame_equal(listed[legacy.columns], legacy, check_exact=True)

    def test_legacy_attrs_follow_bands(self):
        iso = IsochroneModel({**_POP, "survey": "des", "bands": ["i", "g", "z"]})
        assert iso.default_bands == {"des": ("i", "g", "z")}
        assert iso.survey_bands is iso.default_bands
        assert (iso.band_1, iso.band_2) == ("i", "g")
        single = IsochroneModel({**_POP, "survey": "euclid", "bands": ["VIS"]})
        assert (single.band_1, single.band_2) == ("VIS", "VIS")

    def test_no_bands_defaults_to_whole_filter_set(self, stream_config_with_distance):
        """No band configured -> every band of ugali's list for that survey."""
        from ugali.isochrone import Marigo2017

        model = StreamModel(_with_isochrone(stream_config_with_distance, _des_euclid()))
        des, euclid = Marigo2017.band_names["des"], Marigo2017.band_names["euclid"]
        assert model.isochrone.default_bands == {
            "des": tuple(des),
            "euclid": tuple(euclid),
        }
        out = model.sample(100, seed=1)
        expected = [f"des_{b}_true" for b in des] + [f"euclid_{b}_true" for b in euclid]
        assert set(expected).issubset(out.columns)
        assert out[expected].notna().all().all()

    def test_available_bands_cover_the_filter_set(self):
        iso = IsochroneModel(_des_euclid(des_bands=["g", "r"]))
        assert set(iso.available_bands("euclid")) >= {
            "VIS",
            "Y",
            "Blue",
            "J",
            "Red",
            "H",
        }
        assert set(iso.available_bands("des")) >= {"u", "g", "r", "i", "z", "Y"}
        assert "mass_init" not in iso.available_bands("des")
        with pytest.raises(ValueError, match="Unknown isochrone survey"):
            iso.available_bands("roman")

    def test_surveys_list_form_matches_mapping(self, stream_config_with_distance):
        """`surveys: [des, euclid]` is the mapping with each `survey` spelled out."""
        as_list = {**_POP, "surveys": ["des", "euclid"]}
        out_list = StreamModel(
            _with_isochrone(stream_config_with_distance, as_list)
        ).sample(300, seed=2)
        out_map = StreamModel(
            _with_isochrone(stream_config_with_distance, _des_euclid())
        ).sample(300, seed=2)
        pd.testing.assert_frame_equal(out_list, out_map, check_exact=True)

    def test_survey_defaults_to_key_name_not_ugali_des(self):
        """Without `survey`, the filter set is the key's survey name, release
        dropped (ugali itself would silently fall back to DES)."""
        iso = IsochroneModel(
            {**_POP, "surveys": {"lsst_yr4": {"bands": ["g", "r"]}, "roman": {}}}
        )
        assert iso.isos["lsst_yr4"].survey == "lsst"
        assert iso.isos["roman"].survey == "roman"
        explicit = IsochroneModel(
            {**_POP, "surveys": {"custom": {"survey": "des", "bands": ["g"]}}}
        )
        assert explicit.isos["custom"].survey == "des"

    def test_bands_and_legacy_pair_together_raise(self):
        with pytest.raises(ValueError, match="either `bands` or `band_1`/`band_2`"):
            IsochroneModel({**_POP, "survey": "des", "bands": ["g"], "band_1": "g"})

    @pytest.mark.parametrize("bands", [["F158"], ["VIS", "H", "F158"]])
    def test_unavailable_band_in_config_lists_available_bands(self, bands):
        """A Roman band configured for Euclid is rejected, naming the valid ones."""
        with pytest.raises(ValueError, match=r"Available bands.*VIS"):
            IsochroneModel({**_POP, "survey": "euclid", "bands": bands})

    def test_unavailable_band_in_sample_lists_available_bands(self):
        iso = IsochroneModel(_des_euclid(euclid_bands=["VIS"]))
        with pytest.raises(ValueError, match=r"Available bands.*VIS"):
            iso.sample(
                10, 16.0, rng=np.random.default_rng(0), bands={"euclid": ["F158"]}
            )

    # -- sampling on demand ---------------------------------------------------

    def test_sample_any_available_band_on_request(self):
        iso = IsochroneModel(_des_euclid(des_bands=["g", "r"], euclid_bands=["VIS"]))
        mags, _ = iso.sample(
            20,
            16.0,
            rng=np.random.default_rng(0),
            bands={"des": ["z"], "euclid": ["J"]},
        )
        assert set(mags) == {("des", "z"), ("euclid", "J")}

    def test_complete_catalog_accepts_non_default_band(
        self, stream_config_with_distance
    ):
        model = StreamModel(
            _with_isochrone(
                stream_config_with_distance,
                {**_POP, "survey": "des", "release": "yr6", "bands": ["g", "r"]},
            )
        )
        out = model.complete_catalog(
            catalog=None,
            size=30,
            columns_to_add=["phi1", "dist", "des_Y_true"],
            verbose=False,
        )
        assert out["des_Y_true"].notna().all()
        assert "des_g_true" not in out.columns, "only the requested band is added"

    # -- reproducibility (R1-R3) ----------------------------------------------

    def test_one_imf_draw_whatever_the_bands(self):
        """R1: sampling draws the legacy masses (rng.choice over the mass grid)
        and leaves the generator in the same state, for any surveys/bands."""
        iso = IsochroneModel(_des_euclid(des_bands=["g", "r", "i", "z"]))
        grid = iso.iso.sample(mass_min=iso._MASS_MIN, mass_steps=iso._MASS_STEPS)
        pdf = grid[1] / grid[1].sum()
        rng, twin = np.random.default_rng(11), np.random.default_rng(11)
        _, masses = iso.sample(1000, 16.0, rng=rng)
        assert np.array_equal(masses, twin.choice(grid[0], size=1000, p=pdf))
        assert rng.bit_generator.state == twin.bit_generator.state
        assert np.array_equal(
            iso.sample_masses(1000, rng=np.random.default_rng(11)), masses
        )

    @pytest.mark.parametrize(
        "isochrone",
        [
            _des_euclid(des_bands=["g", "r", "i", "z"]),  # more bands
            _des_euclid(des_bands=["g", "r"], euclid_bands=["VIS", "H"]),  # + survey
        ],
        ids=["more_bands", "more_surveys"],
    )
    def test_more_bands_or_surveys_leave_other_columns_unchanged(
        self, stream_config_with_distance, isochrone
    ):
        """R2: same seed -> the shared columns are bit-identical."""
        ref = StreamModel(
            _with_isochrone(
                stream_config_with_distance,
                {**_POP, "surveys": {"des": {"survey": "des", "bands": ["g", "r"]}}},
            )
        ).sample(3000, seed=3)
        out = StreamModel(
            _with_isochrone(stream_config_with_distance, isochrone)
        ).sample(3000, seed=3)
        pd.testing.assert_frame_equal(out[ref.columns], ref, check_exact=True)

    def test_completing_bands_later_matches_sampling_them_at_once(
        self, stream_config_with_distance
    ):
        """R3: complete `i` in a second call (reusing `mass`) -> identical to
        sampling g, r, i together."""
        n = 20000  # enough stars for a few on the horizontal branch

        def model(bands):
            return StreamModel(
                _with_isochrone(
                    stream_config_with_distance,
                    {**_POP, "survey": "des", "bands": bands},
                )
            )

        first = model(["g", "r"]).complete_catalog(
            catalog=None, size=n, rng=np.random.default_rng(5), verbose=False
        )
        hb_offset = model(["g", "r"]).isochrone._hb_offset(first["mass"])
        assert np.count_nonzero(hb_offset) > 0, "sample has no dispersed HB star"
        later = model(["g", "r"]).complete_catalog(
            catalog=first, columns_to_add=["des_i_true"], verbose=False
        )
        at_once = model(["g", "r", "i"]).complete_catalog(
            catalog=None, size=n, rng=np.random.default_rng(5), verbose=False
        )
        pd.testing.assert_frame_equal(later[at_once.columns], at_once, check_exact=True)

    # -- horizontal-branch spread ---------------------------------------------

    def test_magnitudes_depend_on_mass_alone(self):
        """The drawn masses, supplied back, give the same magnitudes, HB stars
        included: the masses alone reproduce a sample."""
        iso = IsochroneModel(_des_euclid(des_bands=["g", "r"], euclid_bands=["VIS"]))
        mags, masses = iso.sample(50000, 16.0, rng=np.random.default_rng(1))
        assert np.count_nonzero(iso._hb_offset(masses)) > 0, "no dispersed HB star"
        again, _ = iso.sample(50000, 16.0, masses=masses)
        for key in mags:
            assert np.array_equal(again[key], mags[key]), key

    def test_hb_spread_is_shared_by_all_bands_and_surveys(self):
        """HB stars are offset by one of ugali's dispersion values, the same in
        every band of every survey, so their colours are the isochrone's."""
        iso = IsochroneModel(
            _des_euclid(des_bands=["g", "r"], euclid_bands=["VIS", "H"])
        )
        mags, masses = iso.sample(50000, 16.0, rng=np.random.default_rng(0))
        offset = None
        for name in ("des", "euclid"):
            undispersed = IsochroneModel._absolute_mags(
                iso.isos[name], masses, iso.default_bands[name]
            )
            for band, mag in undispersed.items():
                this = mags[(name, band)] - 16.0 - mag
                offset = this if offset is None else offset
                assert np.allclose(this, offset, rtol=0, atol=1e-9), (name, band)
        offset = np.round(offset, 6)
        # ugali disperses the HB by up to +/- hb_spread in 0.025 steps.
        steps = offset / 0.025
        assert np.allclose(steps, np.round(steps), rtol=0, atol=1e-6)
        assert np.all(np.abs(offset) <= iso.iso.hb_spread + 1e-9)
        assert np.count_nonzero(offset) > 0
        # Only horizontal-branch masses are dispersed.
        hb_masses = iso.iso.mass_init[iso.iso.stage == iso.iso.hb_stage]
        on_hb = (masses > hb_masses.min()) & (masses < hb_masses.max())
        assert np.all(on_hb[offset != 0])

    def test_hb_offsets_cover_ugali_dispersion_evenly(self):
        """Every one of ugali's HB offsets is used by as many HB masses."""
        iso = IsochroneModel({**_POP, "survey": "des", "bands": ["g"]})
        values, counts = np.unique(iso._hb_mass_offsets, return_counts=True)
        assert np.allclose(values, np.linspace(-0.1, 0.1, 9))
        assert np.all(np.abs(counts / counts.mean() - 1) < 0.02)
