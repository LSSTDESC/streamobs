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
            "lsst_g_true",
            "lsst_r_true",
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
        data_without_mag = stream_catalog.drop(columns=["lsst_g_true", "lsst_r_true"])
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

    def test_mag_sampling_reproducibility(self, mock_injector, seed):
        """Test that magnitude sampling is reproducible with the same random seed.
        Created for debug, but kept since it is a useful test of the sampling function.
        """
        apparent_mag = np.array([20.0, 21.0, 22.0])
        mag_err = np.array([0.1, 0.2, 0.3])
        sample_with_rng1 = mock_injector.sample_measured_magnitudes(
            apparent_mag,
            mag_err,
            rng=np.random.default_rng(seed),
        )
        sample_with_rng1 = pd.to_numeric(sample_with_rng1, errors="coerce")

        sample_with_seed1 = mock_injector.sample_measured_magnitudes(
            apparent_mag,
            mag_err,
            seed=seed,
        )
        sample_with_seed1 = pd.to_numeric(sample_with_seed1, errors="coerce")

        sample_with_rng2 = mock_injector.sample_measured_magnitudes(
            apparent_mag,
            mag_err,
            rng=np.random.default_rng(seed),
        )
        sample_with_rng2 = pd.to_numeric(sample_with_rng2, errors="coerce")

        sample_with_seed2 = mock_injector.sample_measured_magnitudes(
            apparent_mag,
            mag_err,
            seed=seed,
        )
        sample_with_seed2 = pd.to_numeric(sample_with_seed2, errors="coerce")

        assert np.allclose(
            sample_with_rng1, sample_with_seed1, equal_nan=True
        ), "Samples with rng and seed should be the same"
        assert np.allclose(
            sample_with_rng1, sample_with_rng2, equal_nan=True
        ), "Samples with the same rng should be the same"
        assert np.allclose(
            sample_with_seed1, sample_with_seed2, equal_nan=True
        ), "Samples with the same seed should be the same"

    def test_injection_reproducibility(
        self, mock_injector, stream_catalog, seed, verbose
    ):
        """Test that injection is reproducible with the same random seed."""
        injected_catalog_1 = mock_injector.inject(
            stream_catalog, perfect_galstarsep=True, seed=seed, verbose=verbose
        )
        injected_catalog_2 = mock_injector.inject(
            stream_catalog, perfect_galstarsep=True, seed=seed, verbose=verbose
        )

        rng = np.random.default_rng(seed)
        injected_catalog_rng1 = mock_injector.inject(
            stream_catalog, perfect_galstarsep=True, rng=rng, verbose=verbose
        )
        rng = np.random.default_rng(seed)
        injected_catalog_rng2 = mock_injector.inject(
            stream_catalog, perfect_galstarsep=True, rng=rng, verbose=verbose
        )

        for cols in injected_catalog_1.columns:
            assert (
                cols in injected_catalog_2.columns
            ), f"Column {cols} missing in second injection"
            assert (
                cols in injected_catalog_rng1.columns
            ), f"Column {cols} missing in rng1 injection"
            assert (
                cols in injected_catalog_rng2.columns
            ), f"Column {cols} missing in rng2 injection"
            col1 = pd.to_numeric(injected_catalog_1[cols], errors="coerce")
            col2 = pd.to_numeric(injected_catalog_2[cols], errors="coerce")
            col_rng1 = pd.to_numeric(injected_catalog_rng1[cols], errors="coerce")
            col_rng2 = pd.to_numeric(injected_catalog_rng2[cols], errors="coerce")
            if (col1 is not None and col2 is not None) and cols not in [
                "mu1",
                "mu2",
                "rv",
            ]:  # velocities not implemented
                assert np.allclose(
                    col1, col2, equal_nan=True
                ), f"Column {cols} values differ between injections with the same seed"
            if (col1 is not None and col_rng1 is not None) and cols not in [
                "mu1",
                "mu2",
                "rv",
            ]:  # velocities not implemented
                assert np.allclose(
                    col1, col_rng1, equal_nan=True
                ), f"Column {cols} values differ between injections with seed and rng with the same seed"
            if (col_rng1 is not None and col_rng2 is not None) and cols not in [
                "mu1",
                "mu2",
                "rv",
            ]:  # velocities not implemented
                assert np.allclose(
                    col_rng1, col_rng2, equal_nan=True
                ), f"Column {cols} values differ between injections with the same rng"


# ---------------------------------------------------------------------------
# complete_data + bands/survey API
# ---------------------------------------------------------------------------


@pytest.mark.observed
class TestCompleteDataAndAPI:
    """complete_data (single + multi survey) and the unified bands/survey API."""

    def test_complete_data_single_survey(
        self, mock_injector, stream_config_with_distance
    ):
        """complete_data fills ra/dec + true mags and preserves existing columns."""
        df = pd.DataFrame({"phi1": [-3.0, 0.0, 3.0], "phi2": [0.0, 0.0, 0.0]})
        out = mock_injector.complete_data(
            df, bands=["g", "r"], stream_config=stream_config_with_distance, seed=1
        )
        for col in ["ra", "dec", "lsst_g_true", "lsst_r_true"]:
            assert col in out.columns, f"missing {col}"
        assert out[["lsst_g_true", "lsst_r_true"]].notna().all().all()
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
            "lsst_g_true",
            "lsst_r_true",
        ]:
            assert col in out.columns, f"missing {col}"

        # Verify that we have a single column for each true magnitude, not one per survey.
        assert out[["lsst_g_true", "lsst_r_true"]].notna().all().all()

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


# ---------------------------------------------------------------------------
# Selection-draw coupling (classification_coupling) + perfect-separation nesting
# ---------------------------------------------------------------------------


@pytest.mark.observed
class TestSelectionCoupling:
    """The per-object selection uniform: nesting within a survey, coupling across.

    Two independent properties live here.

    **Nesting (within one survey).** ``get_detection_efficiency`` is >=
    ``get_completeness`` pointwise — perfect star/galaxy separation can only help
    — so the realistic selection MUST be a subset of the perfect-separation one.
    That holds only if both ``detect_flag`` calls share one per-object uniform;
    with two independent draws an object can fail the easier cut while passing
    the harder one, which is unphysical.

    **Coupling (across surveys).** ``classification_coupling`` sets
    ``P(selected in both)``: ``p_1 * p_2`` when independent (the default),
    ``min(p_1, p_2)`` when shared.

    The coupling tests deliberately use a **homogeneous** catalog — every object
    at the same sky position and the same true magnitude — so that ``p_k`` is a
    single number per survey and ``p_1 * p_2`` is the correct independent
    prediction. On a realistic catalog it is NOT: per-object probabilities vary
    and are correlated across surveys through the shared true magnitude, so the
    joint rate exceeds the product of the marginals even under independent draws
    (the same heterogeneity effect that
    ``scripts/roman/measure_joint_misclassification.py`` corrects for with 2-D
    cells). Only the reference band is injected, so no non-reference-band S/N cut
    can add per-survey failures on top of the selection draw.
    """

    NS = ("lsst_yr4", "lsst_yr5")
    BANDS = {"lsst_yr4": ["r"], "lsst_yr5": ["r"]}
    N = 20000

    @staticmethod
    def _homogeneous_catalog(injector, n=N, dm_offset=-1.0):
        """n identical objects on one pixel covered by both surveys.

        ``dm_offset`` is measured from the shallower survey's r-band limit
        (extinction folded in), so the selection probability sits near 0.5 and
        the test is sensitive to a mis-specified coupling.
        """
        import healpy as hp

        s4, s5 = injector.surveys["lsst_yr4"], injector.surveys["lsst_yr5"]
        m4, m5 = s4.get_maglim("r"), s5.get_maglim("r")
        nside = hp.get_nside(m4)
        ok = (
            np.isfinite(m4)
            & (m4 > 0)
            & (m4 != hp.UNSEEN)
            & np.isfinite(m5)
            & (m5 > 0)
            & (m5 != hp.UNSEEN)
        )
        good = np.flatnonzero(ok)
        pix = int(good[good.size // 2])
        ra, dec = hp.pix2ang(nside, pix, lonlat=True)
        ebv_pix = hp.ang2pix(hp.get_nside(s4.ebv_map), ra, dec, lonlat=True)
        extinction = float(s4.get_extinction("r", pixel=ebv_pix))
        mag_true = float(m4[pix]) - extinction + dm_offset
        return pd.DataFrame(
            {
                "ra": np.full(n, ra),
                "dec": np.full(n, dec),
                "lsst_r_true": np.full(n, mag_true),
            }
        )

    def _rates(self, injector, catalog, **kw):
        out = injector.inject(
            catalog.copy(), bands=self.BANDS, seed=11, verbose=False, **kw
        )
        assert not [
            c for c in out.columns if c.startswith("_u_select")
        ], "the private coupling column must be dropped before inject() returns"
        f4 = out[f"{self.NS[0]}_flag_observed"].to_numpy().astype(bool)
        f5 = out[f"{self.NS[1]}_flag_observed"].to_numpy().astype(bool)
        p4, p5 = float(f4.mean()), float(f5.mean())
        # guard: an all-or-nothing marginal would make every assertion vacuous
        assert 0.1 < p4 < 0.9 and 0.1 < p5 < 0.9, (
            f"marginals p4={p4:.3f} p5={p5:.3f} are not intermediate; the survey "
            "products moved and _homogeneous_catalog's dm_offset needs retuning"
        )
        return p4, p5, float((f4 & f5).mean())

    @staticmethod
    def _mc_tol(p, n, nsigma=5.0):
        return nsigma * np.sqrt(max(p * (1.0 - p), 1e-12) / n)

    def test_perfect_galstarsep_flag_is_superset(
        self, mock_injector, stream_catalog, verbose
    ):
        """flag_observed must nest inside flag_perfect_galstarsep, with no exceptions."""
        out = mock_injector.inject(
            stream_catalog.copy(), perfect_galstarsep=True, verbose=verbose, seed=3
        )
        realistic = out["lsst_yr4_flag_observed"].to_numpy().astype(bool)
        perfect = out["lsst_yr4_flag_perfect_galstarsep"].to_numpy().astype(bool)
        n_violating = int((realistic & ~perfect).sum())
        assert n_violating == 0, (
            f"{n_violating} objects pass the realistic selection but fail the "
            "perfect star/galaxy-separation one; both detect_flag calls must "
            "share a single per-object uniform."
        )
        # not vacuous: perfect separation has to gain something
        assert perfect.sum() > realistic.sum()

    def test_default_coupling_is_independent(self, mock_multisurvey_injector):
        """The default must stay `independent`: joint == p4 * p5 to MC error."""
        cat = self._homogeneous_catalog(mock_multisurvey_injector)
        p4, p5, joint = self._rates(mock_multisurvey_injector, cat)
        assert abs(joint - p4 * p5) < self._mc_tol(p4 * p5, self.N), (
            f"joint={joint:.4f} is not the independent product {p4 * p5:.4f} "
            f"(p4={p4:.4f}, p5={p5:.4f})"
        )
        assert joint < min(p4, p5) - 0.05, "independent must sit well below min()"

    def test_shared_coupling_gives_min_of_marginals(self, mock_multisurvey_injector):
        """One shared uniform collapses `selected in both` onto min(p1, p2) exactly."""
        cat = self._homogeneous_catalog(mock_multisurvey_injector)
        p4, p5, joint = self._rates(
            mock_multisurvey_injector, cat, classification_coupling="shared"
        )
        # exact, not statistical: identical objects with p4 <= p5 make survey 4's
        # selection a strict subset of survey 5's under a shared u.
        assert joint == pytest.approx(min(p4, p5), abs=1e-12), (
            f"joint={joint:.6f} != min(p4, p5)={min(p4, p5):.6f} under shared coupling"
        )

    def test_rho_zero_reproduces_independent(self, mock_multisurvey_injector):
        """The Gaussian copula at rho=0 reproduces the independent joint rate."""
        cat = self._homogeneous_catalog(mock_multisurvey_injector)
        p4, p5, joint = self._rates(
            mock_multisurvey_injector, cat, classification_coupling=0.0
        )
        assert abs(joint - p4 * p5) < self._mc_tol(p4 * p5, self.N)

    def test_rho_one_reproduces_shared(self, mock_multisurvey_injector):
        """rho=1 is the shared (comonotonic) limit."""
        cat = self._homogeneous_catalog(mock_multisurvey_injector)
        p4, p5, joint = self._rates(
            mock_multisurvey_injector, cat, classification_coupling=1.0
        )
        assert joint == pytest.approx(min(p4, p5), abs=1e-12)

    def test_intermediate_rho_interpolates(self, mock_multisurvey_injector):
        """0 < rho < 1 must sit strictly between the independent and shared rates."""
        cat = self._homogeneous_catalog(mock_multisurvey_injector)
        _, _, j_ind = self._rates(
            mock_multisurvey_injector, cat, classification_coupling="independent"
        )
        p4, p5, j_rho = self._rates(
            mock_multisurvey_injector, cat, classification_coupling=0.6
        )
        _, _, j_sh = self._rates(
            mock_multisurvey_injector, cat, classification_coupling="shared"
        )
        tol = self._mc_tol(p4 * p5, self.N)
        assert j_ind - tol < j_rho < j_sh + tol, (
            f"rho=0.6 joint={j_rho:.4f} is not between independent={j_ind:.4f} "
            f"and shared={j_sh:.4f}"
        )

    def test_coupling_dict_resolves_by_source_type(self, mock_multisurvey_injector):
        """A {source_type: value} dict picks the entry for this call's source_type."""
        cat = self._homogeneous_catalog(mock_multisurvey_injector)
        p4, p5, joint = self._rates(
            mock_multisurvey_injector,
            cat,
            source_type="stars",
            classification_coupling={"stars": "shared", "galaxies": "independent"},
        )
        assert joint == pytest.approx(min(p4, p5), abs=1e-12)

    @pytest.mark.parametrize("bad", ["comonotonic", 1.5, -2.0, {"galaxies": "shared"}])
    def test_bad_coupling_rejected(self, mock_multisurvey_injector, bad):
        """Unknown strings, out-of-range rho, and a dict missing source_type raise."""
        cat = self._homogeneous_catalog(mock_multisurvey_injector, n=10)
        with pytest.raises(ValueError):
            mock_multisurvey_injector.inject(
                cat,
                bands=self.BANDS,
                seed=11,
                verbose=False,
                classification_coupling=bad,
            )

    def test_u_select_length_is_checked(self, mock_injector, mock_survey):
        """A mis-sized u_select is a programming error, not silently broadcast."""
        with pytest.raises(ValueError):
            mock_injector.detect_flag(
                0,
                survey=mock_survey,
                mag=np.array([22.0, 23.0]),
                band="r",
                u_select=np.array([0.5]),
            )

    def test_u_select_is_used_verbatim(self, mock_injector, mock_survey):
        """u_select fully determines the outcome: u=0 always selects, u=1 never does."""
        mag = np.full(20, 22.0)
        assert mock_injector.detect_flag(
            0, survey=mock_survey, mag=mag, band="r", u_select=np.zeros(20)
        ).all()
        assert not mock_injector.detect_flag(
            0, survey=mock_survey, mag=mag, band="r", u_select=np.ones(20)
        ).any()
