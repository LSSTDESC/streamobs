"""
tests/test_euclid.py
====================
Implementation behaviour tests for the Euclid releases (``euclid/q1``,
``euclid/dr1``).

These complement the generic contract checks every registered survey gets in
``tests/test_surveys.py`` with the Euclid-specific behaviour:

- the ``euclid/q1`` Survey loads with VIS as the reference band and Y/J/H as
  forced-photometry bands, with the extinction and saturation values the
  config documents;
- the ugali ``euclid`` isochrone set is already in AB, so the Roman Vega->AB
  offset must be a no-op for Euclid bands and the isochrone must produce AB
  magnitudes in the band names the survey uses (``VIS``, ``Y``, ``J``, ``H``);
- ``inject()`` produces ``euclid_q1``-namespaced columns, including a
  forced-photometry band that resolves to the ``_nocut`` error curves;
- ``euclid/dr1`` is uniform at the Q1-anchored depth and shares Q1's tables.

Tests that need the data files under ``data/surveys/euclid_q1/`` skip when
they are absent.
"""

import os

import numpy as np
import pytest

_DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "surveys")


def _has_products(release):
    d = os.path.join(_DATA, f"euclid_{release}")
    return os.path.isdir(d) and any(
        f.endswith((".fits.gz", ".csv")) for f in os.listdir(d)
    )


_skip_no_q1 = pytest.mark.skipif(
    not _has_products("q1"), reason="euclid_q1 products not present under data/surveys/"
)
_skip_no_dr1 = pytest.mark.skipif(
    not _has_products("dr1"), reason="euclid_dr1 products not present under data/surveys/"
)


@pytest.fixture(scope="module")
def euclid_q1():
    if not _has_products("q1"):
        pytest.skip("euclid_q1 products not present")
    from streamobs import surveys

    return surveys.Survey.load("euclid", release="q1", verbose=False)


@pytest.fixture(scope="module")
def euclid_dr1():
    if not _has_products("dr1"):
        pytest.skip("euclid_dr1 products not present")
    from streamobs import surveys

    return surveys.Survey.load("euclid", release="dr1", verbose=False)


@pytest.fixture(scope="module")
def euclid_q1_injector(euclid_q1):
    from streamobs.observed import StreamInjector

    return StreamInjector(survey=euclid_q1, verbose=False)


def _valid(m):
    """Mask of on-footprint pixels: hp.read_map leaves hp.UNSEEN in place."""
    import healpy as hp

    return np.isfinite(m) & (m != hp.UNSEEN) & (m > 0)


def _covered_positions(survey, band, n, rng):
    """ra/dec of n random covered pixels of ``band``'s maglim map."""
    import healpy as hp

    m = survey.maglim_maps[band]
    pix = np.where(np.isfinite(m) & (m > survey.saturation[band]))[0]
    pick = rng.choice(pix, size=n, replace=True)
    ra, dec = hp.pix2ang(hp.npix2nside(m.size), pick, lonlat=True)
    return ra, dec


# ---------------------------------------------------------------------------
# Survey loading
# ---------------------------------------------------------------------------
@pytest.mark.surveys
class TestEuclidQ1Survey:
    @_skip_no_q1
    def test_namespace_and_bands(self, euclid_q1):
        assert euclid_q1.namespace == "euclid_q1"
        assert set(euclid_q1.bands) == {"VIS", "Y", "J", "H"}
        assert euclid_q1.completeness_band == "VIS"

    @_skip_no_q1
    def test_extinction_coefficients(self, euclid_q1):
        # F99 (R_V = 3.1) at the SVO pivot wavelengths x 0.86 (SF11), see config
        expected = {"VIS": 1.790, "Y": 0.879, "J": 0.595, "H": 0.400}
        for band, coeff in expected.items():
            assert np.isclose(euclid_q1.coeff_extinc[band], coeff, atol=1e-3)
        # ordering with wavelength is the sanity check that survives any
        # revision of the absolute values
        assert (
            euclid_q1.coeff_extinc["VIS"]
            > euclid_q1.coeff_extinc["Y"]
            > euclid_q1.coeff_extinc["J"]
            > euclid_q1.coeff_extinc["H"]
        )

    @_skip_no_q1
    def test_saturation(self, euclid_q1):
        assert euclid_q1.saturation["VIS"] == 17.5
        for b in ("Y", "J", "H"):
            assert euclid_q1.saturation[b] == 15.5

    @_skip_no_q1
    def test_anchored_depths(self, euclid_q1):
        """Map medians sit on the Wide-survey 5-sigma point-source depths."""
        for band, m5 in (("VIS", 26.2), ("Y", 24.5), ("J", 24.5), ("H", 24.5)):
            m = euclid_q1.maglim_maps[band]
            assert np.isclose(np.median(m[_valid(m)]), m5, atol=0.05)

    @_skip_no_q1
    def test_forced_bands_use_nocut_curves(self, euclid_q1):
        for band in ("Y", "J", "H"):
            for kind in ("catalog", "sample"):
                fn = euclid_q1._resolve_log_photo_error(kind, band=band)
                assert fn is getattr(euclid_q1, f"log_photo_error_{kind}_nocut")

    @_skip_no_q1
    def test_sample_curve_is_inflated_catalog_curve(self, euclid_q1):
        """sample = catalog x f with f > 1: the anchor says the reported errors
        are optimistic at the depth, so the noise draw must exceed `magerr`."""
        grid = np.arange(-6.0, 0.0, 0.25)
        ratio = 10 ** (
            euclid_q1.log_photo_error_sample(grid) - euclid_q1.log_photo_error_catalog(grid)
        )
        assert np.all(ratio > 1.5) and np.all(ratio < 2.5)
        assert np.allclose(ratio, ratio[0], rtol=1e-3), "constant factor by construction"

    @_skip_no_q1
    def test_classification_collapses_faintward(self, euclid_q1):
        """POINT_LIKE_PROB loses the stars near the depth; detection does not."""
        maglim = 26.2
        det = euclid_q1.get_detection_efficiency("VIS", np.array([22.0, 25.25]), maglim)
        cls = euclid_q1.get_classification_efficiency("VIS", np.array([22.0, 25.25]), maglim)
        assert det[0] > 0.9 and det[1] > 0.85
        assert cls[0] > 0.9 and cls[1] < 0.3


# ---------------------------------------------------------------------------
# Isochrones: the ugali `euclid` set is AB and uses the survey's band names
# ---------------------------------------------------------------------------
@pytest.mark.model
class TestEuclidIsochrone:
    def test_vega_to_ab_is_noop_for_euclid_bands(self):
        from streamobs.model import ROMAN_VEGA_TO_AB, IsochroneModel

        for band in ("VIS", "Y", "J", "H"):
            assert band not in ROMAN_VEGA_TO_AB
            assert IsochroneModel._to_ab(None, band, 20.0) == 20.0

    def test_euclid_isochrone_samples_ab_magnitudes(self):
        """A 12 Gyr metal-poor isochrone in VIS/H gives finite AB magnitudes
        with the colour range of an old main sequence + RGB."""
        from streamobs.model import IsochroneModel

        try:
            iso = IsochroneModel(
                {
                    "name": "Marigo2017",
                    "survey": "euclid",
                    "release": "q1",
                    "age": 12.0,
                    "z": 0.0006,
                    "band_1": "VIS",
                    "band_2": "H",
                }
            )
        except Exception as exc:  # pragma: no cover - depends on ~/.ugali
            pytest.skip(f"ugali euclid isochrones unavailable: {exc}")
        assert iso.surveys == ["euclid_q1"]
        assert iso.survey_bands["euclid_q1"] == ("VIS", "H")
        # sample() returns ({(namespace, band): apparent_mag}, masses)
        mags, masses = iso.sample(500, 16.8, rng=np.random.default_rng(1))
        assert len(masses) == 500
        vis = np.asarray(mags[("euclid_q1", "VIS")], dtype=float)
        h = np.asarray(mags[("euclid_q1", "H")], dtype=float)
        assert np.all(np.isfinite(vis)) and np.all(np.isfinite(h))
        colour = vis - h
        assert -0.5 < np.median(colour) < 2.0, "VIS - H out of range for an old population"
        assert np.min(vis) > 16.8 - 4.0 and np.max(vis) < 16.8 + 14.0


# ---------------------------------------------------------------------------
# Injection — namespaced columns, forced band
# ---------------------------------------------------------------------------
@pytest.mark.surveys
class TestEuclidQ1Injection:
    @_skip_no_q1
    def test_inject_produces_namespaced_columns(self, euclid_q1, euclid_q1_injector):
        import pandas as pd

        from streamobs.columns import err_col, flag_col, obs_col, true_col

        rng = np.random.default_rng(3)
        ra, dec = _covered_positions(euclid_q1, "VIS", 40, rng)
        df = pd.DataFrame(
            {
                "ra": ra,
                "dec": dec,
                true_col("VIS", "euclid_q1"): rng.uniform(19.0, 25.0, 40),
                true_col("H", "euclid_q1"): rng.uniform(18.0, 24.0, 40),
            }
        )
        out = euclid_q1_injector.inject(df, bands=["VIS", "H"], verbose=False)
        for col in (
            true_col("VIS", "euclid_q1"),  # euclid_VIS_true  (name-only)
            obs_col("VIS", "euclid_q1"),  # euclid_q1_VIS_obs
            err_col("VIS", "euclid_q1"),
            obs_col("H", "euclid_q1"),  # forced band: resolves to the _nocut curves
            err_col("H", "euclid_q1"),
            flag_col("euclid_q1"),
        ):
            assert col in out.columns, f"missing {col} after inject()"
        assert np.allclose(out[true_col("VIS", "euclid_q1")].values, df[true_col("VIS", "euclid_q1")].values)
        errs = out[err_col("H", "euclid_q1")].to_numpy(dtype=float)
        assert np.any(np.isfinite(errs)) and np.all(errs[np.isfinite(errs)] > 0)


# ---------------------------------------------------------------------------
# DR1
# ---------------------------------------------------------------------------
@pytest.mark.surveys
class TestEuclidDR1Survey:
    @_skip_no_dr1
    def test_uniform_depth_at_anchor(self, euclid_dr1):
        for band, m5 in (("VIS", 26.2), ("Y", 24.5), ("J", 24.5), ("H", 24.5)):
            m = euclid_dr1.maglim_maps[band]
            v = m[_valid(m)]
            assert v.size > 0
            assert np.allclose(v, m5, atol=1e-3), f"{band} map is not uniform at {m5}"

    @_skip_no_dr1
    def test_footprint_area(self, euclid_dr1):
        import healpy as hp

        m = euclid_dr1.maglim_maps["VIS"]
        area = _valid(m).sum() * hp.nside2pixarea(hp.npix2nside(m.size), degrees=True)
        assert 2000 < area < 2300, f"DR1 footprint {area:.0f} deg^2, expected ~2100"

    @_skip_no_dr1
    @_skip_no_q1
    def test_shares_q1_tables(self, euclid_q1, euclid_dr1):
        grid = np.arange(-6.0, 0.5, 0.5)
        assert np.allclose(euclid_q1.completeness(grid), euclid_dr1.completeness(grid), equal_nan=True)
        assert np.allclose(
            euclid_q1.log_photo_error_catalog(grid), euclid_dr1.log_photo_error_catalog(grid), equal_nan=True
        )
