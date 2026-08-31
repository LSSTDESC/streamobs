"""Unit tests for the DES/DECam Balrog reducer helpers.

The reducer is a standalone script (it must run on a machine with no streamobs
install), so it is loaded by path rather than imported as a package module.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "des" / "balrog_selection_function.py"


@pytest.fixture(scope="module")
def bsf():
    if not SCRIPT.exists():  # pragma: no cover
        pytest.skip(f"{SCRIPT} not present")
    spec = importlib.util.spec_from_file_location("balrog_selection_function", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _conf(mag, a, b, n=10_000):
    return pd.DataFrame({"mag_g": mag, "a": a, "b": b,
                         "n_pos": [n] * len(mag), "n_neg": [n] * len(mag)})


class TestDeconvolution:
    def test_perfect_surrogate_is_a_no_op(self, bsf):
        """a=1, b=0 means S IS the classifier, so nothing should change."""
        mag = np.arange(18.0, 25.0, 0.25)
        eff = np.linspace(0.95, 0.30, mag.size)
        conf = _conf(mag, np.ones(mag.size), np.zeros(mag.size))
        out, info = bsf.deconvolve_classification_eff(eff, conf, mag)
        np.testing.assert_allclose(out, eff, atol=1e-12)
        assert info["max_abs_change"] == pytest.approx(0.0, abs=1e-12)

    def test_recovers_a_known_efficiency(self, bsf):
        """Round-trip: build eff_S from a known eff_X, invert, get eff_X back."""
        mag = np.arange(18.0, 25.0, 0.25)
        eff_x = np.linspace(0.98, 0.40, mag.size)
        a = np.linspace(0.99, 0.80, mag.size)
        b = np.linspace(0.001, 0.02, mag.size)
        eff_s = a * eff_x + b * (1 - eff_x)
        out, info = bsf.deconvolve_classification_eff(eff_s, _conf(mag, a, b), mag)
        np.testing.assert_allclose(out, eff_x, atol=1e-10)
        assert info["n_bins_corrected"] == mag.size

    def test_undoes_the_faint_end_underestimate(self, bsf):
        """A surrogate that loses recall reads low; the inversion corrects up."""
        # a/b as actually measured around g ~ 24.9, held flat across a few bins
        # so there is something to interpolate between.
        mag = np.arange(24.375, 25.5, 0.25)
        a = np.full(mag.size, 0.829)
        b = np.full(mag.size, 0.010)
        eff_x = np.full(mag.size, 0.75)
        eff_s = a * eff_x + b * (1 - eff_x)
        assert np.all(eff_s < eff_x)  # the bias this exists to remove
        assert eff_s[0] == pytest.approx(0.6243, abs=1e-4)  # ~17% low, uncorrected
        out, _ = bsf.deconvolve_classification_eff(eff_s, _conf(mag, a, b), mag)
        np.testing.assert_allclose(out, eff_x, atol=1e-10)

    def test_poorly_separated_bins_are_left_alone(self, bsf):
        """As a -> b the inversion blows up, so those bins must pass through."""
        mag = np.arange(18.0, 22.0, 0.25)
        eff = np.full(mag.size, 0.6)
        conf = _conf(mag, np.full(mag.size, 0.30), np.full(mag.size, 0.25))
        out, info = bsf.deconvolve_classification_eff(eff, conf, mag)
        np.testing.assert_allclose(out, eff)
        assert info["n_bins_corrected"] == 0

    def test_result_stays_a_probability(self, bsf):
        """Noise must not push the corrected efficiency outside [0, 1]."""
        mag = np.arange(18.0, 25.0, 0.25)
        rng = np.random.default_rng(0)
        eff = np.clip(rng.normal(0.9, 0.3, mag.size), 0, 1)
        conf = _conf(mag, np.full(mag.size, 0.85), np.full(mag.size, 0.02))
        out, _ = bsf.deconvolve_classification_eff(eff, conf, mag)
        assert np.all((out >= 0.0) & (out <= 1.0))

    def test_low_count_bins_are_not_used(self, bsf):
        """Confusion bins with too few objects must not drive a correction."""
        mag = np.arange(18.0, 22.0, 0.25)
        eff = np.full(mag.size, 0.6)
        conf = _conf(mag, np.full(mag.size, 0.9), np.zeros(mag.size), n=5)
        out, info = bsf.deconvolve_classification_eff(eff, conf, mag)
        np.testing.assert_allclose(out, eff)
        assert info["n_bins_corrected"] == 0


class TestClassifierVendoring:
    def test_bdf_classifier_separates_stars_from_galaxies(self, bsf):
        """Point-like (small T) at high S/N -> star classes; large T -> galaxy."""
        s2n = np.full(4, 100.0)
        star = bsf.bdf_extended_class_dr3gold(np.full(4, 0.0), s2n)
        gal = bsf.bdf_extended_class_dr3gold(np.full(4, 2.0), s2n)
        assert np.all(star <= 1)
        assert np.all(gal >= 3)

    def test_sentinels_return_the_no_data_class(self, bsf):
        ext = bsf.bdf_extended_class_dr3gold(
            np.array([bsf.SENTINEL, 0.1, 0.1]),
            np.array([100.0, bsf.SENTINEL, -1.0]),
        )
        assert np.all(ext == -9)


class TestSchemaRegistry:
    def test_both_decam_surveys_are_registered(self, bsf):
        assert set(bsf.SCHEMAS) == {"delve", "des_y6"}

    def test_per_survey_magnitude_defaults(self, bsf):
        """des_y6 must default to bdf: its truth side is the deep-field BDF mag,
        so comparing against meas_bdf_mag keeps (obs - true) apples-to-apples."""
        assert bsf.SCHEMAS["delve"](mag_kind=None).mag_kind == "psf"
        assert bsf.SCHEMAS["des_y6"](mag_kind=None).mag_kind == "bdf"

    def test_des_truth_mag_default_includes_extinction(self, bsf):
        """`raw` has the smaller median offset but a 12x wider per-tile spread,
        because extinction varies tile to tile -- so it must not be the default."""
        assert bsf.SCHEMAS["des_y6"]().truth_mag_kind == "deredden"

    def test_des_requires_its_extra_inputs(self, bsf):
        """A des_y6 run without the measured file / labels / surrogate cannot be
        correct, so it must fail loudly rather than silently degrade."""
        import argparse

        schema = bsf.SCHEMAS["des_y6"]()
        args = argparse.Namespace(catalog="a.h5", measured=None,
                                  truth_labels=None, surrogate=None)
        with pytest.raises(SystemExit, match="--measured"):
            schema.open(args)


class TestProductColumnContract:
    """The reducer's CSV headers must be exactly what streamobs reads back.

    SurveyFactory loads the misclassification product inside a bare
    ``except: pass``, so a mis-named column does not raise -- the product just
    silently becomes None. That is the pre-existing des_yr6 failure mode, so it
    has to be pinned by a round-trip rather than by eyeballing the header.
    """

    @staticmethod
    def _write(tmp_path, header, cols):
        import numpy as np

        path = tmp_path / "curve.csv"
        np.savetxt(path, np.column_stack(cols), delimiter=",",
                   header=header, fmt="%.6f", comments="")
        return str(path)

    def test_misclass_header_round_trips_through_streamobs(self, tmp_path):
        from streamobs import surveys

        delta = np.arange(-8.0, 1.0, 0.25)
        rate = np.clip(0.01 + 0.02 * (delta + 8) / 9, 0, 1)
        path = self._write(tmp_path, "mag_g,delta_mag,missclassification_eff",
                           [delta + 25.0, delta, rate])
        fn = surveys.SurveyFactory.set_completeness(
            path, delta_saturation=-9.0, selection="missclassified")
        assert fn(-4.0) == pytest.approx(np.interp(-4.0, delta, rate), abs=1e-6)

    def test_efficiency_header_round_trips_through_streamobs(self, tmp_path):
        from streamobs import surveys

        delta = np.arange(-8.0, 1.0, 0.25)
        det = np.clip(1.0 - 0.05 * (delta + 8) / 9, 0, 1)
        cls_ = np.full(delta.size, 0.92)
        path = self._write(
            tmp_path,
            "mag_g,delta_mag,detection_eff,classification_eff,classification_detection_eff",
            [delta + 25.0, delta, det, cls_, det * cls_])
        for sel, want in (("detected", det), ("classified", cls_),
                          ("both", det * cls_)):
            fn = surveys.SurveyFactory.set_completeness(
                path, delta_saturation=-9.0, selection=sel)
            assert fn(-4.0) == pytest.approx(np.interp(-4.0, delta, want), abs=1e-6)

    def test_reducer_emits_the_expected_header(self, bsf):
        """Guard the literal the reducer writes, since the loader will not."""
        src = SCRIPT.read_text()
        assert "missclassification_eff" in src, (
            "the misclassification product must use streamobs' column name")
        assert "misclass_rate" not in src, (
            "misclass_rate loads as None -- streamobs reads missclassification_eff")


class TestSigmaConvention:
    def test_sig_sn5_is_the_5sigma_magnitude_error(self, bsf):
        assert bsf.SIG_SN5 == pytest.approx(0.21715, abs=1e-5)
