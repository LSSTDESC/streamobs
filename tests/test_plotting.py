"""
tests/test_plotting.py
======================
Smoke test for ``streamobs.plotting.plot_inject``.

Guards against the regression where ``plot_inject`` read bare, non-namespaced
columns (``flag_observed``, ``r_obs``, ...) and so raised "Missing required
columns" on the injector's now-namespaced (``{name}_{release}``) output.
"""

import matplotlib

matplotlib.use("Agg")  # headless backend for tests

import pytest

from streamobs.plotting import plot_inject


@pytest.mark.observed
class TestPlotInject:
    def test_plot_inject_namespaced_columns(
        self, mock_injector, stream_catalog, verbose
    ):
        """plot_inject consumes the namespaced injector output without error."""
        cat = mock_injector.inject(stream_catalog, bands=["g", "r"], verbose=verbose)
        # mock_injector is LSST/yr4 -> namespace lsst_yr4; columns are lsst_yr4_*.
        fig, ax = plot_inject(cat, mock_injector.survey, bands=["g", "r"])
        assert fig is not None
        assert len(ax) == 3
        matplotlib.pyplot.close(fig)
