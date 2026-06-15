"""
Column-name helpers for injected catalogs.

These centralize the naming convention so the injector is not hard-coded to
specific bands (``mag_r_obs``, ``magerr_g``, ...). Two conventions are
supported, selected by the ``survey`` argument:

- **Legacy / single-survey** (``survey=None``): ``mag_<band>`` (true),
  ``mag_<band>_obs`` (observed), ``magerr_<band>`` (error), ``flag_observed``.
  This is what :class:`~streamobs.observed.StreamInjector` emits and what
  downstream consumers already read, so it is preserved unchanged.
- **Multi-survey** (``survey="roman"``, ``"lsst"``, ...): ``<survey>_<band>_true``,
  ``<survey>_<band>_obs``, ``<survey>_<band>_err``, ``<survey>_flag_observed``.
  Used by the (future) ``MultiSurveyInjector`` so each band's columns are
  namespaced by survey.
"""


def true_col(band, survey=None):
    """Column holding the *true* (noiseless) apparent magnitude for ``band``."""
    return f"{survey}_{band}_true" if survey else f"mag_{band}"


def obs_col(band, survey=None):
    """Column holding the *observed* (noisy) magnitude for ``band``."""
    return f"{survey}_{band}_obs" if survey else f"mag_{band}_obs"


def err_col(band, survey=None):
    """Column holding the reported magnitude error for ``band``."""
    return f"{survey}_{band}_err" if survey else f"magerr_{band}"


def flag_col(survey=None):
    """Column holding the detection flag (band-independent)."""
    return f"{survey}_flag_observed" if survey else "flag_observed"
