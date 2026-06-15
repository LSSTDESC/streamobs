"""
Column-name helpers for injected catalogs.

These centralize the naming convention so the injector is not hard-coded to
specific bands. A single, uniform ``<band>_true`` / ``<band>_obs`` /
``<band>_err`` convention is used everywhere; the only difference between the
single- and multi-survey cases is an optional survey prefix:

- **Single-survey** (``survey=None``): ``<band>_true`` (true / noiseless),
  ``<band>_obs`` (observed / noisy), ``<band>_err`` (reported error),
  ``flag_observed``.
- **Multi-survey** (``survey="roman"``, ``"lsst"``, ...): ``<survey>_<band>_true``,
  ``<survey>_<band>_obs``, ``<survey>_<band>_err``, ``<survey>_flag_observed``.
  Used by :class:`~streamobs.observed.MultiSurveyInjector` so each band's
  columns are namespaced by survey.

.. note::
   This convention intentionally **drops** the historical ``mag_<band>`` /
   ``mag_<band>_obs`` / ``magerr_<band>`` names — it is not backward compatible
   with catalogs written by older ``streamobs`` versions.
"""


def true_col(band, survey=None):
    """Column holding the *true* (noiseless) apparent magnitude for ``band``."""
    return f"{survey}_{band}_true" if survey else f"{band}_true"


def obs_col(band, survey=None):
    """Column holding the *observed* (noisy) magnitude for ``band``."""
    return f"{survey}_{band}_obs" if survey else f"{band}_obs"


def err_col(band, survey=None):
    """Column holding the reported magnitude error for ``band``."""
    return f"{survey}_{band}_err" if survey else f"{band}_err"


def flag_col(survey=None):
    """Column holding the detection flag (band-independent)."""
    return f"{survey}_flag_observed" if survey else "flag_observed"


def perfect_flag_col(survey=None):
    """Column holding the perfect star/galaxy-separation flag (band-independent)."""
    return f"{survey}_flag_perfect_galstarsep" if survey else "flag_perfect_galstarsep"
