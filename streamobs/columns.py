"""
Column-name helpers for injected catalogs.

These centralize the naming convention so the injector is not hard-coded to
specific bands. Injected catalogs are **always** survey-namespaced —
``<namespace>_<band>_true`` (true / noiseless), ``<namespace>_<band>_obs``
(observed / noisy), ``<namespace>_<band>_err`` (reported error), and
``<namespace>_flag_observed`` — produced by
:class:`~streamobs.observed.StreamInjector` whether it serves one survey or
several. The namespace is the survey's :attr:`~streamobs.surveys.Survey.namespace`
(``{name}_{release}``), so it includes the release on every column kind
(e.g. ``lsst_yr5_r_obs``, ``roman_dc2_F158_obs``) and the same survey at two
releases never collides.

The ``survey`` argument therefore identifies the namespace. ``survey=None`` is
retained only as a low-level fallback that yields the bare ``<band>_…`` /
``flag_observed`` names; the injector itself never uses it.

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
