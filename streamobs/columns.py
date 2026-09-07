"""
Column-name helpers for injected catalogs.

These centralize the naming convention so the injector is not hard-coded to
specific bands. Injected catalogs are **always** survey-namespaced —
``<name>_<band>_true`` (true / noiseless), ``<namespace>_<band>_obs``
(observed / noisy), ``<namespace>_<band>_err`` (reported error), and
``<namespace>_flag_observed`` — produced by
:class:`~streamobs.observed.StreamInjector` whether it serves one survey or
several. The namespace is the survey's :attr:`~streamobs.surveys.Survey.namespace`
(``{name}_{release}``).

**Observed, error, and flag columns carry the full namespace**, release included
(e.g. ``lsst_yr5_r_obs``, ``roman_dc2_F158_obs``), so the same survey at two
releases never collides.

**True-magnitude columns are the deliberate exception**: a true (noiseless)
magnitude does not depend on which survey release observed it, so these are keyed
on the survey **name only** and the release is dropped (e.g. ``roman_F158_true``,
shared across ``roman_dc2`` and ``roman_hlwas_*``; ``lsst_r_true``, shared across
every LSST release). :func:`true_col` therefore takes the leading
``{name}`` component of whatever namespace it is handed. See
``docs/source/column_convention.md`` for the full table.

The ``survey`` argument therefore identifies the namespace. ``survey=None`` is
retained only as a low-level fallback that yields the bare ``<band>_…`` /
``flag_observed`` names; the injector itself never uses it.

.. note::
   This convention intentionally **drops** the historical ``mag_<band>`` /
   ``mag_<band>_obs`` / ``magerr_<band>`` names — it is not backward compatible
   with catalogs written by older ``streamobs`` versions.
"""


def true_col(band, survey_namespace=None):
    """Column holding the *true* (noiseless) apparent magnitude for ``band``.

    Unlike :func:`obs_col` / :func:`err_col` / :func:`flag_col`, this drops the
    release and keys on the survey **name** only — a true magnitude is
    release-independent, so ``lsst_yr1`` and ``lsst_yr5`` both yield
    ``lsst_<band>_true``.
    """

    # Split survey_namespace "{name}_{release}" into survey and release if needed
    if isinstance(survey_namespace, str):
        survey_name = survey_namespace.split("_")[0]
    else:
        survey_name = None

    return f"{survey_name}_{band}_true" if survey_name else f"{band}_true"


def obs_col(band, survey_namespace=None):
    """Column holding the *observed* (noisy) magnitude for ``band``."""
    return f"{survey_namespace}_{band}_obs" if survey_namespace else f"{band}_obs"


def err_col(band, survey_namespace=None):
    """Column holding the reported magnitude error for ``band``."""
    return f"{survey_namespace}_{band}_err" if survey_namespace else f"{band}_err"


def flag_col(survey_namespace=None):
    """Column holding the detection flag (band-independent)."""
    return f"{survey_namespace}_flag_observed" if survey_namespace else "flag_observed"


def perfect_flag_col(survey_namespace=None):
    """Column holding the perfect star/galaxy-separation flag (band-independent)."""
    return (
        f"{survey_namespace}_flag_perfect_galstarsep"
        if survey_namespace
        else "flag_perfect_galstarsep"
    )


def detected_flag_col(survey_namespace=None):
    """Column holding the DETECTION-ONLY flag (band-independent).

    ``P(detected)`` with no star/galaxy classification term, for BOTH source types --
    stars via the stellar detection efficiency, galaxies via the galaxy detection
    efficiency. Emitted only when ``inject(detection_only_flag=True)``.

    Distinct from :func:`perfect_flag_col`: "perfect star/galaxy separation" means a
    galaxy is correctly REJECTED, so that flag is stars-only and would be False for every
    galaxy. This one asks whether the object was detected at all, which is what a joint
    selection needs from the survey that is not supplying morphology.
    """
    return (
        f"{survey_namespace}_flag_detected" if survey_namespace else "flag_detected"
    )


def coupling_col(survey_namespace=None):
    """Internal column holding this survey's pre-drawn selection uniform.

    Written by :meth:`~streamobs.observed.StreamInjector.inject` only when
    ``classification_coupling`` is not ``"independent"``, consumed by
    ``_inject_one_survey``, and dropped before ``inject`` returns — it is a
    private scratch column, not part of the output convention.
    """
    return (
        f"_u_select_{survey_namespace}" if survey_namespace else "_u_select"
    )
