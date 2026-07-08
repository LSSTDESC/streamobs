from __future__ import annotations

from typing import Callable, Tuple

import numpy as np
import pandas as pd
import scipy.interpolate
from matplotlib.path import Path
from ugali.analysis.isochrone import factory as isochrone_factory

from streamobs.model import ROMAN_VEGA_TO_AB

"""
Match filter module for stellar stream analysis.

This module provides functions to create isochrone-based matched-filter
polygons and select stars that fall within these polygons in color-magnitude
space.

The matched-filter technique works by defining a polygon in color-magnitude
space (CMD) that encompasses the expected stellar locus for a given stellar
population (defined by age, metallicity, and distance). Stars falling within
this polygon are selected as candidate members.

The polygon edges are constructed via **spline interpolation** of the
isochrone locus, which guarantees bounded, well-behaved polygon vertices
even at the faint end where the raw exponential error model would otherwise
blow up.

Functions
---------
build_match_filter
    Construct a CMD polygon from an isochrone model using spline boundaries.
is_in_match_filter
    Select objects inside the matched-filter polygon.
build_filter_splines
    Low-level: build near/far color boundary splines from isochrone arrays.
"""


# =============================================================================
# Default isochrone and filter parameters
# =============================================================================
# From DES 2018: https://arxiv.org/pdf/1801.03097

DEFAULT_AGE_GYR = 13.0  # Stellar population age in Gigayears
DEFAULT_METALLICITY = 0.0002  # Metallicity (Z), typical for old halo populations
DEFAULT_DISTANCE_MODULUS_SPREAD = 0.5  # Half-width of distance modulus range [mag]
DEFAULT_COLOR_SPREAD = [
    0.05,
    0.05,
]  # Color padding [blue_side, red_side] [mag]
DEFAULT_ERROR_MULTIPLIER = [2.0, 2.0]  # Error scaling [blue_side, red_side]

# Absolute magnitude of Main Sequence Turn-Off for Marigo2017 isochrone
MSTO_ABSOLUTE_MAG = 3.5
DEFAULT_RGB_CLIP_MAG = 0.2  # Clip RGB at (MSTO - 0.2) mag

default_errors_des2018 = {
    "baseline_error": 0.001,
    "exp_pivot": 27.09,
    "exp_scale": 1.09,
}
default_errors = {
    "baseline_error": 0.004775486092612673,
    "exp_pivot": 28.421419633248796,
    "exp_scale": 1.0011829218659076,
}  # fitted on DC2 mag err

# Per-survey defaults for color_min/color_max clipping applied to spline edges.
# These keep the polygon bounded regardless of the exponential error model.
# For LSST g-r the MS+RGB locus lives in [~0, ~0.7]; we allow [-0.5, 1.5].
# For Roman NIR bands (e.g. F106-F158 in AB) the locus is [-0.15, ~0.35] after
# Vega→AB conversion; we clip at [-0.5, 1.5] to prevent error blow-up.
_SURVEY_COLOR_LIMITS = {
    "lsst": (-0.5, 1.5),
    "des": (-0.5, 1.5),
    "sdss": (-0.5, 1.5),
    "roman": (-0.5, 1.5),
    "default": (-0.5, 1.5),
}

# Per-survey bright absolute magnitude cut (points brighter are noisier/RGB clutter)
_SURVEY_ABS_MAG_MIN = {
    "lsst": 2.9,
    "des": 2.9,
    "sdss": 2.9,
    "roman": 2.9,
    "default": 2.9,
}

# Default faint apparent-magnitude limit used when sampling the polygon grid
_DEFAULT_APP_MAG_MAX = 30.0

# Number of magnitude grid points when sampling splines into polygon vertices
_POLYGON_GRID_POINTS = 200


def error_model(
    magnitude,
    baseline_error=default_errors["baseline_error"],
    exp_pivot=default_errors["exp_pivot"],
    exp_scale=default_errors["exp_scale"],
    verbose=False,
):
    """Compute the median photometric error as a function of magnitude.

    Uses an exponential error model calibrated for typical survey data.

    Parameters
    ----------
    magnitude : float or array-like
        Apparent magnitude(s) for which to compute the error.
    baseline_error : float, optional
        Constant floor of the error curve.
    exp_pivot : float, optional
        Pivot magnitude of the exponential growth.
    exp_scale : float, optional
        Scale length of the exponential.
    verbose : bool, optional
        Print parameter values.

    Returns
    -------
    error : float or array-like
        Photometric error(s) corresponding to the input magnitude(s).
    """
    if verbose:
        print(
            "Using following parameters values:", baseline_error, exp_pivot, exp_scale
        )
    return baseline_error + np.exp((magnitude - exp_pivot) / exp_scale)


def build_filter_splines(
    iso_color: np.ndarray,
    iso_mag: np.ndarray,
    *,
    mu: float = 0.0,
    mag_is_absolute: bool = True,
    abs_mag_min: float | None = 2.9,
    app_mag_max: float | None = None,
    color_min: float = -0.5,
    color_max: float = 1.5,
    dmu: float = 0.5,
    C: tuple | list = (0.05, 0.05),
    E: tuple | list = (2.0, 2.0),
    err: Callable[[np.ndarray], np.ndarray] | None = None,
) -> Tuple[Callable, Callable]:
    """Build matched-filter boundary splines from precomputed isochrone arrays.

    Constructs two interpolating functions — ``spline_near`` (blue/near edge)
    and ``spline_far`` (red/far edge) — that give the color boundaries of the
    matched-filter region as a function of **apparent magnitude**.

    The splines are built from the isochrone locus after applying magnitude
    cuts and color clipping, so they remain bounded even where the raw
    exponential error model would diverge.

    Parameters
    ----------
    iso_color : array_like
        Isochrone color array (e.g. g−r) for the locus.
    iso_mag : array_like
        Isochrone magnitude array.  Interpreted as absolute magnitude when
        ``mag_is_absolute=True``, apparent magnitude otherwise.
    mu : float, optional
        Distance modulus (used only when ``mag_is_absolute=True``).
    mag_is_absolute : bool, optional
        Whether ``iso_mag`` is in absolute magnitudes.
    abs_mag_min : float or None, optional
        Bright absolute-magnitude cut applied before spline fitting.
    app_mag_max : float or None, optional
        Faint apparent-magnitude cut applied before spline fitting.
    color_min, color_max : float, optional
        Hard clip bounds applied to the near/far color envelopes before
        spline fitting.  These prevent runaway exponential values from
        entering the polygon.
    dmu : float, optional
        Half-width of the distance-modulus spread used to shift magnitudes
        for the near/far error evaluation.
    C : tuple/list [near, far], optional
        Additive color padding on the near (blue) and far (red) sides.
    E : tuple/list [near, far], optional
        Multiplicative error scaling for near and far sides.
    err : callable or None, optional
        Function ``err(apparent_mag) -> mag_error``.  When ``None`` the
        default exponential model is used.

    Returns
    -------
    spline_near, spline_far : callable
        Functions mapping apparent magnitude → color boundary.  They return
        ``nan`` outside the valid magnitude range so that callers can
        detect out-of-bounds regions without needing explicit range checks.
    """
    if err is None:
        err = error_model

    # Allow E to be scalar (backward compat) or [near, far]
    if np.isscalar(E):
        E_near = E_far = float(E)
    else:
        E_near, E_far = float(E[0]), float(E[1])

    if np.isscalar(C):
        C_near = C_far = float(C)
    else:
        C_near, C_far = float(C[0]), float(C[1])

    iso_color = np.asarray(iso_color, dtype=float)
    iso_mag = np.asarray(iso_mag, dtype=float)
    if iso_color.shape != iso_mag.shape:
        raise ValueError("iso_color and iso_mag must have the same shape")

    if mag_is_absolute:
        abs_mag = iso_mag
        app_mag = iso_mag + float(mu)
    else:
        app_mag = iso_mag
        abs_mag = iso_mag - float(mu)

    # Keep only finite points and apply magnitude cuts
    sel = np.isfinite(iso_color) & np.isfinite(app_mag)
    if abs_mag_min is not None:
        sel &= abs_mag > float(abs_mag_min)
    if app_mag_max is not None:
        sel &= app_mag < float(app_mag_max)

    iso_color = iso_color[sel]
    app_mag = app_mag[sel]

    if iso_color.size < 2:
        raise ValueError("Not enough valid points to build filter splines")

    # Compute near (blue) and far (red) magnitude offsets from dmu spread
    mnear = app_mag - dmu / 2.0
    mfar = app_mag + dmu / 2.0

    # Build color envelopes; clip to [color_min, color_max] to suppress
    # exponential blow-up at the faint end
    color_far = np.clip(
        iso_color + E_far * err(mfar) + C_far,
        color_min,
        color_max,
    )
    color_near = np.clip(
        iso_color - E_near * err(mnear) - C_near,
        color_min,
        color_max,
    )

    # interp1d requires strictly increasing x
    order = np.argsort(app_mag)
    x = app_mag[order]
    y_near = color_near[order]
    y_far = color_far[order]

    # Drop duplicate x values (keep first occurrence)
    x_unique, unique_idx = np.unique(x, return_index=True)
    y_near = y_near[unique_idx]
    y_far = y_far[unique_idx]

    if x_unique.size < 2:
        raise ValueError("Not enough unique magnitude points to build filter splines")

    spline_near = scipy.interpolate.interp1d(
        x_unique, y_near, bounds_error=False, fill_value=np.nan
    )
    spline_far = scipy.interpolate.interp1d(
        x_unique, y_far, bounds_error=False, fill_value=np.nan
    )

    return spline_near, spline_far


def build_match_filter(
    distance_modulus,
    age=DEFAULT_AGE_GYR,
    metallicity=DEFAULT_METALLICITY,
    distance_modulus_spread=DEFAULT_DISTANCE_MODULUS_SPREAD,
    color_spread=DEFAULT_COLOR_SPREAD,
    error_multiplier=DEFAULT_ERROR_MULTIPLIER,
    rgb_clip_mag=DEFAULT_RGB_CLIP_MAG,
    color_cut=True,
    red_cap_mode="error_shrunk",
    verbose=False,
    error_kwargs={},
    survey="lsst",
    isochrone_model="Marigo2017",
    band_1="g",
    band_2="r",
):
    """Build an isochrone matched-filter polygon in color-magnitude space.

    Creates a closed polygon that encompasses the expected locus of stars
    at a given distance, accounting for photometric errors and distance spread.

    **Implementation note:** polygon edges are constructed via spline
    interpolation of the isochrone locus rather than raw concatenation of
    error-broadened arrays.  This keeps all polygon vertices strictly bounded
    (within the ``color_min``/``color_max`` clip derived from the survey),
    preventing the exponential error model from producing degenerate colors
    at the faint end (the old approach produced blue-edge colors of ~ -113
    at distance_modulus ≈ 16.8 with default LSST parameters).

    Parameters
    ----------
    distance_modulus : float
        Distance modulus (m − M) of the stellar population in magnitudes.
    age : float, optional
        Isochrone age in Gigayears (Gyr).
    metallicity : float, optional
        Isochrone metallicity as mass fraction Z.
    distance_modulus_spread : float, optional
        Half-width of the distance modulus range to consider [mag].
    color_spread : list of float, optional
        Additive color padding as [blue_side, red_side] in magnitudes.
    error_multiplier : list of float, optional
        Multiplicative scaling for photometric errors as [blue_side, red_side].
    rgb_clip_mag : float or None, optional
        If provided, clips the Red Giant Branch at this magnitude offset
        from the Main Sequence Turn-Off (MSTO).
    color_cut : bool, optional
        If True, the red spline edge is additionally capped so the filter
        does not extend to very red stars (see ``red_cap_mode``).
    red_cap_mode : {'error_shrunk', 'constant'}, optional
        How the red-edge cap behaves when ``color_cut=True``:

        - ``'error_shrunk'`` (default): cap at
          ``(color_max - 0.25) - error(mag_far)``, so the cap moves blueward
          with the photometric error and the filter narrows (eventually
          closes) at the faint end — the historical pre-spline behaviour,
          now applied to the spline-sampled edges.  The far edge is floored
          at the near edge so the polygon stays simple instead of
          self-crossing.
        - ``'constant'``: cap at ``color_max - 0.25`` at all magnitudes.
    verbose : bool, optional
        Print diagnostic information.
    error_kwargs : dict, optional
        Extra keyword arguments forwarded to :func:`error_model`.
    survey : str, optional
        Survey name passed to the ugali isochrone factory.  Also controls
        per-survey color clip limits and absolute-magnitude cuts.
    isochrone_model : str, optional
        Ugali isochrone model name (e.g. ``'Marigo2017'``).
    band_1, band_2 : str, optional
        Bands defining the CMD: color = band_1 − band_2, magnitude = band_1.
        Defaults (``'g'``, ``'r'``) preserve the historical LSST behaviour.
        For Roman pass e.g. ``band_1='F106'``, ``band_2='F158'``.
        Roman isochrone magnitudes are returned by ugali in Vega and are
        converted to AB here (via the ``ROMAN_VEGA_TO_AB`` table that the
        injection path uses), so the filter selects on the same photometric
        system as the injected catalogs.

    Returns
    -------
    polygon_vertices : ndarray, shape (N, 2)
        Vertices of the closed polygon as (color, magnitude) pairs.
        Column 0: color (band_1 − band_2)
        Column 1: apparent band_1 magnitude
        Suitable for use with ``matplotlib.path.Path.contains_points()``.

    Notes
    -----
    - Uses Marigo2017 isochrone models from the ugali package by default.
    - The Vega→AB conversion for Roman bands is applied before any spline
      fitting, so the polygon is in the same photometric system as the data.
    """
    survey_lower = survey.lower() if isinstance(survey, str) else "default"

    # Look up per-survey color clip limits
    color_min, color_max = _SURVEY_COLOR_LIMITS.get(
        survey_lower, _SURVEY_COLOR_LIMITS["default"]
    )
    abs_mag_min = _SURVEY_ABS_MAG_MIN.get(survey_lower, _SURVEY_ABS_MAG_MIN["default"])

    # --- Generate isochrone model ---
    isochrone = isochrone_factory(
        isochrone_model,
        age=age,
        z=metallicity,
        survey=survey,
        band_1=band_1,
        band_2=band_2,
    )

    # ugali returns Roman isochrones in Vega; our catalogs/maglim maps are AB.
    # AB = Vega + offset (a no-op for non-Roman bands) — same convention as
    # StreamInjector._to_ab, so the filter and the injected stream share one
    # photometric system.
    ab_1 = ROMAN_VEGA_TO_AB.get(band_1, 0.0)
    ab_2 = ROMAN_VEGA_TO_AB.get(band_2, 0.0)
    isochrone_color = isochrone.color + (ab_1 - ab_2)  # band_1 - band_2
    isochrone_absolute_mag = isochrone.mag + ab_1  # absolute band_1 magnitude

    # --- Optional: Clip the Red Giant Branch ---
    if rgb_clip_mag is not None:
        rgb_bright_limit = MSTO_ABSOLUTE_MAG - rgb_clip_mag
        if verbose:
            print(f"Clipping RGB at absolute mag = {rgb_bright_limit:.2f}")
        valid_stars_mask = isochrone_absolute_mag > rgb_bright_limit
        isochrone_color = isochrone_color[valid_stars_mask]
        isochrone_absolute_mag = isochrone_absolute_mag[valid_stars_mask]

    # --- Build spline-based color boundaries ---
    # The err callable accepts apparent magnitudes; we wrap error_model with
    # any user-supplied error_kwargs.
    def _err(m):
        return error_model(m, **error_kwargs)

    # Apparent magnitude range for polygon sampling
    app_mag_min_data = isochrone_absolute_mag.min() + distance_modulus
    app_mag_max_data = isochrone_absolute_mag.max() + distance_modulus
    # Clip to a sane range
    app_mag_min_data = max(app_mag_min_data, distance_modulus + abs_mag_min)
    app_mag_max_data = min(app_mag_max_data, _DEFAULT_APP_MAG_MAX)

    spline_near, spline_far = build_filter_splines(
        iso_color=isochrone_color,
        iso_mag=isochrone_absolute_mag,
        mu=distance_modulus,
        mag_is_absolute=True,
        abs_mag_min=abs_mag_min,
        app_mag_max=None,  # no faint hard cut; let the spline domain govern
        color_min=color_min,
        color_max=color_max,
        dmu=distance_modulus_spread,
        C=color_spread,
        E=error_multiplier,
        err=_err,
    )

    # --- Sample splines on a fine magnitude grid to get polygon vertices ---
    mag_grid = np.linspace(app_mag_min_data, app_mag_max_data, _POLYGON_GRID_POINTS)

    near_colors = spline_near(mag_grid)
    far_colors = spline_far(mag_grid)

    # Remove grid points where either spline is nan (outside interpolation range)
    valid = np.isfinite(near_colors) & np.isfinite(far_colors)
    mag_grid = mag_grid[valid]
    near_colors = near_colors[valid]
    far_colors = far_colors[valid]

    if len(mag_grid) < 2:
        raise ValueError(
            "Too few valid spline points to build a polygon. "
            "Check isochrone parameters and distance_modulus."
        )

    # Optional color cut: cap the red (far) edge so the polygon does not
    # extend to very red stars.
    if color_cut:
        red_cap = color_max - 0.25
        if red_cap_mode == "constant":
            far_colors = np.minimum(far_colors, red_cap)
        elif red_cap_mode == "error_shrunk":
            cap = red_cap - _err(mag_grid + distance_modulus_spread / 2.0)
            far_colors = np.minimum(far_colors, cap)
            # floor at the near edge: the filter closes rather than
            # self-crossing where the shrinking cap dips below it
            far_colors = np.maximum(far_colors, near_colors)
        else:
            raise ValueError(
                f"red_cap_mode must be 'constant' or 'error_shrunk', got {red_cap_mode!r}"
            )

    if verbose:
        print(f"Polygon magnitude range: {mag_grid.min():.2f} – {mag_grid.max():.2f}")
        print(f"Color range: [{near_colors.min():.3f}, {far_colors.max():.3f}]")

    # --- Concatenate forward (far/red) and backward (near/blue) edges ---
    # Forward: red edge traces increasing magnitude
    # Backward: blue edge traces decreasing magnitude (closes the polygon)
    polygon_color = np.concatenate([far_colors, near_colors[::-1]])
    polygon_mag = np.concatenate([mag_grid, mag_grid[::-1]])

    return np.column_stack([polygon_color, polygon_mag])


def is_in_match_filter(
    mag_1, mag_2, polygon_vertices=None, match_filter_params=None, verbose=False
):
    """Select objects inside the matched-filter polygon.

    ``mag_1``/``mag_2`` are the same two bands the polygon was built with
    (color = mag_1 − mag_2, magnitude axis = mag_1) — ``'g'``/``'r'`` for
    LSST, e.g. ``'F106'``/``'F158'`` for Roman.

    Parameters
    ----------
    mag_1 : array_like
        Apparent magnitude in band_1.
    mag_2 : array_like
        Apparent magnitude in band_2.
    polygon_vertices : ndarray, shape (N, 2) or None, optional
        Pre-built polygon from :func:`build_match_filter`.  If ``None``,
        ``match_filter_params`` must be given.
    match_filter_params : dict or None, optional
        Keyword arguments forwarded to :func:`build_match_filter` to build
        the polygon on the fly.
    verbose : bool, optional
        Print selection statistics.

    Returns
    -------
    selection_mask : ndarray of bool
        ``True`` for objects inside the polygon.
    """
    mag_1 = pd.to_numeric(mag_1, errors="coerce")
    mag_2 = pd.to_numeric(mag_2, errors="coerce")

    # Explicitly exclude NaN values (from 'BAD_MAG' or invalid data)
    valid_mask = ~(np.isnan(mag_1) | np.isnan(mag_2))

    color_mag_coords = np.column_stack([mag_1 - mag_2, mag_1])

    if polygon_vertices is None:
        if match_filter_params is None:
            raise ValueError(
                "Either polygon_vertices or match_filter_params must be provided."
            )
        polygon_vertices = build_match_filter(**match_filter_params)

    polygon_path = Path(polygon_vertices)
    # Only check valid (non-NaN) points
    selection_mask = np.zeros(len(mag_1), dtype=bool)
    selection_mask[valid_mask] = polygon_path.contains_points(
        color_mag_coords[valid_mask]
    )

    if verbose:
        n_valid = np.sum(valid_mask)
        n_selected = np.sum(selection_mask)
        selection_fraction = n_selected / len(mag_1) * 100 if len(mag_1) > 0 else 0
        print(
            f"Match filter selects {selection_fraction:.2f}% of all stars "
            f"({n_selected}/{len(mag_1)}, {len(mag_1)-n_valid} stars with NaN)"
        )

    return selection_mask
