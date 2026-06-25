import pandas as pd
import numpy as np
from matplotlib.path import Path
from ugali.analysis.isochrone import factory as isochrone_factory


"""
Match filter module for stellar stream analysis.

This module provides functions to create isochrone-based matched-filter 
polygons and select stars that fall within these polygons in color-magnitude space.

The matched-filter technique works by defining a polygon in color-magnitude 
space (CMD) that encompasses the expected stellar locus for a given stellar 
population (defined by age, metallicity, and distance). Stars falling within 
this polygon are selected as candidate members.

Functions
---------
compute_photometric_error
    Compute magnitude-dependent photometric errors.
build_isochrone_polygon
    Construct a CMD polygon from an isochrone model.
select_stars_by_isochrone
    Select stars within an isochrone-based CMD polygon.
select_stars_with_error_model
    Select stars with a custom photometric error model.

Notes
-----
This module uses the Dotter2008 isochrone models from the ugali package.
"""


# =============================================================================
# Default isochrone and filter parameters
# =============================================================================
# From DES 2018: https://arxiv.org/pdf/1801.03097

DEFAULT_AGE_GYR = 13.0  # Stellar population age in Gigayears
DEFAULT_METALLICITY = 0.0002  # Metallicity (Z), typical for old halo populations
DEFAULT_DISTANCE_MODULUS_SPREAD = 0.5  # Half-width of distance modulus range [mag]
DEFAULT_COLOR_SPREAD = [0.05, 0.05]  # Color padding [blue_side, red_side] [mag] # note: smaller than DES2018, symmetric ideal case
DEFAULT_ERROR_MULTIPLIER = [2.0, 2.0]  # Error scaling [blue_side, red_side]

# Absolute magnitude of Main Sequence Turn-Off for Dotter2008 isochrone
# WARNING: This value is model-dependent and hardcoded for Dotter2008
MSTO_ABSOLUTE_MAG = 3.5 # it is the same for Marigo2017
DEFAULT_RGB_CLIP_MAG = 0.2  # Clip RGB at (MSTO - 0.2) mag

default_errors_des2018 =  {'baseline_error': 0.001, 'exp_pivot': 27.09, 'exp_scale': 1.09}
default_errors =  {'baseline_error': 0.004775486092612673, 'exp_pivot': 28.421419633248796, 'exp_scale': 1.0011829218659076} # fitted on DC2 mag err

def error_model(magnitude, baseline_error=default_errors['baseline_error'], exp_pivot=default_errors['exp_pivot'], exp_scale=default_errors['exp_scale'], verbose=False):
    """
    Compute the median photometric error as a function of magnitude.
    
    Uses an exponential error model calibrated for typical survey data.
    
    Parameters
    ----------
    magnitude : float or array-like
        Apparent magnitude(s) for which to compute the error.
    
    Returns
    -------
    error : float or array-like
        Photometric error(s) corresponding to the input magnitude(s).
    """
    if verbose:
        print("Using following parameters values:",baseline_error, exp_pivot, exp_scale )
    return baseline_error + np.exp((magnitude - exp_pivot) / exp_scale)


def build_match_filter(
    distance_modulus,
    age=DEFAULT_AGE_GYR,
    metallicity=DEFAULT_METALLICITY,
    distance_modulus_spread=DEFAULT_DISTANCE_MODULUS_SPREAD,
    color_spread=DEFAULT_COLOR_SPREAD,
    error_multiplier=DEFAULT_ERROR_MULTIPLIER,
    rgb_clip_mag=DEFAULT_RGB_CLIP_MAG,
    color_cut=True,
    verbose=False,
    error_kwargs = {},
    survey = "lsst",
    isochrone_model = 'Marigo2017',
):
    """
    Build an isochrone matched-filter polygon in color-magnitude space.
    
    Creates a closed polygon that encompasses the expected locus of stars
    at a given distance, accounting for photometric errors and distance spread.
    The polygon is constructed by tracing the isochrone at both the near and
    far distance limits, with additional padding based on photometric errors.

    The polygon structure in CMD space:
    
        Color (g-r) -->
        
        ^                    * (MSTO region, faint end)
        |                   / \\
        |   Blue edge -->  /   \\  <-- Red edge
        |                 /     \\
        Mag (g)          /       \\
        |               /         \\
        |              *-----------*  (RGB, bright end)
        v
    
    The blue edge is offset blueward (negative color) from the isochrone,
    and the red edge is offset redward (positive color). Both edges include
    magnitude-dependent photometric error corrections.

    Parameters
    ----------
    distance_modulus : float
        Distance modulus (m - M) of the stellar population in magnitudes.
        Relates apparent magnitude (m) to absolute magnitude (M).
    age : float, optional
        Isochrone age in Gigayears (Gyr). Default is 12.0 Gyr,
        appropriate for old halo populations.
    metallicity : float, optional
        Isochrone metallicity as mass fraction Z (not [Fe/H]).
        Default is 0.0006 (~[Fe/H] = -2.0), typical for metal-poor halo stars.
    distance_modulus_spread : float, optional
        Half-width of the distance modulus range to consider, in magnitudes.
        Accounts for depth along the line of sight. Default is 0.5 mag.
    color_spread : list of float, optional
        Additive color padding as [blue_side, red_side] in magnitudes.
        Adds a fixed offset to broaden the polygon. Default is [0.1, 0.1].
    error_multiplier : list of float, optional
        Multiplicative scaling for photometric errors as [blue_side, red_side].
        Controls how much the polygon expands based on magnitude-dependent
        errors. Default is [5.0, 1.0] (more expansion on blue side).
    rgb_clip_mag : float or None, optional
        If provided, clips the Red Giant Branch at this magnitude offset
        from the Main Sequence Turn-Off (MSTO). Stars brighter than
        (MSTO_ABSOLUTE_MAG - rgb_clip_mag) are excluded. Useful to avoid
        contamination from foreground giants.
    color_cut : bool, optional
        If True, restricts the match filter to stars with g-r < 1.25 - error(g),
        where error(g) is the magnitude-dependent photometric error.
        This removes very red stars from the filter. Default is True.

    Returns
    -------
    polygon_vertices : ndarray, shape (N, 2)
        Vertices of the closed polygon as (color, magnitude) pairs.
        Column 0: g-r color
        Column 1: apparent g-band magnitude
        Suitable for use with matplotlib.path.Path.contains_points().

    Examples
    --------
    >>> vertices = build_isochrone_polygon(distance_modulus=16.5, age=12.0)
    >>> polygon = Path(vertices)
    >>> inside = polygon.contains_points(star_coords)
    
    Notes
    -----
    - Uses the Marigo2017 isochrone models from the ugali package.
    - The MSTO absolute magnitude is hardcoded for Dotter2008 isochrones.
    - Photometric errors increase exponentially toward fainter magnitudes,
      resulting in a wider polygon at the faint end.
    """
    # --- Generate isochrone model ---
    isochrone = isochrone_factory(isochrone_model, age=age, z=metallicity, survey = survey)
    isochrone_color = isochrone.color  # Intrinsic g-r color
    isochrone_absolute_mag = isochrone.mag  # Absolute g-band magnitude

    # --- Optional: Clip the Red Giant Branch ---
    if rgb_clip_mag is not None:
        rgb_bright_limit = MSTO_ABSOLUTE_MAG - rgb_clip_mag
        if verbose:
            print(f"Clipping RGB at absolute mag = {rgb_bright_limit:.2f}")
        valid_stars_mask = isochrone_absolute_mag > rgb_bright_limit
        isochrone_color = isochrone_color[valid_stars_mask]
        isochrone_absolute_mag = isochrone_absolute_mag[valid_stars_mask]

    # --- Compute apparent magnitudes at distance boundaries ---
    half_spread = distance_modulus_spread / 2.0
    apparent_mag_near = isochrone_absolute_mag + distance_modulus - half_spread
    apparent_mag_far = isochrone_absolute_mag + distance_modulus + half_spread

    # --- Unpack padding and error multiplier parameters ---
    blue_color_padding, red_color_padding = color_spread
    blue_error_multiplier, red_error_multiplier = error_multiplier

    # --- Compute magnitude-dependent photometric errors ---
    photometric_error_far = error_model(apparent_mag_far,**error_kwargs)
    photometric_error_near = error_model(apparent_mag_near, **error_kwargs)

    # --- Build polygon edges ---
    # Red edge: trace isochrone forward, offset redward
    red_edge_color = (
        isochrone_color
        + red_error_multiplier * photometric_error_far
        + red_color_padding
    )

    # Blue edge: trace isochrone backward (reversed), offset blueward
    blue_edge_color = (
        isochrone_color[::-1]
        - blue_error_multiplier * photometric_error_near[::-1]
        - blue_color_padding
    )
    
    # --- Apply color cut if requested ---
    if color_cut:
        # Restrict to g - r < 1.25 - error(g) to remove very red stars
        color_base = 1.25
        color_limit_far = color_base - photometric_error_far
        red_edge_color = np.minimum(red_edge_color, color_limit_far)

    # --- Concatenate edges to form closed polygon ---
    polygon_color = np.concatenate([red_edge_color, blue_edge_color])
    polygon_apparent_mag = np.concatenate([
        isochrone_absolute_mag,
        isochrone_absolute_mag[::-1]
    ]) + distance_modulus

    return np.column_stack([polygon_color, polygon_apparent_mag])


def is_in_match_filter(mag_g, mag_r, polygon_vertices=None, match_filter_params=None, verbose = False):
    mag_g = pd.to_numeric(mag_g, errors="coerce")
    mag_r = pd.to_numeric(mag_r, errors="coerce")
    
    # Explicitly exclude NaN values (from 'BAD_MAG' or invalid data)
    valid_mask = ~(np.isnan(mag_g) | np.isnan(mag_r))
    
    color_mag_coords = np.column_stack([
        mag_g - mag_r,
        mag_g
    ])

    if polygon_vertices is None:
        if match_filter_params is None:
            raise ValueError("Either polygon_vertices or match_filter_params must be provided.")
        polygon_vertices = build_match_filter(**match_filter_params)
    
    polygon_path = Path(polygon_vertices)
    # Only check valid (non-NaN) points
    selection_mask = np.zeros(len(mag_g), dtype=bool)
    selection_mask[valid_mask] = polygon_path.contains_points(color_mag_coords[valid_mask])

    if verbose:
        n_valid = np.sum(valid_mask)
        n_selected = np.sum(selection_mask)
        selection_fraction = n_selected / len(mag_g) * 100 if len(mag_g) > 0 else 0
        print(f'Match filter selects {selection_fraction:.2f}% of all stars ({n_selected}/{len(mag_g)}, {len(mag_g)-n_valid} stars with NaN)') 
    
    return selection_mask

