"""
Shared pytest fixtures for the streamobs test suite.

All heavy objects (stream configs, survey mocks, sample DataFrames) are
defined once here and injected where needed via pytest fixtures, so every
test module can stay focused on behaviour rather than setup.
"""

import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def seed():
    """Fixed random seed for reproducibility."""
    return 42

@pytest.fixture(scope="session")
def rng(seed):
    """Random number generator initialized with a fixed seed."""
    return np.random.default_rng(seed)


# ---------------------------------------------------------------------------
# Minimal stream model config (no isochrone, no velocity — pure geometry)
# ---------------------------------------------------------------------------

@pytest.fixture
def minimal_stream_config():
    """Minimal StreamModel config: density + track only (no isochrone / velocity)."""
    return {
        "density": {
            "type": "uniform",
            "xmin": -10.0,
            "xmax": 10.0,
        },
        "track": {
            "center": {
                "type": "polynomial",
                "coeffs": [0.0],          # flat track at phi2 = 0
            },
            "spread": {
                "type": "polynomial",
                "coeffs": [0.1],          # constant width σ = 0.1 deg
            },
            "sampler": "gaussian",
        },
    }


@pytest.fixture
def stream_config_with_distance(minimal_stream_config):
    """StreamModel config that also produces a distance modulus column."""
    cfg = dict(minimal_stream_config)
    cfg["distance_modulus"] = {
        "center": {
            "type": "polynomial",
            "coeffs": [16.0],   # ~16 kpc
        },
        "spread": {
            "type": "polynomial",
            "coeffs": [0.0],
        },
        "sampler": "gaussian",
    }
    return cfg


# ---------------------------------------------------------------------------
# Tiny sample DataFrame used for injection tests
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_catalog_radec():
    """
    Small catalog with ra, dec, mag_g, mag_r columns — ready for injection
    without needing coordinate conversion or a stream model.
    """
    rng = np.random.default_rng(0)
    n = 50
    return pd.DataFrame({
        "ra":    rng.uniform(30.0, 60.0, n),
        "dec":   rng.uniform(-20.0, 0.0, n),
        "mag_g": rng.uniform(22.0, 25.0, n),
        "mag_r": rng.uniform(21.5, 24.5, n),
    })


@pytest.fixture
def sample_catalog_phi():
    """Small catalog with phi1/phi2 stream coordinates (no ra/dec yet)."""
    rng = np.random.default_rng(1)
    n = 30
    phi1 = rng.uniform(-8.0, 8.0, n)
    phi2 = rng.normal(0.0, 0.1, n)
    return pd.DataFrame({
        "phi1":  phi1,
        "phi2":  phi2,
        "mag_g": rng.uniform(22.0, 25.0, n),
        "mag_r": rng.uniform(21.5, 24.5, n),
    })


# ---------------------------------------------------------------------------
# Mock Survey — avoids any file I/O for unit tests
# ---------------------------------------------------------------------------

@pytest.fixture
def base_maglim():
    """Base magnitude limit for the mock survey."""
    return 26.0

def _make_mock_survey(name="mock", release=None, nside=64):
    """Build a synthetic Survey object backed by uniform HEALPix maps."""
    import healpy as hp
    import scipy.interpolate
    from streamobs.surveys import Survey

    npix = hp.nside2npix(nside)

    # Uniform magnitude-limit maps at 25 mag
    maglim_val = 25.0
    maglim_map = np.full(npix, maglim_val)

    # Uniform E(B-V) = 0.05 (low dust)
    ebv_map = np.full(npix, 0.05)

    # Simple completeness function: 1 for delta_mag < 0, 0 otherwise
    delta_arr = np.linspace(-12, 2, 200)
    eff_arr = np.where(delta_arr < 0.0, 1.0, 0.0)
    completeness_fn = scipy.interpolate.interp1d(
        delta_arr, eff_arr, bounds_error=False, fill_value=0.0
    )

    # Simple log-photo-error: constant log10(err) = -2  (err = 0.01 mag)
    log_err_arr = np.full_like(delta_arr, -2.0)
    log_photo_error_fn = scipy.interpolate.interp1d(
        delta_arr, log_err_arr, bounds_error=False, fill_value=1.0
    )

    survey = Survey(
        name=name,
        release=release,
        bands=["g", "r"],
        maglim_maps={"g": maglim_map.copy(), "r": maglim_map.copy()},
        coeff_extinc={"g": 3.303, "r": 2.285},
        saturation={"g": 16.0, "r": 16.0},
        sys_error={"g": 0.005, "r": 0.005},
        ebv_map=ebv_map,
        coverage=(maglim_map > 0).astype(float),
        completeness=completeness_fn,
        completeness_band="r",
        delta_saturation=-10.4,
        log_photo_error=log_photo_error_fn,
    )
    return survey


@pytest.fixture
def mock_survey():
    """Synthetic Survey with uniform maps — no file I/O required."""
    return _make_mock_survey()


@pytest.fixture
def mock_injector(mock_survey):
    """StreamInjector pre-loaded with the mock survey."""
    from streamobs.observed import StreamInjector
    return StreamInjector(mock_survey)



# ---------------------------------------------------------------------------
# Magnitudes samples
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def saturation_magnitudes():
    """Array of magnitudes spanning the saturation threshold."""
    return np.linspace(10.0, 17.0, 10)

@pytest.fixture(scope="session")
def bright_magnitudes(rng):
    """Array of magnitudes well above the saturation threshold."""
    return rng.uniform(18.0, 21.0, 100)

@pytest.fixture(scope="session")
def faint_magnitudes(rng):
    """Array of magnitudes well below the saturation threshold."""
    return rng.uniform(27.0, 30.0, 100)



