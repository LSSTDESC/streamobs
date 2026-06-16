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

@pytest.fixture(scope="session")
def verbose():
    """Control verbosity of test output."""
    return True

# ---------------------------------------------------------------------------
# Minimal stream model config (no isochrone, no velocity — pure geometry)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def minimal_stream_config():
    """Minimal StreamModel config: density + track only (no isochrone / velocity)."""
    return {
        'density': {'type': 'Uniform', 'xmin': -9.0, 'xmax': 9.0}, 

        # Track model
        'track': {'center': {'type': 'Constant', 'value': 0.0}, # center line of the stream in degrees
                    'spread': {'type': 'Constant', 'value': 0.2}, # spread of the stream in degrees
                    'sampler': 'Gaussian'}, # how to sample across the stream

        # Isochrone model
        'isochrone': {'name': 'Marigo2017', # isochrone set name
                        'survey': 'lsst', # survey for filter set
                        'age': 12.0, # Age in Gyr of the population
                        'z': 0.0006, # Metallicity of the population
                        'band_1': 'g', # first band for color-magnitude
                        'band_2': 'r', # second band for color-magnitude
                        'band_1_detection': True}, 
    }

@pytest.fixture(scope="session")
def stream_config_with_distance(minimal_stream_config):
    """StreamModel config that also produces a distance modulus column."""
    cfg = dict(minimal_stream_config)
    cfg["distance_modulus"] ={'center':  {'type': 'Constant', 'value': 16.8}, 
                                 'spread': {'type': 'Constant', 'value': 0.0}, 
                                }
    return cfg

# ---------------------------------------------------------------------------
# Tiny sample DataFrame used for injection tests
# ---------------------------------------------------------------------------

@pytest.fixture
def sample_catalog_radec(rng):
    """
    Small catalog with ra, dec, mag_g, mag_r columns — ready for injection
    without needing coordinate conversion or a stream model.
    """
    n = 50
    return pd.DataFrame({
        "ra":    rng.uniform(30.0, 60.0, n),
        "dec":   rng.uniform(-20.0, 0.0, n),
    })


@pytest.fixture
def sample_catalog_phi(rng):
    """Small catalog with phi1/phi2 stream coordinates (no ra/dec yet)."""
    n = 50
    phi1 = rng.uniform(-8.0, 8.0, n)
    phi2 = rng.normal(0.0, 0.1, n)
    return pd.DataFrame({
        "phi1":  phi1,
        "phi2":  phi2,
    })


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

# ---------------------------------------------------------------------------
# Default survey config for injection tests
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def base_maglim():
    """Base magnitude limit for the mock survey."""
    return 26

@pytest.fixture(scope="session")
def mock_survey(verbose):
    from streamobs import surveys
    return surveys.Survey.load(
        survey="lsst",
        release="yr4",
        verbose=verbose,
    )

@pytest.fixture(scope="session")
def mock_injector(mock_survey, verbose):
    from streamobs.observed import StreamInjector
    return StreamInjector(survey=mock_survey, verbose=verbose)


@pytest.fixture(scope="session")
def stream_catalog(stream_config_with_distance):
    from streamobs.model import StreamModel
    model = StreamModel(stream_config_with_distance)
    samples = model.sample(1000)
    return samples