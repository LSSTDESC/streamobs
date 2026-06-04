# Quickstart Guide

This page provides a quick introduction to get you started with StreamObs.

First, make sure you have StreamObs installed (see [Installation](installation.md) for details).


## Generating Mock Stream Data

StreamObs can generate parametric stellar stream models using configuration files or directly in Python.

### Example 1: Using Command-Line Tools

The simplest way to generate a stream is using the provided command-line scripts:

```bash
# Generate a simple linear stream
./bin/generate_stream.py config/toy1_config.yaml -o toy1_stream.csv --plot

# Generate a sinusoidal stream
./bin/generate_stream.py config/toy2_config.yaml -o toy2_stream.csv --plot

# Generate Pal 5 stream from interpolation
./bin/generate_stream.py config/pal5_config.yaml -o pal5_stream.csv --plot
```

These commands will:
- Read the stream parameters from the configuration file
- Generate mock stellar positions (phi1, phi2) and magnitudes
- Save the results to a CSV file
- Create visualization plots (with `--plot` flag)

### Example 2: Generating Streams in Python

You can also generate streams programmatically:

```python
import numpy as np
import pandas as pd
from streamobs.model import StreamModel
from streamobs.utils import parse_config

# Load configuration
config = parse_config('config/toy1_config.yaml')

# Create stream model and generate stars
stream = StreamModel(config['stream'])
stream_df = stream.sample(config['stream']['nstars'])

# The dataframe contains: phi1, phi2, distance, magnitudes, etc.
print(stream_df.head())
```

## Converting to Observable Quantities

Once you have mock stream data (either generated or from simulations), you can convert it to realistic observations using the `StreamInjector` class.

### Example 3: Applying Survey Effects

```python
import numpy as np
import pandas as pd
from streamobs import surveys, observed

# Load a survey (e.g., LSST Year 1)
lsst_survey = surveys.Survey.load(survey='lsst', release='yr1')

# Create the stream injector
injector = observed.StreamInjector(lsst_survey)

# Create or load your mock stream data
# Here we create a simple test dataset
# This could contain (ra, dec) coordinates if you want to skip coordinate transformation
rng = np.random.default_rng(42)
mock_data = pd.DataFrame({
    'phi1': rng.uniform(-5, 5, 1000),      # Stream longitude
    'phi2': rng.uniform(-1, 1, 1000),      # Stream latitude
    'mag_g': rng.uniform(18, 28, 1000),    # g-band magnitude
    'mag_r': rng.uniform(18, 28, 1000),    # r-band magnitude
})

# Apply survey effects: footprint, extinction, photometric errors
observed_data = injector.inject(
    mock_data, 
    seed=42,
    bands = ['r', 'g'] # Choose the bands to use
    mask_type=['footprint', 'ebv'],  # Restrict injection to a mask created from the footprind and low dust area
    verbose=True
)

print(f"Input stars: {len(mock_data)}")
print(f"Detected stars: {len(observed_data[observed_data['flag_observed']==1])}")
```

### What the Injector Does

The `StreamInjector` applies several observational effects:

1. **Coordinate conversion**: Converts stream coordinates (phi1, phi2) to sky coordinates (RA, Dec)
2. **Extinction**: Applies Galactic dust extinction corrections
3. **Photometric errors**: Adds realistic magnitude uncertainties
4. **Detection completeness**: Applies magnitude-dependent detection probability


The output dataframe includes:
- `ra`, `dec`: Sky coordinates
- `mag_g_obs`, `mag_r_obs`: Observed magnitudes with errors
- `magerr_g`, `magerr_r`: Photometric uncertainties
- `flag_observed`: Detection and clasification flag (1=detected & classified as a star, 0=not detected of not classified as a star)


## Next Steps

Now that you've seen the basics, you can:

- **Learn more about configuration files**: See the `config/` directory for examples of stream and survey configurations
- **Explore the API**: Check the [API Reference](modules.rst) for detailed documentation of all classes and functions
- **View detailed examples**: Visit the `notebooks/` directory for Jupyter notebooks with more complex workflows
- **Customize surveys**: Learn how to define your own survey parameters in YAML files
- **Try different stream models**: Experiment with different parametric stream descriptions

## Key Concepts

- **Mock generation**: Create idealized stream data with known properties
- **Survey injection**: Apply realistic observational effects to mock data
- **Configuration-driven**: Use YAML files to define stream and survey parameters
- **Modular design**: Use individual components (generation, photometry, observation) independently

## Getting Help

- **Documentation**: Browse the full documentation for detailed information
- **Examples**: Check the `notebooks/` directory for worked examples
- **Issues**: Report bugs or request features on [GitHub](https://github.com/LSSTDESC/streamobs/issues)
