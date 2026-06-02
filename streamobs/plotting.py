#!/usr/bin/env python
"""
Plotting functions.
"""

import numpy as np
import pylab as plt
import os
import pandas as pd

try:
    import healpy as hp
    import skyproj
except ImportError:
    healpy = None
    skyproj = None

def draw_stream(phi1, phi2):
    """ Create 2d histogram and draw stellar distribution.

    Parameters
    ----------
    phi1, phi2 : coordinates of stars (deg)

    Returns:
    --------
    num, xedges, yedges, im : output of imshow
    """

    phi1_bins = np.linspace(phi1.min(), phi1.max(), 100)
    phi2_bins = np.linspace(phi2.min(), phi2.max(), 100)

    hist, xedges, yedges = np.histogram2d(phi1, phi2,
                                           bins=[phi1_bins, phi2_bins],
    )

    ax = plt.gca()
    ax.set_title("Mock Stream")
    ax.set_xlabel("$\phi_1$ [deg]")
    ax.set_ylabel("$\phi_2$ [deg]")

    image = ax.imshow(
        hist.T,
        extent=[xedges.min(), xedges.max(), yedges.min(), yedges.max()],
        aspect="auto",
        vmin=np.percentile(hist, 0.1),
        vmax=np.percentile(hist, 99.9),
        origin="lower",
        cmap="gray_r",
    )
    return hist, xedges, yedges, image

def plot_stream(phi1, phi2):
    """Plot binned histogram of stream.

    Parameters
    ----------
    phi1, phi2 : coordinates of stars (deg)

    Returns
    -------
    fig, ax : the figure and axis
    """
    fig = plt.figure(dpi=175, figsize=(10, 2.5))
    gs = fig.add_gridspec(1, 1, wspace=0, hspace=0)
    plt.subplots_adjust(bottom=0.2)

    ax = fig.add_subplot(gs[0, 0])
    hist, xedges, yedges, image = draw_stream(phi1, phi2)
    return fig, ax


###############################################################################
####################### Injection plotting functions ##########################
###############################################################################


def plot_inject(data, survey, bands=None,range=None, **kwargs):
    """
    Plot injection results showing observed vs unobserved stars.
    
    Creates three panels:
    - Left: Injected stars in the survey footprint (using magnitude limit map)
    - Center: HR diagram using true magnitudes
    - Right: HR diagram using measured magnitudes with photometric errors
    
    Parameters
    ----------
    data : pandas.DataFrame or dict-like
        Data containing the injected stream. Must have columns:
        - 'ra', 'dec': Sky coordinates in degrees
        - 'flag_observed': Boolean flag for detected stars
        - 'mag_<band>': True magnitudes for each band
        - 'mag_<band>_obs': Observed magnitudes for each band
    survey : Survey
        Survey object containing magnitude limit maps and other properties.
    bands : list of str, optional
        Bands to use for HR diagram. Default is ['g', 'r'].
        Will use first two bands if more are provided.
    **kwargs : dict
        Additional keyword arguments:
        - save : bool, if True saves the figure
        - folder : str, path to save folder (default: ../data/outputs/)
        - survey_name : str, survey name for filename (default: uses survey object)
    
    Returns
    -------
    fig, ax : tuple
        The matplotlib figure and axes array.
    
    Examples
    --------
    >>> fig, ax = plot_inject(observed_data, survey, bands=['g', 'r'], save=True)
    """
    # Convert to DataFrame if needed
    if not isinstance(data, pd.DataFrame):
        data = pd.DataFrame(data)
    
    # Default to g and r bands
    if bands is None:
        bands = ['g', 'r']
    
    # Use first two bands for color-magnitude diagram
    if len(bands) < 2:
        raise ValueError("Need at least 2 bands for HR diagram")
    band1, band2 = bands[0], bands[1]
    
    # Check required columns
    required_cols = ['ra', 'dec', 'flag_observed', 
                     f'mag_{band1}', f'mag_{band2}',
                     f'mag_{band1}_obs', f'mag_{band2}_obs']
    missing_cols = [col for col in required_cols if col not in data.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Get detection flags
    sel = data["flag_observed"].astype(bool)
    
    # Create figure
    fig, ax = plt.subplots(1, 3, figsize=(14, 6))
    plt.subplots_adjust(left=0.075, right=0.98)
    
    # --- Panel 1: Spatial distribution with magnitude limit map ---
    ax[0].set_title("Injection of the stream into the footprint")
    
    # Get magnitude limit map for the first band
    try:
        maglim_map = survey.maglim_maps[band1]
        nside = hp.get_nside(maglim_map)

        if range is not None:
            ra_min, ra_max = range[0], range[1]
            dec_min, dec_max = range[2], range[3]
        else:
            # Determine RA/Dec range from data with some padding
            ra_min, ra_max = data['ra'].min() - 5, data['ra'].max() + 5
            dec_min, dec_max = data['dec'].min() - 5, data['dec'].max() + 5
        
        # Create grid for the map
        x = np.linspace(ra_min, ra_max, 200)
        y = np.linspace(dec_min, dec_max, 200)
        XX, YY = np.meshgrid(x, y)
        
        # Get magnitude limits at each grid point
        ZZ = maglim_map[hp.ang2pix(nside, XX, YY, lonlat=True)]
        ZZ[ZZ < 0] = np.nan  # Mask invalid pixels
        
        # Plot magnitude limit map
        # PF: I would consider using skyproj for this spatial plot to get a proper sky projection
        mesh = ax[0].pcolormesh(XX, YY, ZZ, shading="auto", cmap='viridis')
        fig.colorbar(mesh, ax=ax[0], label=f"Mag Limit ({band1})")
        
    except (KeyError, AttributeError) as e:
        # If map not available, just show a simple plot
        print(f"Warning: Could not plot magnitude limit map: {e}")
    
    # Plot stars
    ax[0].scatter(
        data["ra"], data["dec"], s=2, alpha=0.5, color="gray", label="Unobserved"
    )
    ax[0].scatter(
        data["ra"][sel],
        data["dec"][sel],
        s=4,
        alpha=1.0,
        color="red",
        label="Observed",
    )
    ax[0].invert_xaxis()
    ax[0].set_xlabel("RA (deg)")
    ax[0].set_ylabel("Dec (deg)")
    ax[0].legend()
    
    # --- Panel 2: HR diagram with true magnitudes ---
    ax[1].set_title("HR diagram using True magnitudes")
    
    color_true = data[f"mag_{band1}"] - data[f"mag_{band2}"]
    mag_true = data[f"mag_{band1}"]
    
    ax[1].scatter(
        color_true,
        mag_true,
        s=2,
        alpha=0.5,
        color="gray",
        label="Unobserved",
    )
    ax[1].scatter(
        color_true[sel],
        mag_true[sel],
        s=4,
        alpha=1.0,
        color="red",
        label="Observed",
    )
    ax[1].set_xlim(-0.5, 1.5)
    ax[1].set_ylim(30, 16)
    ax[1].set_xlabel(f"({band1}-{band2})")
    ax[1].set_ylabel(band1)
    ax[1].legend()
    
    # --- Panel 3: HR diagram with measured magnitudes ---
    ax[2].set_title("HR diagram using sampled observed magnitudes")
    
    # Convert measured magnitudes to numeric, handling "BAD_MAG" strings
    mag1_obs = pd.to_numeric(data[f"mag_{band1}_obs"], errors="coerce")
    mag2_obs = pd.to_numeric(data[f"mag_{band2}_obs"], errors="coerce")
    
    # Mask out bad measurements
    mask_good = (~mag1_obs.isna()) & (~mag2_obs.isna())
    
    if mask_good.sum() > 0:
        color_obs = mag1_obs - mag2_obs
        
        ax[2].scatter(
            color_obs[mask_good],
            mag1_obs[mask_good],
            s=2,
            alpha=0.5,
            color="gray",
            label="Unobserved",
        )
        ax[2].scatter(
            color_obs[sel & mask_good],
            mag1_obs[sel & mask_good],
            s=4,
            alpha=1.0,
            color="red",
            label="Observed",
        )
    else:
        print("Warning: No valid measured magnitudes to plot")
    
    ax[2].set_xlim(-0.5, 1.5)
    ax[2].set_ylim(30, 16)
    ax[2].set_xlabel(f"({band1}-{band2})")
    ax[2].set_ylabel(band1)
    ax[2].legend()
    
    # --- Save figure if requested ---
    if kwargs.pop("save", False):
        current_dir = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(current_dir, "..", "data", "outputs")
        folder = kwargs.pop("folder", path)
        
        if not os.path.exists(folder):
            os.makedirs(folder)
        
        # Get survey name for filename
        survey_name = kwargs.pop("survey_name", getattr(survey, 'name', 'unknown'))
        
        # Find next available file number
        existing_files = os.listdir(folder)
        num_fichier = len(
            [f for f in existing_files if f.startswith(f"injection_{survey_name}_")]
        )
        file_name = f"injection_{survey_name}_{num_fichier:03d}"
        filename = os.path.join(folder, file_name + ".png")
        
        print(f"Saving figure: {filename}...")
        plt.savefig(filename, dpi=150, bbox_inches='tight')
    
    return fig, ax


def plot_stream_in_mask(ra,dec,mask, nest=False, output_folder = None):
    """Plot stream in mask using healpy and skyproj.

    Parameters
    ----------
    ra, dec : coordinates of stars (deg)
    mask : boolean array, True for masked points

    Returns
    -------
    fig, ax : the figure and axis
    """
    if skyproj is None:
        raise ImportError("healpy or skyproj not installed, cannot plot stream in mask.")
    
    if mask is None:
        mask = np.ones(hp.nside2npix(32), dtype=bool)  # Default mask if none provided

    mask = mask.astype(float)
    fig, ax = plt.subplots(figsize=(8, 5))
    sp = skyproj.McBrydeSkyproj(ax=ax)
    sp.draw_hpxmap(mask, nest=nest)
    sp.ax.scatter(ra, dec, color='red', s=10, label='Stream', alpha=0.5)
    sp.ax.legend(loc='lower right')
    plt.colorbar()
    fig.tight_layout()

    if output_folder is not None:
        os.makedirs(output_folder, exist_ok=True)
        filename = os.path.join(output_folder, 'stream_in_mask.png')
        if filename is not None:
            print(f"Writing figure: {filename}...")
            plt.savefig(filename)
        
    return fig, ax


def plot_healsparse_map(map):
    if isinstance(map, np.ndarray): # healpix map
        healpix_map = map
    elif isinstance(map, str): # filename
        healpix_map = hp.read_map(map, verbose=False)
    elif 'healsparse' in globals():
       import healsparse as hsp
       if isinstance(map, hsp.HealSparseMap):
           healsparse_map = map
           nside_sparse = healsparse_map.nside_sparse 
           nest = False  
           healpix_map = healsparse_map.generate_healpix_map(nside=nside_sparse,nest=nest)
    else:
        raise ValueError("map must be a healpix map array, filename, or HealSparseMap object")
         
    fig, ax = plt.subplots(figsize=(8, 5))
    sp = skyproj.McBrydeSkyproj(ax=ax)
    sp.draw_hpxmap(healpix_map, nest=False)
    plt.colorbar()   
    fig.tight_layout()

    return fig, ax


def plot_maglim(maglim_map):
    fig, ax = plot_healsparse_map(maglim_map)
    return fig, ax