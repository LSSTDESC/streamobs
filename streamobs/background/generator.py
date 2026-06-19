"""
Fast per-pixel background generation from precomputed CMD grids.
"""

import numpy as np

from ..surveys import Survey
from .storage import BackgroundStorage


class LightBackgroundGenerator:
    """
    Generate background catalogs rapidly from precomputed CMD histogram grids.

    For each HEALPix pixel in the requested sky region the generator:

    1. Retrieves the local magnitude limits from the survey's HEALPix maps.
    2. Selects the nearest precomputed CMD histogram in ``(maglim_ref, delta)``
       space.
    3. Scales the reference object count to the pixel area via a Poisson draw.
    4. Samples ``(color, mag)`` pairs from the 2-D histogram.
    5. Samples sky positions uniformly within the pixel.

    Stars and galaxies are handled independently so that users can build
    custom mixed background models.

    Parameters
    ----------
    storage : BackgroundStorage
        Storage backend holding the precomputed CMD grids.
    survey : Survey
        Survey with real (non-uniform) HEALPix magnitude limit maps used to
        look up per-pixel depth.
    bands : tuple of str, optional
        ``(band1, band2)`` — must match the bands used when building resources.
        Default ``('g', 'r')``.
    **kwargs
        Reserved for future use.

    Examples
    --------
    >>> gen = LightBackgroundGenerator(storage, survey, bands=('g', 'r'))
    >>> catalog = gen.generate(
    ...     phi1_limits=(-20, 20),
    ...     phi2_limits=(-2, 2),
    ...     gc_frame=frame,
    ...     nside=4096,
    ...     source_type='both',
    ... )
    """

    def __init__(
        self,
        storage: BackgroundStorage,
        survey: Survey,
        bands=("g", "r"),
        **kwargs,
    ):
        self.storage = storage
        self.survey = survey
        self.bands = bands
        # Lazy cache: {source_type: grid_dict}
        self._resources: dict = {}

    def generate(
        self,
        phi1_limits,
        phi2_limits,
        gc_frame,
        nside=4096,
        source_type="both",
        **kwargs,
    ) -> "pd.DataFrame":
        """
        Generate a background catalog for the given sky region.

        Parameters
        ----------
        phi1_limits : tuple of float
            ``(phi1_min, phi1_max)`` in degrees.
        phi2_limits : tuple of float
            ``(phi2_min, phi2_max)`` in degrees.
        gc_frame : gala.coordinates.GreatCircleICRSFrame
            Great-circle frame mapping ``(phi1, phi2)`` to ``(RA, Dec)``.
        nside : int, optional
            HEALPix resolution. Default ``4096``.
        source_type : str, optional
            ``'stars'``, ``'galaxies'``, or ``'both'``. Default ``'both'``.
        **kwargs
            rng : numpy.random.Generator, optional
            seed : int, optional

        Returns
        -------
        pd.DataFrame
            Background catalog with columns ``phi1``, ``phi2``, ``color``,
            ``mag`` (reference band), and ``source_type``.
        """
        ...

    def _generate_one_type(
        self,
        phi1_limits,
        phi2_limits,
        gc_frame,
        source_type: str,
        nside: int,
        **kwargs,
    ) -> "pd.DataFrame":
        """
        Generate objects of a single source type pixel by pixel.

        Parameters
        ----------
        phi1_limits, phi2_limits : tuple of float
        gc_frame : gala.coordinates.GreatCircleICRSFrame
        source_type : str
            ``'stars'`` or ``'galaxies'``.
        nside : int
        **kwargs

        Returns
        -------
        pd.DataFrame
        """
        ...

    def _get_footprint_pixels(
        self,
        phi1_limits,
        phi2_limits,
        gc_frame,
        nside: int,
    ) -> np.ndarray:
        """
        Return HEALPix pixel indices that cover the given ``(phi1, phi2)`` box.

        Parameters
        ----------
        phi1_limits, phi2_limits : tuple of float
        gc_frame : gala.coordinates.GreatCircleICRSFrame
        nside : int

        Returns
        -------
        np.ndarray of int
            Pixel indices.
        """
        ...

    def _get_maglim_pixels(
        self,
        pixels: np.ndarray,
        band: str,
        nside: int,
    ) -> np.ndarray:
        """
        Retrieve magnitude limits for a set of pixels from the survey map.

        Parameters
        ----------
        pixels : np.ndarray of int
        band : str
        nside : int
            Resolution of ``pixels``; the survey map may have a different nside
            and will be re-pixelized if needed.

        Returns
        -------
        np.ndarray of float
        """
        ...

    def _select_cmd_distribution(
        self,
        maglim_ref: float,
        delta: float,
        source_type: str,
    ) -> dict:
        """
        Nearest-neighbour lookup in the precomputed CMD grid.

        Parameters
        ----------
        maglim_ref : float
            Reference band magnitude limit for this pixel.
        delta : float
            ``maglim_band1 - maglim_ref`` for this pixel.
        source_type : str
            ``'stars'`` or ``'galaxies'``.

        Returns
        -------
        dict
            ``{'cmd_hist', 'color_edges', 'mag_edges', 'n_ref',
            'area_ref_deg2'}`` for the nearest grid point.
        """
        ...

    def _scale_n_objects(
        self,
        n_ref: int,
        area_ref_deg2: float,
        pixel_area_deg2: float,
    ) -> int:
        """
        Draw the expected number of objects for a pixel via Poisson scaling.

        Parameters
        ----------
        n_ref : int
            Number of objects in the reference simulation.
        area_ref_deg2 : float
            Sky area of the reference simulation in deg².
        pixel_area_deg2 : float
            Sky area of the target pixel in deg².

        Returns
        -------
        int
            Poisson draw around ``n_ref * pixel_area_deg2 / area_ref_deg2``.
        """
        ...

    def _sample_from_cmd(
        self,
        cmd_hist: np.ndarray,
        color_edges: np.ndarray,
        mag_edges: np.ndarray,
        n_objects: int,
        rng: np.random.Generator,
    ) -> "pd.DataFrame":
        """
        Draw ``n_objects`` (color, mag) pairs from the 2-D CMD histogram.

        Parameters
        ----------
        cmd_hist : np.ndarray
            2-D array of counts with shape ``(n_color, n_mag)``.
        color_edges : np.ndarray
        mag_edges : np.ndarray
        n_objects : int
        rng : numpy.random.Generator

        Returns
        -------
        pd.DataFrame
            Columns ``color`` and ``mag``.
        """
        ...

    def _sample_positions(
        self,
        n_objects: int,
        pixel: int,
        nside: int,
        gc_frame,
        rng: np.random.Generator,
    ) -> "pd.DataFrame":
        """
        Sample ``(phi1, phi2)`` positions uniformly within a HEALPix pixel.

        Parameters
        ----------
        n_objects : int
        pixel : int
            HEALPix pixel index.
        nside : int
        gc_frame : gala.coordinates.GreatCircleICRSFrame
            Used to convert ``(RA, Dec)`` corners back to ``(phi1, phi2)``.
        rng : numpy.random.Generator

        Returns
        -------
        pd.DataFrame
            Columns ``phi1`` and ``phi2``.
        """
        ...
