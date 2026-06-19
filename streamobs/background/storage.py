"""
Persistence layer for precomputed background CMD histogram grids.
"""

import os


class BackgroundStorage:
    """
    Save and load precomputed color–magnitude diagram (CMD) histogram grids.

    One compressed parquet file is stored per ``(source_type, bands)``
    combination, e.g. ``stars_gr.parquet``.  All files live under
    ``data/background/`` (which is excluded from version control).

    **File format** — long-format parquet with one row per
    ``(maglim_r, maglim_g, color_center, mag_center)`` cell::

        maglim_r | maglim_g | color_center | mag_center | count | n_ref | area_ref_deg2

    Parquet's columnar compression handles repeated ``color_center`` and
    ``mag_center`` values efficiently.

    Parameters
    ----------
    base_path : str, optional
        Root directory for resource files.  Defaults to
        ``{package_root}/data/background/``.
    survey_name : str, optional
        Survey identifier (e.g. ``'lsst'``).  Used in the file path.

    Examples
    --------
    >>> storage = BackgroundStorage(survey_name='lsst')
    >>> path = storage.get_path('stars', ('g', 'r'))
    >>> # -> .../data/background/lsst/stars_gr.parquet
    """

    def __init__(self, base_path=None, survey_name="lsst", **kwargs):
        if base_path is None:
            _pkg_root = os.path.join(os.path.dirname(__file__), "..", "..", "data")
            base_path = os.path.join(_pkg_root, "background")
        self.base_path = base_path
        self.survey_name = survey_name

    def get_path(self, source_type: str, bands: tuple) -> str:
        """
        Build the file path for a given ``(source_type, bands)`` combination.

        Parameters
        ----------
        source_type : str
            ``'stars'`` or ``'galaxies'``.
        bands : tuple of str
            Band names in order, e.g. ``('g', 'r')``.

        Returns
        -------
        str
            Absolute path to the parquet file, e.g.
            ``{base_path}/lsst/stars_gr.parquet``.
        """
        bands_str = "".join(bands)
        filename = f"{source_type}_{bands_str}.parquet"
        return os.path.join(self.base_path, self.survey_name, filename)

    def save_data(
        self,
        data: dict,
        source_type: str,
        bands: tuple,
        **kwargs,
    ):
        """
        Persist the CMD histogram grid to a compressed parquet file.

        The ``data`` dict is keyed by ``(maglim_r, maglim_g)`` and each value
        is a dict with keys ``cmd_hist``, ``color_edges``, ``mag_edges``,
        ``n_ref``, and ``area_ref_deg2`` (as returned by
        :meth:`~streamobs.background.resource_builder.BackgroundResourceBuilder._build_one_config`).
        The grid is flattened to long-format before writing.

        Parameters
        ----------
        data : dict
            CMD histogram grid keyed by ``(maglim_r, maglim_g)``.
        source_type : str
            ``'stars'`` or ``'galaxies'``.
        bands : tuple of str
            Band names, e.g. ``('g', 'r')``.
        **kwargs
            compression : str, optional
                Parquet compression codec. Default is ``'zstd'``.
        """
        ...

    def load_data(
        self,
        source_type: str,
        bands: tuple,
        **kwargs,
    ) -> dict:
        """
        Load a CMD histogram grid from the parquet file for this combination.

        Reads the file returned by :meth:`get_path` and reconstructs the
        nested dict keyed by ``(maglim_r, maglim_g)``.

        Parameters
        ----------
        source_type : str
            ``'stars'`` or ``'galaxies'``.
        bands : tuple of str
            Band names, e.g. ``('g', 'r')``.

        Returns
        -------
        dict
            CMD grid keyed by ``(maglim_r, maglim_g)`` → ``{'cmd_hist',
            'color_edges', 'mag_edges', 'n_ref', 'area_ref_deg2'}``.
        """
        ...

    def exists(self, source_type: str, bands: tuple) -> bool:
        """
        Check whether the resource file for this combination exists on disk.

        Parameters
        ----------
        source_type : str
            ``'stars'`` or ``'galaxies'``.
        bands : tuple of str
            Band names, e.g. ``('g', 'r')``.

        Returns
        -------
        bool
        """
        return os.path.exists(self.get_path(source_type, bands))
