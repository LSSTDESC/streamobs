"""
Persistence layer for precomputed background CMD histogram grids.
"""

import os

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


class BackgroundStorage:
    """
    Save and load precomputed color–magnitude diagram (CMD) histogram grids.

    One compressed parquet file per ``(source_type, bands)`` combination,
    e.g. ``stars_gr.parquet``.  Inside that file there is one row per
    ``(maglim_b2, maglim_b1)`` grid point, where ``b1 = bands[0]`` (color
    band) and ``b2 = bands[1]`` (reference/magnitude band).

    **File format** — one row per ``(maglim_b2, maglim_b1)`` pair::

        maglim_b2 | maglim_b1 | n_ref | area_ref_deg2
        | color_edge_min | color_edge_max | n_color
        | mag_edge_min   | mag_edge_max   | n_mag
        | counts  (list of n_color × n_mag floats, row-major)

    Bin edges are derived from ``(edge_min, edge_max, n_bins)`` on load.
    Bin centers are not stored; compute them from edges when needed.

    When reading for a specific ``(maglim_b2, maglim_b1)`` pair, pyarrow
    predicate pushdown is used so that only the relevant row groups are
    read from disk — the full file is never loaded into memory.

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
        Build the file path for a ``(source_type, bands)`` combination.

        Parameters
        ----------
        source_type : str
            ``'stars'`` or ``'galaxies'``.
        bands : tuple of str
            Band names in order, e.g. ``('g', 'r')``.

        Returns
        -------
        str
            Absolute path, e.g. ``{base_path}/lsst/stars_gr.parquet``.
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
        Persist the full CMD histogram grid to a single parquet file.

        Parameters
        ----------
        data : dict
            Full grid keyed by ``(maglim_b2, maglim_b1)``, each value being
            a dict with keys ``cmd_hist``, ``color_edges``, ``mag_edges``,
            ``n_ref``, ``area_ref_deg2``.  ``b1 = bands[0]``, ``b2 = bands[1]``.
        source_type : str
            ``'stars'`` or ``'galaxies'``.
        bands : tuple of str
            Band names, e.g. ``('g', 'r')``.
        **kwargs
            compression : str, optional
                Parquet compression codec. Default ``'zstd'``.
        """
        compression = kwargs.get("compression", "zstd")

        rows = []
        for (maglim_b2, maglim_b1), d in data.items():
            color_edges = d["color_edges"]
            mag_edges = d["mag_edges"]
            rows.append(
                {
                    "maglim_b2": round(float(maglim_b2), 4),
                    "maglim_b1": round(float(maglim_b1), 4),
                    "n_ref": int(d["n_ref"]),
                    "area_ref_deg2": float(d["area_ref_deg2"]),
                    "color_edge_min": float(color_edges[0]),
                    "color_edge_max": float(color_edges[-1]),
                    "n_color": int(len(color_edges) - 1),
                    "mag_edge_min": float(mag_edges[0]),
                    "mag_edge_max": float(mag_edges[-1]),
                    "n_mag": int(len(mag_edges) - 1),
                    "counts": d["cmd_hist"].ravel().tolist(),
                }
            )

        df = pd.DataFrame(rows)
        path = self.get_path(source_type, bands)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if os.path.exists(path):
            os.remove(path)
        df.to_parquet(path, compression=compression, index=False)

    def load_data(
        self,
        source_type: str,
        bands: tuple,
        maglim_b2: float,
        maglim_b1: float,
    ) -> dict:
        """
        Load the CMD histogram for a specific ``(maglim_b2, maglim_b1)`` pair.

        Uses pyarrow predicate pushdown — only the relevant row groups are
        read from disk.

        Parameters
        ----------
        source_type : str
            ``'stars'`` or ``'galaxies'``.
        bands : tuple of str
            Band names, e.g. ``('g', 'r')``.
        maglim_b2 : float
            Magnitude limit for ``bands[1]`` (reference band).
        maglim_b1 : float
            Magnitude limit for ``bands[0]`` (color band).

        Returns
        -------
        dict
            ``{'cmd_hist', 'color_edges', 'mag_edges', 'n_ref', 'area_ref_deg2'}``.
        """
        row = self._load_table(source_type, bands, maglim_b2, maglim_b1).to_pandas().iloc[0]
        return self._row_to_dict(row)

    def load_all(self, source_type: str, bands: tuple) -> dict:
        """
        Load the full CMD histogram grid from the parquet file.

        Returns
        -------
        dict
            ``{(maglim_b2, maglim_b1): {'cmd_hist', 'color_edges', 'mag_edges',
            'n_ref', 'area_ref_deg2'}}`` where ``b1 = bands[0]``, ``b2 = bands[1]``.
        """
        df = self._load_table(source_type, bands).to_pandas()
        return {
            (row["maglim_b2"], row["maglim_b1"]): self._row_to_dict(row)
            for _, row in df.iterrows()
        }

    def exists(self, source_type: str, bands: tuple) -> bool:
        """
        Check whether the resource file for this ``(source_type, bands)``
        combination exists on disk.

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

    def _load_table(
        self,
        source_type: str,
        bands: tuple,
        maglim_b2: float = None,
        maglim_b1: float = None,
    ):
        """Read the parquet file, optionally filtering to a single row.

        When ``maglim_b2`` and ``maglim_b1`` are given, pyarrow predicate
        pushdown is applied so only the matching row groups are read.
        When both are ``None``, the full file is returned.
        """
        path = self.get_path(source_type, bands)
        filters = None
        if maglim_b2 is not None and maglim_b1 is not None:
            filters = [
                ("maglim_b2", "=", round(float(maglim_b2), 4)),
                ("maglim_b1", "=", round(float(maglim_b1), 4)),
            ]
        return pq.read_table(path, filters=filters)

    @staticmethod
    def _row_to_dict(row) -> dict:
        """Reconstruct a histogram dict from a single DataFrame row."""
        n_color = int(row["n_color"])
        n_mag = int(row["n_mag"])
        color_edges = np.linspace(row["color_edge_min"], row["color_edge_max"], n_color + 1)
        mag_edges = np.linspace(row["mag_edge_min"], row["mag_edge_max"], n_mag + 1)
        counts = row["counts"]
        if hasattr(counts, "as_py"):
            counts = counts.as_py()
        return {
            "cmd_hist": np.array(counts).reshape(n_color, n_mag),
            "color_edges": color_edges,
            "mag_edges": mag_edges,
            "n_ref": int(row["n_ref"]),
            "area_ref_deg2": float(row["area_ref_deg2"]),
        }
