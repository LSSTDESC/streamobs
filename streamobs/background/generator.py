"""
Fast per-pixel background generation from precomputed CMD grids.
"""

import warnings

import astropy.coordinates as coord
import astropy.units as u
import healpy as hp
import numpy as np
import pandas as pd

from ..columns import obs_col
from ..surveys import Survey
from .storage import BackgroundStorage


class LightBackgroundGenerator:
    """
    Generate background catalogs rapidly from precomputed CMD histogram grids.

    For each HEALPix pixel in the requested sky region the generator:

    1. Retrieves the effective magnitude limit (observed maglim minus dust extinction).
    2. Bilinearly interpolates the CMD histogram grid at that effective maglim
       (taking into account the dust absorption in each band).
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
        ``(band_g, band_r)`` — must match the bands used when building resources.
        Default ``('g', 'r')``.
    **kwargs
        Reserved for future use.

    Examples
    --------
    >>> gen = LightBackgroundGenerator(storage, survey, bands=('g', 'r'))
    >>> catalog, meta = gen.generate(
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
    ):
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
        tuple of (pd.DataFrame, dict)
            Background catalog with columns ``phi1``, ``phi2``,
            ``mag_{band_r}`` (e.g. ``mag_r``), ``mag_{band_g}`` (e.g. ``mag_g``),
            and ``source_type``; plus a metadata dict with keys ``nside``,
            ``color_edges``, ``mag_edges``, ``band1``, ``band2``.
        """
        rng = kwargs.get("rng")
        if rng is None:
            rng = np.random.default_rng(kwargs.get("seed"))

        # Cap nside to the resolution of the magnitude-limit maps.
        # There's no benefit to sampling positions at sub-pixel resolution relative to the depth maps.
        maglim_nsides = [
            hp.get_nside(self.survey.maglim_maps[b])
            for b in self.bands
            if self.survey.maglim_maps.get(b) is not None
        ]
        if maglim_nsides:
            maglim_nside = min(maglim_nsides)
            if nside > maglim_nside:
                warnings.warn(
                    f"Requested nside={nside} exceeds the magnitude-limit map "
                    f"resolution (nside={maglim_nside}). "
                    f"Capping to nside={maglim_nside}.",
                    UserWarning,
                    stacklevel=2,
                )
                nside = maglim_nside

        active = ["stars", "galaxies"] if source_type == "both" else [source_type]

        catalogs = []
        for st in active:
            df = self._generate_one_type(
                phi1_limits, phi2_limits, gc_frame, st, nside, rng
            )
            if len(df) > 0:
                catalogs.append(df)

        namespace = self.survey.namespace
        b1, b2 = self.bands[0], self.bands[1]
        col1, col2 = obs_col(b1, namespace), obs_col(b2, namespace)

        catalog = (
            pd.concat(catalogs, ignore_index=True)
            if catalogs
            else pd.DataFrame(
                columns=["ra", "dec", "phi1", "phi2", col1, col2, "source_type"]
            )
        )

        meta = {
            "nside": nside,
            "namespace": namespace,
            "band1": b1,
            "band2": b2,
        }
        for st in active:
            if self._resources.get(st):
                first = next(iter(self._resources[st].values()))
                meta["color_edges"] = first["color_edges"]
                meta["mag_edges"] = first["mag_edges"]
                break

        return catalog, meta

    def _load_resources(self, source_type: str):
        """Load CMD grid for *source_type* from storage and cache it."""
        if source_type not in self._resources:
            self._resources[source_type] = self.storage.load_all(
                source_type, self.bands
            )

    def _generate_one_type(
        self,
        phi1_limits,
        phi2_limits,
        gc_frame,
        source_type: str,
        nside: int,
        rng: np.random.Generator,
    ) -> pd.DataFrame:
        """Generate objects of a single source type, pixel by pixel.

        All pixel positions are collected in (ra, dec) space first, then a
        single batched gala coordinate transform converts them to (phi1, phi2).
        This avoids per-pixel transform overhead.
        """
        self._load_resources(source_type)

        pixels = self._get_footprint_pixels(phi1_limits, phi2_limits, gc_frame, nside)
        namespace = self.survey.namespace
        b1, b2 = self.bands[0], self.bands[1]  # b1 = color band, b2 = reference band
        col1, col2 = obs_col(b1, namespace), obs_col(b2, namespace)
        empty = pd.DataFrame(
            columns=["ra", "dec", "phi1", "phi2", col1, col2, "source_type"]
        )
        if len(pixels) == 0:
            return empty

        pixel_area_deg2 = hp.nside2pixarea(nside, degrees=True)
        # b2 is the reference (magnitude) band; b1 is the color (secondary) band.
        # Grid keys are (maglim_b2, maglim_b1), matching how storage was built.
        maglim_b2_eff = self._get_effective_maglim(pixels, b2, nside)
        maglim_b1_eff = self._get_effective_maglim(pixels, b1, nside)

        # Accumulate samples; defer gala transform to a single batch call
        all_ra: list = []
        all_dec: list = []
        all_col2: list = []  # reference-band magnitudes
        all_col1: list = []  # color-band magnitudes

        for i, pixel in enumerate(pixels):
            if maglim_b2_eff[i] is None or maglim_b1_eff[i] is None:
                continue
            if maglim_b2_eff[i] <= 0 or maglim_b1_eff[i] <= 0:
                continue
            if np.isnan(maglim_b2_eff[i]) or np.isnan(maglim_b1_eff[i]):
                continue

            cmd = self._interpolate_cmd(
                float(maglim_b2_eff[i]), float(maglim_b1_eff[i]), source_type
            )
            if cmd is None or cmd["cmd_hist"].sum() == 0:
                continue

            # Scale the pixel's expected object count from the reference area to the pixel area.
            n_objects = self._scale_n_objects(
                cmd["cmd_hist"].sum(), cmd["area_ref_deg2"], pixel_area_deg2, rng
            )
            if n_objects == 0:
                continue

            # Sample (color, mag) pairs from the 2-D CMD histogram and sample positions within the pixel.
            df_cmd = self._sample_from_cmd(
                cmd["cmd_hist"], cmd["color_edges"], cmd["mag_edges"], n_objects, rng
            )
            ra_pix, dec_pix = self._sample_positions(n_objects, pixel, nside, rng)

            all_ra.append(ra_pix)
            all_dec.append(dec_pix)
            all_col2.append(df_cmd["mag"].values)
            all_col1.append((df_cmd["mag"] + df_cmd["color"]).values)

        if not all_ra:
            return empty

        ra_all = np.concatenate(all_ra)
        dec_all = np.concatenate(all_dec)

        # Single batched coordinate transform for all pixels combined
        sky = coord.SkyCoord(ra=ra_all * u.deg, dec=dec_all * u.deg, frame="icrs")
        gc_coords = sky.transform_to(gc_frame)
        phi1_all = gc_coords.phi1.deg
        phi2_all = gc_coords.phi2.deg

        # hp.query_polygon with inclusive=True returns pixels that merely touch
        # the box boundary, so sampled positions can fall outside the limits.
        # Filter to the exact requested region here.
        phi1_min, phi1_max = phi1_limits
        phi2_min, phi2_max = phi2_limits
        in_box = (
            (phi1_all >= phi1_min)
            & (phi1_all <= phi1_max)
            & (phi2_all >= phi2_min)
            & (phi2_all <= phi2_max)
        )

        col2_all = np.concatenate(all_col2)
        col1_all = np.concatenate(all_col1)

        if not in_box.any():
            return empty

        return pd.DataFrame(
            {
                "ra": ra_all[in_box],
                "dec": dec_all[in_box],
                "phi1": phi1_all[in_box],
                "phi2": phi2_all[in_box],
                col2: col2_all[in_box],
                col1: col1_all[in_box],
                "source_type": source_type,
            }
        )

    def _get_footprint_pixels(
        self,
        phi1_limits,
        phi2_limits,
        gc_frame,
        nside: int,
    ) -> np.ndarray:
        """Return HEALPix pixel indices covering the ``(phi1, phi2)`` bounding box."""
        phi1_min, phi1_max = phi1_limits
        phi2_min, phi2_max = phi2_limits

        # 4 corners of the (phi1, phi2) box, listed counterclockwise
        phi1_corners = np.array([phi1_min, phi1_max, phi1_max, phi1_min]) * u.deg
        phi2_corners = np.array([phi2_min, phi2_min, phi2_max, phi2_max]) * u.deg
        stream_coords = coord.SkyCoord(
            phi1=phi1_corners, phi2=phi2_corners, frame=gc_frame
        )
        icrs = stream_coords.icrs
        ra_corners = icrs.ra.deg
        dec_corners = icrs.dec.deg

        # hp.ang2vec returns (N, 3) — shape expected by query_polygon
        vecs = hp.ang2vec(ra_corners, dec_corners, lonlat=True)
        return hp.query_polygon(nside, vecs, inclusive=True)

    def _get_effective_maglim(
        self,
        pixels: np.ndarray,
        band: str,
        nside: int,
    ) -> np.ndarray:
        """
        Compute effective magnitude limit for each pixel.

        ``maglim_eff = maglim_obs - A_band``  (dust reduces the effective depth).
        """
        ra, dec = hp.pix2ang(nside, pixels, lonlat=True)

        # Maglim — re-pixelise to the survey map's nside if different
        maglim_map = self.survey.maglim_maps[band]
        nside_maglim = hp.get_nside(maglim_map)
        pix_maglim = hp.ang2pix(nside_maglim, ra, dec, lonlat=True)
        maglim = maglim_map[pix_maglim].astype(float)

        # Extinction
        if self.survey.ebv_map is not None and band in self.survey.coeff_extinc:
            nside_ebv = hp.get_nside(self.survey.ebv_map)
            pix_ebv = hp.ang2pix(nside_ebv, ra, dec, lonlat=True)
            extinction = self.survey.coeff_extinc[band] * self.survey.ebv_map[pix_ebv]
        else:
            extinction = np.zeros_like(maglim)

        return maglim - extinction

    def _interpolate_cmd(
        self,
        maglim_b2: float,
        maglim_b1: float,
        source_type: str,
    ) -> dict:
        """
        Bilinear interpolation of CMD histogram in ``(maglim_b2, maglim_b1)`` space.

        ``maglim_b2`` is the effective limit for the reference band (``bands[1]``);
        ``maglim_b1`` for the color band (``bands[0]``).  Grid keys follow the
        same convention as :class:`BackgroundStorage`: ``(maglim_b2, maglim_b1)``.

        All grid points share the same bin edges, so per-bin linear weighting
        is valid. The result is clipped to ≥ 0.
        """
        grid = self._resources[source_type]
        if not grid:
            return None

        keys = list(grid.keys())
        r_vals = sorted(set(mr for (mr, _) in keys))
        g_vals = sorted(set(mg for (_, mg) in keys))

        # Clamp to grid range
        maglim_b2 = float(np.clip(maglim_b2, r_vals[0], r_vals[-1]))
        maglim_b1 = float(np.clip(maglim_b1, g_vals[0], g_vals[-1]))

        # Bracket in b2 (first key axis)
        r1_idx = max(0, int(np.searchsorted(r_vals, maglim_b2, side="right")) - 1)
        r2_idx = min(len(r_vals) - 1, r1_idx + 1)
        r1, r2 = r_vals[r1_idx], r_vals[r2_idx]
        wr = (maglim_b2 - r1) / (r2 - r1) if r2 != r1 else 0.0

        # Bracket in b1 (second key axis)
        g1_idx = max(0, int(np.searchsorted(g_vals, maglim_b1, side="right")) - 1)
        g2_idx = min(len(g_vals) - 1, g1_idx + 1)
        g1, g2 = g_vals[g1_idx], g_vals[g2_idx]
        wg = (maglim_b1 - g1) / (g2 - g1) if g2 != g1 else 0.0

        # Bilinear weights for 4 corners
        corners = [
            ((1 - wr) * (1 - wg), (r1, g1)),
            ((1 - wr) * wg, (r1, g2)),
            (wr * (1 - wg), (r2, g1)),
            (wr * wg, (r2, g2)),
        ]

        H_interp = None
        n_ref_interp = 0.0
        area_interp = 0.0
        w_total = 0.0
        first_key = None

        for w, key in corners:
            if key not in grid or w == 0:
                continue
            if H_interp is None:
                H_interp = grid[key]["cmd_hist"] * w
                first_key = key
            else:
                H_interp = H_interp + grid[key]["cmd_hist"] * w
            n_ref_interp += w * grid[key]["n_ref"]
            area_interp += w * grid[key]["area_ref_deg2"]
            w_total += w

        if H_interp is None:
            return None

        # Re-normalise weights if some corners were missing
        if w_total > 0 and abs(w_total - 1.0) > 1e-9:
            H_interp /= w_total
            n_ref_interp /= w_total
            area_interp /= w_total

        return {
            "cmd_hist": np.clip(H_interp, 0, None),
            "color_edges": grid[first_key]["color_edges"],
            "mag_edges": grid[first_key]["mag_edges"],
            "n_ref": n_ref_interp,
            "area_ref_deg2": area_interp,
        }

    def _scale_n_objects(
        self,
        n_detected: float,
        area_ref_deg2: float,
        pixel_area_deg2: float,
        rng: np.random.Generator,
    ) -> int:
        """Poisson draw for the expected object count in one pixel.

        ``n_detected`` is ``cmd_hist.sum()`` — the number of sources actually
        detected in the reference area at this effective maglim, not the number
        injected (``n_ref``).  Using detected counts means the rate naturally
        reflects survey depth: a shallower CMD has fewer counts and therefore
        produces fewer objects per pixel.
        """
        lam = n_detected * pixel_area_deg2 / area_ref_deg2
        return int(rng.poisson(lam))

    def _sample_from_cmd(
        self,
        cmd_hist: np.ndarray,
        color_edges: np.ndarray,
        mag_edges: np.ndarray,
        n_objects: int,
        rng: np.random.Generator,
    ) -> pd.DataFrame:
        """
        Draw ``n_objects`` (color, mag) pairs from the 2-D CMD histogram.

        Returns a DataFrame with columns ``color`` and ``mag``.
        """
        total = cmd_hist.sum()
        if total == 0 or n_objects == 0:
            return pd.DataFrame({"color": np.array([]), "mag": np.array([])})

        p = cmd_hist.ravel() / total
        flat_idx = rng.choice(len(p), size=n_objects, p=p)
        i_color, i_mag = np.unravel_index(flat_idx, cmd_hist.shape)

        # Uniform sample within each selected bin cell
        dc = color_edges[i_color + 1] - color_edges[i_color]
        dm = mag_edges[i_mag + 1] - mag_edges[i_mag]
        color = color_edges[i_color] + rng.uniform(0.0, 1.0, n_objects) * dc
        mag = mag_edges[i_mag] + rng.uniform(0.0, 1.0, n_objects) * dm

        return pd.DataFrame({"color": color, "mag": mag})

    def _sample_positions(
        self,
        n_objects: int,
        pixel: int,
        nside: int,
        rng: np.random.Generator,
    ):
        """
        Sample ``n_objects`` (RA, Dec) positions uniformly within a HEALPix pixel.

        Uses rejection sampling in the pixel's bounding box with sin(dec) weighting
        for uniform-on-sphere distribution.

        Returns
        -------
        ra, dec : np.ndarray
            Arrays of length ``n_objects`` in degrees.
        """
        vecs = hp.boundaries(nside, pixel, step=1)  # (3, 4)
        corner_ra, corner_dec = hp.vec2ang(vecs.T, lonlat=True)

        if corner_ra.max() - corner_ra.min() > 180:
            corner_ra = np.where(corner_ra > 180, corner_ra - 360, corner_ra)

        ra_min, ra_max = float(corner_ra.min()), float(corner_ra.max())
        dec_min, dec_max = float(corner_dec.min()), float(corner_dec.max())
        sin_dec_min = np.sin(np.deg2rad(dec_min))
        sin_dec_max = np.sin(np.deg2rad(dec_max))

        accepted_ra: list = []
        accepted_dec: list = []
        batch_size = max(n_objects * 8, 64)

        while len(accepted_ra) < n_objects:
            batch_ra = rng.uniform(ra_min, ra_max, batch_size) % 360.0
            sin_dec = rng.uniform(sin_dec_min, sin_dec_max, batch_size)
            batch_dec = np.rad2deg(np.arcsin(np.clip(sin_dec, -1.0, 1.0)))
            pix_test = hp.ang2pix(nside, batch_ra, batch_dec, lonlat=True)
            mask = pix_test == pixel
            accepted_ra.extend(batch_ra[mask].tolist())
            accepted_dec.extend(batch_dec[mask].tolist())

        return np.array(accepted_ra[:n_objects]), np.array(accepted_dec[:n_objects])
