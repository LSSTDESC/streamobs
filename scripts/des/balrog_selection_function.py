#!/usr/bin/env python
"""Build streamobs selection-function products from a Balrog SSI catalog.

Runs on the machine that holds the Balrog catalog and emits only small files
(~10 KB of CSVs plus, optionally, a HEALPix depth map).  Nothing here imports
streamobs; the only dependencies are numpy, h5py, pandas and (for depth maps)
healpy / healsparse.

Products, following ``streamobs/docs/source/selection_function_methodology.md``:

  <out>/<tag>_stellar_efficiency_cut<band>.csv
        mag_<band>, delta_mag, detection_eff, classification_eff,
        classification_detection_eff
  <out>/<tag>_photoerror_<band>.csv            SAMPLE  (truth scatter -> noise draw)
  <out>/<tag>_photoerror_<band>_catalog.csv    CATALOG (reported magerr -> S/N cut)
  <out>/<tag>_photoerror_<band>_nocut.csv          SAMPLE, no S/N cut  (forced bands)
  <out>/<tag>_photoerror_<band>_catalog_nocut.csv  CATALOG, no S/N cut (forced bands)
  <out>/<tag>_photoerror_<band>{,_catalog}_raw.csv    provenance, pre-afterburner
  <out>/<tag>_galaxy_misclass_cut<band>.csv        stellar contamination
  <out>/<tag>_galaxy_misclass_cut<band>_raw.csv    provenance, pre-afterburner
  <out>/<tag>_maglim_<band>_nside<N>.fits.gz   only with --write-maglim
  <out>/<tag>_audit.json                       counts + anchors for validation

Everything is keyed to ``delta_mag = mag_true - maglim(pixel)`` with maglim on a
5 sigma scale.  The S/N > 5 detection cut is *baked into* the efficiency curves;
the injector must not re-apply it.

Depth handling
--------------
The input depth map supplies the *spatial structure*; the *absolute scale* is
truth-anchored from the injections themselves — each band's map is shifted so
that its median equals the magnitude at which the truth-based scatter of
(observed - true) reaches sigma = 2.5/ln10/5 = 0.2171.  That makes
``delta_mag = 0 <-> sigma = 0.217`` by construction, so a map delivered on an
S/N = 10 scale (or any other) is corrected automatically rather than by
assuming a +2.5*log10(2) = 0.753 mag offset.  ``--no-anchor`` disables this and
trusts the map as delivered.

Two passes are made over the catalog.  Pass 1 accumulates the truth-anchor
statistics (needs no depth zero point); pass 2 accumulates every delta_mag-keyed
histogram using the anchored depth.  Both stream in chunks.

Which population defines the depth (``--anchor-sample``)
--------------------------------------------------------
The anchor sample **drops the reference-band S/N cut** (``nosnr``, the default).
The alternative -- anchoring on the same detected+classified sample the
photo-error curve describes, which would put sigma = 0.2171 at delta_mag = 0 by
construction -- was measured on the full DES Y6 catalog and rejected:

    band   nosnr shift   detected shift
    g         -0.297         +0.395
    r         -0.280         -0.637
    i         -0.180         -0.207
    z         -0.124         -0.116

Under ``detected`` the g and r depths move in opposite directions by more than a
magnitude in total, which cannot be real depth.  The reason is structural: the
S/N cut is applied in the *reference band only*, so the g anchor is truncated by
its own cut while the r/i/z anchor samples are selected on *g* S/N -- a different
and band-dependent distortion of each scatter curve.  ``nosnr`` gives smooth,
same-sign shifts that order correctly with wavelength.

This does mean the sample photo-error curve reaches sigma = 0.2171 somewhat
faintward of delta_mag = 0 rather than exactly at it.  That is already true of
every truth-anchored streamobs release for a related reason -- the maglim is
anchored to the SAMPLE curve while ``get_photo_error`` returns the CATALOG curve
-- which is why they all set ``skip_snr_maglim_check`` in the test registry.
Roman's own anchor sample does include the detection cut, but its methodology
footnote records the difference there as <= 0.01 mag, i.e. Roman's anchor is
numerically indistinguishable from the no-cut one.  DES diverges only because its
S/N cut truncates much harder relative to its error inflation (1.46 vs ~2), so
``nosnr`` is the choice that actually matches Roman's *effective* convention.

Usage
-----
    python balrog_selection_function.py \
        --survey delve \
        --catalog /path/to/BalrogOfTheStars_Catalog_V4.hdf5 \
        --maglim-map g=/path/delve_maglim_g.hsp r=/path/delve_maglim_r.hsp \
        --out ./delve_dr3_gold_products

    python balrog_selection_function.py \
        --survey des_y6 \
        --catalog /path/to/fiducial_injected_sof.hdf5 \
        --measured /path/to/fiducial_matched_measured_sof.hdf5 \
        --maglim-map g=/path/des_y6_5_sig_maglim_band_g_nside_512.hsp \
        --out ./des_y6_products
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Conventions shared with streamobs
# ---------------------------------------------------------------------------
SNR_DEPTH = 5.0
SIG_SN5 = 2.5 / np.log(10) / SNR_DEPTH  # 0.21715 mag
DET_EFF_DELTA_MAX = 1.0  # zero detection_eff faintward of this delta_mag
EFF_DELTA_MIN = -11.0  # drop efficiency rows brighter than this
MIN_COUNT_EFF = 20  # drop efficiency bins with fewer true stars
# Balrog injects each parent source many times (DES: 2.8M deep-field objects ->
# 145.7M injections, ~51 each), so rows in a magnitude bin are NOT independent --
# the effective sample size is the number of DISTINCT parent sources.  Bright
# bins are the danger: a 5.9 deg^2 deep field holds few bright stars, and those
# are saturated in the (much deeper) deep-field imaging, so their injected
# morphology is corrupt.  Left unguarded this produced a classification_eff of
# 0.685 at mag_g = 18.125 between neighbours of 1.000 and 0.947, and a run of
# implausibly exact 1.0000 bins brightward of it.  Schemas that expose a parent
# id set `src_id`; those that do not fall back to the row count alone.
# 100 distinct sources gives a binomial precision of ~2% at p ~ 0.95, and on the
# DES deep fields it lands the threshold at g ~ 18.5, just faintward of where the
# parent star sample stops being usable.  Distinct deep-field stars per 0.25 mag:
#   g 16.0-17.75:  1,1,3,0,4,0,2      (nothing at all)
#   g 17.75-18.00: 13
#   g 18.00-18.25: 49    <- gave classification_eff 0.685, 8.5 sigma below its
#                           neighbours: not noise, but a real property of 49
#                           saturated-in-the-deep-field objects, and flatly
#                           contradicted by the 98.0% the Y6 Gold paper reports
#                           over 17.5 <= i <= 22.5 (Table A.3)
#   g 18.25-18.50: 76
#   g 18.50-18.75: 123   <- first bin worth trusting
MIN_UNIQUE_SRC = 100
MIN_COUNT_PE = 20  # drop photo-error bins with fewer stars
MAG_BINS = np.arange(15.0, 29.0 + 1e-6, 0.25)
MAG_MID = 0.5 * (MAG_BINS[1:] + MAG_BINS[:-1])
DELTA_BINS = np.arange(-12.0, 2.0 + 1e-6, 0.12)
DELTA_MID = 0.5 * (DELTA_BINS[1:] + DELTA_BINS[:-1])
ANCHOR_BINS = np.arange(16.0, 28.0 + 1e-6, 0.25)
ANCHOR_MID = 0.5 * (ANCHOR_BINS[1:] + ANCHOR_BINS[:-1])
# quantile grid used to accumulate (obs - true) distributions per bin without
# holding the rows in memory
RESID_BINS = np.arange(-3.0, 3.0 + 1e-9, 0.005)
RESID_MID = 0.5 * (RESID_BINS[1:] + RESID_BINS[:-1])

BAD = -9.0e8  # anything below this is a sentinel, not a magnitude
SENTINEL = -9.999e9  # desqr's sentinel for failed BDF measurements


# ---------------------------------------------------------------------------
# Star/galaxy classification
#
# Vendored from desqr/desqr/classify.py :: bdf_interp_array_y6gold /
# bdf_extended_class_dr3gold (A. Drlica-Wagner), so the remote machine needs no
# desqr checkout.  This is the *official* DELVE DR3 Gold classifier and it
# reuses the DES Y6 Gold interpolation nodes, so DES Y6 and DELVE DR3 are
# classified identically — which is what makes the two releases comparable in
# the overlap region.
#
# It depends only on BDF_T and BDF_S2N, both of which Balrog carries.  (The
# newer EXT_XGB cannot be used here: three of its six features — CONC,
# WAVG_SPREAD_MODEL_I, WAVG_SPREADERR_MODEL_I — are not measured for injected
# sources.  See arXiv:2501.05739 Appendix A.2.)
#
# Classes:  0 high-confidence star, 1 likely star, 2 likely galaxy,
#           3 high-confidence galaxy, 4 high-purity galaxy, -9 failure.
# ---------------------------------------------------------------------------
_BDF_X = np.array(
    [-3.0, 0.79891862, 0.90845217, 0.98558583, 1.05791208, 1.13603715,
     1.22479487, 1.33572223, 1.48983602, 1.74124395, 2.43187589, 6.0]
)
_BDF_Y = np.array(
    [
        [0.028, 0.028, 0.008, 0.0, 0.004, 0.012, 0.012, 0.004, 0.012, 0.024, 0.04, 0.04],
        [-0.028, -0.028, -0.028, -0.028, -0.028, -0.028, -0.028, -0.012, 0.005, 0.022, 0.04, 0.04],
        [-0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.012, 0.008, 0.016, 0.016],
        [0.252, 0.252, 0.188, 0.14, 0.096, 0.104, 0.052, 0.048, 0.04, 0.052, 0.088, 0.088],
    ]
)


def bdf_extended_class_dr3gold(bdf_T, bdf_s2n):
    """DELVE DR3 / DES Y6 Gold BDF extended classifier.

    ``scipy.interpolate.interp1d`` in the original is linear, so ``np.interp``
    reproduces it exactly while keeping this numpy-only.
    """
    bdf_T = np.asarray(bdf_T, dtype=float)
    bdf_s2n = np.asarray(bdf_s2n, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        x = np.log10(bdf_s2n)

    sentinel = (
        np.isclose(bdf_T, SENTINEL)
        | np.isclose(bdf_s2n, SENTINEL)
        | (bdf_s2n <= 0.0)
        | (x < _BDF_X[0])
        | (x > _BDF_X[-1])
        | ~np.isfinite(x)
    )
    xc = np.clip(np.where(np.isfinite(x), x, _BDF_X[0]), _BDF_X[0], _BDF_X[-1])

    ext = np.zeros(xc.size, dtype=int)
    for y in _BDF_Y:
        ext += (bdf_T > np.interp(xc, _BDF_X, y)).astype(int)
    ext[sentinel] = -9
    return ext


# ---------------------------------------------------------------------------
# Survey schema adapters
# ---------------------------------------------------------------------------
class BalrogSchema:
    """Maps one Balrog release's column layout onto the quantities we need.

    Subclasses fill in ``read_chunk`` and return a dict with these keys, each a
    1-D array of length n (rows in the chunk):

        ra, dec            truth position, degrees
        is_star            bool, true star (truth label)
        detected           bool
        classified         bool, detected AND classified as a point source
        in_footprint       bool, sky is not masked here (applied to BOTH the
                           efficiency numerator and denominator)
        tile               bytes/str, tile name (for the zero-point audit)
        true_mag[band]     float, injected apparent magnitude **including
                           extinction**, i.e. what the survey should have
                           observed.  streamobs keys delta_mag to the reddened
                           true magnitude (``surveys.py`` docstring: "True
                           apparent magnitude(s) including extinction"), so any
                           catalog storing dereddened truth must add A_band back
                           here.
        obs_mag[band]      float, measured magnitude (nan where invalid)
        obs_magerr[band]   float, reported magnitude error (nan where invalid)
    """

    name = "base"
    bands: tuple = ()

    def open(self, args):
        raise NotImplementedError

    def n_rows(self) -> int:
        raise NotImplementedError

    def read_chunk(self, lo: int, hi: int) -> dict:
        raise NotImplementedError

    def close(self):
        pass


class DelveBalrog(BalrogSchema):
    """DELVE ``BalrogOfTheStars_Catalog_V*.hdf5``.

    Flat top-level datasets, one row per injection, ``meas_detected`` marking
    recovery.  The (N, 4) measured arrays are ordered g, r, i, z.
    """

    name = "delve"
    bands = ("g", "r", "i", "z")
    BAND_INDEX = {"g": 0, "r": 1, "i": 2, "z": 3}  # verified against truth fluxes
    ZP = 30.0
    # truth_FLUX_* is DEREDDENED (verified: adding A_sfd98 moves the median
    # (obs - true) residual from +0.129 to +0.019 in g, +0.099 to +0.023 in r),
    # so the reddened truth magnitude is  ZP - 2.5*log10(flux) + A_band.
    EXT_COL = {"g": "Ag_sfd98", "r": "Ar_sfd98", "i": "Ai_sfd98", "z": "Az_sfd98"}

    def __init__(self, ext_max=1, mag_kind=None, badpix_max=None,
                 ref_band="g", snr_detect=SNR_DEPTH, truth_mag_kind=None):
        self.ext_max = ext_max  # star iff 0 <= EXT <= ext_max
        self.mag_kind = mag_kind or "psf"  # "psf" or "bdf"
        # DELVE's truth_FLUX_* is dereddened and A_band is always added back
        # below, so this knob is DES-only; accepted and ignored here.
        self.truth_mag_kind = truth_mag_kind
        self.badpix_max = badpix_max
        self.ref_band = ref_band
        self.snr_detect = snr_detect

    def open(self, args):
        self.f = h5py.File(args.catalog, "r")
        return self

    def n_rows(self) -> int:
        return int(self.f["meas_detected"].shape[0])

    def read_chunk(self, lo, hi):
        f = self.f
        sl = slice(lo, hi)
        detected = f["meas_detected"][sl] > 0
        out = {
            "ra": f["truth_ra"][sl],
            "dec": f["truth_dec"][sl],
            "is_star": f["truth_STAR"][sl].astype(bool),
            "tile": f["tilename"][sl],
        }

        # Two DIFFERENT kinds of cut, which must not be conflated:
        #
        #  * ``in_footprint`` is a property of the SKY, not of the object -
        #    FLAGS_FOREGROUND masks bright stars, the Magellanic Clouds, etc.
        #    It is applied to the numerator AND the denominator, because a star
        #    that was never going to be measured in a masked region says nothing
        #    about detection efficiency.  streamobs handles masked sky through
        #    the footprint/maglim map; folding it into the efficiency curve as
        #    well would double-count it (and it does: leaving it in the
        #    numerator only drags the bright-end plateau down to ~0.46).
        #
        #  * ``quality`` IS part of the selection function - the per-object
        #    flags an observer applies to get a clean sample.  It stays in the
        #    numerator only, and is the origin of a bright-end plateau slightly
        #    below unity (cf. the Roman products, plateau ~0.91).
        #
        # NOTE: undetected rows carry meas_* == 0.0 (not NaN, not the sentinel),
        # so they pass every ``== 0`` test - the detection gate is what keeps
        # them out.  ``meas_bdf_deblend_flags`` is a broken constant in V4
        # (never 0) and ``meas_mask_flags`` is identically 0; neither is used.
        out["in_footprint"] = f["FLAGS_FOREGROUND"][sl] == 0
        quality = (f["meas_flags"][sl] == 0) & (f["meas_bdf_flags"][sl] == 0)
        if self.badpix_max is not None:
            quality &= f["meas_badpix_frac"][sl] < self.badpix_max
        # Reference-band S/N cut.  streamobs requires the S/N > 5 cut to be
        # BAKED INTO the efficiency curves (the injector must not re-apply it),
        # so it belongs here.  Measured effect on DELVE V4: no change at all for
        # delta_mag < -0.5 (ratio 1.000), then 0.79 -> 0.53 at delta_mag +0.25
        # and 0.71 -> 0.19 at +0.75 -- i.e. it only bites in the last magnitude,
        # where it is what makes detection cross 50% at the 5 sigma depth.
        quality_nosnr = quality.copy()
        if self.snr_detect:
            j = self.BAND_INDEX[self.ref_band]
            quality &= f["meas_psf_flux_s2n"][sl][:, j] > self.snr_detect
        out["detected"] = detected & quality
        out["detected_nosnr"] = detected & quality_nosnr

        # star classification: official DR3 Gold / Y6 Gold BDF classifier
        ext = bdf_extended_class_dr3gold(f["meas_bdf_T"][sl], f["meas_bdf_s2n"][sl])
        is_ptsrc = (ext >= 0) & (ext <= self.ext_max)
        out["ext_class"] = ext
        out["classified"] = out["detected"] & is_ptsrc
        # Anchor population: same selection WITHOUT the S/N cut.  The truth
        # anchor asks "at what magnitude does the full population's truth-based
        # scatter reach sigma=0.217"; conditioning that sample on having passed
        # S/N>5 keeps only the well-measured objects, shrinking the scatter and
        # pushing the anchor artificially deep (measured: +0.23 mag on DELVE
        # V4).  cf. the anchor-sample footnote in
        # selection_function_methodology.md, where the same effect is <0.01 mag
        # for Roman and was therefore deferred.
        out["classified_nosnr"] = out["detected_nosnr"] & is_ptsrc

        true_mag, obs_mag, obs_err = {}, {}, {}
        magcol = "meas_psf_mag" if self.mag_kind == "psf" else "meas_bdf_mag"
        fluxcol = "meas_psf_flux" if self.mag_kind == "psf" else "meas_bdf_flux"
        errcol = fluxcol + "_err"
        mm = f[magcol][sl]
        fl = f[fluxcol][sl]
        fe = f[errcol][sl]
        with np.errstate(all="ignore"):
            for b in self.bands:
                j = self.BAND_INDEX[b]
                tf = f[f"truth_FLUX_{b.upper()}"][sl]
                tm = np.where(tf > 0, self.ZP - 2.5 * np.log10(np.abs(tf)), np.nan)
                true_mag[b] = tm + f[self.EXT_COL[b]][sl]  # deredden -> reddened
                om = mm[:, j].astype(float)
                om[(om < BAD) | ~detected] = np.nan
                obs_mag[b] = om
                oe = 1.0857362 * fe[:, j] / fl[:, j]
                oe[~np.isfinite(oe) | (oe <= 0) | ~detected] = np.nan
                obs_err[b] = oe
        out["true_mag"], out["obs_mag"], out["obs_magerr"] = true_mag, obs_mag, obs_err
        return out

    def close(self):
        self.f.close()


class DesY6Balrog(BalrogSchema):
    """DES Y6 Balrog (Anbajagane et al. 2025, arXiv:2501.05683).

    Two files, not one:

      A = ``fiducial_injected_sof.hdf5``          145,724,947 rows, truth only
      B = ``fiducial_matched_measured_sof.hdf5``  119,667,013 rows, ``meas_*``

    **B is exactly A masked by ``detected_0.5arcsec_match_radius``, in the same
    row order** -- verified element-wise on ``id`` at two independent offsets.
    So no hash join is needed, just a lockstep walk.  But ``chunks()`` is
    iterated TWICE by ``main`` (zero-point audit, then products), so a stateful
    cursor into B would desynchronise on the second pass.  ``open()`` therefore
    precomputes the cumulative detected count, making ``read_chunk`` idempotent.

    Truth star/galaxy label
    -----------------------
    Neither file carries one.  The paper's own definition is the *colour-based*
    ``KNN_CLASS`` of the parent Y3 deep-field catalogue, joined on
    ``id`` == deep-field ``COADD_OBJECT_ID`` (matches 100% of injections).
    Build the lookup with ``scripts/des/build_des_truth_labels.py``:

        KNN_CLASS 2 -> star, 1 -> galaxy, 0 -> UNCLASSIFIED (22.5% of
        injections; no usable NIR photometry, so the colour classifier could
        not run).  Class 0 is excluded from ``in_footprint``, which removes it
        from the numerator AND denominator of every curve -- it is neither a
        star nor a known galaxy, so it cannot inform either.

    A ``|bdf_T| < 0.02`` morphological proxy is NOT used: it is only 36.8% pure.

    Star classification
    -------------------
    ``EXT_XGB`` is not computable for injections (arXiv:2501.05739 App. A.2), so
    ``classified`` uses the surrogate from
    ``scripts/des/build_des_xgb_surrogate.py``, and ``main`` deconvolves the
    resulting efficiency with that script's per-magnitude confusion table.
    """

    name = "des_y6"
    bands = ("g", "r", "i", "z")
    BAND_INDEX = {"g": 0, "r": 1, "i": 2, "z": 3}  # verified via ext_mags ordering
    KNN_STAR, KNN_GALAXY, KNN_UNCLASSIFIED = 2, 1, 0

    def __init__(self, ext_max=1, mag_kind=None, badpix_max=None,
                 ref_band="g", snr_detect=SNR_DEPTH, truth_mag_kind="deredden"):
        self.ext_max = ext_max
        # bdf, not psf: the truth side is the deep-field BDF magnitude, so
        # comparing against meas_bdf_mag keeps (obs - true) apples-to-apples.
        self.mag_kind = mag_kind or "bdf"
        self.badpix_max = badpix_max
        self.ref_band = ref_band
        self.snr_detect = snr_detect
        self.truth_mag_kind = truth_mag_kind
        self._surrogate = None
        self._features = None
        self._truth = None
        self._fg_mask = None
        self._fg_nside = 0

    def open(self, args):
        if not getattr(args, "measured", None):
            raise SystemExit("--measured is required for --survey des_y6")
        if not getattr(args, "truth_labels", None):
            raise SystemExit("--truth-labels is required for --survey des_y6")
        if not getattr(args, "surrogate", None):
            raise SystemExit("--surrogate is required for --survey des_y6")

        self.fa = h5py.File(args.catalog, "r")
        self.fb = h5py.File(args.measured, "r")

        det = self.fa["detected_0.5arcsec_match_radius"][:] > 0
        n_b = int(self.fb["meas_id"].shape[0])
        if int(det.sum()) != n_b:
            raise SystemExit(
                f"lockstep broken: A has {int(det.sum()):,} detected rows but B "
                f"has {n_b:,}. The two files are not a matched pair."
            )
        # cum[i] = number of detected rows strictly before i  ->  the B offset
        # for a chunk starting at i.  Idempotent across repeated passes.
        self._det = det
        self._cum = np.concatenate(([0], np.cumsum(det)))
        print(f"  des_y6: A={det.size:,} rows, B={n_b:,} detected (lockstep verified)")

        self._fg_mask = self._build_foreground_mask(args)

        truth = pd.read_parquet(args.truth_labels, columns=["ID", "KNN_CLASS"])
        self._truth = pd.Series(truth["KNN_CLASS"].to_numpy(),
                                index=truth["ID"].to_numpy())
        print(f"  des_y6: truth labels for {len(self._truth):,} deep-field objects")

        import xgboost as xgb

        self._surrogate = xgb.XGBClassifier()
        self._surrogate.load_model(args.surrogate)
        self._features = Path(args.surrogate).with_name(
            "des_y6_xgb_features.txt").read_text().split()
        self._thresh = args.surrogate_threshold
        print(f"  des_y6: surrogate with {len(self._features)} features, "
              f"threshold {self._thresh}")
        return self

    def _build_foreground_mask(self, args):
        """Turn ``meas_FLAGS_FOREGROUND`` into a *positional* sky mask.

        The flag lives only in the measured file, so it is unknown for
        undetected injections.  Cutting on it directly would drop masked objects
        from the efficiency numerator while leaving them in the denominator --
        and it is not a rare flag (7.5% of detected rows are non-zero), so that
        asymmetry would bias ``detection_eff`` low by several percent.

        But foreground masking is a property of the SKY (bright stars, the
        Magellanic Clouds), not of the object, so it can be recovered by
        position: a HEALPix pixel is masked if most detections in it are
        flagged.  That mask then applies to detected and undetected rows alike,
        which is the convention the curves require (cf. the DELVE case, where
        the flag ships for every row and leaving it in the numerator alone drags
        the bright-end plateau from ~0.95 to ~0.46).

        Returns ``None`` when disabled, which falls back to the depth map's own
        footprint as the only sky cut.
        """
        if args.foreground_nside <= 0:
            print("  des_y6: positional foreground mask DISABLED")
            return None

        import healpy as hp

        nside = args.foreground_nside
        npix = hp.nside2npix(nside)
        n_flag = np.zeros(npix, dtype=np.int32)
        n_tot = np.zeros(npix, dtype=np.int32)
        n_b = int(self.fb["meas_id"].shape[0])
        step = 10_000_000
        for lo in range(0, n_b, step):
            hi = min(lo + step, n_b)
            ra = self.fb["meas_ra"][lo:hi]
            dec = self.fb["meas_dec"][lo:hi]
            fg = self.fb["meas_FLAGS_FOREGROUND"][lo:hi]
            ok = np.isfinite(ra) & np.isfinite(dec)
            pix = hp.ang2pix(nside, ra[ok], dec[ok], lonlat=True)
            np.add.at(n_tot, pix, 1)
            np.add.at(n_flag, pix, (fg[ok] != 0).astype(np.int32))
            print(f"    foreground mask {hi:,}/{n_b:,}", end="\r")
        print()
        seen = n_tot > 0
        masked = np.zeros(npix, dtype=bool)
        masked[seen] = (n_flag[seen] / n_tot[seen]) > 0.5
        frac = masked[seen].sum() / max(seen.sum(), 1)
        print(f"  des_y6: foreground mask nside={nside}, "
              f"{seen.sum():,} populated pixels, {frac:.4f} masked")
        self._fg_nside = nside
        return masked

    def n_rows(self) -> int:
        return int(self._det.size)

    def _features_frame(self, m: dict) -> np.ndarray:
        """Assemble the surrogate feature matrix in the trained column order.

        Must mirror ``build_des_xgb_surrogate.derive`` exactly; the names are the
        real-catalogue ones and the mapping to ``meas_*`` is the contract
        documented in that script.
        """
        with np.errstate(all="ignore"):
            col = {
                "BDF_T": m["bdf_T"],
                "BDF_T_ERR": m["bdf_T_err"],
                "BDF_T_RATIO": m["bdf_T_ratio"],
                "log_bdf_s2n": np.log10(np.where(m["bdf_s2n"] > 0, m["bdf_s2n"], np.nan)),
                "PSF_T": m["psf_T"],
                "psf_bdf_T": m["psf_T"] - m["bdf_T"],
                "conc_gap": m["gap_mag"][:, 2] - m["bdf_mag"][:, 2],
                "conc_gap_r": m["gap_mag"][:, 1] - m["bdf_mag"][:, 1],
                "conc_ap8": m["ap8_mag"][:, 2] - m["bdf_mag"][:, 2],
                "conc_ap8_r": m["ap8_mag"][:, 1] - m["bdf_mag"][:, 1],
                "BDF_MAG_I": m["bdf_mag"][:, 2],
                "BDF_MAG_G": m["bdf_mag"][:, 0],
            }
        return np.column_stack([col[f] for f in self._features])

    def read_chunk(self, lo, hi):
        a, b = self.fa, self.fb
        sl = slice(lo, hi)
        det = self._det[sl]
        b0, b1 = int(self._cum[lo]), int(self._cum[hi])
        bsl = slice(b0, b1)
        n = hi - lo

        def scatter(vals, dtype=float, fill=np.nan):
            """Place a detected-only (B) column back on the full (A) grid."""
            out = np.full(n, fill, dtype=dtype)
            out[det] = vals
            return out

        ids = a["id"][sl]
        knn = self._truth.reindex(ids).to_numpy()

        out = {
            "ra": a["ra"][sl],
            "dec": a["dec"][sl],
            "is_star": knn == self.KNN_STAR,
            "tile": a["wide_tilename"][sl],  # NOT `tilename` (deep-field source tile)
            # Parent deep-field object. Each is injected ~51 times, so this is
            # what makes the effective sample size countable (see MIN_UNIQUE_SRC).
            "src_id": ids,
        }

        # Sky mask AND label availability gate the denominator.  KNN_CLASS == 0
        # is unclassified, not galaxy: keeping it would contaminate the galaxy
        # misclassification denominator and dilute the star denominator.
        in_fp = np.isin(knn, (self.KNN_STAR, self.KNN_GALAXY))
        if self._fg_mask is not None:
            import healpy as hp

            pix = hp.ang2pix(self._fg_nside, out["ra"], out["dec"], lonlat=True)
            in_fp &= ~self._fg_mask[pix]  # positional: applies to ALL rows
        out["in_footprint"] = in_fp

        quality = det.copy()
        quality &= scatter(b["meas_flags"][bsl], dtype=float, fill=1.0) == 0
        quality &= scatter(b["meas_bdf_flags"][bsl], dtype=float, fill=1.0) == 0
        # <0.2% of injections land within 1.5" of a real object; the paper
        # removes them from every analysis (Sec. 3.4).
        quality &= scatter(b["meas_BALROG_FLAG_BLEND"][bsl], dtype=float, fill=1.0) == 0
        if self.badpix_max is not None:
            quality &= scatter(b["meas_badpix_frac"][bsl], dtype=float,
                               fill=np.inf) < self.badpix_max

        quality_nosnr = quality.copy()
        if self.snr_detect:
            j = self.BAND_INDEX[self.ref_band]
            # Gate on the S/N of the SAME photometry the magnitudes and errors
            # come from.  meas_psf_flux_s2n would be inconsistent under
            # --mag-kind bdf (the DES default), and cutting on one photometry
            # family while measuring the scatter of another distorts the
            # surviving population's error distribution.
            if self.mag_kind == "bdf":
                fl = np.asarray(b["meas_bdf_flux"][bsl][:, j], dtype=float)
                fe = np.asarray(b["meas_bdf_flux_err"][bsl][:, j], dtype=float)
                with np.errstate(all="ignore"):
                    s2n_det = np.where(fe > 0, fl / fe, 0.0)
                s2n = scatter(s2n_det, fill=0.0)
            else:
                s2n = scatter(b["meas_psf_flux_s2n"][bsl][:, j], fill=0.0)
            quality &= s2n > self.snr_detect
        out["detected"] = quality
        out["detected_nosnr"] = quality_nosnr

        # --- surrogate classification -------------------------------------
        meas = {
            "bdf_T": b["meas_bdf_T"][bsl],
            "bdf_T_err": b["meas_bdf_T_err"][bsl],
            "bdf_T_ratio": b["meas_bdf_T_ratio"][bsl],
            "bdf_s2n": b["meas_bdf_s2n"][bsl],
            "psf_T": b["meas_psf_T"][bsl],
            "bdf_mag": b["meas_bdf_mag"][bsl],
            "gap_mag": b["meas_gap_mag"][bsl],
            "ap8_mag": b["meas_psf_mag_aper8"][bsl],
        }
        for k, v in meas.items():
            v = np.asarray(v, dtype=float)
            v[v < BAD] = np.nan
            if k.endswith("mag"):
                v[v > 37.0] = np.nan
            meas[k] = v

        X = self._features_frame(meas)
        is_ptsrc_det = np.zeros(b1 - b0, dtype=bool)
        usable = np.isfinite(X).all(axis=1)
        if usable.any():
            prob = self._surrogate.predict_proba(X[usable])[:, 1]
            is_ptsrc_det[usable] = prob >= self._thresh
        is_ptsrc = scatter(is_ptsrc_det, dtype=bool, fill=False)

        out["classified"] = out["detected"] & is_ptsrc
        out["classified_nosnr"] = out["detected_nosnr"] & is_ptsrc

        # --- magnitudes ----------------------------------------------------
        # Truth side.  streamobs keys delta_mag to the *reddened* true apparent
        # magnitude -- what the survey should have observed (surveys.py:617).
        # Like DELVE, DES needs extinction handled explicitly, and the deciding
        # diagnostic is the SPREAD of the per-tile (obs - true) offsets, not
        # their median.  Measured on 4M injections:
        #
        #   truth_mag_kind          ref offset   per-tile spread(16-84)/2
        #   raw       bdf_mag           -0.0125            0.0300
        #   reddened  bdf_mag+ext       -0.0625            0.0025
        #   deredden  bdf_mag_dered+ext -0.0275            0.0025
        #
        # Extinction varies tile to tile, so leaving it out stamps every tile
        # with its own offset -- hence `raw`'s 12x wider spread.  `deredden` is
        # both the physically correct construction (de-redden at the deep field,
        # re-redden at the injection position) and the smallest residual, so it
        # is the default.  A global median alone is misleading here: `raw` has
        # the smallest |median| and is nonetheless wrong.
        true_mag, obs_mag, obs_err = {}, {}, {}
        tcol = "bdf_mag_deredden" if self.truth_mag_kind == "deredden" else "bdf_mag"
        tm_all = np.asarray(a[tcol][sl], dtype=float)
        tm_all[tm_all < BAD] = np.nan
        tm_all[tm_all > 37.0] = np.nan
        if self.truth_mag_kind == "raw":
            ext = np.zeros_like(tm_all)
        else:
            ext = scatter_2d(b["ext_mags"][bsl], det, n)

        magcol = "meas_bdf_mag" if self.mag_kind == "bdf" else "meas_psf_mag"
        fluxcol = "meas_bdf_flux" if self.mag_kind == "bdf" else "meas_psf_flux"
        mm = np.asarray(b[magcol][bsl], dtype=float)
        fl = np.asarray(b[fluxcol][bsl], dtype=float)
        fe = np.asarray(b[fluxcol + "_err"][bsl], dtype=float)
        with np.errstate(all="ignore"):
            for band in self.bands:
                j = self.BAND_INDEX[band]
                true_mag[band] = tm_all[:, j] + ext[:, j]
                om = scatter(mm[:, j])
                om[(om < BAD) | (om > 37.0)] = np.nan
                obs_mag[band] = om
                # No per-band magnitude error column exists in either file.
                oe = scatter(1.0857362 * fe[:, j] / fl[:, j])
                oe[~np.isfinite(oe) | (oe <= 0)] = np.nan
                obs_err[band] = oe
        out["true_mag"], out["obs_mag"], out["obs_magerr"] = true_mag, obs_mag, obs_err
        return out

    def close(self):
        self.fa.close()
        self.fb.close()


def scatter_2d(vals, det, n):
    """Place a detected-only (N,4) block back on the full row grid."""
    out = np.full((n, vals.shape[1]), np.nan)
    v = np.asarray(vals, dtype=float)
    v[v < BAD] = np.nan
    out[det] = v
    return out


SCHEMAS = {"delve": DelveBalrog, "des_y6": DesY6Balrog}


# ---------------------------------------------------------------------------
# Depth maps
# ---------------------------------------------------------------------------
VALID_MAGLIM = (15.0, 30.0)  # anything outside this is junk, not a depth
MEDIAN_NSIDE = 2048  # resolution used only for the scalar reference median


class MaglimMosaic:
    """One or more depth maps covering complementary sky, sampled as a unit.

    The DELVE DR3 depth is split across two disjoint files —
    ``delve_dr32_*_maglim_wmean.hsp`` and
    ``delve_dr311+dr312_*_maglim_Nov28th.hsp`` — which between them cover
    99.76% of the Balrog injections while overlapping on exactly 0% (measured).
    Either one alone throws away roughly half the footprint, so they are
    mosaicked here.  Where maps do overlap, the first listed wins.
    """

    def __init__(self, paths):
        self.parts = [MaglimMap(p) for p in paths]
        self.nside = max(p.nside for p in self.parts)
        self._offset = 0.0

    @property
    def offset(self):
        return self._offset

    @offset.setter
    def offset(self, v):
        self._offset = v
        for p in self.parts:
            p.offset = v

    @property
    def median(self):
        # area-weighted across the parts, then shifted
        w = np.array([p.area_deg2 for p in self.parts])
        m = np.array([p._raw_median for p in self.parts])
        good = np.isfinite(m) & (w > 0)
        return float(np.sum(w[good] * m[good]) / np.sum(w[good])) + self._offset

    @property
    def area_deg2(self):
        return float(sum(p.area_deg2 for p in self.parts))

    def sample(self, ra, dec):
        out = self.parts[0].sample(ra, dec)
        for p in self.parts[1:]:
            miss = ~np.isfinite(out)
            if not miss.any():
                break
            out[miss] = p.sample(np.asarray(ra)[miss], np.asarray(dec)[miss])
        return out

    def to_healpix(self, nside_out):
        out = self.parts[0].to_healpix(nside_out)
        for p in self.parts[1:]:
            d = p.to_healpix(nside_out)
            miss = ~np.isfinite(out)
            out[miss] = d[miss]
        return out

    def describe(self):
        return " + ".join(
            f"{Path(p.path).name}(nside={p.nside}, {p.area_deg2:.0f}deg2, "
            f"med={p._raw_median:.3f})"
            for p in self.parts
        )


class MaglimMap:
    """A depth map sampled at its **native** resolution.

    HealSparse maps are kept sparse and queried with ``get_values_pos``, so a
    6 GB DELVE DR3.2 map is never densified — we get full native resolution for
    the per-object depth without paying for a full-sky dense array.  A dense
    map is only ever built on demand, for ``--write-maglim``.

    ``offset`` carries the truth-anchor shift and is applied at sample time, so
    anchoring never rewrites the map.

    Values outside ``VALID_MAGLIM`` are treated as missing: the DES Y6 maps
    ship with junk outliers (r/i max 37.5/34.1, y-band min -91) that would
    otherwise poison both the median and the per-object depth.
    """

    def __init__(self, path):
        import healpy as hp

        self.path = str(path)
        self.offset = 0.0
        self.sparse = None
        self.dense = None
        if self.path.endswith((".hsp", ".hs")):
            import healsparse as hsp

            self.sparse = hsp.HealSparseMap.read(self.path)
            self.nside = self.sparse.nside_sparse
            # The reference median is a single scalar, so take it from a
            # degraded copy: materialising valid_pixels at nside 16384 would
            # mean billions of int64 indices for no gain in the median.
            ref = self.sparse
            if ref.nside_sparse > MEDIAN_NSIDE:
                ref = ref.degrade(MEDIAN_NSIDE, reduction="mean")
            vals = ref[ref.valid_pixels].astype(float)
            del ref
        else:
            d = np.asarray(hp.read_map(self.path), dtype=float)
            d[(d == hp.UNSEEN) | ~np.isfinite(d)] = np.nan
            self.dense = d
            self.nside = hp.npix2nside(d.size)
            vals = d[np.isfinite(d)]
        good = np.isfinite(vals) & (vals >= VALID_MAGLIM[0]) & (vals <= VALID_MAGLIM[1])
        self.n_valid = int(good.sum())
        self.n_clipped = int(vals.size - good.sum())
        self._raw_median = float(np.median(vals[good])) if good.any() else np.nan
        # n_valid is counted at MEDIAN_NSIDE for HealSparse inputs, so derive
        # the area from that resolution rather than from nside_sparse
        nside_stat = min(self.nside, MEDIAN_NSIDE) if self.sparse is not None else self.nside
        self.area_deg2 = self.n_valid * hp.nside2pixarea(nside_stat, degrees=True)

    @property
    def median(self):
        return self._raw_median + self.offset

    def sample(self, ra, dec):
        """Depth at each position, nan off-footprint or where the value is junk."""
        if self.sparse is not None:
            v = np.asarray(
                self.sparse.get_values_pos(np.asarray(ra), np.asarray(dec), lonlat=True),
                dtype=float,
            )
        else:
            import healpy as hp

            pix = hp.ang2pix(self.nside, np.asarray(ra), np.asarray(dec), lonlat=True)
            v = self.dense[pix]
        v = np.where(
            np.isfinite(v) & (v >= VALID_MAGLIM[0]) & (v <= VALID_MAGLIM[1]), v, np.nan
        )
        return v + self.offset

    def to_healpix(self, nside_out):
        """Dense RING map at ``nside_out``, with the anchor offset applied."""
        import healpy as hp

        if self.sparse is not None:
            m = self.sparse
            if m.nside_sparse > nside_out:
                m = m.degrade(nside_out, reduction="mean")
            d = np.asarray(m.generate_healpix_map(nest=False), dtype=float)
        else:
            d = self.dense.copy()
            if self.nside > nside_out:
                d = hp.ud_grade(d, nside_out, order_in="RING")
        d[(d == hp.UNSEEN) | ~np.isfinite(d)] = np.nan
        d = np.where(
            (d >= VALID_MAGLIM[0]) & (d <= VALID_MAGLIM[1]), d + self.offset, np.nan
        )
        return d


# ---------------------------------------------------------------------------
# Streaming accumulators
# ---------------------------------------------------------------------------
class ResidualHist:
    """2-D histogram of (obs - true) residuals in bins of some key variable.

    Lets us recover per-bin percentiles (and hence the truth-based scatter)
    without ever holding the rows in memory.
    """

    def __init__(self, key_bins):
        self.key_bins = key_bins
        self.h = np.zeros((key_bins.size - 1, RESID_BINS.size - 1), dtype=np.int64)

    def add(self, key, resid):
        m = np.isfinite(key) & np.isfinite(resid)
        if not m.any():
            return
        self.h += np.histogram2d(key[m], resid[m], bins=[self.key_bins, RESID_BINS])[
            0
        ].astype(np.int64)

    def counts(self):
        return self.h.sum(axis=1)

    def percentile(self, q):
        """Per-key-bin percentile of the residual distribution (nan if empty)."""
        out = np.full(self.h.shape[0], np.nan)
        tot = self.h.sum(axis=1)
        cum = np.cumsum(self.h, axis=1)
        for i in np.where(tot > 0)[0]:
            target = q / 100.0 * tot[i]
            j = int(np.searchsorted(cum[i], target))
            j = min(j, RESID_MID.size - 1)
            out[i] = RESID_MID[j]
        return out

    def scatter(self):
        """(p84 - p16) / 2 per key bin."""
        return (self.percentile(84) - self.percentile(16)) / 2.0


class TileOffsets:
    """Per-tile median of (observed - reddened true) for bright true stars.

    A well-calibrated tile sits at a common small offset.  In DELVE Balrog V4 a
    subset of tiles (concentrated at RA >~ 310 and RA <~ 80) instead sit up to
    -0.72 mag away, in discrete per-tile steps — a truth-side zero-point
    problem, not photon noise.  Left uncorrected these tiles inflate the
    apparent photometric scatter by more than an order of magnitude and would
    wreck both the error model and the truth anchor.
    """

    OFF_BINS = np.arange(-1.5, 1.5 + 1e-9, 0.005)
    OFF_MID = 0.5 * (OFF_BINS[1:] + OFF_BINS[:-1])

    def __init__(self):
        self.index = {}
        self.rows = []

    def add(self, tiles, resid):
        m = np.isfinite(resid)
        if not m.any():
            return
        tiles, resid = np.asarray(tiles)[m], resid[m]
        for t in np.unique(tiles):
            key = t.decode() if isinstance(t, bytes) else str(t)
            if key not in self.index:
                self.index[key] = len(self.rows)
                self.rows.append(np.zeros(self.OFF_MID.size, dtype=np.int64))
            self.rows[self.index[key]] += np.histogram(
                resid[tiles == t], bins=self.OFF_BINS
            )[0].astype(np.int64)

    def medians(self, min_count=30):
        out = {}
        for key, i in self.index.items():
            h = self.rows[i]
            tot = h.sum()
            if tot < min_count:
                continue
            j = int(np.searchsorted(np.cumsum(h), 0.5 * tot))
            out[key] = float(self.OFF_MID[min(j, self.OFF_MID.size - 1)])
        return out


class MedianHist:
    """Histogram of log10(reported magerr) in bins of a key variable."""

    LOG_BINS = np.arange(-4.0, 1.0 + 1e-9, 0.005)
    LOG_MID = 0.5 * (LOG_BINS[1:] + LOG_BINS[:-1])

    def __init__(self, key_bins):
        self.key_bins = key_bins
        self.h = np.zeros((key_bins.size - 1, self.LOG_BINS.size - 1), dtype=np.int64)

    def add(self, key, value):
        m = np.isfinite(key) & np.isfinite(value)
        if not m.any():
            return
        self.h += np.histogram2d(
            key[m], value[m], bins=[self.key_bins, self.LOG_BINS]
        )[0].astype(np.int64)

    def counts(self):
        return self.h.sum(axis=1)

    def median(self):
        out = np.full(self.h.shape[0], np.nan)
        tot = self.h.sum(axis=1)
        cum = np.cumsum(self.h, axis=1)
        for i in np.where(tot > 0)[0]:
            j = int(np.searchsorted(cum[i], 0.5 * tot[i]))
            j = min(j, self.LOG_MID.size - 1)
            out[i] = self.LOG_MID[j]
        return out


# ---------------------------------------------------------------------------
# Afterburner (mirrors the Roman/LSST generators)
# ---------------------------------------------------------------------------
def apply_photoerr_corrections(tab, curve_id, corrections_path):
    if corrections_path is None or not Path(corrections_path).exists():
        return tab.copy()
    import yaml

    with open(corrections_path) as fh:
        corr = yaml.safe_load(fh) or {}
    rules = corr.get(curve_id, {})
    out = tab.copy()
    for rule, params in rules.items():
        if rule == "clamp_faint":
            m = out["delta_mag"] >= float(params["delta_mag_min"])
            out.loc[m, "log_mag_err"] = np.maximum(
                out.loc[m, "log_mag_err"], float(params["value"])
            )
            print(f"  [afterburner] {curve_id}: clamp_faint on {int(m.sum())} bins")
        elif rule == "cut_bright":
            n = int((out["delta_mag"] < float(params["delta_mag_min"])).sum())
            out = out[out["delta_mag"] >= float(params["delta_mag_min"])].copy()
            print(f"  [afterburner] {curve_id}: cut_bright dropped {n} bins")
        else:
            print(f"  [afterburner] WARNING unknown rule {rule!r} for {curve_id}")
    return out


def reapply_corrections(out_dir, tag, band, corrections_path):
    """Re-run the afterburner from the saved ``*_raw.csv`` provenance files.

    The raw curves exist precisely so the manual cleanup can be revised without
    re-deriving anything -- the reduction itself is a multi-hour pass over
    hundreds of GB, and tuning a YAML threshold should not cost that.
    """
    out = Path(out_dir)
    pe_header = "delta_mag,log_mag_err"
    mis_header = f"mag_{band},delta_mag,missclassification_eff"
    n = 0
    # each curve carries its own header: the misclassification table has three
    # columns, not the photo-error pair's two, and writing the wrong one leaves
    # a file SurveyFactory loads as None through its bare except
    for name, cid, header in [
        (f"{tag}_photoerror_{band}", f"{band}_sample", pe_header),
        (f"{tag}_photoerror_{band}_catalog", f"{band}_catalog", pe_header),
        (f"{tag}_photoerror_{band}_nocut", f"{band}_sample_nocut", pe_header),
        (f"{tag}_photoerror_{band}_catalog_nocut", f"{band}_catalog_nocut", pe_header),
        (f"{tag}_galaxy_misclass_cut{band}", f"{band}_misclass", mis_header),
    ]:
        raw = out / f"{name}_raw.csv"
        if not raw.exists():
            print(f"  missing {raw}, skipping")
            continue
        tab = pd.read_csv(raw)
        if list(tab.columns) != header.split(","):
            raise SystemExit(
                f"{raw.name}: columns {list(tab.columns)} do not match the "
                f"expected header {header!r}"
            )
        cleaned = apply_photoerr_corrections(tab, cid, corrections_path)
        write_csv(out / f"{name}.csv", cleaned, header)
        n += 1
    if not n:
        raise SystemExit(f"no *_raw.csv found in {out}")
    print(f"re-applied {corrections_path} to {n} curves in {out}")


def deconvolve_classification_eff(eff_cls, conf, mag_mid, min_sep=0.2,
                                  min_count=200):
    """Invert a surrogate's selection back onto the classifier it approximates.

    The surrogate ``S`` is not the real classifier ``X``, so what the injections
    measure is ``eff_S = P(S=1|star)`` rather than ``eff_X = P(X=1|star)``.
    Given the per-magnitude confusion measured on the real catalogue,

        a = P(S=1 | X=1)        b = P(S=1 | X=0)
        eff_S = a*eff_X + b*(1 - eff_X)   ->   eff_X = (eff_S - b)/(a - b)

    Only bins where the two classes are actually separated are corrected: as
    ``a -> b`` the inversion amplifies noise without bound, so ``min_sep``
    gates it and everything else is left as measured.

    Returns ``(eff_out, info)``.
    """
    good = (conf["n_pos"] > min_count) & (conf["n_neg"] > min_count)
    if good.sum() < 2:
        # Nothing well-measured to interpolate between: leave the curve as
        # measured rather than crashing or inventing a correction.
        return np.asarray(eff_cls), {
            "applied": False, "reason": "too few usable confusion bins",
            "n_bins_corrected": 0, "median_abs_change": 0.0,
            "max_abs_change": 0.0,
        }
    a_i = np.interp(mag_mid, conf.loc[good, "mag_g"], conf.loc[good, "a"],
                    left=np.nan, right=np.nan)
    b_i = np.interp(mag_mid, conf.loc[good, "mag_g"], conf.loc[good, "b"],
                    left=np.nan, right=np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        sep = a_i - b_i
        corrected = (eff_cls - b_i) / sep
    ok = np.isfinite(corrected) & (sep > min_sep)
    out = np.where(ok, np.clip(corrected, 0.0, 1.0), eff_cls)
    changed = np.abs(np.nan_to_num(out) - np.nan_to_num(eff_cls))
    info = {
        "applied": True,
        "n_bins_corrected": int(ok.sum()),
        "median_abs_change": float(np.nanmedian(changed[ok])) if ok.any() else 0.0,
        "max_abs_change": float(np.nanmax(changed)) if ok.any() else 0.0,
    }
    return out, info


def write_csv(path, df, header):
    np.savetxt(path, df.values, delimiter=",", header=header, fmt="%.6f", comments="")
    print(f"  wrote {path} ({len(df)} rows)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    tag = args.tag or args.survey

    if args.reapply_corrections:
        return reapply_corrections(out, tag, args.band, args.corrections)

    # Required for a real reduction, but not when only re-cleaning curves, so
    # they are validated here rather than by argparse.
    for flag, val in (("--catalog", args.catalog), ("--maglim-map", args.maglim_map)):
        if not val:
            raise SystemExit(f"{flag} is required (omit only with --reapply-corrections)")

    schema = SCHEMAS[args.survey](
        ext_max=args.ext_max,
        mag_kind=args.mag_kind,
        badpix_max=args.badpix_max,
        ref_band=args.band,
        snr_detect=args.snr_detect,
        truth_mag_kind=args.truth_mag_kind,
    ).open(args)
    n = schema.n_rows()
    row0 = max(0, args.row_start)
    if args.max_rows:
        n = min(n, row0 + args.max_rows)
    print(f"{tag}: rows {row0:,}-{n:,}, bands {schema.bands}, ref band {args.band}")

    maps = {}
    for spec in args.maglim_map:
        b, path = spec.split("=", 1)
        maps[b] = MaglimMosaic(path.split(","))
        mm = maps[b]
        print(f"  depth map {b}: {mm.describe()}")
        print(f"      -> {mm.area_deg2:.0f} deg2, median {mm.median:.3f}")
    if args.band not in maps:
        raise SystemExit(f"no depth map given for reference band {args.band!r}")

    def chunks():
        for lo in range(row0, n, args.chunk):
            yield lo, min(lo + args.chunk, n)

    ref = args.band

    # ---- pass 1: per-tile zero-point audit ------------------------------
    tiles = TileOffsets()
    print("\npass 1/2: per-tile zero-point audit")
    for lo, hi in chunks():
        c = schema.read_chunk(lo, hi)
        sel = c["is_star"] & c["classified"]
        bright = (
            sel
            & c["in_footprint"]
            & (c["true_mag"][ref] > args.zp_mag_min)
            & (c["true_mag"][ref] < args.zp_mag_max)
        )
        tiles.add(c["tile"][bright], (c["obs_mag"][ref] - c["true_mag"][ref])[bright])
        print(f"  {hi:,}/{n:,}", end="\r")
    print()

    tile_off = tiles.medians()
    bad_tiles = {}
    if tile_off:
        vals = np.array(list(tile_off.values()))
        ref_off = float(np.median(vals))
        bad_tiles = {
            k: v for k, v in tile_off.items() if abs(v - ref_off) > args.tile_zp_max
        }
        print(
            f"  tile zero-points: {len(tile_off)} tiles measured, reference offset "
            f"{ref_off:+.4f}, spread(16-84)/2 "
            f"{(np.percentile(vals, 84) - np.percentile(vals, 16)) / 2:.4f}"
        )
        print(
            f"  {len(bad_tiles)} tiles deviate by > {args.tile_zp_max} mag "
            f"-> {args.tile_zp}"
        )
    else:
        ref_off = 0.0
        print("  tile zero-points: no tile had enough bright stars; skipping")

    def tile_shift(tile_arr):
        """Per-row zero-point correction to subtract from the observed mag."""
        corr = np.zeros(tile_arr.size)
        drop = np.zeros(tile_arr.size, dtype=bool)
        if args.tile_zp == "none" or not tile_off:
            return corr, drop
        for t in np.unique(tile_arr):
            key = t.decode() if isinstance(t, bytes) else str(t)
            m = tile_arr == t
            if key in bad_tiles:
                if args.tile_zp == "reject":
                    drop |= m
                    continue
                corr[m] = bad_tiles[key] - ref_off
            elif key in tile_off:
                corr[m] = tile_off[key] - ref_off
        return corr, drop

    # ---- pass 2: truth anchor + delta_mag-keyed products ----------------
    # The anchor must be measured *after* the tile correction, but the anchor
    # only shifts each depth map by a constant, and a constant map shift simply
    # translates the delta_mag axis.  So delta is accumulated against the raw
    # map here and the axis is relabelled once the shift is known — exact, and
    # it keeps this to two passes over the catalog.
    print("\npass 2/2: truth anchor + efficiency + photo-error + misclassification")
    anchor = {b: ResidualHist(ANCHOR_BINS) for b in maps}
    # distinct parent sources per magnitude bin, packed as bin*2**32 + src_id
    uniq_src = set()       # true stars
    uniq_src_gal = set()   # true galaxies (the misclassification denominator)
    n_all = np.zeros(MAG_BINS.size - 1, dtype=np.int64)  # all injected true stars
    n_det = np.zeros_like(n_all)
    n_cls = np.zeros_like(n_all)
    n_gal_det = np.zeros_like(n_all)  # detected true galaxies
    n_gal_cls = np.zeros_like(n_all)  # ... misclassified as stars
    pe_sample = ResidualHist(DELTA_BINS)
    pe_catalog = MedianHist(DELTA_BINS)
    # Forced-photometry counterparts: the same two curves measured on the
    # population that has NOT had the reference-band S/N cut applied.  Only
    # the reference band's photometry is conditioned on its own detection;
    # every other band is forced, so it must not inherit that conditioning.
    pe_sample_nc = ResidualHist(DELTA_BINS)
    pe_catalog_nc = MedianHist(DELTA_BINS)
    n_dropped = 0

    med_raw = {b: maps[b].median for b in maps}  # offsets are still zero here

    for lo, hi in chunks():
        c = schema.read_chunk(lo, hi)
        corr, drop = tile_shift(c["tile"])
        n_dropped += int(drop.sum())
        keep = ~drop

        ml = maps[ref].sample(c["ra"], c["dec"])
        tm_ref = c["true_mag"][ref]
        ok = np.isfinite(ml) & np.isfinite(tm_ref) & keep & c["in_footprint"]
        delta = tm_ref - ml  # against the RAW map; axis relabelled below

        star = c["is_star"] & ok
        gal = (~c["is_star"]) & ok
        n_all += np.histogram(tm_ref[star], MAG_BINS)[0]
        if "src_id" in c:
            for sel_mask, acc in ((star, uniq_src), (gal, uniq_src_gal)):
                ib = np.digitize(tm_ref[sel_mask], MAG_BINS) - 1
                good_b = (ib >= 0) & (ib < MAG_BINS.size - 1)
                acc.update(
                    np.unique(
                        ib[good_b].astype(np.int64) * (1 << 32)
                        + c["src_id"][sel_mask][good_b].astype(np.int64)
                    ).tolist()
                )
        n_det += np.histogram(tm_ref[star & c["detected"]], MAG_BINS)[0]
        n_cls += np.histogram(tm_ref[star & c["classified"]], MAG_BINS)[0]
        n_gal_det += np.histogram(tm_ref[gal & c["detected"]], MAG_BINS)[0]
        n_gal_cls += np.histogram(tm_ref[gal & c["classified"]], MAG_BINS)[0]

        sel = star & c["classified"]
        resid_ref = c["obs_mag"][ref] - corr - tm_ref
        pe_sample.add(delta[sel], resid_ref[sel])
        sel_nc = star & c["classified_nosnr"]
        pe_sample_nc.add(delta[sel_nc], resid_ref[sel_nc])
        with np.errstate(all="ignore"):
            log_magerr = np.log10(c["obs_magerr"][ref])
            pe_catalog.add(delta[sel], log_magerr[sel])
            pe_catalog_nc.add(delta[sel_nc], log_magerr[sel_nc])
        for b in maps:
            # Which population defines "the depth" -- see --anchor-sample.
            cls_key = "classified" if args.anchor_sample == "detected" else "classified_nosnr"
            selb = c["is_star"] & c[cls_key] & keep & c["in_footprint"]
            anchor[b].add(
                c["true_mag"][b][selb],
                (c["obs_mag"][b] - corr - c["true_mag"][b])[selb],
            )
        print(f"  {hi:,}/{n:,}", end="\r")
    print()
    schema.close()
    if n_dropped:
        print(f"  dropped {n_dropped:,} rows in rejected tiles")

    # ---- truth anchor ----------------------------------------------------
    m5, shift = {}, {}
    for b in maps:
        if args.no_anchor:
            m5[b], shift[b] = med_raw[b], 0.0
        else:
            sc_a = anchor[b].scatter()
            cnt_a = anchor[b].counts()
            good = np.isfinite(sc_a) & (sc_a > 0) & (cnt_a >= 100)
            if good.sum() < 3:
                raise SystemExit(f"band {b}: too few bins to truth-anchor")
            m5[b] = float(
                np.interp(np.log10(SIG_SN5), np.log10(sc_a[good]), ANCHOR_MID[good])
            )
            shift[b] = m5[b] - med_raw[b]
            maps[b].offset = shift[b]
        print(
            f"  {b}: map median {med_raw[b]:.3f} -> truth-anchored {m5[b]:.3f} "
            f"(shift {shift[b]:+.3f})"
        )

    maglim_ref = med_raw[ref] + shift[ref]
    delta_mid = DELTA_MID - shift[ref]  # relabel the raw-map delta axis
    print(f"reference maglim ({ref}) = {maglim_ref:.3f}")

    # ---- efficiency ------------------------------------------------------
    with np.errstate(invalid="ignore", divide="ignore"):
        eff_det = n_det / n_all
        # classification_eff = cls/det, so it needs its OWN count guard: at the
        # faint end n_det collapses and the ratio turns up unphysically on a
        # handful of objects (seen: 0.29 -> 0.42 faintward of delta_mag ~1.75).
        eff_cls = np.where(n_det >= MIN_COUNT_EFF, n_cls / np.maximum(n_det, 1), np.nan)

    # ---- surrogate deconvolution ----------------------------------------
    # The surrogate S is not EXT_XGB, so the measured eff_S = P(S=1|star) is a
    # biased estimate of what we want, eff_X = P(EXT_XGB<=1|star).  With
    #     a = P(S=1 | EXT_XGB<=1),  b = P(S=1 | EXT_XGB>1)
    # measured on the real catalogue per magnitude bin,
    #     eff_S = a*eff_X + b*(1-eff_X)   ->   eff_X = (eff_S - b)/(a - b).
    # Assumes a,b are the same for true stars as for the mixed population they
    # were measured on -- second-order, and checked externally against SPLASH.
    deconv = {"applied": False}
    if args.confusion:
        eff_cls, deconv = deconvolve_classification_eff(
            eff_cls, pd.read_csv(args.confusion), MAG_MID)
        deconv["confusion"] = str(args.confusion)
        print(f"\n  deconvolved classification_eff in "
              f"{deconv['n_bins_corrected']} bins (median |change| "
              f"{deconv['median_abs_change']:.4f}, "
              f"max {deconv['max_abs_change']:.4f})")
    elif args.survey == "des_y6":
        print("\n  WARNING: no --confusion given; classification_eff describes "
              "the SURROGATE, not EXT_XGB")

    with np.errstate(invalid="ignore", divide="ignore"):
        # classification_detection_eff must stay consistent with the (possibly
        # deconvolved) classification_eff rather than being recomputed from the
        # raw counts, or the two curves would describe different classifiers.
        eff_both = eff_det * np.nan_to_num(eff_cls)
    eff = pd.DataFrame(
        {
            f"mag_{ref}": MAG_MID,
            "delta_mag": MAG_MID - maglim_ref,
            "detection_eff": eff_det,
            "classification_eff": eff_cls,
            "classification_detection_eff": eff_both,
        }
    )
    enough = n_all >= MIN_COUNT_EFF
    n_src = None
    if uniq_src:
        n_src = np.zeros(MAG_BINS.size - 1, dtype=np.int64)
        for packed in uniq_src:
            n_src[packed >> 32] += 1
        thin = enough & (n_src < MIN_UNIQUE_SRC)
        if thin.any():
            # CLAMP, do not drop.  streamobs's injector holds the loaded curve
            # flat at its first (brightest) row for anything brighter than
            # that -- so whatever value survives here becomes literally the
            # bright-edge value applied to every brighter magnitude.
            # Truncating the bright end instead of clamping would let a thin,
            # noisy bin become that held-flat value, i.e. assert that bright,
            # well-detected stars are nearly unrecoverable.  Instead hold
            # classification_eff at the brightest well-sampled value: DES does
            # not classify a g = 17 star worse than a g = 18.5 one.
            # detection_eff is left as measured; it sits at ~1.0 and does not
            # depend on morphology, so it is not affected by the corrupt
            # bright-end injection profiles.
            reliable = np.where(enough & (n_src >= MIN_UNIQUE_SRC))[0]
            if reliable.size:
                i0 = int(reliable[0])
                bright = np.zeros_like(enough)
                bright[:i0] = enough[:i0]
                eff_cls = eff_cls.copy()
                eff_cls[bright] = eff_cls[i0]
                eff["classification_eff"] = eff_cls
                eff["classification_detection_eff"] = eff_det * np.nan_to_num(eff_cls)
                print(f"  clamped classification_eff in {int(bright.sum())} bright "
                      f"bins (mag_{ref} <= {MAG_MID[i0 - 1]:.2f}) to the brightest "
                      f"well-sampled value {eff_cls[i0]:.4f}: fewer than "
                      f"{MIN_UNIQUE_SRC} distinct parent sources there, and the few "
                      f"bright deep-field stars are saturated in the deep imaging, "
                      f"so their injected morphology is not representative")
    eff = eff[enough].fillna(0.0)
    eff = eff[eff["delta_mag"] >= EFF_DELTA_MIN]
    faint = eff["delta_mag"] > DET_EFF_DELTA_MAX
    eff.loc[faint, ["detection_eff", "classification_detection_eff"]] = 0.0
    print(f"  zeroed detection_eff in {int(faint.sum())} bins (delta_mag > {DET_EFF_DELTA_MAX})")
    write_csv(
        out / f"{tag}_stellar_efficiency_cut{ref}.csv",
        eff,
        f"mag_{ref},delta_mag,detection_eff,classification_eff,classification_detection_eff",
    )

    # ---- photo-error, two curves ----------------------------------------
    sc = pe_sample.scatter()
    cnt = pe_sample.counts()
    med_rep = pe_catalog.median()
    keep = np.isfinite(sc) & (sc > 0) & (cnt >= MIN_COUNT_PE)
    sample_tab = pd.DataFrame(
        {"delta_mag": delta_mid[keep], "log_mag_err": np.log10(sc[keep])}
    )
    catalog_tab = pd.DataFrame(
        {"delta_mag": delta_mid[keep], "log_mag_err": med_rep[keep]}
    )
    near = (sample_tab["delta_mag"] > -3) & (sample_tab["delta_mag"] < 0.5)
    factor = 10 ** (
        sample_tab["log_mag_err"].values - catalog_tab["log_mag_err"].values
    )
    inflation = (
        float(np.nanmedian(factor[near.values])) if near.any() else float("nan")
    )
    print(f"\n  ERROR-INFLATION FACTOR (truth scatter / reported): {inflation:.2f}")
    print("  (~1 = reported errors are well calibrated; Roman DC2 was ~2)\n")

    # The forced-photometry pair, built exactly the same way off the no-S/N-cut
    # population.  streamobs applies these to every band that is not the
    # reference band; see Survey._resolve_log_photo_error.
    sc_nc = pe_sample_nc.scatter()
    cnt_nc = pe_sample_nc.counts()
    med_rep_nc = pe_catalog_nc.median()
    keep_nc = np.isfinite(sc_nc) & (sc_nc > 0) & (cnt_nc >= MIN_COUNT_PE)
    sample_nc_tab = pd.DataFrame(
        {"delta_mag": delta_mid[keep_nc], "log_mag_err": np.log10(sc_nc[keep_nc])}
    )
    catalog_nc_tab = pd.DataFrame(
        {"delta_mag": delta_mid[keep_nc], "log_mag_err": med_rep_nc[keep_nc]}
    )

    for name, tab, cid in [
        (f"{tag}_photoerror_{ref}", sample_tab, f"{ref}_sample"),
        (f"{tag}_photoerror_{ref}_catalog", catalog_tab, f"{ref}_catalog"),
        (f"{tag}_photoerror_{ref}_nocut", sample_nc_tab, f"{ref}_sample_nocut"),
        (f"{tag}_photoerror_{ref}_catalog_nocut", catalog_nc_tab, f"{ref}_catalog_nocut"),
    ]:
        write_csv(out / f"{name}_raw.csv", tab, "delta_mag,log_mag_err")
        write_csv(
            out / f"{name}.csv",
            apply_photoerr_corrections(tab, cid, args.corrections),
            "delta_mag,log_mag_err",
        )

    # ---- galaxy misclassification ---------------------------------------
    with np.errstate(invalid="ignore", divide="ignore"):
        misclass = np.where(n_gal_det > 0, n_gal_cls / np.maximum(n_gal_det, 1), np.nan)
    # Same thin-parent-sample hazard as the stellar curve, and just as visible:
    # unguarded, the bright end read 0.63 and 0.45 at mag_g ~ 15.4-15.6, i.e. half
    # of bright galaxies called stars.  Clamp (do not drop) for the same reason --
    # streamobs's injector holds the loaded curve flat at its first surviving
    # row for anything brighter.  Bright galaxies are well resolved, so holding
    # the rate at the brightest well-sampled value is also the physically
    # sensible shape.
    if uniq_src_gal:
        n_src_g = np.zeros(MAG_BINS.size - 1, dtype=np.int64)
        for packed in uniq_src_gal:
            n_src_g[packed >> 32] += 1
        ok_g = (n_gal_det >= MIN_COUNT_EFF) & (n_src_g >= MIN_UNIQUE_SRC)
        rel_g = np.where(ok_g)[0]
        if rel_g.size:
            j0 = int(rel_g[0])
            thin_g = (n_gal_det >= MIN_COUNT_EFF) & (np.arange(misclass.size) < j0)
            if thin_g.any():
                misclass = misclass.copy()
                misclass[thin_g] = misclass[j0]
                print(f"  clamped missclassification_eff in {int(thin_g.sum())} "
                      f"bright bins (mag_{ref} <= {MAG_MID[j0 - 1]:.2f}) to "
                      f"{misclass[j0]:.4f}: fewer than {MIN_UNIQUE_SRC} distinct "
                      f"parent galaxies there")
    mis = pd.DataFrame(
        {
            f"mag_{ref}": MAG_MID,
            "delta_mag": MAG_MID - maglim_ref,
            # MUST be `missclassification_eff` (sic).  That is the column
            # streamobs' set_completeness(selection="missclassified") reads, and
            # SurveyFactory wraps the load in a bare `except: pass`, so any other
            # spelling makes the product load as None with no error at all --
            # which is exactly the pre-existing des_yr6 failure mode this
            # product is meant to fix.
            "missclassification_eff": misclass,
        }
    )
    mis = mis[n_gal_det >= MIN_COUNT_EFF].fillna(0.0)
    # Untouched measurement first, then the tracked bright-end cut: brightward of
    # it the true-galaxy counts per bin get small and the rate swings between
    # zero and tens of percent, so it stops being a measurement. Same afterburner
    # YAML and same cut_bright rule as the photo-error curves; see
    # scripts/apply_misclass_cuts.py, which applies it to the releases whose
    # generators are too expensive to re-run.
    write_csv(
        out / f"{tag}_galaxy_misclass_cut{ref}_raw.csv",
        mis,
        f"mag_{ref},delta_mag,missclassification_eff",
    )
    write_csv(
        out / f"{tag}_galaxy_misclass_cut{ref}.csv",
        apply_photoerr_corrections(mis, f"{ref}_misclass", args.corrections),
        f"mag_{ref},delta_mag,missclassification_eff",
    )

    # ---- anchored depth maps --------------------------------------------
    if args.write_maglim and not args.no_anchor:
        import healpy as hp

        for b, mm in maps.items():
            nso = args.maglim_nside or mm.nside
            d = mm.to_healpix(nso)
            # Mask junk pixels inherited from the input depth map.  A handful
            # survive the (15, 30) sanity window while being physically
            # impossible -- before this, r reached 28.99 (4 mag deeper than its
            # own median) and 15.56.  They are rare (0.01-0.08% per band) but an
            # injection landing on one gets a nonsense delta_mag, so they are
            # dropped rather than shipped.
            if args.maglim_clip > 0:
                fin = np.isfinite(d)
                if fin.any():
                    med_b = float(np.median(d[fin]))
                    bad = fin & (np.abs(d - med_b) > args.maglim_clip)
                    if bad.any():
                        print(f"  {b}: masked {int(bad.sum())} pixels deviating "
                              f"> {args.maglim_clip} mag from the median "
                              f"{med_b:.3f} ({100 * bad.sum() / fin.sum():.4f}%)")
                        d = np.where(bad, np.nan, d)
            arr = np.where(np.isfinite(d), d, hp.UNSEEN)
            fn = out / f"{tag}_maglim_{b}_nside{nso}.fits.gz"
            hp.write_map(str(fn), arr, overwrite=True, dtype=np.float32)
            print(f"  wrote {fn} (median {np.nanmedian(d):.3f})")

    # ---- audit record ----------------------------------------------------
    audit = {
        "tag": tag,
        "survey": args.survey,
        "catalog": str(args.catalog),
        "rows_read": int(n),
        "ref_band": ref,
        "classifier": (
            "ext_xgb_surrogate" if args.survey == "des_y6"
            else "bdf_extended_class_dr3gold"
        ),
        "deconvolution": deconv,
        "ext_max": args.ext_max,
        "mag_kind": args.mag_kind,
        "snr_detect": args.snr_detect,
        "truth_anchored": not args.no_anchor,
        "anchor_sample": args.anchor_sample,
        "m5_truth_anchored": {b: float(v) for b, v in m5.items()},
        "maglim_ref_median": maglim_ref,
        "error_inflation_factor": inflation,
        "n_true_stars_binned": int(n_all.sum()),
        "n_detected": int(n_det.sum()),
        "n_classified": int(n_cls.sum()),
    }
    with open(out / f"{tag}_audit.json", "w") as fh:
        json.dump(audit, fh, indent=2)
    print(f"  wrote {out / f'{tag}_audit.json'}")
    print("\ndone — copy the CSVs (and the audit json) back; they are tiny.")


def build_parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--survey", required=True, choices=sorted(SCHEMAS))
    p.add_argument("--catalog", default=None, help="Balrog HDF5 catalog")
    p.add_argument(
        "--measured",
        default=None,
        help="des_y6 only: fiducial_matched_measured_sof.hdf5 (the meas_* file). "
        "DELVE ships one combined table and does not need this",
    )
    p.add_argument(
        "--truth-labels",
        default=None,
        help="des_y6 only: parquet with ID + KNN_CLASS from "
        "scripts/des/build_des_truth_labels.py",
    )
    p.add_argument(
        "--surrogate",
        default=None,
        help="des_y6 only: EXT_XGB surrogate json from "
        "scripts/des/build_des_xgb_surrogate.py (its _features.txt must sit "
        "alongside it)",
    )
    p.add_argument("--surrogate-threshold", type=float, default=0.5)
    p.add_argument(
        "--foreground-nside",
        type=int,
        default=4096,
        help="des_y6 only: nside of the POSITIONAL foreground mask built from "
        "meas_FLAGS_FOREGROUND (a pixel is masked if >50%% of its detections "
        "are flagged). The flag exists only for detected rows, so cutting on "
        "it directly would bias detection_eff low; 0 disables the mask and "
        "leaves the depth-map footprint as the only sky cut",
    )
    p.add_argument(
        "--confusion",
        default=None,
        help="des_y6 only: per-magnitude confusion CSV (mag_g,a,b) used to "
        "deconvolve the surrogate's selection back onto EXT_XGB. Without it "
        "classification_eff describes the SURROGATE, not EXT_XGB",
    )
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--tag", default=None, help="output filename prefix (default: --survey)")
    p.add_argument(
        "--maglim-map",
        nargs="+",
        default=None,
        metavar="BAND=PATH[,PATH...]",
        help="depth map(s) per band, e.g. g=/path/maglim_g.hsp. Give several "
        "comma-separated paths to mosaic complementary footprints (DELVE DR3 "
        "needs both dr32 and dr311+dr312); first listed wins on overlap",
    )
    p.add_argument("--band", default="g", help="reference band for the curves")
    p.add_argument(
        "--maglim-nside",
        type=int,
        default=1024,
        help="nside of the depth maps WRITTEN OUT by --write-maglim. Input maps "
        "are always sampled at their native resolution (HealSparse maps are "
        "queried sparsely, never densified); 0 writes at native nside",
    )
    p.add_argument("--chunk", type=int, default=2_000_000, help="rows per chunk")
    p.add_argument("--max-rows", type=int, default=0, help="stop after N rows (testing)")
    p.add_argument("--row-start", type=int, default=0, help="first row to read (testing)")
    p.add_argument(
        "--ext-max",
        type=int,
        default=1,
        choices=[0, 1, 2],
        help="star iff 0 <= EXT_BDF <= this. 0 = pure stellar sample, "
        "1 = complete stellar sample (the DES Y6 default), 2 = looser",
    )
    p.add_argument(
        "--mag-kind",
        default=None,
        choices=["psf", "bdf"],
        help="measured magnitude to use (default: per-survey -- psf for delve, "
        "bdf for des_y6, which matches its deep-field BDF truth magnitude)",
    )
    p.add_argument(
        "--truth-mag-kind",
        default="deredden",
        choices=["raw", "reddened", "deredden"],
        help="des_y6 only: how to build the true apparent magnitude. "
        "'deredden' = bdf_mag_deredden + ext_mags (default: de-redden at the "
        "deep field, re-redden at the injection position), "
        "'reddened' = bdf_mag + ext_mags, 'raw' = bdf_mag as delivered. "
        "The per-tile SPREAD of (obs-true) is the diagnostic, not its median: "
        "'raw' has the smallest median offset yet a 12x wider spread, because "
        "extinction varies tile to tile",
    )
    p.add_argument("--corrections", default=None, help="afterburner YAML")
    p.add_argument(
        "--reapply-corrections",
        action="store_true",
        help="do NOT re-derive: just re-run the afterburner over the *_raw.csv "
        "already in --out and rewrite the cleaned curves. Editing a cleanup "
        "threshold should not cost a multi-hour pass over the catalogue",
    )
    p.add_argument(
        "--snr-detect",
        type=float,
        default=SNR_DEPTH,
        help="require reference-band S/N above this for an object to count as "
        "detected (default 5, the streamobs convention -- the curves own the "
        "S/N cut and the injector must not re-apply it). 0 disables, keeping "
        "Balrog's raw detection flag",
    )
    p.add_argument(
        "--badpix-max",
        type=float,
        default=None,
        help="also require meas_badpix_frac < this (DELVE V4: median 0.14, so a "
        "cut at 0.1 keeps only ~24%% of detections -- use with care)",
    )
    p.add_argument(
        "--tile-zp",
        default="reject",
        choices=["correct", "reject", "none"],
        help="what to do with tiles whose zero point is off: drop them "
        "(default -- matches DES Y6 Balrog, which ships an explicit "
        "flags_bad_zp on 11.4%% of injections rather than correcting), correct "
        "them, or ignore the audit",
    )
    p.add_argument(
        "--tile-zp-max",
        type=float,
        default=0.05,
        help="a tile is 'bad' if its median (obs-true) deviates from the "
        "footprint median by more than this (mag)",
    )
    p.add_argument("--zp-mag-min", type=float, default=19.0, help="bright end of the ZP audit")
    p.add_argument("--zp-mag-max", type=float, default=22.0, help="faint end of the ZP audit")
    p.add_argument(
        "--anchor-sample",
        default="nosnr",
        choices=["nosnr", "detected"],
        help="which population defines the depth. 'nosnr' (default) drops the "
        "reference-band S/N cut from the anchor sample only. 'detected' anchors "
        "on the same sample the photo-error curve describes, which would make "
        "sigma = 0.2171 at delta_mag = 0 hold by construction -- but MEASURED ON "
        "DES IT IS UNSTABLE ACROSS BANDS (see the module docstring); use it only "
        "for diagnostics",
    )
    p.add_argument(
        "--no-anchor",
        action="store_true",
        help="trust the input depth map's absolute scale (skip pass 1)",
    )
    p.add_argument("--write-maglim", action="store_true", help="emit anchored depth maps")
    p.add_argument(
        "--maglim-clip",
        type=float,
        default=1.5,
        help="mask written-out depth pixels deviating more than this many mag "
        "from their band median. The input healsparse maps carry a small number "
        "of physically impossible pixels that pass the (15, 30) window; 0 "
        "disables the mask and ships them",
    )
    return p


if __name__ == "__main__":
    main(build_parser().parse_args())
