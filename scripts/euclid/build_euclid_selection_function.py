#!/usr/bin/env python
"""Build streamobs selection-function products for Euclid from the Q1 MER catalogue.

Two inputs, both already on disk here:

  * the Euclid Q1 MER Final Catalog as a HATS collection
    (/astro/store/shire/hats/catalogs/euclid_q1, 29.8M rows, 63 deg^2 of the
    three Euclid Deep Fields observed to nominal Wide-survey depth), read with
    lsdb and reduced to the point-like sample;
  * the ECDFS spectroscopic truth compilation of Gatto et al. 2026 (A&A 709,
    A79; ``catalog.dat.gz``: 15,611 sources, 4,361 stars / 11,250 galaxies,
    with LSST DP1 ugrizy photometry), which supplies the star/galaxy labels.

Products (``streamobs/docs/source/selection_function_methodology.md``):

  <out>/euclid_q1_maglim_<band>_nside128.fits.gz   depth maps, VIS Y J H
  <out>/euclid_q1_stellar_efficiency_cutvis.csv
        mag_VIS, delta_mag, detection_eff, classification_eff,
        classification_detection_eff
  <out>/euclid_q1_photoerror_vis_catalog.csv       CATALOG (reported magerr)
  <out>/euclid_q1_photoerror_vis.csv               SAMPLE  (noise draw)
  <out>/euclid_q1_photoerror_vis_catalog_nocut.csv CATALOG, no S/N cut
  <out>/euclid_q1_photoerror_vis_nocut.csv         SAMPLE,  no S/N cut
  <out>/euclid_q1_galaxy_misclass_cutvis.csv       stellar contamination
  <out>/euclid_q1_audit.json                       counts, anchors, factors

Photometry
----------
The reference band is VIS (I_E), measured with ``FLUX_VIS_PSF``: it is the
point-source estimator, it is total (the 2xFWHM aperture flux is NOT
aperture-corrected and sits 0.09 mag fainter for stars), and it is what the
Euclid x DES stellar work in this group already uses. The NISP bands use
``FLUX_{Y,J,H}_TEMPLFIT``, the template-fit flux measured with the VIS
detection as prior -- i.e. forced photometry at the VIS position, which is
exactly the role streamobs assigns to every non-reference band. The aperture
NIR fluxes are also uncorrected (0.12-0.15 mag fainter than TEMPLFIT) and
their reported S/N = 5 lands ~1.2 mag shallower, so they are not used.

Depth convention: externally anchored
-------------------------------------
The catalogue's own reported PSF errors put the median S/N = 5 point-source
depth at I_E ~ 27.2, a full magnitude deeper than the 26.2 the Euclid Wide
Survey delivers for 5 sigma point sources (Scaramella et al. 2022), and the
reported-error slope is 0.27 dex/mag where background-limited photometry gives
0.4 -- the signature of errors that are increasingly underestimated faintward,
the same situation as the Roman DC2 mock (methodology doc, *Depth maps*).
There the fix is to truth-anchor the map to the scatter of (obs - true). Euclid
has no injection catalogue and the truth compilation carries no Euclid-band
photometry (a colour transform from DP1/DP2 ugriz floors at ~0.1 mag, useless
for this), so the map is anchored to the *published* depth instead: the
reported-error map supplies the spatial structure, exactly as the DES
healsparse maps do, and its median is shifted onto the Wide-survey 5 sigma
point-source depth per band. The implied error-inflation factor
``f = 10**(0.4 * (native - anchored))`` plays the role Roman's truth-measured
factor plays:

  * "detected" means reported S/N > 5 f, i.e. true S/N > 5 under the anchor;
  * the SAMPLE photo-error curve is the CATALOG curve scaled by f, so that
    sigma_sample reaches 2.5/ln10/5 = 0.2171 at delta_mag = 0 -- the property
    truth anchoring enforces by construction. A constant factor is the
    minimal assumption; it is not a measurement, and the audit records it.

``--no-anchor`` ships the native reported-error scale instead (f = 1).

Truth-based curves
------------------
Stars and galaxies from the ECDFS compilation are matched to Q1 at 1". A truth
star is *detected* when it has a counterpart with a valid PSF flux, passes
SPURIOUS_FLAG == 0 and has reported S/N > 5 f; it is *classified* when
POINT_LIKE_PROB > 0.5, the MER point-source probability (96% complete / 96%
pure on this sample, and the best of the Euclid classifiers tested here).
DET_QUALITY_FLAG is deliberately not used: it is set on 94% of stars brighter
than I_E = 19 (bright-star bits), so it is not an observer's point-source
quality cut. The magnitude axis is the observed PSF I_E; the few truth stars
without one (undetected, or no PSF flux) fall back to a DP1 (r, i, z) colour
transform fitted on the matched stars, so they stay in the denominator.

Run with the streamobs env (lsdb 0.9.2 / hats 0.9.2):
    /astro/store/shiren/conda-envs/stream_team/envs/streamobs/bin/python \
        scripts/euclid/build_euclid_selection_function.py --figdir docs/source/_static/euclid_q1
"""

from __future__ import annotations

import argparse
import json
import pathlib
import time

import healpy as hp
import numpy as np
import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[2]
OUT_DEFAULT = REPO / "data/surveys/euclid_q1"
HATS_EUCLID = "/astro/store/shire/hats/catalogs/euclid_q1"
TRUTH_DEFAULT = (
    "/astro/store/shire/pferguso/projects/euclid/euclid_x_des/catalog.dat.gz"
)
TAG = "euclid_q1"
REF = "VIS"
BANDS = ["VIS", "Y", "J", "H"]
FLUX = {
    "VIS": "FLUX_VIS_PSF",
    "Y": "FLUX_Y_TEMPLFIT",
    "J": "FLUX_J_TEMPLFIT",
    "H": "FLUX_H_TEMPLFIT",
}
# Euclid Wide Survey 5 sigma point-source depths (Scaramella et al. 2022,
# A&A 662, A112). Q1 observed the EDFs to nominal Wide depth (Euclid Q1 overview,
# arXiv:2503.15302), so these are the fiducial depths for Q1 and for DR1.
EXTERNAL_M5 = {"VIS": 26.2, "Y": 24.5, "J": 24.5, "H": 24.5}
# VIS: POINT_LIKE_PROB is undefined for 96% of sources brighter than 17.0 and
# the point-like fraction collapses 0.90 -> 0.64 -> 0.03 across 17.5-17.0-16.5.
# NISP: TEMPLFIT S/N still rises down to 16.0 (no saturation signature); the
# brightest point-like stars in the catalogue sit at ~15.4, so 15.5 is a floor
# rather than a measurement.
SATURATION = {"VIS": 17.5, "Y": 15.5, "J": 15.5, "H": 15.5}
# pivot wavelengths (Angstrom) from the SVO Filter Profile Service,
# Euclid/VIS.vis and Euclid/NISP.{Y,J,H}
PIVOT_AA = {"VIS": 7103.37, "Y": 10785.39, "J": 13620.63, "H": 17648.80}
# Schlafly & Finkbeiner 2011 recalibration of SFD98 E(B-V); the des/yr6 and
# delve/dr3_gold coefficients are on this footing (their Table 6 = 0.86 x F99)
SF11_FACTOR = 0.86

SNR_DEPTH = 5.0
SIG_SN5 = 2.5 / np.log(10) / SNR_DEPTH  # 0.21715 mag
INT_NULL = 32767  # MER int16 null
NSIDE = 128
MIN_PIX_STARS = 20  # point-like sources needed for a pixel depth
MAGLIM_CLIP = 1.5  # mask pixels this far from the band median
MATCH_RADIUS = 1.0  # arcsec
MIN_COUNT_EFF = 20
MIN_COUNT_MIS = 100
MIN_COUNT_PE = 20
DET_EFF_DELTA_MAX = 1.0
MAG_BINS = np.arange(17.5, 27.0 + 1e-6, 0.5)
MAG_MID = 0.5 * (MAG_BINS[1:] + MAG_BINS[:-1])
DELTA_BINS = np.arange(-10.0, 2.0 + 1e-6, 0.12)
DELTA_MID = 0.5 * (DELTA_BINS[1:] + DELTA_BINS[:-1])
# truth-catalogue layout: RA, Dec, then (mag, err) pairs in ugrizy order,
# then is_star, is_galaxy, originating survey
T_RA, T_DEC, T_R, T_RERR, T_I, T_IERR, T_Z, T_ZERR, T_STAR = (
    0,
    1,
    6,
    7,
    8,
    9,
    10,
    11,
    14,
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def num(df, col):
    return pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)


def flux_to_mag(flux):
    """MER fluxes are in uJy: m_AB = -2.5 log10(f) + 23.9."""
    flux = np.asarray(flux, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(
            flux > 0, -2.5 * np.log10(np.where(flux > 0, flux, 1)) + 23.9, np.nan
        )


def mag_and_err(df, band):
    f, e = num(df, FLUX[band]), num(df, FLUX[band].replace("FLUX_", "FLUXERR_"))
    ok = np.isfinite(f) & (f > 0) & np.isfinite(e) & (e > 0) & (f < INT_NULL)
    m = np.where(ok, flux_to_mag(np.where(ok, f, 1)), np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        err = np.where(ok, 1.0857362 * e / np.where(ok, f, 1), np.nan)
    return m, err


def f99_a_over_ebv(lambda_aa, rv=3.1):
    """Fitzpatrick (1999) A_lambda / E(B-V) for R_V = 3.1 at lambda (Angstrom).

    The optical/IR part of F99 is a cubic spline through fixed anchor points
    in x = 1/lambda (um^-1); these are the Table 3 anchors for R_V = 3.1.
    Valid for lambda > 3700 A, which covers every Euclid band.
    """
    from scipy.interpolate import CubicSpline

    if rv != 3.1:
        raise ValueError("only the R_V = 3.1 anchors are tabulated here")
    x_anchor = np.array([0.0, 0.377, 0.820, 1.667, 1.828, 2.141, 2.433])
    k_anchor = np.array([0.0, 0.265, 0.829, 2.688, 3.055, 3.806, 4.315])
    spl = CubicSpline(x_anchor, k_anchor)
    x = 1e4 / np.asarray(lambda_aa, dtype=float)
    if np.any(x > x_anchor[-1]):
        raise ValueError("wavelength outside the optical/IR spline domain")
    return float(spl(x)) if np.ndim(x) == 0 else spl(x)


def wilson(k, n, z=1.0):
    k, n = np.asarray(k, float), np.asarray(n, float)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = k / n
        denom = 1 + z**2 / n
        centre = (p + z**2 / (2 * n)) / denom
        half = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2)) / denom
    return p, np.clip(centre - half, 0, 1), np.clip(centre + half, 0, 1)


def write_csv(path, df, header):
    np.savetxt(path, df.values, delimiter=",", header=header, fmt="%.6f", comments="")
    print(f"  wrote {path.name} ({len(df)} rows)")


def truth_cone(ra, dec, pad_deg=0.15):
    ra_c, dec_c = 0.5 * (ra.min() + ra.max()), 0.5 * (dec.min() + dec.max())
    r = np.hypot((ra - ra_c) * np.cos(np.radians(dec)), dec - dec_c).max()
    return ra_c, dec_c, r + pad_deg


# ---------------------------------------------------------------------------
# depth maps
# ---------------------------------------------------------------------------
def error_slope(mag, err, sat, snr_lo=10.0, snr_hi=100.0):
    """Measured slope of log10(magerr) vs mag over 10 < S/N < 100, for the audit.

    Background-limited photometry gives 0.4 dex/mag. The Q1 PSF errors give
    ~0.2-0.27 here, i.e. the reported errors grow more slowly than the sky
    noise allows -- increasingly underestimated faintward. That is why the
    depth is NOT obtained by extrapolating along this slope (doing so puts the
    5 sigma depth at I_E ~ 29), and why the absolute scale is anchored
    externally instead.
    """
    snr = 1.0857362 / err
    sel = (
        np.isfinite(mag)
        & np.isfinite(err)
        & (snr > snr_lo)
        & (snr < snr_hi)
        & (mag > sat + 0.5)
    )
    bins = np.arange(np.floor(mag[sel].min() * 10) / 10, mag[sel].max() + 0.1, 0.1)
    idx = np.digitize(mag[sel], bins) - 1
    xs, ys = [], []
    for i in range(len(bins) - 1):
        m = idx == i
        if m.sum() >= 50:
            xs.append(0.5 * (bins[i] + bins[i + 1]))
            ys.append(np.median(np.log10(err[sel][m])))
    slope, _ = np.polyfit(xs, ys, 1)
    return float(slope), int(sel.sum())


def depth_map(mag, err, ra, dec, sat, nside=NSIDE, snr_lo=20.0, snr_hi=60.0):
    """Reported-error S/N = 5 depth per RING pixel, background-limited scaling.

    Each object at 20 < S/N < 60 is carried to S/N = 5 with the background-
    limited relation m5 = m + 2.5 log10(S/N / 5), and the pixel takes the
    median. The S/N window keeps the extrapolation short (1.5-2.7 mag) and
    stays where the point-like selection is still pure (I_E < 25); the
    catalogue's own faint-end error slope is not trusted, see error_slope.
    This supplies the *spatial structure* of the map; its absolute scale is
    replaced by the anchor in main().
    """
    snr = 1.0857362 / err
    ok = (
        np.isfinite(mag)
        & np.isfinite(err)
        & (snr > snr_lo)
        & (snr < snr_hi)
        & (mag > sat + 1.0)
    )
    m5 = mag[ok] + 2.5 * np.log10(snr[ok] / SNR_DEPTH)
    pix = hp.ang2pix(nside, ra[ok], dec[ok], lonlat=True)
    npix = hp.nside2npix(nside)
    order = np.argsort(pix)
    pix_s, m5_s = pix[order], m5[order]
    uniq, start, count = np.unique(pix_s, return_index=True, return_counts=True)
    out = np.full(npix, np.nan)
    for p, s, c in zip(uniq, start, count):
        if c >= MIN_PIX_STARS:
            out[p] = np.median(m5_s[s : s + c])
    return out, int(ok.sum())


def write_map(path, dense):
    arr = np.where(np.isfinite(dense), dense, hp.UNSEEN).astype(np.float32)
    hp.write_map(str(path), arr, overwrite=True, dtype=np.float32)
    good = np.isfinite(dense)
    print(
        f"  wrote {path.name}: median {np.nanmedian(dense):.3f}, "
        f"{good.sum()} pixels = {good.sum() * hp.nside2pixarea(hp.npix2nside(dense.size), degrees=True):.1f} deg2"
    )


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main(args):
    import lsdb

    t0 = time.time()
    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    figdir = pathlib.Path(args.figdir) if args.figdir else None
    if figdir:
        figdir.mkdir(parents=True, exist_ok=True)

    # ---- 1. the point-like Q1 sample -------------------------------------
    cols = ["RIGHT_ASCENSION", "DECLINATION", "POINT_LIKE_PROB", "SPURIOUS_FLAG"]
    for b in BANDS:
        cols += [FLUX[b], FLUX[b].replace("FLUX_", "FLUXERR_")]
    cat = lsdb.open_catalog(args.hats, columns=cols)
    pl = cat.query(f"POINT_LIKE_PROB > 0.5 and POINT_LIKE_PROB < {INT_NULL}").compute()
    pl = pl[num(pl, "SPURIOUS_FLAG") == 0].reset_index(drop=True)
    ra_pl, dec_pl = num(pl, "RIGHT_ASCENSION"), num(pl, "DECLINATION")
    print(
        f"point-like, SPURIOUS_FLAG == 0 sample: {len(pl):,} rows  [{time.time() - t0:.0f} s]"
    )

    # ---- 2. depth maps + anchor -----------------------------------------
    maps_native, maps_anchored, depth = {}, {}, {}
    for b in BANDS:
        m, e = mag_and_err(pl, b)
        slope, n_slope = error_slope(m, e, SATURATION[b])
        native, n_fit = depth_map(m, e, ra_pl, dec_pl, SATURATION[b])
        med_native = float(np.nanmedian(native))
        if args.no_anchor:
            shift, m5 = 0.0, med_native
        else:
            m5 = EXTERNAL_M5[b]
            shift = m5 - med_native
        anchored = native + shift
        fin = np.isfinite(anchored)
        bad = fin & (np.abs(anchored - np.nanmedian(anchored)) > MAGLIM_CLIP)
        anchored[bad] = np.nan
        f_bg = 10 ** (0.4 * (med_native - m5))
        maps_native[b], maps_anchored[b] = native, anchored
        depth[b] = dict(
            reported_error_slope_dex_per_mag=slope,
            n_slope=n_slope,
            n_depth=n_fit,
            native_median_bg_extrapolation=med_native,
            anchored_median=float(np.nanmedian(anchored)),
            external_m5=EXTERNAL_M5[b],
            shift=shift,
            inflation_factor_bg_extrapolation=f_bg,
            n_pixels=int(fin.sum()),
            n_pixels_clipped=int(bad.sum()),
            area_deg2=float(
                np.isfinite(anchored).sum() * hp.nside2pixarea(NSIDE, degrees=True)
            ),
        )
        print(
            f"  {b:<3} error slope {slope:.3f} dex/mag; native median ({n_fit:,} objects) "
            f"{med_native:.3f} -> anchored {m5:.3f} (shift {shift:+.3f}, bg-extrapolation f = {f_bg:.2f}), "
            f"{int(fin.sum())} pixels, {int(bad.sum())} clipped"
        )
        write_map(out / f"{TAG}_maglim_{b.lower()}_nside{NSIDE}.fits.gz", anchored)
    # ---- 3. photo-error curves (reference band, point-like sample) --------
    m_ref, e_ref = mag_and_err(pl, REF)
    pix = hp.ang2pix(NSIDE, ra_pl, dec_pl, lonlat=True)
    ml = maps_anchored[REF][pix]
    ok = (
        np.isfinite(m_ref)
        & np.isfinite(e_ref)
        & np.isfinite(ml)
        & (m_ref >= SATURATION[REF])
    )
    delta = m_ref - ml
    snr = 1.0857362 / e_ref
    log_err = np.log10(e_ref)

    def median_curve(sel):
        idx = np.digitize(delta[sel], DELTA_BINS) - 1
        good = (idx >= 0) & (idx < DELTA_MID.size)
        med = np.full(DELTA_MID.size, np.nan)
        cnt = np.bincount(idx[good], minlength=DELTA_MID.size)
        for i in np.where(cnt >= MIN_COUNT_PE)[0]:
            med[i] = np.median(log_err[sel][good][idx[good] == i])
        return med, cnt

    # The inflation factor is read off the no-cut curve AT the anchored depth:
    # f = 0.2171 / median reported sigma of point-like sources at delta_mag = 0.
    # It is what the anchor implies about the reported errors there, with no
    # assumption about how they scale with magnitude (the catalogue's own
    # error slope is 0.2 dex/mag, not the background-limited 0.4, so an
    # extrapolated "native" depth depends on where one extrapolates from --
    # 26.6 from S/N 20-60, 27.2 from the reported S/N = 5 crossing).
    nc_med, nc_cnt = median_curve(ok)
    keep_nc = np.isfinite(nc_med)
    sig_rep0 = float(10 ** np.interp(0.0, DELTA_MID[keep_nc], nc_med[keep_nc]))
    f_ref = 1.0 if args.no_anchor else SIG_SN5 / sig_rep0
    snr_cut = SNR_DEPTH * f_ref
    depth[REF]["inflation_factor"] = f_ref
    depth[REF]["reported_sigma_at_anchored_depth"] = sig_rep0
    print(
        f"reference band {REF}: median reported sigma at the anchored depth {sig_rep0:.4f} -> "
        f"f = {f_ref:.2f}; detection = reported S/N > {snr_cut:.2f}"
    )
    det = ok & (snr > snr_cut)
    cat_med, cat_cnt = median_curve(det)
    keep = np.isfinite(cat_med)
    logf = np.log10(f_ref)
    catalog_tab = pd.DataFrame(
        {"delta_mag": DELTA_MID[keep], "log_mag_err": cat_med[keep]}
    )
    sample_tab = catalog_tab.assign(log_mag_err=catalog_tab["log_mag_err"] + logf)
    catalog_nc_tab = pd.DataFrame(
        {"delta_mag": DELTA_MID[keep_nc], "log_mag_err": nc_med[keep_nc]}
    )
    sample_nc_tab = catalog_nc_tab.assign(
        log_mag_err=catalog_nc_tab["log_mag_err"] + logf
    )
    for name, tab in [
        (f"{TAG}_photoerror_{REF.lower()}_catalog", catalog_tab),
        (f"{TAG}_photoerror_{REF.lower()}", sample_tab),
        (f"{TAG}_photoerror_{REF.lower()}_catalog_nocut", catalog_nc_tab),
        (f"{TAG}_photoerror_{REF.lower()}_nocut", sample_nc_tab),
    ]:
        write_csv(out / f"{name}.csv", tab, "delta_mag,log_mag_err")
    # where the reported-error (catalog) curve crosses S/N = 10 -- the number
    # tests/test_surveys.py needs as snr10_delta_mag
    cc = catalog_tab.sort_values("log_mag_err")
    snr10_delta = float(
        np.interp(np.log10(1.0857362 / 10.0), cc["log_mag_err"], cc["delta_mag"])
    )
    sig0_cat = float(
        10 ** np.interp(0.0, catalog_tab["delta_mag"], catalog_tab["log_mag_err"])
    )
    sig0_sample = float(
        10 ** np.interp(0.0, sample_tab["delta_mag"], sample_tab["log_mag_err"])
    )
    print(
        f"  catalog curve: sigma({sig0_cat:.4f}) at delta_mag = 0, S/N = 10 at delta_mag {snr10_delta:+.3f}; "
        f"sample curve sigma({sig0_sample:.4f}) at delta_mag = 0"
    )

    # ---- 4. truth-based efficiency + misclassification --------------------
    raw = pd.read_fwf(args.truth, header=None, colspecs="infer", compression="gzip")
    t_ra, t_dec = raw[T_RA].to_numpy(float), raw[T_DEC].to_numpy(float)
    is_star = raw[T_STAR].to_numpy().astype(bool)
    r, i_, z = (raw[k].to_numpy(float) for k in (T_R, T_I, T_Z))
    rerr, ierr, zerr = (raw[k].to_numpy(float) for k in (T_RERR, T_IERR, T_ZERR))
    ra_c, dec_c, radius = truth_cone(t_ra, t_dec)
    full = (
        lsdb.open_catalog(
            args.hats,
            columns=[
                "RIGHT_ASCENSION",
                "DECLINATION",
                "POINT_LIKE_PROB",
                "SPURIOUS_FLAG",
                FLUX[REF],
                FLUX[REF].replace("FLUX_", "FLUXERR_"),
            ],
        )
        .cone_search(ra_c, dec_c, radius * 3600)
        .compute()
    )
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    idx, sep, _ = SkyCoord(t_ra * u.deg, t_dec * u.deg).match_to_catalog_sky(
        SkyCoord(num(full, "RIGHT_ASCENSION") * u.deg, num(full, "DECLINATION") * u.deg)
    )
    hit = sep.arcsec < MATCH_RADIUS
    mt = full.iloc[idx].reset_index(drop=True)
    m_obs, e_obs = mag_and_err(mt, REF)
    m_obs, e_obs = np.where(hit, m_obs, np.nan), np.where(hit, e_obs, np.nan)
    plp = num(mt, "POINT_LIKE_PROB")
    plp = np.where(hit & (plp < INT_NULL), plp, np.nan)
    quality = hit & (num(mt, "SPURIOUS_FLAG") == 0)
    snr_obs = 1.0857362 / e_obs
    detected = quality & np.isfinite(m_obs) & (snr_obs > snr_cut)
    classified = detected & (plp > 0.5)
    print(
        f"truth: {len(t_ra):,} sources ({is_star.sum():,} stars); matched {hit.mean():.3f}; "
        f"stars detected {detected[is_star].mean():.3f}, classified {classified[is_star].mean():.3f}; "
        f"Q1 rows in cone {len(full):,}"
    )

    # DP1 colour transform for the truth stars without a Euclid PSF magnitude
    fit = is_star & detected & (snr_obs > 50) & (m_obs > SATURATION[REF] + 0.5)
    for arr in (r, i_, z):
        fit &= np.isfinite(arr)
    for arr in (rerr, ierr, zerr):
        fit &= np.isfinite(arr) & (arr < 0.05)
    X = np.column_stack([np.ones(fit.sum()), (r - i_)[fit], (i_ - z)[fit]])
    coef, *_ = np.linalg.lstsq(X, (m_obs - i_)[fit], rcond=None)
    pred = i_ + coef[0] + coef[1] * (r - i_) + coef[2] * (i_ - z)
    resid = (m_obs - pred)[fit]
    transform_scatter = float((np.percentile(resid, 84) - np.percentile(resid, 16)) / 2)
    mag_axis = np.where(np.isfinite(m_obs) & hit, m_obs, pred)
    n_fallback = int((is_star & ~(np.isfinite(m_obs) & hit) & np.isfinite(pred)).sum())
    n_lost = int((is_star & ~np.isfinite(mag_axis)).sum())
    print(
        f"  DP1 transform I_E - i = {coef[0]:+.3f} {coef[1]:+.3f}(r-i) {coef[2]:+.3f}(i-z), "
        f"scatter {transform_scatter:.3f} mag on {fit.sum()} stars; used for {n_fallback} "
        f"undetected stars, {n_lost} stars have no magnitude at all and are dropped"
    )

    # local depth: the anchored VIS map over the truth cone
    vec = hp.ang2vec(ra_c, dec_c, lonlat=True)
    cone_pix = hp.query_disc(NSIDE, vec, np.radians(radius), inclusive=True)
    maglim_ref = float(np.nanmedian(maps_anchored[REF][cone_pix]))
    print(f"  anchored {REF} depth over the ECDFS cone: {maglim_ref:.3f}")

    star = is_star & np.isfinite(mag_axis)
    n_all = np.histogram(mag_axis[star], MAG_BINS)[0]
    n_det = np.histogram(mag_axis[star & detected], MAG_BINS)[0]
    n_cls = np.histogram(mag_axis[star & classified], MAG_BINS)[0]
    with np.errstate(invalid="ignore", divide="ignore"):
        eff_det = n_det / n_all
        eff_cls = np.where(n_det >= MIN_COUNT_EFF, n_cls / np.maximum(n_det, 1), np.nan)
        eff_both = eff_det * np.nan_to_num(eff_cls)
    eff = pd.DataFrame(
        {
            f"mag_{REF}": MAG_MID,
            "delta_mag": MAG_MID - maglim_ref,
            "detection_eff": eff_det,
            "classification_eff": eff_cls,
            "classification_detection_eff": eff_both,
            "n_true_stars": n_all,
        }
    )
    eff = eff[(n_all >= MIN_COUNT_EFF) & (MAG_MID >= SATURATION[REF])].fillna(0.0)
    # Faint-end extension. The truth compilation runs out of stars at about the
    # depth itself (fewer than MIN_COUNT_EFF per bin faintward of ~I_E 26), and
    # streamobs fills zero beyond a table's last row, which would turn a
    # measured ~0.9 into 0 across one bin. The bins between the last measured
    # one and delta_mag = +1 are therefore filled with what a hard S/N > 5 cut
    # does to Gaussian flux noise, P = Phi(5 (10^(-0.4 delta) - 1)), scaled by
    # the measured detection plateau; classification_eff is held at its last
    # measured value there. These rows are a MODEL, not a measurement -- the
    # audit records how many -- and they are gone as soon as a deeper truth
    # sample exists.
    from scipy.stats import norm

    plateau = float(
        np.median(eff["detection_eff"].to_numpy()[eff["delta_mag"].to_numpy() < -1.0])
    )
    last_mag = float(eff[f"mag_{REF}"].max())
    ext_mid = MAG_MID[
        (MAG_MID > last_mag) & (MAG_MID - maglim_ref <= DET_EFF_DELTA_MAX)
    ]
    n_model = int(ext_mid.size)
    if n_model:
        d = ext_mid - maglim_ref
        p_det = plateau * norm.cdf(SNR_DEPTH * (10 ** (-0.4 * d) - 1.0))
        cls_last = float(eff["classification_eff"].iloc[-1])
        ext = pd.DataFrame(
            {
                f"mag_{REF}": ext_mid,
                "delta_mag": d,
                "detection_eff": p_det,
                "classification_eff": cls_last,
                "classification_detection_eff": p_det * cls_last,
                "n_true_stars": 0,
            }
        )
        eff = pd.concat([eff, ext], ignore_index=True)
        print(
            f"  detection curve: {n_model} faint bin(s) beyond mag_{REF} = {last_mag:.2f} are the "
            f"S/N-threshold model (plateau {plateau:.3f}, classification held at {cls_last:.3f})"
        )
    faint = eff["delta_mag"] > DET_EFF_DELTA_MAX
    eff.loc[faint, ["detection_eff", "classification_detection_eff"]] = 0.0
    eff_out = eff.drop(columns="n_true_stars")
    write_csv(
        out / f"{TAG}_stellar_efficiency_cut{REF.lower()}.csv",
        eff_out,
        f"mag_{REF},delta_mag,detection_eff,classification_eff,classification_detection_eff",
    )

    gal = ~is_star & detected
    n_gal_det = np.histogram(m_obs[gal], MAG_BINS)[0]
    n_gal_cls = np.histogram(m_obs[gal & (plp > 0.5)], MAG_BINS)[0]
    with np.errstate(invalid="ignore", divide="ignore"):
        misclass = n_gal_cls / np.maximum(n_gal_det, 1)
    mis = pd.DataFrame(
        {
            f"mag_{REF}": MAG_MID,
            "delta_mag": MAG_MID - maglim_ref,
            # spelled with two s's: the column SurveyFactory reads
            "missclassification_eff": misclass,
            "n_true_galaxies": n_gal_det,
        }
    )
    mis = mis[n_gal_det >= MIN_COUNT_MIS]
    write_csv(
        out / f"{TAG}_galaxy_misclass_cut{REF.lower()}.csv",
        mis.drop(columns="n_true_galaxies"),
        f"mag_{REF},delta_mag,missclassification_eff",
    )

    # ---- 5. extinction coefficients -------------------------------------
    ext = {}
    for b in BANDS:
        raw_f99 = f99_a_over_ebv(PIVOT_AA[b])
        ext[b] = dict(
            pivot_aa=PIVOT_AA[b], f99_rv31=raw_f99, adopted_sf11=raw_f99 * SF11_FACTOR
        )
        print(
            f"  A_{b}/E(B-V): F99 at {PIVOT_AA[b]:.0f} A = {raw_f99:.3f}, x0.86 = {raw_f99 * SF11_FACTOR:.3f}"
        )

    # ---- 6. audit ---------------------------------------------------------
    audit = {
        "tag": TAG,
        "catalog": args.hats,
        "truth": args.truth,
        "ref_band": REF,
        "bands": BANDS,
        "flux_columns": FLUX,
        "quality_cut": "SPURIOUS_FLAG == 0",
        "classifier": "POINT_LIKE_PROB > 0.5",
        "depth_convention": (
            "native" if args.no_anchor else "external (Wide 5 sigma point source)"
        ),
        "detection_snr_cut_reported": snr_cut,
        "n_pointlike": int(len(pl)),
        "depth": depth,
        "m5_anchored": {b: depth[b]["anchored_median"] for b in BANDS},
        "saturation": SATURATION,
        "maglim_ref_truth_cone": maglim_ref,
        "truth_cone": dict(ra=ra_c, dec=dec_c, radius_deg=radius),
        "n_truth": int(len(t_ra)),
        "n_true_stars": int(is_star.sum()),
        "n_true_stars_binned": int(eff["n_true_stars"].sum()),
        "n_model_extended_bins": n_model,
        "detection_plateau": plateau,
        "n_stars_detected": int((star & detected).sum()),
        "n_stars_classified": int((star & classified).sum()),
        "n_true_galaxies_detected": int(gal.sum()),
        "dp1_transform": dict(
            coef=[float(c) for c in coef],
            scatter=transform_scatter,
            n_fit=int(fit.sum()),
            n_fallback=n_fallback,
            n_dropped=n_lost,
        ),
        "catalog_sigma_at_delta0": sig0_cat,
        "sample_sigma_at_delta0": sig0_sample,
        "snr10_delta_mag": snr10_delta,
        "extinction": ext,
    }
    with open(out / f"{TAG}_audit.json", "w") as fh:
        json.dump(audit, fh, indent=2)
    print(f"  wrote {TAG}_audit.json")

    # ---- 7. figures -------------------------------------------------------
    if figdir:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # (a) truth-based efficiency and misclassification
        fig, ax = plt.subplots(figsize=(8.5, 5.5))
        x = eff[f"mag_{REF}"].to_numpy()
        for key, lab, col in (
            ("detection_eff", f"detection (reported S/N > {snr_cut:.1f})", "#1f77b4"),
            (
                "classification_eff",
                "classification (POINT_LIKE_PROB > 0.5 | detected)",
                "#d62728",
            ),
            ("classification_detection_eff", "detection x classification", "k"),
        ):
            n = eff["n_true_stars"].to_numpy().astype(float)
            if key == "classification_eff":
                n = n * eff["detection_eff"].to_numpy()
            p = eff[key].to_numpy()
            _, lo, hi = wilson(p * np.maximum(n, 1), np.maximum(n, 1))
            lo, hi = np.where(n > 0, lo, p), np.where(n > 0, hi, p)
            ax.plot(
                x,
                p,
                "o-",
                color=col,
                ms=4,
                lw=2 if key.startswith("classification_d") else 1.5,
                label=lab,
            )
            ax.fill_between(x, lo, hi, color=col, alpha=0.15, lw=0)
        ax.plot(
            mis[f"mag_{REF}"],
            mis["missclassification_eff"],
            "s--",
            color="#2ca02c",
            ms=4,
            label="galaxy misclassification (true galaxies called star)",
        )
        for xi, ni in zip(x, eff["n_true_stars"]):
            ax.annotate(
                f"{int(ni)}" if ni else "model",
                (xi, 0.01),
                fontsize=6.5,
                ha="center",
                color="0.45",
            )
        ax.axvline(maglim_ref, color="0.5", ls=":", lw=1)
        ax.text(
            maglim_ref,
            0.5,
            f" anchored 5$\\sigma$ depth {maglim_ref:.2f}",
            rotation=90,
            va="center",
            fontsize=8,
            color="0.4",
        )
        ax.axvspan(MAG_BINS[0] - 1, SATURATION[REF], color="0.9")
        ax.set_xlim(MAG_BINS[0] - 0.5, MAG_BINS[-1])
        ax.set_ylim(-0.02, 1.05)
        ax.set_xlabel("Euclid $I_E$ (PSF)")
        ax.set_ylabel("efficiency")
        ax.grid(alpha=0.25)
        ax.legend(loc="lower left", fontsize=8.5)
        ax.set_title(
            f"Euclid Q1 stellar selection function on the ECDFS truth compilation\n"
            f"{int(is_star.sum()):,} true stars, {int(gal.sum()):,} detected true galaxies "
            f"(Gatto et al. 2026); grey numbers = true stars per bin",
            fontsize=10.5,
        )
        fig.tight_layout()
        fig.savefig(figdir / f"{TAG}_truth_validation.png", dpi=140)
        plt.close(fig)

        # (b) depth anchor: native vs anchored per band + the error curves
        fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
        ax = axes[0]
        for b, col in zip(BANDS, ("k", "#ff7f0e", "#d62728", "#9467bd")):
            v = maps_native[b][np.isfinite(maps_native[b])]
            ax.hist(
                v,
                bins=np.arange(23.5, 28.5, 0.05),
                histtype="step",
                color=col,
                lw=1.3,
                label=f"{b}: native {depth[b]['native_median_bg_extrapolation']:.2f} -> {depth[b]['anchored_median']:.2f}",
            )
            ax.axvline(depth[b]["anchored_median"], color=col, ls="--", lw=1)
        ax.set_xlabel("reported-error S/N = 5 depth per nside-128 pixel (native scale)")
        ax.set_ylabel("pixels")
        ax.legend(fontsize=8.5)
        ax.set_title(
            "Q1 depth maps: structure from the catalogue, scale from the Wide-survey depth",
            fontsize=10,
        )
        ax.grid(alpha=0.25)
        ax = axes[1]
        ax.plot(
            catalog_tab["delta_mag"],
            catalog_tab["log_mag_err"],
            "-",
            color="#1f77b4",
            lw=2,
            label="catalog (reported, detected)",
        )
        ax.plot(
            sample_tab["delta_mag"],
            sample_tab["log_mag_err"],
            "-",
            color="#d62728",
            lw=2,
            label=f"sample = catalog x {f_ref:.2f}",
        )
        ax.plot(
            catalog_nc_tab["delta_mag"],
            catalog_nc_tab["log_mag_err"],
            "--",
            color="#1f77b4",
            lw=1.2,
            label="catalog, no S/N cut",
        )
        ax.plot(
            sample_nc_tab["delta_mag"],
            sample_nc_tab["log_mag_err"],
            "--",
            color="#d62728",
            lw=1.2,
            label="sample, no S/N cut",
        )
        ax.axhline(np.log10(SIG_SN5), color="0.5", ls=":", lw=1)
        ax.axvline(0, color="0.5", ls=":", lw=1)
        ax.set_xlabel("delta_mag = $I_E$ - maglim(pixel)")
        ax.set_ylabel("log$_{10}$ $\\sigma_{I_E}$ [mag]")
        ax.set_xlim(-9, 2)
        ax.legend(fontsize=8.5, loc="upper left")
        ax.grid(alpha=0.25)
        ax.set_title(
            f"VIS photo-error curves ({len(pl):,} point-like sources)", fontsize=10
        )
        fig.tight_layout()
        fig.savefig(figdir / f"{TAG}_depth_anchor.png", dpi=140)
        plt.close(fig)
        print(f"  wrote figures to {figdir}")

    print(f"done in {time.time() - t0:.0f} s")


def build_parser():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--hats", default=HATS_EUCLID, help="Euclid Q1 HATS collection")
    p.add_argument(
        "--truth",
        default=TRUTH_DEFAULT,
        help="ECDFS spectroscopic compilation (catalog.dat.gz)",
    )
    p.add_argument("--out", default=str(OUT_DEFAULT), help="output directory")
    p.add_argument("--figdir", default=None, help="write derivation figures here")
    p.add_argument(
        "--no-anchor",
        action="store_true",
        help="ship the native reported-error depth scale (f = 1) instead of anchoring "
        "to the Wide-survey point-source depths",
    )
    return p


if __name__ == "__main__":
    main(build_parser().parse_args())
