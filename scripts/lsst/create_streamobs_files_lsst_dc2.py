#!/usr/bin/env python
# coding: utf-8
"""LSST DC2 streamobs selection-function products (truth-anchored, r-referenced).

LSST analog of ``scripts/roman/create_streamobs_files_hlwas.py``, built from the
DC2 object + ``truth_merged_summary`` per-tract skims.  Produces, in the same
formats as the Roman products:

  data/surveys/lsst_dc2/lsst_dc2_photoerror_r.csv          # SAMPLE  (truth-based scatter -> noise draw)
  data/surveys/lsst_dc2/lsst_dc2_photoerror_r_catalog.csv  # CATALOG (reported magerr -> S/N cut)
  data/surveys/lsst_dc2/lsst_dc2_stellar_efficiency_cutr.csv
  data/surveys/lsst_dc2/lsst_dc2_maglim_{r,g}_nside1024.fits.gz   # truth-anchored depth maps
  (+ *_raw.csv provenance for the photo-error curves)

Same truth-anchored methodology as Roman (delta_mag = true mag - maglim; S/N>5
baked into the efficiency curves; two-curve photo-error + afterburner; delta_mag>1
detection clamp; corrected ``classification_eff`` spelling).  Key differences:
star/galaxy separation is the LSST ``extendedness`` flag (not a size envelope, so
no classifier module), mags are AB (no Vega->AB), and truth-matching is a direct
single-survey positional match.

Depth maps: the DC2 ``supreme_dc2_dr6d_v3_{b}_maglim_psf_wmean.hs`` HealSparse maps
provide the spatial structure + footprint; we degrade them to nside=1024 and
truth-anchor the absolute scale (shift each so its median = the band's truth-based
S/N=5 depth), so ``delta_mag=0 <-> sigma=0.217`` by construction.

Run with the streamobs env:
  conda activate streamobs
  python scripts/lsst/create_streamobs_files_lsst_dc2.py [--tracts N] [--refresh]
"""

import argparse
import re
from glob import glob
from pathlib import Path

import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from astropy.coordinates import SkyCoord
from astropy.io import fits

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
REPO = Path(__file__).resolve().parents[2]
DATA_DIR = Path("/astro/store/shire/stream_team/stream_finding/data/lsst_dc2")
OUT_DIR = REPO / "data/surveys/lsst_dc2"
FIG_DIR = REPO / "lsst_dc2_scratch/figs"
CORRECTIONS_FILE = (  # generation-time only; lives beside this script, not in config/
    Path(__file__).resolve().parent / "lsst_photoerror_corrections.yaml"
)
CACHE = OUT_DIR / "_cache_matched_stars.parquet"  # matched true-star object rows
CACHE_TS = OUT_DIR / "_cache_truth_stars.parquet"  # truth-star denominator

REF_BAND = "r"  # curves keyed to r (applied to g via delta_mag)
BANDS = ["r", "g"]  # depth maps built for these
SNR_DEPTH = 5  # S/N cut and reference depth
EXT_CUT = 0.5  # extendedness < cut  -> classified as point source (star)
TRUTH_STAR, TRUTH_GAL = 2, 1  # truth_type: 1=galaxy, 2=star, 3=SN
MATCH_RADIUS = 1.0  # arcsec, object -> nearest truth
NSIDE = 1024  # output maglim-map resolution (matches Roman)
MAG_BINS = np.arange(15.0, 29.0 + 1e-6, 0.25)
MAG_MID = 0.5 * (MAG_BINS[1:] + MAG_BINS[:-1])
SIG_SN5 = 2.5 / np.log(10) / SNR_DEPTH  # magerr at S/N=5 = 0.2171
DET_EFF_DELTA_MAX = 1.0  # zero detection_eff for delta_mag > this (faint-tail clamp)
EFF_DELTA_MIN = -11.0  # drop curve rows with delta_mag < this (bright/saturation cut)

OBJ_COLS = (
    ["ra", "dec", "extendedness", "clean", "blendedness"]
    + [f"mag_{b}" for b in BANDS]
    + [f"magerr_{b}" for b in BANDS]
    + [f"psFlux_{REF_BAND}", f"psFluxErr_{REF_BAND}"]
)
TRU_COLS = ["ra", "dec", "truth_type", "id"] + [f"mag_{b}" for b in BANDS]

plt.rcParams.update({"figure.dpi": 110, "font.size": 11})


def load_cols(path, cols):
    """Load selected FITS columns into a DataFrame, converting big-endian FITS
    buffers to native byte order (pandas rejects big-endian arrays)."""
    with fits.open(path) as f:
        d = f[1].data
        out = {}
        for c in cols:
            a = np.asarray(d[c])
            if a.dtype.byteorder == ">":
                a = a.astype(a.dtype.newbyteorder("="))
            out[c] = a
        return pd.DataFrame(out)


def build_matched_catalog(n_tracts=0):
    """Loop tracts: match objects -> nearest truth (<1"), keep matched true-star
    object rows + the all-true-star denominator."""
    obj_frames, truthstar_frames = [], []
    obj_files = sorted(glob(str(DATA_DIR / "dc2_object_run2.2i_dr6_skim_tract_*.fits")))
    if n_tracts:
        obj_files = obj_files[:n_tracts]
    for of in obj_files:
        tract = re.search(r"tract_(\d+)\.fits", of).group(1)
        tf = DATA_DIR / f"dc2_run2.2i_truth_merged_summary_skim_tract_{tract}.fits"
        if not tf.exists():
            print(f"tract {tract}: truth skim not transferred yet -- skipped")
            continue
        obj = load_cols(of, OBJ_COLS)
        tru = load_cols(tf, TRU_COLS)
        c_obj = SkyCoord(ra=obj["ra"].values * u.deg, dec=obj["dec"].values * u.deg)
        c_tru = SkyCoord(ra=tru["ra"].values * u.deg, dec=tru["dec"].values * u.deg)
        idx, d2d, _ = c_obj.match_to_catalog_sky(c_tru)
        obj["match_sep"] = d2d.arcsec
        obj["matched"] = d2d.arcsec < MATCH_RADIUS
        for b in BANDS:
            obj[f"truth_mag_{b}"] = tru[f"mag_{b}"].values[idx]
        obj["truth_type"] = tru["truth_type"].values[idx]
        obj["truth_id"] = tru["id"].values[idx]
        obj_star = obj[obj["matched"] & (obj["truth_type"] == TRUTH_STAR)].copy()
        obj_frames.append(obj_star)
        truthstar_frames.append(
            tru.loc[tru["truth_type"] == TRUTH_STAR, ["id", f"mag_{REF_BAND}"]].copy()
        )
        print(
            f"tract {tract}: {len(obj):,} objs, {obj.matched.mean():.1%} matched; "
            f"{len(obj_star):,} matched true stars"
        )
    cat = pd.concat(obj_frames, ignore_index=True)
    truth_stars = pd.concat(truthstar_frames, ignore_index=True).drop_duplicates("id")
    return cat, truth_stars


def truth_anchor_m5(cat, b, mask):
    """Mag where the TRUTH-BASED scatter of (obs - true) reaches S/N=5.

    DIAGNOSTIC ONLY -- no longer applied to the maps.  The depth scale is the
    pipeline's native reported-error one (see ``native_maglim_map``); this function
    quantifies how much more the real scatter is than the reported errors claim
    (~1.4x for LSST DC2), which the SAMPLE photo-error curve carries instead.
    """
    sub = cat.loc[mask, [f"truth_mag_{b}", f"mag_{b}"]].dropna()
    mt = sub[f"truth_mag_{b}"].values
    dmv = sub[f"mag_{b}"].values - mt
    bins = np.arange(22.0, 28.2, 0.25)
    mid = 0.5 * (bins[1:] + bins[:-1])
    scat = np.full(mid.size, np.nan)
    ib = np.digitize(mt, bins) - 1
    for i in range(mid.size):
        v = dmv[ib == i]
        if v.size >= 100:
            scat[i] = (np.percentile(v, 84) - np.percentile(v, 16)) / 2
    g = np.isfinite(scat)
    return float(np.interp(np.log10(SIG_SN5), np.log10(scat[g]), mid[g]))


def native_maglim_map(band, m5_truth=None):
    """Load the supreme DC2 HealSparse maglim map and degrade to NSIDE, UNSHIFTED.

    The map is the pipeline's own reported-error S/N=5 depth, and that is the scale
    every delta_mag-keyed product uses: ``delta_mag = 0`` is where the *reported*
    magerr reaches ``SIG_SN5`` = 0.2171.  This matches what is available for a real
    survey, where the pipeline's quoted depth is all there is -- so sims and data are
    treated identically.

    The truth-based scatter at that magnitude is LARGER than 0.2171 (reported errors
    are optimistic).  That excess is real and is carried by the separate SAMPLE
    photo-error curve, NOT folded into the depth zero point.

    ``m5_truth`` is only reported, as a diagnostic of the reported-vs-truth gap.
    Returns a dense RING healpix map (hp.UNSEEN off-footprint).
    """
    import healsparse as hsp

    src = OUT_DIR / f"supreme_dc2_dr6d_v3_{band}_maglim_psf_wmean.hs"
    m = hsp.HealSparseMap.read(str(src)).degrade(NSIDE, reduction="mean")
    hpmap = m.generate_healpix_map(nest=False)  # RING, hp.UNSEEN off-footprint
    cov = (hpmap != hp.UNSEEN) & np.isfinite(hpmap)
    raw_med = float(np.median(hpmap[cov]))
    msg = f"  {band}: supreme median={raw_med:.3f} (native reported-error scale, unshifted)"
    if m5_truth is not None:
        msg += (
            f"; truth-based S/N=5 sits at {m5_truth:.3f} "
            f"({m5_truth - raw_med:+.3f} mag) -- DIAGNOSTIC ONLY, not applied"
        )
    print(f"{msg}; {cov.sum():,} pixels @ nside={NSIDE}")
    return hpmap


def main(n_tracts=0, refresh=False):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # ---- 1. matched catalog (cached) ------------------------------------
    if CACHE.exists() and CACHE_TS.exists() and not refresh and not n_tracts:
        print(
            f"loading cached match from {CACHE.name}, {CACHE_TS.name} (use --refresh to rebuild)"
        )
        cat = pd.read_parquet(CACHE)
        truth_stars = pd.read_parquet(CACHE_TS)
    else:
        cat, truth_stars = build_matched_catalog(n_tracts)
        if not n_tracts:
            cat.to_parquet(CACHE)
            truth_stars.to_parquet(CACHE_TS)
            print(f"cached match -> {CACHE.name}, {CACHE_TS.name}")
    print(
        f"\n{len(cat):,} matched true-star object rows; {len(truth_stars):,} unique true stars"
    )

    sn = cat[f"psFlux_{REF_BAND}"] / cat[f"psFluxErr_{REF_BAND}"]
    clean = cat["clean"].astype(bool)
    det_ok = (sn > SNR_DEPTH) & clean
    is_ptsrc = cat["extendedness"] < EXT_CUT
    star_ptsrc = is_ptsrc & clean  # matched true stars are all rows in cat

    # ---- 2. depth maps (r, g): NATIVE pipeline scale, not truth-anchored -
    print("\nbuilding maglim maps (native reported-error S/N=5 scale):")
    m5_truth = {b: truth_anchor_m5(cat, b, star_ptsrc) for b in BANDS}  # diagnostic
    maglim_maps = {b: native_maglim_map(b, m5_truth[b]) for b in BANDS}
    for b in BANDS:
        mlm = maglim_maps[b]
        cov = mlm != hp.UNSEEN
        fname = OUT_DIR / f"lsst_dc2_maglim_{b}_nside{NSIDE}.fits.gz"
        hp.write_map(fname, mlm, overwrite=True, dtype=np.float32)
        print(f"  wrote {fname.name} (median {np.median(mlm[cov]):.2f})")

    mlm_r = maglim_maps[REF_BAND]
    covered = mlm_r != hp.UNSEEN
    MAGLIM_REF = float(np.median(mlm_r[covered]))  # curves keyed to the r map median
    print(f"reference maglim = r map median = {MAGLIM_REF:.3f}")

    # ---- 3. detection & classification efficiency (true stars) ----------
    star_dets = (
        cat.loc[det_ok, ["truth_id", "match_sep", "extendedness"]]
        .sort_values("match_sep")
        .drop_duplicates("truth_id")
    )
    detected_ids = set(star_dets["truth_id"].dropna().astype("i8"))
    classified_ids = set(
        star_dets.loc[star_dets["extendedness"] < EXT_CUT, "truth_id"]
        .dropna()
        .astype("i8")
    )
    ts = truth_stars.copy()
    ts["is_det"] = ts["id"].astype("i8").isin(detected_ids)
    ts["is_cls"] = ts["id"].astype("i8").isin(classified_ids)
    mag = ts[f"mag_{REF_BAND}"].values
    ok = np.isfinite(mag)
    n_all = np.histogram(mag[ok], MAG_BINS)[0]
    n_det = np.histogram(mag[ok & ts["is_det"].values], MAG_BINS)[0]
    n_cls = np.histogram(mag[ok & ts["is_cls"].values], MAG_BINS)[0]
    with np.errstate(invalid="ignore"):
        eff_det = n_det / n_all
        eff_cls = np.where(n_det > 0, n_cls / np.maximum(n_det, 1), np.nan)
        eff_both = n_cls / n_all
    print(
        f"efficiency sample: {ok.sum():,} true stars "
        f"({len(detected_ids):,} detected, {len(classified_ids):,} classified)"
    )

    # ---- 4. two-curve photo-error (delta keyed to the per-pixel r map) ---
    # TWO SELECTIONS, two pairs of curves.
    #
    # det_ok (reported S/N > 5) -> the REFERENCE band.  The injector draws the
    #   detection flag and the noise independently, so the reference-band curve must
    #   describe the population it is applied to: the detected one.  See
    #   docs/source/selection_function_methodology.md, "Why the curves are measured on
    #   the *detected* population (validated 2026-07)".
    # no cut -> every OTHER band.  Those bands are FORCED photometry: an object
    #   detected in r is measured in g whatever its g flux, so the g population is not
    #   conditioned on g-band detection and the det_ok curve understates its errors by
    #   up to 2.6x at delta_mag = +1.  Validated against the directly-measured g curve
    #   on the r-detected sample: the r-derived no-cut curve reproduces it to <=5%
    #   across the full range (notebooks/detok_photoerr_comparison.ipynb).
    def _photoerr_tables(mask, label):
        pe = cat.loc[
            mask,
            ["ra", "dec", f"truth_mag_{REF_BAND}", f"mag_{REF_BAND}", f"magerr_{REF_BAND}"],
        ].dropna()
        pix = hp.ang2pix(NSIDE, pe["ra"].values, pe["dec"].values, lonlat=True)
        ml_local = mlm_r[pix]
        good = ml_local != hp.UNSEEN
        delta = pe[f"truth_mag_{REF_BAND}"].values[good] - ml_local[good]
        dm_obs = (
            pe[f"mag_{REF_BAND}"].values - pe[f"truth_mag_{REF_BAND}"].values
        )[good]
        logerr_reported = np.log10(pe[f"magerr_{REF_BAND}"].values[good])

        dbins = np.arange(np.floor(delta.min() * 10) / 10, 1.5 + 1e-6, 0.12)
        dmid = 0.5 * (dbins[1:] + dbins[:-1])
        log_scatter = np.full(dmid.size, np.nan)
        med_logerr_rep = np.full(dmid.size, np.nan)
        ib = np.digitize(delta, dbins) - 1
        for i in range(dmid.size):
            v = dm_obs[ib == i]
            if v.size >= 20:
                log_scatter[i] = np.log10(
                    (np.percentile(v, 84) - np.percentile(v, 16)) / 2
                )
                med_logerr_rep[i] = np.median(logerr_reported[ib == i])
        keep = np.isfinite(log_scatter)
        factor = 10 ** (log_scatter[keep] - med_logerr_rep[keep])
        near = (dmid[keep] > -3) & (dmid[keep] < 0.5)
        print(
            f"  [{label:<34}] n={len(delta):>9,}  rows={int(keep.sum()):>4}  "
            f"inflation={np.nanmedian(factor[near]):.2f}"
        )
        return (
            pd.DataFrame({"delta_mag": dmid[keep], "log_mag_err": log_scatter[keep]}),
            pd.DataFrame({"delta_mag": dmid[keep], "log_mag_err": med_logerr_rep[keep]}),
        )

    print("\n" + "=" * 74)
    print("PHOTO-ERROR CURVES (inflation = truth scatter / reported magerr)")
    photoerr_tab, catalog_tab = _photoerr_tables(
        star_ptsrc & det_ok, "det_ok -> reference band"
    )
    photoerr_nc, catalog_nc = _photoerr_tables(
        star_ptsrc, "no cut -> forced-photometry bands"
    )
    print("=" * 74 + "\n")

    def _write_curve(tab, stem, rule):
        np.savetxt(
            OUT_DIR / f"{stem}_raw.csv", tab.values, delimiter=",",
            header="delta_mag,log_mag_err", fmt="%.6f",
        )
        clean = _apply_photoerr_corrections(tab, rule, CORRECTIONS_FILE)
        np.savetxt(
            OUT_DIR / f"{stem}.csv", clean.values, delimiter=",",
            header="delta_mag,log_mag_err", fmt="%.6f",
        )
        return clean

    photoerr_clean = _write_curve(photoerr_tab, "lsst_dc2_photoerror_r", "r_sample")
    catalog_clean = _write_curve(
        catalog_tab, "lsst_dc2_photoerror_r_catalog", "r_catalog"
    )
    _write_curve(photoerr_nc, "lsst_dc2_photoerror_r_nocut", "r_sample")
    _write_curve(catalog_nc, "lsst_dc2_photoerror_r_catalog_nocut", "r_catalog")
    print(
        f"wrote 4 photo-error curves + *_raw: sample/catalog (det_ok, reference band) "
        f"and sample/catalog _nocut (forced-photometry bands); "
        f"{len(photoerr_clean)} / {len(photoerr_nc)} rows"
    )

    # ---- 5. stellar efficiency table ------------------------------------
    eff_tab = pd.DataFrame(
        {
            "mag_r": MAG_MID,
            "delta_mag": MAG_MID - MAGLIM_REF,
            "detection_eff": eff_det,
            "classification_eff": eff_cls,
            "classification_detection_eff": eff_both,
        }
    )
    eff_tab = eff_tab[n_all >= 20].fillna(0.0)
    # bright cut: drop rows entirely so the injector's bright-edge hold
    # (flat at the table's first remaining row, below the physical
    # saturation floor) governs brighter magnitudes instead of noisy bins
    _bright = eff_tab["delta_mag"] < EFF_DELTA_MIN
    eff_tab = eff_tab[~_bright]
    print(
        f"  dropped {int(_bright.sum())} bins with delta_mag < {EFF_DELTA_MIN} (bright cut)"
    )
    _faint = eff_tab["delta_mag"] > DET_EFF_DELTA_MAX
    eff_tab.loc[_faint, "detection_eff"] = 0.0
    eff_tab.loc[_faint, "classification_detection_eff"] = 0.0
    print(
        f"  zeroed detection_eff for {int(_faint.sum())} bins with delta_mag > {DET_EFF_DELTA_MAX}"
    )
    np.savetxt(
        OUT_DIR / "lsst_dc2_stellar_efficiency_cutr.csv",
        eff_tab.values,
        delimiter=",",
        header="mag_r,delta_mag,detection_eff,classification_eff,classification_detection_eff",
        fmt="%.6f",
    )
    print(f"wrote lsst_dc2_stellar_efficiency_cutr.csv ({len(eff_tab)} rows)")

    # ---- 6. sanity overlay vs existing external LSST tables -------------
    ext_pe = np.genfromtxt(
        REPO / "data/others/lsst_photoerror_r.csv", delimiter=",", names=True
    )
    ext_eff = np.genfromtxt(
        REPO / "data/others/lsst_stellar_efficiency_cutr.csv", delimiter=",", names=True
    )
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].plot(
        ext_pe["delta_mag"],
        ext_pe["log_mag_err"],
        "-",
        color="0.6",
        label="existing (Tsiane+25)",
    )
    axes[0].plot(
        photoerr_clean.delta_mag,
        photoerr_clean.log_mag_err,
        "C3o-",
        ms=3,
        label="DC2 sample (truth)",
    )
    axes[0].plot(
        catalog_clean.delta_mag,
        catalog_clean.log_mag_err,
        "C0--",
        label="DC2 catalog (reported)",
    )
    axes[0].set(
        xlabel=r"$\Delta$mag", ylabel=r"$\log_{10}\sigma$ [mag]", title="photo-error"
    )
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)
    axes[1].plot(
        ext_eff["delta_mag"],
        ext_eff["classification_detection_eff"],
        "-",
        color="0.6",
        label="existing",
    )
    axes[1].plot(
        eff_tab.delta_mag,
        eff_tab.classification_detection_eff,
        "C3o-",
        ms=3,
        label="DC2 (truth-anchored)",
    )
    axes[1].set(
        xlabel=r"$\Delta$mag", ylabel="classification_detection_eff", title="efficiency"
    )
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    fig.suptitle(
        f"LSST DC2 r-band products (maglim_r={MAGLIM_REF:.2f}, S/N>{SNR_DEPTH})"
    )
    fig.tight_layout()
    fig.savefig(FIG_DIR / "lsst_dc2_products_r.png", dpi=130, bbox_inches="tight")
    print(f"wrote {(FIG_DIR / 'lsst_dc2_products_r.png').relative_to(REPO)}")

    print(
        '\nTODO: galaxy misclassification curve (cosmoDC2 size_true<0.3") + lsst_dc2.yaml wiring'
    )


def _apply_photoerr_corrections(tab, curve_id, corrections_path):
    """Apply afterburner corrections from a YAML file (mirrors the Roman generator)."""
    if not corrections_path.exists():
        print(
            f"  [afterburner] corrections file not found: {corrections_path} — skipping"
        )
        return tab.copy()
    with open(corrections_path) as fh:
        corr = yaml.safe_load(fh) or {}
    rules = corr.get(curve_id, {})
    if not rules:
        return tab.copy()
    out = tab.copy()
    for rule_name, params in rules.items():
        if rule_name == "clamp_faint":
            d_min = float(params["delta_mag_min"])
            v = float(params["value"])
            mask = out["delta_mag"] >= d_min
            out.loc[mask, "log_mag_err"] = np.maximum(out.loc[mask, "log_mag_err"], v)
            print(
                f"  [afterburner] {curve_id}: clamp_faint delta_mag>={d_min:.3f} "
                f"-> log_mag_err=max(raw,{v:.4f}) ({int(mask.sum())} bins)"
            )
        elif rule_name == "cut_bright":
            d_min = float(params["delta_mag_min"])
            n_drop = int((out["delta_mag"] < d_min).sum())
            out = out[out["delta_mag"] >= d_min].copy()
            print(
                f"  [afterburner] {curve_id}: cut_bright dropped {n_drop} bins "
                f"with delta_mag < {d_min:.3f}"
            )
        else:
            print(
                f"  [afterburner] WARNING: unknown rule '{rule_name}' for '{curve_id}' — skipped"
            )
    return out


def reapply_corrections():
    """Rewrite the cleaned photo-error CSVs from their *_raw.csv provenance.

    Re-deriving the curves needs the DC2 object and truth skims (~120 GB), which
    is disproportionate when only the afterburner YAML has changed. The raw
    files are the untouched measurement, so re-cleaning them reproduces exactly
    what a full run would write. Mirrors --reapply-corrections in
    scripts/des/balrog_selection_function.py.

    Note the _nocut curves share the `r_sample` / `r_catalog` rule ids with the
    reference-band pair, so a rule change applies to all four and they stay on a
    common delta_mag grid.
    """
    n = 0
    for stem, rule in [
        ("lsst_dc2_photoerror_r", "r_sample"),
        ("lsst_dc2_photoerror_r_catalog", "r_catalog"),
        ("lsst_dc2_photoerror_r_nocut", "r_sample"),
        ("lsst_dc2_photoerror_r_catalog_nocut", "r_catalog"),
    ]:
        raw = OUT_DIR / f"{stem}_raw.csv"
        if not raw.exists():
            print(f"  [afterburner] missing {raw.name} — skipped")
            continue
        tab = pd.DataFrame(
            np.genfromtxt(raw, delimiter=",", comments="#"),
            columns=["delta_mag", "log_mag_err"],
        )
        clean = _apply_photoerr_corrections(tab, rule, CORRECTIONS_FILE)
        np.savetxt(
            OUT_DIR / f"{stem}.csv", clean.values, delimiter=",",
            header="delta_mag,log_mag_err", fmt="%.6f",
        )
        print(f"  wrote {stem}.csv  ({len(tab)} raw -> {len(clean)} rows)")
        n += 1
    print(f"re-applied {CORRECTIONS_FILE.name} to {n} curves in {OUT_DIR}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--tracts",
        type=int,
        default=0,
        help="limit to the first N tracts (0 = all 79; also bypasses the cache)",
    )
    ap.add_argument("--refresh", action="store_true", help="rebuild the match cache")
    ap.add_argument(
        "--reapply-corrections",
        action="store_true",
        help="only re-clean the photo-error curves from their *_raw.csv files, "
        "without re-deriving anything (no DC2 skims needed)",
    )
    args = ap.parse_args()
    if args.reapply_corrections:
        reapply_corrections()
    else:
        main(n_tracts=args.tracts, refresh=args.refresh)
