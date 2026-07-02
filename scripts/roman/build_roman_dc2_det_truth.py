#!/usr/bin/env python
"""Build a DETECTION-centric (det -> truth) catalog for the Roman-Rubin DC2 mock
(Troxel et al. 2023, arXiv:2209.06829) for catalog-level injection / selection work.

Per coadd tile: load the SExtractor detection catalog + the per-object truth index,
cut detections to S/N>5 (KEEP ALL FLAGS), and match each DETECTION to a truth source
following the paper's recipe: among truth objects within 1", take the up to 3 nearest on
the sky, then assign the one closest in MAGNITUDE (mag-tiebreak). EVERY detection is kept
-- unmatched ones (no truth within 1") get matched=False and NaN truth columns, so the
spurious / false-detection rate is preserved. One row per detection.
Concatenate all tiles into one parquet. REPLACES the old truth->det product.

Mag-tiebreak magnitude: the detection image is the median of the 4 single-band coadds, so
the detection `mag_auto` is compared to a truth broadband mag = mean flux over the 4 truth
bands (truth mag==0.0 means "no flux in band" -> ignored).

Decisions (per user): det->truth direction, mag-tiebreak, KEEP unmatched (flagged),
KEEP all flags, replace existing catalog.
Resource bounds: <=24 worker processes, each holds one small tile.

Run with the pferguso_hats env python (needs pyarrow):
  /astro/users/pferguso/.conda/envs/pferguso_hats/bin/python
"""
import os
# keep every worker single-threaded so N_proc ~= N_cores (stay under the core budget)
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
          "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[v] = "1"

import gzip, io, glob, argparse, time, traceback
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord, search_around_sky
import astropy.units as u
import pyarrow.parquet as pq
from multiprocessing import Pool

DATA = "/astro/store/shire/stream_team/stream_finding/data/roman_mock"
OUTDIR = "/astro/store/shire/pferguso/software/streamobs/data/surveys/roman_dc2"
TILEDIR = os.path.join(OUTDIR, "_tiles_det")
FINAL = os.path.join(OUTDIR, "roman_dc2_det_truth.parquet")

BANDS = ["Y106", "J129", "H158", "F184"]
MATCH_RADIUS_ARCSEC = 1.0   # paper's truth-match radius
SN_MIN = 5.0                # paper's detection threshold
N_NEAREST = 3               # paper: disambiguate among up to 3 nearest within radius

# SExtractor detection columns (53), in file order -- kept verbatim
DET_COLS = ['number', 'flux_auto', 'fluxerr_auto', 'mag_auto', 'magerr_auto', 'kron_radius',
            'background', 'isoareaf_image', 'xwin_image', 'ywin_image', 'alphawin_j2000',
            'deltawin_j2000', 'x2win_image', 'y2win_image', 'xywin_image', 'x2win_world',
            'y2win_world', 'xywin_world', 'awin_world', 'bwin_world', 'thetawin_world',
            'mu_threshold', 'mu_max', 'flags', 'class_star',
            'mag_auto_Y106', 'mag_auto_J129', 'mag_auto_H158', 'mag_auto_F184',
            'magerr_auto_Y106', 'magerr_auto_J129', 'magerr_auto_H158', 'magerr_auto_F184',
            'flux_auto_Y106', 'flux_auto_J129', 'flux_auto_H158', 'flux_auto_F184',
            'fluxerr_auto_Y106', 'fluxerr_auto_J129', 'fluxerr_auto_H158', 'fluxerr_auto_F184',
            'x2win_world_Y106', 'y2win_world_Y106', 'xywin_world_Y106',
            'x2win_world_J129', 'y2win_world_J129', 'xywin_world_J129',
            'x2win_world_H158', 'y2win_world_H158', 'xywin_world_H158',
            'x2win_world_F184', 'y2win_world_F184', 'xywin_world_F184']

TRUTH_MAG = [f"mag_{b}" for b in BANDS]
TRUTH_DERED = [f"dered_{b}" for b in BANDS]
TRUTH_OUT = (["truth_ind", "truth_ra", "truth_dec", "truth_gal_star"] +
             [f"truth_{m}" for m in TRUTH_MAG] + [f"truth_{m}" for m in TRUTH_DERED] +
             ["truth_bb_mag"])

# final column order: tile, all detection cols, the matched-truth payload
COLS = ["tile"] + DET_COLS + ["det_sn", "matched", "match_sep_arcsec"] + TRUTH_OUT


def load_gz_fits(path):
    with gzip.open(path, "rb") as f:
        raw = f.read()
    with fits.open(io.BytesIO(raw)) as h:
        return Table(h[1].data)


def truth_broadband_mag(truth):
    """Mean-flux magnitude across the 4 truth bands (mag==0.0 = missing -> ignored).
    The detection image is the median of the 4 single-band coadds, so mean band flux is the
    natural truth proxy for the detection's mag_auto used in the mag-tiebreak."""
    m = np.vstack([np.asarray(truth[c], float) for c in TRUTH_MAG]).T   # (n, 4)
    m[m == 0.0] = np.nan
    with np.errstate(over="ignore"):
        flux = 10.0 ** (-0.4 * m)
    with np.errstate(invalid="ignore", divide="ignore"):
        return -2.5 * np.log10(np.nanmean(flux, axis=1))   # NaN if all bands missing


def process_tile(tile):
    """Worker entry point: never raises, so the pool keeps going on a bad tile."""
    try:
        return _process_tile_impl(tile)
    except Exception:
        return ("err", tile, traceback.format_exc())


def _process_tile_impl(tile):
    """Return ('ok'|'skip', tile, parquet_path)."""
    out_path = os.path.join(TILEDIR, f"{tile}.parquet")
    if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
        return ("skip", tile, out_path)

    det = load_gz_fits(os.path.join(DATA, "det", f"dc2_det_{tile}.fits.gz"))
    truth = load_gz_fits(os.path.join(DATA, "truth", f"dc2_index_{tile}.fits.gz"))

    # S/N cut on detections (no flags cut, per request)
    with np.errstate(divide="ignore", invalid="ignore"):
        det_sn = np.asarray(det["flux_auto"], float) / np.asarray(det["fluxerr_auto"], float)
    keep = np.isfinite(det_sn) & (det_sn > SN_MIN)
    det = det[keep]
    det_sn = det_sn[keep]
    ndet = len(det)

    # base output: every surviving detection, all SExtractor columns, truth cols = NaN
    df = det[DET_COLS].to_pandas()
    df.insert(0, "tile", tile)
    df["det_sn"] = det_sn
    df["matched"] = False
    df["match_sep_arcsec"] = np.nan
    for c in TRUTH_OUT:
        df[c] = np.nan

    if ndet > 0 and len(truth) > 0:
        det_mag = np.asarray(det["mag_auto"], float)
        tbb = truth_broadband_mag(truth)
        cd = SkyCoord(np.asarray(det["alphawin_j2000"], float) * u.deg,
                      np.asarray(det["deltawin_j2000"], float) * u.deg)
        ct = SkyCoord(np.asarray(truth["ra"], float) * u.deg,
                      np.asarray(truth["dec"], float) * u.deg)
        # all (det, truth) pairs within the match radius
        i_d, i_t, sep2d, _ = search_around_sky(cd, ct, MATCH_RADIUS_ARCSEC * u.arcsec)
        if len(i_d) > 0:
            pairs = pd.DataFrame({"d": i_d, "t": i_t, "sep": sep2d.arcsec,
                                  "dmag": np.abs(tbb[i_t] - det_mag[i_d])})
            # keep up to N_NEAREST nearest truths per detection ...
            pairs = pairs.sort_values(["d", "sep"], kind="mergesort")
            pairs["rank"] = pairs.groupby("d", sort=False).cumcount()
            pairs = pairs[pairs["rank"] < N_NEAREST]
            # ... then assign the one closest in magnitude (NaN dmag sorts last)
            pairs = pairs.sort_values(["d", "dmag"], kind="mergesort", na_position="last")
            best = pairs.drop_duplicates("d", keep="first")
            di = best["d"].to_numpy()
            ti = best["t"].to_numpy()
            df.loc[di, "matched"] = True
            df.loc[di, "match_sep_arcsec"] = best["sep"].to_numpy()
            df.loc[di, "truth_ind"] = np.asarray(truth["ind"], float)[ti]
            df.loc[di, "truth_ra"] = np.asarray(truth["ra"], float)[ti]
            df.loc[di, "truth_dec"] = np.asarray(truth["dec"], float)[ti]
            df.loc[di, "truth_gal_star"] = np.asarray(truth["gal_star"], float)[ti]
            for m in TRUTH_MAG + TRUTH_DERED:
                v = np.asarray(truth[m], float)
                if m in TRUTH_MAG:                 # truth mag 0.0 -> NaN (no flux in band)
                    v = np.where(v == 0.0, np.nan, v)
                df.loc[di, f"truth_{m}"] = v[ti]
            df.loc[di, "truth_bb_mag"] = tbb[ti]

    df = df[COLS]
    df.to_parquet(out_path, engine="pyarrow", index=False, compression="snappy")
    return ("ok", tile, out_path)


def discover_tiles():
    tiles = []
    for p in glob.glob(os.path.join(DATA, "truth", "dc2_index_[0-9]*.fits.gz")):
        if os.path.getsize(p) <= 1000:               # 53-byte empty placeholders
            continue
        t = os.path.basename(p)[len("dc2_index_"):-len(".fits.gz")]
        if os.path.exists(os.path.join(DATA, "det", f"dc2_det_{t}.fits.gz")):
            tiles.append(t)
    return sorted(tiles)


def combine(tiles):
    """Stream per-tile parquet files into one final parquet (low memory)."""
    writer = None
    nrows = 0
    for t in tiles:
        p = os.path.join(TILEDIR, f"{t}.parquet")
        tbl = pq.read_table(p)
        if writer is None:
            writer = pq.ParquetWriter(FINAL, tbl.schema, compression="snappy")
        else:
            tbl = tbl.cast(writer.schema)
        writer.write_table(tbl)
        nrows += tbl.num_rows
    if writer is not None:
        writer.close()
    return nrows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nproc", type=int, default=24)
    ap.add_argument("--limit", type=int, default=0, help="process only first N tiles (smoke test)")
    ap.add_argument("--combine-only", action="store_true")
    args = ap.parse_args()
    args.nproc = max(1, min(args.nproc, 24))             # hard cap at 24 cores

    os.makedirs(TILEDIR, exist_ok=True)
    tiles = discover_tiles()
    if args.limit:
        tiles = tiles[:args.limit]
    print(f"[{time.strftime('%H:%M:%S')}] tiles to process: {len(tiles)}  nproc={args.nproc}", flush=True)

    if not args.combine_only:
        t0 = time.time()
        done = errs = 0
        with Pool(args.nproc) as pool:
            for status, tile, info in pool.imap_unordered(process_tile, tiles, chunksize=2):
                done += 1
                if status == "err":
                    errs += 1
                    print(f"  ERROR {tile}: {info}", flush=True)
                if done % 100 == 0 or done == len(tiles):
                    print(f"[{time.strftime('%H:%M:%S')}] {done}/{len(tiles)} "
                          f"({time.time()-t0:.0f}s, errs={errs})", flush=True)
        print(f"[{time.strftime('%H:%M:%S')}] per-tile done: {done} ok-ish, {errs} errors", flush=True)

    print(f"[{time.strftime('%H:%M:%S')}] combining into {FINAL} ...", flush=True)
    nrows = combine(tiles)
    size_gb = os.path.getsize(FINAL) / 1e9
    print(f"[{time.strftime('%H:%M:%S')}] DONE: {FINAL}  rows={nrows:,}  size={size_gb:.2f} GB", flush=True)


if __name__ == "__main__":
    main()
