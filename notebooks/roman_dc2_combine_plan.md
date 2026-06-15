# Roman DC2 mock: plan to build a single truth catalog with observations

**Status:** PLAN ONLY — not executed yet.
**Source:** Troxel et al. 2023, *A Joint Roman + Rubin Synthetic Wide-Field Imaging Survey*
([arXiv:2209.06829](https://arxiv.org/abs/2209.06829)).
**Data:** `/astro/store/shire/stream_team/stream_finding/data/roman_mock`
(symlinked as `data/surveys/roman_dc2/`). ~20 deg², **1039 coadd tiles**, tiles named by
center `{ra}_{dec}` in degrees (e.g. `50.93_-38.8`).
**Python:** `/astro/store/shiren/conda-envs/stream_team/envs/streamobs/bin/python`
(astropy 7.2, numpy, pandas, healpy; no fitsio).

---

## 1. What the files actually are (verified by inspection)

| File(s) | Count | Granularity | Key columns | Notes |
|---|---|---|---|---|
| `det/dc2_det_{tile}.fits.gz` | 2078 | per coadd, **one row per detection** | `number`, `alphawin_j2000`,`deltawin_j2000` (RA/Dec **deg**), `flux_auto`,`fluxerr_auto`, `mag_auto`, `mag_auto_{Y106,J129,H158,F184}` (+magerr,flux,fluxerr), `flags`, `class_star`, shape moments | SExtractor on the coadd. ~21k rows/tile. |
| `truth/dc2_index_{tile}.fits.gz` | 1479 (1039 real + 440 empty) | per coadd, **one row per object** | `ind` (truth ID), `ra`,`dec` (**deg**), `mag_{Y106,J129,H158,F184}`, `dered_{...}`, `gal_star` (**1=star, 0=gal**), `sca`,`dither`,`x`,`y`,`stamp`,`start_row` | **THE truth product to use.** 440 files are 53-byte empty placeholders (tiles with no objects). |
| `truth/coadd/dc2_index_{tile}.fits.gz` | 1481 | — | — | **Byte-identical** to `truth/dc2_index_{tile}` (confirmed via `filecmp`). Ignore one of them. |
| `truth/dc2_index_star.fits.gz` | 1 | **per (object × SCA exposure)** | `ind,sca,dither,x,y,ra,dec,mag,stamp,xmin,xmax,ymin,ymax,dudx,dudy,dvdx,dvdy,start_row` | 414 MB, ~6.36M rows. RA/Dec in **RADIANS**. Single-epoch index; **not needed** for coadd matching. |
| `truth/dc2_index.fits.gz` | 1 | **per (object × SCA exposure)** | same 18 cols as star index | **78 GB.** All objects (gal+star). Header `NAXIS2=0` (streaming FITS — astropy reads 0 rows; see §5). **Not needed for this task.** |
| `truth/dc2_coaddlist.fits.gz` | 1 | per coadd tile | `tilename,coadd_i,coadd_j,coadd_ra,coadd_dec,d_ra,d_dec,input_list` | 1838 rows; `input_list` = contributing exposures. Optional metadata. |

### The key realization
**We do NOT need the 78 GB `dc2_index.fits` (nor the star index) for this.**
The per-coadd `truth/dc2_index_{tile}.fits.gz` files already carry, per object:
truth RA/Dec (deg), truth mags in all 4 Roman bands, dereddening, and the star/galaxy
flag. They are tile-for-tile aligned with the `det/dc2_det_{tile}.fits.gz` files. So the
whole job reduces to **per-tile positional matching** — tractable and parallel.

The big index is only needed if you later want single-epoch (per-SCA) info: which exposures
saw an object, its pixel position/postage stamp, etc.

---

## 2. The combine = per-tile positional match (truth ← det)

Following the paper's matching recipe:

For each of the **1039 real coadd tiles** (intersection of det-tile names and non-empty truth-tile names):

1. **Load truth tile** → table of objects (ind, ra, dec, mag_*, dered_*, gal_star, ...).
2. **Load det tile** → detections (alphawin_j2000, deltawin_j2000, mag_auto_*, flux_auto, fluxerr_auto, flags, class_star).
3. **Cut det**: `S/N = flux_auto / fluxerr_auto > 5` (paper's threshold; removes only ~0.05%).
   Optionally also `flags == 0` (stricter; the paper removes ~32% via flags — leave this as a toggle, default OFF).
4. **Match** truth → det by sky position with **1.0 arcsec** radius
   (`astropy.coordinates.SkyCoord.match_to_catalog_sky`, or a KD-tree on a local tangent plane).
   - Paper's rule: among matches within 1″, when ambiguous take the closest **in magnitude**
     among the up-to-3 nearest. v1 can just take the **nearest neighbor within 1″**; add the
     mag-tiebreak as a refinement.
5. **Left-join on truth**: every truth row is kept. Attach the matched det columns
   (prefixed `det_`) when a match exists; else NaN. Add `detected` (bool) and `match_sep_arcsec`.
6. **Tag** each row with `tile` (the `{ra}_{dec}` string).
7. **Concatenate** all tiles → one catalog. Write to **parquet** (fast, typed) — and/or FITS.

Output is then a single **truth-complete** catalog: one row per true object across the
footprint, with observed (detected) quantities where they exist. From it you can trivially
derive detection efficiency / completeness vs. magnitude, color, `gal_star`, etc., which is
exactly what stream-injection / selection-function work needs.

### Output schema (proposed)
```
ind, tile, ra, dec, gal_star,
mag_Y106, mag_J129, mag_H158, mag_F184,           # truth
dered_Y106, dered_J129, dered_H158, dered_F184,   # truth
detected (bool), match_sep_arcsec,
det_number, det_ra, det_dec,
det_mag_auto, det_mag_auto_Y106..F184,
det_flux_auto, det_fluxerr_auto, det_sn,
det_flags, det_class_star
```
(`mag_* == 0.0` in truth means "no flux in that band" → treat as NaN on load.)

---

## 3. Column subset for the big index (`dc2_index.fits`, 78 GB) — if you still want it

You offered to pre-subset it. **For the coadd truth+det catalog you can skip it entirely**
(use the per-tile truth files). If you do subset it (e.g. for per-SCA / single-epoch studies),
the 18 columns are:

| keep? | column | type | meaning |
|---|---|---|---|
| ✅ | `ind` | int64 | object ID — links to input truth catalog (join key across files) |
| ✅ | `ra`,`dec` | float64 | sky position (**radians** in this file) |
| ✅ | `mag` | float64 | true magnitude (in that exposure's band) |
| ✅ | `sca` | int64 | Roman detector (1–18) |
| ✅ | `dither` | int64 | exposure / pointing ID |
| ➖ | `x`,`y` | float64 | pixel position on the SCA |
| ➖ | `stamp`,`start_row` | int64 | postage-stamp pointers into the image-stamp file |
| ❌ | `xmin,xmax,ymin,ymax` | int64 | stamp bounding box |
| ❌ | `dudx,dudy,dvdx,dvdy` | float64 | local WCS Jacobian |

Minimal useful subset: **`ind, sca, dither, ra, dec, mag`** (6 of 18 → ~⅓ the size).
Add `x,y,stamp,start_row` if you need to find/cut postage stamps. Note: a single object
appears in **many rows** here (once per SCA exposure it lands on), so you must group by `ind`.

---

## 4. Decisions to confirm before coding

1. **Stars only, or all objects?** streamobs is stellar streams → likely `gal_star==1`.
   But keeping galaxies lets us model the contaminating background. → *propose: keep all,
   carry `gal_star`, filter downstream.*
2. **flags cut?** S/N>5 only (your call), or also `flags==0` (drops ~32%, removes
   blends/edges). → *propose: S/N>5 default, `flags==0` as an option.*
3. **Match radius / tie-break:** 1.0″ nearest-neighbor (v1) vs. paper's mag-tiebreak FoF. → *propose: 1.0″ NN first, refine if needed.*
4. **Output format/location:** parquet under `data/surveys/roman_dc2/` (+ optional FITS).
5. **Galaxy de-blending:** the paper notes 20–30% of Rubin objects split into multiple
   Roman objects; with a 1″ radius a det can match several truths — we keep the truth-centric
   left join (each truth → its nearest det), so this is handled implicitly. Flag if needed.

---

## 5. Implementation notes / gotchas

- **Read gzipped FITS**: `with gzip.open(path,'rb') as f: Table(astropy.io.fits.open(io.BytesIO(f.read()))[1].data)`.
  (astropy's direct `.fits.gz` open chokes on the big NAXIS2=0 files; the per-tile files open fine.)
- **`NAXIS2=0` streaming files** (star/master index only): astropy returns 0 rows. To read:
  skip the primary + extension headers (2880-byte blocks to `END`), then
  `np.frombuffer(rest, dtype=<144-byte big-endian dtype>)`; n_rows = data_bytes // 144.
- **Units**: per-tile truth + det are **degrees**; star/master index are **radians**.
- **Missing mags**: truth `mag_* == 0.0` and det un-matched → use NaN.
- **Parallelism**: 1039 independent tiles. Natural fan-out — one worker per tile, then concat.
  (Good candidate for a Workflow run, or simple multiprocessing / joblib.)
- **Sanity check tile** `50.93_-38.8`: 34378 truth objects (175 stars, 34203 gals);
  21382 detections, 21372 with S/N>5; truth & det RA/Dec ranges overlap as expected.
