# Selection functions from Balrog (DECam surveys)

This technote covers how the streamobs selection-function products are derived
for the DECam surveys — **DES Y6 Gold** and **DELVE DR3 Gold** — from
Balrog synthetic-source-injection (SSI) catalogues, and how to re-run the
derivation on a machine that holds the injections.

It is the DECam counterpart to {doc}`selection_function_methodology`, which
defines the products, the `delta_mag` convention and the truth-anchoring recipe
survey-agnostically. **Read that page first**; this one only covers what is
specific to Balrog, and the one problem that is unique to DES/DELVE: the star
classifier cannot be evaluated on injected sources.

Everything here is driven by a single script,
`scripts/des/balrog_selection_function.py`, which carries a schema adapter per
survey (`--survey delve`, `--survey des_y6`). It streams the HDF5 in chunks,
imports nothing from streamobs, and depends only on numpy / h5py / pandas /
healpy (plus healsparse for `.hsp` depth maps and xgboost for the DES
surrogate) — so it runs unmodified on the machine that holds the catalogues and
emits only ~10 KB of CSVs plus optional depth maps.

## Why Balrog needs special handling

A Balrog catalogue is an *injection* catalogue: known sources are painted into
real survey images, the images are re-reduced, and the recovered measurements
are matched back to what went in. That gives exactly what a selection function
needs — a truth table and a measured table for the same objects, in the real
survey's own noise, masking and crowding.

Two things make it awkward in practice:

1. **The star/galaxy classifier the survey actually publishes cannot be
   recomputed on the injections.** See *The EXT_XGB problem* below.
2. **The truth star/galaxy label is not in the injection files.** For DES it has
   to be joined in from the parent deep-field catalogue.

Both are solved, and both cost accuracy in ways that are quantified rather than
hidden.

## The EXT_XGB problem

DES Y6 Gold's best star/galaxy classifier is `EXT_XGB`, an XGBoost model over
six features (Bechtol et al. 2025, [arXiv:2501.05739](https://arxiv.org/abs/2501.05739),
Table A.2): `CONC`, `BDF_T`, `BDF_T_ERR`, `log10(BDF_S2N)`,
`WAVG_SPREAD_MODEL_I`, `WAVG_SPREADERR_MODEL_I`. Its integer classes come from
thresholding the continuous output (their Eq. A5):

```
EXT_XGB = (XGB_PRED < 0.865) + (XGB_PRED < 0.110)
        + (XGB_PRED < 0.045) + (XGB_PRED < 0.015)
```

so `EXT_XGB = 0` is the purest stellar class and `0 <= EXT_XGB <= 1` is the
"complete" stellar selection. `EXT_XGB = -9` means no classification.

**It cannot be evaluated on injections.** Three of the six features are never
measured for injected sources. The Gold paper says so directly (App. A.2):

> the XGBoost classification has not been implemented on
> simulation–injection–recovery tests due to the fact that the CONC parameter
> was not measured for the simulated object samples

and the Balrog paper (Anbajagane et al. 2025,
[arXiv:2501.05683](https://arxiv.org/abs/2501.05683), Sec. 4.4) reaffirms it and
falls back to `EXT_MASH` — which that same paper shows carries **>10% galaxy
contamination** when used to select stars. Verified directly against both
released Balrog files: no `EXT_XGB`, no `XGB_PRED`, no `SPREAD_MODEL`, no
`CONC`.

Neither paper offers a workaround. streamobs needs one, because the quantity the
injector consumes is

```
classification_eff(delta_mag) = P(EXT_XGB <= 1 | true star, detected)
```

for the cut that will actually be applied to the real catalogue.

### What we do instead: surrogate + deconvolution

The key observation is that Balrog has the truth label but not `EXT_XGB`, while
the real Y6 Gold catalogue has `EXT_XGB` but no truth label — and the two
overlap in *feature* space.

**Step 1 — train a surrogate on the real catalogue.**
`scripts/des/build_des_xgb_surrogate.py` trains a classifier `S` on real
`des_y6_gold` rows against the label `EXT_XGB <= 1`, using only features that
can be constructed *identically* on both sides. All are fitvd/SOF quantities,
and the Balrog files are the `_sof` products, so there is no cross-pipeline
mismatch:

| feature | real Y6 Gold | Balrog matched file |
|---|---|---|
| `BDF_T`, `BDF_T_ERR`, `BDF_T_RATIO` | as named | `meas_bdf_T`, `_err`, `_ratio` |
| `log_bdf_s2n` | `log10(BDF_S2N)` | `log10(meas_bdf_s2n)` |
| `conc_gap` | `GAP_MAG_I - BDF_MAG_I` | `meas_gap_mag[:,2] - meas_bdf_mag[:,2]` |
| `conc_ap8` | `PSF_MAG_APER_8_I - BDF_MAG_I` | `meas_psf_mag_aper8[:,2] - meas_bdf_mag[:,2]` |
| `conc_ap8_r` | `PSF_MAG_APER_8_R - BDF_MAG_R` | `meas_psf_mag_aper8[:,1] - meas_bdf_mag[:,1]` |
| `BDF_MAG_I`, `BDF_MAG_G` | as named | `meas_bdf_mag[:,2]`, `[:,0]` |

Balrog's `(N,4)` arrays are ordered g,r,i,z.

The **GAp** (Gaussian-aperture) concentration matters and is not
interchangeable with the aperture one: the paper's `CONC` is a Gaussian-weighted
flux *dilation* ratio (their Sec. III.4 — flux through a PSF-sized Gaussian
versus one 5% wider, calibrated on the PSF model), not an aperture-magnitude
difference, and the paper explicitly flags GAp fluxes as the Balrog-appropriate
tool. Measured feature importance bears this out: `conc_gap` is the single most
important feature.

**Step 2 — deconvolve.** `S` is not `EXT_XGB`, so applying it to Balrog measures
`eff_S = P(S=1 | star)`, not `eff_X`. The confusion between them is measurable
on the real catalogue, per magnitude bin:

```
a(m) = P(S=1 | EXT_XGB <= 1)        b(m) = P(S=1 | EXT_XGB > 1)
eff_S = a*eff_X + b*(1 - eff_X)   ->   eff_X = (eff_S - b) / (a - b)
```

The reducer applies this inversion with `--confusion`, only where the classes
are well separated (`a - b > 0.2`; as `a -> b` the inversion amplifies noise
without bound). This removes the *first-order bias* rather than quoting the
surrogate's disagreement rate as an irreducible systematic.

> **Without the deconvolution `classification_eff` describes the surrogate, not
> `EXT_XGB`.** The reducer warns if `--confusion` is omitted.

### How good is it, and where does it break

Two checks establish this, both run before anything downstream was built.

**`CONC` reconstruction** (`--precheck`). The real catalogue carries `CONC`
itself, so the headroom is directly measurable — regress `CONC` on the
Balrog-available features and see what is recoverable:

| i-band | R² | residual scatter | `CONC` intrinsic spread |
|---|---|---|---|
| 18–19 | 0.851 | 0.0003 | 0.0060 |
| 20–21 | 0.912 | 0.0006 | 0.0063 |
| 22–23 | 0.878 | 0.0011 | 0.0045 |
| 23–24 | 0.689 | 0.0020 | 0.0043 |
| 24–25 | 0.393 | 0.0034 | 0.0049 |

`CONC` is essentially recoverable brightward of i ≈ 23 (residual ~20× below its
own spread) and degrades faintward, which is where the deconvolution does its
work.

**Domain shift.** A surrogate trained on one catalogue and applied to another is
only valid if the feature distributions agree. Comparing real versus Balrog at
matched magnitude flagged one feature hard: **`PSF_T` is offset by ~1.4σ at
every magnitude** (real median 1.035–1.045, Balrog 0.882–0.886). A constant
offset independent of magnitude is a pipeline difference, not a population one
(cf. the PSFEx/Piff mismatch noted in the Balrog paper, Sec. 3.4). `PSF_T` and
the derived `psf_bdf_T` are therefore **excluded** from the feature set, along
with `conc_gap_r`, which is redundant — fitvd fits one shape jointly across
bands with per-band fluxes, so `GAP - BDF` is band-independent by construction.

Dropping all three changes AUC, agreement and recall by less than 1e-4, so they
carried no unique information. Final performance on 16.8M real rows:

- **AUC 0.9946**, selection agreement **0.9802**, recall 0.8813, precision 0.9491
- `a(m)` falls from 1.000 at g ≈ 17 to 0.642 at g ≈ 27, while `b(m)` stays at
  0.002–0.011 — the surrogate is very pure but loses recall faintward, which is
  exactly what the deconvolution corrects. Left uncorrected,
  `classification_eff` would be **~17% low by g ≈ 25**.

**The residual systematic** is that `a` and `b` are measured on the real
catalogue's mixed star+galaxy population, while what is needed is the same
quantity conditioned on true stars. That dependence is second-order, and it is
checked externally rather than assumed away (see *Validation* below).

## The truth star/galaxy label

**DES Y6 Balrog has no truth label column.** It does inject real stars — the
Balrog paper is explicit that 1/3 of injections (`weight_scheme == Y3_weight`,
confirmed at 33.3% in the file) apply no star/galaxy selection at all, and that
Y6 uses *real* deep-field stars rather than Y3's simulated ones. But the paper's
own star definition is the **colour-based** classifier of the parent Y3
deep-field catalogue (Hartley & Choi et al. 2021), explicitly "classified as
stars using photometric colors (and not morphology)".

That label is recoverable because the injected `id` **is** the deep-field
`COADD_OBJECT_ID` — verified to match 100% of a 5M-injection sample.
`scripts/des/build_des_truth_labels.py` downloads the four deep-field parquets
and builds the lookup. The `KNN_CLASS` encoding is not documented in the
released JSON sidecars (Oracle type info only), so it was established
empirically on bright deep-field objects where morphology is unambiguous:

| `KNN_CLASS` | meaning | median `BDF_T` | frac point-like | median `J-Ks` |
|---|---|---|---|---|
| 2 | **star** | −0.004 | 0.924 | −0.21 |
| 1 | galaxy | 0.891 | 0.015 | +0.50 |
| 0 | **unclassified** | 1.324 | 0.149 | +0.15 |

`KNN_CLASS = 0` is *not* "galaxy" — only 21% of those objects have any NIR
photometry, i.e. the colour classifier had nothing to run on. They are 22.5% of
injections (every injection has `in_VHS_footprint == 1`, so this means "too
faint for usable NIR", not "outside NIR coverage"), and the reducer drops them
from the numerator **and** the denominator of every curve: they are neither
known stars nor known galaxies, so they cannot inform either.

> **Do not substitute a morphological proxy.** A `|bdf_T| < 0.02` cut on the
> injected deep-field morphology selects 7.64% of injections but is only
> **36.8% pure** — per 5M injections it picks up ~141k true stars against ~140k
> galaxies and ~102k unclassified. That is the difference between a stellar
> selection function and a mostly-galaxy one.

DELVE's Balrog-of-the-Stars, by contrast, ships an explicit `truth_STAR` boolean
at roughly 50/50, so none of this applies there.

## Which population defines the depth

The truth anchor asks "at what magnitude does the truth-based scatter of
(obs − true) reach σ = 0.2171?" — but *whose* scatter? The reducer's
`--anchor-sample` selects between the population that has passed the
reference-band S/N cut and the one that has not, and for DES the two answers are
a magnitude apart. Measured on the full catalogue:

| band | `nosnr` shift | `detected` shift |
|---|---|---|
| g | −0.297 | **+0.395** |
| r | −0.280 | **−0.637** |
| i | −0.180 | −0.207 |
| z | −0.124 | −0.116 |

`detected` moves g deeper by 0.40 while moving r shallower by 0.64 — a 1.03 mag
swing between adjacent bands, which cannot be real depth. The cause is
structural: the S/N cut is applied in the **reference band only**, so the g
anchor is truncated by its own cut (self-referential — asking where the scatter
of objects selected for having small scatter reaches a threshold that selection
enforces), while the r/i/z anchor samples are selected on *g* S/N, distorting
each differently. `nosnr` gives smooth, same-sign shifts ordered correctly with
wavelength, and is the default.

The trade-off is that the sample photo-error curve then reaches σ = 0.2171
somewhat faintward of `delta_mag = 0` rather than exactly at it. That is already
true of every truth-anchored streamobs release for a related reason — the maglim
is anchored to the *sample* curve while `get_photo_error` returns the *catalog*
curve — which is why they all carry `skip_snr_maglim_check` in the test registry.

Roman's anchor sample does formally include the detection cut, but its own
methodology footnote records the difference there as ≲ 0.01 mag, i.e. Roman's
anchor is numerically indistinguishable from the no-cut one. DES diverges only
because its S/N cut truncates much harder relative to its error inflation
(1.46 versus Roman's ~2), leaving far less headroom between "reported error
passes S/N > 5" and "truth scatter reaches 0.217". So `nosnr` is the choice that
matches Roman's *effective* convention, even though `detected` matches its code.

> Sanity check to keep: **the per-band anchor shifts must be smooth and
> same-sign.** A band whose shift departs sharply from its neighbours is
> signalling a photometry, selection or sentinel-masking problem in that band,
> not a real depth feature. That check is what caught this.

## Balrog gotchas that cost real time

- **DES ships two files, DELVE one.** `fiducial_matched_measured_sof.hdf5`
  (119,667,013 rows) is exactly `fiducial_injected_sof.hdf5` (145,724,947)
  masked by `detected_0.5arcsec_match_radius`, **in the same row order** — a
  lockstep walk, no hash join. The reducer precomputes the cumulative detected
  count so its chunk reader is idempotent, because the catalogue is traversed
  twice (zero-point audit, then products).
- **Use `wide_tilename`, not `tilename`,** for anything per-tile in DES.
  `tilename` is the *deep-field source* tile — only 181 values (`COSMOS_C*`,
  `SN-*_C*`). DELVE has only `tilename` and it *is* the wide tile.
- **Undetected rows are 0.0, not NaN and not the sentinel,** so they pass every
  `== 0` flag test. The detection gate is what keeps them out.
- **Failed measurements use `-9.999e9`,** coincident with `meas_bdf_flags != 0`.
  Real-catalogue sentinels are heterogeneous and are *not* nulls: `-9999`
  (`CONC`), `-99` (`WAVG_SPREAD_MODEL`), `+99`/`+37.5` (`MAG_AUTO`/capped mags),
  `EXT_XGB = -9`. Value-based masking is mandatory; `isna` silently misses them.
- **`meas_bdf_deblend_flags` is a broken constant** in both surveys; AND-ing it
  zeroes the catalogue. `meas_mask_flags` is identically 0. Use neither.
- **Foreground masking is a property of the sky, not the object.** In DES the
  flag exists only in the measured file, so 7.5% of detected rows carry it while
  undetected rows have no value at all. Cutting on it directly would drop masked
  objects from the numerator but keep them in the denominator. The reducer
  instead rebuilds it as a positional HEALPix mask (nside 4096; a pixel is
  masked if >50% of its detections are flagged, masking 8.9% of populated
  pixels) and applies it to detected and undetected rows alike. In DELVE
  `FLAGS_FOREGROUND` ships for every row — and putting it in the numerator alone
  there drags the bright-end plateau from ~0.95 down to ~0.46.
- **Blends:** the Balrog paper flags injections within 1.5″ of a real object
  (<0.2% of the catalogue) and removes them from every analysis. The reducer
  does the same via `meas_BALROG_FLAG_BLEND`.
- **`flags_bad_zp` does not mean what its name suggests** in DES. It is set on
  11.1% of rows but shows no correlation with per-tile photometric offset
  (Pearson r = −0.002 on an 8% sample), and neither paper defines it. Do not use
  it as a QA cut.
- **Extinction.** DELVE's `truth_FLUX_*` is *dereddened* and needs `A_band`
  added back. DES needs `bdf_mag_deredden + ext_mags` — de-redden at the deep
  field, re-redden at the injection position. The deciding diagnostic is the
  per-tile **spread** of (obs − true), not its median: raw `bdf_mag` has the
  *smallest* median offset yet a 12× wider spread (0.0300 vs 0.0025), because
  extinction varies tile to tile. With the correct choice, 0 of 538 DES tiles
  deviate by more than 0.05 mag.
- **Per-tile zero points differ wildly between the two surveys.** DES is clean;
  DELVE V4 has ~22% of tiles off by >0.05 mag, some by −0.72 (close to
  2.5·log₁₀2, suggesting a factor-of-two flux convention issue on the truth side
  for those tiles). Left uncorrected the DELVE truth anchor returns g ≈ 19.8
  instead of ≈24.5. `--tile-zp reject` is the default and is the treatment DES
  itself applies to this class of artifact.
- **Rejected tiles are a sample cut, not a footprint cut.** The curves are
  `delta_mag`-keyed and the imaging is homogeneous within a survey, so a curve
  measured on good-ZP tiles applies to the whole footprint. Ship the maglim map
  unmasked. Do *not* extrapolate this across surveys.

## Running it

### DES Y6

Three steps; the first two are one-offs.

```bash
PY=/astro/store/shiren/conda-envs/stream_team/envs/streamobs/bin/python

# 1. truth labels (downloads ~5.3 GB of deep-field parquet, emits ~36 MB)
$PY scripts/des/build_des_truth_labels.py --data-dir /path/to/y3_deep_fields

# 2. EXT_XGB surrogate + confusion table (reads the real des_y6_gold catalogue)
$PY scripts/des/build_des_xgb_surrogate.py --out artifacts/des_y6 --precheck

# 3. the reduction itself
$PY scripts/des/balrog_selection_function.py \
    --survey des_y6 --band g \
    --catalog       /path/to/fiducial_injected_sof.hdf5 \
    --measured      /path/to/fiducial_matched_measured_sof.hdf5 \
    --truth-labels  /path/to/des_y6_deepfield_truth.parquet \
    --surrogate     artifacts/des_y6/des_y6_xgb_surrogate.json \
    --confusion     artifacts/des_y6/des_y6_xgb_confusion.csv \
    --maglim-map g=data/surveys/des_yr6/des_y6_5_sig_maglim_band_g_nside_512.hsp \
                 r=data/surveys/des_yr6/des_y6_5_sig_maglim_band_r_nside_512.hsp \
                 i=data/surveys/des_yr6/des_y6_5_sig_maglim_band_i_nside_512.hsp \
                 z=data/surveys/des_yr6/des_y6_5_sig_maglim_band_z_nside_512.hsp \
    --write-maglim --maglim-nside 1024 \
    --corrections config/surveys/des_photoerror_corrections.yaml \
    --out ./des_y6_products --tag des_y6
```

The `des_y6_5_sig_*` maps supply the *spatial structure* only; their absolute
scale is truth-anchored from the injections, so a map delivered on any S/N
convention is corrected automatically rather than by assuming a
+2.5·log₁₀2 = 0.753 mag shift. (Those maps are genuinely 5σ — medians g 25.32 /
r 25.13 / i 24.58 / z 23.88 — despite a long-stale comment in the config that
quoted the 10σ number.)

### DELVE DR3 Gold — runbook for the machine holding the injections

DELVE's products are **not shipped with streamobs**; derive them where the
catalogue lives. The DELVE schema needs no surrogate and no truth-label join —
`truth_STAR` and `FLAGS_FOREGROUND` ship for every row, and the classifier is
`bdf_extended_class_dr3gold`, vendored into the reducer and verified
bit-identical to desqr on 1.5M rows. It needs only `BDF_T` and `BDF_S2N`, both
of which Balrog measures, and it reuses the DES Y6 Gold interpolation nodes, so
DES and DELVE come out classified identically — which is what makes the two
releases comparable.

```bash
D=/path/to/decam/depth/maps
$PY scripts/des/balrog_selection_function.py \
    --survey delve --band g \
    --catalog /path/to/BalrogOfTheStars_Catalog_V4.hdf5 \
    --maglim-map g=$D/delve_dr32_g_maglim_wmean.hsp,$D/delve_dr311+dr312_g_maglim_Nov28th.hsp \
                 r=$D/delve_dr32_r_maglim_wmean.hsp,$D/delve_dr311+dr312_r_maglim_Nov28th.hsp \
    --tile-zp reject --write-maglim --maglim-nside 1024 \
    --corrections config/surveys/delve_photoerror_corrections.yaml \
    --out ./delve_dr3_gold_products --tag delve_dr3_gold --chunk 4000000
```

Roughly 9 minutes and ~38 GB of RAM on the full V4 catalogue.

**Both depth maps are required.** The DR3.2 and DR3.1.1+3.1.2 map sets are
*exactly disjoint* — measured on a 1M-injection sample spread over the whole
catalogue, DR3.2 covers 42.2% and DR3.1.1+3.1.2 covers 57.6%, with **0.00%**
overlap and a 99.76% union. They are complementary halves of the footprint, not
competing vintages. Passing only one silently drops ~half the injections, which
in testing produced an all-zero efficiency table. Comma-separate them; first
listed wins on overlap. They are HealSparse at nside 16384 and are queried
sparsely, never densified — `--maglim-nside` only sets the resolution of maps
*written out*.

Reference band is **g** for both releases; colour is carried by the per-band
depth maps, not by separate per-band efficiency tables.

Check the audit JSON before shipping: the per-tile ZP spread, the truth-anchored
m5 per band (within ~0.1 mag of the map median after the shift), the
error-inflation factor (order 1–2, not 47) and the applied map shift.

Copy back only the CSVs, the depth maps and the audit JSON — they total a few MB.

## Validation

- **Against the published integrated numbers.** Bechtol et al. Table A.3 gives
  stellar efficiency/contamination over two magnitude ranges, which the derived
  curves must reproduce: `0 <= EXT_XGB <= 1` → 98.0% / 4.0% for
  17.5 ≤ i ≤ 22.5 and 94.3% / 12.5% for 16.5 ≤ i ≤ 23.5; `EXT_XGB = 0` →
  92.1% / 1.0% and 79.3% / 1.5%. These constrain the *integral* of the curve;
  the paper publishes per-magnitude performance only graphically (their Fig. 3).
- **Against an external deep truth sample — done, and it passes.**
  `scripts/des/validate_des_xgb_with_splash.py` matches the real Y6 Gold
  catalogue at 0.5″ against SPLASH-SXDF (Mehta et al. 2018) in the same SXDF
  field DES used for their own faint-domain validation. On 151,580
  SPLASH-classified matches the stellar completeness reproduces Table A.3:

  | selection | range | this work | paper |
  |---|---|---|---|
  | `EXT_XGB ≤ 1` | 17.5–22.5 | 0.988 | 0.980 |
  | `EXT_XGB ≤ 1` | 16.5–23.5 | 0.964 | 0.943 |
  | `EXT_XGB = 0` | 17.5–22.5 | 0.947 | 0.921 |
  | `EXT_XGB = 0` | 16.5–23.5 | 0.811 | 0.793 |

  This is the external check on the deconvolution's residual assumption, and it
  is the quantity streamobs actually ships.

  > **Contamination is *not* measurable from this comparison** (0.223 vs the
  > paper's 0.040) and the discrepancy is a property of the truth sample.
  > SPLASH's `STAR_FLAG` is pure but incomplete: flag = 1 sits in a tight
  > point-source locus (median HSC `FLUX_RADIUS` 3.07 vs 6.76 for flag = 0), yet
  > the DES-selected objects it calls "galaxy" have median `FLUX_RADIUS` 3.351
  > with p16 = 3.082 — over half of them lie *on* the stellar locus. They are
  > stars that failed SPLASH's IRAC-dependent star criterion, not DES
  > misclassifications. A pure-but-incomplete star label yields an unbiased
  > `P(selected | star)` and an inflated `P(not star | selected)`. The same
  > limitation means this field cannot validate the galaxy-misclassification
  > product; that needs a *complete* galaxy label, e.g. the HSC PDR3
  > concentration `i_psfflux_mag − i_cmodel_mag` DES used for their Fig. 3.
- **DES against DELVE in `delta_mag` space.** The two footprints overlap over
  only 369 deg² (7.1% of DES) and that overlap is edge-dominated, so a sky-based
  cross-check is weak. Comparing the curves in `delta_mag` needs no sky overlap
  at all, and is precisely what the `delta_mag` keying is for.

> **Bright-end plateau.** In the Roman and LSST products the per-object flag cut
> sits in the numerator and pulls the bright-end plateau down to ~0.91 — a
> property of the adopted cut, not the instrument. **DES Y6 Balrog does not
> behave that way**: its flagged fraction is tiny (0.056% for
> `meas_flags`/`meas_bdf_flags`, plus 0.2% blend-flagged), because injections are
> placed on a sparse 20″ hexagonal grid specifically to avoid injection-injection
> blending. So `detection_eff` plateaus at ~1.0, and that is expected rather than
> a missing cut. Do not "fix" it by reaching for a harsher flag selection, and do
> not read the agreement with the legacy DES product's 1.00 plateau as
> validation — that product reached it for a different reason.
