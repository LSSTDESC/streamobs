# DES

**DES** is supported by **StreamObs**.

# Available releases

## DES Y6 Gold

The survey dataset is described in
[Bechtol et al. 2025](https://arxiv.org/abs/2501.05739) and the catalogues are
documented and publicly available from
[DESDM](https://des.ncsa.illinois.edu/releases). Load it with:

```python
des_yr6 = surveys.Survey.load(survey='des', release='yr6')
```

### Products

All of the selection-function products are **derived in-repo** from the DES Y6
Balrog synthetic-source-injection catalogue
([Anbajagane et al. 2025](https://arxiv.org/abs/2501.05683)) by
`scripts/des/balrog_selection_function.py`. The method, the evidence behind each
choice and the runbook are in :doc:`../balrog_selection_functions`; this page is
just the per-release summary.

| File | Contents |
|---|---|
| `des_yr6_maglim_{g,r,i,z}_nside1024.fits.gz` | truth-anchored S/N = 5 depth maps |
| `des_yr6_stellar_efficiency_cutg.csv` | stellar detection + classification efficiency vs `delta_mag` |
| `des_yr6_photoerror_g.csv` | **sample** photo-error curve — truth scatter, drives the noise draw |
| `des_yr6_photoerror_g_catalog.csv` | **catalog** photo-error curve — reported `magerr`, drives the S/N cut |
| `des_yr6_galaxy_misclass_cutg.csv` | fraction of detected true galaxies classified as point sources |
| `des_yr6_photoerror_g{,_catalog}_raw.csv` | pre-afterburner provenance |
| `des_yr6_audit.json` | counts, anchors and convention flags for the run |

Reference band is **g**. The completeness and photo-error curves are keyed to
`delta_mag = mag_true − maglim(pixel)` and applied band-independently, so colour
is carried by the per-band depth maps rather than by separate per-band tables.

### Depth

Truth-anchored medians, and the shift applied to each input map:

| band | input map median | truth-anchored | shift |
|---|---|---|---|
| g | 25.321 | **25.025** | −0.297 |
| r | 25.130 | **24.850** | −0.280 |
| i | 24.582 | **24.402** | −0.180 |
| z | 23.878 | **23.754** | −0.124 |

The input `des_y6_5_sig_maglim_band_*_nside_512.hsp` maps supply the spatial
structure and footprint (5220 deg²); the absolute scale comes from the
injections. Those input maps are **genuinely 5σ**, contrary to a long-standing
stale comment in the config that quoted the 10σ number — the measured medians
sit ~0.75 mag (= 2.5·log₁₀2) fainter than the published Y6 10σ depths, exactly
as the conversion predicts.

The shifts are smooth, same-sign and ordered with wavelength. That coherence is
a validation check, not a coincidence — see the technote for the anchor
convention it rules out.

### Bands: griz only

**Y is not supported.** The Y6 Balrog injects *griz*, so a Y depth map cannot be
truth-anchored on the same footing and no Y efficiency or photo-error curve can
be derived at all. The legacy Y-band healsparse map was removed rather than
shipped on an unanchored scale. `u` was never part of DES.

### Star/galaxy classification, and the EXT_XGB caveat

`classification_eff` describes the **`0 ≤ EXT_XGB ≤ 1`** ("complete") stellar
selection — the cut you would apply to the real Y6 Gold catalogue.

`EXT_XGB` **cannot be evaluated on injected sources**: three of its six features
(`CONC`, `WAVG_SPREAD_MODEL_I`, `WAVG_SPREADERR_MODEL_I`) are never measured for
them, which Bechtol et al. state outright (App. A.2) and Anbajagane et al.
reaffirm (Sec. 4.4). Neither paper offers a workaround; the Balrog paper falls
back to `EXT_MASH`, which it shows carries >10% galaxy contamination for stellar
selection.

streamobs instead trains a **surrogate** for `EXT_XGB` on the real `des_y6_gold`
catalogue using only Balrog-computable features, then **deconvolves** the
Balrog-measured efficiency with the surrogate's per-magnitude confusion so the
shipped curve describes `EXT_XGB` rather than the surrogate. Surrogate
performance on 16.8M real rows: AUC 0.9946, selection agreement 0.9802.

> **What this means for you.** The curve is calibrated for `0 ≤ EXT_XGB ≤ 1`. If
> you apply a *different* stellar cut to real data — `EXT_MASH`, `EXT_FITVD`, or
> `EXT_XGB = 0` — the shipped `classification_eff` does not describe your
> selection and will bias any completeness correction. Re-derive with the
> matching `--ext-max`, or use a different classifier consistently on both sides.

The residual systematic is that the confusion terms are measured on a mixed
star+galaxy population but applied to true stars. It is validated externally
against SPLASH-SXDF in the same field DES used for their own classifier
validation; see the technote.

### Photometric errors

The two-curve model matters here: the truth-based scatter runs **~1.46×** the
reported `magerr`, so the noise draw must not use the reported errors. Brightward
of `delta_mag ≈ −7` (g ≈ 18) the sample curve is dropped as quantisation noise —
see `config/surveys/des_photoerror_corrections.yaml` for the measured evidence.

`sys_error = 0.005` is retained and contributes 3.1% in quadrature against the
curve's 0.020 floor.

### Known limitations

- **Bright end is thin.** Balrog injects each of 2.8M deep-field objects ~51
  times, so the effective sample size in a magnitude bin is the number of
  *distinct* parent sources, not rows. Bright bins are backed by few distinct
  deep-field stars — which are additionally saturated in the much deeper
  deep-field imaging, corrupting their injected morphology. Bins with too few
  distinct sources are dropped, so the curve simply starts where it is
  statistically meaningful.
- **`detection_eff` plateaus near 1.0**, unlike Roman/LSST (~0.91). This is
  expected: the DES Balrog flagged fraction is tiny because injections sit on a
  sparse 20″ grid chosen to avoid injection blending. It is not a missing cut.
- **The truth star label is colour-based**, from the parent Y3 deep-field
  catalogue, not morphological. A morphological proxy on the same data is only
  ~37% pure.
- **`flags_bad_zp` is not usable** as a QA cut — it shows no correlation with
  per-tile photometric offset and neither paper defines it.

Questions about these files can be addressed to Peter Ferguson.
