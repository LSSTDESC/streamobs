# Roman HLWAS Survey Files

This page documents how the Roman High Latitude Wide Area Survey (HLWAS) selection
function shipped with streamobs was derived. The full, executable derivation is in
[`notebooks/create_streamobs_files_hlwas.ipynb`](https://github.com/LSSTDESC/streamobs/blob/main/notebooks/create_streamobs_files_hlwas.ipynb);
the survey is configured by `config/surveys/roman_hlwas.yaml` and used like any other
survey:

```python
from streamobs.surveys import SurveyFactory
survey = SurveyFactory.create_survey("roman", release="hlwas")
```

## Input data

Everything is measured from the **Roman–Rubin DC2 synthetic survey** of
[Troxel et al. (2023)](https://arxiv.org/abs/2209.06829): ~20 deg² of image-level
simulations of the Roman reference High Latitude Imaging Survey at full depth, in the
four bands Y106/J129/H158/F184 (the mock's `H158` is the Roman **F158** filter).
Detection is SExtractor run on a median Y106+J129+H158+F184 coadd; stars come from
Galfast, galaxies from cosmoDC2.

From the released per-tile detection and truth catalogs we built a single
**detection-centric matched catalog** (`roman_dc2_det_truth.parquet`, built by
`notebooks/build_roman_dc2_det_truth.py`): every S/N>5 detection is matched to a true
source within 1″, breaking ties among the up-to-3 nearest neighbours by taking the
closest in magnitude — the same recipe as the paper.

## Selections (paper-exact)

Following Troxel et al. §"object selection", all products apply:

1. **`flags == 0`** — removes 32% of detections (we reproduce this fraction exactly);
2. **S/N > 5** in the detection image;
3. **matched to a true object** within 1″.

Two data subtleties discovered during the derivation, both handled in the notebook:

- **Duplicate truth entries.** ~26% of truth stars have a second truth entry at the
  *exact same position* under a different index, carrying only the J129 magnitude.
  The match can assign the detection to the duplicate, which would mimic a flat ~12%
  detection-efficiency deficit. Truth stars are therefore collapsed to unique
  positions, and a star counts as detected if *any* of its entries matched.
- **Tile-margin geometry.** Each detection coadd extends ~30″ beyond its truth tile's
  footprint, so per-tile matching leaves ~21% of detections (the margin band)
  unmatched — they are duplicates of detections in the neighbouring tile. Inside the
  truth footprint the unmatched rate is 0.1–5%, consistent with the paper's 1.7%.
  The matched-only selection removes these margin duplicates.

## Products

All files live in `data/surveys/roman_hlwas/`.

### Stellar completeness — `roman_stellar_efficiency_cutf158.csv`

Columns `mag_f158, delta_mag, detection_eff, classifiction_eff,
classification_detection_eff` (same format as the LSST table, including the
`classifiction_eff` spelling that the loader expects). For true stars binned in true
F158 magnitude:

- `detection_eff` — fraction with a clean (`flags==0`) S/N>5 detection matched to them;
- `classifiction_eff` — fraction of those whose detection has SExtractor
  `class_star > 0.5`;
- `classification_detection_eff` — the product (detected *and* classified).

The bright plateau is ~0.90–0.95 (not 1.0): blended stars whose only detection is
flagged are excluded by the paper's `flags==0` cut. The 50% point of the combined
efficiency is at F158 ≈ 27.2 in the simulation.

### Photometric errors — `roman_photoerror_f158.csv`

Columns `delta_mag, log_mag_err`: the binned median log10 of the observed F158
`magerr_auto` of star-classified objects, against `delta_mag = true mag − maglim` at
each object's position.

### Depth maps — `roman_dc2_maglim_f{106,129,158,184}_nside1024[_5sigps].fits.gz`

HEALPix (nside=1024, ring) magnitude-limit maps over the DC2 footprint, computed with
the recipe of [desqr/depth.py](https://github.com/kadrlica/desqr/blob/main/desqr/depth.py):
estimate the slope of log(magerr) vs mag from nearest-neighbour pairs, extrapolate each
object to the magnitude where S/N=5, and take the median per pixel.

Two conventions are provided:

- the raw **measured** maps (median F158 depth 27.11 — the catalog's SExtractor errors
  run ~0.2 mag past the official depths, plausibly because correlated noise in the
  resampled median coadds is underestimated);
- the **`_5sigps`** maps, renormalized so the median equals the official 5σ
  point-source depth of the simulated survey (26.9 in F106/F129/F158, 26.2 in F184).

## The depth convention (important)

The `delta_mag` axes of both tables are keyed to the **official 5σ point-source
depth** convention, *not* to the raw measured depth. This means the tables must be
paired with maglim maps in the same convention — the `_5sigps` maps above, or, for the
real HLWAS footprint, an exposure-time-scaled map normalized to the
[STScI community-defined HLWAS median 5σ point-source depths](https://roman-docs.stsci.edu/roman-community-defined-surveys/high-latitude-wide-area-survey):

| Tier | Area | 5σ point-source total depth (AB) |
|---|---|---|
| Wide | ~2700 deg² | F158 26.2 |
| Medium | ~2400 deg² | F106 26.5, F129 26.4, F158 26.4 |
| Deep | ~19.2 deg² | F106 27.7, F129 27.6, F158 27.5, F184 27.0, F213 25.9 |

With this convention, swapping in a shallower (wide-tier) maglim map automatically
shifts the completeness and photo-error curves to the correct depths, under the
assumption that the selection function depends on magnitude only through
`mag − maglim`.

## Caveats

- The simulation is the *reference* HLIS design (deeper than the current
  community-defined wide tier); the wide-tier behaviour is obtained by translation in
  `delta_mag`, not by an independent simulation.
- Star/galaxy classification is SExtractor `class_star` on the detection coadd; the
  false-star rate for compact (<0.3″) galaxies is <1% brighter than F158 ≈ 25.5,
  rising to ~15–17% at the faint end.
- The detection-centric match assigns a blend to its dominant source, so the
  detection efficiency is slightly conservative for stars blended with brighter
  neighbours.
- The truth catalog only extends to ~15 AB and the simulation does not model
  saturation; the configured `saturation: 15.0` marks the validity limit of the
  table rather than the physical WFI saturation (~17–18 AB for point sources).
- Extinction coefficients are CCM89 (R_V=3.1) evaluated at the filter effective
  wavelengths.
