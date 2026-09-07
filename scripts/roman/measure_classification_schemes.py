#!/usr/bin/env python
"""Compare star/galaxy CLASSIFICATION SCHEMES for the joint Roman x Rubin channel.

WHY
---
With real data the natural selection is "does ROMAN call this a star?", ignoring LSST
morphology entirely and using LSST only for the g photometry that forms the g-F158
colour. The pipeline as written instead requires point-like classification in BOTH
surveys (``fig4_metricC2_ingredients._cmd_cols``, ``which == "joint"``, ANDs the two
``flag_observed`` columns, each of which already folds in classification).

The two schemes differ in three ways that matter, and this script measures all three on
the same truth-centric sample:

1. **Contamination.** How many compact galaxies leak into the star sample.
2. **Completeness.** How many stream stars survive.
3. **Coupling.** Whether the per-survey selection draws can be treated as independent.
   Under Roman-only classification only ONE survey classifies, so the
   ``p_lsst * p_roman`` product that the coupling question is about never arises for the
   classification step -- the residual coupling is only
   ``det_LSST x (det & class)_Roman``.

The headline trade is that keeping both classifications buys sensitivity but obliges the
injector to model a mild coupling (``classification_coupling``), while Roman-only
classification makes independence exactly right but costs S/sqrt(B).

WHAT IS MEASURED
----------------
Truth-centric: every true source in the Roman x LSST overlap, whether or not either
survey detected it -- so detection-side differences between the schemes are captured
(a scheme that ignores LSST classification still admits objects LSST detected but would
never have classified). Per scheme, per source type, per Delta-m selection:

    p_lsst   = P(the LSST-side requirement passes)
    p_roman  = P(the Roman-side requirement passes)
    p_joint  = P(both pass)  <- the selection probability the injector must reproduce
    pred_indep / pred_comono, r_indep / r_comono, fitted Gaussian-copula rho

The schemes, as (LSST-side requirement, Roman-side requirement):

* ``both_classify``    -- (detected AND classified, detected AND classified). Current code.
* ``roman_only_class`` -- (detected only,          detected AND classified). The scheme
  under consideration: LSST contributes photometry, not morphology.
* ``lsst_only_channel``   -- LSST-only channel, single survey, for reference.
* ``roman_only_channel``  -- Roman-only channel, single survey, for reference.

The two single-survey channels have one classifier each and therefore no coupling term
at all; they are reported so the joint schemes can be compared against them.

The Roman side always applies ``flags == 0``, because BOTH shipped Roman products do --
the misclassification curve (``build_roman_galaxy_misclass.py``) and the detection
efficiency (``create_streamobs_files_hlwas.py``, section 2: "and passed the paper's
``flags == 0`` cut"). See ``measure_joint_detection.py`` for why that choice materially
changes the galaxy answer.

Reuses the cached truth catalogs and estimators built by ``measure_joint_detection.py``
and ``measure_joint_misclassification.py``; run those first (the caches make this cheap).

Run (streamobs env):
    /astro/store/shiren/conda-envs/stream_team/envs/streamobs/bin/python \
        scripts/roman/measure_classification_schemes.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ART = REPO / "artifacts"
OUT_CSV = ART / "classification_scheme_comparison.csv"
OUT_META = ART / "classification_scheme_comparison.meta.json"

sys.path.insert(0, str(Path(__file__).resolve().parent))
import measure_joint_detection as D  # noqa: E402
from measure_joint_misclassification import cell_predictions, wilson  # noqa: E402

CROSS_RADIUS = 0.2  # arcsec; see measure_joint_detection.cross_match for the histogram
N_JOINT_MIN = 300   # min joint successes before r_indep/rho are trustworthy


def main():
    ART.mkdir(exist_ok=True)
    t0 = time.time()

    roman = D.build_roman_truth_catalog()
    lsst = D.build_lsst_truth_catalog()
    sample = D.cross_match(lsst, roman, CROSS_RADIUS)

    ml_f158 = D.ref_maglim(D.ROMAN_MAGLIM)[0]
    ml_r = D.ref_maglim(D.LSST_MAGLIM)[0]
    if "dm_roman" not in sample.columns:
        sample["dm_roman"] = sample["mag_H158_true"].to_numpy() - ml_f158
        sample["dm_lsst"] = sample["mag_r"].to_numpy() - ml_r
    dm_l = sample["dm_lsst"].to_numpy()
    dm_r = sample["dm_roman"].to_numpy()
    # Delta-m is the conditioning variable on both axes; a source with no truth flux in
    # a band cannot be placed in a cell (see cell_predictions, which refuses non-finite).
    finite = np.isfinite(dm_l) & np.isfinite(dm_r)

    rgs = sample["roman_gal_star"].to_numpy()
    ltt = sample["truth_type"].to_numpy()
    is_gal = (rgs == 0) & (ltt == D.TRUTH_GAL) & finite
    is_star = (rgs == 1) & (ltt == D.TRUTH_STAR) & finite
    print(
        f"truth-agreeing: {int(is_gal.sum()):,} galaxies, {int(is_star.sum()):,} stars "
        f"({int((~finite).sum()):,} dropped for non-finite Delta-m)"
    )

    size = D.join_true_sizes_full(
        sample["cosmodc2_id"].to_numpy(), sample["cosmodc2_hp"].to_numpy()
    )
    compact = np.isfinite(size) & (size < D.GAL_SIZE_MAX)
    print(f"compact fraction among true galaxies: {compact[is_gal].mean():.4f}")

    ld = sample["lsst_detected"].to_numpy()
    lc = sample["lsst_star_class"].to_numpy()
    rd = sample["roman_detected_flagcut"].to_numpy()
    rc = sample["roman_star_class_flagcut"].to_numpy()
    magH = sample["mag_H158_true"].to_numpy()

    # (LSST-side requirement, Roman-side requirement, has a cross-survey coupling term)
    SCHEMES = {
        "both_classify": (ld & lc, rd & rc, True),
        "roman_only_class": (ld, rd & rc, True),
        "lsst_only_channel": (ld & lc, None, False),
        "roman_only_channel": (None, rd & rc, False),
    }
    POPS = {
        "all_galaxies": is_gal,
        "compact_galaxies": is_gal & compact,
        "stars": is_star,
    }
    # Scan the Roman bright cut: Roman's misclassifications are concentrated faintward,
    # so how much LSST's classifier still adds is a strong function of where you cut.
    F158_CUTS = (26.5, 26.0, 25.5, 25.0, 24.5, 24.0)
    SELS = {"ALL": np.ones(len(sample), bool)}
    SELS.update({f"F158_lt_{c}": magH < c for c in F158_CUTS})

    rows = []
    for pop, popmask in POPS.items():
        for scheme, (fl, fr, coupled) in SCHEMES.items():
            for sel, selmask in SELS.items():
                m = popmask & selmask
                n = int(m.sum())
                if n < 200:
                    continue
                if not coupled:  # single-survey channel: the selection IS the one flag
                    one = (fl if fl is not None else fr)[m]
                    k = int(one.sum())
                    pj = k / n
                    lo, hi = wilson(k, n)
                    rows.append(
                        dict(
                            source_type=pop, scheme=scheme, selection=sel, n=n,
                            p_lsst=float(fl[m].mean()) if fl is not None else np.nan,
                            p_roman=float(fr[m].mean()) if fr is not None else np.nan,
                            p_select=pj, p_select_lo=float(lo), p_select_hi=float(hi),
                            pred_indep=np.nan, pred_comono=np.nan,
                            r_indep=np.nan, r_comono=np.nan, rho_copula=np.nan,
                            n_cells=0, coupled=False,
                        )
                    )
                    continue
                a, b = fl[m], fr[m]
                k = int((a & b).sum())
                pj = k / n
                cp = cell_predictions(dm_l[m], dm_r[m], a, b)
                # The cell-weighted null needs enough JOINT successes to be meaningful:
                # with a handful of events spread over ~1000 cells, sum(w * p_l * p_r) is
                # dominated by empty cells and r_indep becomes noise (measured: r_indep
                # swings to 0.39 with 10 events at F158<24.5). Blank the coupling stats
                # rather than emit a number that invites over-reading; the rates are still
                # reported. See N_JOINT_MIN.
                if k < N_JOINT_MIN:
                    cp = dict(pred_indep=np.nan, pred_comono=np.nan, rho=np.nan,
                              n_cells=cp["n_cells"])
                lo, hi = wilson(k, n)
                rows.append(
                    dict(
                        source_type=pop, scheme=scheme, selection=sel, n=n,
                        n_joint=k,
                        p_lsst=float(a.mean()), p_roman=float(b.mean()),
                        p_select=pj, p_select_lo=float(lo), p_select_hi=float(hi),
                        pred_indep=cp["pred_indep"], pred_comono=cp["pred_comono"],
                        r_indep=pj / cp["pred_indep"] if cp["pred_indep"] > 0 else np.nan,
                        r_comono=pj / cp["pred_comono"] if cp["pred_comono"] > 0 else np.nan,
                        rho_copula=cp["rho"], n_cells=cp["n_cells"], coupled=True,
                    )
                )

    tab = pd.DataFrame(rows)

    # ---- derived sensitivity comparison, joint schemes only ---------------------- #
    # S/sqrt(B): signal = stream stars kept, background = galaxy contaminants admitted.
    # Ratios are taken against `both_classify` so the cost of dropping LSST's
    # classification is explicit. Uses all_galaxies (the population the paper's joint
    # background actually injects: gal_idx = ~is_star, no compact cut).
    deriv = []
    for sel in SELS:
        def rate(pop, scheme):
            q = tab[(tab.source_type == pop) & (tab.scheme == scheme) & (tab.selection == sel)]
            return float(q.p_select.iloc[0]) if len(q) else np.nan

        ref_s, ref_b = rate("stars", "both_classify"), rate("all_galaxies", "both_classify")
        for scheme in ("both_classify", "roman_only_class"):
            s, b = rate("stars", scheme), rate("all_galaxies", scheme)
            snr = (s / ref_s) / np.sqrt(b / ref_b)
            deriv.append(
                dict(
                    selection=sel, scheme=scheme,
                    stars_rel=s / ref_s, contaminants_rel=b / ref_b, snr_rel=snr,
                    delta_sb_mag=-2.5 * np.log10(snr),
                    r_indep_galaxies=float(
                        tab[(tab.source_type == "all_galaxies") & (tab.scheme == scheme)
                            & (tab.selection == sel)].r_indep.iloc[0]),
                    rho_galaxies=float(
                        tab[(tab.source_type == "all_galaxies") & (tab.scheme == scheme)
                            & (tab.selection == sel)].rho_copula.iloc[0]),
                    r_indep_stars=float(
                        tab[(tab.source_type == "stars") & (tab.scheme == scheme)
                            & (tab.selection == sel)].r_indep.iloc[0]),
                    rho_stars=float(
                        tab[(tab.source_type == "stars") & (tab.scheme == scheme)
                            & (tab.selection == sel)].rho_copula.iloc[0]),
                )
            )
    dtab = pd.DataFrame(deriv)
    dtab["kind"] = "derived_sensitivity"
    tab["kind"] = "scheme_rate"
    out = pd.concat([tab, dtab], ignore_index=True)
    out.to_csv(OUT_CSV, index=False, float_format="%.6g")
    print(f"\nwrote {OUT_CSV.relative_to(REPO)} ({len(out)} rows)")

    OUT_META.write_text(json.dumps(dict(
        script=str(Path(__file__).relative_to(REPO)),
        n_cross_matched=int(len(sample)),
        cross_match_radius_arcsec=CROSS_RADIUS,
        ref_maglim_f158=ml_f158, ref_maglim_r=ml_r,
        compact_fraction_true_galaxies=float(compact[is_gal].mean()),
        roman_side="detected AND flags==0 (both shipped Roman products apply it)",
        lsst_classifier=f"extendedness < {D.EXT_CUT}",
        schemes={k: ("LSST det+class" if v[0] is not None and k != "roman_only_class"
                     else "LSST det only" if v[0] is not None else "n/a",
                     "Roman det+class" if v[1] is not None else "n/a")
                 for k, v in SCHEMES.items()},
        snr_definition="signal = stream stars kept; background = all-galaxy contaminants "
                       "admitted; ratios relative to both_classify",
    ), indent=2))

    pd.set_option("display.width", 200, "display.max_columns", 40)
    print("\n===== scheme rates =====")
    print(tab[["source_type", "scheme", "selection", "n", "p_lsst", "p_roman",
               "p_select", "pred_indep", "r_indep", "r_comono", "rho_copula"]]
          .to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print("\n===== derived sensitivity (relative to both_classify) =====")
    print(dtab[["selection", "scheme", "stars_rel", "contaminants_rel", "snr_rel",
                "delta_sb_mag", "r_indep_galaxies", "rho_galaxies"]]
          .to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print(f"\ntotal {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
