#!/usr/bin/env python
"""Train a Balrog-computable surrogate for the DES Y6 Gold EXT_XGB star selection.

Why this exists
---------------
The Y6 Gold star/galaxy classifier ``EXT_XGB`` cannot be evaluated on Balrog
injections.  Three of its six input features (``CONC``, ``WAVG_SPREAD_MODEL_I``,
``WAVG_SPREADERR_MODEL_I``) are never measured for injected sources -- stated
outright in Bechtol et al. 2025 (arXiv:2501.05739) App. A.2, *"the XGBoost
classification has not been implemented on simulation-injection-recovery tests
due to the fact that the CONC parameter was not measured for the simulated
object samples"*, and reaffirmed in Anbajagane et al. 2025 (arXiv:2501.05683)
Sec. 4.4.  Neither paper offers a workaround; the Balrog paper simply falls back
to ``EXT_MASH``, which it shows carries >10% galaxy contamination when selecting
stars.

But the selection function needs
``classification_eff = P(EXT_XGB <= 1 | true star, detected)`` for the cut that
will actually be applied to the real catalogue.  Balrog has the truth label but
not ``EXT_XGB``; the real catalogue has ``EXT_XGB`` but no truth label.  They
overlap in *feature* space, which is what makes this solvable:

  1. **Surrogate.**  Train ``S`` on the real Y6 Gold catalogue against the label
     ``EXT_XGB <= 1``, using only features that are constructible identically on
     both sides (all from the same fitvd/SOF pipeline -- the Balrog files are
     ``_sof``).
  2. **Deconvolution.**  ``S`` applied to Balrog gives ``eff_S = P(S=1|star)``,
     not ``eff_X``.  Measure the confusion on the real catalogue per magnitude
     bin, ``a(m) = P(S=1 | EXT_XGB<=1)`` and ``b(m) = P(S=1 | EXT_XGB>1)``, and
     invert:  ``eff_X = (eff_S - b) / (a - b)``.
     This removes the first-order bias rather than quoting the surrogate's
     disagreement rate as an irreducible systematic.

``CONC`` reconstruction pre-check (``--precheck``) establishes the headroom: on
10M real rows, ``CONC`` is recovered with R^2 = 0.85-0.94 brightward of i ~ 23
(residual scatter ~20x below its intrinsic spread) falling to 0.39 by i ~ 24.5.
The GAp-based concentration dominates the feature importance, consistent with
the paper's definition of ``CONC`` as a Gaussian-weighted flux *dilation* ratio
(Sec. III.4) rather than an aperture-magnitude difference.

Feature contract (must be identical on both sides)
--------------------------------------------------
    feature          real Y6 Gold                     Balrog matched file
    BDF_T            BDF_T                            meas_bdf_T
    BDF_T_ERR        BDF_T_ERR                        meas_bdf_T_err
    BDF_T_RATIO      BDF_T_RATIO                      meas_bdf_T_ratio
    log_bdf_s2n      log10(BDF_S2N)                   log10(meas_bdf_s2n)
    PSF_T            PSF_T                            meas_psf_T
    psf_bdf_T        PSF_T - BDF_T                    meas_psf_T - meas_bdf_T
    conc_gap         GAP_MAG_I - BDF_MAG_I            meas_gap_mag[:,2] - meas_bdf_mag[:,2]
    conc_gap_r       GAP_MAG_R - BDF_MAG_R            meas_gap_mag[:,1] - meas_bdf_mag[:,1]
    conc_ap8         PSF_MAG_APER_8_I - BDF_MAG_I     meas_psf_mag_aper8[:,2] - meas_bdf_mag[:,2]
    conc_ap8_r       PSF_MAG_APER_8_R - BDF_MAG_R     meas_psf_mag_aper8[:,1] - meas_bdf_mag[:,1]
    BDF_MAG_I        BDF_MAG_I                        meas_bdf_mag[:,2]
    BDF_MAG_G        BDF_MAG_G                        meas_bdf_mag[:,0]

Balrog (N,4) arrays are ordered g,r,i,z -- so index 0=g, 1=r, 2=i.

Sentinels are heterogeneous and are NOT nulls, so value-based masking is
mandatory: -9.999e9 (BDF/GAP/PSF), -9999 (CONC), -99 (WAVG_SPREAD_MODEL),
+99 / +37.5 (MAG_AUTO / capped mags), EXT_XGB = -9 (no data).

Outputs
-------
    <out>/des_y6_xgb_surrogate.json        the trained booster
    <out>/des_y6_xgb_features.txt          ordered feature list (the contract)
    <out>/des_y6_xgb_confusion.csv         mag_g, a, b, n_pos, n_neg  (deconvolution)
    <out>/des_y6_xgb_surrogate_audit.json  metrics
"""

from __future__ import annotations

import argparse
import glob
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

GOLD = "/astro/store/shire/hats/catalogs/des/des_y6_gold/des_y6_gold/dataset"

BAD = -9.0e8      # anything below this in a BDF/GAP/PSF column is a sentinel
CONC_BAD = -9990.0
MAG_CAP = 37.0    # 37.5 is the "capped" magnitude value; 99 is MAG_AUTO's null

RAW_COLS = [
    "BDF_T", "BDF_T_ERR", "BDF_T_RATIO", "BDF_S2N", "PSF_T",
    "GAP_MAG_I", "GAP_MAG_R", "PSF_MAG_APER_8_I", "PSF_MAG_APER_8_R",
    "BDF_MAG_I", "BDF_MAG_R", "BDF_MAG_G",
    "EXT_XGB", "XGB_PRED", "CONC", "CONC_FLAGS", "MAG_AUTO_I",
    "FLAGS_FOREGROUND", "FLAGS_GOLD",
]

# Two feature sets.  "full" is everything constructible on both sides; "robust"
# drops the two that the domain-shift check flags and one that is redundant:
#
#   PSF_T      real median 1.035-1.045 vs Balrog 0.882-0.886 -- a constant
#              ~1.4 sigma offset at EVERY magnitude, i.e. a pipeline difference
#              (cf. the PSFEx/Piff mismatch noted in arXiv:2501.05683 Sec. 3.4),
#              not a population difference.  Training on a feature whose
#              distribution is displaced in the apply-domain is unsafe.
#   psf_bdf_T  inherits the same shift (0.70 sigma by i ~ 24.5).
#   conc_gap_r redundant with conc_gap: fitvd's BDF fits ONE shape jointly
#              across bands with per-band fluxes, so (GAP - BDF) is a purely
#              morphological quantity and is band-independent by construction
#              (verified identical wherever both are finite).
#
# "robust" is the default; --feature-set full is kept so the cost of dropping
# them stays measurable.
FEATURES_FULL = [
    "BDF_T", "BDF_T_ERR", "BDF_T_RATIO", "log_bdf_s2n", "PSF_T", "psf_bdf_T",
    "conc_gap", "conc_gap_r", "conc_ap8", "conc_ap8_r",
    "BDF_MAG_I", "BDF_MAG_G",
]
FEATURES_ROBUST = [
    "BDF_T", "BDF_T_ERR", "BDF_T_RATIO", "log_bdf_s2n",
    "conc_gap", "conc_ap8", "conc_ap8_r",
    "BDF_MAG_I", "BDF_MAG_G",
]
FEATURES = FEATURES_ROBUST  # rebound in main() from --feature-set

# Confusion table binning -- g band, to match the reducer's reference band and
# the `mag_g` column of the stellar-efficiency CSV.
MAG_BINS = np.arange(16.0, 27.0 + 1e-9, 0.25)


def derive(d: pd.DataFrame) -> pd.DataFrame:
    """Mask sentinels and build the derived features."""
    for c in RAW_COLS:
        if c in ("EXT_XGB", "XGB_PRED", "CONC", "CONC_FLAGS",
                 "FLAGS_FOREGROUND", "FLAGS_GOLD"):
            continue
        d.loc[d[c] < BAD, c] = np.nan
    for c in ("GAP_MAG_I", "GAP_MAG_R", "PSF_MAG_APER_8_I", "PSF_MAG_APER_8_R",
              "BDF_MAG_I", "BDF_MAG_R", "BDF_MAG_G"):
        d.loc[d[c] > MAG_CAP, c] = np.nan
    d.loc[d["CONC"] < CONC_BAD, "CONC"] = np.nan
    d.loc[d["MAG_AUTO_I"] > 90, "MAG_AUTO_I"] = np.nan

    with np.errstate(all="ignore"):
        d["log_bdf_s2n"] = np.log10(d["BDF_S2N"].where(d["BDF_S2N"] > 0))
    d["psf_bdf_T"] = d["PSF_T"] - d["BDF_T"]
    d["conc_gap"] = d["GAP_MAG_I"] - d["BDF_MAG_I"]
    d["conc_gap_r"] = d["GAP_MAG_R"] - d["BDF_MAG_R"]
    d["conc_ap8"] = d["PSF_MAG_APER_8_I"] - d["BDF_MAG_I"]
    d["conc_ap8_r"] = d["PSF_MAG_APER_8_R"] - d["BDF_MAG_R"]
    return d


def load_gold(n_part: int, seed: int) -> pd.DataFrame:
    files = sorted(glob.glob(f"{GOLD}/**/*.parquet", recursive=True))
    rng = np.random.default_rng(seed)
    pick = rng.choice(len(files), size=min(n_part, len(files)), replace=False)
    frames = []
    for i in pick:
        try:
            frames.append(pq.ParquetFile(files[i]).read(columns=RAW_COLS).to_pandas())
        except Exception as exc:  # a few partitions have odd schemas
            print(f"  skip {Path(files[i]).name}: {type(exc).__name__}")
    d = pd.concat(frames, ignore_index=True)
    print(f"  read {len(d):,} rows from {len(frames)}/{len(files)} partitions")
    return derive(d)


def usable(d: pd.DataFrame) -> np.ndarray:
    ok = (d["EXT_XGB"] >= 0).to_numpy()  # -9 = no classification
    for f in FEATURES:
        ok &= np.isfinite(d[f]).to_numpy()
    return ok


def precheck(d: pd.DataFrame) -> dict:
    """How much of the real CONC is recoverable from Balrog-available features?"""
    import xgboost as xgb
    from sklearn.model_selection import train_test_split

    ok = usable(d) & np.isfinite(d["CONC"]).to_numpy() & (d["CONC_FLAGS"] == 0).to_numpy()
    dd = d[ok]
    print(f"\nCONC pre-check on {len(dd):,} rows")
    X, y, m = dd[FEATURES].to_numpy(), dd["CONC"].to_numpy(), dd["BDF_MAG_I"].to_numpy()
    Xtr, Xte, ytr, yte, _, mte = train_test_split(X, y, m, test_size=0.3, random_state=42)
    reg = xgb.XGBRegressor(n_estimators=400, max_depth=7, learning_rate=0.08,
                           subsample=0.8, colsample_bytree=0.8,
                           tree_method="hist", n_jobs=16)
    reg.fit(Xtr, ytr)
    p = reg.predict(Xte)
    r2 = 1 - np.sum((yte - p) ** 2) / np.sum((yte - yte.mean()) ** 2)
    print(f"  global R2 = {r2:.4f}")
    per = {}
    for lo in np.arange(18, 25, 1.0):
        s = (mte >= lo) & (mte < lo + 1)
        if s.sum() < 500:
            continue
        r = yte[s] - p[s]
        per[f"{lo:.0f}-{lo+1:.0f}"] = {
            "n": int(s.sum()),
            "r2": float(1 - np.sum(r**2) / np.sum((yte[s] - yte[s].mean()) ** 2)),
            "resid": float((np.percentile(r, 84) - np.percentile(r, 16)) / 2),
            "conc_spread": float(
                (np.percentile(yte[s], 84) - np.percentile(yte[s], 16)) / 2),
        }
        v = per[f"{lo:.0f}-{lo+1:.0f}"]
        print(f"  i {lo:.0f}-{lo+1:.0f}: n={v['n']:>7d} R2={v['r2']:+.3f} "
              f"resid={v['resid']:.4f} spread={v['conc_spread']:.4f}")
    return {"global_r2": float(r2), "per_mag": per}


def train(d: pd.DataFrame, ext_max: int, seed: int):
    import xgboost as xgb
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import train_test_split

    ok = usable(d)
    dd = d[ok].reset_index(drop=True)
    y = (dd["EXT_XGB"] <= ext_max).to_numpy().astype(int)
    print(f"\ntraining on {len(dd):,} rows; positives (EXT_XGB<={ext_max}): "
          f"{y.mean():.4f}")

    idx = np.arange(len(dd))
    itr, ite = train_test_split(idx, test_size=0.3, random_state=seed, stratify=y)
    clf = xgb.XGBClassifier(
        n_estimators=600, max_depth=8, learning_rate=0.06,
        subsample=0.8, colsample_bytree=0.8, min_child_weight=10,
        eval_metric="logloss", tree_method="hist", n_jobs=16, random_state=seed,
    )
    clf.fit(dd.loc[itr, FEATURES].to_numpy(), y[itr])

    Xte = dd.loc[ite, FEATURES].to_numpy()
    prob = clf.predict_proba(Xte)[:, 1]
    pred = (prob >= 0.5).astype(int)
    auc = roc_auc_score(y[ite], prob)
    agree = float((pred == y[ite]).mean())
    tp = int(((pred == 1) & (y[ite] == 1)).sum())
    recall = tp / max(int((y[ite] == 1).sum()), 1)
    precision = tp / max(int((pred == 1).sum()), 1)
    print(f"  AUC={auc:.4f}  agreement={agree:.4f}  "
          f"recall={recall:.4f}  precision={precision:.4f}")

    return clf, dd, ite, y, prob, {
        "auc": float(auc), "agreement": agree,
        "recall": float(recall), "precision": float(precision),
        "n_train": int(len(itr)), "n_test": int(len(ite)),
        "positive_rate": float(y.mean()),
    }


def confusion_table(dd, ite, y, prob, thresh: float) -> pd.DataFrame:
    """a(m) = P(S=1 | EXT_XGB<=max), b(m) = P(S=1 | EXT_XGB>max), per mag_g bin."""
    mag = dd.loc[ite, "BDF_MAG_G"].to_numpy()
    s = (prob >= thresh).astype(int)
    yt = y[ite]
    rows = []
    for lo, hi in zip(MAG_BINS[:-1], MAG_BINS[1:]):
        m = (mag >= lo) & (mag < hi)
        pos, neg = m & (yt == 1), m & (yt == 0)
        rows.append({
            "mag_g": 0.5 * (lo + hi),
            "a": float(s[pos].mean()) if pos.sum() else np.nan,
            "b": float(s[neg].mean()) if neg.sum() else np.nan,
            "n_pos": int(pos.sum()),
            "n_neg": int(neg.sum()),
        })
    return pd.DataFrame(rows)


def domain_check(dd: pd.DataFrame, balrog: str, n_rows: int) -> dict:
    """Compare real-catalogue vs Balrog feature distributions at matched mag.

    A shift here invalidates applying the surrogate across the domain boundary,
    so it is measured before anything downstream is built.
    """
    import h5py

    print(f"\ndomain-shift check against {Path(balrog).name}")
    with h5py.File(balrog, "r") as f:
        sl = slice(0, n_rows)
        gap = f["meas_gap_mag"][sl]
        bdf = f["meas_bdf_mag"][sl]
        ap8 = f["meas_psf_mag_aper8"][sl]
        b = pd.DataFrame({
            "BDF_T": f["meas_bdf_T"][sl],
            "BDF_T_ERR": f["meas_bdf_T_err"][sl],
            "BDF_T_RATIO": f["meas_bdf_T_ratio"][sl],
            "BDF_S2N": f["meas_bdf_s2n"][sl],
            "PSF_T": f["meas_psf_T"][sl],
            "BDF_MAG_G": bdf[:, 0], "BDF_MAG_R": bdf[:, 1], "BDF_MAG_I": bdf[:, 2],
            "GAP_MAG_R": gap[:, 1], "GAP_MAG_I": gap[:, 2],
            "PSF_MAG_APER_8_R": ap8[:, 1], "PSF_MAG_APER_8_I": ap8[:, 2],
            "flags": f["meas_flags"][sl], "bdf_flags": f["meas_bdf_flags"][sl],
        })
    for c in b.columns:
        if c in ("flags", "bdf_flags"):
            continue
        b.loc[b[c] < BAD, c] = np.nan
    for c in ("GAP_MAG_I", "GAP_MAG_R", "PSF_MAG_APER_8_I", "PSF_MAG_APER_8_R",
              "BDF_MAG_I", "BDF_MAG_R", "BDF_MAG_G"):
        b.loc[b[c] > MAG_CAP, c] = np.nan
    with np.errstate(all="ignore"):
        b["log_bdf_s2n"] = np.log10(b["BDF_S2N"].where(b["BDF_S2N"] > 0))
    b["psf_bdf_T"] = b["PSF_T"] - b["BDF_T"]
    b["conc_gap"] = b["GAP_MAG_I"] - b["BDF_MAG_I"]
    b["conc_gap_r"] = b["GAP_MAG_R"] - b["BDF_MAG_R"]
    b["conc_ap8"] = b["PSF_MAG_APER_8_I"] - b["BDF_MAG_I"]
    b["conc_ap8_r"] = b["PSF_MAG_APER_8_R"] - b["BDF_MAG_R"]
    b = b[(b["flags"] == 0) & (b["bdf_flags"] == 0)]

    out = {}
    print(f"  {'feature':16s} {'i-bin':>10s} {'real med':>10s} {'balrog med':>11s} "
          f"{'d/spread':>9s}")
    for lo in (20.0, 22.0, 23.0, 24.0):
        rs = (dd["BDF_MAG_I"] >= lo) & (dd["BDF_MAG_I"] < lo + 1)
        bs = (b["BDF_MAG_I"] >= lo) & (b["BDF_MAG_I"] < lo + 1)
        if rs.sum() < 200 or bs.sum() < 200:
            continue
        for f in FEATURES:
            if f in ("BDF_MAG_I", "BDF_MAG_G"):
                continue
            r_v = dd.loc[rs, f].to_numpy()
            b_v = b.loc[bs, f].to_numpy()
            b_v = b_v[np.isfinite(b_v)]
            if b_v.size < 200:
                continue
            r_med, b_med = np.nanmedian(r_v), np.nanmedian(b_v)
            spread = (np.nanpercentile(r_v, 84) - np.nanpercentile(r_v, 16)) / 2
            norm = abs(b_med - r_med) / spread if spread > 0 else np.nan
            out[f"{f}@{lo:.0f}"] = {
                "real_median": float(r_med), "balrog_median": float(b_med),
                "shift_over_spread": float(norm),
            }
            flag = "  <-- SHIFT" if norm > 0.5 else ""
            print(f"  {f:16s} {lo:5.0f}-{lo+1:.0f} {r_med:10.4f} {b_med:11.4f} "
                  f"{norm:9.3f}{flag}")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="artifacts/des_y6")
    ap.add_argument("--n-partitions", type=int, default=40,
                    help="how many des_y6_gold HATS partitions to sample")
    ap.add_argument("--ext-max", type=int, default=1,
                    help="star iff 0 <= EXT_XGB <= this (1 = the 'complete' sample)")
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--feature-set", choices=["robust", "full"], default="robust",
                    help="'robust' drops PSF_T/psf_bdf_T (domain-shifted) and "
                         "conc_gap_r (redundant); see the module comment")
    ap.add_argument("--precheck", action="store_true",
                    help="also run the CONC reconstruction pre-check")
    ap.add_argument(
        "--balrog",
        default="/astro/store/shire/pferguso/des_y6_balrog/fiducial_matched_measured_sof.hdf5",
        help="matched Balrog file for the domain-shift check ('' to skip)")
    ap.add_argument("--balrog-rows", type=int, default=2_000_000)
    args = ap.parse_args()

    global FEATURES
    FEATURES = FEATURES_FULL if args.feature_set == "full" else FEATURES_ROBUST

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    print(f"feature set: {args.feature_set} ({len(FEATURES)} features)")

    print("loading des_y6_gold:")
    d = load_gold(args.n_partitions, args.seed)

    audit: dict = {"ext_max": args.ext_max, "threshold": args.threshold,
                   "feature_set": args.feature_set, "features": FEATURES}
    if args.precheck:
        audit["conc_precheck"] = precheck(d)

    clf, dd, ite, y, prob, metrics = train(d, args.ext_max, args.seed)
    audit["surrogate"] = metrics

    conf = confusion_table(dd, ite, y, prob, args.threshold)
    conf.to_csv(out / "des_y6_xgb_confusion.csv", index=False)
    good = conf[(conf["n_pos"] > 200) & (conf["n_neg"] > 200)]
    print(f"\nconfusion table ({len(good)} usable bins):")
    print(good.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    clf.save_model(out / "des_y6_xgb_surrogate.json")
    (out / "des_y6_xgb_features.txt").write_text("\n".join(FEATURES) + "\n")

    if args.balrog:
        audit["domain_shift"] = domain_check(dd, args.balrog, args.balrog_rows)

    (out / "des_y6_xgb_surrogate_audit.json").write_text(json.dumps(audit, indent=2))
    print(f"\nwrote {out}/des_y6_xgb_surrogate.json, _features.txt, "
          f"_confusion.csv, _audit.json")


if __name__ == "__main__":
    main()
