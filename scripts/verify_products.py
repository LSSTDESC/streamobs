#!/usr/bin/env python
"""Verify the streamobs selection-function products for DES Y6 and DELVE DR3 Gold.

Runs the checks that the product contract in
``docs/source/selection_function_methodology.md`` implies but that no single
test asserts end to end, and emits a manifest suitable for a Zenodo record.

Five groups of checks:

  MANIFEST     size + sha256 of every shipped product
  CONTRACT     column names, value ranges, the faint clamp, curve orderings
  DEPTH        each maglim map loads, and its median matches the audit JSON
  PHYSICS      anchor shifts, error-inflation factor, bright plateau, the
               delta_mag at which combined efficiency crosses 50%
  CROSS        DES vs DELVE in delta_mag space -- the comparison the
               methodology calls the main cross-check, since in the DES
               footprint DELVE DR3 Gold *is* DES Y6 and the two curves must
               agree there

Exit status is 0 only if every check passes.  Figures are written to --figdir.

Usage
-----
    python verify_products.py \
        --release des_yr6=/path/to/data/surveys/des_yr6 \
        --release delve_dr3_gold=/path/to/products \
        --figdir figs/verification
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

DET_EFF_DELTA_MAX = 1.0  # faint clamp from the product contract
SIG_SN5 = 2.5 / np.log(10) / 5.0  # 0.21715 mag

results = []  # (group, name, passed, detail)


def check(group, name, passed, detail=""):
    results.append((group, name, bool(passed), detail))
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    return passed


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
def find_products(d: Path, tag: str):
    """Map logical product -> path, tolerating the two naming tags in use."""
    def one(pattern):
        hits = sorted(d.glob(pattern))
        return hits[0] if hits else None

    return {
        "efficiency": one(f"{tag}_stellar_efficiency_cut?.csv"),
        "photoerror_sample": one(f"{tag}_photoerror_?.csv"),
        "photoerror_catalog": one(f"{tag}_photoerror_?_catalog.csv"),
        "photoerror_sample_nocut": one(f"{tag}_photoerror_?_nocut.csv"),
        "photoerror_catalog_nocut": one(f"{tag}_photoerror_?_catalog_nocut.csv"),
        "misclass": one(f"{tag}_galaxy_misclass_cut?.csv"),
        "audit": one(f"{tag}_audit.json"),
        "maglim": sorted(d.glob(f"{tag}_maglim_*_nside*.fits.gz")),
    }


def verify_contract(tag, prod):
    """Column names, ranges, the faint clamp, and curve orderings."""
    eff = pd.read_csv(prod["efficiency"])

    ref = None
    for c in eff.columns:
        if c.startswith("mag_"):
            ref = c.split("_", 1)[1]
    check("CONTRACT", f"{tag}: efficiency has a mag_<band> column", ref is not None,
          f"band = {ref}")

    want = {f"mag_{ref}", "delta_mag", "detection_eff", "classification_eff",
            "classification_detection_eff"}
    check("CONTRACT", f"{tag}: efficiency columns exact", set(eff.columns) == want,
          f"{sorted(eff.columns)}")
    # the legacy misspelling must not have crept back in
    check("CONTRACT", f"{tag}: 'classification_eff' spelled correctly",
          "classifiction_eff" not in eff.columns)

    for c in ["detection_eff", "classification_eff", "classification_detection_eff"]:
        v = eff[c].to_numpy()
        check("CONTRACT", f"{tag}: {c} within [0,1]",
              np.all((v >= 0) & (v <= 1)), f"min {v.min():.4f} max {v.max():.4f}")
        check("CONTRACT", f"{tag}: {c} finite", np.all(np.isfinite(v)))

    # classification_detection_eff = det*cls, so it can never exceed detection_eff
    bad = int(np.sum(eff["classification_detection_eff"] > eff["detection_eff"] + 1e-9))
    check("CONTRACT", f"{tag}: classification_detection_eff <= detection_eff", bad == 0,
          f"{bad} violating rows")

    # the faint clamp
    faint = eff["delta_mag"] > DET_EFF_DELTA_MAX
    if faint.any():
        z = eff.loc[faint, ["detection_eff", "classification_detection_eff"]].to_numpy()
        check("CONTRACT", f"{tag}: faint clamp zeroes delta_mag > {DET_EFF_DELTA_MAX}",
              np.all(z == 0), f"{int(faint.sum())} clamped rows")
    else:
        check("CONTRACT", f"{tag}: efficiency grid stops at the clamp", True,
              "no rows faintward of the clamp")

    for key in ("photoerror_sample", "photoerror_catalog",
                "photoerror_sample_nocut", "photoerror_catalog_nocut"):
        pe = pd.read_csv(prod[key])
        check("CONTRACT", f"{tag}: {key} columns exact",
              list(pe.columns) == ["delta_mag", "log_mag_err"], f"{list(pe.columns)}")
        check("CONTRACT", f"{tag}: {key} finite and sorted",
              np.all(np.isfinite(pe.to_numpy())) and pe["delta_mag"].is_monotonic_increasing)

    mis = pd.read_csv(prod["misclass"])
    check("CONTRACT", f"{tag}: misclass uses the key streamobs reads",
          "missclassification_eff" in mis.columns, f"{list(mis.columns)}")

    # the two photo-error curves must share a delta_mag grid so they can be
    # interpolated against each other
    a = pd.read_csv(prod["photoerror_sample"])["delta_mag"].to_numpy()
    b = pd.read_csv(prod["photoerror_catalog"])["delta_mag"].to_numpy()
    check("CONTRACT", f"{tag}: sample and catalog share a delta_mag grid",
          a.shape == b.shape and np.allclose(a, b))

    # The sample (truth) scatter must exceed the reported error -- the whole
    # reason the two-curve model exists.  Only assert it brightward of the
    # depth: faintward, the detected population is truncated (only objects that
    # scattered bright are recovered), so the *measured* truth scatter
    # compresses while the reported magerr keeps rising, and the two curves can
    # invert.  That is the "curves are measured on the detected population"
    # effect in the methodology doc, not a product defect -- so it is reported
    # rather than failed.
    s = pd.read_csv(prod["photoerror_sample"])["log_mag_err"].to_numpy()
    c = pd.read_csv(prod["photoerror_catalog"])["log_mag_err"].to_numpy()
    bright = a <= 0.0
    frac = float(np.mean(s[bright] >= c[bright] - 1e-9))
    check("CONTRACT", f"{tag}: truth scatter >= reported error (delta_mag <= 0)",
          frac > 0.99, f"{100*frac:.1f}% of {int(bright.sum())} bins")

    # The forced-photometry pair must exist for a multi-band release: streamobs
    # raises rather than applying the detected-population curve to a band whose
    # photometry was forced.
    nc_s = pd.read_csv(prod["photoerror_sample_nocut"])
    nc_c = pd.read_csv(prod["photoerror_catalog_nocut"])
    check("CONTRACT", f"{tag}: _nocut curves share the cut curves' delta_mag grid",
          len(nc_s) == len(nc_c) and np.allclose(nc_s["delta_mag"], nc_c["delta_mag"]))

    # They must differ from the cut pair, and only faintward: the reference-band
    # S/N cut cannot change the scatter of objects well above it.
    both = pd.merge(pd.read_csv(prod["photoerror_sample"]), nc_s,
                    on="delta_mag", suffixes=("_cut", "_nc"))
    d = both["log_mag_err_nc"] - both["log_mag_err_cut"]
    check("CONTRACT", f"{tag}: _nocut differs from the detected-population curve",
          not np.allclose(d, 0.0), f"max |diff| {np.abs(d).max():.4f} dex")
    bright = both["delta_mag"] < -0.5
    check("CONTRACT", f"{tag}: _nocut agrees brightward of the depth",
          bool(np.allclose(d[bright], 0.0, atol=1e-6)),
          f"{int(bright.sum())} bins with delta_mag < -0.5")
    # Faintward the S/N cut truncates the detected sample, so its measured
    # scatter is the *smaller* of the two.
    faint = both["delta_mag"] > 0.5
    if faint.any():
        check("CONTRACT", f"{tag}: _nocut scatter >= cut scatter faintward",
              bool((d[faint] >= -1e-6).all()),
              f"min diff {d[faint].min():+.4f} dex over {int(faint.sum())} bins")

    inv = a[(s < c - 1e-9)]
    if inv.size:
        print(f"       note: sample < catalog in {inv.size} bin(s), "
              f"delta_mag {inv.min():+.2f} to {inv.max():+.2f} "
              f"(faint truncation; see methodology doc)")
    return ref, eff


def verify_depth(tag, prod, audit):
    import healpy as hp

    for path in prod["maglim"]:
        band = path.name.split("_maglim_")[1].split("_")[0]
        m = hp.read_map(str(path), dtype=np.float64)
        nside = hp.npix2nside(m.size)
        mask = np.isfinite(m) & (m > 0)
        good = m[mask]
        med = float(np.median(good))
        area = mask.sum() * hp.nside2pixarea(nside, degrees=True)
        check("DEPTH", f"{tag}: maglim {band} loads", good.size > 0,
              f"nside {nside}, {good.size:,} valid px, {area:,.0f} deg^2, median {med:.3f}")

        # the filename advertises a resolution; the map must actually have it.
        # to_healpix only ever *degrades*, so requesting --maglim-nside above
        # the input map's own nside silently leaves the data coarser than the
        # name claims.
        claimed = path.name.split("_nside")[1].split(".")[0]
        check("DEPTH", f"{tag}: maglim {band} nside matches its filename",
              str(nside) == claimed, f"file says nside{claimed}, map is nside {nside}")
        ref = audit.get("m5_truth_anchored", {}).get(band)
        if ref is not None:
            # the written map is masked at >1.5 mag from the median, so it moves
            # slightly off the anchor; 0.05 mag is the tolerance the run reports
            check("DEPTH", f"{tag}: maglim {band} median matches the anchor",
                  abs(med - ref) < 0.05, f"map {med:.3f} vs anchor {ref:.3f}")


def verify_physics(tag, eff, ref, audit):
    infl = audit.get("error_inflation_factor")
    check("PHYSICS", f"{tag}: error-inflation factor is physical",
          infl is not None and 1.0 <= infl <= 3.0, f"{infl:.3f}")

    check("PHYSICS", f"{tag}: depth is truth-anchored", audit.get("truth_anchored") is True,
          f"anchor sample = {audit.get('anchor_sample', 'n/a')}")

    bright = eff[eff["delta_mag"] < -3.0]
    plateau = float(bright["detection_eff"].median()) if len(bright) else float("nan")
    check("PHYSICS", f"{tag}: bright-end detection plateau is high",
          plateau > 0.85, f"{plateau:.3f}")

    # combined efficiency should cross 50% within a magnitude of the 5-sigma depth
    d = eff["delta_mag"].to_numpy()
    y = eff["classification_detection_eff"].to_numpy()
    order = np.argsort(d)
    d, y = d[order], y[order]
    cross = np.nan
    below = np.where(y < 0.5)[0]
    if below.size and below[0] > 0:
        i = below[0]
        cross = float(np.interp(0.5, [y[i], y[i - 1]], [d[i], d[i - 1]]))
    check("PHYSICS", f"{tag}: combined efficiency crosses 50% near delta_mag 0",
          np.isfinite(cross) and -1.0 < cross < 1.0, f"crossing at delta_mag = {cross:+.3f}")
    return plateau, cross


def verify_cross(rel, figdir):
    """DES vs DELVE on a common delta_mag axis."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not {"des_yr6", "delve_dr3_gold"} <= set(rel):
        check("CROSS", "both releases present for the cross-check", False,
              "need des_yr6 and delve_dr3_gold")
        return

    grid = np.arange(-6.0, 1.01, 0.05)
    curves, pe = {}, {}
    for tag in ("des_yr6", "delve_dr3_gold"):
        eff = rel[tag]["eff"]
        o = np.argsort(eff["delta_mag"].to_numpy())
        curves[tag] = np.interp(grid, eff["delta_mag"].to_numpy()[o],
                                eff["classification_detection_eff"].to_numpy()[o])
        p = pd.read_csv(rel[tag]["prod"]["photoerror_sample"])
        pe[tag] = np.interp(grid, p["delta_mag"], p["log_mag_err"],
                            left=np.nan, right=np.nan)

    band = (grid > -4.0) & (grid < 0.0)
    d_eff = np.abs(curves["des_yr6"] - curves["delve_dr3_gold"])[band]
    check("CROSS", "efficiency curves agree over -4 < delta_mag < 0",
          float(np.median(d_eff)) < 0.15,
          f"median |diff| {np.median(d_eff):.3f}, max {np.max(d_eff):.3f}")

    m = band & np.isfinite(pe["des_yr6"]) & np.isfinite(pe["delve_dr3_gold"])
    d_pe = np.abs(pe["des_yr6"] - pe["delve_dr3_gold"])[m]
    check("CROSS", "photo-error curves agree over -4 < delta_mag < 0",
          d_pe.size > 0 and float(np.median(d_pe)) < 0.3,
          f"median |diff| {np.median(d_pe):.3f} dex" if d_pe.size else "no overlap")

    figdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for tag, c in (("des_yr6", "C0"), ("delve_dr3_gold", "C1")):
        ax[0].plot(grid, curves[tag], c, lw=1.8, label=tag)
        ax[1].plot(grid, pe[tag], c, lw=1.8, label=tag)
    ax[0].axhline(0.5, color="0.6", lw=0.8, ls=":")
    ax[0].axvline(0.0, color="0.6", lw=0.8, ls=":")
    ax[0].set(xlabel=r"$\Delta$mag", ylabel="classification $\\times$ detection eff",
              title="Combined stellar efficiency")
    ax[1].axhline(np.log10(SIG_SN5), color="0.6", lw=0.8, ls=":")
    ax[1].axvline(0.0, color="0.6", lw=0.8, ls=":")
    ax[1].set(xlabel=r"$\Delta$mag", ylabel=r"$\log_{10}\sigma_{\rm mag}$ (truth scatter)",
              title="Photometric error (sample curve)")
    for a in ax:
        a.legend(frameon=False, fontsize=9)
        a.grid(alpha=0.25)
    fig.suptitle("DES Y6 vs DELVE DR3 Gold in $\\Delta$mag space", y=1.02)
    fig.tight_layout()
    out = figdir / "des_vs_delve_delta_mag.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  wrote {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--release", action="append", required=True,
                    metavar="TAG=DIR", help="repeatable, e.g. des_yr6=/path/to/dir")
    ap.add_argument("--figdir", default="figs/verification")
    ap.add_argument("--manifest", default=None, help="write the manifest as JSON here")
    args = ap.parse_args()

    rel, manifest = {}, {}
    for spec in args.release:
        tag, _, d = spec.partition("=")
        d = Path(d)
        prod = find_products(d, tag)
        missing = [k for k, v in prod.items() if not v]
        print(f"\n=== {tag}  ({d}) ===")
        if not check("MANIFEST", f"{tag}: all products present", not missing,
                     f"missing {missing}" if missing else ""):
            continue

        files = [v for k, v in prod.items() if k != "maglim"] + list(prod["maglim"])
        manifest[tag] = {}
        for f in sorted(files):
            manifest[tag][f.name] = {"bytes": f.stat().st_size, "sha256": sha256(f)}
        print(f"  {len(files)} files, "
              f"{sum(v['bytes'] for v in manifest[tag].values())/1e6:.2f} MB")

        audit = json.loads(prod["audit"].read_text())
        ref, eff = verify_contract(tag, prod)
        verify_depth(tag, prod, audit)
        verify_physics(tag, eff, ref, audit)
        rel[tag] = {"prod": prod, "eff": eff, "audit": audit, "ref": ref}

    print("\n=== cross-survey ===")
    verify_cross(rel, Path(args.figdir))

    if args.manifest:
        Path(args.manifest).write_text(json.dumps(manifest, indent=2, sort_keys=True))
        print(f"  wrote {args.manifest}")

    n_fail = sum(1 for *_, ok, _ in [(g, n, p, d) for g, n, p, d in results] if not ok)
    print(f"\n{'=' * 60}\n{len(results) - n_fail}/{len(results)} checks passed")
    if n_fail:
        print("\nFAILURES:")
        for g, n, ok, d in results:
            if not ok:
                print(f"  [{g}] {n} -- {d}")
    return 1 if n_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
