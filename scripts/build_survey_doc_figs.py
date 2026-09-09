#!/usr/bin/env python
"""Build the standard survey-documentation figures for any release.

Every survey page carries the same three figures, so a reader can compare
releases at a glance without learning a new layout each time:

  1. <TAG>_efficiency.png  - stellar detection, classification and their
                             product, with the galaxy misclassification rate on
                             the same axes
  2. <TAG>_photoerror.png  - the four photo-error curves: the detected-population
                             sample/catalog pair used for the reference band, and
                             the _nocut pair used for the forced-photometry bands
  3. <TAG>_depth.png       - the reference-band depth on the sky next to its
                             distribution, i.e. footprint and depth in one figure

Everything is driven from ``config/surveys/<release>.yaml``: the products are
located by the same keys streamobs itself reads, so a new release picks up
figures as soon as it has a config, with no per-survey code.

Survey-specific extras live in their own scripts and are not replaced by this
one -- scripts/des/build_des_survey_doc_figs.py (EXT_XGB surrogate confusion,
SPLASH validation), scripts/roman/build_roman_survey_doc_figs.py.

Usage
-----
  python scripts/build_survey_doc_figs.py                    # every release
  python scripts/build_survey_doc_figs.py des_yr6 lsst_dp2   # named releases
"""

import pathlib
import sys

import matplotlib as mpl
import numpy as np
import yaml

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parents[1]
CONFIG = REPO / "config/surveys"
DATA = REPO / "data/surveys"
STATIC = REPO / "docs/source/_static"
OTHERS = DATA.parent / "others"

SIG_SN5 = 2.5 / np.log(10) / 5.0  # 0.21715

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 130, "font.size": 10,
    "axes.grid": True, "grid.alpha": 0.25, "legend.frameon": False,
})

C_DET, C_CLS, C_BOTH, C_MIS = "#1f77b4", "#d62728", "k", "#2ca02c"
C_SAMPLE, C_CATALOG = "#d62728", "#1f77b4"


def read_curve(path):
    """Read a product CSV; the last '#'-prefixed line (if any) is the header."""
    path = pathlib.Path(path)
    header, n_comment = None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                header, n_comment = line.lstrip("#").strip(), n_comment + 1
            else:
                break
    if header is None:
        return np.genfromtxt(path, delimiter=",", names=True)
    return np.genfromtxt(path, delimiter=",", names=header.split(","),
                         skip_header=n_comment)


def load_release(release):
    """Resolve a release's products from its config. Returns None if unusable."""
    cfg_path = CONFIG / f"{release}.yaml"
    if not cfg_path.exists():
        return None
    cfg = yaml.safe_load(open(cfg_path))
    files = cfg.get("survey_files", {})
    props = cfg.get("survey_properties", {})
    ddir = DATA / release

    ref = files.get("completeness_band")
    if ref is None:
        return None

    def find(key):
        name = files.get(key)
        if not name:
            return None
        for cand in (ddir / name, OTHERS / name):
            if cand.exists():
                return cand
        return None

    return {
        "release": release,
        "name": cfg.get("name", release),
        "ref": ref,
        "bands": props.get("bands", []),
        "eff": find("completeness"),
        "mis": find("gal_misclassification"),
        "pe_sample": find("log_photo_error_sample"),
        "pe_catalog": find("log_photo_error_catalog"),
        "pe_sample_nc": find("log_photo_error_sample_nocut"),
        "pe_catalog_nc": find("log_photo_error_catalog_nocut"),
        "maglim_ref": find(f"maglim_map_{ref}"),
        "out": STATIC / release,
    }


# ---------------------------------------------------------------------------
def fig_efficiency(r):
    if r["eff"] is None:
        print("    skip efficiency: no completeness table")
        return
    eff = read_curve(r["eff"])
    ref = r["ref"]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))

    ax.plot(eff["delta_mag"], eff["detection_eff"], lw=2, color=C_DET,
            label="detection")
    ax.plot(eff["delta_mag"], eff["classification_eff"], lw=2, color=C_CLS,
            label=r"classification $|$ detected")
    ax.plot(eff["delta_mag"], eff["classification_detection_eff"], lw=2.4,
            color=C_BOTH, label=r"detection $\times$ classification")
    if r["mis"] is not None:
        mis = read_curve(r["mis"])
        ax.plot(mis["delta_mag"], mis["missclassification_eff"], lw=1.8,
                color=C_MIS, ls="--", label="galaxy misclassified as star")

    # where the combined efficiency crosses 50%
    d, y = eff["delta_mag"], eff["classification_detection_eff"]
    o = np.argsort(d)
    d, y = d[o], y[o]
    below = np.where(y < 0.5)[0]
    if below.size and below[0] > 0:
        i = below[0]
        x50 = np.interp(0.5, [y[i], y[i - 1]], [d[i], d[i - 1]])
        ax.plot([x50], [0.5], "o", ms=5, color=C_BOTH, zorder=5)
        ax.annotate(f"50% at $\\Delta$mag {x50:+.2f}", xy=(x50, 0.5),
                    xytext=(0.60, 0.72), textcoords=ax.transAxes, fontsize=8,
                    color="0.3",
                    arrowprops=dict(arrowstyle="->", color="0.5", lw=0.8))

    ax.axvline(0.0, color="0.5", lw=0.9, ls=":")
    ax.text(0.02, 0.30, "maglim", fontsize=8, color="0.45", rotation=90)
    ax.set_xlabel(rf"$\Delta$mag = mag$_{ref}^{{\rm true}}$ $-$ maglim$_{ref}$")
    ax.set_ylabel("efficiency")
    ax.set_ylim(-0.02, 1.12)
    ax.set_xlim(max(eff["delta_mag"].min(), -8.5), min(eff["delta_mag"].max(), 2.5))
    ax.set_title(f"{r['release']} — stellar efficiency and galaxy "
                 f"misclassification ({ref})", fontsize=10)
    ax.legend(loc="lower left", fontsize=8.5)
    fig.tight_layout()
    out = r["out"] / f"{r['release']}_efficiency.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"    wrote {out.name}")


def fig_photoerror(r):
    if r["pe_sample"] is None or r["pe_catalog"] is None:
        print("    skip photoerror: no curves")
        return
    ref = r["ref"]
    s, c = read_curve(r["pe_sample"]), read_curve(r["pe_catalog"])
    has_nc = r["pe_sample_nc"] is not None and r["pe_catalog_nc"] is not None

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.semilogy(s["delta_mag"], 10 ** s["log_mag_err"], lw=2, color=C_SAMPLE,
                label=f"sample — truth scatter ({ref}, reference band)")
    ax.semilogy(c["delta_mag"], 10 ** c["log_mag_err"], lw=2, color=C_CATALOG,
                label=f"catalog — reported magerr ({ref}, reference band)")
    if has_nc:
        snc, cnc = read_curve(r["pe_sample_nc"]), read_curve(r["pe_catalog_nc"])
        ax.semilogy(snc["delta_mag"], 10 ** snc["log_mag_err"], lw=1.5,
                    color=C_SAMPLE, ls="--", label="sample, no S/N cut — forced bands")
        ax.semilogy(cnc["delta_mag"], 10 ** cnc["log_mag_err"], lw=1.5,
                    color=C_CATALOG, ls="--", label="catalog, no S/N cut — forced bands")

    ax.axhline(SIG_SN5, color="0.5", lw=0.9, ls=":")
    ax.text(0.5, 0.94, r"$\sigma = 0.217$ (S/N $= 5$)", transform=ax.transAxes,
            ha="center", va="top", fontsize=8, color="0.4")
    ax.axvline(0.0, color="0.5", lw=0.9, ls=":")
    ax.set_xlabel(rf"$\Delta$mag = mag$_{ref}^{{\rm true}}$ $-$ maglim$_{ref}$")
    ax.set_ylabel(r"$\sigma_{\rm mag}$")
    sub = ("solid: reference band, conditioned on detection.  "
           "dashed: forced-photometry bands" if has_nc else
           "no _nocut curves shipped for this release")
    ax.set_title(f"{r['release']} — photometric error model\n{sub}", fontsize=10)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    out = r["out"] / f"{r['release']}_photoerror.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"    wrote {out.name}")


def fig_depth(r):
    if r["maglim_ref"] is None:
        print("    skip depth: no reference-band maglim map")
        return
    import healpy as hp

    f = str(r["maglim_ref"])
    if f.endswith(".hsp"):
        import healsparse as hsp
        m = hsp.HealSparseMap.read(f)
        nside = m.nside_sparse
        dense = m.generate_healpix_map(nside=nside, nest=False)
    else:
        dense = hp.read_map(f, dtype=np.float64)
        nside = hp.npix2nside(dense.size)

    good = np.isfinite(dense) & (dense > 0) & (dense != hp.UNSEEN)
    v = dense[good]
    med = float(np.median(v))
    area = good.sum() * hp.nside2pixarea(nside, degrees=True)
    lo, hi = np.round(np.percentile(v, [2, 98]), 2)
    masked = np.where(good, dense, hp.UNSEEN)

    fig = plt.figure(figsize=(10.4, 4.0))
    hp.mollview(masked, fig=fig.number, sub=(1, 2, 1), cmap="viridis",
                min=lo, max=hi, unit="mag", title="footprint and depth",
                cbar=True)
    hp.graticule(dpar=30, dmer=60, color="0.7", lw=0.4, verbose=False)

    ax = fig.add_subplot(1, 2, 2)
    ax.hist(v, bins=70, color=C_DET, alpha=0.85)
    ax.axvline(med, color="k", lw=1.6, label=f"median {med:.3f}")
    ax.set_xlabel(f"5$\\sigma$ {r['ref']}-band depth")
    ax.set_yticks([])
    ax.set_title(f"{area:,.0f} deg², nside {nside}", fontsize=9)
    ax.legend(fontsize=8.5, loc="upper left")
    ax.grid(alpha=0.25)

    fig.suptitle(f"{r['release']} — {r['ref']}-band depth", fontsize=11, y=1.02)
    fig.tight_layout()
    out = r["out"] / f"{r['release']}_depth.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"    wrote {out.name}  ({area:,.0f} deg², median {med:.3f}, nside {nside})")


def main(argv):
    releases = argv or sorted(p.stem for p in CONFIG.glob("*.yaml")
                              if "corrections" not in p.stem)
    built, skipped = 0, []
    for rel in releases:
        r = load_release(rel)
        if r is None:
            skipped.append((rel, "no config or no completeness_band"))
            continue
        if r["eff"] is None and r["maglim_ref"] is None:
            skipped.append((rel, "no products on disk"))
            continue
        r["out"].mkdir(parents=True, exist_ok=True)
        print(f"  {rel}  (ref band {r['ref']})")
        fig_efficiency(r)
        fig_photoerror(r)
        fig_depth(r)
        built += 1
    print(f"\nbuilt figures for {built} release(s) -> {STATIC}")
    for rel, why in skipped:
        print(f"  skipped {rel}: {why}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
