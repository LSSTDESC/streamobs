#!/usr/bin/env python
"""Build the standard survey-documentation figures for any release.

Every release gets a depth figure; the releases that own their selection
function also get efficiency and photo-error figures:

  1. <TAG>_depth.png       - the depth of each mapped band on the sky, with the
                             distributions below. Small-footprint releases are
                             plotted zoomed rather than all-sky, since a 21 deg^2
                             field is invisible in a Mollweide projection.
  2. <TAG>_efficiency.png  - stellar detection, classification and their
                             product, with the galaxy misclassification rate on
                             the same axes
  3. <TAG>_photoerror.png  - the four photo-error curves: the detected-population
                             sample/catalog pair used for the reference band, and
                             the _nocut pair used for the forced-photometry bands

The lsst_yr*, lsst_dp2 and roman_hlwas_* releases *symlink* their curves from
lsst_dc2 and roman_dc2, so their efficiency and photo-error figures would be
byte-identical duplicates. They are skipped automatically: a release is taken to
own its curves only when the resolved parent directory of its completeness table
is its own resolved directory.

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

plt.rcParams.update(
    {
        "figure.dpi": 130,
        "savefig.dpi": 130,
        "font.size": 10,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "legend.frameon": False,
    }
)

# Bands shown on the depth figure, per survey. The selection function is
# derived in one reference band, so the extra band is there to show how the
# footprint and depth vary with wavelength, not to imply per-band curves.
DEPTH_BANDS = {
    "lsst": ["g", "r"],
    "des": ["g", "r"],
    "delve": ["g", "r"],
    "roman": ["F158", "F106"],
}
# Footprints smaller than this across are plotted zoomed instead of all-sky.
ZOOM_MAX_RADIUS_DEG = 15.0

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
    return np.genfromtxt(
        path, delimiter=",", names=header.split(","), skip_header=n_comment
    )


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

    eff = find("completeness")
    owns_curves = bool(eff) and eff.resolve().parent == ddir.resolve()

    survey = cfg.get("name", release)
    want = DEPTH_BANDS.get(survey, [ref])
    depth_bands = [b for b in want if find(f"maglim_map_{b}")]
    if not depth_bands and find(f"maglim_map_{ref}"):
        depth_bands = [ref]

    return {
        "release": release,
        "name": survey,
        "ref": ref,
        "owns_curves": owns_curves,
        "depth_bands": depth_bands,
        "bands": props.get("bands", []),
        "eff": eff,
        "mis": find("gal_misclassification"),
        "pe_sample": find("log_photo_error_sample"),
        "pe_catalog": find("log_photo_error_catalog"),
        "pe_sample_nc": find("log_photo_error_sample_nocut"),
        "pe_catalog_nc": find("log_photo_error_catalog_nocut"),
        "maglim_ref": find(f"maglim_map_{ref}"),
        **{f"maglim_{b}": find(f"maglim_map_{b}") for b in depth_bands},
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

    ax.plot(
        eff["delta_mag"], eff["detection_eff"], lw=2, color=C_DET, label="detection"
    )
    ax.plot(
        eff["delta_mag"],
        eff["classification_eff"],
        lw=2,
        color=C_CLS,
        label=r"classification $|$ detected",
    )
    ax.plot(
        eff["delta_mag"],
        eff["classification_detection_eff"],
        lw=2.4,
        color=C_BOTH,
        label=r"detection $\times$ classification",
    )
    if r["mis"] is not None:
        mis = read_curve(r["mis"])
        ax.plot(
            mis["delta_mag"],
            mis["missclassification_eff"],
            lw=1.8,
            color=C_MIS,
            ls="--",
            label="galaxy misclassified as star",
        )

    # where the combined efficiency crosses 50%
    d, y = eff["delta_mag"], eff["classification_detection_eff"]
    o = np.argsort(d)
    d, y = d[o], y[o]
    below = np.where(y < 0.5)[0]
    if below.size and below[0] > 0:
        i = below[0]
        x50 = np.interp(0.5, [y[i], y[i - 1]], [d[i], d[i - 1]])
        ax.plot([x50], [0.5], "o", ms=5, color=C_BOTH, zorder=5)
        ax.annotate(
            f"50% at $\\Delta$mag {x50:+.2f}",
            xy=(x50, 0.5),
            xytext=(0.60, 0.72),
            textcoords=ax.transAxes,
            fontsize=8,
            color="0.3",
            arrowprops=dict(arrowstyle="->", color="0.5", lw=0.8),
        )

    ax.axvline(0.0, color="0.5", lw=0.9, ls=":")
    ax.text(0.02, 0.30, "maglim", fontsize=8, color="0.45", rotation=90)
    ax.set_xlabel(rf"$\Delta$mag = mag$_{ref}^{{\rm true}}$ $-$ maglim$_{ref}$")
    ax.set_ylabel("efficiency")
    ax.set_ylim(-0.02, 1.12)
    ax.set_xlim(max(eff["delta_mag"].min(), -8.5), min(eff["delta_mag"].max(), 2.5))
    ax.set_title(
        f"{r['release']} — stellar efficiency and galaxy " f"misclassification ({ref})",
        fontsize=10,
    )
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
    ax.semilogy(
        s["delta_mag"],
        10 ** s["log_mag_err"],
        lw=2,
        color=C_SAMPLE,
        label=f"sample — truth scatter ({ref}, reference band)",
    )
    ax.semilogy(
        c["delta_mag"],
        10 ** c["log_mag_err"],
        lw=2,
        color=C_CATALOG,
        label=f"catalog — reported magerr ({ref}, reference band)",
    )
    if has_nc:
        snc, cnc = read_curve(r["pe_sample_nc"]), read_curve(r["pe_catalog_nc"])
        ax.semilogy(
            snc["delta_mag"],
            10 ** snc["log_mag_err"],
            lw=1.5,
            color=C_SAMPLE,
            ls="--",
            label="sample, no S/N cut — forced bands",
        )
        ax.semilogy(
            cnc["delta_mag"],
            10 ** cnc["log_mag_err"],
            lw=1.5,
            color=C_CATALOG,
            ls="--",
            label="catalog, no S/N cut — forced bands",
        )

    ax.axhline(SIG_SN5, color="0.5", lw=0.9, ls=":")
    ax.text(
        0.5,
        0.94,
        r"$\sigma = 0.217$ (S/N $= 5$)",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=8,
        color="0.4",
    )
    ax.axvline(0.0, color="0.5", lw=0.9, ls=":")
    ax.set_xlabel(rf"$\Delta$mag = mag$_{ref}^{{\rm true}}$ $-$ maglim$_{ref}$")
    ax.set_ylabel(r"$\sigma_{\rm mag}$")
    sub = (
        "solid: reference band, conditioned on detection.  "
        "dashed: forced-photometry bands"
        if has_nc
        else "no _nocut curves shipped for this release"
    )
    ax.set_title(f"{r['release']} — photometric error model\n{sub}", fontsize=10)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    out = r["out"] / f"{r['release']}_photoerror.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"    wrote {out.name}")


def _read_dense(path):
    """Return (dense RING map, nside) for either a .hsp or a .fits.gz map."""
    import healpy as hp

    f = str(path)
    if f.endswith(".hsp"):
        import healsparse as hsp

        m = hsp.HealSparseMap.read(f)
        return m.generate_healpix_map(nside=m.nside_sparse, nest=False), m.nside_sparse
    dense = hp.read_map(f, dtype=np.float64)
    return dense, hp.npix2nside(dense.size)


def _footprint_extent(nside, good_pix):
    """Centre and angular radius of a footprint, safe across the RA=0 wrap.

    Uses the mean unit vector rather than a mean of angles, so a footprint
    spanning RA 0->360 (the HLWAS wide tier) does not report a meaningless
    centre.
    """
    import healpy as hp

    vec = np.array(hp.pix2vec(nside, good_pix))
    mean = vec.mean(axis=1)
    norm = np.linalg.norm(mean)
    if norm == 0:  # antipodally symmetric: no useful centre
        return None, 180.0
    mean = mean / norm
    cosang = np.clip(vec.T @ mean, -1.0, 1.0)
    radius = float(np.degrees(np.arccos(cosang).max()))
    lon, lat = hp.vec2ang(mean, lonlat=True)
    return (float(lon[0]), float(lat[0])), radius


def fig_depth(r):
    """Depth on the sky for each mapped band, with the distributions below."""
    import healpy as hp

    bands = [b for b in r["depth_bands"] if r.get(f"maglim_{b}") is not None]
    if not bands:
        print("    skip depth: no maglim map")
        return

    maps = {}
    for b in bands:
        dense, nside = _read_dense(r[f"maglim_{b}"])
        good = np.isfinite(dense) & (dense > 0) & (dense != hp.UNSEEN)
        if not good.any():
            continue
        maps[b] = {
            "dense": np.where(good, dense, hp.UNSEEN),
            "nside": nside,
            "good": good,
            "v": dense[good],
            "area": good.sum() * hp.nside2pixarea(nside, degrees=True),
        }
    if not maps:
        print("    skip depth: no valid pixels")
        return
    bands = list(maps)

    # zoom when the footprint is compact; an all-sky projection wastes the frame
    b0 = bands[0]
    centre, radius = _footprint_extent(maps[b0]["nside"], np.where(maps[b0]["good"])[0])
    zoom = centre is not None and radius < ZOOM_MAX_RADIUS_DEG

    n = len(bands)
    fig = plt.figure(figsize=(5.4 * n, 6.6))

    for i, b in enumerate(bands):
        d = maps[b]
        lo, hi = np.round(np.percentile(d["v"], [2, 98]), 2)
        title = f"{b}   {d['area']:,.0f} deg²   median {np.median(d['v']):.3f}"
        if zoom:
            # gnomview reso is arcmin/pixel, so size the frame to the footprint
            xsize = 480
            reso = max(2.4 * radius * 60.0 / xsize, 0.2)
            hp.gnomview(
                d["dense"],
                fig=fig.number,
                sub=(2, n, i + 1),
                rot=centre,
                xsize=xsize,
                reso=reso,
                min=lo,
                max=hi,
                cmap="viridis",
                title=title,
                unit="mag",
                notext=True,
            )
            hp.graticule(dpar=2, dmer=2, color="0.75", lw=0.4, verbose=False)
        else:
            hp.mollview(
                d["dense"],
                fig=fig.number,
                sub=(2, n, i + 1),
                min=lo,
                max=hi,
                cmap="viridis",
                title=title,
                unit="mag",
            )
            hp.graticule(dpar=30, dmer=60, color="0.7", lw=0.4, verbose=False)

    ax = fig.add_subplot(2, 1, 2)
    for b, c in zip(bands, (C_DET, C_CLS, C_MIS, "#9467bd")):
        ax.hist(maps[b]["v"], bins=70, histtype="step", lw=1.8, color=c, label=b)
        ax.axvline(np.median(maps[b]["v"]), color=c, lw=1.0, ls=":")
    ax.set_xlabel(r"5$\sigma$ point-source depth")
    ax.set_yticks([])
    ax.legend(fontsize=9, title="band", title_fontsize=8)
    ax.grid(alpha=0.25)

    proj = f"zoomed, {radius:.1f}° radius" if zoom else "all-sky"
    fig.suptitle(
        f"{r['release']} — depth by band ({proj}, nside " f"{maps[bands[0]]['nside']})",
        fontsize=11,
        y=1.045,
    )
    out = r["out"] / f"{r['release']}_depth.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"    wrote {out.name}  [{', '.join(bands)}] {proj}")


def main(argv):
    releases = argv or sorted(
        p.stem for p in CONFIG.glob("*.yaml") if "corrections" not in p.stem
    )
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
        tag = "" if r["owns_curves"] else "   [inherits curves -> depth only]"
        print(f"  {rel}  (ref band {r['ref']}){tag}")
        fig_depth(r)
        if r["owns_curves"]:
            fig_efficiency(r)
            fig_photoerror(r)
        built += 1
    print(f"\nbuilt figures for {built} release(s) -> {STATIC}")
    for rel, why in skipped:
        print(f"  skipped {rel}: {why}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
