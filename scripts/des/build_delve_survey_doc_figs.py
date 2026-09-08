"""
build_delve_survey_doc_figs.py
------------------------------
Generate the figures for docs/source/surveys/DELVE.md from the shipped
delve/dr3_gold products. Mirrors scripts/des/build_des_survey_doc_figs.py.

  1. DELVE_efficiencies.png - detection_eff, classification_eff,
                              classification_detection_eff + the galaxy
                              missclassification_eff on the same axes
  2. DELVE_errors.png       - all four photo-error curves: the sample/catalog
                              pair for the reference band and the _nocut pair
                              for the forced-photometry bands, with the S/N=5
                              line and the measured inflation
  3. DELVE_depth.png        - griz truth-anchored depth histograms, with the
                              input-map medians marked so the anchor shift is
                              visible
  4. DELVE_depth_map.png    - the g-band depth map on the sky

There is deliberately no surrogate-confusion figure and no external-validation
figure: DELVE needs no EXT_XGB surrogate (its classifier is computable on the
injections), and no SPLASH-equivalent truth catalogue overlaps the footprint.

All go to docs/source/_static/delve_dr3_gold/

Usage
-----
  python scripts/des/build_delve_survey_doc_figs.py
"""

import pathlib

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parents[2]
DATA = REPO / "data/surveys/delve_dr3_gold"
OUT = REPO / "docs/source/_static/delve_dr3_gold"
OUT.mkdir(parents=True, exist_ok=True)

BANDS = ("g", "r", "i", "z")
# input healsparse map medians (the two disjoint DR3 halves, mosaicked), for the
# anchor-shift annotation
INPUT_MEDIAN = {"g": 24.178, "r": 23.675, "i": 23.204, "z": 22.560}
SIG_SN5 = 2.5 / np.log(10) / 5.0  # 0.21715

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 130, "font.size": 10,
    "axes.grid": True, "grid.alpha": 0.25, "legend.frameon": False,
})


def read_csv(path):
    """Read a product CSV; the last '#'-prefixed line (if any) is the header."""
    path = pathlib.Path(path)
    header = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                header = line.lstrip("#").strip()
            else:
                break
    if header is None:
        return np.genfromtxt(path, delimiter=",", names=True)
    return np.genfromtxt(path, delimiter=",", names=header.split(","),
                         skip_header=sum(1 for _ in open(path)
                                         if _.startswith("#")))


def fig_efficiencies():
    eff = read_csv(DATA / "delve_dr3_gold_stellar_efficiency_cutg.csv")
    mis = read_csv(DATA / "delve_dr3_gold_galaxy_misclass_cutg.csv")

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(eff["delta_mag"], eff["detection_eff"], lw=2, color="#1f77b4",
            label="detection")
    ax.plot(eff["delta_mag"], eff["classification_eff"], lw=2, color="#d62728",
            label=r"classification $|$ detected")
    ax.plot(eff["delta_mag"], eff["classification_detection_eff"], lw=2.4,
            color="k", label="detection $\\times$ classification")
    ax.plot(mis["delta_mag"], mis["missclassification_eff"], lw=1.8,
            color="#2ca02c", ls="--", label="galaxy misclassified as star")

    # the misclassification curve is noise-dominated brightward of -4
    ax.axvspan(mis["delta_mag"].min(), -4.0, color="0.85", alpha=0.5, zorder=0)
    ax.annotate("misclassification\nnoise-dominated:\nfew true galaxies per bin",
                xy=(-4.5, 0.10), xytext=(0.03, 0.30), textcoords=ax.transAxes,
                fontsize=7.5, color="0.3", ha="left",
                arrowprops=dict(arrowstyle="->", color="0.5", lw=0.8))

    # the bright-end plateau sits below unity by construction
    plateau = float(np.nanmedian(eff["detection_eff"][eff["delta_mag"] < -3.0]))
    ax.axhline(plateau, color="#1f77b4", lw=0.8, ls=":", alpha=0.7)
    ax.text(0.985, 0.86, f"bright-end plateau {plateau:.3f}\nper-object quality "
            "flags are in the\nnumerator (Roman/LSST convention)",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=7.5, color="#1f77b4")

    ax.axvline(0.0, color="0.5", lw=0.9, ls=":")
    ax.text(0.06, 0.30, "maglim", fontsize=8, color="0.45", rotation=90)
    ax.set_xlabel(r"$\Delta$mag = mag$_g^{\rm true}$ $-$ maglim$_g$")
    ax.set_ylabel("efficiency")
    ax.set_ylim(-0.02, 1.12)
    ax.set_xlim(eff["delta_mag"].min(), 2.0)
    ax.set_title("DELVE DR3 Gold stellar efficiency and galaxy "
                 "misclassification\n"
                 r"($0 \leq$ bdf_extended_class_dr3gold $\leq 1$)", fontsize=10)
    ax.legend(loc="center left", fontsize=8.5, bbox_to_anchor=(0.02, 0.52))
    fig.tight_layout()
    fig.savefig(OUT / "DELVE_efficiencies.png")
    plt.close(fig)
    print(f"  wrote {OUT/'DELVE_efficiencies.png'}")


def fig_errors():
    s = read_csv(DATA / "delve_dr3_gold_photoerror_g.csv")
    c = read_csv(DATA / "delve_dr3_gold_photoerror_g_catalog.csv")
    s_nc = read_csv(DATA / "delve_dr3_gold_photoerror_g_nocut.csv")
    c_nc = read_csv(DATA / "delve_dr3_gold_photoerror_g_catalog_nocut.csv")
    sig_s, sig_c = 10 ** s["log_mag_err"], 10 ** c["log_mag_err"]
    sig_s_nc, sig_c_nc = 10 ** s_nc["log_mag_err"], 10 ** c_nc["log_mag_err"]

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7.2, 6.0), sharex=True,
                                  gridspec_kw={"height_ratios": [2.4, 1]})
    ax.semilogy(s["delta_mag"], sig_s, lw=2, color="#d62728",
                label="sample (truth scatter) $\\rightarrow$ noise draw")
    ax.semilogy(c["delta_mag"], sig_c, lw=2, color="#1f77b4",
                label="catalog (reported magerr) $\\rightarrow$ S/N cut")
    ax.semilogy(s_nc["delta_mag"], sig_s_nc, lw=1.5, color="#d62728", ls="--",
                label="sample, no S/N cut (forced bands)")
    ax.semilogy(c_nc["delta_mag"], sig_c_nc, lw=1.5, color="#1f77b4", ls="--",
                label="catalog, no S/N cut (forced bands)")
    ax.axhline(SIG_SN5, color="0.5", lw=0.9, ls=":")
    ax.text(0.985, 0.42, r"$\sigma = 0.217$ (S/N $=5$)", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8, color="0.4")
    ax.axvline(0.0, color="0.5", lw=0.9, ls=":")
    ax.set_ylabel(r"$\sigma_{\rm mag}$")
    ax.set_title("DELVE DR3 Gold photometric error model (g)\n"
                 "solid: reference band, conditioned on detection.  "
                 "dashed: forced-photometry bands", fontsize=10)
    ax.legend(loc="upper left", fontsize=8)
    ax.text(0.985, 0.03,
            "the two pairs agree brightward of the depth and separate only\n"
            "faintward, where the S/N cut truncates the detected sample and\n"
            "its measured scatter turns over instead of continuing to rise",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
            color="0.35")

    n = min(sig_s.size, sig_c.size)
    ax2.plot(s["delta_mag"][:n], sig_s[:n] / sig_c[:n], lw=2, color="k",
             label="reference band")
    m = min(sig_s_nc.size, sig_c_nc.size)
    ax2.plot(s_nc["delta_mag"][:m], sig_s_nc[:m] / sig_c_nc[:m], lw=1.5,
             color="0.45", ls="--", label="forced bands")
    ax2.axhline(1.0, color="0.5", lw=0.9, ls=":")
    ax2.set_ylabel("sample / catalog")
    ax2.set_xlabel(r"$\Delta$mag = mag$_g^{\rm true}$ $-$ maglim$_g$")
    ax2.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(OUT / "DELVE_errors.png")
    plt.close(fig)
    print(f"  wrote {OUT/'DELVE_errors.png'}")


def fig_depth():
    import healpy as hp

    fig, axes = plt.subplots(2, 2, figsize=(8.2, 5.4))
    for ax, b in zip(axes.ravel(), BANDS):
        m = hp.read_map(str(DATA / f"delve_dr3_gold_maglim_{b}_nside512.fits.gz"))
        good = np.isfinite(m) & (m > 0) & (m != hp.UNSEEN)
        v = m[good]
        med = np.median(v)
        ax.hist(v, bins=80, color="#1f77b4", alpha=0.85)
        ax.axvline(med, color="k", lw=1.6, label=f"map median {med:.3f}")
        ax.axvline(INPUT_MEDIAN[b], color="#d62728", lw=1.4, ls="--",
                   label=f"input map {INPUT_MEDIAN[b]:.3f}")
        ax.set_title(f"{b}  (shift {med - INPUT_MEDIAN[b]:+.3f})", fontsize=9)
        ax.set_xlabel("truth-anchored 5$\\sigma$ depth")
        ax.set_yticks([])
        ax.legend(fontsize=7.5, loc="upper left")
    fig.suptitle("DELVE DR3 Gold truth-anchored depth, nside 512 — all four "
                 "shifts share a sign\n(medians are of the written map, which is "
                 "masked >1.5 mag from centre, so they sit ~0.02 below the anchor)",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "DELVE_depth.png")
    plt.close(fig)
    print(f"  wrote {OUT/'DELVE_depth.png'}")

    m = hp.read_map(str(DATA / "delve_dr3_gold_maglim_g_nside512.fits.gz"))
    m = np.where(np.isfinite(m) & (m > 0), m, hp.UNSEEN)
    v = m[m != hp.UNSEEN]
    lo, hi = np.percentile(v, [2, 98])
    fig = plt.figure(figsize=(7.6, 4.4))
    hp.mollview(m, fig=fig.number,
                title="DELVE DR3 Gold truth-anchored g depth (nside 512)",
                unit="mag", min=lo, max=hi, cmap="viridis")
    hp.graticule(dpar=30, dmer=60, color="0.7", lw=0.4)
    fig.savefig(OUT / "DELVE_depth_map.png", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {OUT/'DELVE_depth_map.png'}")


def main():
    print(f"building DELVE survey doc figures from {DATA}")
    fig_efficiencies()
    fig_errors()
    fig_depth()
    print(f"done -> {OUT}")


if __name__ == "__main__":
    main()
