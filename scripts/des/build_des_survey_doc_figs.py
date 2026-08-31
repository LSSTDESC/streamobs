"""
build_des_survey_doc_figs.py
----------------------------
Generate the figures for docs/source/surveys/DES.md from the shipped des/yr6
products. Mirrors scripts/roman/build_roman_survey_doc_figs.py.

  1. DES_efficiencies.png  - detection_eff, classification_eff,
                             classification_detection_eff + the galaxy
                             missclassification_eff on the same axes
  2. DES_errors.png        - the two photo-error curves (sample vs catalog),
                             with the S/N=5 line and the measured inflation
  3. DES_depth.png         - griz truth-anchored depth histograms, with the
                             input-map medians marked so the anchor shift is
                             visible
  4. DES_depth_map.png     - the g-band depth map on the sky
  5. DES_validation.png    - SPLASH-SXDF stellar completeness vs magnitude
                             against the Bechtol et al. Table A.3 benchmarks
                             (skipped if the validation CSVs are absent)

All go to docs/source/_static/des_yr6/

Usage
-----
  python scripts/des/build_des_survey_doc_figs.py
"""

import pathlib

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parents[2]
DATA = REPO / "data/surveys/des_yr6"
ART = REPO / "artifacts/des_y6"
OUT = REPO / "docs/source/_static/des_yr6"
OUT.mkdir(parents=True, exist_ok=True)

BANDS = ("g", "r", "i", "z")
# input healsparse map medians, for the anchor-shift annotation
INPUT_MEDIAN = {"g": 25.321, "r": 25.130, "i": 24.582, "z": 23.878}
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
    eff = read_csv(DATA / "des_yr6_stellar_efficiency_cutg.csv")
    mis = read_csv(DATA / "des_yr6_galaxy_misclass_cutg.csv")

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(eff["delta_mag"], eff["detection_eff"], lw=2, color="#1f77b4",
            label="detection")
    ax.plot(eff["delta_mag"], eff["classification_eff"], lw=2, color="#d62728",
            label=r"classification $|$ detected")
    ax.plot(eff["delta_mag"], eff["classification_detection_eff"], lw=2.4,
            color="k", label="detection $\\times$ classification")
    ax.plot(mis["delta_mag"], mis["missclassification_eff"], lw=1.8,
            color="#2ca02c", ls="--",
            label="galaxy misclassified as star")

    # the bright end is clamped, not measured -- say so on the figure
    clamp_hi = -6.65
    ax.axvspan(eff["delta_mag"].min(), clamp_hi, color="0.85", alpha=0.5, zorder=0)
    ax.annotate("classification clamped:\n<100 distinct parent sources",
                xy=(-7.5, 0.62), xytext=(-6.3, 0.42), fontsize=7.5, color="0.3",
                ha="left", arrowprops=dict(arrowstyle="->", color="0.5", lw=0.8))

    ax.axvline(0.0, color="0.5", lw=0.9, ls=":")
    ax.text(0.06, 0.30, "maglim", fontsize=8, color="0.45", rotation=90)
    ax.set_xlabel(r"$\Delta$mag = mag$_g^{\rm true}$ $-$ maglim$_g$")
    ax.set_ylabel("efficiency")
    ax.set_ylim(-0.02, 1.12)
    ax.set_xlim(eff["delta_mag"].min(), 2.0)
    ax.set_title("DES Y6 stellar efficiency and galaxy misclassification "
                 r"($0 \leq$ EXT_XGB $\leq 1$)", fontsize=10)
    ax.legend(loc="lower left", fontsize=8.5, ncol=2,
              bbox_to_anchor=(0.0, -0.02))
    fig.tight_layout()
    fig.savefig(OUT / "DES_efficiencies.png")
    plt.close(fig)
    print(f"  wrote {OUT/'DES_efficiencies.png'}")


def fig_errors():
    s = read_csv(DATA / "des_yr6_photoerror_g.csv")
    c = read_csv(DATA / "des_yr6_photoerror_g_catalog.csv")
    sig_s, sig_c = 10 ** s["log_mag_err"], 10 ** c["log_mag_err"]

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7.2, 5.6), sharex=True,
                                  gridspec_kw={"height_ratios": [2.4, 1]})
    ax.semilogy(s["delta_mag"], sig_s, lw=2, color="#d62728",
                label="sample (truth scatter) $\\rightarrow$ noise draw")
    ax.semilogy(c["delta_mag"], sig_c, lw=2, color="#1f77b4",
                label="catalog (reported magerr) $\\rightarrow$ S/N cut")
    ax.axhline(SIG_SN5, color="0.5", lw=0.9, ls=":")
    ax.text(0.5, SIG_SN5 * 1.25, r"$\sigma = 0.217$ (S/N $=5$)",
            fontsize=8, color="0.4")
    ax.axvline(0.0, color="0.5", lw=0.9, ls=":")
    ax.set_ylabel(r"$\sigma_{\rm mag}$")
    ax.set_ylim(3e-4, 1.2)
    ax.set_title("DES Y6 two-curve photometric error model (g)", fontsize=10)
    ax.legend(loc="lower right", fontsize=8.5)

    n = min(sig_s.size, sig_c.size)
    ax2.plot(s["delta_mag"][:n], sig_s[:n] / sig_c[:n], lw=2, color="k")
    ax2.axhline(1.0, color="0.5", lw=0.9, ls=":")
    ax2.set_ylabel("sample / catalog")
    ax2.set_xlabel(r"$\Delta$mag = mag$_g^{\rm true}$ $-$ maglim$_g$")
    ax2.text(0.98, 0.92,
             "inflation is strongly magnitude-dependent:\n"
             r"$\sim$1.5 near the limit, $\sim$35 at the bright end, where the"
             "\nreported error is statistical only and misses the ~0.02 mag\n"
             "systematic floor the truth scatter actually sees",
             transform=ax2.transAxes, ha="right", va="top", fontsize=7.5,
             color="0.35")
    fig.tight_layout()
    fig.savefig(OUT / "DES_errors.png")
    plt.close(fig)
    print(f"  wrote {OUT/'DES_errors.png'}")


def fig_depth():
    import healpy as hp

    fig, axes = plt.subplots(2, 2, figsize=(8.2, 5.4))
    for ax, b in zip(axes.ravel(), BANDS):
        f = DATA / f"des_yr6_maglim_{b}_nside1024.fits.gz"
        m = hp.read_map(str(f), verbose=False) if "verbose" in \
            hp.read_map.__code__.co_varnames else hp.read_map(str(f))
        good = np.isfinite(m) & (m > 0) & (m != hp.UNSEEN)
        v = m[good]
        med = np.median(v)
        ax.hist(v, bins=80, color="#1f77b4", alpha=0.85)
        ax.axvline(med, color="k", lw=1.6, label=f"anchored {med:.3f}")
        ax.axvline(INPUT_MEDIAN[b], color="#d62728", lw=1.4, ls="--",
                   label=f"input map {INPUT_MEDIAN[b]:.3f}")
        ax.set_title(f"{b}  (shift {med - INPUT_MEDIAN[b]:+.3f})", fontsize=9)
        ax.set_xlabel("truth-anchored 5$\\sigma$ depth")
        ax.set_yticks([])
        ax.legend(fontsize=7.5, loc="upper left")
    fig.suptitle("DES Y6 truth-anchored depth, nside 1024 "
                 "(shifts are smooth and ordered with wavelength)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "DES_depth.png")
    plt.close(fig)
    print(f"  wrote {OUT/'DES_depth.png'}")

    m = hp.read_map(str(DATA / "des_yr6_maglim_g_nside1024.fits.gz"))
    m = np.where(np.isfinite(m) & (m > 0), m, hp.UNSEEN)
    fig = plt.figure(figsize=(7.6, 4.4))
    hp.mollview(m, fig=fig.number, title="DES Y6 truth-anchored g depth "
                "(nside 1024)", unit="mag", min=24.4, max=25.6, cmap="viridis")
    hp.graticule(dpar=30, dmer=60, color="0.7", lw=0.4)
    fig.savefig(OUT / "DES_depth_map.png", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {OUT/'DES_depth_map.png'}")


def fig_validation():
    f1 = ART / "des_y6_splash_validation_ext1.csv"
    f0 = ART / "des_y6_splash_validation_ext0.csv"
    if not (f1.exists() and f0.exists()):
        print("  (no SPLASH validation CSVs; run "
              "scripts/des/validate_des_xgb_with_splash.py first)")
        return
    import pandas as pd

    d1, d0 = pd.read_csv(f1), pd.read_csv(f0)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for d, lab, col in ((d1, r"$0 \leq$ EXT_XGB $\leq 1$ (complete)", "#1f77b4"),
                        (d0, r"EXT_XGB $=0$ (pure)", "#d62728")):
        g = d[d["n_star"] > 20]
        ax.plot(g["mag"], g["completeness"], "o-", color=col, lw=2, ms=4, label=lab)

    # Bechtol et al. Table A.3 integrated benchmarks
    for lo, hi, val, col in ((17.5, 22.5, 0.980, "#1f77b4"),
                             (16.5, 23.5, 0.943, "#1f77b4"),
                             (17.5, 22.5, 0.921, "#d62728"),
                             (16.5, 23.5, 0.793, "#d62728")):
        ax.hlines(val, lo, hi, color=col, ls=":", lw=1.6, alpha=0.8)
    ax.plot([], [], ls=":", color="0.4", label="Bechtol et al. Table A.3 (integrated)")

    ax.set_xlabel(r"DES MAG_AUTO_$i$")
    ax.set_ylabel("stellar completeness")
    ax.set_ylim(0, 1.05)
    ax.set_title("DES Y6 EXT_XGB completeness vs SPLASH-SXDF truth", fontsize=10)
    ax.legend(loc="lower left", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(OUT / "DES_validation.png")
    plt.close(fig)
    print(f"  wrote {OUT/'DES_validation.png'}")


def fig_surrogate():
    f = ART / "des_y6_xgb_confusion.csv"
    if not f.exists():
        print("  (no confusion table; run build_des_xgb_surrogate.py first)")
        return
    import pandas as pd

    c = pd.read_csv(f)
    g = c[(c["n_pos"] > 200) & (c["n_neg"] > 200)]
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(g["mag_g"], g["a"], "o-", color="#1f77b4", lw=2, ms=3.5,
            label=r"$a = P(S{=}1 \,|\, {\rm EXT\_XGB} \leq 1)$")
    ax.plot(g["mag_g"], g["b"], "o-", color="#d62728", lw=2, ms=3.5,
            label=r"$b = P(S{=}1 \,|\, {\rm EXT\_XGB} > 1)$")
    ax.set_xlabel(r"BDF_MAG_$g$")
    ax.set_ylabel("probability")
    ax.set_ylim(-0.02, 1.05)
    ax.set_title("EXT_XGB surrogate confusion, the deconvolution input\n"
                 r"eff$_X$ = (eff$_S$ $-$ b) / (a $-$ b)", fontsize=10)
    ax.legend(loc="center left", fontsize=9)
    ax.text(0.97, 0.42, "recall falls faintward: without the deconvolution\n"
            "classification_eff would read ~17% low by $g\\approx25$",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color="0.35")
    fig.tight_layout()
    fig.savefig(OUT / "DES_surrogate_confusion.png")
    plt.close(fig)
    print(f"  wrote {OUT/'DES_surrogate_confusion.png'}")


if __name__ == "__main__":
    print("building DES Y6 doc figures:")
    fig_efficiencies()
    fig_errors()
    fig_depth()
    fig_surrogate()
    fig_validation()
    print(f"done -> {OUT}")
