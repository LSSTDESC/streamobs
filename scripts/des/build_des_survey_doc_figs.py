"""
build_des_survey_doc_figs.py
----------------------------
Generate the DES-SPECIFIC figures for docs/source/surveys/DES.md.

The three standard survey figures -- efficiency, photo-error and
footprint/depth -- are built for every release by
scripts/build_survey_doc_figs.py. Only the two figures with no counterpart in
any other survey live here:

  1. DES_surrogate_confusion.png - the EXT_XGB surrogate's confusion against the
                                   real catalogue, per magnitude bin
  2. DES_validation.png          - SPLASH-SXDF stellar completeness vs magnitude
                                   against the Bechtol et al. Table A.3
                                   benchmarks (skipped if the CSVs are absent)

Both go to docs/source/_static/des_yr6/

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
    return np.genfromtxt(
        path,
        delimiter=",",
        names=header.split(","),
        skip_header=sum(1 for _ in open(path) if _.startswith("#")),
    )


def fig_validation():
    f1 = ART / "des_y6_splash_validation_ext1.csv"
    f0 = ART / "des_y6_splash_validation_ext0.csv"
    if not (f1.exists() and f0.exists()):
        print(
            "  (no SPLASH validation CSVs; run "
            "scripts/des/validate_des_xgb_with_splash.py first)"
        )
        return
    import pandas as pd

    d1, d0 = pd.read_csv(f1), pd.read_csv(f0)
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for d, lab, col in (
        (d1, r"$0 \leq$ EXT_XGB $\leq 1$ (complete)", "#1f77b4"),
        (d0, r"EXT_XGB $=0$ (pure)", "#d62728"),
    ):
        g = d[d["n_star"] > 20]
        ax.plot(g["mag"], g["completeness"], "o-", color=col, lw=2, ms=4, label=lab)

    # Bechtol et al. Table A.3 integrated benchmarks
    for lo, hi, val, col in (
        (17.5, 22.5, 0.980, "#1f77b4"),
        (16.5, 23.5, 0.943, "#1f77b4"),
        (17.5, 22.5, 0.921, "#d62728"),
        (16.5, 23.5, 0.793, "#d62728"),
    ):
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
    ax.plot(
        g["mag_g"],
        g["a"],
        "o-",
        color="#1f77b4",
        lw=2,
        ms=3.5,
        label=r"$a = P(S{=}1 \,|\, {\rm EXT\_XGB} \leq 1)$",
    )
    ax.plot(
        g["mag_g"],
        g["b"],
        "o-",
        color="#d62728",
        lw=2,
        ms=3.5,
        label=r"$b = P(S{=}1 \,|\, {\rm EXT\_XGB} > 1)$",
    )
    ax.set_xlabel(r"BDF_MAG_$g$")
    ax.set_ylabel("probability")
    ax.set_ylim(-0.02, 1.05)
    ax.set_title(
        "EXT_XGB surrogate confusion, the deconvolution input\n"
        r"eff$_X$ = (eff$_S$ $-$ b) / (a $-$ b)",
        fontsize=10,
    )
    ax.legend(loc="center left", fontsize=9)
    ax.text(
        0.97,
        0.42,
        "recall falls faintward: without the deconvolution\n"
        "classification_eff would read ~17% low by $g\\approx25$",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        color="0.35",
    )
    fig.tight_layout()
    fig.savefig(OUT / "DES_surrogate_confusion.png")
    plt.close(fig)
    print(f"  wrote {OUT/'DES_surrogate_confusion.png'}")


if __name__ == "__main__":
    print("building DES Y6 doc figures:")
    fig_surrogate()
    fig_validation()
    print(f"done -> {OUT}")
