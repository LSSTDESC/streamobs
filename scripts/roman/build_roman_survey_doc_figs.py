"""
build_roman_survey_doc_figs.py
-------------------------------
Generate the two survey-page figures for docs/source/surveys/Roman.md:

  1. Roman_efficiencies.png  – detection_eff, classification_eff,
                               classification_detection_eff (stellar efficiency CSV)
                               + missclassification_eff (galaxy misclass CSV)
  2. Roman_errors.png        – 10^log_mag_err vs delta_mag for the sample and
                               catalog curves (log y-axis)

Both go to  docs/source/_static/roman_dc2/

Usage
-----
  conda activate streamobs
  python scripts/roman/build_roman_survey_doc_figs.py
"""

import pathlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO = pathlib.Path(__file__).resolve().parents[2]
DATA = REPO / "data/surveys/roman_dc2"
OUT  = REPO / "docs/source/_static/roman_dc2"
OUT.mkdir(parents=True, exist_ok=True)

EFF_CSV      = DATA / "roman_stellar_efficiency_cutf158.csv"
MISCLASS_CSV = DATA / "roman_galaxy_misclass_cutf158.csv"
ERR_SAMP_CSV = DATA / "roman_photoerror_f158.csv"
ERR_CAT_CSV  = DATA / "roman_photoerror_f158_catalog.csv"

# ---------------------------------------------------------------------------
# Helper: read a multi-comment-header CSV (last comment line = column names)
# ---------------------------------------------------------------------------
def read_csv(path):
    """Read a CSV where '#'-prefixed lines precede data; last '#' line = header."""
    path = pathlib.Path(path)
    header_line = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                header_line = line.lstrip("#").strip()
            else:
                break
    if header_line is None:
        raise ValueError(f"No comment header found in {path}")
    cols = [c.strip() for c in header_line.split(",")]
    data = np.genfromtxt(path, delimiter=",", comments="#")
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return {c: data[:, i] for i, c in enumerate(cols)}


# ---------------------------------------------------------------------------
# Style matching existing LSST figures: clean, seaborn-ish
# ---------------------------------------------------------------------------
mpl.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 12,
    "axes.labelsize": 12,
    "legend.fontsize": 10,
    "figure.dpi": 150,
})

# ---------------------------------------------------------------------------
# Figure 1 – efficiencies
# ---------------------------------------------------------------------------
eff  = read_csv(EFF_CSV)
mis  = read_csv(MISCLASS_CSV)

fig, ax = plt.subplots(figsize=(7, 5))

ax.plot(eff["delta_mag"], eff["detection_eff"],
        label="Detection efficiency", lw=2, color="#1f77b4")
ax.plot(eff["delta_mag"], eff["classification_eff"],
        label="Classification efficiency", lw=2, color="#ff7f0e")
ax.plot(eff["delta_mag"], eff["classification_detection_eff"],
        label="Combined efficiency\n(detection × classification)", lw=2,
        color="#2ca02c")
ax.plot(mis["delta_mag"], mis["missclassification_eff"],
        label="Galaxy misclassification efficiency\n(compact galaxies)", lw=2,
        ls="--", color="#d62728")

ax.set_xlabel(r"$\Delta m = m - m_{\rm lim}$ [mag]")
ax.set_ylabel("Efficiency")
ax.set_xlim(left=eff["delta_mag"].min() - 0.1)
ax.set_ylim(-0.03, 1.08)
ax.axhline(0.5, color="k", lw=0.7, ls=":", alpha=0.5)
ax.axvline(0.0, color="k", lw=0.7, ls=":", alpha=0.5, label="Magnitude limit")
ax.legend(loc="upper right", framealpha=0.85)
ax.set_title("Roman F158 stellar efficiencies and galaxy contamination")
fig.tight_layout()

out_eff = OUT / "Roman_efficiencies.png"
fig.savefig(out_eff, dpi=150)
plt.close(fig)
print(f"Written: {out_eff}  ({out_eff.stat().st_size // 1024} kB)")


# ---------------------------------------------------------------------------
# Figure 2 – photometric errors (log y-axis)
# ---------------------------------------------------------------------------
samp = read_csv(ERR_SAMP_CSV)
cat  = read_csv(ERR_CAT_CSV)

fig, ax = plt.subplots(figsize=(7, 5))

ax.semilogy(samp["delta_mag"], 10 ** samp["log_mag_err"],
            label="Sample curve (truth-based scatter, drives noise draw)",
            lw=2, color="#1f77b4")
ax.semilogy(cat["delta_mag"], 10 ** cat["log_mag_err"],
            label="Catalog curve (reported magerr, drives S/N cut)",
            lw=2, color="#ff7f0e", ls="--")

ax.set_xlabel(r"$\Delta m = m - m_{\rm lim}$ [mag]")
ax.set_ylabel(r"Photometric uncertainty $\sigma_m$ [mag]")
ax.axvline(0.0, color="k", lw=0.7, ls=":", alpha=0.5, label="Magnitude limit")
ax.legend(loc="upper left", framealpha=0.85)
ax.set_title("Roman F158 photometric error model")
fig.tight_layout()

out_err = OUT / "Roman_errors.png"
fig.savefig(out_err, dpi=150)
plt.close(fig)
print(f"Written: {out_err}  ({out_err.stat().st_size // 1024} kB)")


# ---------------------------------------------------------------------------
# Figure 3 – HLWAS "all" all-sky depth maps (F158 + F106)
# ---------------------------------------------------------------------------
import healpy as hp

HLWAS_ALL_DIR = REPO / "data/surveys/roman_hlwas_all"

# Band configuration: (fits_stem, label, panel_title)
BANDS = [
    ("roman_hlwas_all_maglim_f158_nside1024.fits.gz", "F158", "HLWAS all  F158"),
    ("roman_hlwas_all_maglim_f106_nside1024.fits.gz", "F106", "HLWAS all  F106"),
]

# Collect maps and derive shared colour scale
maps_allsky = {}
for fname, band, title in BANDS:
    path = HLWAS_ALL_DIR / fname
    if path.exists():
        m = hp.read_map(str(path), nest=False, verbose=False)
        valid = m > hp.UNSEEN + 1
        maps_allsky[band] = {"arr": m, "valid": valid, "title": title, "fname": fname}

n_panels = len(maps_allsky)

# Shared colour limits: 1st–99th percentile of all covered pixels, rounded to 0.05
all_valid_vals = np.concatenate([v["arr"][v["valid"]] for v in maps_allsky.values()])
vmin_raw = np.percentile(all_valid_vals, 1)
vmax_raw = np.percentile(all_valid_vals, 99)
# Round to nearest 0.05 for clean colorbar ticks
vmin = round(float(vmin_raw) / 0.05) * 0.05
vmax = round(float(vmax_raw) / 0.05) * 0.05

# healpy mollview draws into its own figure; we use sub= to layout panels
fig_allsky = plt.figure(figsize=(12 * n_panels, 6))

for i, (band, info) in enumerate(maps_allsky.items()):
    m_plot = info["arr"].copy()
    m_plot[~info["valid"]] = hp.UNSEEN
    med = np.median(info["arr"][info["valid"]])
    title_str = f"{info['title']}\nmedian = {med:.3f} AB"
    hp.mollview(
        m_plot,
        title=title_str,
        min=vmin,
        max=vmax,
        cmap="viridis",
        unit="5σ maglim [AB]",
        badcolor="white",
        bgcolor="white",
        fig=fig_allsky.number,
        sub=(1, n_panels, i + 1),
        notext=True,
    )

fig_allsky.suptitle(
    "Roman HLWAS 'all' tier — exposure-time-scaled depth maps (nside=1024)",
    fontsize=14,
    y=1.01,
)
fig_allsky.tight_layout()

out_allsky = OUT / "Roman_depth_allsky.png"
fig_allsky.savefig(out_allsky, dpi=150, bbox_inches="tight")
plt.close(fig_allsky)
print(f"Written: {out_allsky}  ({out_allsky.stat().st_size // 1024} kB)")
