#!/usr/bin/env python
"""Diagnose the failing PR #67 S/N=10 assertion in test_log_photo_error_behavior.

The test does:
    error = get_photo_error(band, maglim - 0.75, maglim, kind="catalog")
    snr   = 1 / error
    assert np.isclose(snr, 10.0, atol=0.5)

It fails for EVERY band (the message names 'g' only because bands sort first).
This figure shows why: `1/magerr` is not S/N. S/N = 2.5/ln(10)/magerr, so S/N=10
means magerr = 0.10857, not 0.1000. The curve is right; the conversion is not.

Local, gitignored helper (do not commit). Run with:
  python scripts/plot_snr10_test_diagnosis.py
"""
import matplotlib.pyplot as plt
import numpy as np

import streamobs.surveys as S

K = 2.5 / np.log(10)  # 1.08574 -- magerr = K / SNR
PROBE = -0.75  # delta_mag the test evaluates at
OUT = "snr10_test_diagnosis.png"

sv = S.Survey.load(survey="lsst", release="dc2", verbose=False)
ml = float(np.median([v for v in [26.846]]))
ref = sv.completeness_band

d = np.arange(-2.0, 0.05, 0.005)
curves = {}
for b in ("r", "g"):
    sv.sys_error[b] = 0.0
    curves[b] = np.array(
        [
            float(np.atleast_1d(sv.get_photo_error(b, ml + x, ml, kind="catalog"))[0])
            for x in d
        ]
    )

e_r = float(np.atleast_1d(sv.get_photo_error("r", ml + PROBE, ml, kind="catalog"))[0])
e_g = float(np.atleast_1d(sv.get_photo_error("g", ml + PROBE, ml, kind="catalog"))[0])

fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))

# ---- panel 1: the curves and the two candidate targets --------------------
a = ax[0]
a.plot(d, curves["r"], color="C0", lw=2.2, label=f"r  (reference band -> det_ok curve)")
a.plot(
    d,
    curves["g"],
    color="C3",
    lw=2.2,
    ls="--",
    label="g  (forced photometry -> _nocut curve)",
)
a.axhline(0.1000, color="k", ls=":", lw=1.4)
a.axhline(K / 10, color="green", ls="-.", lw=1.4)
a.axvline(PROBE, color="0.5", lw=1.0)
a.annotate(
    "magerr = 0.1000\n(what the test's 1/err == 10 implies)",
    xy=(-1.95, 0.1000),
    xytext=(-1.93, 0.113),
    fontsize=8,
    color="k",
)
a.annotate(
    f"magerr = {K/10:.5f}\n(TRUE S/N = 10)",
    xy=(-1.95, K / 10),
    xytext=(-1.93, 0.083),
    fontsize=8,
    color="green",
)
a.plot([PROBE], [e_r], "o", color="C0", ms=7, zorder=5)
a.plot([PROBE], [e_g], "o", color="C3", ms=7, zorder=5)
a.set(
    xlabel=r"$\Delta m = m - {\rm maglim}$",
    ylabel=r"catalog $\sigma$ (mag)",
    xlim=(-2, 0),
    ylim=(0.02, 0.16),
    title="Both curves land on TRUE S/N = 10 at the probe",
)
a.legend(fontsize=8, loc="upper left")
a.grid(alpha=0.3)

# ---- panel 2: S/N under both conventions ---------------------------------
b_ = ax[1]
b_.axhspan(9.5, 10.5, color="green", alpha=0.12, label="test tolerance (10 $\\pm$ 0.5)")
b_.plot(
    d, K / curves["r"], color="C0", lw=2.2, label="r: TRUE S/N $=2.5/\\ln 10/\\sigma$"
)
b_.plot(d, K / curves["g"], color="C3", lw=2.2, ls="--", label="g: TRUE S/N")
b_.plot(
    d,
    1 / curves["r"],
    color="C0",
    lw=1.6,
    alpha=0.55,
    ls=":",
    label="r: what the test computes ($1/\\sigma$)",
)
b_.plot(
    d,
    1 / curves["g"],
    color="C3",
    lw=1.6,
    alpha=0.55,
    ls=":",
    label="g: what the test computes ($1/\\sigma$)",
)
b_.axhline(10, color="k", lw=1.0)
b_.axvline(PROBE, color="0.5", lw=1.0)
for e, c, nm in ((e_r, "C0", "r"), (e_g, "C3", "g")):
    b_.plot([PROBE], [K / e], "o", color=c, ms=7, zorder=5)
    b_.plot([PROBE], [1 / e], "s", color=c, ms=6, zorder=5, alpha=0.7)
b_.annotate(
    f"TRUE S/N: r={K/e_r:.2f}, g={K/e_g:.2f}  PASS",
    xy=(PROBE, K / e_r),
    xytext=(-1.95, 11.6),
    fontsize=9,
    color="green",
)
b_.annotate(
    f"1/$\\sigma$: r={1/e_r:.3f}, g={1/e_g:.3f}  FAIL",
    xy=(PROBE, 1 / e_g),
    xytext=(-1.95, 8.1),
    fontsize=9,
    color="firebrick",
)
b_.set(
    xlabel=r"$\Delta m = m - {\rm maglim}$",
    ylabel="S/N",
    xlim=(-2, 0),
    ylim=(4, 14),
    title="The 1.0857 factor is the whole failure",
)
b_.legend(fontsize=7.5, loc="lower left")
b_.grid(alpha=0.3)

fig.suptitle(
    "PR #67 S/N=10 assertion: the curves are correct, the S/N conversion is not "
    f"(det_ok vs _nocut differ by only {100*(e_g/e_r-1):.1f}% here)",
    fontsize=11,
)
fig.tight_layout()
fig.savefig(OUT, dpi=140, bbox_inches="tight")
print(f"wrote {OUT}")
print(f"  probe delta_mag = {PROBE}")
print(f"  r (det_ok): sigma={e_r:.5f}  1/sigma={1/e_r:.3f}  TRUE S/N={K/e_r:.3f}")
print(f"  g (_nocut): sigma={e_g:.5f}  1/sigma={1/e_g:.3f}  TRUE S/N={K/e_g:.3f}")
print(f"  det_ok vs _nocut difference at probe: {100*(e_g/e_r-1):.2f}%")
