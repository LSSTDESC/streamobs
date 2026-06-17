#!/usr/bin/env python
"""
Example of generating a matched filter map for pal 5
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate


def inverse_transfom_sample(vals, pdf, size):
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]
    fn = scipy.interpolate.interp1d(cdf, list(range(0, len(cdf))))
    x_new = np.random.uniform(size=np.rint(size).astype(int))
    x_new[x_new < 1e-3] = (
        1e-3  # running into error that values close to 0 are flagged as being
    )
    # below interp range
    index = np.rint(fn(x_new)).astype(int)
    return vals[index]


if __name__ == "__main__":

    config = {
        "phi1_min": -10,
        "phi1_max": 20,
        "phi2_min": -9,
        "phi2_max": 3,
        "bin_size_x": 0.1,
        "bin_size_y": 0.1,
        "maglim_g": 26.5,
        "bg": 2.5,
        "n_stream_stars": 5e3,
    }
    nbinx = int(
        np.round((config["phi1_max"] - config["phi1_min"]) / config["bin_size_x"], 0)
    )
    nbiny = int(
        np.round((config["phi2_max"] - config["phi2_min"]) / config["bin_size_y"], 0)
    )

    input_dat = pd.read_csv("./data/erkal_2016_pal_5_input.csv")
    phi2_interp = scipy.interpolate.interp1d(
        input_dat["phi1"].values, input_dat["phi2"].values
    )
    width_interp = scipy.interpolate.interp1d(
        input_dat["phi1"].values, input_dat["width"].values
    )
    # create stream
    stream = pd.DataFrame(
        {
            "phi1": inverse_transfom_sample(
                input_dat["phi1"].values,
                input_dat["nstars"].values,
                config["n_stream_stars"],
            )
        }
    )
    stream["phi2"] = np.random.normal(
        phi2_interp(stream["phi1"]), width_interp(stream["phi1"])
    )
    # create matched filter map
    streamMap, xe, ye = np.histogram2d(
        stream["phi1"],
        stream["phi2"],
        bins=[nbinx, nbiny],
        normed=False,
        range=[
            [config["phi1_min"], config["phi1_max"]],
            [config["phi2_min"], config["phi2_max"]],
        ],
    )
    bgMap = np.ones_like(streamMap) * config["bg"]
    mfMap = streamMap + np.random.poisson(bgMap)

    # create figure
    fig = plt.figure(dpi=175, figsize=(10, 2.5))
    gs = fig.add_gridspec(1, 1, wspace=0, hspace=0)

    ax = fig.add_subplot(gs[0, 0])
    ax.set_title("Mock Pal5 (fit from erkal 2016)")
    ax.set_xticks(np.arange(-20, 20, 5))
    ax.set_xlabel("$\phi_1$ [deg]")
    ax.set_ylabel("$\phi_2$ [deg]")

    bit = ax.imshow(
        (mfMap).T,
        extent=[xe.min(), xe.max(), ye.min(), ye.max()],
        aspect="auto",
        vmax=np.percentile(mfMap, 99.9),
        origin="lower",
        cmap="gray_r",
    )
    plt.savefig("./example_mock_pal5_matched_filter_map.png", facecolor="white")
