#!/usr/bin/env python
"""
More modular stream generation example.
"""
import os
import copy
from importlib import reload

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate


import streamobs.model; reload(streamobs.model)
from streamobs.model import StreamModel, BackgroundModel
import streamobs.plotting; reload(streamobs.plotting)
from streamobs.plotting import plot_stream
from streamobs.utils import parse_config

def generate_stream(config):
    """ Generate the simulated stream.

    Parameters
    ----------
    config : configuration file

    Return
    ------
    stars_df : output data frame
    """

    print("Generating stream...")
    stream = StreamModel(config['stream'])
    stream_df = stream.sample(config['stream']['nstars'] )
    print(f"  generated {len(stream_df)} stream stars.")

    print("Generating background...")
    bkg = BackgroundModel(config['background'])
    bkg_df = bkg.sample(config['background']['nstars'])
    print(f"  generated {len(bkg_df)} background stars.")

    print("Combining stream and background.")
    out = pd.concat([stream_df, bkg_df])
    out['flag'] = np.hstack([np.ones(len(stream_df),dtype=int),
                            np.zeros(len(bkg_df),dtype=int)])

    return out


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config',
                        help='configuration file')
    parser.add_argument('-o','--outfile', default='stream_out.csv',
                        help='output file')
    parser.add_argument('-p','--plot', action='store_true',
                        help='plot stream and background')
    args = parser.parse_args()

    print(f"Reading config: {args.config}")
    config = parse_config(args.config)

    # Generate the stream and background
    stars = generate_stream(config)

    print(f"Writing {args.outfile}")
    stars.to_csv(args.outfile, index=False)

    if args.plot:
        print(f"Plotting stars...")
        fig,ax = plot_stream(stars['phi1'], stars['phi2'])
        pngfile = args.outfile.replace('.csv','.png')
        print(f"Writing {pngfile}")
        plt.savefig(pngfile, facecolor="white")
