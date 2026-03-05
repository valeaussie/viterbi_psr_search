#!/usr/bin/env python

import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy.signal as sig
import scipy.fftpack as fft
import argparse
from numpy.random import exponential
from viterbi import viterbi,backtrace
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib

#to use LaTeX for text rendering in figures
matplotlib.rcParams["text.latex.preamble"] += (
    r'\usepackage[dvips]{graphicx}'
    r'\usepackage{amsmath}'
    r'\usepackage{amssymb}'
)
#figure size and style
plt.rcParams["figure.figsize"] = [10,8]
rc('text', usetex=True)
rc('font', size=16.0)
rc('font',**{'family':'serif'})
plt.rcParams['figure.facecolor'] = 'white'

def make_fake_spec(NT, Nbin):
    """
    Generate a fake spectrogram (emission matrix) for pure noise.

    Parameters
    ----------
    NT : int
        Number of time steps (coherent segments).

    Nbin : int
        Number of frequency bins (HMM states).

    Returns
    -------
    spec : ndarray of shape (Nbin, NT)
        Matrix of exponentially distributed random variables.

    Explanation
    -----------
    Under Gaussian time-domain noise, periodogram power in each
    frequency bin follows an exponential distribution (chi-square
    with 2 degrees of freedom). This function simulates that.
    """
    return exponential(size=(Nbin, NT))

if __name__ == '__main__':
    # --------------------------------------------
    # Parse command-line arguments
    # --------------------------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument('--NT', type=int, help="Number of timesteps")
    parser.add_argument('--Nbin', type=int, help="Number of frequency bins")
    parser.add_argument('--Nreal', type=int, help="Number of realisations")
    parser.add_argument('--ll_file', type=str, default=None, help="File to save loglikes to")
    parser.add_argument('--ll_thresh', type=float, default=-1e7, help="Threshold above which to save loglikes")  
    args = parser.parse_args()

    # --------------------------------------------
    # Allocate storage for results
    # --------------------------------------------
    # loglikes will store final Viterbi scores for:
    #   each frequency bin (rows)
    #   each realisation (columns)
    loglikes = np.zeros((args.Nbin, args.Nreal))

    # --------------------------------------------
    # Monte Carlo loop
    # --------------------------------------------
    for i in range(args.Nreal):
        print(f"Doing realisation {i+1} of {args.Nreal}")

        # Generate synthetic noise-only spectrogram
        # Shape = (frequency bin, time step)
        spec = make_fake_spec(args.NT, args.Nbin)

        # Run Viterbi dynamic programming
        # delta[f,t] = best cumulative score ending at frequency f at time t
        delta, backptrs = viterbi(spec)

        # Store the final cumulative scores at the last time step
        # delta[:, -1] gives best score for each possible ending frequency
        loglikes[:, i] =  delta[:, -1]

    # --------------------------------------------
    # Flatten results across bins and realisations
    # --------------------------------------------
    # This gives us a single array of log-likelihoods (Viterbi scores) for all bins and realisations.
    loglikes = loglikes.flatten()
    print(np.mean(loglikes))
    print(np.std(loglikes))

    # --------------------------------------------
    # Optionally save high-score values
    # --------------------------------------------
    if args.ll_file:
        np.savetxt(args.ll_file, loglikes[loglikes > args.ll_thresh])
