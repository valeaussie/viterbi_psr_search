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

matplotlib.rcParams["text.latex.preamble"] += r'\usepackage[dvips]{graphicx}\usepackage{amsmath}\usepackage{amssymb}'
plt.rcParams["figure.figsize"] = [10,8]
rc('text', usetex=True)
rc('font', size=16.0)
rc('font',**{'family':'serif'})
plt.rcParams['figure.facecolor'] = 'white'

def make_fake_spec(NT, Nbin):
    return exponential(size=(Nbin, NT))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--NT', type=int, help="Number of timesteps")
    parser.add_argument('--Nbin', type=int, help="Number of frequency bins")
    parser.add_argument('--Nreal', type=int, help="Number of realisations")
    parser.add_argument('--ll_file', type=str, default=None, help="File to save loglikes to")
    parser.add_argument('--ll_thresh', type=float, default=-1e7, help="Threshold above which to save loglikes")
    
    args = parser.parse_args()

    loglikes = np.zeros((args.Nbin, args.Nreal))

    for i in range(args.Nreal):
        print(f"Doing realisation {i+1} of {args.Nreal}")
        spec = make_fake_spec(args.NT, args.Nbin)

        delta, backptrs = viterbi(spec)

        loglikes[:, i] =  delta[:, -1]

    loglikes = loglikes.flatten()
    print(np.mean(loglikes))
    print(np.std(loglikes)) 
    if args.ll_file:
        np.savetxt(args.ll_file, loglikes[loglikes > args.ll_thresh])
