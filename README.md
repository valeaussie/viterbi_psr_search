# Viterbi HMM Pulsar Search

## Overview

This repository implements a semi-coherent pulsar search algorithm designed to detect radio pulsars in compact binary systems.  

The method combines:

- A **Hidden Markov Model (HMM)** to model time-varying pulse frequency  
- The **Viterbi algorithm** to recover the most likely frequency evolution  
- A **matched filtering approach** (Schuster periodogram) for per-segment detection  

This approach is designed to handle strong Doppler modulation caused by short orbital period binaries.

---

## Method Summary

1. The input time series is divided into short segments.
2. A periodogram (matched filter) is computed for each segment.
3. An HMM models how the pulsar frequency evolves between segments.
4. The Viterbi algorithm finds the most probable frequency track.
5. Detection statistics are computed from the recovered path.

This allows semi-coherent (coherent in chunks, incoherently combined between chunks) detection of compact binary pulsars under realistic survey conditions.

---

## Repository Structure

### Core Code

- `code/viterbi.py`  
  Core HMM and Viterbi implementation.

- `code/process_timeseries.py`  
  Preprocessing and segmentation of time series data.

- `code/process_timeseries_het.py`  
  Heterodyned time series processing.

- `code/do_noise.py`  
  Noise generation and handling (for null test).

---

### Simulation & Experiments

- `generate_data.py` (multiple directories)  
  Scripts to generate simulated pulsar signals.

- `scan_*.py`  
  Search pipelines for parameter scans.

- `Porb/`  
  Orbital period experiments.

- `flux/`  
  Flux sensitivity studies.

- `harmonics/`  
  Harmonic analysis experiments.

- `scalloping/`  
  Scalloping loss investigations.

- `sample_datasets/`  
  Example datasets for testing.

---

## Environment

Recommended:

You will have your environment
source /fred/oz022/$USER/software/envrmnts/setup_viterbi_psr.sh

```bash
pip install numpy scipy matplotlib

## Run code - From code/experiments/your_experiment
Running a Search

From inside your experiment directory:

python ../../process_timeseries_het.py --params search.params


You may override any parameter from the command line:

python ../../process_timeseries_het.py --params search.params --Nsft 32


Command-line arguments override values in search.params.

Required Parameters

These must be provided in search.params (or via CLI):

infile — Path to float32 .dat time series

tsamp — Sampling interval (seconds)

Tsft — Length of each coherent segment (seconds)

f0 — Bottom of search frequency band (Hz)

bw — Width of search frequency band (Hz)

Optional Parameters

Nsft (default: use full dataset)

padding_factor (default: 1)

num_harm (default: 1)

top_paths (default: 1)

out_prefix (default: "search")

plot_path (True/False)

save_delta (True/False)

spec_flo, spec_fhi (plot limits)

Example search.params
infile = "data/yourfile.dat"
tsamp = 6.4e-5
Tsft = 440
Nsft = 16
f0 = 150
bw = 50
padding_factor = 1
out_prefix = "exp_01_first_test"
plot_path = True
num_harm = 1


