# Viterbi HMM Pulsar Search

## Overview

This repository implements a semi-coherent pulsar search algorithm designed to detect radio pulsars in compact binary systems.  

The method combines:

- A **Hidden Markov Model (HMM)** to model time-varying pulse frequency  
- The **Viterbi algorithm** to recover the most likely frequency evolution  
- A **matched filtering approach** (Schuster periodogram) for per-segment detection  

This approach is designed to handle Doppler modulations caused by short orbital period binaries.

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

- `code/process_timeseries_het.py`  
  Heterodyned time series processing.

- `code/do_noise.py`  
  Noise generation and handling (for null test).

## Environment

Recommended:

You will have your environment
source setup_your_env.sh

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

```

---

## End-to-End Search Pipeline (Real Data)

This section describes how to run a full search on a real observation, starting from a `.fits` filterbank file and ending with a folded candidate plot. Replace `path/to/...` with the actual paths on your system.

The pipeline has four stages:

1. **Convert** SIGPROC/PSRFITS `.fits` to filterbank `.fil` (if needed).
2. **Dedisperse** the filterbank at the trial DM(s) to produce a dedispersed time series.
3. **Search** the dedispersed time series with the Viterbi HMM tracker.
4. **Fold** the original `.fil` using the spin parameters (`f0`, `f1`, `f2`) recovered by the search.

External tools required (installed and on `PATH`):

- `digifil` (from DSPSR)
- `dedisperse_all_fil` (from `sigpyproc` / your dedispersion package)
- `psrfold_fil` (from PulsarX)
- Python with `numpy`, `scipy`, `matplotlib` (for the Viterbi step)

---

### Step 1: Convert `.fits` to `.fil`

```bash
digifil -threads 8 -b 8 -o path/to/output.fil path/to/input.fits
```

Flags:

- `-threads 8` parallel threads.
- `-b 8` output 8-bit samples.
- `-o` output filename.

**Note:** Record the start MJD of the observation from the `.fits` header (e.g. via `psrfits_dump` or `fitshead`). You will need it later to convert the search reference time to a `pepoch` for folding.

---

### Step 2: Dedisperse

```bash
dedisperse_all_fil \
    --dms 24.356441 \
    --ddm 0.1 \
    --ndm 1 \
    --rfi kadaneF 8 4 zdot \
    --fillPatch rand \
    --format presto \
    -v \
    -f path/to/output.fil
```

Flags:

- `--dms 24.356441` central trial DM in pc cm$^{-3}$.
- `--ddm 0.1` step size between trial DMs.
- `--ndm 1` number of trial DMs (set `--ndm 1` for a single DM, or larger for a small DM grid).
- `--rfi kadaneF 8 4 zdot` RFI mitigation chain.
- `--fillPatch rand` fill masked samples with random noise.
- `--format presto` write PRESTO-format output (`.dat` plus `.inf` sidecar).

This produces, for each trial DM, a 32-bit float dedispersed time series `output_DM<value>.dat` together with an `.inf` file that records the sampling interval `tsamp`, the start MJD, and the number of samples. Open the `.inf` file and note the value of `Width of each time series bin (sec)`. This is the `tsamp` you must pass to the Viterbi search.

---

### Step 3: Run the Viterbi Search

The Viterbi script reads a raw 32-bit float time series. The PRESTO `.dat` written by `dedisperse_all_fil` is exactly that, so it can be passed directly as `infile`.

Set up an experiment directory and a parameter file:

```
code/experiments/example_run/
    search.params
```

Example `search.params`:

```ini
infile = "path/to/output_DM24.356441.dat"
tsamp = 6.4e-5
Tsft = 30
Nsft = 64
f0 = 370
bw = 20
padding_factor = 1
num_harm = 1
top_paths = 1
out_prefix = "example_run"
plot_path = True
```

Notes on parameter choices:

- `tsamp` must match the `Width of each time series bin (sec)` in the `.inf` file produced by Step 2.
- `Tsft` is the coherent segment length in seconds. It must be short enough that the signal frequency drift within one segment is much less than $1/T_\mathrm{sft}$. For a binary with orbital period $P_\mathrm{b}$, a common choice is $T_\mathrm{sft} \lesssim P_\mathrm{b}/30$.
- `Nsft` $\times$ `Tsft` must not exceed the duration of the dedispersed time series. Set `Nsft = -1` (or omit it) to use the full dataset.
- `f0` and `bw` define the search band. `f0` is the **lower** edge in Hz, and the searched band is $[f_0,\ f_0 + \mathrm{bw}]$.
- `num_harm` controls incoherent harmonic summing. Set to $1$ for the fundamental only.

Run the search from the experiment directory:

```bash
cd code/experiments/example_run
python ../../process_timeseries_het.py --params search.params
```

CLI flags override values in `search.params`, e.g.:

```bash
python ../../process_timeseries_het.py --params search.params --Nsft 32 --num-harm 2
```

Outputs (written to the current directory):

- `example_run_spectrogram_4_harmonics.png` heteroscedastic spectrogram with the recovered Viterbi path overlaid (if `plot_path = True`).
- `example_run_loglikes_4_harmonics.png` log-likelihood vs frequency at the final time step.
- `example_run_paths.dat` top-N paths and their scores.
- `example_run_track.dat` columns: `time_s`, `frequency_hz` (Viterbi path), `fitted_frequency_hz` (sinusoidal fit).
- `example_run_sinusoid_fit.dat` the recovered spin parameters at the reference time `tref`. Key fields:
  - `tref_s` reference time in seconds **from the start of the dedispersed time series**.
  - `f0_at_tref_hz` spin frequency at `tref`.
  - `f1_at_tref_hz_per_s` first frequency derivative at `tref`.
  - `f2_at_tref_hz_per_s2` second frequency derivative at `tref`.

The terminal will also print these values, e.g.:

```
Recovered f0 at tref=945.000 s: 378.327569311425 Hz
Recovered f1: -1.624881583924e-06 Hz/s
Recovered f2: -2.340076538380e-09 Hz/s^2
```

---

### Step 4: Convert `tref` to `pepoch` (MJD)

`psrfold_fil` requires `--pepoch` as a Modified Julian Date, but the search returns `tref` in seconds since the start of the dedispersed time series. Convert with:

$$
\mathrm{pepoch} \;=\; \mathrm{MJD}_\mathrm{start} \;+\; \frac{t_\mathrm{ref}}{86400}
$$

where $\mathrm{MJD}_\mathrm{start}$ is the start MJD recorded in the `.inf` file produced by Step 2 (field `Epoch of observation (MJD)`), and $86400$ is the number of seconds in a day. For example, with $\mathrm{MJD}_\mathrm{start} = 59391.19581$ and $t_\mathrm{ref} = 945.0$ s:

$$
\mathrm{pepoch} = 59391.19581 + \frac{945.0}{86400} = 59391.20675
$$

> **Important:** the start MJD in the `.inf` file refers to the start of the **dedispersed** time series after any cropping or padding applied by `dedisperse_all_fil`. If you instead read the start MJD from the original `.fits` header, make sure it corresponds to the same time origin used by the Viterbi search. If in doubt, use the `.inf` value.

---

### Step 5: Fold

Fold the original `.fil` (not the dedispersed `.dat`) at the recovered parameters:

```bash
psrfold_fil \
    -v \
    -t 2 \
    --f0 378.327569311425 \
    --f1 -1.624881583924e-06 \
    --f2 -2.340076538380e-09 \
    --pepoch 59391.20675 \
    --dm 24.356441 \
    --clfd 2.0 \
    --rfi kadaneF 8 4 zdot \
    --fillPatch rand \
    -n 256 \
    -b 64 \
    -L 30 \
    --template path/to/PulsarX/include/template/meerkat_fold.template \
    -f path/to/output.fil
```

Flags:

- `--f0 --f1 --f2` taken from `example_run_sinusoid_fit.dat`.
- `--pepoch` MJD computed in Step 4.
- `--dm` **must equal the DM used for dedispersion in Step 2** (here, `24.356441`).
- `--clfd 2.0` time-domain RFI cleaning threshold.
- `--rfi kadaneF 8 4 zdot` matched RFI chain to the search.
- `-n 256` number of phase bins.
- `-b 64` number of frequency channels in the output archive.
- `-L 30` sub-integration length in seconds.
- `--template` PulsarX standard fold template (path system-dependent).
- `-f` input filterbank.


If you are on OzStar 

---

### Pitfalls and Sanity Checks

- **`tsamp` mismatch.** If the `tsamp` in `search.params` does not match the `.inf` file, all recovered frequencies will be scaled incorrectly. Check the `.inf` file every time.
- **DM consistency.** The fold DM in Step 5 must match the dedispersion DM in Step 2. Mixing values from different runs (as is easy to do when copying commands) will produce a smeared profile or a non-detection.
- **Search band.** Make sure the true spin frequency lies inside $[f_0,\ f_0 + \mathrm{bw}]$. If folding gives no profile, widen the band and re-run, or inspect `example_run_loglikes_4_harmonics.png` to see where the peak likelihood actually falls.
- **Pepoch convention.** A wrong `--pepoch` will rotate the apparent profile in phase but, more importantly, will mis-apply `f1` and `f2`, smearing the profile across sub-integrations even when `f0` is correct. Verify by checking that the folded time-vs-phase plot is straight, not curved.
- **Sinusoidal fit limitations.** The script fits a single sinusoid to the recovered Viterbi path with bounded amplitude (≤ $1$ Hz) and angular frequency (≤ $1$ rad/s). For wide-orbit binaries or highly non-sinusoidal drifts, the fit may fail or return poor `f1`, `f2`. In that case, fit a polynomial to `example_run_track.dat` directly and substitute those values into Step 5.