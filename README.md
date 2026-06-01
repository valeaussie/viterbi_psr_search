# Viterbi HMM Pulsar Search

## Overview

This repository implements a semi-coherent blind search for radio pulsars in compact binary systems, based on a hidden Markov model (HMM) framework combined with the Viterbi algorithm.

The method models the time-varying apparent spin frequency of a pulsar as a hidden Markov chain. Each coherent segment of the dedispersed time series contributes a normalised Fourier power (Schuster periodogram) observation. The Viterbi algorithm then finds the most probable frequency track through the resulting time-frequency spectrogram, and a log-likelihood statistic is formed to assess significance.

This approach is designed to handle Doppler modulations from short-period binary companions, where the apparent spin frequency drifts by more than one frequency bin over the total observation, making standard FFT-based searches insensitive.

The method is described in:

> O'Leary, Dunn & Melatos (2026), *Discovering pulsars in compact binaries with a hidden Markov model*, arXiv:2601.00500

---

## Repository Structure

```
code/
    viterbi.py                      # Core HMM and Viterbi algorithm
    viterbi_pipeline.py             # Main search script (lean blind mode + full fit mode)
    viterbi_search_and_fit.py       # Full Viterbi search + polynomial + Kepler fit
    do_noise.py                     # Noise generation for null tests
    experiments/
        <experiment_name>/
            jobs/                   # SLURM job scripts for each pipeline stage
            data/                   # Input data and params files
            <run_name>/             # Stage outputs (stage0_inspect/, stage1_clean/, ...)
```

---

## Environment

Activate your environment before running any stage:

```bash
source setup_your_env.sh
```

The Python dependencies are:

```bash
pip install numpy scipy matplotlib sigpyproc psrqpy pyyaml
```

External tools required (installed and on `PATH` or via Apptainer/Singularity):

- `filtool` (PulsarX, via Apptainer SIF)
- `prepdata` (PRESTO)
- `psrfold_fil` (PulsarX, via Apptainer SIF)
- `psrstat` (PSRCHIVE)

---

## Pipeline Overview

The full pipeline runs in six stages. Each stage is implemented as a SLURM job script (or short Python script) inside `code/experiments/<experiment_name>/jobs/`. Stage outputs are written to a versioned experiment directory, e.g. `47Tuc_blind_search_v1/`.

```
Stage 0    Inspect filterbank metadata, query ATNF catalogue
Stage 1    RFI cleaning (filtool) + bandpass diagnostics
Stage 2    Incoherent dedispersion (prepdata, barycentred and topocentric)
Stage 3    Blind Viterbi sweep over subband and Tsft grid (lean mode)
Stage 3a   Aggregate loglike curves, peak-find, deduplicate candidates
Stage 3b   Full Viterbi search + polynomial + Kepler fit on each candidate
Stage 4    Fold each candidate with psrfold_fil (polynomial and Kepler rows)
Stage 5    Collate S/N from folded archives, produce ranked candidate CSV
```

---

## Stage 0: Inspect Filterbank

**Script:** `jobs/stage0_inspect.py`

Reads the filterbank header using `sigpyproc`, computes derived quantities (Nyquist frequency, maximum coherent timescale per known binary), queries the ATNF pulsar catalogue for known pulsars at the target position, and writes observation metadata to YAML.

```bash
python stage0_inspect.py \
    --fil  path/to/observation.fil \
    --exp-dir path/to/experiment_dir
```

Use `--offline` on compute nodes without network access (skips the ATNF query).

**Outputs** (inside `<exp-dir>/`):

```
config/observation.yaml
config/known_pulsars.yaml
stage0_inspect/inspect.log
stage0_inspect/header_raw.txt
provenance/stage0_runinfo.txt
```

---

## Stage 1: RFI Cleaning

**Scripts:** `jobs/stage1_filtool_only.sh`, `jobs/stage1b_bandpass.sh`, `jobs/stage1_clean_postprocess.py`

Stage 1 runs `filtool` (PulsarX) inside an Apptainer container to clean RFI from the filterbank. Stage 1b then computes before/after bandpass diagnostics and writes RFI statistics to YAML.

```bash
sbatch jobs/stage1_filtool_only.sh
sbatch jobs/stage1b_bandpass.sh     # after stage1_filtool_only completes
```

The RFI mitigation chain used in the 47 Tuc example is:

```
-z zdot kadaneF 8 4 kadaneT 8 4 zap 925 960 zap 1525 1612 zap 1675 1720 --fillPatch rand
```

Adapt the zap ranges and algorithms to match your observing band and known RFI environment.

**Outputs** (inside `stage1_clean/`):

```
<rootname>_01.fil          # Cleaned filterbank
filtool.log
bandpass.png               # Before/after bandpass comparison
rfi_stats.yaml
```

---

## Stage 2: Dedispersion

**Script:** `jobs/stage2_dedisp.sh`

Runs PRESTO `prepdata` to produce a dedispersed time series at a single trial DM. Two versions are produced: barycentred (default) and topocentric (`-nobary`), allowing the effect of barycentring on the search to be compared.

```bash
sbatch jobs/stage2_dedisp.sh
```

Edit `DM` and `DM_TAG` in the script to match your target. The output `.inf` file records `tsamp`, the start MJD, and the number of real (non-padded) samples. **Read `tsamp` from the `.inf` file** and use it in all downstream stages.

**Outputs** (inside `stage2_dedisp/<DM_TAG>_bary/` and `stage2_dedisp/<DM_TAG>_topo/`):

```
<rootname>_bary.dat   +   <rootname>_bary.inf
<rootname>_topo.dat   +   <rootname>_topo.inf
```

---

## Stage 3: Blind Viterbi Sweep

**Script:** `jobs/stage3_blind.sh`

Runs `viterbi_pipeline.py` in `--lean` mode over a grid of frequency subbands and coherent timescales $T_\mathrm{sft}$. Each Slurm array task covers one (subband, $T_\mathrm{sft}$) pair. Lean mode saves only the log-likelihood-vs-frequency curve and the top path; it does not compute fits, produce plots, or write candfiles.

The grid in the 47 Tuc example is:

- Subbands: 100 to 800 Hz, width 10 Hz, step 9 Hz (1 Hz overlap between adjacent subbands).
- $N_T$ grid: {16, 32, 64, 128, 256}, giving $T_\mathrm{sft} = T_\mathrm{obs} / N_T$.

```bash
sbatch jobs/stage3_blind.sh
```

The array size `--array=0-N` must equal `n_subbands * n_Tsft - 1`. Recompute if you change the grid.

**Outputs** (inside `stage3_viterbi/blind_v1/Nt<N>/f0_<F>/`):

```
blind_Nt<N>_f0_<F>_loglike_curve.dat    # log-likelihood vs frequency
blind_Nt<N>_f0_<F>_paths.dat            # top Viterbi path(s)
blind_Nt<N>_f0_<F>.params               # parameter file used for this run
```

---

## Stage 3a: Aggregate and Deduplicate Candidates

**Script:** `jobs/aggregate_blind.py`

This step is run manually (not as a SLURM job) after stage 3 completes. It scans all `*_loglike_curve.dat` files in the blind output tree, identifies peaks above a robust threshold (default: 8 MAD-sigma above the median), optionally cross-matches them against known pulsars, deduplicates peaks within a frequency tolerance across subbands and $N_T$ values, and writes a ranked candidate list.

```bash
python jobs/aggregate_blind.py \
    --blind-dir  path/to/stage3_viterbi/blind_v1/ \
    --known-yaml path/to/config/known_pulsars.yaml \
    --out-csv    path/to/stage3_viterbi/blind_v1/candidates.csv \
    --n-sigma    8.0 \
    --dedup-tol-hz 0.5
```

**Outputs:**

```
candidates.csv        # All peaks above threshold, sorted by log-likelihood
candidates_dedup.csv  # Deduplicated candidate list; input to stage 3b
```

Each row in `candidates_dedup.csv` records: `peak_freq_hz`, `peak_loglike`, `multiplicity` (number of raw peaks merged), `nt_values` (which $N_T$ values recovered the candidate), `Nt_best`, `subband_f0`, and `known_match`.

---

## Stage 3b: Full Viterbi Fit

**Script:** `jobs/stage3b_fit.sh`

Runs `viterbi_search_and_fit.py` on each deduplicated candidate. Each task re-runs the Viterbi search in a narrow 10 Hz window centred on the candidate frequency, then fits two independent models to the recovered path:

- A polynomial Taylor expansion (default order 5, giving $f_0 \ldots f_5$) evaluated at the path midpoint.
- A circular-orbit Kepler model (giving $f_\mathrm{spin}$, semi-amplitude $K$, orbital period $P_b$, and zero-crossing time $T_0$), from which Taylor coefficients are derived analytically.

A classification flag (`binary` or `isolated`) is written based on whether the Kepler amplitude is significant relative to the frequency resolution.

```bash
sbatch jobs/stage3b_fit.sh
```

The array size `--array=0-N` must equal the number of rows in `candidates_dedup.csv` minus one.

**Outputs** (inside `stage3_viterbi/candidates/cand_NNN/`):

```
cand_NNN.params
cand_NNN_psrfold.candfile        # 6-column psrfold_fil input (poly row 0, Kepler row 1)
cand_NNN_psrfold.candfile.info   # Sidecar: pepoch, F2, Keplerian parameters
cand_NNN_fit_summary.dat         # Full fit results and classification
cand_NNN_track.dat               # Viterbi path + polynomial fit + Kepler fit
cand_NNN_spectrogram.png
cand_NNN_loglikes.png
```

**Note on the candfile:** `psrfold_fil --candfile` accepts only the 6-column schema `(id, dm, acc, F0, F1, S/N)`. F2 and orbital parameters cannot be passed through the candfile. Read `--pepoch` and `--f2` values from the `.info` sidecar file when folding at stage 4.

---

## Stage 4: Fold Candidates

**Script:** `jobs/stage4_fold.sh`

Folds each candidate twice using `psrfold_fil` (PulsarX inside Apptainer): once using the polynomial fit parameters (row 0 of the candfile) and once using the Kepler-derived Taylor parameters (row 1). Each fold runs in its own subdirectory so that PulsarX auto-named output files do not overwrite each other.

The `--pepoch` and `--f2` values for each row are read from the `.info` sidecar file produced by stage 3b.

```bash
sbatch jobs/stage4_fold.sh
```

The array size `--array=0-N` must equal the number of candidate directories in `stage3_viterbi/candidates/` minus one.

**Outputs** (inside `stage4_fold/cand_NNN/`):

```
poly/
    J0000-00_..._00001.ar      # Folded archive (polynomial fit)
    J0000-00_..._00001.png
    J0000-00_..._00001.cands
kepler/
    J0000-00_..._00001.ar      # Folded archive (Kepler fit)
    J0000-00_..._00001.png
    J0000-00_..._00001.cands
```

---

## Stage 5: Collate S/N

**Script:** `jobs/stage5_collate.sh`, `jobs/collate_fold_snr.py`

Runs `psrstat` on all folded `.ar` files to extract S/N, cross-matches results against `candidates_dedup.csv`, selects the best fit (polynomial or Kepler) per candidate by S/N, and writes a ranked CSV.

```bash
sbatch jobs/stage5_collate.sh
```

**Output:**

```
stage4_fold/fold_snr_ranked.csv
```

Each row records: `cand_id`, `freq_hz`, `snr_best`, `best_fit` (poly or kepler), `snr_poly`, `snr_kepler`, `peak_loglike`, `multiplicity`, `nt_values`, and `known_match`.

---

## Running a Single-Target Search (Quick Start)

If you have a dedispersed `.dat` file and want to run a single Viterbi search without the full pipeline infrastructure, use `viterbi_pipeline.py` directly:

```bash
python code/viterbi_pipeline.py --params search.params
```

Example `search.params`:

```ini
infile         = "path/to/observation_DM24.36.dat"
tsamp          = 7.656e-05
Tsft           = 450
Nsft           = -1
f0             = 370.0
bw             = 20.0
padding_factor = 1
num_harm       = 1
top_paths      = 1
out_prefix     = "run_01"
plot_path      = True
dm             = 24.356441
mjd_start      = 59391.125996527684
```

Any parameter can be overridden from the command line:

```bash
python code/viterbi_pipeline.py --params search.params --Tsft 225 --out_prefix run_02
```

Set `Nsft = -1` to use the full length of the time series. Set `--lean` to skip fits and plots (useful when running many subbands).

### Key parameter choices

- `tsamp` must match the `Width of each time series bin (sec)` value in the PRESTO `.inf` file exactly. A mismatch will scale all recovered frequencies incorrectly.
- `Tsft` must be short enough that the signal frequency drift within one segment is much less than $1/T_\mathrm{sft}$. The maximum coherent timescale for a binary with projected semi-major axis $a_1$ and spin frequency $f_0$ is $T_\mathrm{coh,max} \approx (c / (2\pi f_0 a_1))^{1/2} P_b^{1/2}$ (see O'Leary et al. 2026, equation 4).
- `f0` is the **lower** edge of the search band. The search covers $[f_0,\ f_0 + \mathrm{bw}]$.
- `dm` and `mjd_start` are required only if you want the script to write a `psrfold_fil` candfile. Read `mjd_start` from the `.inf` file, not from the original `.fits` header.

---

## Parameter Reference

| Parameter | Required | Description |
|---|---|---|
| `infile` | yes | Path to float32 PRESTO `.dat` time series |
| `tsamp` | yes | Sampling interval (s); must match `.inf` |
| `Tsft` | yes | Coherent segment length (s) |
| `f0` | yes | Lower edge of search band (Hz) |
| `bw` | yes | Search band width (Hz) |
| `Nsft` | no | Number of segments; `-1` uses full data |
| `padding_factor` | no | SFT zero-padding factor (default 1) |
| `num_harm` | no | Number of harmonics to sum incoherently (default 1) |
| `top_paths` | no | Number of top Viterbi paths to save (default 1) |
| `poly_order` | no | Polynomial fit order (default 5, giving $f_0 \ldots f_5$) |
| `out_prefix` | no | Output file prefix (default `search`) |
| `plot_path` | no | Overplot Viterbi path on spectrogram (default False) |
| `save_delta` | no | Save full Viterbi delta matrix (default False) |
| `dm` | no | DM (pc cm$^{-3}$); required for candfile output |
| `mjd_start` | no | Start MJD of dedispersed time series; required for candfile output |
| `spec_flo`, `spec_fhi` | no | Frequency axis limits for spectrogram plot (Hz) |

---

## Output Files

| File | Description |
|---|---|
| `<prefix>_loglike_curve.dat` | Log-likelihood vs terminating frequency (primary detection product) |
| `<prefix>_paths.dat` | Top Viterbi path(s) with log-likelihood |
| `<prefix>_track.dat` | Path + polynomial fit + Kepler fit vs time |
| `<prefix>_fit_summary.dat` | Full polynomial and Kepler fit results, classification |
| `<prefix>_psrfold.candfile` | `psrfold_fil` input (6-column, poly row 0 + Kepler row 1) |
| `<prefix>_psrfold.candfile.info` | Sidecar: pepoch, F2, Keplerian parameters |
| `<prefix>_spectrogram.png` | Time-frequency spectrogram with Viterbi path overlaid |
| `<prefix>_loglikes.png` | Log-likelihood vs frequency at final time step |

---

## Pitfalls and Sanity Checks

**`tsamp` mismatch.** If `tsamp` in the params file does not match the `.inf` file, all recovered frequencies will be scaled by `tsamp_wrong / tsamp_true`. Check the `.inf` file every time, especially after re-running dedispersion.

**DM consistency.** The DM passed to `psrfold_fil` at stage 4 must equal the DM used for dedispersion at stage 2. Mixing values from different runs will smear the folded profile.

**Search band.** Verify that the true spin frequency lies inside $[f_0,\ f_0 + \mathrm{bw}]$. If folding gives no profile, inspect the `_loglike_curve.dat` or `_loglikes.png` to see where the peak likelihood falls and adjust the band.

**Pepoch convention.** The `pepoch` written by `viterbi_pipeline.py` is `mjd_start + t_ref / 86400`, where `t_ref` is the path midpoint in seconds from the start of the dedispersed time series. A wrong `pepoch` mis-applies `f1` and `f2`, smearing the profile across sub-integrations even when `f0` is correct. Verify by checking that the folded time-vs-phase plot is straight, not curved.

**Candfile F2 limitation.** The `psrfold_fil` 6-column candfile schema does not carry F2. Read `--f2` from the `.info` sidecar file and pass it directly on the `psrfold_fil` command line when folding high-spin-down or short-period binary candidates.

**Blind sweep array size.** The `--array` range in `stage3_blind.sh` must equal `n_subbands * n_Tsft - 1` exactly. Recompute after any change to the frequency grid or $N_T$ list.

**Stage 3b array size.** The `--array` range in `stage3b_fit.sh` must equal the number of rows in `candidates_dedup.csv` minus one. Re-check after re-running `aggregate_blind.py` with different thresholds.

---

## Example: 47 Tucanae Blind Search

The experiment `code/experiments/47_Tuc_O_Alessandro_Ridolfi/` contains all job scripts and configuration files for a blind search of a 4-hour MeerKAT L-band observation of 47 Tucanae (observation `47Tuc_22UL_1of2_L.fil`, DM = 24.356441 pc cm$^{-3}$). The search covers 100 to 800 Hz in 78 subbands of width 10 Hz (step 9 Hz) at five coherent timescales ($N_T \in \{16, 32, 64, 128, 256\}$), for a total of 390 Slurm array tasks at stage 3.