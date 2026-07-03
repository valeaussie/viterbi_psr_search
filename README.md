# Viterbi HMM Pulsar Search

## Overview

This repository implements a semi-coherent blind search for radio pulsars in compact binary systems, based on a hidden Markov model (HMM) solved with the Viterbi algorithm.

The method models the time-varying apparent spin frequency of a pulsar as a hidden Markov chain. The observation for each coherent time segment is the normalised Fourier power (Schuster periodogram) of the dedispersed, heterodyned, and downsampled radio intensity time series. The Viterbi algorithm finds the most probable frequency track through the resulting time-frequency spectrogram and returns a log-likelihood statistic used to assess significance.

The approach is designed to handle the Doppler modulation from a short-period binary companion, where the apparent spin frequency drifts across many Fourier bins over the total observation, making standard FFT-based searches insensitive. Because the transition model is an unbiased random walk (equal probability of moving by -1, 0, or +1 frequency bin per coherent segment), the pipeline does not search over orbital parameters explicitly. The coherent segment length is instead chosen short enough that the frequency drift per segment is at most one bin, for the widest binary orbits of interest.

The method is described in:

> [O'Leary, Dunn & Melatos (2026)](https://iopscience.iop.org/article/10.3847/1538-4357/ae3288)

---

## Repository Structure

```
code/
    viterbi.py                      # Core HMM: Viterbi recursion, numpy implementation
    viterbi_pipeline.py             # Main search script: SFT computation, lean blind mode
    viterbi_search_and_fit.py       # Full Viterbi search + polynomial + Kepler orbital fit
    do_noise.py                     # Noise-only null test helper
    experiments/
        <experiment_name>/
            jobs/                   # SLURM job scripts and helper Python scripts
                stage0_inspect.py
                stage1a_inspect_dynamic_spectrum.py
                stage1b_filtool_only.sh
                stage1c_bandpass.sh (calling stage1_clean_postprocess.py)
                stage2_dedisp.sh
                stage3_blind.sh
                stage3_viterbi.sh
                stage3b_fit.sh
                stage4_fold.sh
                stage5_collate.sh
                aggregate_blind.py      # Stage 3a candidate aggregation
                collate_fold_snr.py     # Stage 5 S/N extraction
                make_candidate_gallery_base64.py
                gen_search_params.py
                patch_known_pulsars.py
            data/                   # Input filterbank, params files
            config/                 # observation.yaml, known_pulsars.yaml
            <run_name>/             # Stage outputs (stage0_inspect/, stage1_clean/, ...)
```

---

## Environment

project environment:

```bash
source /fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh
```

Python dependencies:

```
numpy  scipy  matplotlib  sigpyproc  psrqpy  pyyaml
```

External tools required (available via the setup script or Apptainer SIF):

- `filtool` (PulsarX)
- `prepdata` (PRESTO)
- `psrfold_fil` (PulsarX)
- `psrstat` (PSRCHIVE)

---

## Pipeline Overview

The pipeline runs in six stages. Each stage is a SLURM job script inside `jobs/` (except stage 0 which is just an inspection). Outputs go to a versioned experiment directory, e.g. `47Tuc_blind_search_v1/`.

### Stage 0: Inspect and configure

```bash
python stage0_inspect.py --fil <path>.fil --exp-dir <experiment-dir>
```

Reads the filterbank header using `sigpyproc`, computes derived quantities (Nyquist frequency, Chandler ceiling per pulsar), queries the ATNF catalogue for known pulsars in the field, and writes:

- `config/observation.yaml`: telescope and observation parameters
- `config/known_47tuc_pulsars.yaml`: known pulsars with F0, DM, and binary parameters
- `stage0_inspect/inspect.log`

Use `--offline` on compute nodes without network access.

> **Note:** The ATNF catalogue at the moment lags behind and some pulsars in the 47 TUC are missing. I have added them to
> `known_47tuc_pulsars.yaml` running:
>
> ```bash
> python jobs/patch_known_pulsars.py \
>     --yaml config/known_47tuc_pulsars.yaml \
>     --out  config/known_47tuc_pulsars.yaml
> ```
>
> The script currently adds the 15 new 47 Tuc pulsars from
> [Chen & Risbud et al. (2026)](https://www.aanda.org/articles/aa/pdf/2026/06/aa59650-26.pdf).
> Use `--dry-run` to preview changes without writing.

### Stage 1: RFI

Here we check if there is RFI and remove it if needed.

To check whether the data actually has RFI
`stage1_inspect_dynamic_spectrum.py` running:
```bash
python stage1a_inspect_dynamic_spectrum.py \
    --fil /fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/data/47Tuc_22UL_1of2_L.fil \
    --flo 880 \
    --fhi 970 \
    --nsamples 10000 \
    --outdir /fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1/stage1_clean
```
This produces `dynspec_rfi_ratio.png` and `dynspec_stats.txt` in the output directory.
Check `dynspec_rfi_ratio.png`: if any channel's ratio exceeds 1.05 (the red dotted line),
that frequency range is RFI-elevated and should be hard-zapped. The exact elevation values
are in `dynspec_stats.txt`. The `--flo` and `--fhi` values above are examples for the
47 Tuc L-band observation; adjust them to cover the known RFI environment of your need.
If the ratio test shows elevation above baseline, run `stage1b_filtool_only.sh` to clean 
then `stage1c_bandpass.sh` to verify

### Stage 2: Dedispersion

```bash
sbatch --array=0-<N_DM-1> stage2_dedisp.sh
```

Runs `prepdata` (PRESTO) to produce a dedispersed and barycentred time series for each DM trial. The DM grid is defined inside `stage2_dedisp.sh`. Each array task handles one DM value and writes a `.dat` and `.inf` file to `stage2_dedisp/DM<XX.XX>_bary/`.

### Stage 3: Blind Viterbi search

```bash
sbatch --array=0-<N_jobs-1> stage3_blind.sh
```

The array size must equal `N_subbands * N_Nt_values * N_DM - 1`. Each task calls `viterbi_pipeline.py` in lean blind mode for one (DM, Nt, subband) combination, writing a `*_loglike_curve.dat` and a sibling `*.params` file to `stage3_viterbi/blind_v1/DM<XX.XX>/Nt<N>/f0_<F>/`.

The `*.params` file records `Tsft`, `out_prefix`, and other search parameters. It is read by `aggregate_blind.py` to infer `T_obs` without any additional inputs.

### Stage 3a: Candidate aggregation

```bash
python aggregate_blind.py \
    --blind-dir  <exp-dir>/stage3_viterbi/blind_v1 \
    --known-yaml <exp-dir>/config/known_47tuc_pulsars.yaml \
    --out-dir    <exp-dir>/stage3_viterbi/blind_v1 \
    --false-alarm-rate 0.1
```

Reads every `*_loglike_curve.dat` file, groups them by `Nt`, and computes a separate detection threshold for each `Nt` group following the exponential-tail method of O'Leary et al. (2026), Equations 13-17 and Appendix D.

**Threshold method.** The Viterbi log-likelihood scales with `Nt` (it is a cumulative sum over `Nt` coherent segments), so a single global threshold across all `Nt` values is not meaningful. Within each `Nt` group, all loglike values from all subbands and DM trials are pooled. One can choose the percentile and the tail above that percentile is fitted as an exponential using the maximum-likelihood estimator for the rate parameter (in out run we chose 95% percentile). The threshold is then set by inverting the desired false alarm probability per subband (`--false-alarm-rate`, default 0.1). The number of independent Viterbi runs used to build the noise distribution is taken to be the number of curve files in that `Nt` group (the role of `N_real` in the paper).

**Deduplication.** Candidates are deduplicated in frequency with a tolerance derived automatically from the coarsest bin width in the `Nt` grid: `tol = 2 * N_T_max / T_obs`. This ensures that the same pulsar detected at different `Nt` values (whose peak frequencies can differ by up to one coarse bin due to quantisation) is correctly merged into a single candidate, without over-merging distinct pulsars.

Outputs:

- `candidates_raw.csv`: one row per above-threshold peak, with DM, Nt, subband, frequency, loglike, and known-pulsar cross-match.
- `candidates_dedup.csv`: deduplicated candidate list ranked by loglike, with columns `dm_count` and `multiplicity` indicating how many DM trials and (Nt, subband) combinations recovered the candidate.
- `threshold_calibration.txt`: full calibration record (lambda, L_tail, L_th, sigma_L_th) for each Nt group.

### Stage 3b: Full fit

```bash
sbatch --array=0-<N_candidates-1> stage3b_fit.sh
```

Runs `viterbi_search_and_fit.py` on each candidate from `candidates_dedup.csv`. Performs a refined Viterbi search followed by a polynomial frequency evolution fit and optionally a Kepler orbital fit. Writes per-candidate result files to `stage3_viterbi/candidates/`.

### Stage 4: Folding

```bash
sbatch --array=0-<N_candidates-1> stage4_fold.sh
```

Runs `psrfold_fil` (PulsarX, inside the Apptainer SIF) to fold the cleaned filterbank at each candidate frequency and DM. Writes `.ar` files to `stage4_fold/`.

### Stage 5: Collate S/N

```bash
sbatch stage5_collate.sh
```

Runs `collate_fold_snr.py`, which calls `psrstat` (PSRCHIVE) on all folded `.ar` files to extract S/N and produce a ranked CSV at `stage4_fold/fold_snr_ranked.csv`.

---

## Key Parameters

| Parameter | Typical value | Notes |
|-----------|--------------|-------|
| Subband width | 10 Hz | `bw` in `stage3_blind.sh` |
| Search band | 50 to 1050 Hz | Covers all known MSPs |
| `Nt` grid | 16, 32, 64, 128, 256 | Coherent segment counts |
| DM grid | 21 trials | Centred on cluster DM |
| False alarm rate | 0.1 per subband | `--false-alarm-rate` in Stage 3a |
| Dedup tolerance | auto (2 x coarsest bin) | `--dedup-bin-widths 2` |

---

## Application to 47 Tucanae

### Data

The pipeline has been applied to a Parkes radio telescope observation of the globular cluster 47 Tucanae (47 Tuc), observed at L-band (centre frequency ~1374 MHz) using the Murriyang multibeam system. The filterbank records Stokes I intensity with 96 frequency channels across a 288 MHz bandwidth, 1-bit digitisation, and a sampling time of approximately 64 microseconds. The total observation spans roughly 4 hours.

The data were processed with the call path described in Stages 0-5 above, using the experiment directory `47Tuc_blind_search_v1`.

### Search configuration

- DM grid: 21 trials spanning approximately 23.4 to 25.4 pc/cm^3, centred on the cluster DM of 24.36 pc/cm^3.
- `Nt` grid: 16, 32, 64, 128, 256 (corresponding to coherent timescales from ~56 s to ~900 s).
- Subband width: 10 Hz; total band: 50-1050 Hz; 100 subbands.
- Total Viterbi jobs: 21 DM x 5 Nt x 100 subbands = 10,500 per polarisation.

### Known pulsars in 47 Tucanae

47 Tuc is the most pulsar-rich globular cluster in the Southern sky. The table below lists all 42 known pulsars as of early 2026, drawn from the ATNF pulsar catalogue and from Chen & Risbud et al. (2026, A&A), who reported 15 new millisecond pulsars discovered with MeerKAT.

| Name | Period (ms) | F0 (Hz) | DM (pc/cm^3) | Binary type |
|------|------------|---------|-------------|-------------|
| J0024-7204C  |   5.757 |  173.71 | 24.591 | isolated |
| J0024-7204D  |   5.358 |  186.65 | 24.741 | isolated |
| J0024-7204E  |   3.536 |  282.78 | 24.240 | DD |
| J0024-7204F  |   2.624 |  381.16 | 24.384 | isolated |
| J0024-7204G  |   4.040 |  247.50 | 24.434 | isolated |
| J0024-7204H  |   3.210 |  311.49 | 24.375 | DD |
| J0024-7204I  |   3.485 |  286.94 | 24.430 | ELL1 |
| J0024-7204J  |   2.101 |  476.05 | 24.594 | BTX |
| J0024-7204L  |   4.346 |  230.09 | 24.399 | isolated |
| J0024-7204M  |   3.677 |  271.99 | 24.426 | isolated |
| J0024-7204N  |   3.054 |  327.44 | 24.557 | isolated |
| J0024-7204O  |   2.643 |  378.31 | 24.358 | BTX |
| J0024-7204P  |   3.643 |  274.50 | 24.290 | BT |
| J0024-7204Q  |   4.033 |  247.94 | 24.279 | ELL1 |
| J0024-7204R  |   3.480 |  287.32 | 24.361 | ELL1 |
| J0024-7204S  |   2.830 |  353.31 | 24.382 | ELL1 |
| J0024-7204T  |   7.588 |  131.78 | 24.421 | ELL1 |
| J0024-7204U  |   4.343 |  230.26 | 24.340 | ELL1 |
| J0024-7204V  |   4.810 |  207.89 | 24.105 | BTX |
| J0024-7204W  |   2.352 |  425.11 | 24.370 | BTX |
| J0024-7204X  |   4.772 |  209.58 | 24.538 | ELL1 |
| J0024-7204Y  |   2.197 |  455.24 | 24.475 | ELL1 |
| J0024-7204Z  |   4.554 |  219.57 | 24.450 | isolated |
| J0024-7204aa |   1.840 |  543.48 | 24.921 | isolated |
| J0024-7204ab |   3.705 |  269.93 | 24.326 | isolated |
| J0024-7204ac |   2.740 |  364.96 | 24.460 | BT |
| J0024-7204ad |   3.740 |  267.38 | 24.410 | BT |
| J0024-7204ae |   3.870 |  258.40 | 24.340 | He WD |
| J0024-7204af |   2.990 |  334.45 | 24.340 | BW/RB |
| J0024-7204ag |   9.760 |  102.46 | 24.410 | binary |
| J0024-7204ah |   3.070 |  325.73 | 24.360 | binary |
| J0024-7204ai |  13.030 |   76.75 | 24.470 | binary (ecc) |
| J0024-7204aj |   6.360 |  157.23 | 24.380 | isolated |
| J0024-7204ak |   3.520 |  284.09 | 23.910 | binary |
| J0024-7204al |   2.670 |  374.53 | 24.110 | BW |
| J0024-7204am |   4.160 |  240.38 | 24.550 | binary |
| J0024-7204an |   2.610 |  383.14 | 24.120 | binary |
| J0024-7204ao |   1.880 |  531.91 | 23.650 | isolated |
| J0024-7204ap |   5.110 |  195.69 | 24.360 | isolated |
| J0024-7204aq |   3.040 |  328.95 | 23.630 | binary |
| J0024-7204ar |   9.760 |  102.46 | 24.160 | binary |
| J0024-7204as |   4.020 |  248.76 | 24.660 | binary |

Pulsars ae through as (15 entries, source column `Chen+2026`) were discovered by Chen & Risbud et al. (2026) using MeerKAT and are not yet in the ATNF catalogue. They are included in `config/known_47tuc_pulsars.yaml` via `patch_known_pulsars.py`.

Binary classification codes: DD = double neutron star; ELL1 = low-eccentricity binary (ELL1 timing model); BTX = black-widow/redback with complex orbital model; BT = Blandford-Teukolsky; BW = black widow; RB = redback; He WD = helium white dwarf companion.

### Notes on the threshold calibration

The exponential-tail threshold used in Stage 3a was calibrated empirically from the 47 Tuc observation itself rather than from off-source data or Monte Carlo noise simulations (which are not available for this dataset). Because 47 Tuc contains many bright pulsars, the extreme tail of the pooled loglike distribution is contaminated by real signals, which can cause the fitted exponential rate to be underestimated and the threshold to be inflated. The threshold should therefore be treated as conservative. If a candidate's loglike is only marginally above the threshold, it should not be dismissed: follow it up with folding regardless.

---

## Candidate Vetting

After Stage 5, the ranked CSV at `stage4_fold/fold_snr_ranked.csv` is the primary output. To assess whether a candidate is a real new pulsar:

1. Check `known_match` in `candidates_dedup.csv`. If it matches a known pulsar, the detection is a re-detection, which is still useful for verifying pipeline sensitivity.
2. Check `dm_count`. A real pulsar should peak in significance at or near the cluster DM (~24.36 pc/cm^3) and be recovered across multiple DM trials. Detections at only one DM trial are suspect.
3. Check `multiplicity`. A real pulsar should appear in multiple subbands and ideally multiple `Nt` values.
4. Inspect the folded `.ar` file. A real pulsar shows a coherent pulse profile and dispersion sweep consistent with the cluster DM.
5. Check for harmonics. If a candidate at frequency `f` has a corresponding detection at `f/2` or `2f`, one is likely a harmonic of the other.