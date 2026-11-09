# Reproducing the 47 Tuc analysis

This repository contains the analysis programs for *Applying a Hidden Markov Model Search for Compact Binary Pulsars to a MeerKAT Observation of 47 Tucanae*. This guide follows the manuscript inspected on 9 November 2026. The previous repository guide is preserved in [README_old.md](README_old.md); use this guide for the paper's current selections and original-filterbank folds.

The analysis has several distinct parts:

| Part | Search input and time coordinate | Role in the paper |
|---|---|---|
| Original targeted recovery of O | PulsarX time series at DM 24.356441, topocentric | Tracks, parameter estimates, four folding models and a coherence-timescale scan |
| Self-blinded HMM search | 21 PRESTO barycentric time series from statistically cleaned data | Original selection: D, F, I, J and O |
| Supplementary selection | Saved blind HMM endpoint scores at 16 segments | Adds C; E and third-harmonic L remain provisional |
| PRESTO comparison | The same 21 cleaned-data barycentric time series | Direct one-harmonic acceleration searches; seven secure recoveries, with H and L provisional |
| Catalogue-informed HMM follow-up | Separate narrow searches of the barycentric data | Additional folding diagnostics, including C/L and later J/O follow-up; does not change blind outcomes |
| Final comparison folds | Original channelised filterbank | Native PRESTO reference and candidate folds used in the current comparison table |

No new pulsar is confirmed. HMM likelihoods are ranking scores, without a calibrated search-wide false-alarm probability. PRESTO profile significances quoted here are nominal folded-profile statistics. Differences in selection, interval length and follow-up prevent interpreting the comparison as matched sensitivity.

## What is supplied, and what is still needed

Git contains the retained scientific programs, required source copies, fixed catalogue configuration and the tested reference timing models. It excludes observations, containers, Python environments, search grids and results, large model caches, native folding products, scheduler logs and archived runs. Generated files named below are expected outputs or inputs to a later stage, not promises that those files are present in a fresh clone.

Obtain the original MeerKAT observation from the observation owners, including permission to use it. No verified public download for this exact filterbank is supplied here. The expected file is `data/47Tuc_22UL_1of2_L.fil` inside the experiment. It is an 8-bit SIGPROC filterbank with one intensity stream, 128 channels, `fch1=1708.23828125` MHz, `foff=-6.6875` MHz, sampling time approximately `7.65607476635514e-5` s and start MJD `59391.125996527684`. It contains 188,091,904 samples, lasts 14,400.4568 s and occupies 24,075,764,064 bytes including its 352-byte header. This identifies the input; a similarly named observation is not interchangeable.

**One preparation step remains unresolved.** The original targeted O search used `redetection_of_47TucO/47TucO_redetect_cfbf00000_Plan1_1_DM24.36.dat` and its matching `.inf`. The retained `.inf` records DM **24.356441**, despite the rounded filename, no barycentring, 188,026,156 samples and sampling time `0.00007656074766355` s. The paper specifies PulsarX `dedisperse_all_fil`, Kadane filtering and zero-DM matched filtering (`--kadaneF 8 4 --zdot`). The exact original launcher and dedispersion-plan file have not been located among the inspected retained scripts. Obtain that invocation/plan, or the original prepared time series, from the authors before claiming an exact rebuild of this experiment from the raw filterbank. The smaller sample count must not be replaced silently with the original filterbank's length. The commands below reproduce the downstream targeted analysis once this input is available.

Several follow-up programs are preserved replay drivers. They require generated manifests, original selected tracks and model summaries, and sometimes the preceding inspection-fold products which they protect against overwriting. These are operational inputs, not alternative final results. Recreate them in the dependency order below, or obtain the corresponding scientific-input bundle from the authors. A clone alone is not a self-contained archive of those inputs. Do not reuse a historical manifest with different files or timestamps: prepare a new manifest for a fresh run. No new orchestration or analysis code is introduced by this documentation update.

## Software and local paths

The HMM needs Python, NumPy, SciPy, Matplotlib and PyYAML. Saved timing measurements used Python 3.10.4, NumPy 1.23.5, SciPy 1.15.3 and Matplotlib 3.10.8. The PRESTO container recorded Python 3.10.12, NumPy 1.23.5, SciPy 1.15.2 and Matplotlib 3.10.0. Preserve the environment used for a reproduction: small numerical changes can alter tie-breaking, nonlinear fits and candidate ordering.

Other requirements are:

- PulsarX `filtool` and `dedisperse_all_fil`; `psrfold_fil` and Tempo2 are needed if reproducing the PulsarX timing-predictor procedures. See the [PulsarX source and installation instructions](https://github.com/ypmen/PulsarX). The retained cleaning script names the `pulsarx_20241219.sif` image. Obtain that image separately for the original environment. The upstream Docker image or a source build provides a usable alternative, but is not identified here as the same build.
- PRESTO `prepdata`, `realfft`, `rednoise`, `accelsearch`, `prepfold` and its Python bindings, including `presto.presto`, `presto.sigproc` and `presto.prepfold`. The production comparison ran inside `pulsar-miner-2.0-prestozl-22feb2025_latest.sif`. This filename identifies the retained environment, not a fully pinned public software release. Obtain the exact image from the authors for numerical replay. For rebuilding, use the [PRESTO source](https://github.com/scottransom/presto) and its installation instructions at the historical revision used by that image. Current PRESTO has changed timing dependencies; installing today's release is not evidence of reproducing this build.
- Classic pulsar-timing **TEMPO**, its runtime files, Solar-system ephemerides, clock corrections and MeerKAT observatory coordinates. The native folding helpers invoke `tempo -z`, use MeerKAT site `m`, barycentric site `@`, and require the PRESTO MeerKAT code `MK`. Candidate time conversion uses DE405 and `UTC(NIST)`. Reference ephemerides retain their supplied timing conventions and may require other ephemerides, including DE421. See the [TEMPO reference manual](https://tempo.sourceforge.net/reference_manual.html).
- Ghostscript and the plot conversion programs invoked by the retained helpers, including `ps2pdf`, `pdfinfo` and `pdftoppm`. Slurm and `srun` are required by the saved batch launchers and some timing wrappers. PSRCHIVE is used by older scripts, not for the final comparison-table statistics.

The image also contains PULSAR_MINER, whose [upstream repository](https://github.com/alex88ridolfi/PULSAR_MINER) describes its own pipeline. **The comparison reported here did not call that wrapper.** Checking `pulsar_miner -version` during software inspection was not a PULSAR_MINER search. The unrelated `data/data.cfg` is not the configuration of the reported PRESTO comparison.

For portable command-line programs, set paths in your shell:

```bash
export REPO=/absolute/path/to/your/viterbi_psr_search
export EXP="$REPO/code/experiments/47_Tuc_O_Alessandro_Ridolfi"
export BLIND="$EXP/47Tuc_blind_search_v1"
export JOBS="$EXP/jobs"
export REDET="$EXP/redetection_of_47TucO"
export PROD="$EXP/psr_miner_search/hmm_presto_v2_20261027T060742Z/production"
export FOLLOW="$PROD/followup"
export RAW="$EXP/data/47Tuc_22UL_1of2_L.fil"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg
```

`REPO` is your own checkout. Do not copy the author's environment directory. Install the dependencies in your own environment and make their executables and Python modules available.

Some retained drivers embed the original OzSTAR repository prefix. An alternative to changing their source is to map your checkout into that prefix inside a container. [Apptainer bind mounts](https://apptainer.org/docs/user/main/bind_paths_and_mounts.html) support a different host and container path:

```bash
export PAPER_PREFIX=/fred/oz022/vdimarco/software/install/viterbi_psr_search
export PRESTO_SIF=/absolute/path/to/pulsar-miner-2.0-prestozl-22feb2025_latest.sif
apptainer exec --bind "$REPO:$PAPER_PREFIX" "$PRESTO_SIF" bash
# Inside this shell, set REPO to PAPER_PREFIX and set EXP, BLIND, etc. again.
```

This needs no OzSTAR directory on the host. If the observation is outside the checkout, bind it explicitly onto the expected `data/47Tuc_22UL_1of2_L.fil` path too. Use the Python entry points inside an allocation with the dependencies active. The original shell launchers additionally contain `module load`, setup-script, image, working-directory and log paths; treat these as site-specific submission templates. Configure those paths in a personal launcher outside Git, or bind a local setup script onto their expected setup location. Do not assume exported variables override shell variables that a retained script assigns unconditionally. Nested container launches need their image and bind paths configured separately.

All execution commands below are instructions for a subsequent reproduction. This documentation task did not run them. On OzSTAR, read the current official documentation and project/account/QoS limits before submission. Use measured representative runtimes with a modest margin for wall-time requests, and let Slurm enforce permitted concurrency. The wall times printed in old templates are not recommendations for a new run.

## 1. Original targeted recovery of O

### Track and coherence scan

Keep this experiment under `redetection_of_47TucO`, separate from `47Tuc_blind_search_v1/redetect_O_J_matched_20261105`. Use the topocentric PulsarX `.dat` and `.inf` identified above. The search band is 376.31--380.31 Hz, with one harmonic, no spectral padding, a fifth-order frequency fit and a three-choice frequency transition: minus one, zero or plus one Fourier bin. The main track has 32 segments.

The exact retained duration-based scan is in `grid_search.sh`. Its segment counts and requested durations are:

| Segments | Duration (s) | Segments | Duration (s) | Segments | Duration (s) |
|---:|---:|---:|---:|---:|---:|
| 16 | 899.713885 | 31 | 464.368420 | 64 | 224.928433 |
| 20 | 719.771093 | 32 | 449.856943 | 128 | 112.464217 |
| 24 | 599.809282 | 40 | 359.885508 | 256 | 56.232108 |
| 28 | 514.122187 | 48 | 299.904603 | 512 | 28.116016 |
| 30 | 479.847395 | 56 | 257.061055 | 1024 | 14.058008 |

For example, run this in the prepared environment and an appropriate compute allocation:

```bash
cd "$REDET"
python "$REPO/code/viterbi_search_and_fit.py" \
  --infile "$REDET/47TucO_redetect_cfbf00000_Plan1_1_DM24.36.dat" \
  --tsamp 0.00007656074766355 --Tsft 449.856943 --Nsft 32 \
  --f0 376.31 --bw 4 --poly-order 5 --plot-path \
  --dm 24.356441 --mjd-start 59391.125996527684038 \
  --known-label J0024-7204O --out_prefix 47TucO_redetect_Nt32
```

Repeat with the table's `--Tsft` and `--Nsft` pairs and corresponding output prefixes. Do not substitute the later `--segments` mode if reproducing the original duration-based scan. The outputs include `_track.dat`, `_loglikes.dat`, `_fit_summary.dat`, `_paths.dat`, spectrogram and likelihood PDFs, and sampling metadata.

Each segment has its mean removed, is heterodyned to the band centre and Fourier-resampled. The code truncates an odd resampled length to an even length. It normalises Fourier power by the resampled segment variance times its sample count and labels frequencies using the actual resampling interval. Keep these implementation details when comparing track bins.

The expected main result is a quasi-sinusoidal track with peak score approximately 564.9. The central track sample is the reference epoch, at index `N_T // 2`, rather than the arithmetic midpoint between the two central samples. The polynomial and circular fits provide all six local frequency coefficients through `f5`. The circular implementation uses a sine convention; converting its epoch to the paper's cosine/ascending-node convention adds one quarter of the fitted orbital period. Adding that offset to the start MJD does not barycentre an epoch.

### Targeted noise reference

Use `compute_Lth.py` after each search, with the actual number of bins in its likelihood file. `stage_Lth.sh` records 10,000 noise realisations, nominal per-subband alpha 0.1, 200 bootstrap samples and seed equal to the segment count. The targeted calculation uses the default tail percentile 99.99; this differs from the blind pooled 95th-percentile selection.

```bash
python "$REDET/compute_Lth.py" --NT 32 --Nbin 1798 --Nreal 10000 \
  --alpha-prime 0.1 --nboot 200 --seed 32 --Tcoh 449.856943 \
  --out /absolute/path/to/new/Lth_Nt32.csv
```

`collate_Lth.py` combines the per-timescale CSVs. The paper quotes a reference level of 95.4 +/- 0.1 for the main targeted case. It is not a calibrated false-alarm threshold for this real-data search.

### Fits and native folds

Use `folding_with_PRESTO/make_polyco_fold.py` for both Taylor-based models. It integrates the six frequency coefficients into a sixth-order phase polynomial, retaining the central reference epoch. Its polycos are topocentric and must not receive another barycentric correction.

```bash
python "$REDET/folding_with_PRESTO/make_polyco_fold.py" \
  --fit "$REDET/47TucO_redetect_Nt32_fit_summary.dat" \
  --inf "$REDET/47TucO_redetect_cfbf00000_Plan1_1_DM24.36.inf" \
  --fil "$RAW" --outdir /absolute/path/to/new/O_Nt32_folds --name 47TucO_Nt32
```

Run the two `prepfold` commands written to `prepfold_polyco_commands.txt` in the intended output directory. They use `-filterbank -topo`, the generated polycos, and local refinement with `-npfact 2 -ndmfact 3 -searchfdd`. The script's default radio-subband setting is retained for this original experiment; do not silently replace it with the later comparison's explicit `-nsub 32`.

`make_prepfold_inputs.py` also writes commands for the full circular orbit, converting `K` to `x=K P_b/(2 pi f_spin)` and the sine epoch to the ascending-node convention. Use its **kepler_orbit** command for that panel. Its direct `-f/-fd/-fdd` polynomial and Taylor commands truncate higher terms and are not the final six-coefficient Taylor folds. `poly_grid_fold.sh` is likewise not the PRESTO six-coefficient fold sequence used in the current scan table.

Fold the independent observation-epoch O ephemeris separately. Its retained input is `J0024-7204O_obsepoch.par`, also preserved as the `paper_O_obsepoch_BT` reference model in the reference-fold package below. A barycentric ephemeris must be evaluated at the topocentric sample times with the timing software; it is not a topocentric fitted-frequency model.

The current four-panel result is 59.4, 26.1, 11.9 and 11.3 nominal PRESTO sigma for the independent reference, polynomial, Keplerian Taylor and full-orbit folds. `make_ntscan_folds.sh` generates polynomial polycos and commands for all 15 scan points. The current scan's strongest fold is 26.1 at 32 segments. Do not use old PSRCHIVE S/N summaries as these PRESTO profile statistics.

## 2. Prepare the self-blinded input

### Statistical RFI cleaning

`jobs/stage1_rfi_cleaning.sh` reads the original filterbank and runs PulsarX `filtool`. Its scientific command, with local paths substituted, is:

```bash
mkdir -p "$BLIND/stage1_clean"
cd "$BLIND/stage1_clean"
filtool --threads 8 --ra 00:24:04.6500 --dec=-72:04:53.8000 \
  --telescope MeerKAT --rootname 47Tuc_22UL_1of2_L_RFIcleaned \
  --source_name J0024-7204O \
  -z zdot kadaneF 8 4 kadaneT 8 4 --fillPatch rand -f "$RAW"
```

The required output is `47Tuc_22UL_1of2_L_RFIcleaned_01.fil`. This is statistical mitigation without predefined frequency masks or hard-zapped bands. The searched data were produced by job 13920756 on 6 July 2026. The older June hard-zapped preparation and its `rfi_stats.yaml` do not describe this input. The cleaning uses random replacement; the retained launcher does not specify a random seed. Replaying that command is therefore not a promise of bit-for-bit identical cleaned samples. Obtain the retained cleaned input if exact numerical replay is required.

### Barycentric dedispersion

`jobs/stage2_dedisp_new.sh` invokes PRESTO `prepdata` at the 21 DMs in `jobs/dm_list.txt`, from 23.50 to 25.50 in steps of 0.10 pc cm^-3. There is no `-nobary`: barycentring is enabled.

```bash
while IFS= read -r dm; do
  case "$dm" in ''|\#*) continue ;; esac
  tag=$(printf 'DM%05.2f' "$dm")
  out="$BLIND/stage2_dedisp/${tag}_bary"
  mkdir -p "$out"
  (cd "$out" && prepdata -dm "$dm" -o "47Tuc_${tag}_bary" \
    "$BLIND/stage1_clean/47Tuc_22UL_1of2_L_RFIcleaned_01.fil")
done < "$JOBS/dm_list.txt"
```

Run through the corresponding Slurm array for production; the loop shows the exact per-DM operation. Each directory contains `.dat`, `.inf` and the saved dedispersion log when launched through the retained script. All 21 retained logs identify the July cleaned filterbank. Keep each `.inf` beside its `.dat`: it supplies the epoch, DM, sampling interval, barycentric flag and valid-sample range. Read the On/Off pairs, not just the padded `Number of bins` field. The stored series have 192,000,000 samples; their trailing padding is not observed signal. Barycentric epochs are DM-dependent and differ from the raw topocentric epoch.

## 3. Self-blinded HMM grid and original selection

The grid has 112 bands with starts `range(50,1050,9)` Hz, width 10 Hz, and segment counts `16,32,64,128,256,512,1024`. Its effective frequency coverage is 50--1059 Hz and it contains 16,464 configurations. `num_harm=1`, `padding_factor=1`, `top_paths=1` and lean output are used.

The current code's `--segments N` mode divides the valid sample count by `N` using integer division, constructs exactly `N` equal segments and leaves fewer than `N` trailing samples unused. It records the actual sampling in `_sampling.json`. It does not derive the segment duration from PRESTO's padded length. Uniform starting-state weights and the equal three-choice transition link the normalised per-segment powers; harmonic powers are not summed.

The saved production data live under `stage3_viterbi/blind_v2`. Choose a fresh output location when reproducing; the measured launcher refuses to overwrite completed measurements.

```bash
export SEARCH_CODE_DIR="$PROD/frozen/code"
export SEARCH_OUTPUT_BASE="$BLIND/stage3_viterbi/blind_v2"
export SEARCH_GRID_MANIFEST="$PROD/grid_manifest.json"
export SEARCH_GRID_SIZE=16464
export SEARCH_RUN_ID=47tuc-reproduction
export SEARCH_CANDIDATE_DIR="$BLIND/stage3_viterbi/candidates"
python "$JOBS/make_search_grid.py" --jobs-dir "$JOBS" --out "$SEARCH_GRID_MANIFEST"
```

`stage3_blind.sh` maps a global index to DM, segment count and subband, writes a `.params` file and runs `measure_search.py` around `viterbi_pipeline.py --lean` on one CPU. `stage3_batch.sh` groups configurations into scheduler tasks; `SEARCH_CONFIGS_PER_WORKER` controls work per allocation, not simultaneous-job concurrency. Set the array size from the grid and current Slurm limits. Per-configuration outputs are `_loglike_curve.dat`, `_paths.dat`, `_sampling.json`, `_timing.json` and `.params` under `DMxx.xx/NtN/f0_FFF.FFF/`.

After every configuration has completed, run the original selection with explicit settings:

```bash
python "$JOBS/aggregate_blind.py" \
  --blind-dir "$SEARCH_OUTPUT_BASE" --known-yaml "$BLIND/config/known_47tuc_pulsars.yaml" \
  --out-dir "$SEARCH_CANDIDATE_DIR" --expected-curves 16464 \
  --expected-grid "$SEARCH_GRID_MANIFEST" --run-id "$SEARCH_RUN_ID" \
  --tail-percentile 95 --false-alarm-rate 0.1 --dedup-tol-hz 0.035
python "$JOBS/classify_known.py" \
  --cands "$SEARCH_CANDIDATE_DIR/candidates_timescale_alternatives.csv" \
  --known-yaml "$BLIND/config/known_47tuc_pulsars.yaml" \
  --out "$SEARCH_CANDIDATE_DIR/classified_timescale_alternatives.csv"
python "$JOBS/flag_rfi.py" \
  --classified-csv "$SEARCH_CANDIDATE_DIR/classified_timescale_alternatives.csv" \
  --known-yaml "$BLIND/config/known_47tuc_pulsars.yaml" \
  --cand-dir "$SEARCH_CANDIDATE_DIR" --out-csv "$SEARCH_CANDIDATE_DIR/rfi_flags.csv"
```

`stage3_aggregate.sh` wraps this sequence and its existing completeness check. Do not use the aggregator's default tail percentile in place of the explicit 95 above.

At each segment count, the selector pools every terminating-bin score over all DMs and subbands, fits an exponential tail above the 95th percentile and applies the nominal alpha=0.1 prescription. It records every above-cut interior local maximum. Peaks within 0.035 Hz are grouped at a fixed segment count, retaining their source memberships; alternatives at different segment counts remain separate. Expected totals are 4,397 raw peaks and 415 groups. Raw accumulated scores must not be used to choose between different segment counts.

Catalogue association allows integer harmonics through eight and rational multiples with numerator and denominator at most three, within 0.15 Hz and 0.3 pc cm^-3. The lowest harmonic order wins, followed by frequency and DM residuals. These are inspection leads, not automatic identifications. The expected automatic labels comprise 27 fundamental and 51 harmonic/rational associations. The secure original recoveries are D, F, I, J and O.

RFI labels are warnings and retain every candidate: 990--1010 Hz for the crowded region, within 0.05 Hz of multiples of 25 Hz, and within 0.4 Hz of 216.68, 433.357, 650.03 or 866.71 Hz. The original groups include 254, 27 and 47 flags respectively. A DM coincidence or absence of a warning does not establish an astrophysical signal.

## 4. Direct PRESTO comparison and search costs

Use `production/benchmark/run_paired_production.sbatch` with `SEARCH_PACKAGE` set to the `production` directory. Its `benchmark/controlled_benchmark.py` imports shared I/O and timing functions from `presto_pilot.py`; the latter is required even though this is not a new pilot study. `benchmark/hmm_sources` holds the exact source copies used for the stage measurements.

The production launcher runs two HMM measurements at DM 24.40 for each segment count, then runs the complete PRESTO comparison in the SIF image. The Python entry points are:

```bash
python "$PROD/benchmark/controlled_benchmark.py" hmm --root /absolute/path/to/new/paired_run --repeats 2
python "$PROD/benchmark/controlled_benchmark.py" presto --root /absolute/path/to/new/paired_run --all-dms
```

Run the second command in the PRESTO image and both inside a real Slurm allocation; the driver checks for one. Its embedded data root must resolve through the configured bind mapping. It creates the method subdirectories exclusively, so use a new run root. The production launcher also calls the retained `verify_paired.py`; this is an existing runtime dependency, not a new validation program added for this release.

At each DM the driver extracts the full observation, eight contiguous sections and sixteen contiguous sections: `21*(1+8+16)=525` searches. For valid sample count `V`, section `i` of `p` covers `[V*i//p, V*(i+1)//p)`. It subtracts the section mean and adds zero padding to the stored FFT length `N_stored//p`. It advances the section `.inf` epoch by its first valid sample, using Decimal arithmetic. A 30- or 15-minute section therefore has an approximate duration; boundaries use actual valid samples, not a nominal 1,800 or 900 seconds.

For every section the scientific commands are:

```bash
realfft -fwd section.dat
rednoise section.fft
accelsearch -ncpus 1 -numharm 1 -zmax 200 -wmax 0 -sigma 2 \
  -flo 50 -fhi 1059 section_red.fft
```

The driver retains candidate text tables, binary `.cand` files, section metadata, stage logs, measurements, `candidates_all.csv` and summaries. It deletes only its transient section `.dat` and Fourier arrays after completing each section. Preserve the binary candidates: their full-precision `r,z,w` parameters are required for folding; rounded text columns are insufficient.

After the search, classify `candidates_all.csv` with the same catalogue and association tolerances and apply advisory RFI flags, producing `classified_all.csv` and `rfi_flags.csv` in the PRESTO run directory. There are 5,538 retained PRESTO candidates in the saved run. `followup/export_presto.py --production "$PROD" --out "$FOLLOW/presto_parameters_fullprecision.json"` exports the native binary parameters. Its retained run name is `paired_19348703`; a fresh replay must map its new run directory to that expected name, or invoke the export operations with the appropriate root in a personal launcher. This path is an execution record, not a claim that a new Slurm job will have that number.

The export follows PRESTO's phase convention: `z0=z-w/2`, `r0=r-z0/2-w/6`, then `f_start=r0/T`, `fdot_start=z0/T^2`, `fddot=w/T^3`, where `T` is the padded FFT duration. There is no jerk search here (`wmax=0`). The same `zmax` gives different physical acceleration coverage at different frequencies and interval lengths.

The complete HMM grid consumes approximately 174.6 CPU-hours; all 525 PRESTO intervals consume approximately 30.5 CPU-minutes. These totals exclude common dedispersion, candidate follow-up and folding. The HMM per-task `_timing.json` records and PRESTO `presto_timings.csv` distinguish CPU time from elapsed time. Do not substitute elapsed allocation duration or allocated CPU-hours for measured process CPU time. `production/finalization/compare_results.py` combines the retained production accounting and measurements; its generated inputs are not shipped in Git.

For the internal-stage table, use the DM24.40, 16-segment HMM measurements, one 10-Hz band, two repeats, and PRESTO's full-observation section over 50--1059 Hz. On one AMD EPYC 7543 core, preparation/tracking were 29.26--29.38/0.74--0.77 s; PRESTO FFT/whitening/acceleration were 3.21/3.55/16.86 s. These cover different bandwidths and are not a direct equal-band timing comparison. The leading-work table uses the manuscript's stated `n log2(n)` FFT estimates and sparse three-predecessor state-update counts, not another timing experiment.

## 5. Recover selected tracks and prepare folding models

The selected local endpoint is not always the global maximum in its source subband. Preserve `source_curve`, `source_params` and `source_bin` from the selection. A narrow re-search centred on a candidate can select another path and is not a substitute for the original blind candidate.

`followup/prepare_remote.py` prepares `track_manifest.json` from the selected, classified and flagged tables. It reads `production/submitted_production.json` to locate the original search and selection. For a new run, a local run record needs `candidate_output` and `hmm_output` pointing to the newly generated directories. This preparation script also submits recovery jobs; use it only after configuring its scheduler paths and resources. The actual scientific recovery is:

```bash
python "$FOLLOW/recover_tracks.py" recover --manifest "$FOLLOW/track_manifest.json" --index SOURCE_INDEX
```

Repeat over the manifest's source groups (166 in the retained run). It reruns only the relevant source configurations to save their accumulated-score arrays, checks agreement with the selected saved scores, then backtraces each selected endpoint using the original transition rule and tie-breaking. Outputs include each original track, `candidate.json`, the selected curve and the source-state cache. `common_fold.py` fits the fifth-order polynomial and full circular frequency model to this original track and writes `models.json` and `track_with_fits.dat`. These fits supply the candidate models; no catalogue timing parameters enter them.

The preserved preparation sequence uses `prepare_folds.py`, `common_fold.py`, and the final summary serializer in `restored_all_415/code/restore_all.py`. `prepare_folds.py` creates `fold_manifest.json` and also launches inspection pilots, so inspect its local submission settings before invoking it. Its `.dat` inspection folds provide protected intermediate files expected by later replay drivers, **not the statistics in the current comparison table**. The summary serializer's `summary(...)` function writes the saved polynomial and circular parameters at their original barycentric reference epoch. `add_kepler_taylor_block.py` appends all six analytically derived circular-model Taylor coefficients; its candidate index and placement records are generated bookkeeping inputs and must refer to the current selection. Do not use its historical archive-restoration machinery to change candidate membership.

The native HMM preparation entry point is `restored_all_415/code/native_hmm_prepfold_v7.py --prepare`. It expects the 415 candidate summaries with both `poly_f0..f5` and `kepler_f0..f5`, their original `.inf` files and protected intermediate products. It writes `restored_all_415/native_presto_hmm_v1/manifest.json`. Preparing that manifest does not require running its superseded cleaned-filterbank native folds.

For PRESTO, `native_presto_results_v1/code/native_presto_results_v2.py prepare` reads the exported binary parameters and preceding inspection products, converts barycentric interval boundaries onto the filterbank time axis, and writes its native model manifest. Again, only its model preparation is needed before the final original-filterbank sequence below.

These preservation-oriented preparation scripts have more input requirements than the underlying numerical fits. If starting from the saved scientific-input bundle, keep its summaries, candidate identities and native model manifests together and configure the recorded paths. Do not fabricate completion markers or substitute different tracks merely to pass their preservation checks.

## 6. Catalogue-informed narrow searches and supplementary selection

### The earlier narrow follow-up grid

The retained `redetect_<letter>/run_viterbi_target.sh` scripts search two-Hz windows centred on catalogue fundamental frequencies, across the following grids. C also has `run_targeted_array.sh`. Their different valid-sample/duration handling is preserved; use the actual retained target script and saved model epochs rather than replacing them all with the later exact-segment mode.

| Target | DMs (step 0.10) | Requested segment counts |
|---|---|---|
| C | 24.30--24.80 | 2, 4, 8, 16, 32 |
| E | 24.00--24.50 | 12, 16, 24, 32, 48, 64, 96, 128 |
| H | 24.10--24.60 | 12, 16, 24, 32, 48, 64, 96 |
| L | 24.20--24.70 | 2, 4, 8, 16, 32 |
| M | 24.20--24.70 | 2, 4, 8, 16, 32 |
| R | 24.20--24.50 | 52, 58, 64, 80, 96, 128 |
| T | 24.20--24.70 | 13, 16, 24, 32, 48, 64, 96 |
| Y | 24.20--24.70 | 36, 48, 64, 96, 128, 192 |

Use `F0_hz` for the named pulsar in the retained `merged_identity_reference.json`, and pass `TARGET`, `F0_HZ`, `NT`, `DMS` and `BW_HZ=2` to the corresponding script. Each array index selects a DM. Its `--show` mode displays those settings without submitting anything. Outputs are the likelihood, track, fits and diagnostics under `redetect_<letter>/viterbi_out/DMxx.xx_NtN/`.

`targeted_native_v1/code/targeted_workflow_v2.py prepare` reads those saved products, groups likelihoods by target and segment count, and selects the strongest saved DM at each setting. It uses the original HMM native manifest's filterbank metadata and settings and the retained catalogue identity reference. Run this before adding the later J/O directories: it scans `redetect_*` and was written for these eight named targets. Do not broaden its glob over later follow-up directories without adapting a personal launcher. It writes `manifest.json`, `selected_runs.csv` and diagnostic thresholds. These are targeted follow-up settings, not new blind selection criteria.

### Supplementary 16-segment selection

This is a reselection of the **existing** blind scores. `lower_threshold_audit_v1/code/audit.py` and `supplement.py` are retained because they prepare inputs used by the reported supplementary selection. Their directory name does not make their exploratory reports part of the paper's final results or this Git update.

The existing dependency sequence is:

```bash
python "$FOLLOW/lower_threshold_audit_v1/code/audit.py" initialize
python "$FOLLOW/lower_threshold_audit_v1/code/audit.py" worker --nt 16
# Repeat worker --nt N for N=32,64,128,256,512,1024 in suitable allocations.
python "$FOLLOW/lower_threshold_audit_v1/code/audit.py" report
python "$FOLLOW/lower_threshold_audit_v1/code/supplement.py"
python "$FOLLOW/supplementary_Nt16_v1/code/workflow.py" prepare
```

`initialize` reads the original grid, candidate tables, known-source identity reference and earlier targeted native manifest. Each worker reads saved endpoint curves, not the observations. The supplementary workflow consumes `inputs.json`, the Nt16 curve inventory and peak arrays, and `Nt16_target_inspection.csv`. It applies `L>58` uniformly over all DMs and subbands at 16 segments, with the same 0.035-Hz grouping. The expected total is 4,021 raw peaks and 126 groups: 38 original and **88 additional** groups. Catalogue controls affect review priority and candidate ordering, not admission to that list.

Run `workflow.py recover --index SOURCE_INDEX` for its source groups, then `workflow.py manifest`. The retained `recover_helper.py` backtraces the original endpoint; `native_helper.py` supplies the same time-transfer/folding functions. The resulting `supplementary_candidates.csv`, `supplementary_members.csv`, `selection.json`, exact tracks and six-coefficient models feed the final native folds. Do not replace the supplementary cut with a target-centred search.

C's fundamental score is 59.94, E's 83.73 and L's third-harmonic score 58.76, all below the original 16-segment threshold of 84.95. C is a secure supplementary recovery; E and L remain provisional. For the current table, `L_cut=58` at 16 segments and the original per-segment-count threshold elsewhere. The quoted margin must belong to the same original or supplementary candidate as the fold outside parentheses.

## 7. Final native folds on the original filterbank

Use `psr_miner_search/paper_results` to locate the retained products on the authors' filesystem. It is a view of existing outputs; most of its links and tables are not shipped in Git. Superseded runs under `old_runs` are not the current results.

The final candidate engine is `47Tuc_original_filterbank_candidate_folds_20261104/code/candidate_folds_v3.py`, with its own `native_helper.py`. It reads the four prepared native model manifests for blind HMM, PRESTO, earlier targeted HMM and supplementary HMM. `prepare` constructs new `row_batches`, output records and the original-input manifest; `fold --index I` executes one row; `report` combines the records.

```bash
export CAND_RAW="$EXP/psr_miner_search/47Tuc_original_filterbank_candidate_folds_20261104"
python "$CAND_RAW/code/candidate_folds_v3.py" prepare
python "$CAND_RAW/code/candidate_folds_v3.py" fold --index I
python "$CAND_RAW/code/candidate_folds_v3.py" report
```

Replace `I` with every required manifest row; `worker_v3.sbatch` is the retained one-CPU worker template. Prepare only in a fresh output tree, after the earlier model manifests exist. The final preparation stores raw-input metadata and recalculates PRESTO start/end fractions against the **original** filterbank duration. Reusing cleaned-filterbank fractions would fold different intervals.

For HMM, TEMPO transfers each saved barycentric six-coefficient model onto the topocentric original data. It checks the barycentric polynomial transfer separately, then generates MeerKAT polycos at the observing frequency and candidate DM. It preserves the reference epoch and all coefficients through `f5`. For PRESTO, it uses the binary-candidate phase model and corresponding section epoch and boundaries. Invalid or non-positive saved Taylor models are recorded as invalid; they are not silently truncated or replaced.

The final candidate `prepfold` command uses `-filterbank -topo -noxwin -ncpus 1 -nsub 32 -npart 64`, the transferred polycos, correct start/end fractions, and `-npfact 2 -ndmfact 3 -searchfdd`. Higher HMM coefficients remain fixed during refinement of DM and residual frequency terms through `f2`. Profile length is PRESTO's native default for the sampling time and model frequency. Outputs include `.pfd`, `.bestprof`, PostScript/PDF/PNG diagnostics, polycos, timing-model inputs, logs and completion records, with `presto_original_filterbank_` prefixes. Read sigma from the native best-profile statistic; it is not peak/off-pulse S/N.

### Independent reference column

`47Tuc_original_filterbank_best_ephemeris_folds_20261104/ephemerides_original` contains the 42 tested timing models for 23 pulsars, including alternatives. `source_manifest.json` identifies their origin, values and hashes. The public originals are linked from [Paulo Freire's 47 Tuc timing page](https://www3.mpifr-bonn.mpg.de/staff/pfreire/47Tuc/); the supplied observation-epoch O model is explicitly distinguished. The retained files, rather than whichever model is newest when this is read, define this comparison.

The final program is `code/fold_raw.py`, importing `fold_engine.py` and `verified_native_helper.py`:

```bash
export REF_RAW="$EXP/psr_miner_search/47Tuc_original_filterbank_best_ephemeris_folds_20261104"
python "$REF_RAW/code/fold_raw.py" prepare
python "$REF_RAW/code/fold_raw.py" fold --index I
```

Run all 42 indices in suitable allocations. `prepare` uses the retained source manifest and original timing files; the older inventory/download programs are not necessary to recreate that fixed input set. It changes only the internal PSR label for TEMPO compatibility. Reference folds cover the full original observation, retain the independent timing model, and use the same 64 subintegrations, 32 subbands and standard local PRESTO refinement. The program compares 12- and 15-coefficient TEMPO representations before folding. `paper_table_reference.py` reads saved results for the reference column; use the resulting fold diagnostics as well as the nominal statistic when reviewing an association.

For the table's 14 pulsars, the selected reference models are `freire_earlier_00` for D, J, G, Q, M and L; `paper_O_obsepoch_BT` for O; and `freire_latest` for C, E, F, I, H, N and U. Here `latest` is the saved manifest's label. It is not an instruction to fetch a replacement.

## 8. Later J and O complete-orbit follow-up

These experiments use the cleaned barycentric search series and the **original** filterbank for every native fold. They remain separate from the original topocentric O experiment.

### J

`47Tuc_blind_search_v1/redetect_J/code/workflow.py` prepares a grid at DMs 24.50, 24.60 and 24.70 and segment counts `32,40,48,52,56,60,62,63,64,65,66,68,72,76,80,88,96,112,128`: 57 searches. The four-Hz band is 474.04685844061--478.04685844061 Hz. It derives equal segment lengths from valid samples. Its modes are `prepare`, `search --index I`, `fold-prepare`, `fold --index I` and `status`; folding indices cover both polynomial and circular-Taylor alternatives. Its worker, measurement wrapper, helper and search-source copies are retained. Use the scalar execution modes in allocations; the `pilot`, `production` and `fold-submit` modes also submit jobs and require local scheduler configuration.

`redetect_J/full_orbit_v1/code/full_orbit.py prepare` selects the saved DM24.50, 40-segment circular fit. `build` integrates its full periodic frequency curve and constructs short polycos; `standard_refinement` performs the reported native fold. The full periodic model gives 62.4 nominal sigma in the current table. Its orbital parameters come from the HMM track, without a published J orbital model.

### O

`redetect_O_J_matched_20261105/code/workflow.py` uses DMs 24.30, 24.40 and 24.50 and segment counts `16,20,24,26,28,30,31,32,33,34,36,38,40,44,48,56,64`: 51 searches and 102 polynomial/Taylor folds. Its four-Hz window is centred on the original blind O candidate's saved ending frequency. Use the same scalar preparation/search/fold modes as J. `selected_candidate.json` fixes the reviewed choice at DM24.30 and 36 segments; it is a selection input, not a catalogue ephemeris.

For the selected track, run the retained complete-model sequence:

```bash
export OFULL="$BLIND/redetect_O_J_matched_20261105/full_orbit_v1"
python "$OFULL/code/full_orbit.py" prepare
python "$OFULL/code/full_orbit.py" build
python "$OFULL/code/full_orbit.py" fixed
python "$OFULL/code/full_orbit.py" standard_refinement
python "$OFULL/code/full_orbit.py" frequency_refinement
```

The model integrates `f_spin-K*sin[2*pi*(t-T0)/P_b]` analytically in barycentric time. A synthetic 1,024-Hz clock provides the TEMPO time mapping onto the raw sample times. Ten-minute overlapping polyco blocks encode the full sinusoid; an independent shorter representation checks interpolation and phase transfer. The orbital amplitude, period and phase remain those fitted to the track. Standard refinement changes DM and the residual terms through `f2`, not the orbit.

Expected O statistics are 22.3 for this track's Keplerian Taylor fold, 20.8 for the fixed complete orbit, 38.4 for standard refinement and 26.4 for the wider frequency/DM-only variant. `code/spin_prepare.py` then derives the measured spin correction from the latter PFD, preserving all orbital terms, and prepares `spin_refined/manifest.json`. Run `spin_refined/code/full_orbit.py build` and its `frequency_refinement` mode for the corrected-spin raw refold, which gives 27.0. The strongest separate polynomial scan fold is 30.0 at DM24.40 and 30 segments; it is a different track and is not the starting point for the 22.3-to-38.4 comparison.

These additional folds do not replace the original blind O likelihood or establish another blind recovery. Some residual phase drift remains. The earlier topocentric O polynomial result of 26.1 belongs to the first experiment, not this scan. Unreported pulse-alignment refinements are outside this guide's reproduction sequence.

## 9. Mapping the current paper to outputs

Use manuscript labels and filenames as well as figure numbers, which can change during revision.

| Manuscript item | Source and expected output |
|---|---|
| Main O spectrogram and likelihood figures (`fig:track`, `fig:loglike`) | Original topocentric Nt32 track and likelihood products; paper filenames `47TucO_spectrogram.pdf` and `47TucO_loglikes.pdf` |
| Fold-parameter table (`tab:foldparams`) | Nt32 `_fit_summary.dat`; both local frequency expansions retain `f0..f5`; convert the circular epoch convention explicitly |
| Four O fold panels (`fig:folds`) | Independent observation-epoch reference, full-coefficient polynomial/Taylor polycos, and full circular-orbit fold: 59.4, 26.1, 11.9, 11.3 |
| Coherence table (`tab:ntscan`) | Fifteen original duration-based searches and polynomial polyco folds; use actual frequency-bin counts and native PRESTO sigma |
| Scan-track panels (`fig:tracks`) | Original Nt16, Nt28 and Nt1024 spectrograms |
| Scan-likelihood panels (`fig:loglikes`) | Original Nt28 and Nt256 likelihood curves; scores and noise reference levels are timescale-dependent |
| Blind comparison (`tab:blind_search_results`) | Reviewed candidates listed below; original-filterbank native folds and separately identified targeted values |
| Leading computational work (`tab:computational_work`) | Manuscript's FFT-work expressions and sparse state-update count for one DM, one HMM timescale and one full PRESTO interval |
| Internal stage timings (`tab:measured_cost`) | Paired DM24.40 measurements: Nt16 HMM repeat records and full-interval PRESTO stage records |

The blind HMM identities outside parentheses are:

| Pulsar | Selected candidate | Applied-cut margin | Nominal fold sigma |
|---|---|---:|---:|
| D | `cand_000007_DM24.80_Nt32_f186.651067Hz_sb185.000_bary` | 827.49 | 36.9 |
| J | `cand_000073_DM24.60_Nt64_f475.955529Hz_sb473.000_bary` | 706.13 | 17.4 |
| O | `cand_000052_DM24.40_Nt32_f378.304435Hz_sb374.000_bary` | 363.30 | 21.4 |
| F | `cand_000055_DM24.40_Nt16_f381.157806Hz_sb374.000_bary` | 233.49 | 35.5 |
| I | `cand_000037_DM24.40_Nt16_f286.944418Hz_sb284.000_bary` | 42.45 | 16.7 |
| C | `supp16_002_DM24.50_f173.707800Hz` | 1.94 | 7.0 |
| E | `supp16_000_DM24.30_f282.771147Hz` | 25.73 | 6.7 |
| L | `supp16_003_DM24.40_f690.263285Hz` (third harmonic) | 0.76 | 5.6 |

All these HMM folds cover the full observation. Within the reviewed alternatives for that selected candidate and harmonic, use the stronger valid model fold. Parenthesised targeted values are J=62.4, O=38.4, C=15.0 and L=10.6. C and L use fundamental-frequency targeted searches at two segments. Their supplementary candidates use 16 segments, and supplementary L is at the third harmonic. None of these parenthesised values inherits a likelihood margin from its targeted track.

The PRESTO identities and folded intervals are:

| Pulsar | Selected candidate | Interval | Nominal fold sigma |
|---|---|---|---:|
| D | `presto_DM24.60_parts1_section00_cand000002` | Full | 96.2 |
| J | `presto_DM24.50_parts8_section00_cand000002` | Approximately 30 min | 46.4 |
| O | `presto_DM24.40_parts8_section00_cand000003` | Approximately 30 min | 49.9 |
| F | `presto_DM24.30_parts1_section00_cand000005` | Full | 39.3 |
| I | `presto_DM24.20_parts8_section01_cand000014` | Approximately 30 min, fundamental | 6.8 |
| C | `presto_DM24.30_parts1_section00_cand000041` | Full | 15.2 |
| E | `presto_DM24.50_parts8_section00_cand000020` | Approximately 30 min | 8.1 |
| L | `presto_DM24.50_parts1_section00_cand000063` | Full, third harmonic | 8.6 |
| H | `presto_DM24.30_parts16_section01_cand000012` | Approximately 15 min, second harmonic | 5.7 |

PRESTO I's reviewed second-harmonic alternative reaches 8.2, but the main table uses its fundamental. Reference sigmas for D,J,O,F,I,C,E,L,G,Q,M,N,H,U are respectively 97.0,69.2,59.4,39.7,25.5,15.0,12.8,10.9,5.6,5.5,5.4,5.3,5.2,5.2. The automatic candidate labels and largest printed fold statistics alone do not determine the paper's outcomes. HMM H/M/R/T/Y targeted products do not establish additional secure recoveries; in particular, the strongest historical H/M tracks cannot automatically be assigned to those pulsars.

**The retained `build_paper_table_original_data_with_targeted_O.py` predates the current manuscript's J/O entries and table layout. It is not an up-to-date table generator.** The mappings above follow the current manuscript. This documentation release leaves that script unchanged and excludes it from the newly added reproduction sources. Assemble the current table from the reviewed identities and final fold records, keeping the main blind values and parenthesised follow-up separate.

The appendix's feature near 938.45 Hz is a retained historical illustration, numbered 193 in the old plots and corresponding to corrected group `cand_000146` (DM23.60, 16 segments, endpoint 938.454463 Hz). Its panels do not generate the corrected candidate totals or comparison table. No historical inspection pipeline or archive is required for the current selection. Obtain those illustrative panels separately if rebuilding the complete manuscript figure package. The phase-error appendices give analytic/idealised calculations, not another measured search or a calibration of the reported fold losses.

## Reproduction limits to retain in any report

An exact numerical replay requires the original observation, retained cleaned samples or their random-cleaning provenance, the original targeted PulsarX launcher/plan or prepared series, and the historical software/timing environment. The staged sources document the scientific steps and expose their saved-input requirements; they do not establish bit-for-bit reproducibility from a fresh clone alone. No analyses were rerun to prepare this guide. Keep the two time coordinates, original and supplementary selections, targeted follow-up, and original-data final folds explicit when reporting a reproduction.
