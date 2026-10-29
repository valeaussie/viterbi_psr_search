# Viterbi HMM Pulsar Search

## Overview

This repository implements a semi-coherent blind search for radio pulsars in compact binary systems, based on a hidden Markov model (HMM) solved with the Viterbi algorithm.

The method models the time-varying apparent spin frequency of a pulsar as a hidden Markov chain. The observation for each coherent time segment is the normalised Fourier power (Schuster periodogram) of the dedispersed, heterodyned, and downsampled radio intensity time series. The Viterbi algorithm finds the most probable frequency track through the resulting time-frequency spectrogram and returns a log-likelihood statistic used to assess significance.

The approach is designed to handle the Doppler modulation from a short-period binary companion, where the apparent spin frequency drifts across many Fourier bins over the total observation, making standard FFT-based searches insensitive. Because the transition model is an unbiased random walk (equal probability of moving by -1, 0, or +1 frequency bin per coherent segment), the pipeline does not search over orbital parameters explicitly. The coherent segment length is instead chosen short enough that the frequency drift per segment is at most one bin, for the widest binary orbits of interest.

The method is described in:

> [O'Leary, Dunn & Melatos (2026)](https://iopscience.iop.org/article/10.3847/1538-4357/ae3288)

---

## Candidate diagnostic plots

The following PDFs contain the 415 selected HMM candidates from the saved
`47Tuc_blind_search_v1/stage4_fold/comparison_v1/hmm/` comparison. Each PDF
has one page per candidate, in the same order, with bookmarks for navigation.
These are snapshots restored from the saved search arrays and model fits;
adding them here does not rerun the search or update the fits. They are
separate from the corrected `blind_v2` production rerun.

| Plot collection | Contents |
|---|---|
| [Spectrograms with fits](docs/plots/hmm_candidate_review/01_spectrograms_with_fits.pdf) | Time-frequency maps with the Viterbi track, polynomial fit and circular Kepler fit. |
| [Likelihoods with thresholds](docs/plots/hmm_candidate_review/02_likelihoods_with_thresholds.pdf) | Full source-subband ending-state score curves, selection thresholds, selected endpoints and local peak insets. |
| [Plain likelihoods over 2 Hz](docs/plots/hmm_candidate_review/03_plain_likelihoods_2Hz.pdf) | Score curves within 1 Hz of each selected endpoint, clipped to the original subband, without thresholds, markers or insets. |

The model overlays are diagnostic fits. In particular, a circular Kepler
curve is not a validated orbital solution.

### RFI flags in the plot titles

The RFI flags are automatic frequency-based warnings, not confirmed RFI
classifications. The thresholds below use the candidate's selected
`peak_freq_hz`, which is the final Viterbi-track frequency in these plots.
These are barycentric search frequencies in Hz, not receiver subband
frequencies in MHz. The rules do not assess the fold or its DM dependence.

| Flag | What triggers it | Why it was flagged |
|---|---|---|
| `rfi_swarm` | Frequency between 990 and 1010 Hz, inclusive. | A crowded region around 1 kHz where interference generated many similar candidates. |
| `rfi_comb` | Within 0.05 Hz of a multiple of 25 Hz. | Possible association with regularly spaced interference tones, such as 75, 100 or 125 Hz. |
| `rfi_family` | Within 0.4 Hz of 216.68, 433.357, 650.03 or 866.71 Hz. | Near a previously identified interference family and its harmonics. |
| `none` | None of these frequency rules matched. | No warning from these rules; this does not establish that the candidate is free of RFI. |

The rules are applied in the order swarm, comb, then family, and the first
match supplies the label. **Retained** means the candidate was kept for
inspection and folding despite its flag. Every selected candidate was
retained, including flagged candidates and known pulsars. A flag alone is
not a reason to reject a candidate; DM-zero comparisons and folding provide
additional evidence.

The thresholds and retention policy are defined in
[`flag_rfi.py`](code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/flag_rfi.py).

For this plotted candidate set, RFI flags are advisory and all candidates
are retained. References later in this README to moving folders or excluding
flagged candidates describe the earlier Stage 3c policy.

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
                stage1a_plots_inspect_dynamic_spectrum.py
                stage1_rfi_cleaning.sh
                stage1b_bandpass_diagnostic.sh (calling stage1_clean_postprocess.py)
                stage2_dedisp_new.sh
                stage3_blind.sh
                stage3a_viterbi_for_47Tuco.sh
                stage3b_fit.sh
                stage4_fold.sh
                stage5_collate.sh
                aggregate_blind.py      # Stage 3a candidate aggregation
                classify_known.py       # Cross-match candidates against known pulsars
                flag_rfi.py             # Stage 3c RFI separation
                get_known_label.py      # Format one classified row as a plot legend label
                parse_dedup_row.py      # Read one candidate row for the Stage 3b launcher
                collate_fold_snr.py     # Stage 5 S/N extraction
                make_candidate_gallery_base64.py
                gen_search_params.py
                patch_known_pulsars.py
            data/                   # Input filterbank, params files
            config/                 # observation.yaml, known_47tuc_pulsars.yaml
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

Check whether the data contains RFI, then clean it.

To inspect the dynamic spectrum:

```bash
python stage1a_plots_inspect_dynamic_spectrum.py \
    --fil <exp-dir>/data/47Tuc_22UL_1of2_L.fil \
    --flo 880 \
    --fhi 970 \
    --nsamples 10000 \
    --outdir <exp-dir>/stage1_clean
```

This produces `dynspec_rfi_ratio.png` and `dynspec_stats.txt` in the output directory. The `--flo` and `--fhi` values above cover one satellite band; adjust them to the RFI environment of your own observation.

**Hard zapping is not used in this pipeline.** For the 47 Tuc MeerKAT data, the FBFUSE beamformer pre-normalises the filterbank, and per-channel mean and baseline ratios across the known satellite RFI bands were measured to be within 0.3 per cent of unity. There is therefore no persistent narrowband RFI to remove by hand, and hard-zapping named frequency ranges only discards clean data. Run the ratio test as a diagnostic, not as a trigger for zapping.

Cleaning is performed by `filtool` (PulsarX) using statistical filters only:

```bash
sbatch stage1_rfi_cleaning.sh
```

The RFI configuration is a single line inside that script:

```
RFI_FLAGS="-z zdot kadaneF 8 4 kadaneT 8 4 --fillPatch rand"
```

`zdot` removes the correlated zero-DM component, and `kadaneF` and `kadaneT` apply the Kadane maximum-subarray filter in frequency and in time respectively. No frequency ranges are zapped. The cleaned filterbank is written to `stage1_clean/` and is the input to Stage 2.

To verify the bandpass after cleaning:

```bash
sbatch stage1b_bandpass_diagnostic.sh
```

To generate cleaning statistics:

```bash
python stage1c_clean_postprocess.py
```

This writes `rfi_stats.yaml`, recording how many channels were affected and the flags used. Regenerate it whenever the cleaning is rerun, since a stale copy will describe a filterbank that no longer exists.

### Stage 2: Dedispersion

```bash
sbatch --array=0-<N_DM-1> stage2_dedisp_new.sh
```

Runs `prepdata` (PRESTO) to produce a dedispersed and barycentred time series for each DM trial. The DM grid is defined inside `stage2_dedisp_new.sh`. Each array task handles one DM value and writes a `.dat` and `.inf` file to `stage2_dedisp/DM<XX.XX>_bary/`.

### Stage 3: Blind Viterbi search

Each task calls `viterbi_pipeline.py` in lean blind mode for one
(DM, Nt, subband) combination, writing a `*_loglike_curve.dat` and a sibling
`*.params` file to `stage3_viterbi/blind_v1/DM<XX.XX>/Nt<N>/f0_<F>/`.

The `*.params` file records `Tsft`, `out_prefix`, and other search
parameters. It is read by `aggregate_blind.py` to infer `T_obs` without any
additional inputs.

The total number of tasks is `N_DM * N_Nt * N_subbands`. For the 47 Tuc run
that is `21 * 7 * 112 = 16464`, so array indices run from 0 to 16463.

**The array will not fit in one submission.** OzSTAR sets `MaxArraySize` to
2048, so the highest permitted array index is 2047 and any submission asking
for more is rejected outright. Confirm the limit before sizing a run:

```bash
scontrol show config | grep -i MaxArraySize
```

`stage3_blind.sh` handles this by offsetting the task index rather than by
splitting the work list. Line 60 reads:

```bash
IDX=$(( SLURM_ARRAY_TASK_ID + ${IDX_OFFSET:-0} ))
```

so each submission covers a 2048-wide slice of the global index range, and
`IDX_OFFSET` selects which slice. The full 16,464 tasks therefore take nine
submissions: eight complete chunks plus a remainder of 80.

```bash
for k in $(seq 0 7); do
    sbatch --array=0-2047 --export=ALL,IDX_OFFSET=$((k * 2048)) stage3_blind.sh
done
sbatch --array=0-79 --export=ALL,IDX_OFFSET=16384 stage3_blind.sh
```

**`--export=ALL,IDX_OFFSET=N` is mandatory, and the `ALL` is not optional.**
Writing `--export=IDX_OFFSET=N` alone sets `SBATCH_EXPORT` to that single
variable and suppresses the rest of the submitting environment. The job then
starts without the module and path setup it depends on. Worse, `IDX_OFFSET`
itself defaults to zero if it fails to arrive, so a broken submission does
not error: it silently reruns chunk zero and overwrites its output. Every
chunk appears to succeed while eight ninths of the search never happens.

After the last chunk finishes, verify the count before aggregating:

```bash
find <exp-dir>/stage3_viterbi/blind_v1 -name '*_loglike_curve.dat' | wc -l
```

The answer must equal `N_DM * N_Nt * N_subbands` exactly. Anything less
means a chunk failed or was overwritten, and `aggregate_blind.py` will
happily build a candidate list from the incomplete set without complaining.

### Stage 3a: Candidate aggregation

```bash
python aggregate_blind.py \
    --blind-dir  <exp-dir>/stage3_viterbi/blind_v1 \
    --known-yaml <exp-dir>/config/known_47tuc_pulsars.yaml \
    --out-dir    <exp-dir>/stage3_viterbi/blind_v1 \
    --false-alarm-rate 0.1 \
    --tail-percentile 95 \
    --dedup-tol-hz 0.035
```

Reads every `*_loglike_curve.dat` file, groups them by `Nt`, and computes a separate detection threshold for each `Nt` group following the exponential-tail method of O'Leary et al. (2026), Equations 13-17 and Appendix D.

**Threshold method.** The Viterbi log-likelihood scales with `Nt` (it is a cumulative sum over `Nt` coherent segments), so a single global threshold across all `Nt` values is not meaningful. Within each `Nt` group, all loglike values from all subbands and DM trials are pooled. One can choose the percentile and the tail above that percentile is fitted as an exponential using the maximum-likelihood estimator for the rate parameter (in our run we chose 95% percentile). The threshold is then set by inverting the desired false alarm probability per subband (`--false-alarm-rate`, default 0.1). The number of independent Viterbi runs used to build the noise distribution is taken to be the number of curve files in that `Nt` group (the role of `N_real` in the paper).

**Deduplication.** Candidates are deduplicated in frequency with a tolerance that, by default, is derived automatically from the coarsest bin width in the `Nt` grid: `tol = 2 * N_T_max / T_obs`. The intent is that the same pulsar detected at different `Nt` values, whose peak frequencies can differ by up to one coarse bin due to quantisation, is merged into a single candidate.
Warning: this default couples the dedup tolerance to the largest `Nt` in the grid, and it over-merges. Because the tolerance is set by `N_T_max`, extending the `Nt` grid widens the merge window for every candidate in the search, including finely resolved ones that are genuinely distinct. With `T_obs` of about `14700s`, `N_T_max = 256` gives a tolerance of `0.035 Hz`, while `N_T_max = 1024` gives `0.139 Hz`, a factor of four, purely as a side effect of adding coarse `Nt` values for the tight-orbit regime. At the wider tolerance the 47 Tuc candidate list contains 88 entries; at the narrower one it contains 194. The difference is entirely in the merging, since the loglike curves, the tail percentile and the per-`Nt` thresholds are unaffected by this setting. The 106 additional candidates are absorbed into neighbours, not rejected on significance.
Recommendation. Set the tolerance explicitly with `--dedup-tol-hz` rather than relying on the auto value, so that a later change to the `Nt` grid cannot silently move it. For the 47 Tuc run the tolerance was fixed at 0.035 Hz, which corresponds to the coarsest bin of the validated part of the grid `(N_T_max = 256)`. Detections at `Nt = 512` and `Nt = 1024` are then matched into that window rather than widening it.
Note that a fine tolerance splits a single bright pulsar into several near-F0 entries, since its per-trial peak frequencies scatter by more than the merge window. This is expected and is the safer failure direction: over-splitting duplicates information and can be collapsed afterwards, whereas over-merging destroys it and can only be recovered by re-running aggregation.

Outputs:

- `candidates_raw.csv`: one row per above-threshold peak, with DM, Nt, subband, frequency, loglike, and known-pulsar cross-match.
- `candidates_dedup.csv`: deduplicated candidate list ranked by loglike, with columns `dm_count` and `multiplicity` indicating how many DM trials and (Nt, subband) combinations recovered the candidate.
- `threshold_calibration.txt`: full calibration record (lambda, L_tail, L_th, sigma_L_th) for each Nt group.

### Stage 3b: Full fit

**Prerequisite: run `classify_known.py` first.** Stage 3b looks up each candidate's known-pulsar classification in `classified.csv` in order to annotate the plot legends. If that file does not exist, the lookup returns an empty label and **every plot is produced unannotated, with no error and no warning**. Confirm the file is present before launching:

```bash
ls <exp-dir>/stage3_viterbi/candidates/classified.csv
```

The full order is therefore: Stage 3a aggregation, then `classify_known.py`
(see [Known-pulsar classification and plot annotation](#known-pulsar-classification-and-plot-annotation)),
then Stage 3b, then Stage 3c RFI separation. Stage 3c must come last of the
four, because it moves the `cand_NNN` folders that Stage 3b creates.

Set the array upper bound to the number of candidates minus one, that is, the number of data rows in `candidates_dedup.csv` (its line count minus one for the header). For the 47 Tuc run with 194 candidates the range is `0-193`.

```bash
sbatch --array=0-<N_candidates-1> stage3b_fit.sh
```

Runs `viterbi_search_and_fit.py` on each candidate from `candidates_dedup.csv`. Performs a refined Viterbi search followed by a polynomial frequency evolution fit and optionally a Kepler orbital fit. Writes per-candidate result files to `stage3_viterbi/candidates/`.

**Note on `Nt_best` versus the actual segment count.** `Nt_best`, as recorded in
`candidates_dedup.csv` and used throughout Stage 3b, is a nominal label from the
coherent-timescale grid, not a guarantee of how many real segments a candidate's
refit will use.

`Tsft` (equivalently `T_coh`) for Stage 3b is computed as `Tobs / Nt_best`. 
Both `stage3_blind.sh` and `stage3b_fit.sh` derive `Tobs` at run time as `NSAMP * TSAMP`, 
reading both values from the matching `.inf` file. `NSAMP` is the `Number of bins in the time series` field, 
which counts PRESTO's zero-padded length, so the resulting `Tobs` is longer than the real, non-padded observation.

`viterbi_search_and_fit.py` separately trims the loaded `.dat` data down to the
real sample count recorded in the `.inf` file before building the spectrogram,
matching the trim already applied in `viterbi_pipeline.py`. As a result, the
number of full coherent segments that actually fit into the real data, `Nsft`,
shown in the Viterbi path legend as `N_seg`, can be noticeably lower than
`Nt_best`, particularly for short-`T_coh` candidates, where the fixed padded
tail (roughly 300 s for this dataset) is a large fraction of a single segment.

For example, a candidate with `Nt_best = 1024` has `T_coh = 14.36` s. At that
`T_coh`, roughly 21 of the segments implied by the padded `Tobs` fall entirely
inside the padded tail and are excluded, leaving `Nsft = 1003` real segments.
This is expected, and does not indicate a problem with the candidate or the
fit, `Nt_best` is a grid label, while `N_seg` in the plot is the actual number
of segments the refit used.


### Stage 3c: RFI separation

```bash
python flag_rfi.py \
    --classified-csv <exp-dir>/stage3_viterbi/candidates/classified.csv \
    --known-yaml     <exp-dir>/config/known_47tuc_pulsars.yaml \
    --cand-dir       <exp-dir>/stage3_viterbi/candidates \
    --out-csv        <exp-dir>/stage3_viterbi/candidates/rfi_flags.csv
```

Identifies interference among the candidates that `classify_known.py` left
`unclassified`, and moves their `cand_NNN` folders into `candidates_rfi/`.
Nothing is ever deleted, and candidates labelled `known_main` or
`known_harmonic` are excluded from consideration by construction.

Without `--move` the script only reports. Add `--move` once the dry-run
totals look right.

**Three rules**, applied in order. The first to match claims the candidate,
so the recorded reason reflects rule order, not physical priority.

| Rule | Test | Rationale |
|------|------|-----------|
| `rfi_swarm` | frequency in 990 to 1010 Hz | One broad interference structure fragmented by the dedup tolerance into many nominal candidates |
| `rfi_comb` | within 0.05 Hz of a multiple of 25 Hz | Round-number terrestrial tones |
| `rfi_family` | within 0.4 Hz of 216.68, 433.357, 650.03 or 866.71 Hz | A signal confirmed terrestrial by DM-zero folding, plus its harmonics |

**Safety check.** Before anything is reported or moved, all three rules are
applied to the candidates already labelled as known pulsars. If any rule
would catch one of them, the script aborts and changes nothing. The known
pulsars are the control group, and no rule that damages them is acceptable.
This is how the family tolerance was set: at 0.6 Hz the rule catches the 2/3
harmonic of J0024-7204ah at 217.155 Hz and the check fails, so the tolerance
was tightened to 0.4 Hz, which still covers every observed family member.

**Override.** A candidate lying within 0.1 Hz of any `n/m` harmonic of a
catalogue pulsar, for `n` and `m` up to 4, is never vetoed. It is labelled
`hold_near_known` and left in place for inspection. Veto windows are blunt,
and the swarm window alone spans 20 Hz, so without this a genuine harmonic
falling inside would vanish with no trace.

**Outputs.**

- `rfi_flags.csv`: every unclassified candidate, with verdict, reason, and
  its nearest catalogue harmonic and residual.
- `candidates_rfi/rfi_manifest.csv`: written on `--move`, so the reason for
  each removal travels with the folders rather than living only in a log.

**Result on the 47 Tuc run.** Of 151 unclassified candidates, 144 were moved
(104 swarm, 13 comb, 27 family), one was held, and six were kept. The six
comprise five apparent fragments of known pulsars that `classify_known.py`
missed on too tight a frequency tolerance, all offset by about 0.16 Hz from
their parent, and one candidate at 938.45 Hz whose nearest catalogue
harmonic is almost 4 Hz away.

**The parameters are specific to this dataset.** The swarm window is the
measured extent of one interference feature in this observation, not a
physical scale, and it must be re-derived for other data. The rules encode
what was found here, not a general RFI taxonomy.

**Two notes on the rules as written.**

The 650 Hz entry in the family list is redundant and slightly wrong. Fitted
frequencies for those candidates sit at 649.982 Hz, which is 650.000 with
the barycentric shift applied, matching the comb rather than the third
harmonic of 216.68 Hz that the family list predicts. The comb rule fires
first in any case, so the outcome is unaffected.

The strongest argument against the swarm is sharpness, not count. Genuine
comb lines are sharp: one candidate each, with fitted frequencies scattering
by well under one part in a hundred thousand. The 1 kHz feature spans
0.12 Hz across 32 fitted candidates, running continuously through the round
value with no gap. A rotating object has one frequency. That smear does not.


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
| `Nt` grid | 16, 32, 64, 128, 256, 512, 1024 | Coherent segment counts |
| DM grid | 21 trials | Centred on cluster DM |
| False alarm rate | 0.1 per subband | `--false-alarm-rate` in Stage 3a |
| Dedup tolerance | 0.035 Hz (set explicitly) | `--dedup-tol-hz`; do not rely on the auto default, see Stage 3a warning |
| Tail percentile | 95 | `--tail-percentile`; script default is 99.99 |
| Swarm window | 990 to 1010 Hz | `SWARM_LO`, `SWARM_HI` in `flag_rfi.py`; dataset-specific |
| Comb tolerance | 0.05 Hz on multiples of 25 Hz | `COMB_TOL`, `COMB_SPACING` in `flag_rfi.py` |
| Family tolerance | 0.4 Hz | `FAM_TOL` in `flag_rfi.py`; set by the control-group safety check |

---

## Application to 47 Tucanae

### Data

The pipeline has been applied to a MeerKAT observation of the globular cluster 47 Tucanae (47 Tuc), observed at L-band with the MeerKAT L-band receiver, centred at 1283.582 MHz. The filterbank records total-intensity data with 128 frequency channels across an 856 MHz bandwidth, 8-bit digitisation, and a sampling time of 76.56 microseconds. The total observation spans 14400 s (4 hours).

The data were processed with the call path described in Stages 0-5 above, using the experiment directory `47Tuc_blind_search_v1`.

### Search configuration

- DM grid: 21 trials spanning 23.50 to 25.50 pc/cm^3 in steps of 0.10 pc/cm^3, bracketing the cluster DM of 24.36 pc/cm^3.
- Nt grid: 16, 32, 64, 128, 256, 512, 1024, corresponding to coherent timescales T_coh = Tobs/Nt from approximately 14s to approximately 919s.
- Subband width: 10 Hz; total band: 50-1050 Hz; 112 subbands (9 Hz step, 1 Hz overlap).
- Total Viterbi jobs: 21 DM x 7 Nt x 112 subbands = 16,464 per polarisation.

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

After Stage 5, the ranked CSV at `stage4_fold/fold_snr_ranked.csv` is the primary output of the fold.
The `cand_NNN/` folders produced by Stage 3b correspond to every row of
`candidates_dedup.csv`, in row order. No filtering is applied at that point:
the set includes RFI, re-detections of known pulsars, and harmonics of known
pulsars, as well as any genuinely new candidates. Stage 3b is deliberately
unselective, so that nothing is discarded before it has a diagnostic plot.

Filtering happens afterwards, in Stage 3c. `flag_rfi.py` moves candidates
matching the interference rules into `candidates_rfi/` and records why in
`rfi_flags.csv`. Nothing is deleted, so a veto can be reversed by moving the
folder back. After Stage 3c, `candidates/` holds only the known pulsars, any
candidate held near a catalogue harmonic, and the genuinely unexplained.
That reduced set is what goes forward to folding.

Because the numbering follows row order in `candidates_dedup.csv`, the `cand_NNN` indices are not stable across aggregation runs. Re-running Stage 3a with different settings renumbers every candidate. Refer to candidates by frequency and parent pulsar in any writing that needs to outlive a single run.
To assess whether a candidate is a real new pulsar:

1. Check `classification` in `classified.csv` (see below). Do not rely on the known_match column of `candidates_dedup.csv`: it flags only exact frequency coincidences with a known pulsar and misses harmonics entirely. A candidate at twice a known pulsar's `F0` has an empty known_match and will otherwise pass straight through as if it were new.
2. Check `verdict` in `rfi_flags.csv`. A candidate marked `hold_near_known`
   sits inside an RFI veto window but close to a catalogue harmonic, and was
   spared for exactly that reason. It needs a decision rather than a default.
   Candidates already moved to `candidates_rfi/` carry their reason in
   `rfi_manifest.csv`.
3. Check `dm_count`. A real pulsar should peak in significance at or near the cluster DM (~24.36 pc/cm^3) and be recovered across multiple DM trials. Detections at only one DM trial are suspect.
4. Check `multiplicity`. A real pulsar should appear in multiple subbands and ideally multiple `Nt` values.
5. Inspect the folded `.ar` file. A real pulsar shows a coherent pulse profile and dispersion sweep consistent with the cluster DM.
6. Check for harmonics of other candidates. If a candidate at frequency `f` has a corresponding detection at `f/2` or `2f`, one is likely a harmonic of the other. Harmonics of known pulsars are handled by `classify_known.py` and should already be labelled. Note that harmonic relations are not restricted to integer ratios: a candidate can sit at `3/2` of another, so integer-only checking is insufficient.


### Known-pulsar classification and plot annotation

47 Tuc contains 42 known pulsars, several of them bright, so a large fraction of any candidate list consists of re-detections and their harmonics. Two helper scripts identify these automatically and label the diagnostic plots.

`classify_known.py` reads `candidates_dedup.csv` and cross-matches every candidate against `config/known_47tuc_pulsars.yaml`, writing `classified.csv` alongside the candidate list:

```bash
python classify_known.py \
    --cands      <exp-dir>/stage3_viterbi/candidates/candidates_dedup.csv \
    --known-yaml <exp-dir>/config/known_47tuc_pulsars.yaml \
    --out        <exp-dir>/stage3_viterbi/candidates/classified.csv
```

Each candidate is labelled `known_main` (the pulsar itself), `known_harmonic` (a rational multiple of a known `F0`), or `unclassified` (RFI, or genuinely new).

A match is accepted only if the frequency lies within `--tol-hz` of `(n/m) * F0` **and** the candidate DM lies within `--dm-tol` of that pulsar's DM. The DM gate is what prevents RFI from being mislabelled: terrestrial signals do not sit at a cluster pulsar's DM, so they fail the second condition regardless of how well the frequency happens to line up.

Allowed ratios are integers up to `--max-int` plus fractions `p/q` with `p,q` up to `--max-frac-pq`. Restricting the fractions matters, because the catalogue is dense enough that unrestricted ratio searching will find a spurious match for almost any frequency. Where several matches are acceptable, the lowest harmonic order wins, then the smallest frequency residual, then the smallest DM residual.

The output records the parent pulsar, the ratio, the parent `F0` and DM, the predicted frequency, both residuals, and an `n_ambiguous` count flagging candidates that more than one pulsar could explain.

On the 47 Tuc run, 194 candidates resolve to 17 `known_main` and 26 `known_harmonic` entries, drawn from five distinct pulsars (D, F, I, J and O). The harmonics include the full ladder of D at 2x, 3x, 4x and 5x, and the second harmonics of J, O and F. The candidates at 952.10 Hz and 756.61 Hz are the second harmonics of J0024-7204J and J0024-7204O respectively. Both parents are themselves binaries, which accounts for the smooth Doppler-modulated Viterbi tracks at those frequencies. A smooth sinusoidal track is therefore not, on its own, evidence of a new binary: always check the harmonic classification before treating one as a discovery.

`get_known_label.py` formats one classified row into a short caption for the plot legend. Stage 3b calls it once per candidate and passes the result to `viterbi_search_and_fit.py` through the per-candidate params file as `known_label`, which appears as the first legend entry on the spectrogram plot. Candidates classified as `unclassified` produce an empty label and their plots are left unannotated. This is intentional: RFI plots are retained as reference examples of what interference looks like in this pipeline.

Stage 3b writes both diagnostic plots as PDF rather than PNG, so that the curves, axes and labels remain vector and publication quality. The grey spectrogram background is an embedded raster at the native resolution of the array, which no dpi setting will change.