# 47 Tuc: HMM search reproducibility

This repository contains the code, small fixed inputs and three candidate
diagnostic PDFs for the targeted redetection of 47 Tuc O and the blind
HMM/Viterbi search of the MeerKAT observation. The observation itself, other
generated results, exploratory runs, software regression tests and PRESTO
benchmarking work remain on the cluster.

The search code includes the corrections made in October 2026: exact segment
counts after removing PRESTO padding, consistent FFT frequency labels, and
candidate alternatives retained separately at each timescale. Corrected runs use a new `blind_v2` output tree and frozen source copies.
Existing paper figures and candidate counts describe the historical run until
regenerated and verified.

## What is included

Paths below are relative to
`code/experiments/47_Tuc_O_Alessandro_Ridolfi/` unless stated otherwise.

| Files | Purpose |
|---|---|
| `code/viterbi.py` | Frequency-state recursion and backtracking |
| `code/viterbi_pipeline.py` | Blind search, including lean mode |
| `code/viterbi_search_and_fit.py` | Search, track fitting and diagnostic products |
| `jobs/stage0_inspect.py` and `jobs/patch_known_pulsars.py` | Observation inspection and catalogue preparation |
| `jobs/stage1*` | RFI cleaning and its diagnostics |
| `jobs/stage2_dedisp_new.sh` | Blind-search dedispersion and barycentring |
| `jobs/dm_list.txt` and `jobs/nt_list.txt` | Fixed DM and segment-count grids |
| `jobs/stage3_blind.sh`, `aggregate_blind.py`, `classify_known.py`, `stage3b_fit.sh` and their helpers | Blind search, selection, classification and follow-up |
| `jobs/stage4*`, `stage5_collate.sh` and result-inspection helpers | Folding, collation and candidate inspection |
| `47Tuc_blind_search_v1/config/` | Fixed observation metadata and known-pulsar catalogue snapshots |
| `redetection_of_47TucO/` selected scripts | Targeted search, threshold calculation, significance summaries, figures and folding |
| `redetection_of_47TucO/47TucO_psrcat.par` and `J0024-7204O_obsepoch.par` | Small reference ephemeris snapshots |

The Monte Carlo calculation in `redetection_of_47TucO/compute_Lth.py` is part
of the scientific threshold calculation and is included. Exploratory tests
and the separate synthetic software regression suite are kept locally only.

## Inputs and environment

The observation is `47Tuc_22UL_1of2_L.fil`. Raw and cleaned filterbanks,
dedispersed `.dat`/`.inf` pairs and container images must be supplied separately.
The committed YAML files are snapshots used to describe the observation and
catalogue; their absolute paths and file-existence fields describe the original
cluster environment, not a new checkout.

The cluster environment is loaded with:

```sh
source /fred/oz022/vdimarco/software/envrmnts/setup_viterbi_psr.sh
```

The scripts use Python with NumPy, SciPy, Matplotlib, PyYAML, sigpyproc and
psrqpy. External processing uses PRESTO, PulsarX and PSRCHIVE as appropriate.
Some plotting scripts use LaTeX. The environment setup and container images
are external to this repository. Record their versions and the source commit
with every production run; this repository does not yet provide a fully
pinned environment.

## Targeted redetection of 47 Tuc O

`redetection_of_47TucO/grid_search.sh` records the targeted frequency band,
DM, observation epoch and coherent-duration grid. Its input is the targeted
topocentric dedispersed time series named in that launcher, with its matching
`.inf`. This is a separate input from the cleaned, barycentred blind-search
DM grid. The raw data preparation and the actual DM recorded in the `.inf`
must be checked before rebuilding that targeted input.

`stage_Lth.sh` and `compute_Lth.py` calculate noise thresholds for the targeted
grid. `collate_Lth.py`, `norm_viterbi_significance.py`,
`plot_ntscan_summary.py` and `plot_ntscan_foldsnr.py` produce the associated
summaries and figures. `fold_on_track.py`, `poly_grid_fold.sh` and
`kepler_fit_to_par.py` provide the track/model folding and parameter conversion.
Some historical summaries contain fixed measured values, so they are not a
substitute for regenerating measurements after a new search.

## Blind search

The scripts in `jobs/` record inspection, cleaning, dedispersion, the blind
search, follow-up and folding. The committed grids contain 21 DMs from 23.50
to 25.50 in steps of 0.10, and segment counts
16, 32, 64, 128, 256, 512 and 1024.

With the corrected `segments` option, the Python programs remove padding
using the matching `.inf`, then divide valid samples into the requested
number of equal integer-length segments. They save the actual duration and
sample counts in `_sampling.json`. The historical `Tsft`/`Nsft` mode remains
available for reproducing older settings.

`aggregate_blind.py` writes `candidates_raw.csv` and
`candidates_timescale_alternatives.csv`. It merges nearby DM/subband peaks
within a timescale and retains alternatives across timescales. `Nt_best`
identifies a row's timescale; it is not a comparison of raw scores between
timescales. The threshold formula is unchanged, and its nominal false alarm
probabilities in correlated observation data still need validation.

Every selected candidate has an informative `candidate_id`, for example
`cand_000123_DM24.40_Nt32_f378.307411Hz_sb374.000_bary`. The ordinal is
fixed within the frozen candidate manifest; it is not a cross-run identifier.
`candidate_members.csv` maps every above-threshold local peak to its group,
including source curve, frequency bin and parameters. The representative's
DM/Nt/subband are retained. Nearby DM/subband peaks may still be grouped
within Nt; the complete raw table and group membership remain available.
Local-peak separation is one bin, avoiding the historical five-bin suppression.
The outermost bin of each curve is not treated as an interior local maximum;
the overlapping subbands cover interior band boundaries.

`make_search_grid.py` freezes the expected grid. Its `--validate-results`
mode checks every timing and sampling record before aggregation.
`aggregate_blind.py` requires `--expected-curves`, `--expected-grid` and
`--run-id`, fails on incomplete/unreadable curves or failed thresholds, and
refuses to overwrite a frozen candidate manifest. The observation rerun uses
`--tail-percentile 95 --false-alarm-rate 0.1 --dedup-tol-hz 0.035`, retaining
the historical threshold prescription without claiming a calibrated false
alarm probability for correlated trials.

Classify the new table with `classify_known.py` before follow-up. It preserves
all candidate IDs, provenance and frequency precision. Its classification
file must correspond to that exact table.
`stage3b_fit.sh` uses `classified_timescale_alternatives.csv` in the
`stage3_viterbi/candidates/` folder. Do not reuse the old 194-row classification or array
bounds. Set submission array sizes from the current grid/table and respect
the cluster's array-size limit.

## Candidate diagnostic plots

The following PDFs contain the 415 selected HMM candidates from the saved
`47Tuc_blind_search_v1/stage4_fold/comparison_v1/hmm/` comparison. Each PDF
has one page per candidate, in the same order, with bookmarks for navigation.
These are snapshots restored from the saved search arrays and model fits;
adding them here does not rerun the search or update the fits. They are
separate from the corrected `blind_v2` production rerun described above.

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

## Before a new production run

No job is submitted by preparing or committing this repository. Review the
site-specific paths and the scientific run configuration before submission:

- Select new output directories rather than overwriting historical results.
- Check that targeted inputs are topocentric and blind-search inputs are the
  intended cleaned, barycentred DM series.
- Set the same corrected candidate directory and table in the later folding,
  collation and inspection steps. The defaults now use
  `stage3_viterbi/candidates/` and its frozen table throughout.
- `stage1b_bandpass_diagnostic.sh` currently names
  `stage1_clean_postprocess.py`; the included implementation is
  `stage1c_clean_postprocess.py`. Update that path before using the wrapper,
  or invoke the included Python script directly with its input arguments.
- Validate representative corrected searches, then regenerate thresholds,
  candidate classifications, fits and folds as required. Do not mix old and
  corrected results under one run label.

These are outstanding run-configuration checks. Repository cleanup alone does
not validate the corrected observation analysis or regenerate the paper.

## What Git ignores

The root `.gitignore` is an exact allowlist: all paths are ignored unless
explicitly included. This applies to new Python and shell files too. Local
tests, alternate ephemerides, `simulate_dataset/`, `code/tests/`,
`psr_miner_search/`, logs, backup packages, observation files and generated
outputs stay on the cluster and are excluded from new commits by default.
The three candidate diagnostic PDFs linked above are explicit exceptions;
other files under `docs/` remain ignored unless individually allowed.

To include another file, add its exact path and any missing parent-directory
exceptions to `.gitignore`, then inspect `git status` and the staged diff.
Avoid `git add -f`, which bypasses this policy. Files removed from tracking
by this cleanup are still on disk. Earlier public commits retain their old
contents; no published history has been rewritten.

## Retention and timing in the corrected rerun

`flag_rfi.py` reports every selected candidate, including known pulsars, and
records advisory frequency flags. It has no move action. All selected
candidates remain eligible for fitting and folding. These frequency
heuristics alone do not establish that a signal is terrestrial.
`collate_fold_snr.py` joins by candidate ID and retains missing/failed folds
with an explicit status. `select_coherent_candidates.py` writes every
selected candidate, including missing or unreadable tracks. Path variability
is descriptive; it is not by itself a binary classification. There is no
cross-timescale raw-likelihood ranking or top-30 exclusion.

Freeze the source and grid before submitting `stage3_batch.sh`. Set
`SEARCH_CODE_DIR`, `SEARCH_OUTPUT_BASE`, `SEARCH_GRID_SIZE` and optionally
`SEARCH_CONFIGS_PER_WORKER` (default 16). For 16,464 configurations this gives
1,029 workers, submitted as `--array=0-1028%16`. Each worker runs its
configurations serially on one CPU. This stays below the site's array and
pending-job limits. `stage3_blind.sh` remains the one-configuration launcher.
Both write a separate `_timing.json` per search through `measure_search.py`.
The record includes complete command wall/CPU time, peak memory, CPU model,
affinity, thread settings, return status and source/parameter hashes.
Dedispersion and candidate follow-up are excluded from this search timer.
Sum CPU times over the complete grid and report scheduler elapsed time
separately. Retain and report any hardware differences.
