# 47 Tuc: HMM search reproducibility

This repository contains the code and small fixed inputs for the targeted
redetection of 47 Tuc O and the blind HMM/Viterbi search of the MeerKAT
observation. It does not contain the observation itself, generated results,
exploratory runs, software regression tests or PRESTO benchmarking work.
Those files remain on the cluster.

The search code includes the corrections made in October 2026: exact segment
counts after removing PRESTO padding, consistent FFT frequency labels, and
candidate alternatives retained separately at each timescale. The corrected
observation search has not yet been rerun. Existing paper figures and candidate
counts must therefore be distinguished from results of the corrected code.

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

Classify the new table with `classify_known.py` before follow-up. Its
row-indexed classification file must correspond to that exact table.
`stage3b_fit.sh` uses `classified_timescale_alternatives.csv` in its new
follow-up output folder. Do not reuse the old 194-row classification or array
bounds. Set submission array sizes from the current grid/table and respect
the cluster's array-size limit.

## Before a new production run

No job is submitted by preparing or committing this repository. Review the
site-specific paths and the scientific run configuration before submission:

- Select new output directories rather than overwriting historical results.
- Check that targeted inputs are topocentric and blind-search inputs are the
  intended cleaned, barycentred DM series.
- Set the same corrected candidate directory and table in the later folding,
  collation and inspection steps. Some launcher defaults still refer to the
  historical `candidates/` directory.
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

To include another file, add its exact path and any missing parent-directory
exceptions to `.gitignore`, then inspect `git status` and the staged diff.
Avoid `git add -f`, which bypasses this policy. Files removed from tracking
by this cleanup are still on disk. Earlier public commits retain their old
contents; no published history has been rewritten.
