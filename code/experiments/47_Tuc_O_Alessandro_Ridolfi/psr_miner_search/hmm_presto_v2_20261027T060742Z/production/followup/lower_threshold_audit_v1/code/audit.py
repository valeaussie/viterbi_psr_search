"""Audit saved HMM endpoint scores. Never rerun a search, track fit, or fold.

All writes are confined to a new audit directory. Input hashes and original
selection reproduction make the analysis reviewable. Frequency proximity is
an inspection lead, never a confirmed pulsar association.
"""
from pathlib import Path
from collections import Counter
import argparse
import bisect
import csv
import hashlib
import io
import json
import math
import os
import resource
import subprocess
import time
from datetime import datetime, timezone

import numpy as np

E = Path('/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi')
P = E / 'psr_miner_search/hmm_presto_v2_20261027T060742Z/production'
W = P / 'followup/lower_threshold_audit_v1'
BLIND = E / '47Tuc_blind_search_v1/stage3_viterbi/blind_v2'
ORIGINAL = E / '47Tuc_blind_search_v1/stage3_viterbi/candidates'
SEGMENTS = [16, 32, 64, 128, 256, 512, 1024]
FRACTIONS = [1., .9, .8, .7, .6, .5, .4, .3, .2, .1, 0.]
TARGETS = ['C', 'H', 'L', 'M', 'E']
TOL = .035


def stamp():
    return datetime.now(timezone.utc).isoformat()


def sha_bytes(data):
    return hashlib.sha256(data).hexdigest()


def write_json(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n')


def write_csv(path, rows, fields=None):
    if fields is None:
        fields = list(dict.fromkeys(k for r in rows for k in r))
    with path.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def config():
    return json.loads((W / 'inputs.json').read_text())


def initialize():
    """Read small metadata files only, before any compute job is submitted."""
    grid = json.loads((P / 'grid_manifest.json').read_text())
    assert len(grid['tasks']) == 16464
    thresholds = {}
    old_rows = list(csv.DictReader((ORIGINAL / 'candidates_raw.csv').open()))
    for row in old_rows:
        nt = int(row['Nt'])
        threshold = float(row['threshold'])
        if str(nt) in thresholds:
            assert thresholds[str(nt)] == threshold
        thresholds[str(nt)] = threshold
    assert sorted(map(int, thresholds)) == SEGMENTS
    assert len(old_rows) == 4397
    cat_path = P / 'followup/pipeline_assessment_v1/inputs/merged_identity_reference.json'
    catalog = json.loads(cat_path.read_text())
    selected = [r for r in catalog if r['name'] in ['J0024-7204' + x for x in TARGETS]]
    assert len(selected) == 5
    controls = []
    targeted = P / 'followup/targeted_native_v1/manifest.json'
    if targeted.exists():
        for row in json.loads(targeted.read_text())['selected_runs']:
            if row['target'] in TARGETS:
                controls.append({k: row.get(k) for k in [
                    'target', 'dm', 'segments', 'curve', 'track', 'summary',
                    'peak_frequency_hz', 'peak_loglike', 'nq', 'segment_duration_s',
                    'catalogue_frequency_hz', 'threshold', 'selection_scope']})
    frozen = [P / 'grid_manifest.json', cat_path,
              ORIGINAL / 'candidates_raw.csv',
              ORIGINAL / 'candidates_timescale_alternatives.csv',
              ORIGINAL / 'threshold_calibration.txt']
    inputs = dict(created_utc=stamp(), thresholds=thresholds, targets=selected,
                  tasks=grid['tasks'], fractions=FRACTIONS, tolerance_hz=TOL,
                  source_root=str(BLIND), target_controls=controls,
                  source_hashes={str(p): sha_bytes(p.read_bytes()) for p in frozen},
                  code_sha256=sha_bytes(Path(__file__).read_bytes()),
                  rule='L_cut(Nt,r) = median(all endpoint scores at Nt) + r * (original_cut(Nt) - median)',
                  notes=[
                      'Fractions are operating points, not false-alarm probabilities or calibrated significances.',
                      'Cuts are applied identically at every searched frequency for each Nt.',
                      'Greedy 0.035-Hz grouping reproduces the original frequency-only grouping within Nt; it is not a source count.',
                      'All DM/subband appearances are correlated; multiplicity is not independent evidence.',
                      'Timescale alternatives are retained; raw likelihoods are not compared between Nt.',
                      'Window: conservative orbital Doppler envelope + 0.001 Hz per harmonic + two bins.',
                      'Harmonic windows are location probes, not harmonic summing or recovery confirmation.',
                      'The median reference tier is deliberately permissive and may be overwhelmed by noise.'
                  ])
    write_json(W / 'inputs.json', inputs)
    write_csv(W / 'targeted_controls.csv', controls)
    print(json.dumps(dict(tasks=len(grid['tasks']), targets=[r['name'] for r in selected],
                          thresholds=thresholds, controls=len(controls))), flush=True)


def load_curve(task):
    path = BLIND / task['curve']
    before = path.stat()
    start = time.perf_counter()
    raw = path.read_bytes()
    array = np.loadtxt(io.BytesIO(raw))
    after = path.stat()
    assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
    assert array.ndim == 2 and array.shape[1] == 2 and array.shape[0] >= 3
    assert np.isfinite(array).all() and np.all(np.diff(array[:, 0]) > 0)
    peaks = np.flatnonzero((array[1:-1, 1] >= array[:-2, 1]) &
                           (array[1:-1, 1] > array[2:, 1])) + 1
    assert np.allclose(np.diff(array[:, 0]), np.median(np.diff(array[:, 0])), rtol=1e-5, atol=1e-10)
    info = dict(index=task['index'], source_curve=task['curve'], Nt=task['Nt'],
                dm=task['dm'], subband_f0=task['subband_f0'], samples=len(array),
                frequency_lo_hz=float(array[0, 0]), frequency_hi_hz=float(array[-1, 0]),
                bin_hz=float(np.median(np.diff(array[:, 0]))), bytes=len(raw),
                sha256=sha_bytes(raw), mtime_ns=before.st_mtime_ns,
                load_and_peak_seconds=time.perf_counter()-start)
    return array, peaks, info


def location_probes(task, array, peaks, inputs):
    freqs, scores = array[:, 0], array[:, 1]
    df = float(np.median(np.diff(freqs)))
    cutoff = inputs['thresholds'][str(task['Nt'])]
    records = []
    for cat in inputs['targets']:
        fspin = float(cat['F0_hz'])
        velocity_fraction = 0.
        if cat.get('PB_days') and cat.get('A1_lts'):
            ecc = float(cat.get('ECC') or 0.)
            velocity_fraction = (2 * math.pi * float(cat['A1_lts']) /
                                 (float(cat['PB_days']) * 86400.) *
                                 (1 + ecc) / math.sqrt(1 - ecc**2))
        for harmonic in range(1, math.floor(1059 / fspin) + 1):
            centre = harmonic * fspin
            halfwidth = centre * velocity_fraction + .001 * harmonic + 2 * df
            if centre + halfwidth < freqs[0] or centre - halfwidth > freqs[-1]:
                continue
            candidates = peaks[(freqs[peaks] >= centre-halfwidth) &
                               (freqs[peaks] <= centre+halfwidth)]
            nearest_bin = int(np.argmin(np.abs(freqs-centre)))
            nearest_peak = int(peaks[np.argmin(np.abs(freqs[peaks]-centre))]) if len(peaks) else None
            highest = int(candidates[np.argmax(scores[candidates])]) if len(candidates) else None
            row = dict(target=cat['name'].replace('J0024-7204', ''), harmonic=harmonic,
                       Nt=task['Nt'], dm=task['dm'], catalogue_dm=cat['DM_pccm3'],
                       within_catalogue_dm_0p3=abs(float(task['dm'])-float(cat['DM_pccm3'])) <= .30000001,
                       subband_f0=task['subband_f0'], source_curve=task['curve'],
                       source_index=task['index'], reference_hz=centre,
                       window_halfwidth_hz=halfwidth, orbital_envelope_hz=centre*velocity_fraction,
                       bin_hz=df, original_cutoff=cutoff,
                       reference_inside_curve=bool(freqs[0] <= centre <= freqs[-1]),
                       reference_bin=nearest_bin, reference_bin_hz=float(freqs[nearest_bin]),
                       reference_bin_score=float(scores[nearest_bin]),
                       window_peak_count=len(candidates))
            for prefix, index in [('nearest_peak', nearest_peak), ('window_peak', highest)]:
                row[prefix+'_bin'] = index
                row[prefix+'_hz'] = float(freqs[index]) if index is not None else None
                row[prefix+'_score'] = float(scores[index]) if index is not None else None
                row[prefix+'_offset_hz'] = float(freqs[index]-centre) if index is not None else None
                row[prefix+'_inside_window'] = bool(abs(freqs[index]-centre) <= halfwidth) if index is not None else False
                row[prefix+'_above_original'] = bool(scores[index] > cutoff) if index is not None else False
            records.append(row)
    return records


def pilot():
    inputs = config()
    timings = []
    for nt in SEGMENTS:
        tasks = [r for r in inputs['tasks'] if r['Nt'] == nt]
        # Sample 24 paths spread across the DM/subband grid for each Nt.
        for ix in np.linspace(0, len(tasks)-1, 24, dtype=int):
            array, peaks, info = load_curve(tasks[int(ix)])
            location_probes(tasks[int(ix)], array, peaks, inputs)
            timings.append(info)
    groups = []
    for nt in SEGMENTS:
        rows = [r for r in timings if r['Nt'] == nt]
        times = [r['load_and_peak_seconds'] for r in rows]
        groups.append(dict(Nt=nt, pilot_curves=len(rows),
                           median_seconds=float(np.median(times)),
                           p95_seconds=float(np.percentile(times, 95)),
                           projected_read_seconds=float(np.percentile(times, 95)*2352)))
    write_json(W / 'pilot.json', dict(utc=stamp(), timings=timings, groups=groups,
                                    max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    print(json.dumps(groups), flush=True)


def worker(nt):
    started = time.perf_counter()
    inputs = config()
    assert sha_bytes(Path(__file__).read_bytes()) == inputs['code_sha256']
    folder = W / ('Nt' + str(nt))
    folder.mkdir(exist_ok=False)
    tasks = [r for r in inputs['tasks'] if r['Nt'] == nt]
    assert len(tasks) == 2352
    all_scores, peak_scores, peak_freqs, peak_bins, peak_curves = [], [], [], [], []
    inventory, probes = [], []
    for ix, task in enumerate(tasks):
        array, peaks, info = load_curve(task)
        all_scores.append(array[:, 1].copy())
        peak_scores.append(array[peaks, 1].copy())
        peak_freqs.append(array[peaks, 0].copy())
        peak_bins.append(peaks.astype(np.int32))
        peak_curves.append(np.full(len(peaks), ix, dtype=np.int32))
        inventory.append(info)
        probes.extend(location_probes(task, array, peaks, inputs))
    combined = np.concatenate(all_scores)
    median = float(np.median(combined))
    old_cut = inputs['thresholds'][str(nt)]
    assert old_cut > median
    score_quantiles = {str(p):float(np.percentile(combined, p)) for p in [50,90,95,99,99.9,99.99]}
    del combined, all_scores, array
    ps, pf, pb, pc = map(np.concatenate, [peak_scores, peak_freqs, peak_bins, peak_curves])
    keep = ps > median
    ps, pf, pb, pc = ps[keep], pf[keep], pb[keep], pc[keep]
    del peak_scores, peak_freqs, peak_bins, peak_curves, keep
    dms = np.array([r['dm'] for r in tasks])
    bands = np.array([r['subband_f0'] for r in tasks])
    order = np.lexsort((pb, bands[pc], dms[pc], pf, -ps))
    ps, pf, pb, pc = ps[order], pf[order], pb[order], pc[order]
    del order
    # A compact ledger preserves every audit peak, without huge Python dict lists.
    with (folder/'peaks_above_median.npz').open('xb') as out:
        np.savez_compressed(out, score=ps, frequency_hz=pf, source_bin=pb,
                            curve_index=pc, global_curve_index=np.array([t['index'] for t in tasks]))
    write_csv(folder/'curve_inventory.csv', inventory)
    for row in probes:
        row['median_at_Nt'] = median
        for prefix in ['nearest_peak', 'window_peak']:
            val = row[prefix+'_score']
            row[prefix+'_excess_fraction'] = (val-median)/(old_cut-median) if val is not None else None
    write_csv(folder/'target_locations.csv', probes)
    # Representatives are created in the same order as the historical greedy
    # algorithm. Choosing the earliest matching representative reproduces it.
    anchors, anchor_ids, families = [], [], []
    counts = []
    cursor = 0
    for fraction in FRACTIONS:
        cut = median + fraction*(old_cut-median) if fraction != 1 else old_cut
        while cursor < len(ps) and ps[cursor] > cut:
            freq, score, ci, bi = float(pf[cursor]), float(ps[cursor]), int(pc[cursor]), int(pb[cursor])
            left, right = bisect.bisect_left(anchors, freq-TOL), bisect.bisect_right(anchors, freq+TOL)
            if left < right:
                family = families[min(anchor_ids[left:right])]
            else:
                fid = len(families)
                task = tasks[ci]
                family = dict(family_id=f'Nt{nt}_family{fid:06d}', Nt=nt,
                              representative_hz=freq, representative_score=score,
                              representative_dm=task['dm'], representative_subband_f0=task['subband_f0'],
                              source_curve=task['curve'], source_bin=bi,
                              excess_fraction=(score-median)/(old_cut-median),
                              originally_selected=score > old_cut, raw_members=0,
                              dm_set=set(), subband_set=set())
                families.append(family)
                where = bisect.bisect_left(anchors, freq)
                anchors.insert(where, freq)
                anchor_ids.insert(where, fid)
            family['raw_members'] += 1
            family['dm_set'].add(float(tasks[ci]['dm']))
            family['subband_set'].add(int(tasks[ci]['subband_f0']))
            cursor += 1
        counts.append(dict(Nt=nt, excess_fraction_cut=fraction, score_cut=cut,
                           median=median, original_cut=old_cut,
                           raw_peaks=cursor, grouped_candidates=len(families),
                           groups_with_2_or_more_DMs=sum(len(f['dm_set'])>=2 for f in families),
                           groups_in_2_or_more_subbands=sum(len(f['subband_set'])>=2 for f in families)))
    family_rows = []
    for f in families:
        row = {k:v for k,v in f.items() if k not in ['dm_set', 'subband_set']}
        row.update(dm_count=len(f['dm_set']), dm_values=';'.join(f'{x:.2f}' for x in sorted(f['dm_set'])),
                   subband_count=len(f['subband_set']), subband_starts=';'.join(map(str,sorted(f['subband_set']))),
                   membership_scope='All peaks above the Nt median; correlated appearances, not independent detections')
        family_rows.append(row)
    write_csv(folder/'families_above_median.csv', family_rows)
    write_csv(folder/'counts.csv', counts)
    old_expected = [r for r in csv.DictReader((ORIGINAL/'candidates_raw.csv').open()) if int(r['Nt']) == nt]
    original_ids = {(r['source_curve'], int(r['source_bin'])) for r in old_expected}
    n_orig = counts[0]['raw_peaks']
    actual_ids = {(tasks[int(pc[j])]['curve'], int(pb[j])) for j in range(n_orig)}
    assert actual_ids == original_ids, (nt, len(actual_ids), len(original_ids))
    original_groups = [r for r in csv.DictReader((ORIGINAL/'candidates_timescale_alternatives.csv').open()) if int(r['Nt_best']) == nt]
    assert counts[0]['grouped_candidates'] == len(original_groups)
    # Check the actual representative locations as well as the counts.
    expected_reps={(r['source_curve'],int(r['source_bin'])) for r in original_groups}
    actual_reps={(r['source_curve'],r['source_bin']) for r in family_rows if r['originally_selected']}
    assert actual_reps == expected_reps
    result = dict(Nt=nt, curves=len(tasks), bins=sum(r['samples'] for r in inventory),
                  score_quantiles=score_quantiles, baseline_exactly_reproduced=True,
                  raw_original=n_orig, groups_original=len(original_groups),
                  max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  wall_s=time.perf_counter()-started, job_id=os.environ.get('SLURM_JOB_ID'),
                  code_sha256=inputs['code_sha256'], completed_utc=stamp())
    write_json(folder/'complete.json', result)
    print(json.dumps(result), flush=True)


def report():
    inputs = config()
    counts, locations, completions = [], [], []
    for nt in SEGMENTS:
        folder = W/('Nt'+str(nt))
        completions.append(json.loads((folder/'complete.json').read_text()))
        counts.extend(list(csv.DictReader((folder/'counts.csv').open())))
        locations.extend(list(csv.DictReader((folder/'target_locations.csv').open())))
    assert sum(r['curves'] for r in completions) == 16464
    assert sum(r['raw_original'] for r in completions) == 4397
    assert sum(r['groups_original'] for r in completions) == 415
    for path, digest in inputs['source_hashes'].items():
        assert sha_bytes(Path(path).read_bytes()) == digest, path
    write_csv(W/'counts_by_timescale.csv', counts)
    totals = []
    for fraction in FRACTIONS:
        rows = [r for r in counts if float(r['excess_fraction_cut']) == fraction]
        totals.append(dict(excess_fraction_cut=fraction,
                           raw_peaks=sum(int(r['raw_peaks']) for r in rows),
                           grouped_timescale_candidates=sum(int(r['grouped_candidates']) for r in rows),
                           additional_groups=sum(int(r['grouped_candidates']) for r in rows)-415,
                           groups_with_2_or_more_DMs=sum(int(r['groups_with_2_or_more_DMs']) for r in rows)))
    write_csv(W/'counts_total.csv', totals)
    write_csv(W/'target_locations_all.csv', locations)
    selected = []
    for target in TARGETS:
        harmonics = sorted({int(r['harmonic']) for r in locations if r['target'] == target})
        for harmonic in harmonics:
            for nt in SEGMENTS:
                for scope in ['any_searched_DM', 'within_0.3_of_catalogue_DM']:
                    rows = [r for r in locations if r['target']==target and int(r['harmonic'])==harmonic and int(r['Nt'])==nt
                            and (scope=='any_searched_DM' or r['within_catalogue_dm_0p3']=='True')]
                    candidates = [r for r in rows if r['window_peak_score']]
                    if candidates:
                        row = max(candidates, key=lambda r:float(r['window_peak_score']))
                        selected.append(dict(row, selection_scope=scope,
                                             interpretation='Highest nearby endpoint peak, not a confirmed target track'))
    write_csv(W/'target_best_by_timescale.csv', selected)
    fundamental = []
    for target in TARGETS:
        rows = [r for r in selected if r['target']==target and int(r['harmonic'])==1 and r['selection_scope']=='within_0.3_of_catalogue_DM']
        if rows:
            fundamental.append(max(rows, key=lambda r:float(r['window_peak_excess_fraction'])))
    write_csv(W/'fundamental_inspection_leads.csv', fundamental)
    # Compact human-readable report, with the denominator defined explicitly.
    lines = ['# Lower-threshold candidate-count audit', '',
             'This is a reselection audit of saved blind-search endpoint scores. No search, track recovery, fit, or fold was run.', '',
             'For each segment count, let m be the median of all saved endpoint scores and L_old the existing threshold. The audit cut is m + r (L_old - m). r=1 reproduces the old selection; r=0 is the deliberately permissive median tier. r is not a significance or false-alarm probability.', '',
             'Grouping exactly reproduces the historical greedy 0.035-Hz grouping within each segment count. DMs and overlapping subbands are correlated. Summing the groups across segment counts retains timescale alternatives and does not count unique physical sources.', '',
             '| r | Raw peaks | Grouped timescale candidates | Additional groups |',
             '|---:|---:|---:|---:|']
    lines += [f"| {r['excess_fraction_cut']:g} | {r['raw_peaks']:,} | {r['grouped_timescale_candidates']:,} | {r['additional_groups']:,} |" for r in totals]
    lines += ['', '## Fundamental-frequency inspection leads', '',
              'The following are the highest relative-to-threshold window peaks within 0.3 of each catalogue DM, chosen across timescales. They are proximity-based leads, not redetections. Every DM, timescale, subband and in-band integer harmonic is in the accompanying CSV files.', '',
              '| Target | Segments | DM | Peak Hz | Score | Old cut | r at peak |',
              '|---|---:|---:|---:|---:|---:|---:|']
    for r in fundamental:
        lines.append(f"| {r['target']} | {r['Nt']} | {float(r['dm']):.2f} | {float(r['window_peak_hz']):.9f} | {float(r['window_peak_score']):.3f} | {float(r['original_cutoff']):.3f} | {float(r['window_peak_excess_fraction']):.3f} |")
    lines += ['', '## Interpretation and limits', '',
              '- Pulsar windows use a conservative orbital Doppler envelope, plus 0.001 Hz per harmonic and two frequency bins. This is not a phase-connected orbital ephemeris.',
              '- The nearest reference-bin score, nearest local peak and highest peak in the declared window are distinct and are recorded separately.',
              '- A Viterbi endpoint near a pulsar frequency can belong to a path that followed noise, interference or another source. Confirmation needs the original selected path and pulse/dispersion evidence.',
              '- Lowering a cut cannot recover a feature that is not an interior local maximum. The original algorithm omits the outermost frequency bins.',
              '- Targeted runs use narrower bands, sometimes different historical frequency grids and segment counts. Their scores are not directly interchangeable with the blind curves.',
              '- Multiplicity and frequency-only groups are workload diagnostics, not independent evidence or a measured false-alarm rate.',
              '- All original 4,397 selected peaks and 415 representative groups were reproduced, including source-curve and source-bin identities.',
              '- All frozen selection/catalogue metadata hashes are unchanged. Every input curve has a recorded SHA-256 and stable before/after file metadata.', '',
              'Files: counts_total.csv, counts_by_timescale.csv, target_best_by_timescale.csv, target_locations_all.csv, targeted_controls.csv. Nt*/families_above_median.csv and Nt*/peaks_above_median.npz preserve the broad audit selection and raw peak identities.', '']
    with (W/'REPORT.md').open('x') as out:
        out.write('\n'.join(lines))
    write_json(W/'complete.json', dict(completed_utc=stamp(), curves=16464,
                                     baseline_reproduced=dict(raw_peaks=4397, groups=415),
                                     original_metadata_hashes_unchanged=True,
                                     workers=completions, totals=totals,
                                     fundamental_inspection_leads=fundamental))
    print(json.dumps(dict(totals=totals, fundamental_inspection_leads=fundamental)), flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('mode', choices=['initialize','pilot','worker','report'])
    parser.add_argument('--nt', type=int)
    a=parser.parse_args()
    if a.mode=='initialize': initialize()
    elif a.mode=='pilot': pilot()
    elif a.mode=='worker': worker(a.nt)
    else: report()
