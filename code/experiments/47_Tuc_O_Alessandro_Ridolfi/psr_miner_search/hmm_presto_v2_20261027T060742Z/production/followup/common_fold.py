"""Inspection folds in barycentric sample time, shared by both blind searches.

These are dedispersed time-series folds, not channelized confirmation folds.
The frequency refinement is follow-up work and is never part of search timing.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time

import numpy as np

NBIN = 64
NFINE = 1024
NSUB = 64
REFINE_CYCLES = np.linspace(-8., 8., 1025)


def phase(t, model):
    if model['kind'] == 'polynomial':
        u = t - model['reference_time_s']
        c = [0.] + [v / math.factorial(k+1) for k, v in enumerate(model['frequency_derivatives'])]
        return np.polynomial.polynomial.polyval(u, c)
    if model['kind'] == 'circular':
        f, K, Pb, T0 = model['parameters']
        w = 2 * np.pi / Pb
        return f*t + K/w * np.cos(w*(t-T0))
    raise ValueError(model['kind'])


def accumulate(x, dt, start, stop, model, df=0., nfine=NFINE):
    """Use sample j at j*dt; fold only the explicitly valid, searched samples."""
    n = stop-start
    sums = np.zeros((nfine, NBIN))
    counts = np.zeros((nfine, NBIN))
    variance_sum = 0.
    for j0 in range(start, stop, 1_000_000):
        j1 = min(j0+1_000_000, stop)
        j = np.arange(j0, j1, dtype=np.int64)
        t = (j-model.get('time_origin_sample', 0))*dt
        phi = phase(t, model) + df*(t-model.get('reference_time_s', 0.))
        bins = np.floor(np.remainder(phi, 1.)*NBIN).astype(np.int64)
        si = np.minimum((j-start)*nfine//n, nfine-1)
        weights = np.asarray(x[j0:j1], dtype=np.float64)
        weights -= weights.mean()
        variance_sum += np.dot(weights, weights)
        flat = si*NBIN+bins
        sums += np.bincount(flat, weights=weights, minlength=nfine*NBIN).reshape(nfine, NBIN)
        counts += np.bincount(flat, minlength=nfine*NBIN).reshape(nfine, NBIN)
    return sums, counts, variance_sum/n


def choose_frequency_offset(sums, counts, duration):
    # Approximate coherence on fine subintegrations proposes an offset; the
    # final profile is always folded again from the original samples.
    means = sums / np.maximum(counts, 1)
    means -= means.mean(axis=1, keepdims=True)
    coeff = np.fft.rfft(means, axis=1)[:, 1:5]
    u = (np.arange(len(sums))+.5)/len(sums)
    powers = []
    for cycles in REFINE_CYCLES:
        ph = np.exp(-2j*np.pi*cycles*u[:, None]*np.arange(1, 5)[None, :])
        powers.append(float(np.sum(np.abs(np.sum(coeff*ph, axis=0))**2)))
    powers = np.array(powers)
    idx = int(np.argmax(powers))
    return float(REFINE_CYCLES[idx]/duration), powers, idx


def profile_products(sums, counts):
    prof = sums.sum(axis=0)/np.maximum(counts.sum(axis=0), 1)
    subsums = sums.reshape(NSUB, -1, NBIN).sum(axis=1)
    subcounts = counts.reshape(NSUB, -1, NBIN).sum(axis=1)
    sub = subsums/np.maximum(subcounts, 1)
    centered = prof-np.median(prof)
    scale = 1.4826*np.median(np.abs(centered))
    score = float(np.max(centered)/scale) if scale>0 else 0.
    return prof, sub, score


def fold_model(x, dt, start, stop, model, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    wall0, cpu0 = time.perf_counter(), time.process_time()
    out.mkdir(parents=True, exist_ok=True)
    sums, counts, variance = accumulate(x, dt, start, stop, model)
    fixed, fixed_sub, fixed_score = profile_products(sums, counts)
    duration = (stop-start)*dt
    df, powers, best = choose_frequency_offset(sums, counts, duration)
    refined_sums, refined_counts, _ = accumulate(x, dt, start, stop, model, df=df)
    refined, refined_sub, refined_score = profile_products(refined_sums, refined_counts)
    np.savez_compressed(out/'fold_arrays.npz', fixed_profile=fixed, fixed_subintegrations=fixed_sub,
                        refined_profile=refined, refined_subintegrations=refined_sub,
                        fixed_sums=sums, fixed_counts=counts, refined_sums=refined_sums,
                        refined_counts=refined_counts, refinement_cycles=REFINE_CYCLES,
                        refinement_coherent_power=powers, sample_variance=variance)
    np.savetxt(out/'profiles.txt', np.c_[(np.arange(NBIN)+.5)/NBIN, fixed, refined],
               header='pulse_phase fixed_profile refined_profile (mean-subtracted input intensity)')
    metadata = dict(model=model, start_sample=start, end_sample_exclusive=stop, tsamp_s=dt,
                    duration_s=duration, bins=NBIN, subintegrations=NSUB,
                    refinement=dict(frequency_offset_hz=df, trials=len(REFINE_CYCLES),
                                    range_cycles=[-8., 8.], best_grid_index=best,
                                    at_boundary=bool(best in (0, len(REFINE_CYCLES)-1)),
                                    harmonics_used=4, final_profile_refolded_from_samples=True),
                    fixed_inspection_score=fixed_score, refined_inspection_score=refined_score,
                    score_warning='Peak/MAD inspection statistic, not a calibrated detection significance; frequency trials and two HMM models must be accounted for.',
                    wall_s=time.perf_counter()-wall0, cpu_s=time.process_time()-cpu0,
                    time_frame='barycentric', sample_time='sample j is at j*tsamp',
                    input_processing='one-million-sample chunks, each chunk mean subtracted, no padding folded')
    (out/'fold.json').write_text(json.dumps(metadata, indent=2)+'\n')
    fig, axes = plt.subplots(2, 2, figsize=(10, 6), constrained_layout=True)
    for col, (prof, sub, title) in enumerate(((fixed, fixed_sub, 'Fixed search-derived model'),
                                             (refined, refined_sub, 'Frequency refinement'))):
        axes[0, col].plot(np.arange(NBIN)/NBIN, prof, color='k', lw=1)
        axes[0, col].set(title=title, xlabel='Pulse phase', ylabel='Mean intensity')
        # Remove each row's baseline for display only.
        display = sub-np.median(sub, axis=1, keepdims=True)
        sd = np.std(display)
        axes[1, col].imshow(display, origin='lower', aspect='auto', extent=[0, 1, 0, duration/3600],
                             cmap='gray_r', vmin=-2*sd, vmax=4*sd)
        axes[1, col].set(xlabel='Pulse phase', ylabel='Time into searched interval (h)')
    fig.suptitle(out.parent.name+'\n'+out.name+'; barycentric dedispersed time series', fontsize=9)
    fig.savefig(out/'fold.png', dpi=130)
    plt.close(fig)
    return metadata


def load_fit_helpers(frozen):
    sys.path.insert(0, str(frozen))
    spec = importlib.util.spec_from_file_location('frozen_fits', frozen/'viterbi_search_and_fit.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def hmm_models(candidate, frozen, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    data = np.loadtxt(candidate/(candidate.name+'_track.dat'))
    t, f = data.T
    info = json.loads((candidate/'candidate.json').read_text())
    cache = Path(info['source_cache'])
    curve = np.loadtxt(cache/'source_loglike_curve.dat')
    df = float(np.median(np.diff(curve[:, 0])))
    helper = load_fit_helpers(frozen)
    poly = helper.fit_polynomial_at_midpoint(t, f, 5, df)
    kep = helper.fit_kepler(t, f, df, info['sampling']['segment_duration_s'])
    models = dict(polynomial=dict(kind='polynomial', reference_time_s=poly[0],
                                  frequency_derivatives=poly[1].tolist(),
                                  frequency_derivative_errors=poly[2].tolist(),
                                  chi2=poly[3], dof=poly[4], reduced_chi2=poly[5], order=poly[7]),
                  circular=dict(kind='circular', parameters=kep[0].tolist(),
                                parameter_names=['f_spin_hz', 'K_hz', 'Pb_s', 'T0_s'],
                                parameter_errors=kep[1].tolist(), chi2=kep[2], dof=kep[3], reduced_chi2=kep[4],
                                warning='Circular apparent-frequency model fitted to quantized track; orbital parameters are candidate estimates.'))
    (out/'models.json').write_text(json.dumps(models, indent=2)+'\n')
    np.savetxt(out/'track_with_fits.dat', np.c_[t, f, poly[6], kep[5]],
               header='time_s original_selected_track_hz polynomial_hz circular_hz')
    fig, ax = plt.subplots(figsize=(9, 3), constrained_layout=True)
    ax.step(t/3600, f, where='mid', color='0.5', label='Original blind track')
    ax.plot(t/3600, poly[6], label='Polynomial, order 5')
    ax.plot(t/3600, kep[5], label='Full circular model')
    ax.set(xlabel='Time from barycentric epoch (h)', ylabel='Frequency (Hz)', title=candidate.name)
    ax.legend(fontsize=8)
    fig.savefig(out/'track_models.png', dpi=140)
    plt.close(fig)
    return info, models


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--method', choices=['hmm', 'presto'], required=True)
    ap.add_argument('--index', type=int, required=True)
    a = ap.parse_args()
    manifest = json.loads(a.manifest.read_text())
    item = manifest[a.method][a.index]
    root = Path(manifest['output_root'])/a.method
    root.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    if a.method == 'hmm':
        candidate = Path(item['candidate_directory'])
        out = root/candidate.name
        if (out/'complete.json').exists():
            print('Already complete:', out)
            return
        out.mkdir(exist_ok=True)
        info, models = hmm_models(candidate, Path(manifest['frozen_code']), out)
        sample = info['sampling']
        source = Path(manifest['dedispersion_root'])/('DM%.2f_bary'%float(info['dm_best']))/('47Tuc_DM%.2f_bary.dat'%float(info['dm_best']))
        x = np.memmap(source, dtype='<f4', mode='r')
        start, stop, dt = 0, int(sample['valid_samples']), float(sample['tsamp_s'])
        result = {}
        for name, model in models.items():
            result[name] = fold_model(x, dt, start, stop, model, out/name)
        record = dict(candidate_id=candidate.name, method='hmm', source=str(source), candidate=info,
                      folds=result, total_wall_s=time.perf_counter()-started)
        (out/'complete.json').write_text(json.dumps(record, indent=2)+'\n')
        print('Completed HMM candidate:', candidate.name, 'wall', record['total_wall_s'])
    else:
        # Each worker owns an original PRESTO section and all its candidates.
        # This keeps the input pages reusable while preserving candidate IDs.
        results = []
        for row in item['candidates']:
            cid = row['candidate_id']
            out = root/cid
            if (out/'complete.json').exists():
                print('Already complete:', cid)
                continue
            out.mkdir(exist_ok=True)
            meta = row['section_metadata']
            source = Path(meta['source'])
            start, stop = int(meta['start_sample']), int(meta['end_sample_exclusive'])
            dt = float(meta['fft_duration_s'])/int(meta['fft_samples'])
            x = np.memmap(source, dtype='<f4', mode='r')
            model = dict(kind='polynomial', reference_time_s=0., time_origin_sample=start,
                         frequency_derivatives=[row['f_start_hz'], row['fdot_start_hz_s'], row['fddot_hz_s2']],
                         origin='Exact PRESTO binary r,z,w converted as installed prepfold.c',
                         source_binary=row['source_binary'], source_binary_sha256=row['source_binary_sha256'])
            result = fold_model(x, dt, start, stop, model, out/'acceleration')
            record = dict(candidate_id=cid, method='presto', source=str(source), candidate=row, folds=dict(acceleration=result))
            (out/'complete.json').write_text(json.dumps(record, indent=2)+'\n')
            results.append(cid)
            print('Completed PRESTO candidate:', cid, flush=True)
        print('Completed PRESTO section:', a.index, 'new candidates', len(results), 'wall', time.perf_counter()-started)


if __name__ == '__main__':
    main()
