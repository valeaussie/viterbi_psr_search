"""Append the saved circular model's local Taylor block to 415 current summaries.

No search, track recovery, refitting, folding, or plot generation is performed.
Only the selected current summaries are changed. Originals are backed up first.
The frozen helper functions are extracted in isolation to avoid importing the
search pipeline or relying on its Python environment.
"""
from pathlib import Path
from decimal import Decimal, getcontext
from types import SimpleNamespace
import argparse
import ast
import csv
import hashlib
import json
import math
import os
import stat
import tempfile

getcontext().prec = 40
ROOT = Path(__file__).resolve().parent.parent
EXPERIMENT = Path('/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi')
CANDIDATES = EXPERIMENT / '47Tuc_blind_search_v1/stage3_viterbi/candidates'
MODELS = EXPERIMENT / '47Tuc_blind_search_v1/stage4_fold/comparison_v1/hmm'
FROZEN = ROOT.parent.parent / 'frozen/code/viterbi_search_and_fit.py'
AUDIT = ROOT / 'kepler_taylor_augmentation'
MARKER = '# Kepler-derived Taylor expansion at the polynomial reference epoch'
DERIVED_KEYS = ['kepler_tref_s', 'kepler_tref_mjd'] + ['kepler_f%d' % k for k in range(6)]

def sha(path):
    result = hashlib.sha256()
    with path.open('rb') as source:
        for chunk in iter(lambda: source.read(1024*1024), b''):
            result.update(chunk)
    return result.hexdigest()

def parse_summary(text):
    values = {}
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        fields = line.split()
        if fields[0] in values:
            raise ValueError('Duplicate summary key: ' + fields[0])
        values[fields[0]] = fields[1:]
    return values

def isolated_frozen_helpers():
    tree = ast.parse(FROZEN.read_text())
    wanted = {'kepler_model', 'kepler_derivs'}
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in wanted]
    assert {n.name for n in nodes} == wanted
    isolated = ast.fix_missing_locations(ast.Module(body=nodes, type_ignores=[]))
    # These two frozen functions use only np.pi, np.sin and np.array. Their
    # scalar operations are supplied by the standard library; no imports or
    # other pipeline functions are executed.
    namespace = {'np': SimpleNamespace(pi=math.pi, sin=math.sin, array=list)}
    exec(compile(isolated, str(FROZEN), 'exec'), namespace)
    return namespace['kepler_model'], namespace['kepler_derivs']

def entry_state(path):
    if path.is_symlink():
        return {'symlink': str(path.readlink())}
    info = path.stat()
    return dict(bytes=info.st_size, sha256=sha(path), mtime_ns=info.st_mtime_ns,
                mode=stat.S_IMODE(info.st_mode))

def other_files_snapshot(ids, summaries):
    result = {}
    for cid in ids:
        for path in sorted((CANDIDATES / cid).rglob('*')):
            if path in summaries:
                continue
            if path.is_symlink() or path.is_file():
                result[str(path)] = entry_state(path)
    return result

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply', action='store_true', help='Otherwise calculate/check only; write no files.')
    args = parser.parse_args()
    model, derivatives = isolated_frozen_helpers()
    with (ROOT / 'candidate_index.csv').open(newline='') as source:
        ids = [r['candidate_id'] for r in csv.DictReader(source)]
    assert len(ids) == len(set(ids)) == 415
    with (ROOT / 'candidate_placement_index.csv').open(newline='') as source:
        summary_paths = {r['candidate_id']: Path(r['target']) for r in csv.DictReader(source)
                         if r['member'].endswith('_fit_summary.dat')}
    assert set(summary_paths) == set(ids)
    plans = []
    max_curve_difference = 0.0
    max_derivative_scaled_difference = 0.0
    for cid in ids:
        target = summary_paths[cid]
        assert target.parent == CANDIDATES / cid and target.is_file() and not target.is_symlink()
        original = target.read_bytes()
        assert original.endswith(b'\n'), 'Expected existing summary to end in a newline'
        values = parse_summary(original.decode())
        assert values['candidate_id'] == [cid] and values['frequency_frame'] == ['barycentric']
        assert values['selected_candidate_retained'] == ['True']
        saved_path = MODELS / cid / 'models.json'
        saved = json.loads(saved_path.read_text())
        parameters = saved['circular']['parameters']
        assert len(parameters) == 4 and all(math.isfinite(x) for x in parameters)
        fspin, amplitude, pb, t0 = parameters
        assert amplitude >= 0 and pb > 0
        for key, value in zip(('kepler_f_spin_hz', 'kepler_K_hz', 'kepler_Pb_s', 'kepler_T0_s'), parameters):
            assert math.isclose(float(values[key][0]), value, rel_tol=3e-15, abs_tol=1e-25), (cid, key)
        tref = saved['polynomial']['reference_time_s']
        assert abs(float(values['poly_tref_s'][0]) - tref) < 1e-9
        epoch = Decimal(values['input_epoch_mjd'][0])
        tref_mjd = epoch + Decimal(str(tref)) / Decimal(86400)
        assert abs(tref_mjd - Decimal(values['poly_tref_mjd'][0])) < Decimal('1e-14')
        coefficients = [float(x) for x in derivatives(tref, *parameters, 5)]
        assert len(coefficients) == 6 and all(math.isfinite(x) for x in coefficients)

        # Independent recurrence check: d_(k+2) = -omega^2 d_k for a sinusoid.
        omega = 2 * math.pi / pb
        phase = omega * (tref - t0)
        modulation = [-amplitude * math.sin(phase), -amplitude * omega * math.cos(phase)]
        for k in range(2, 6):
            modulation.append(-omega * omega * modulation[k-2])
        independent = [fspin + modulation[0]] + modulation[1:]
        for k in range(6):
            scale = max(abs(fspin), 1.) if k == 0 else max(amplitude * omega**k, 1e-300)
            error = abs(coefficients[k] - independent[k]) / scale
            max_derivative_scaled_difference = max(max_derivative_scaled_difference, error)
            assert error < 1e-10, (cid, k, error)
        assert abs(coefficients[0] - model(tref, *parameters)) < 1e-12

        # Compare the saved full fitted curve, not a truncated Taylor curve
        # extrapolated over the whole observation.
        reference_point_found = False
        curve_path = MODELS / cid / 'track_with_fits.dat'
        for line in curve_path.read_text().splitlines():
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            time_s, original_track, polynomial, circular = map(float, line.split())
            difference = abs(model(time_s, *parameters) - circular)
            assert difference < 1e-9, (cid, time_s, difference)
            max_curve_difference = max(max_curve_difference, difference)
            if abs(time_s - tref) < 1e-9:
                reference_point_found = True
                assert abs(coefficients[0] - circular) < 1e-9
        assert reference_point_found

        block_lines = [
            '', MARKER,
            '# Evaluated analytically from the saved circular fit; no refitting.',
            '# f_k is the k-th derivative of apparent frequency at kepler_tref.',
            '# f(t) = sum_k f_k * (t-kepler_tref)^k / k!; units: Hz s^(-k).',
            '# These are point estimates; no derivative uncertainty propagation was performed.',
            '# Figure 3(c) uses f0, f1 and f2; f3 through f5 are also retained.',
            'kepler_tref_s ' + values['poly_tref_s'][0],
            'kepler_tref_mjd ' + values['poly_tref_mjd'][0],
        ]
        block_lines += ['kepler_f%d %.17e' % (k, value) for k, value in enumerate(coefficients)]
        block = ('\n'.join(block_lines) + '\n').encode()
        present = [key in values for key in DERIVED_KEYS]
        assert not any(present) or all(present), 'Partial existing Taylor block: ' + cid
        already = all(present)
        if already:
            assert MARKER in original.decode()
            for k, value in enumerate(coefficients):
                assert float(values['kepler_f%d' % k][0]) == value
            assert values['kepler_tref_s'] == values['poly_tref_s']
            assert values['kepler_tref_mjd'] == values['poly_tref_mjd']
            start = original.index(('\n' + MARKER).encode())
            base = original[:start]
            assert original[start:] == block
        else:
            base = original
        plans.append(dict(candidate_id=cid, path=target, original=original, base=base,
                          block=block, already=already, coefficients=coefficients,
                          model_sha256=sha(saved_path), curve_sha256=sha(curve_path)))
    print(json.dumps(dict(candidates=len(plans), summaries_to_update=sum(not p['already'] for p in plans),
                          existing_complete_blocks=sum(p['already'] for p in plans),
                          max_saved_curve_difference_hz=max_curve_difference,
                          max_derivative_scaled_difference=max_derivative_scaled_difference,
                          frozen_helper_sha256=sha(FROZEN), apply=args.apply)), flush=True)
    if not args.apply:
        return

    AUDIT.mkdir(exist_ok=True)
    backups = AUDIT / 'original_summaries'
    backups.mkdir(exist_ok=True)
    before = other_files_snapshot(ids, set(summary_paths.values()))
    before_path = AUDIT / 'other_candidate_files_before.json'
    if before_path.exists():
        assert json.loads(before_path.read_text()) == before
    else:
        with before_path.open('x') as output:
            json.dump(before, output, indent=2)
            output.write('\n')
    for item in plans:
        backup = backups / item['path'].name
        if backup.exists():
            assert backup.read_bytes() == item['base']
        else:
            with backup.open('xb') as output:
                output.write(item['base'])
    for i, item in enumerate(plans):
        path = item['path']
        assert path.read_bytes() == item['original'], 'Concurrent summary change: ' + str(path)
        if not item['already']:
            mode = stat.S_IMODE(path.stat().st_mode)
            temporary = None
            try:
                with tempfile.NamedTemporaryFile(mode='wb', prefix='.kepler_taylor_', dir=path.parent,
                                                 delete=False) as output:
                    temporary = Path(output.name)
                    output.write(item['original'] + item['block'])
                    output.flush()
                    os.fsync(output.fileno())
                os.chmod(temporary, mode)
                assert path.read_bytes() == item['original']
                os.replace(temporary, path)
                temporary = None
            finally:
                if temporary is not None:
                    temporary.unlink()
        assert path.read_bytes() == item['base'] + item['block']
        updated = parse_summary(path.read_text())
        previous = parse_summary(item['base'].decode())
        assert all(updated[k] == v for k, v in previous.items())
        assert all(k in updated for k in DERIVED_KEYS)
        if (i+1) % 100 == 0:
            print('Verified complete Taylor blocks:', i+1, '/', len(plans), flush=True)
    after = other_files_snapshot(ids, set(summary_paths.values()))
    assert before == after, 'A non-summary candidate file changed'
    for item in plans:
        row = dict(candidate_id=item['candidate_id'], summary=str(item['path']),
                   original_sha256=hashlib.sha256(item['base']).hexdigest(),
                   updated_sha256=sha(item['path']), model_sha256=item['model_sha256'],
                   saved_curve_sha256=item['curve_sha256'],
                   kepler_coefficients=item['coefficients'])
        with (AUDIT / (item['candidate_id'] + '.json')).open('w') as output:
            json.dump(row, output, indent=2)
            output.write('\n')
    report = dict(candidates=415, summaries_complete=415,
                  changed_this_run=sum(not item['already'] for item in plans),
                  coefficients_per_candidate=6, order=5, same_reference_epoch_as_polynomial=True,
                  existing_summary_contents_preserved=True, original_summaries_backed_up=415,
                  other_candidate_files_verified_unchanged=len(before),
                  saved_curve_validation_max_difference_hz=max_curve_difference,
                  derivative_validation_max_scaled_difference=max_derivative_scaled_difference,
                  frozen_helper_sha256=sha(FROZEN), code_path=str(Path(__file__).resolve()),
                  code_sha256=sha(Path(__file__)), search_rerun=False, track_recovery=False,
                  refitted=False, folded=False, plots_regenerated=False,
                  updated_directory=str(CANDIDATES), backups=str(backups),
                  archived_ZIP_and_combined_PDFs_unchanged=True)
    (AUDIT / 'audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)
    print('Example appended block:\n' + plans[0]['block'].decode(), flush=True)

if __name__ == '__main__':
    main()
