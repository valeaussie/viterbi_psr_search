"""Native channelized PRESTO folds of saved six-coefficient HMM models.

Read-only scientific inputs. No refitting, new search or catalogue parameters.
TEMPO receives F0 through F5 at the saved barycentric reference epoch. Its
MeerKAT polycos convert that phase model to the filterbank time axis. A separate
barycentric-site polyco test checks the six-coefficient transfer before folding.
"""
from pathlib import Path
from decimal import Decimal, getcontext
from datetime import datetime, timezone
import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import time
import traceback

getcontext().prec = 40
EXPERIMENT = Path('/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi')
PRODUCTION = EXPERIMENT / 'psr_miner_search/hmm_presto_v2_20261027T060742Z/production'
RESTORED = PRODUCTION / 'followup/restored_all_415'
WORKFLOW = RESTORED / 'native_presto_hmm_v1'
CANDIDATES = EXPERIMENT / '47Tuc_blind_search_v1/stage3_viterbi/candidates'
FOLD_ROOT = EXPERIMENT / '47Tuc_blind_search_v1/stage4_fold/comparison_v1/hmm'
FILTERBANK = EXPERIMENT / '47Tuc_blind_search_v1/stage1_clean/47Tuc_22UL_1of2_L_RFIcleaned_01.fil'
DEDISP = EXPERIMENT / '47Tuc_blind_search_v1/stage2_dedisp'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2) + '\n')


def summary(data):
    result = {}
    for line in data.splitlines():
        fields = line.split()
        if not fields or fields[0].startswith('#'):
            continue
        if fields[0] in result:
            raise ValueError('Duplicate summary field: ' + fields[0])
        result[fields[0]] = fields[1:]
    return result


def inf_fields(path):
    return {line.split('=', 1)[0].strip(): line.split('=', 1)[1].strip()
            for line in Path(path).read_text().splitlines() if '=' in line}


def inf_get(data, prefix):
    found = [v for k, v in data.items() if k.startswith(prefix)]
    if len(found) != 1:
        raise ValueError((prefix, found))
    return found[0]


def old_fold_snapshot(cid):
    root = FOLD_ROOT / cid
    # Includes existing fold arrays, plots, models and completion record.
    paths = [root / name for name in ['complete.json', 'models.json', 'track_with_fits.dat', 'track_models.png']]
    for folder in ['polynomial', 'circular']:
        paths.extend(root / folder / name for name in ['fold_arrays.npz', 'fold.json', 'fold.png', 'profiles.txt'])
    if not all(path.is_file() for path in paths):
        raise ValueError('Expected original fold products: ' + cid)
    return {str(path): dict(sha256=sha(path), size=path.stat().st_size,
                            mtime_ns=path.stat().st_mtime_ns) for path in paths}


def prepare():
    if (WORKFLOW / 'manifest.json').exists():
        raise FileExistsError('Manifest already exists; inspect it instead of silently replacing it')
    from presto import sigproc
    header, header_size = sigproc.read_header(str(FILTERBANK))
    assert header['nchans'] == 128 and header['nbits'] == 8
    dt = float(header['tsamp'])
    nifs = int(header.get('nifs', 1))
    samples = (FILTERBANK.stat().st_size - header_size) * 8 // (128 * 8 * nifs)
    duration = samples * dt
    with (RESTORED / 'candidate_index.csv').open(newline='') as source:
        ids = [row['candidate_id'] for row in csv.DictReader(source)]
    assert len(ids) == len(set(ids)) == 415
    rows = []
    snapshots = {}
    for number, cid in enumerate(ids):
        directory = CANDIDATES / cid
        spath = directory / (cid + '_fit_summary.dat')
        values = summary(spath.read_text())
        candidate = json.loads((directory / 'candidate.json').read_text())
        assert values['frequency_frame'] == ['barycentric']
        assert candidate['selected_candidate_retained']
        assert values['poly_tref_s'] == values['kepler_tref_s']
        assert values['poly_tref_mjd'] == values['kepler_tref_mjd']
        dm_text = values['dm_pc_cm-3'][0]
        inf = DEDISP / ('DM%.2f_bary' % float(dm_text)) / ('47Tuc_DM%.2f_bary.inf' % float(dm_text))
        info = inf_fields(inf)
        assert inf_get(info, 'Barycentered?') == '1'
        assert Decimal(inf_get(info, 'Epoch of observation')) == Decimal(values['input_epoch_mjd'][0])
        assert abs(float(inf_get(info, 'Width of each time series bin')) - dt) < 1e-15
        assert abs(float(inf_get(info, 'Dispersion measure')) - float(dm_text)) < 1e-10
        snapshots[cid] = old_fold_snapshot(cid)
        for model, prefix, folder in [('polynomial', 'poly', 'polynomial'),
                                       ('kepler_taylor', 'kepler', 'kepler_taylor')]:
            coefficients = [values[prefix + '_f' + str(k)][0] for k in range(6)]
            assert all(math.isfinite(float(v)) for v in coefficients)
            epoch = values[prefix + '_tref_mjd'][0]
            exact_epoch = Decimal(values['input_epoch_mjd'][0]) + Decimal(values[prefix + '_tref_s'][0]) / Decimal(86400)
            assert abs(exact_epoch - Decimal(epoch)) < Decimal('1e-14')
            rows.append(dict(index=len(rows), candidate_number=number, candidate_id=cid,
                             candidate_json=str(directory / 'candidate.json'), summary=str(spath),
                             summary_sha256=sha(spath), inf=str(inf), inf_sha256=sha(inf),
                             model=model, coefficient_strings=coefficients,
                             reference_epoch_mjd=epoch, input_epoch_mjd=values['input_epoch_mjd'][0],
                             reference_time_s=values[prefix + '_tref_s'][0], dm=dm_text,
                             ra=inf_get(info, 'J2000 Right Ascension'),
                             dec=inf_get(info, 'J2000 Declination'),
                             source_valid_samples=int(candidate['sampling']['valid_samples']),
                             output_directory=str(FOLD_ROOT / cid / folder),
                             candidate_label='HMM%06d' % number,
                             classification=candidate['classification'], rfi_flag=candidate['rfi_flag']))
        if (number + 1) % 100 == 0:
            print('Prepared', number + 1, '/ 415 candidates', flush=True)
    WORKFLOW.mkdir(exist_ok=True)
    (WORKFLOW / 'logs').mkdir(exist_ok=True)
    write_json(WORKFLOW / 'original_fold_files.json', snapshots)
    manifest = dict(created_utc=datetime.now(timezone.utc).isoformat(), rows=rows,
                    candidates=415, models=2, folds=830, filterbank=str(FILTERBANK),
                    filterbank_metadata=dict(header=header, header_bytes=header_size,
                                             samples=samples, duration_s=duration,
                                             size_bytes=FILTERBANK.stat().st_size,
                                             mtime_ns=FILTERBANK.stat().st_mtime_ns),
                    solar_system_ephemeris='DE405', clock='UTC(NIST)',
                    purpose='Native channelized PRESTO follow-up using saved HMM coefficients; excluded from search timing',
                    settings=dict(ncpus=1, npfact=2, ndmfact=3, searchfdd=True,
                                  npart=64, nsub=32, profile_bins='PRESTO default from sample time and folding frequency'),
                    no_catalogue_parameters=True, no_refitting=True, no_candidate_exclusions=True)
    write_json(WORKFLOW / 'manifest.json', manifest)
    print(json.dumps({key: value for key, value in manifest.items() if key != 'rows'}, indent=2), flush=True)


def par_text(row, dm):
    # TZ mode decodes the first four name characters as integer HHMM. This is
    # only an internal neutral name. RAJ/DECJ still come from the search .inf.
    lines = ['PSR 0000+0000', 'RAJ ' + row['ra'], 'DECJ ' + row['dec'],
             'PEPOCH ' + row['reference_epoch_mjd'], 'POSEPOCH ' + row['reference_epoch_mjd'],
             'DM ' + str(dm), 'EPHEM DE405', 'CLK UTC(NIST)', 'NITS 1',
             'TZRMJD ' + row['reference_epoch_mjd'], 'TZRFRQ 0', 'TZRSITE @']
    lines.extend('F%d %s' % (k, value) for k, value in enumerate(row['coefficient_strings']))
    return '\n'.join(lines) + '\n'


def tempo_polycos(directory, row, site, frequency, dm, start, stop):
    directory.mkdir()
    (directory / 'model.par').write_text(par_text(row, dm))
    (directory / 'tz.in').write_text('%s 12 30 12 430\n\n\n%s 30 12 12 %.9f\n' %
                                    (site, '0000+0000', frequency))
    command = ['tempo', '-z', '-f', 'model.par']
    with (directory / 'tempo.log').open('w') as output:
        result = subprocess.run(command, cwd=directory, input='%d %d\n' % (math.floor(start) - 1, math.ceil(stop)),
                                text=True, stdout=output, stderr=subprocess.STDOUT)
    if result.returncode or not (directory / 'polyco.dat').is_file():
        print((directory / 'tempo.log').read_text(), flush=True)
        raise RuntimeError('TEMPO generation failed: ' + str(directory / 'tempo.log'))
    blocks = read_polycos(directory / 'polyco.dat', start, stop,
                          directory / 'polyco_for_observation.dat')
    assert blocks and all(block['ncoeff'] == 12 and block['span_minutes'] == 30 for block in blocks)
    return blocks


def read_polycos(path, coverage_start=None, coverage_stop=None, selected_path=None):
    lines = iter(Path(path).read_text().splitlines(keepends=True))
    blocks = []
    selected_text = []
    for line in lines:
        if not line.strip():
            continue
        first = line.split()
        second_line = next(lines)
        second = second_line.split()
        assert len(first) >= 5 and len(second) >= 6, (first, second)
        ncoeff = int(second[4])
        assert ncoeff % 3 == 0
        coefficient_lines = [next(lines) for _ in range(ncoeff // 3)]
        midpoint = Decimal(first[3])
        half_span = Decimal(second[3]) / Decimal(2880)
        # TEMPO generates complete calendar days. Far-away blocks can overflow
        # its printed phase field for steep high-order Taylor expansions. They
        # are not used to fold this observation. Preserve the raw file, but pass
        # only overlapping blocks to PRESTO, without changing any numeric field.
        if coverage_start is not None and (midpoint + half_span < Decimal(str(coverage_start))
                                          or midpoint - half_span > Decimal(str(coverage_stop))):
            continue
        coefficients = []
        for coefficient_line in coefficient_lines:
            coefficients.extend(float(x.replace('D', 'e')) for x in coefficient_line.split())
        assert len(coefficients) == ncoeff
        blocks.append(dict(name=first[0], midpoint=midpoint,
                           phase=Decimal(second[0].replace('D', 'E')), frequency=float(second[1].replace('D', 'e')),
                           span_minutes=int(second[3]), ncoeff=ncoeff,
                           coefficients=coefficients, header_line=first))
        selected_text.extend([line, second_line] + coefficient_lines)
    if selected_path is not None:
        Path(selected_path).write_text(''.join(selected_text))
    return blocks


def evaluate(blocks, epoch):
    possible = [block for block in blocks
                if abs(epoch - block['midpoint']) * Decimal(1440) < Decimal(block['span_minutes']) / 2]
    if not possible:
        raise ValueError('No polyco covers ' + str(epoch))
    block = min(possible, key=lambda b: abs(epoch - b['midpoint']))
    minutes = float((epoch - block['midpoint']) * Decimal(1440))
    correction = sum(c * minutes ** k for k, c in enumerate(block['coefficients']))
    phase = float(block['phase'] % 1) + minutes * 60 * block['frequency'] + correction
    frequency = block['frequency'] + sum(k * c * minutes ** (k - 1) / 60
                                        for k, c in enumerate(block['coefficients']) if k)
    return phase, frequency


def verify_bary_transfer(blocks, row, duration):
    f = [float(value) for value in row['coefficient_strings']]
    ref = Decimal(row['reference_epoch_mjd'])
    start = Decimal(row['input_epoch_mjd'])
    zero = evaluate(blocks, ref)[0] % 1
    phase_errors = []
    frequency_errors = []
    for index in range(1001):
        epoch = start + Decimal(str(duration)) * Decimal(index) / Decimal(1000 * 86400)
        seconds = float((epoch - ref) * Decimal(86400))
        phase, frequency = evaluate(blocks, epoch)
        expected_phase = sum(v * seconds ** (k + 1) / math.factorial(k + 1) for k, v in enumerate(f))
        expected_frequency = sum(v * seconds ** k / math.factorial(k) for k, v in enumerate(f))
        phase_errors.append(abs(math.remainder(phase - zero - expected_phase, 1)))
        frequency_errors.append(abs(frequency - expected_frequency))
    result = dict(samples=1001, frequency_coefficients=6,
                  max_phase_error_cycles=max(phase_errors),
                  max_frequency_error_hz=max(frequency_errors),
                  phase_tolerance_cycles=1e-4, frequency_tolerance_hz=1e-7,
                  arbitrary_constant_phase_removed=True)
    result['passed'] = result['max_phase_error_cycles'] < 1e-4 and result['max_frequency_error_hz'] < 1e-7
    if not result['passed']:
        raise ValueError('All-coefficient TEMPO transfer check failed: ' + json.dumps(result))
    return result


def render_ps(path):
    gs = shutil.which('gs')
    if not gs:
        raise RuntimeError('Ghostscript is required to preserve native PRESTO PDF and PNG diagnostics')
    pdf = path.with_suffix('.pdf')
    png = path.with_suffix('.png')
    # PGPLOT's PostScript is landscape content on a portrait coordinate system.
    # Let Ghostscript record its orientation in the PDF, then rasterize the PDF
    # so both formats display the native page upright without raster editing.
    subprocess.run([gs, '-q', '-dNOPAUSE', '-dBATCH', '-dSAFER', '-dEPSCrop',
                    '-dAutoRotatePages=/All', '-sDEVICE=pdfwrite', '-sOutputFile=' + str(pdf), str(path)], check=True)
    subprocess.run([gs, '-q', '-dNOPAUSE', '-dBATCH', '-dSAFER', '-dUseCropBox',
                    '-sDEVICE=png16m', '-r130', '-sOutputFile=' + str(png), str(pdf)], check=True)
    assert pdf.stat().st_size > 1000 and png.stat().st_size > 1000
    return pdf, png


def run(index, validate_only=False):
    if not validate_only and not shutil.which('gs'):
        raise RuntimeError('Ghostscript is required for the native PDF and PNG diagnostics')
    manifest = json.loads((WORKFLOW / 'manifest.json').read_text())
    row = manifest['rows'][index]
    assert row['index'] == index
    assert sha(row['summary']) == row['summary_sha256'], 'Summary changed after manifest preparation'
    assert sha(row['inf']) == row['inf_sha256'], 'Input metadata changed'
    source = Path(manifest['filterbank'])
    assert source.stat().st_size == manifest['filterbank_metadata']['size_bytes']
    assert source.stat().st_mtime_ns == manifest['filterbank_metadata']['mtime_ns']
    target = Path(row['output_directory'])
    prefix = 'presto_' + row['model'] + '_f0_f5'
    marker = target / (prefix + '_complete.json')
    if marker.exists():
        prior = json.loads(marker.read_text())
        assert prior['row'] == row
        assert all(sha(target / name) == value for name, value in prior['product_sha256'].items())
        print('Already complete:', marker, flush=True)
        return
    if not validate_only:
        if any(target.glob(prefix + '*')):
            raise FileExistsError('Existing native products without completion marker: ' + str(target))
    scratch = Path(tempfile.mkdtemp(prefix='hmm_native_%04d_' % index,
                                   dir=os.environ.get('SLURM_TMPDIR', '/tmp')))
    started = time.perf_counter()
    record = dict(row=row, source=str(source), scratch=str(scratch),
                  started_utc=datetime.now(timezone.utc).isoformat(),
                  job_id=os.environ.get('SLURM_JOB_ID'), code_sha256=sha(__file__),
                  prepfold_executable=shutil.which('prepfold'), tempo_executable=shutil.which('tempo'))
    try:
        duration = manifest['filterbank_metadata']['duration_s']
        start_bary = float(row['input_epoch_mjd'])
        stop_bary = start_bary + duration / 86400
        print('Checking all six coefficients:', row['candidate_id'], row['model'], flush=True)
        bary = tempo_polycos(scratch / 'bary_check', row, '@', 0., 0., start_bary, stop_bary)
        record['barycentric_coefficient_transfer'] = verify_bary_transfer(bary, row, duration)
        print(json.dumps(record['barycentric_coefficient_transfer']), flush=True)
        header = manifest['filterbank_metadata']['header']
        frequencies = [float(header['fch1']) + k * float(header['foff']) for k in range(int(header['nchans']))]
        record['polyco_observing_frequency_mhz'] = (min(frequencies) + max(frequencies)) / 2
        start_topo = float(header['tstart'])
        topo = tempo_polycos(scratch / 'topocentric', row, 'm', record['polyco_observing_frequency_mhz'],
                             row['dm'], start_topo, start_topo + duration / 86400)
        # Coverage and finite phase/frequency checks on the actual topocentric time axis.
        topo_frequencies = []
        for sample in range(1001):
            epoch = Decimal(str(start_topo)) + Decimal(str(duration)) * Decimal(sample) / Decimal(1000 * 86400)
            phase, frequency = evaluate(topo, epoch)
            assert math.isfinite(phase) and math.isfinite(frequency)
            topo_frequencies.append(frequency)
        record['topocentric_frequency_range_hz'] = [min(topo_frequencies), max(topo_frequencies)]
        if min(topo_frequencies) <= 0:
            raise ValueError('Saved Taylor expansion has a non-positive frequency within the observation; no truncation or refit applied')
        if validate_only:
            write_json(WORKFLOW / ('transfer_check_%04d.json' % index), record)
            print('Transfer and coverage checks passed:', index, flush=True)
            return
        (scratch / 'input.fil').symlink_to(source)
        # Relabel only the neutral TEMPO identifiers. Retain its phase/frequency
        # coefficients and every other character, and preserve the raw output.
        raw_polycos = (scratch / 'topocentric/polyco_for_observation.dat').read_text()
        relabelled = re.sub(r'(?m)^(\s*)0000\+0000(?=\s)',
                             lambda match: match.group(1) + row['candidate_label'], raw_polycos)
        assert relabelled != raw_polycos
        assert all(block['name'] == '0000+0000' for block in topo)
        (scratch / 'model.polycos').write_text(relabelled)
        command = ['prepfold', '-filterbank', '-topo', '-noxwin', '-ncpus', '1',
                   '-dm', row['dm'], '-npfact', '2', '-ndmfact', '3', '-nsub', '32',
                   '-npart', '64', '-psr', row['candidate_label'], '-polycos', 'model.polycos',
                   '-searchfdd', '-o', prefix, 'input.fil']
        record['command'] = command
        print('Starting native PRESTO fold:', index, flush=True)
        with (scratch / 'prepfold.log').open('w') as log:
            result = subprocess.run(command, cwd=scratch, stdout=log, stderr=subprocess.STDOUT)
        record['returncode'] = result.returncode
        pfd_paths = list(scratch.glob(prefix + '*.pfd'))
        ps_paths = list(scratch.glob(prefix + '*.pfd.ps'))
        if result.returncode or len(pfd_paths) != 1 or len(ps_paths) != 1:
            raise RuntimeError('Native PRESTO fold did not produce one complete PFD and diagnostic: ' + str(scratch))
        from presto.prepfold import pfd
        folded = pfd(str(pfd_paths[0]))
        assert folded.nsub == 32 and folded.npart == 64
        assert math.isfinite(float(folded.profs.sum()))
        record['native_profile_shape'] = list(folded.profs.shape)
        record['native_pfd_sha256'] = sha(pfd_paths[0])
        bestprof = Path(str(pfd_paths[0]) + '.bestprof')
        assert bestprof.is_file()
        record['bestprof_header'] = [line for line in bestprof.read_text().splitlines() if line.startswith('#')]
        render_ps(ps_paths[0])
        record['wall_s'] = time.perf_counter() - started
        record['completed_utc'] = datetime.now(timezone.utc).isoformat()
        target.mkdir(exist_ok=True)
        preserved = {
            prefix + '.par': scratch / 'topocentric/model.par',
            prefix + '.polycos': scratch / 'model.polycos',
            prefix + '_prepfold.log': scratch / 'prepfold.log',
            prefix + '_tempo.log': scratch / 'topocentric/tempo.log',
            prefix + '_bary_check.polycos': scratch / 'bary_check/polyco.dat',
            prefix + '_bary_check_observation.polycos': scratch / 'bary_check/polyco_for_observation.dat',
            prefix + '_bary_check.par': scratch / 'bary_check/model.par',
            prefix + '_bary_check_tempo.log': scratch / 'bary_check/tempo.log',
            prefix + '_tz.in': scratch / 'topocentric/tz.in',
            prefix + '_tempo_raw.polycos': scratch / 'topocentric/polyco.dat',
        }
        preserved.update({path.name: path for path in scratch.glob(prefix + '*.pfd*')})
        assert not any((target / name).exists() for name in preserved)
        for name, path in preserved.items():
            shutil.copyfile(path, target / name)
        record['product_sha256'] = {name: sha(target / name) for name in preserved}
        record['search_models_refitted'] = False
        record['frequency_coefficients_used'] = 6
        record['presto_local_refinement_enabled'] = True
        record['original_fold_files_preserved'] = True
        original = json.loads((WORKFLOW / 'original_fold_files.json').read_text())[row['candidate_id']]
        for path, before in original.items():
            assert sha(path) == before['sha256']
            assert Path(path).stat().st_mtime_ns == before['mtime_ns']
        write_json(marker, record)
        print('Completed:', marker, 'wall_s:', record['wall_s'], flush=True)
        shutil.rmtree(scratch)
    except Exception:
        error_directory = WORKFLOW / ('failure_details_%04d_v7' % index)
        error_directory.mkdir(exist_ok=True)
        for subdir in ['bary_check', 'topocentric']:
            source_directory = scratch / subdir
            if source_directory.exists():
                shutil.copytree(source_directory, error_directory / subdir, dirs_exist_ok=True)
        record['error_details'] = str(error_directory)
        record['wall_s'] = time.perf_counter() - started
        record['error'] = traceback.format_exc()
        write_json(WORKFLOW / ('failure_%04d_v7.json' % index), record)
        print(record['error'], flush=True)
        print('Failure scratch preserved:', scratch, flush=True)
        raise


def report():
    manifest = json.loads((WORKFLOW / 'manifest.json').read_text())
    results, missing = [], []
    for row in manifest['rows']:
        prefix = 'presto_' + row['model'] + '_f0_f5'
        marker = Path(row['output_directory']) / (prefix + '_complete.json')
        if marker.is_file():
            record = json.loads(marker.read_text())
            assert record['row'] == row
            assert all(sha(Path(row['output_directory']) / name) == value
                       for name, value in record['product_sha256'].items())
            results.append(record)
        else:
            missing.append(row['index'])
    summary_report = dict(expected_folds=830, completed_folds=len(results),
                          complete_candidates=sum(all(any(r['row']['index'] == 2 * c + m for r in results)
                                                       for m in [0, 1]) for c in range(415)),
                          missing_indices=missing,
                          total_worker_wall_s=sum(r['wall_s'] for r in results),
                          output_root=str(FOLD_ROOT),
                          all_original_fold_files_preserved=True)
    originals = json.loads((WORKFLOW / 'original_fold_files.json').read_text())
    for paths in originals.values():
        for path, before in paths.items():
            assert sha(path) == before['sha256'] and Path(path).stat().st_mtime_ns == before['mtime_ns']
    write_json(WORKFLOW / 'completion_report.json', summary_report)
    with (WORKFLOW / 'fold_index.csv').open('w', newline='') as output:
        writer = csv.DictWriter(output, fieldnames=['candidate_id', 'model', 'directory', 'pdf', 'png', 'wall_s', 'significance_header'])
        writer.writeheader()
        for record in results:
            row = record['row']
            products = record['product_sha256']
            writer.writerow(dict(candidate_id=row['candidate_id'], model=row['model'],
                                 directory=row['output_directory'],
                                 pdf=next(name for name in products if name.endswith('.pdf')),
                                 png=next(name for name in products if name.endswith('.png')),
                                 wall_s=record['wall_s'], significance_header='; '.join(
                                     line for line in record['bestprof_header'] if 'Prob(Noise)' in line)))
    print(json.dumps(summary_report, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare', action='store_true')
    parser.add_argument('--index', type=int)
    parser.add_argument('--validate-only', action='store_true')
    parser.add_argument('--report', action='store_true')
    args = parser.parse_args()
    if args.prepare:
        prepare()
    elif args.report:
        report()
    elif args.index is not None:
        run(args.index, args.validate_only)
    else:
        parser.error('Choose --prepare, --report or --index')
