"""J-matched targeted O coherence-time scan, with preserved inputs and measured requests."""
from pathlib import Path
from datetime import datetime, timezone
from decimal import Decimal
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
import re
import subprocess
import sys

W = Path(__file__).resolve().parents[1]
CODE = W / 'code'
E = W.parents[1]
SOURCE = E.parents[1]
NT = [16, 20, 24, 26, 28, 30, 31, 32, 33, 34, 36, 38, 40, 44, 48, 56, 64]
DMS = ['24.30', '24.40', '24.50']


def now():
    return datetime.now(timezone.utc).isoformat()


def save(path, data):
    with Path(path).open('x') as out:
        json.dump(data, out, indent=2, allow_nan=False)
        out.write('\n')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def prepare():
    from presto import sigproc
    h = module(CODE / 'native_helper.py', 'J_prepare_helper')
    original = E / 'data/47Tuc_22UL_1of2_L.fil'
    header, size = sigproc.read_header(str(original))
    samples = (original.stat().st_size - size) * 8 // (header['nbits'] * header['nchans'] * header.get('nifs', 1))
    metadata = dict(header=header, header_bytes=size, samples=samples,
                    duration_s=samples * header['tsamp'], size_bytes=original.stat().st_size,
                    mtime_ns=original.stat().st_mtime_ns)
    rows = []
    for dm in DMS:
        dat = E / '47Tuc_blind_search_v1/stage2_dedisp' / ('DM%s_bary' % dm) / ('47Tuc_DM%s_bary.dat' % dm)
        inf = dat.with_suffix('.inf')
        info = h.inf_fields(inf)
        assert dat.is_file() and h.inf_get(info, 'Barycentered?') == '1'
        assert abs(float(h.inf_get(info, 'Dispersion measure')) - float(dm)) < 1e-9
        for nt in NT:
            cid = 'redetect_O_DM%s_Nt%d' % (dm, nt)
            directory = W / 'viterbi_out' / ('DM%s_Nt%d' % (dm, nt))
            directory.mkdir(exist_ok=True)
            assert not any(directory.iterdir()), 'Refuse preparing over existing search products'
            rows.append(dict(index=len(rows), candidate_id=cid, dm=dm, segments=nt,
                dat=str(dat), dat_size_bytes=dat.stat().st_size, dat_mtime_ns=dat.stat().st_mtime_ns,
                inf=str(inf), inf_sha256=sha(inf), tsamp=h.inf_get(info, 'Width of each time series bin'),
                input_epoch_mjd=h.inf_get(info, 'Epoch of observation'),
                ra=h.inf_get(info, 'J2000 Right Ascension'), dec=h.inf_get(info, 'J2000 Declination'),
                prefix=str(directory / cid), directory=str(directory)))
    limits = {}
    commands = {
        'association': ['sacctmgr','-nP','show','assoc','where','account=oz022','format=Account,User,QOS,DefaultQOS,GrpTRES,MaxJobs,MaxSubmitJobs'],
        'qos': ['sacctmgr','-nP','show','qos','where','name=normal','format=Name,MaxWall,MaxTRESPU,MaxTRESPA,GrpTRES,MaxJobsPU,MaxSubmitPU'],
        'slurm_config': ['scontrol','show','config'],
        'partitions': ['scontrol','show','partition'],
    }
    for name, command in commands.items():
        limits[name] = dict(command=command, output=subprocess.check_output(command, text=True))
    save(W / 'limits_at_preparation.json', dict(checked_utc=now(), limits=limits,
        documentation=['https://supercomputing.swin.edu.au/docs/1a-nt_oz/what-s-new.html',
        'https://supercomputing.swin.edu.au/docs/2-ozstar/oz-partition.html',
        'https://supercomputing.swin.edu.au/docs/2-ozstar/oz-slurm-create.html']))
    save(W / 'search_manifest.json', dict(created_utc=now(), rows=rows, segments=NT, dms=DMS,
        lower_frequency_hz=float(h.summary((E/"47Tuc_blind_search_v1/stage3_viterbi/candidates/cand_000052_DM24.40_Nt32_f378.304435Hz_sb374.000_bary/cand_000052_DM24.40_Nt32_f378.304435Hz_sb374.000_bary_fit_summary.dat").read_text())["selected_ending_frequency_hz"][0])-2., bandwidth_hz=4., num_harm=1, poly_order=5,
        search_input_frame='barycentric', search_input_preparation='Existing RFI-cleaned PRESTO stage2 time series; unchanged from the blind J comparison.',
        O_workflow_difference='This matched scan uses the same barycentric blind-search preparation as J. The older topocentric O scan and folds remain separate and unchanged.',
        filterbank=str(original), filterbank_metadata=metadata, target='O', scope='targeted redetection, not a new blind-search result',
        maximum_concurrent_tasks=None, code_sha256={p.name:sha(p) for p in CODE.glob('*.py')}))
    print('Prepared', len(rows), 'targeted search runs', flush=True)


def search(index):
    m = json.loads((W / 'search_manifest.json').read_text()); r = m['rows'][index]
    assert r['index'] == index and sha(r['inf']) == r['inf_sha256']
    assert Path(r['dat']).stat().st_size == r['dat_size_bytes'] and Path(r['dat']).stat().st_mtime_ns == r['dat_mtime_ns']
    assert all(sha(CODE / name) == value for name, value in m['code_sha256'].items())
    prefix = Path(r['prefix']); timing = Path(str(prefix) + '_timing.json')
    if timing.exists():
        assert json.loads(timing.read_text())['returncode'] == 0
        print('Search already complete', index); return
    if list(prefix.parent.glob(prefix.name + '*')):
        raise FileExistsError('Partial outputs require inspection: ' + str(prefix))
    command = [sys.executable, str(CODE / 'viterbi_search_and_fit.py'), '--infile', r['dat'],
        '--tsamp', r['tsamp'], '--segments', str(r['segments']), '--f0', str(m['lower_frequency_hz']),
        '--bw', str(m['bandwidth_hz']), '--poly-order', '5', '--plot-path', '--dm', r['dm'],
        '--mjd-start', r['input_epoch_mjd'], '--known-label', 'J0024-7204O', '--out_prefix', r['prefix']]
    with Path(str(prefix) + '.log').open('x') as out:
        subprocess.run([sys.executable, str(CODE / 'measure_search.py'), '--out', str(timing), '--', *command],
                       stdout=out, stderr=subprocess.STDOUT, check=True)
    sampling = json.loads(Path(str(prefix) + '_sampling.json').read_text())
    assert sampling['mode'] == 'exact_segments' and sampling['actual_segments'] == r['segments']
    assert sampling['frequency_axis_version'] == 2
    import numpy as np
    track = np.loadtxt(str(prefix) + '_track.dat'); curve = np.loadtxt(str(prefix) + '_loglikes.dat')
    assert len(track) == r['segments'] and np.isfinite(track).all() and np.isfinite(curve).all()
    peak = int(np.argmax(curve[:, 1])); assert abs(curve[peak, 0] - track[-1, 1]) < 1e-9
    save(prefix.parent / 'search_complete.json', dict(row=r, sampling=sampling,
        endpoint_frequency_hz=float(track[-1, 1]), peak_loglike=float(curve[peak, 1]),
        product_sha256={p.name:sha(p) for p in prefix.parent.glob(prefix.name + '*') if p.is_file()},
        completed_utc=now(), timing=json.loads(timing.read_text())))
    print('Completed J search', index, r['dm'], r['segments'], flush=True)


def fold(index):
    import numpy as np
    from presto import psr_utils
    # PRESTO's FFT rotation must retain odd profile lengths, here often 27 bins.
    def odd_safe_fft_rotate(arr, bins):
        arr=np.asarray(arr); f=np.arange(arr.size//2+1,dtype=float)
        return np.fft.irfft(np.exp(2j*np.pi*f*bins/arr.size)*np.fft.rfft(arr),n=arr.size)
    psr_utils.fft_rotate = odd_safe_fft_rotate
    m = json.loads((W / 'search_manifest.json').read_text()); r = m['rows'][index // 2]
    native = Path(r['directory']) / 'native_followup'
    marker = native / 'manifest.json'
    if not marker.exists():
        raise FileNotFoundError('Prepare fold manifests after search completion')
    h = module(CODE / 'native_helper.py', 'J_native_helper'); h.WORKFLOW = native
    h.run(index % 2)


def fold_prepare():
    m = json.loads((W / 'search_manifest.json').read_text())
    h = module(CODE / 'native_helper.py', 'J_fold_prepare_helper')
    for r in m['rows']:
        directory = Path(r['directory']); complete = directory / 'search_complete.json'
        if not complete.exists(): continue
        native = directory / 'native_followup'
        if (native / 'manifest.json').exists(): continue
        native.mkdir(exist_ok=True)
        fit = Path(r['prefix'] + '_fit_summary.dat'); values = h.summary(fit.read_text())
        sampling = json.loads(Path(r['prefix'] + '_sampling.json').read_text())
        rows = []
        for model, prefix in [('polynomial', 'poly'), ('kepler_taylor', 'kepler')]:
            coefficients = [values[prefix + '_f' + str(k)][0] for k in range(6)]
            assert all(math.isfinite(float(c)) for c in coefficients)
            target = native / model; target.mkdir()
            tref = values[prefix + '_tref_s'][0]
            rows.append(dict(index=len(rows), candidate_number=r['index'], candidate_id=r['candidate_id'],
                candidate_json='', summary=str(fit), summary_sha256=sha(fit), inf=r['inf'], inf_sha256=r['inf_sha256'],
                model=model, coefficient_strings=coefficients, reference_epoch_mjd=str(Decimal(r['input_epoch_mjd']) + Decimal(tref) / Decimal(86400)),
                input_epoch_mjd=r['input_epoch_mjd'], reference_time_s=tref, dm=r['dm'], ra=r['ra'], dec=r['dec'],
                source_valid_samples=sampling['valid_samples'], output_directory=str(target),
                candidate_label='OTG%06d' % r['index'], classification='targeted_O', rfi_flag='', segments=r['segments']))
        save(native / 'original_fold_files.json', {r['candidate_id']:{str(p):dict(sha256=sha(p),mtime_ns=p.stat().st_mtime_ns)
            for p in [fit, Path(r['prefix'] + '_track.dat'), Path(r['inf'])]}})
        save(native / 'manifest.json', dict(rows=rows, filterbank=m['filterbank'], filterbank_metadata=m['filterbank_metadata']))
    print('Prepared native fold manifests', flush=True)


def submit(mode, indices, minutes, memory_mb, dependency=None, label=None):
    path = W / 'batches' / ('%s_%s.json' % (mode, label or len(list((W/'batches').iterdir()))))
    save(path, indices)
    command = ['sbatch','--parsable','--account=oz022','--qos=normal','--job-name=J-'+mode,
        '--cpus-per-task=1','--ntasks=1','--mem=%dM'%memory_mb,'--time=%02d:%02d:00'%divmod(minutes,60),
        '--array=0-%d'%(len(indices)-1), '--export=ALL,J_MODE=%s,J_BATCH=%s,J_CODE=%s'%(mode,path,CODE),
        '--output='+str(W/'logs'/('%s-%%A_%%a.out'%mode)), '--error='+str(W/'logs'/('%s-%%A_%%a.err'%mode))]
    if dependency: command.append('--dependency=afterok:'+dependency)
    command.append(str(CODE/'worker.sbatch'))
    job = subprocess.check_output(command,text=True).strip().split(';')[0]
    return dict(job_id=job,mode=mode,indices=indices,minutes=minutes,memory_mb=memory_mb,command=command)


def pilot():
    m=json.loads((W/'search_manifest.json').read_text())
    indices=[r['index'] for r in m['rows'] if r['dm']=='24.40' and r['segments'] in [16,32,64]]
    # Comparable existing full searches requested 10 minutes and 4 GiB. These
    # three jobs measure this exact search-plus-fit before production requests.
    job=submit('search',indices,5,4096,label='pilot')
    save(W/'submitted_pilot.json',dict(created_utc=now(),jobs=[job],maximum_concurrent_tasks=None))
    print(json.dumps(job),flush=True)


def production():
    m=json.loads((W/'search_manifest.json').read_text());p=json.loads((W/'submitted_pilot.json').read_text())
    indices=p['jobs'][0]['indices']; measurements=[]
    for i in indices:
        c=json.loads((Path(m['rows'][i]['directory'])/'search_complete.json').read_text())
        assert c['timing']['returncode']==0;measurements.append(c['timing'])
    seconds=max(r['wall_s'] for r in measurements);rss=max(r['max_rss_kib'] for r in measurements)/1024
    minutes=max(2,math.ceil((seconds*1.25+20)/60));memory=max(1536,math.ceil(rss*1.25/256)*256)
    rest=[r['index'] for r in m['rows'] if r['index'] not in indices]
    jobs=[submit('search',rest,minutes,memory,label='production')]
    save(W/'submitted_search_production.json',dict(created_utc=now(),jobs=jobs,maximum_concurrent_tasks=None,
        pilot_max_wall_s=seconds,pilot_max_rss_mb=rss,safety_factor=1.25))
    print(json.dumps(jobs),flush=True)


def fold_submit():
    fold_prepare()
    m=json.loads((W/'search_manifest.json').read_text())
    assert all((Path(r['directory'])/'search_complete.json').exists() for r in m['rows'])
    # Existing J native PRESTO folds, at the same profile length and full
    # refinement settings, measured about 620 seconds. 15 min is >25% margin.
    jobs=[submit('fold',list(range(2*len(m['rows']))),28,2048,label='production')]
    save(W/'submitted_folds.json',dict(created_utc=now(),jobs=jobs,maximum_concurrent_tasks=None,
        wall_basis='Measured O native fold at 35 bins and identical refinement: 1317 s; 28 minutes includes 25 percent plus startup margin.'))
    print(json.dumps(jobs),flush=True)


def status():
    m=json.loads((W/'search_manifest.json').read_text()); rows=[]
    for r in m['rows']:
        d=Path(r['directory']); complete=d/'search_complete.json'
        for model in ['polynomial','kepler_taylor']:
            marker=d/'native_followup'/model/('presto_'+model+'_f0_f5_complete.json')
            row=dict(index=r['index'],dm=r['dm'],segments=r['segments'],search_complete=complete.exists(),model=model,fold_complete=marker.exists())
            if complete.exists():
                c=json.loads(complete.read_text());row.update(coherence_s=c['sampling']['segment_duration_s'],endpoint_hz=c['endpoint_frequency_hz'],score=c['peak_loglike'])
            if marker.exists():
                c=json.loads(marker.read_text()); headers=c['bestprof_header']
                row.update(sigma=next(float(q[1]) for line in headers if (q:=re.search(r'\(~?\s*([\d.eE+-]+)\s+sigma\)',line))),
                    transfer_passed=c['barycentric_coefficient_transfer']['passed'],wall_s=c['wall_s'],
                    refined_dm=next(line.split('=')[1].strip() for line in headers if line.startswith('# Best DM')),
                    png=next(str(marker.parent/name) for name in c['product_sha256'] if name.endswith('.pfd.png')),
                    pdf=next(str(marker.parent/name) for name in c['product_sha256'] if name.endswith('.pfd.pdf')))
            rows.append(row)
    result=dict(checked_utc=now(),search_complete=sum(r['search_complete'] for r in rows)//2,
                search_expected=len(m['rows']),fold_complete=sum(r['fold_complete'] for r in rows),fold_expected=len(rows),rows=rows)
    dest=W/'status.json';tmp=W/'status.json.tmp';tmp.write_text(json.dumps(result,indent=2)+'\n');tmp.replace(dest)
    with (W/'fold_summary.csv').open('w',newline='') as out:
        fields=list(dict.fromkeys(k for r in rows for k in r));wr=csv.DictWriter(out,fieldnames=fields);wr.writeheader();wr.writerows(rows)
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)
    for r in sorted([r for r in rows if r['fold_complete']],key=lambda r:r['sigma'],reverse=True)[:8]:print(json.dumps(r),flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','pilot','production','search','fold-prepare','fold-submit','fold','status']);p.add_argument('--index',type=int)
    a=p.parse_args()
    if a.mode=='prepare':prepare()
    elif a.mode=='pilot':pilot()
    elif a.mode=='production':production()
    elif a.mode=='search':search(a.index)
    elif a.mode=='fold-prepare':fold_prepare()
    elif a.mode=='fold-submit':fold_submit()
    elif a.mode=='fold':fold(a.index)
    else:status()
