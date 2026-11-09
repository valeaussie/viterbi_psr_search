"""Native channelized follow-up of all saved PRESTO acceleration candidates.

Read-only search inputs, exact binary r/z/w convention, original searched
intervals, TEMPO time conversion, full local refinement. Never overwrite old
folds. All scientific code is preserved on OzSTAR before execution.
"""
from pathlib import Path
from decimal import Decimal
from collections import Counter, defaultdict
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import time
import traceback

E = Path('/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi')
P = E / 'psr_miner_search/hmm_presto_v2_20261027T060742Z/production'
W = P / 'followup/native_presto_results_v1'
CODE = W / 'code'
ROOT = E / '47Tuc_blind_search_v1/stage4_fold/comparison_v1/presto'
FIL = E / '47Tuc_blind_search_v1/stage1_clean/47Tuc_22UL_1of2_L_RFIcleaned_01.fil'
INPUT = P / 'followup/presto_parameters_fullprecision.json'
HELPER = P / 'followup/restored_all_415/code/native_hmm_prepfold_v7.py'
PREFIX = 'presto_native_refined'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def save(path, obj):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(obj, indent=2) + '\n')
    temp.replace(path)


def helper():
    spec = importlib.util.spec_from_file_location('verified_hmm_tempo', HELPER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prepare():
    import numpy as np
    from presto import sigproc, psr_utils
    import presto.presto as presto
    assert not (W / 'manifest.json').exists()
    h = helper()
    data = json.loads(INPUT.read_text())['candidates']
    assert len(data) == len({r['candidate_id'] for r in data}) == 5538
    header, header_bytes = sigproc.read_header(str(FIL))
    dt = float(header['tsamp'])
    assert header['nchans'] == 128 and header['nbits'] == 8 and header['nifs'] == 1
    nsamp = (FIL.stat().st_size-header_bytes)//128
    duration = nsamp*dt
    # Same PRESTO/TEMPO barycenter routine and DE405 clock convention as prepdata.
    # Include midpoints to bound interpolation errors for original section limits.
    seconds = np.linspace(0, duration, 1441)
    top = float(header['tstart']) + seconds / 86400
    bary = np.zeros_like(top)
    velocity = np.zeros_like(top)
    presto.barycenter(top, bary, velocity, '00:24:04.6500', '-72:04:53.8000', 'MK', 'DE405')
    elapsed = (bary-bary[0])*86400
    error = float(np.max(np.abs(np.interp(seconds[1::2], seconds[::2], elapsed[::2])-elapsed[1::2])))
    assert error < dt/4, error
    save(W/'time_conversion.json', dict(topocentric_seconds=seconds.tolist(),
         barycentric_elapsed_seconds=elapsed.tolist(), barycentric_start_mjd=float(bary[0]),
         mean_velocity_c=float(np.mean(velocity)), midpoint_interpolation_error_s=error,
         ephemeris='DE405', observatory='MK', clock='UTC(NIST)',
         explanation='DM-dependent constant high-channel delay cancels in elapsed times; origin verified separately against each prepdata .inf'))
    rows = []
    snapshots = {}
    inf_cache = {}
    binary_hashes = {}
    for i, candidate in enumerate(data):
        cid = candidate['candidate_id']
        target = ROOT / cid / 'acceleration'
        assert target.is_dir(), target
        assert not any(target.glob(PREFIX+'*')), target
        old = [target/name for name in ['fold_arrays.npz','fold.json','fold.png','profiles.txt']]
        old.append(target.parent/'complete.json')
        assert all(p.is_file() for p in old)
        snapshots[cid] = {str(p):dict(size=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,
                                     sha256=sha(p) if p.name=='fold.png' else None) for p in old}
        binary = candidate['source_binary']
        if binary not in binary_hashes:
            binary_hashes[binary] = sha(binary)
        assert binary_hashes[binary] == candidate['source_binary_sha256']
        prop = presto.fourierprops()
        presto.get_rzw_cand(binary, int(candidate['presto_candidate_number']), prop)
        meta = candidate['section_metadata']
        T = float(meta['fft_duration_s'])
        assert abs(T-int(meta['fft_samples'])*dt)<1e-8
        z0 = float(prop.z)-float(prop.w)/2
        r0 = float(prop.r)-z0/2-float(prop.w)/6
        coefficients = [r0/T, z0/T**2, float(prop.w)/T**3]
        for actual, key in zip(coefficients,['f_start_hz','fdot_start_hz_s','fddot_hz_s2']):
            assert math.isclose(actual,float(candidate[key]),rel_tol=1e-13,abs_tol=1e-18)
        infpath = str(Path(meta['source']).with_suffix('.inf'))
        if infpath not in inf_cache:
            info = h.inf_fields(infpath)
            assert h.inf_get(info,'Barycentered?') == '1'
            epoch = h.inf_get(info,'Epoch of observation')
            dm = float(h.inf_get(info,'Dispersion measure'))
            # prepdata uses Doppler-shifted high channel and mean observing velocity.
            fhigh = max(float(header['fch1']),float(header['fch1'])+127*float(header['foff']))
            delay = float(psr_utils.delay_from_DM(dm, fhigh*(1+float(np.mean(velocity)))))
            origin_difference = (float(epoch)-float(bary[0]))*86400+delay
            assert abs(origin_difference)<dt/2, (epoch,origin_difference)
            inf_cache[infpath] = dict(epoch=epoch,dm=dm,origin_difference_s=origin_difference,
                                    ra=h.inf_get(info,'J2000 Right Ascension'),
                                    dec=h.inf_get(info,'J2000 Declination'),sha256=sha(infpath))
        info = inf_cache[infpath]
        assert abs(info['dm']-float(candidate['dm_best']))<1e-8
        start_s = int(meta['start_sample'])*dt
        stop_s = int(meta['end_sample_exclusive'])*dt
        start_top = float(np.interp(start_s,elapsed,seconds))
        stop_top = float(np.interp(stop_s,elapsed,seconds))
        epoch = Decimal(info['epoch'])+Decimal(str(start_s))/Decimal(86400)
        # Upper estimate for default bin count; actual native shape is recorded.
        bins = min(64,max(1,math.ceil(1/(coefficients[0]*dt))))
        rows.append(dict(index=i,candidate=candidate,candidate_id=cid,
             output_directory=str(target),candidate_label='PST%06d'%i,
             coefficient_strings=[format(f,'.17e') for f in coefficients]+['0','0','0'],
             reference_epoch_mjd=str(epoch),input_epoch_mjd=str(epoch),reference_time_s='0',
             dm=str(info['dm']),ra=info['ra'],dec=info['dec'],section_duration_s=stop_s-start_s,
             requested_topocentric_start_s=start_top,requested_topocentric_stop_s=stop_top,
             start_fraction=start_top/duration,end_fraction=stop_top/duration,
             estimated_bins=bins,inf_path=infpath,inf_sha256=info['sha256'],
             clipped_stop_s=max(0,stop_s-float(elapsed[-1]))))
    save(W/'original_files.json',snapshots)
    manifest = dict(rows=rows,candidates=len(rows),filterbank=str(FIL),filterbank_header=header,
                    filterbank_size=FIL.stat().st_size,filterbank_mtime_ns=FIL.stat().st_mtime_ns,
                    filterbank_duration_s=duration,input_sha256=sha(INPUT),helper_sha256=sha(HELPER),
                    inf_verification=inf_cache,full_refinement=dict(npfact=2,ndmfact=3,searchfdd=True,npart=64,nsub=32),
                    original_intervals=True,catalogue_timing_parameters=False,candidate_exclusions=0)
    save(W/'manifest.json',manifest)
    # Reduced-grid timing probe measures the large-bin optimizer without hours of pilot delay.
    # It is stored separately and never counted as a production fold.
    small = min(rows,key=lambda r:(abs(r['estimated_bins']-15),-r['section_duration_s']))
    medium = min(rows,key=lambda r:(abs(r['estimated_bins']-35),r['section_duration_s']))
    large = min((r for r in rows if r['estimated_bins']==64),key=lambda r:r['section_duration_s'])
    pilots=[dict(index=small['index'],timing=False),dict(index=medium['index'],timing=False),
            dict(index=large['index'],timing=True)]
    save(W/'pilots.json',pilots)
    print(json.dumps(dict(candidates=len(rows),bin_counts=dict(Counter(r['estimated_bins'] for r in rows)),
                          pilots=pilots,example=rows[0]['output_directory']),indent=2),flush=True)


def run(index,timing=False):
    import numpy as np
    from presto.prepfold import pfd
    h = helper()
    manifest = json.loads((W/'manifest.json').read_text())
    assert sha(HELPER)==manifest['helper_sha256']
    row = manifest['rows'][index]
    target = (W/'timing'/('%04d'%index)) if timing else Path(row['output_directory'])
    marker = target/(PREFIX+'_complete.json')
    if marker.exists():
        prior=json.loads(marker.read_text())
        assert prior['row']==row
        assert all(sha(target/k)==v for k,v in prior['products'].items())
        print('Already complete',index,flush=True)
        return
    target.mkdir(parents=True,exist_ok=True)
    assert not any(target.glob(PREFIX+'*')), 'Refusing to overwrite native products without completion marker'
    original=json.loads((W/'original_files.json').read_text())[row['candidate_id']]
    assert FIL.stat().st_size==manifest['filterbank_size'] and FIL.stat().st_mtime_ns==manifest['filterbank_mtime_ns']
    assert sha(row['candidate']['source_binary'])==row['candidate']['source_binary_sha256']
    started=time.monotonic()
    scratch=Path(tempfile.mkdtemp(prefix='presto_native_%04d_'%index,dir=os.environ.get('SLURM_TMPDIR','/tmp')))
    record=dict(row=row,timing_probe=timing,job_id=os.environ.get('SLURM_JOB_ID'),code_sha256=sha(__file__))
    try:
        duration=row['section_duration_s']
        bary_start=float(row['reference_epoch_mjd'])
        bary=h.tempo_polycos(scratch/'bary_check',row,'@',0.,0.,bary_start,bary_start+duration/86400)
        record['transfer_check']=h.verify_bary_transfer(bary,row,duration)
        header=manifest['filterbank_header']
        obsfreq=(float(header['fch1'])+float(header['fch1'])+127*float(header['foff']))/2
        top_start=float(header['tstart'])+row['requested_topocentric_start_s']/86400
        top_stop=float(header['tstart'])+row['requested_topocentric_stop_s']/86400
        # Include a record-width margin because prepfold rounds requested boundaries.
        topo=h.tempo_polycos(scratch/'topocentric',row,'m',obsfreq,row['dm'],top_start-2/86400,top_stop+2/86400)
        text=(scratch/'topocentric/polyco_for_observation.dat').read_text()
        text=re.sub(r'(?m)^(\s*)0000\+0000(?=\s)',lambda m:m.group(1)+row['candidate_label'],text)
        (scratch/'model.polycos').write_text(text)
        (scratch/'input.fil').symlink_to(FIL)
        command=['prepfold','-filterbank','-topo','-noxwin','-ncpus','1','-dm',row['dm'],
                 '-npfact','1' if timing else '2','-ndmfact','1' if timing else '3',
                 '-nsub','32','-npart','64','-psr',row['candidate_label'],'-polycos','model.polycos',
                 '-start',format(row['start_fraction'],'.17g'),'-end',format(row['end_fraction'],'.17g'),
                 '-searchfdd','-o',PREFIX,'input.fil']
        record['command']=command
        print('Starting',index,row['candidate_id'],'timing_probe',timing,flush=True)
        native_start=time.monotonic()
        optimization_start=None
        logpath=scratch/'prepfold.log'
        with logpath.open('w') as log:
            process=subprocess.Popen(command,cwd=scratch,stdout=log,stderr=subprocess.STDOUT)
            while process.poll() is None:
                if optimization_start is None and 'Optimizing...' in logpath.read_text(errors='replace'):
                    optimization_start=time.monotonic()
                time.sleep(.5)
        native_stop=time.monotonic()
        record['returncode']=process.returncode
        assert process.returncode==0, 'prepfold failed'
        paths=list(scratch.glob(PREFIX+'*.pfd'))
        assert len(paths)==1
        folded=pfd(str(paths[0]))
        assert folded.nsub==32 and folded.npart==64
        assert np.isfinite(folded.profs).all()
        record['profile_shape']=list(folded.profs.shape)
        record['folded_topocentric_epoch_mjd']=float(folded.tepoch)
        record['folded_samples']=int(folded.stats[:,:,0].sum()/folded.nsub)
        record['native_wall_s']=native_stop-native_start
        record['raw_fold_s']=(optimization_start or native_stop)-native_start
        record['optimization_and_render_s']=native_stop-(optimization_start or native_stop)
        record['optimizer_start_detected']=optimization_start is not None
        logtext=logpath.read_text(errors='replace')
        match=re.search(r'Searching (\d+) DMs, (\d+) periods, (\d+) p-dots, and (\d+) p-dotdots',logtext)
        assert match, 'Cannot verify refinement grid'
        record['trial_dimensions']=[int(v) for v in match.groups()]
        record['trial_count']=math.prod(record['trial_dimensions'])
        record['bestprof_header']=[l for l in Path(str(paths[0])+'.bestprof').read_text().splitlines() if l.startswith('#')]
        h.render_ps(Path(str(paths[0])+'.ps'))
        products={PREFIX+'.pfd':paths[0]}
        for suffix in ['.bestprof','.ps','.pdf','.png']:
            products[PREFIX+'.pfd'+suffix]=Path(str(paths[0])+suffix)
        products.update({PREFIX+'.par':scratch/'topocentric/model.par',
             PREFIX+'.polycos':scratch/'model.polycos',PREFIX+'_prepfold.log':logpath,
             PREFIX+'_tempo.log':scratch/'topocentric/tempo.log',
             PREFIX+'_bary_check.polycos':scratch/'bary_check/polyco_for_observation.dat',
             PREFIX+'_tempo_raw.polycos':scratch/'topocentric/polyco.dat'})
        for name,path in products.items():
            with (target/name).open('xb') as dest,path.open('rb') as source:
                shutil.copyfileobj(source,dest)
        for name,before in original.items():
            path=Path(name)
            assert path.stat().st_size==before['size'] and path.stat().st_mtime_ns==before['mtime_ns']
            if before['sha256']:
                assert sha(path)==before['sha256']
        record['products']={name:sha(target/name) for name in products}
        record['wall_s']=time.monotonic()-started
        record['original_files_unchanged']=True
        record['original_search_interval_used']=True
        save(marker,record)
        print('Completed',index,'wall_s',record['wall_s'],'shape',record['profile_shape'],flush=True)
        shutil.rmtree(scratch)
    except Exception:
        record['error']=traceback.format_exc()
        record['scratch']=str(scratch)
        error=W/('failure_%04d%s'%(index,'_timing' if timing else ''))
        error.mkdir(exist_ok=True)
        for path in scratch.rglob('*'):
            if path.is_file() and not path.is_symlink() and path.stat().st_size<20_000_000:
                dest=error/path.relative_to(scratch)
                dest.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(path,dest)
        save(error/'failure.json',record)
        raise


def pilots():
    plan=json.loads((W/'pilots.json').read_text())
    env=os.environ.copy()
    jobs=[]
    for item in plan:
        i=item['index']
        env['PRESTO_INDEX']=str(i)
        env['PRESTO_TIMING']='1' if item['timing'] else '0'
        command=['sbatch','--parsable','--account=oz022','--export=ALL','--time=00:30:00',
                 '--output='+str(W/('logs/pilot-%04d-%%j.out'%i)),
                 '--error='+str(W/('logs/pilot-%04d-%%j.err'%i)),str(CODE/'worker_v2.sbatch')]
        job=subprocess.check_output(command,env=env,text=True).strip().split(';')[0]
        jobs.append(dict(**item,job_id=job,command=command))
        save(W/'submitted_pilots.json',dict(jobs=jobs))
    print(json.dumps(jobs,indent=2))


def submit():
    assert not (W/'submitted_production.json').exists()
    manifest=json.loads((W/'manifest.json').read_text())
    controls=[]
    for item in json.loads((W/'pilots.json').read_text()):
        row=manifest['rows'][item['index']]
        folder=W/'timing'/('%04d'%item['index']) if item['timing'] else Path(row['output_directory'])
        control=json.loads((folder/(PREFIX+'_complete.json')).read_text())
        assert control['transfer_check']['passed'] and control['optimizer_start_detected']
        controls.append(control)
    calibration=max(c['optimization_and_render_s']/(c['trial_count']*c['profile_shape'][-1]) for c in controls)
    raw_rate=max(c['raw_fold_s']/c['row']['section_duration_s'] for c in controls)
    raw_overhead=max(c['wall_s']-c['native_wall_s'] for c in controls)
    groups=defaultdict(list)
    for row in manifest['rows']:
        if (Path(row['output_directory'])/(PREFIX+'_complete.json')).exists():
            continue
        n=row['estimated_bins']
        trials=(6*n+1)*(4*n+1)**3
        expected=raw_overhead+raw_rate*row['section_duration_s']+calibration*trials*n
        minutes=max(10,math.ceil(expected*1.25/300)*5)
        groups[minutes].append(row['index'])
    assoc=subprocess.check_output(['sacctmgr','-nP','show','assoc','where','account=oz022','user=vdimarco',
          'format=Account,User,QOS,DefaultQOS,GrpTRES,MaxJobs,MaxSubmitJobs'],text=True)
    qos=subprocess.check_output(['sacctmgr','-nP','show','qos','where','name=normal',
          'format=Name,MaxTRESPU,MaxJobsPU,MaxSubmitJobsPU,MaxWall'],text=True)
    state=dict(jobs=[],complete=False,calibration_seconds_per_trial_bin=calibration,
               raw_seconds_per_observation_second=raw_rate,overhead_seconds=raw_overhead,
               association=assoc,qos=qos,margin=1.25,maximum_concurrent_tasks=None,
               expected_candidates=5538,old_folds_preserved=True,
               official_docs=['https://supercomputing.swin.edu.au/docs/1a-nt_oz/what-s-new.html',
                              'https://supercomputing.swin.edu.au/docs/3-jobmon/monitoring.html'])
    save(W/'submitted_production.json',state)
    env=os.environ.copy()
    env.pop('PRESTO_INDEX',None)
    env['PRESTO_TIMING']='0'
    # Cluster MaxArraySize is 2048. Remap index through manifest batches;
    # arrays have at most 2000 tasks, with no % concurrency restriction.
    for minutes,indices in sorted(groups.items()):
        for offset in range(0,len(indices),2000):
            batch=indices[offset:offset+2000]
            batchname='batch_%04d.json'%len(state['jobs'])
            save(W/batchname,batch)
            env['PRESTO_BATCH']=str(W/batchname)
            limit='%02d:%02d:00'%divmod(minutes,60)
            command=['sbatch','--parsable','--account=oz022','--export=ALL','--time='+limit,
                 '--array=0-%d'%(len(batch)-1),'--output='+str(W/'logs/full-%A_%a.out'),
                 '--error='+str(W/'logs/full-%A_%a.err'),str(CODE/'worker_v2.sbatch')]
            job=subprocess.check_output(command,env=env,text=True).strip().split(';')[0]
            state['jobs'].append(dict(job_id=job,minutes=minutes,indices=batch,command=command))
            save(W/'submitted_production.json',state)
    deps=':'.join(j['job_id'] for j in state['jobs'])
    command=['sbatch','--parsable','--account=oz022','--export=ALL','--dependency=afterany:'+deps,
             '--kill-on-invalid-dep=yes','--output='+str(W/'logs/report-%j.out'),
             '--error='+str(W/'logs/report-%j.err'),str(CODE/'report_v2.sbatch')]
    state['report_job_id']=subprocess.check_output(command,env=env,text=True).strip().split(';')[0]
    state['complete']=True
    save(W/'submitted_production.json',state)
    print(json.dumps(dict(jobs=[dict(job_id=j['job_id'],tasks=len(j['indices']),minutes=j['minutes']) for j in state['jobs']],
                          report=state['report_job_id'],candidates=5538,maximum_concurrent_tasks=None),indent=2))


def report():
    manifest=json.loads((W/'manifest.json').read_text())
    results=[]
    for row in manifest['rows']:
        marker=Path(row['output_directory'])/(PREFIX+'_complete.json')
        error=W/('failure_%04d'%row['index'])/'failure.json'
        if marker.exists():
            data=json.loads(marker.read_text())
            results.append(dict(index=row['index'],candidate_id=row['candidate_id'],status='complete',
                 wall_s=data['wall_s'],pdf=str(marker.parent/(PREFIX+'.pfd.pdf')),
                 bestprof_header=data['bestprof_header']))
        else:
            results.append(dict(index=row['index'],candidate_id=row['candidate_id'],status='failed' if error.exists() else 'missing',
                                error=json.loads(error.read_text()).get('error') if error.exists() else None))
    old=json.loads((W/'original_files.json').read_text())
    changed=[]
    for files in old.values():
        for name,before in files.items():
            p=Path(name)
            if not p.is_file() or p.stat().st_size!=before['size'] or p.stat().st_mtime_ns!=before['mtime_ns'] or (before['sha256'] and sha(p)!=before['sha256']):
                changed.append(name)
    data=dict(counts=dict(Counter(r['status'] for r in results)),original_files_changed=changed,results=results)
    save(W/'completion_report.json',data)
    print(json.dumps({k:v for k,v in data.items() if k!='results'},indent=2))
    assert not changed


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=['prepare','run','pilots','submit','report'])
    parser.add_argument('--index',type=int)
    parser.add_argument('--timing',action='store_true')
    args=parser.parse_args()
    W.mkdir(parents=True,exist_ok=True)
    (W/'logs').mkdir(exist_ok=True)
    if args.mode=='run':
        index=args.index
        if index is None:
            index=int(os.environ['SLURM_ARRAY_TASK_ID'])
            if os.environ.get('PRESTO_BATCH'):
                index=json.loads(Path(os.environ['PRESTO_BATCH']).read_text())[index]
        run(index,args.timing)
    else:
        globals()[args.mode]()
