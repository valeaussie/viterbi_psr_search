"""Supplementary Nt16 endpoints, exact track recovery and native folds.

The uniform score cut is exploratory. Catalogue probes label review priorities
after selection; they never alter selection or supply folding coefficients.
"""
from pathlib import Path
from collections import defaultdict,Counter
from decimal import Decimal
import argparse,bisect,csv,hashlib,html,importlib.util,json,math,os,shutil,subprocess,sys,time
import numpy as np
W=Path(__file__).resolve().parents[1];P=W.parents[1];E=P.parents[2]
A=P/'followup/lower_threshold_audit_v1';B=E/'47Tuc_blind_search_v1/stage3_viterbi/blind_v2'
F=P/'frozen/code';C=W/'candidates';D=E/'47Tuc_blind_search_v1/stage2_dedisp'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def save(p,x):
    with Path(p).open('x') as out:json.dump(x,out,indent=2);out.write('\n')
def rows(p):return list(csv.DictReader(Path(p).open()))
def table(p,rr):
    with Path(p).open('x',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=list(dict.fromkeys(k for r in rr for k in r)));writer.writeheader();writer.writerows(rr)
def module(path,name):
    sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
def check_inputs():
    for path,digest in load(A/'inputs.json')['source_hashes'].items():assert sha(path)==digest,path
def prepare():
    check_inputs();C.mkdir();(W/'source_cache').mkdir()
    cfg=load(A/'inputs.json');tasks=[t for t in cfg['tasks'] if t['Nt']==16]
    inv={r['source_curve']:r for r in rows(A/'Nt16/curve_inventory.csv')}
    with np.load(A/'Nt16/peaks_above_median.npz') as archive:
        z={k:archive[k] for k in ['score','frequency_hz','curve_index','source_bin']}
    mask=z['score']>58
    assert np.count_nonzero(mask)==4021
    anchors=[];ids=[];groups=[];members=[]
    for j in np.flatnonzero(mask):
        frequency=float(z['frequency_hz'][j]);score=float(z['score'][j]);t=tasks[int(z['curve_index'][j])]
        lo=bisect.bisect_left(anchors,frequency-.035);hi=bisect.bisect_right(anchors,frequency+.035)
        if lo<hi:gi=min(ids[lo:hi])
        else:
            gi=len(groups);where=bisect.bisect_left(anchors,frequency);anchors.insert(where,frequency);ids.insert(where,gi)
            groups.append(dict(group_id=f'Nt16_family{gi:06d}',representative_hz=frequency,representative_score=score,
                dm=t['dm'],subband_f0=t['subband_f0'],source_curve=t['curve'],source_bin=int(z['source_bin'][j]),
                original=score>cfg['thresholds']['16']))
        members.append(dict(group_id=groups[gi]['group_id'],source_curve=t['curve'],source_bin=int(z['source_bin'][j]),
            frequency_hz=frequency,score=score,dm=t['dm'],subband_f0=t['subband_f0']))
    assert len(groups)==126 and sum(g['original'] for g in groups)==38
    selected=[g for g in groups if not g['original']];assert len(selected)==88
    probes=rows(A/'Nt16_target_inspection.csv');priorities={}
    for target,harmonic in [('C',1),('C',2),('E',1),('L',3)]:
        r=next(r for r in probes if r['target']==target and int(r['harmonic'])==harmonic)
        priorities[(r['source_curve'],int(r['window_peak_bin']))]=f'{target}_h{harmonic}'
    selected.sort(key=lambda g:(0 if (g['source_curve'],g['source_bin']) in priorities else 1,-g['representative_score']))
    old=load(P/'followup/track_manifest.json');cached={g['source_curve']:i for i,g in enumerate(old['sources'])}
    output=[];sources={}
    for i,g in enumerate(selected):
        cid=f'supp16_{i:03d}_DM{g["dm"]:.2f}_f{g["representative_hz"]:.6f}Hz'
        subset=[m for m in members if m['group_id']==g['group_id']]
        row=dict(candidate_id=cid,cand_index=i,group_id=g['group_id'],source_curve=g['source_curve'],
            source_params=g['source_curve'].replace('_loglike_curve.dat','.params'),source_bin=g['source_bin'],
            peak_freq_hz=g['representative_hz'],peak_loglike=g['representative_score'],dm_best=g['dm'],
            Nt_best=16,subband_f0=g['subband_f0'],threshold=58.,original_threshold=cfg['thresholds']['16'],
            priority_probe=priorities.get((g['source_curve'],g['source_bin']),''),classification='supplementary_unconfirmed',
            parent_pulsar='',rfi_flag='',multiplicity=len(subset),dm_values=';'.join(f'{v:.2f}' for v in sorted({m['dm'] for m in subset})),
            dm_count=len({m['dm'] for m in subset}))
        output.append(row)
        if row['source_curve'] not in sources:
            sources[row['source_curve']]=dict(source_curve=row['source_curve'],source_params=row['source_params'],
                original_sha256=inv[row['source_curve']]['sha256'],candidates=[],old_cache=(str(P/'followup/source_cache'/f'source_{cached[row["source_curve"]]:03d}') if row['source_curve'] in cached else None))
        sources[row['source_curve']]['candidates'].append(row)
    assert sum(bool(r['priority_probe']) for r in output)==4
    table(W/'all_Nt16_groups_at_58.csv',groups);table(W/'all_Nt16_members_at_58.csv',members)
    table(W/'supplementary_candidates.csv',output)
    table(W/'supplementary_members.csv',[m for m in members if m['group_id'] in {g['group_id'] for g in selected}])
    save(W/'selection.json',dict(cut=58.,Nt=16,all_groups=126,original_groups=38,new_groups=88,raw_peaks=4021,
        source_runs=len(sources),sources=list(sources.values()),candidates=output,original_input_hashes=cfg['source_hashes'],
        selection='Uniform exploratory cut over all searched DM/subbands; original frequency grouping, target-independent selection.',
        alternatives='All source-curve/bin members retained in supplementary_members.csv; multiplicity is correlated.'))
    save(W/'code_manifest.json',{p.name:sha(p) for p in (W/'code').iterdir() if p.is_file()})
    print(json.dumps(dict(new_candidates=88,source_runs=len(sources),reusable_sources=sum(s['old_cache'] is not None for s in sources.values()),priorities=[r['priority_probe'] for r in output[:4]])))
def recover(index):
    started=time.monotonic();sel=load(W/'selection.json');s=sel['sources'][index]
    cache=W/'source_cache'/f'source_{index:03d}';cache.mkdir(exist_ok=True)
    if (cache/'complete.json').exists():print('Already recovered',index);return
    helper=module(W/'code/recover_helper.py','recover_helper')
    assert sha(B/s['source_curve'])==s['original_sha256']
    if s['old_cache']:
        original=Path(s['old_cache']);assert (original/'complete.json').is_file()
        for suffix in ['_loglike_curve.dat','_delta.dat','_sampling.json']:
            shutil.copyfile(original/('source'+suffix),cache/('source'+suffix))
    else:
        values=helper.load_params(B/s['source_params']);values.update(out_prefix=str(cache/'source'),save_delta=True)
        params=cache/'source.params';params.write_text('\n'.join(f'{k} = {v!r}' for k,v in values.items())+'\n')
        with (cache/'recovery.log').open('w') as log:
            subprocess.run([sys.executable,str(F/'viterbi_pipeline.py'),'--params',str(params),'--lean'],stdout=log,stderr=subprocess.STDOUT,check=True)
    curve=np.loadtxt(cache/'source_loglike_curve.dat');old=np.loadtxt(B/s['source_curve']);delta=np.loadtxt(cache/'source_delta.dat')
    np.testing.assert_allclose(curve[:,0],old[:,0],rtol=0,atol=1e-10)
    np.testing.assert_allclose(curve[:,1],old[:,1],rtol=0,atol=1e-6)
    np.testing.assert_allclose(delta[:,-1],old[:,1],rtol=0,atol=1e-6)
    sample=load(cache/'source_sampling.json');assert delta.shape==(len(curve),16)
    tsft=sample['segment_duration_s'];spec=helper.emissions_from_delta(delta);assert np.min(spec)>-1e-7
    np.savez_compressed(cache/'original_search_state.npz',delta=delta,spectrogram=spec,frequency_hz=curve[:,0],segment_duration_s=tsft)
    t=(np.arange(16)+.5)*tsft
    sys.path.insert(0,str(F));fit=module(F/'viterbi_search_and_fit.py','fit_helpers')
    for row in s['candidates']:
        out=C/row['candidate_id'];out.mkdir(exist_ok=True);end=row['source_bin']
        assert abs(curve[end,0]-row['peak_freq_hz'])<1e-9 and abs(delta[end,-1]-row['peak_loglike'])<1e-6
        indices=helper.backtrace_delta(delta,end);freq=curve[indices,0]
        assert np.max(np.abs(np.diff(indices)))<=1
        np.testing.assert_allclose(spec[indices,np.arange(16)].sum(),row['peak_loglike'],rtol=0,atol=1e-6)
        np.savetxt(out/(out.name+'_track.dat'),np.c_[t,freq]);np.save(out/'original_state_indices.npy',indices)
        np.savetxt(out/'selected_source_curve.dat',old)
        save(out/'candidate.json',dict(row,exact_original_path=True,sampling=sample,source_cache=str(cache),selected_candidate_retained=True,
            source_curve_max_score_difference=float(np.max(np.abs(curve[:,1]-old[:,1])))))
        helper.plot_candidate(freq,curve[:,0],spec,tsft,row,out)
        df=float(np.median(np.diff(curve[:,0])));poly=fit.fit_polynomial_at_midpoint(t,freq,5,df);kep=fit.fit_kepler(t,freq,df,tsft)
        kcoeff=fit.kepler_derivs(poly[0],*kep[0],5)
        models=dict(reference_time_s=poly[0],polynomial=poly[1].tolist(),kepler_taylor=np.asarray(kcoeff).tolist(),
            polynomial_order=poly[7],polynomial_chi2=poly[3],kepler_parameters=kep[0].tolist(),kepler_chi2=kep[2],
            caveat='Models fitted only to the original selected blind track; no catalogue frequency/ephemeris supplied.')
        save(out/'models.json',models)
        summary=out/'fit_summary.dat';summary.write_text('\n'.join([f'reference_time_s {poly[0]:.15g}']+[f'{name}_f{k} {v:.17g}' for name in ['polynomial','kepler_taylor'] for k,v in enumerate(models[name])])+'\n')
        np.savetxt(out/'track_with_fits.dat',np.c_[t,freq,poly[6],kep[5]])
    save(cache/'complete.json',dict(source_curve=s['source_curve'],candidate_ids=[r['candidate_id'] for r in s['candidates']],
        max_score_difference=float(np.max(np.abs(curve[:,1]-old[:,1]))),wall_s=time.monotonic()-started,job_id=os.environ.get('SLURM_JOB_ID'),reused=s['old_cache']))
    print('Recovered',index,'wall',time.monotonic()-started,flush=True)
def manifest():
    check_inputs();sel=load(W/'selection.json');old=load(P/'followup/targeted_native_v1/manifest.json')
    h=module(W/'code/native_helper.py','native_prepare');rr=[];snapshots={}
    for n,c in enumerate(sel['candidates']):
        out=C/c['candidate_id'];info=load(out/'candidate.json');models=load(out/'models.json');summary=out/'fit_summary.dat'
        dm=c['dm_best'];inf=D/f'DM{dm:.2f}_bary'/f'47Tuc_DM{dm:.2f}_bary.inf';v=h.inf_fields(inf);epoch=h.inf_get(v,'Epoch of observation')
        assert h.inf_get(v,'Barycentered?')=='1';tref=models['reference_time_s']
        snapshots[out.name]={str(p):dict(sha256=sha(p),mtime_ns=p.stat().st_mtime_ns) for p in [summary,out/'candidate.json',out/(out.name+'_track.dat')]}
        for name in ['polynomial','kepler_taylor']:
            coeff=models[name];finite=all(math.isfinite(x) for x in coeff)
            grid=np.linspace(0,old['filterbank_metadata']['duration_s'],1001)-tref
            freq=sum(x*grid**k/math.factorial(k) for k,x in enumerate(coeff))
            valid=finite and float(np.min(freq))>0
            native=out/'native'/name;native.mkdir(parents=True,exist_ok=True)
            row=dict(index=len(rr),candidate_number=n,candidate_id=out.name,candidate_json=str(out/'candidate.json'),
                summary=str(summary),summary_sha256=sha(summary),inf=str(inf),inf_sha256=sha(inf),model=name,
                coefficient_strings=[f'{x:.17g}' for x in coeff],reference_epoch_mjd=str(Decimal(epoch)+Decimal(str(tref))/Decimal(86400)),
                input_epoch_mjd=epoch,reference_time_s=str(tref),dm=f'{dm:.2f}',ra=h.inf_get(v,'J2000 Right Ascension'),
                dec=h.inf_get(v,'J2000 Declination'),source_valid_samples=info['sampling']['valid_samples'],output_directory=str(native),
                candidate_label=f'SUP{n:06d}',classification=c['classification'],rfi_flag='',priority_probe=c['priority_probe'],
                preflight_positive_frequency=valid,estimated_bins=min(64,max(1,math.ceil(1/((coeff[0] if finite and coeff[0]>0 else c['peak_freq_hz'])*old['filterbank_metadata']['header']['tsamp'])))))
            rr.append(row)
    save(W/'manifest.json',dict(rows=rr,filterbank=old['filterbank'],filterbank_metadata=old['filterbank_metadata'],settings=old['settings'],
        selection_threshold=58,segment_count=16,all_coefficients_from_blind_tracks=True))
    save(W/'original_fold_files.json',snapshots)
    print('Prepared native rows',len(rr),'valid',sum(r['preflight_positive_frequency'] for r in rr),flush=True)
def native(index):
    from presto import psr_utils
    def rotate(arr,bins):
        arr=np.asarray(arr);f=np.arange(arr.size//2+1,dtype=float)
        return np.fft.irfft(np.exp(2j*np.pi*f*bins/arr.size)*np.fft.rfft(arr),n=arr.size)
    psr_utils.fft_rotate=rotate
    h=module(W/'code/native_helper.py','native_run');h.WORKFLOW=W;h.run(index)
def quick(index):
    m=load(W/'selection.json');c=m['candidates'][index];out=C/c['candidate_id'];models=load(out/'models.json');info=load(out/'candidate.json')
    h=module(W/'code/common_fold.py','inspection_fold')
    source=D/f'DM{c["dm_best"]:.2f}_bary'/f'47Tuc_DM{c["dm_best"]:.2f}_bary.dat';x=np.memmap(source,dtype='<f4',mode='r')
    for name in ['polynomial','kepler_taylor']:
        dest=out/'inspection'/name
        if (dest/'fold.json').exists():continue
        if not all(math.isfinite(x) for x in models[name]):
            dest.mkdir(parents=True,exist_ok=True);save(dest/'invalid_model.json',dict(reason='Nonfinite track-fit coefficients; no replacement coefficients supplied.'));continue
        h.fold_model(x,info['sampling']['tsamp_s'],0,info['sampling']['valid_samples'],
            dict(kind='polynomial',reference_time_s=models['reference_time_s'],frequency_derivatives=models[name]),dest)
def submit(mode,indices,minutes,label,dependency=None):
    assert not (W/f'submitted_{label}.json').exists()
    batch=W/f'{label}_batch.json'
    if batch.exists():assert load(batch)==indices
    else:save(batch,indices)
    cmd=[os.environ.get('SUPP_SBATCH','sbatch'),'--parsable','--account=oz022','--qos=normal','--ntasks=1','--cpus-per-task=1','--mem='+('1G' if mode=='native' else '4G'),
        '--time=%02d:%02d:00'%divmod(minutes,60),'--job-name='+label,'--chdir='+str(W),'--array=0-'+str(len(indices)-1),
        '--output='+str(W/f'logs/{label}-%A_%a.out'),'--error='+str(W/f'logs/{label}-%A_%a.err'),
        '--export=ALL,SUPP_ROOT='+str(W)+',SUPP_MODE='+mode+',SUPP_BATCH='+str(batch),str(W/'code/worker.sbatch')]
    if dependency:cmd[-1:-1]=['--dependency='+dependency,'--kill-on-invalid-dep=yes']
    job=subprocess.check_output(cmd,text=True).strip().split(';')[0];record=dict(job_id=job,mode=mode,indices=indices,minutes=minutes,command=cmd)
    save(W/f'submitted_{label}.json',record);print(json.dumps(record),flush=True)
    return job
def launch_recovery(pilot=False):
    sel=load(W/'selection.json');pending=[i for i,s in enumerate(sel['sources']) if not (W/'source_cache'/f'source_{i:03d}'/'complete.json').exists()]
    if pilot:
        i=next(i for i in pending if any(c['priority_probe']=='C_h1' for c in sel['sources'][i]['candidates']))
        # Existing production search timings provide an initial, measured bound.
        p=B/sel['sources'][i]['source_curve'];timing=load(str(p).replace('_loglike_curve.dat','_timing.json'))
        if not (W/'original_timing_pilot.json').exists():save(W/'original_timing_pilot.json',timing)
        submit('recover',[i],5,'recovery-pilot')
    else:
        pilot=next(load(p) for p in (W/'source_cache').glob('*/complete.json'))
        minutes=max(2,math.ceil((pilot['wall_s']*1.25+30)/60))
        return submit('recover',pending,minutes,'recovery-production')
def chain():
    pilot=load(W/'submitted_recovery-pilot-v2.json')['job_id']
    submit('recover-chain',[0],2,'recovery-chain','afterok:'+pilot)
def recover_chain():
    job=launch_recovery(False)
    submit('finish-recovery',[0],3,'prepare-folds','afterok:'+job)
    submit('report',[0],2,'recovery-report','afterany:'+job)
def finish_recovery():
    manifest()
    quickjob=submit('quick',list(range(4)),5,'control-inspection')
    launch_native(True)
    submit('report',[0],2,'inspection-report','afterany:'+quickjob)
def launch_native(controls=False):
    m=load(W/'manifest.json');cal=load(P/'followup/native_presto_results_v1/submitted_production.json');duration=m['filterbank_metadata']['duration_s'];groups=defaultdict(list)
    completed=[]
    if not controls:
        for row in m['rows']:
            if not row['priority_probe']:continue
            marker=Path(row['output_directory'])/('presto_'+row['model']+'_f0_f5_complete.json')
            rec=load(marker);completed.append(dict(bins=rec['native_profile_shape'][-1],wall_s=rec['wall_s']))
        save(W/'remaining_runtime_calibration.json',dict(measurements=completed,safety_margin=1.25,
            rule='Use the nearest measured bin count at or above the requested count; scale optimizer work by trials times bins, retain measured raw-data baseline. Five-minute rounding.',
            original_calibration=str(P/'followup/native_presto_results_v1/submitted_production.json')))
    for r in m['rows']:
        if not r['preflight_positive_frequency'] or bool(r['priority_probe'])!=controls:continue
        n=r['estimated_bins'];trials=(6*n+1)*(4*n+1)**3
        prediction=cal['overhead_seconds']+cal['raw_seconds_per_observation_second']*duration+cal['fixed_optimizer_seconds']+cal['calibration_seconds_per_trial_bin']*trials*n
        if completed:
            reference=min(r['bins'] for r in completed if r['bins']>=n)
            wall=max(r['wall_s'] for r in completed if r['bins']==reference)
            baseline=cal['overhead_seconds']+cal['raw_seconds_per_observation_second']*duration+cal['fixed_optimizer_seconds']
            reference_work=(6*reference+1)*(4*reference+1)**3*reference
            prediction=baseline+max(0,wall-baseline)*(trials*n)/reference_work
        minutes=max(10,math.ceil(prediction*1.25/300)*5);groups[minutes].append(r['index'])
    jobs=[]
    for minutes,indices in sorted(groups.items()):jobs.append(submit('native',indices,minutes,('control' if controls else 'remaining')+f'-native-{minutes}'))
    if not controls:submit('report',[0],5,'final-report','afterany:'+':'.join(jobs))
def report():
    check_inputs();sel=load(W/'selection.json');manifest=load(W/'manifest.json') if (W/'manifest.json').exists() else {'rows':[]}
    rr=[]
    for r in manifest['rows']:
        out=Path(r['output_directory']);prefix='presto_'+r['model']+'_f0_f5';marker=out/(prefix+'_complete.json')
        rec=load(marker) if marker.exists() else None
        failure=W/f'failure_{r["index"]:04d}_v7.json'
        status='complete' if rec else ('invalid_model' if not r['preflight_positive_frequency'] else ('failed' if failure.exists() else 'pending'))
        if rec:
            assert rec['row']==r
            assert all(sha(out/n)==v for n,v in rec['product_sha256'].items())
        rr.append(dict(candidate_id=r['candidate_id'],priority_probe=r['priority_probe'],model=r['model'],status=status,
            png=next((str(out/n) for n in rec['product_sha256'] if n.endswith('.pfd.png')),'') if rec else '',
            pdf=next((str(out/n) for n in rec['product_sha256'] if n.endswith('.pfd.pdf')),'') if rec else '',
            significance_header='; '.join(s for s in rec['bestprof_header'] if 'Prob(Noise)' in s) if rec else ''))
    version=int(time.time());result=dict(candidate_groups=88,recovered_sources=sum((W/'source_cache'/f'source_{i:03d}'/'complete.json').exists() for i in range(len(sel['sources']))),
        expected_sources=len(sel['sources']),native_fold_status=dict(Counter(r['status'] for r in rr)),output_root=str(W),original_inputs_unchanged=True)
    save(W/f'status_{version}.json',result)
    if rr:table(W/f'fold_status_{version}.csv',rr)
    # A persistent, filterable index links original endpoint plots and any folds.
    cards=[]
    for c in sel['candidates']:
        cid=c['candidate_id'];rel='candidates/'+cid+'/'
        links=[f'<a href="{html.escape(rel+cid+"_spectrogram.png")}">Track and score</a>']
        for name in ['polynomial','kepler_taylor']:
            p=C/cid/'inspection'/name/'fold.png'
            if p.exists():links.append(f'<a href="{p.relative_to(W)}">{name} inspection fold</a>')
        for r in rr:
            if r['candidate_id']==cid and r['pdf']:links.append(f'<a href="{Path(r["pdf"]).relative_to(W)}">{r["model"]} native fold</a>')
        cards.append('<article><h3>'+html.escape(cid+' '+c['priority_probe'])+'</h3><p>Score %.3f; DM %.2f; %d correlated members.</p>'%(c['peak_loglike'],c['dm_best'],c['multiplicity'])+' | '.join(links)+'</article>')
    (W/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Nt16 supplementary candidates</title><style>body{font:16px system-ui;max-width:1100px;margin:25px}article{border-top:1px solid #ddd;padding:12px}h3{font-size:16px}</style><h1>88 supplementary Nt16 candidates</h1><p>Uniform exploratory score cut: 58. Catalogue labels are inspection leads. Original blind results preserved.</p>'+''.join(cards))
    print(json.dumps(result),flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('mode');a.add_argument('--index',type=int);args=a.parse_args()
    if args.mode=='prepare':prepare()
    elif args.mode=='recover':recover(args.index)
    elif args.mode=='manifest':manifest()
    elif args.mode=='native':native(args.index)
    elif args.mode=='quick':quick(args.index)
    elif args.mode=='pilot':launch_recovery(True)
    elif args.mode=='retry-pilot':submit('recover',[2],5,'recovery-pilot-v2')
    elif args.mode=='recover-all':launch_recovery(False)
    elif args.mode=='chain':chain()
    elif args.mode=='recover-chain':recover_chain()
    elif args.mode=='finish-recovery':finish_recovery()
    elif args.mode=='controls':launch_native(True)
    elif args.mode=='remaining':launch_native(False)
    elif args.mode=='quick-controls':submit('quick',list(range(4)),5,'control-inspection')
    elif args.mode=='report':report()
    else:raise ValueError(args.mode)
