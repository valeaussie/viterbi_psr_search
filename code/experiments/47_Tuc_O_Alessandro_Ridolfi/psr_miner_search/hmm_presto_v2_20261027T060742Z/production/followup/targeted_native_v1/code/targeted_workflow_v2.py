"""Native targeted HMM follow-up, preserving saved searches and all coefficients.

Likelihood thresholds are diagnostic applications of the historical tail
prescription, not calibrated false-alarm probabilities. No HMM search is rerun.
"""
from pathlib import Path
from decimal import Decimal
from collections import defaultdict,Counter
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
import time

W=Path(__file__).resolve().parents[1];P=W.parents[1];E=P.parents[2]
ROOT=E/'47Tuc_blind_search_v1'
HELPER=P/'followup/restored_all_415/code/native_hmm_prepfold_v7.py'
AGG=P/'frozen/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/aggregate_blind.py'
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')
def loadmodule(path,name):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
def fields(path):
    out={}
    for line in Path(path).read_text().splitlines():
        v=line.split()
        if v and not v[0].startswith('#'):out[v[0]]=v[1:]
    return out
def csvsave(path,rows):
    with Path(path).open('x',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)

def make_plot(row):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    start=time.monotonic();data=np.loadtxt(row['curve'])
    destination=Path(row['output']);destination.mkdir(parents=True,exist_ok=True)
    prefix=destination/'likelihood_with_threshold'
    assert not prefix.with_suffix('.pdf').exists()
    fig,ax=plt.subplots(figsize=(9,5))
    ax.plot(data[:,0],data[:,1],lw=.8,color='tab:blue',label='Saved final-state log likelihood')
    ax.axhline(row['threshold'],color='tab:red',ls='--',lw=1.2,label='Diagnostic threshold (tail prescription)')
    ax.axvline(row['catalogue_frequency_hz'],color='0.4',ls=':',lw=1,label='Catalogue frequency (reference)')
    ax.set(xlabel='Barycentric frequency (Hz)',ylabel='Log likelihood',title='47 Tuc %s: DM %.2f, $N_T=%d$'%(row['target'],row['dm'],row['segments']))
    ax.ticklabel_format(axis='x',style='plain',useOffset=False)
    ax.legend(fontsize=8);fig.text(.5,.015,'Diagnostic threshold; false-alarm probability is not calibrated.',ha='center',fontsize=8)
    fig.tight_layout(rect=[0,.035,1,1]);fig.savefig(prefix.with_suffix('.pdf'));fig.savefig(prefix.with_suffix('.png'),dpi=150);plt.close(fig)
    save(destination/'likelihood_threshold.json',dict(row,plotting_wall_s=time.monotonic()-start,curve_sha256=sha(row['curve'])))
    return time.monotonic()-start

def prepare():
    import numpy as np
    h=loadmodule(HELPER,'hmm_target_helper_prepare');agg=loadmodule(AGG,'original_threshold_prescription')
    old=json.loads((P/'followup/restored_all_415/native_presto_hmm_v1/manifest.json').read_text())
    catalogue={r['name'].replace('J0024-7204',''):r for r in json.loads((P/'followup/pipeline_assessment_v1/inputs/merged_identity_reference.json').read_text())}
    plots=[];groups=defaultdict(list);snapshots={};warnings=[]
    for directory in sorted(ROOT.glob('redetect_*')):
        if not directory.is_dir():continue
        target=directory.name.split('_',1)[1]
        for path in sorted(directory.glob('viterbi_out/**/*_loglikes.dat')):
            base=str(path)[:-len('_loglikes.dat')]
            match=re.search(r'DM([\d.]+)_Nt(\d+)',path.name)
            if not match:continue
            dm=float(match[1]);nt=int(match[2]);data=np.loadtxt(path,ndmin=2)
            assert data.shape[1]==2 and np.isfinite(data).all()
            ipk=int(np.argmax(data[:,1]));track=Path(base+'_track.dat');fit=Path(base+'_fit_summary.dat')
            assert track.is_file() and fit.is_file()
            t=np.loadtxt(track,ndmin=2);assert len(t)==nt
            spacing=float(np.median(np.diff(t[:,0]))) if nt>1 else float(t[0,0]*2)
            out=directory/'native_followup_v1'/('DM%.2f_Nt%d'%(dm,nt))
            record=dict(target=target,dm=dm,segments=nt,curve=str(path),track=str(track),summary=str(fit),
                peak_frequency_hz=float(data[ipk,0]),peak_loglike=float(data[ipk,1]),nq=len(data),
                segment_duration_s=spacing,track_end_s=float(t[-1,0]+spacing/2),
                catalogue_frequency_hz=float(catalogue[target]['F0_hz']),output=str(out),
                selection_scope='Strongest saved final-state likelihood, not a catalogue-forced frequency')
            plots.append(record);groups[(target,nt)].append((record,data[:,1]))
    selected=[];thresholds=[]
    for key,items in sorted(groups.items()):
        values=np.concatenate([d for r,d in items]);nq=int(np.median([len(d) for r,d in items]))
        # Small targeted grids cannot support the blind grid's 99.99th-percentile
        # tail. Use a declared 95th-percentile cutoff for every targeted group.
        tail,rate,threshold,error,n=agg.compute_paper_threshold(values,len(items),nq,.1,tail_percentile=95.)
        tr=dict(target=key[0],segments=key[1],threshold=threshold,l_tail=tail,lambda_hat=rate,
                formal_threshold_error=error,n_tail=n,n_curves=len(items),nq=nq,tail_percentile=95.,nominal_alpha=.1,
                interpretation='Historical exponential-tail prescription; diagnostic, not a calibrated FAP')
        thresholds.append(tr)
        for r,_ in items:r['threshold']=threshold;r['threshold_metadata']=tr
        best=max((r for r,d in items),key=lambda r:(r['peak_loglike'],-abs(r['dm']-catalogue[key[0]]['DM_pccm3'])))
        selected.append(best)
    rows=[];originals={}
    for candidate_number,item in enumerate(selected):
        directory=Path(item['output']);directory.mkdir(parents=True,exist_ok=True)
        values=fields(item['summary']);dm=item['dm'];tag='DM%.2f'%dm
        inf=ROOT/'stage2_dedisp'/(tag+'_bary')/('47Tuc_'+tag+'_bary.inf')
        info=h.inf_fields(inf);epoch=h.inf_get(info,'Epoch of observation')
        assert h.inf_get(info,'Barycentered?')=='1'
        cid='target_%s_DM%.2f_Nt%d'%(item['target'],dm,item['segments'])
        originals[cid]={str(p):dict(sha256=sha(p),size=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in map(Path,[item['summary'],item['track'],item['curve']])}
        if item['track_end_s']>old['filterbank_metadata']['duration_s']+2:
            warnings.append(dict(candidate_id=cid,warning='Historical segment grid extends beyond valid input duration; preserve saved epochs/coefficients, fold only available filterbank data.',track_end_s=item['track_end_s'],filterbank_duration_s=old['filterbank_metadata']['duration_s']))
        for model,prefix in [('polynomial','poly'),('kepler_taylor','kepler')]:
            order=int(values.get('poly_order_used',['5'])[0]) if model=='polynomial' else 5
            coefficients=[('0' if model=='polynomial' and k>order else values.get(prefix+'_f'+str(k),['nan'])[0]) for k in range(6)]
            finite=all(math.isfinite(float(c)) for c in coefficients)
            tref=values[prefix+'_tref_s'][0]
            row=dict(index=len(rows),candidate_number=candidate_number,candidate_id=cid,
                candidate_json='',summary=item['summary'],summary_sha256=sha(item['summary']),inf=str(inf),inf_sha256=sha(inf),
                model=model,coefficient_strings=coefficients,reference_epoch_mjd=str(Decimal(epoch)+Decimal(tref)/Decimal(86400)),
                input_epoch_mjd=epoch,reference_time_s=tref,dm='%.2f'%dm,
                ra=h.inf_get(info,'J2000 Right Ascension'),dec=h.inf_get(info,'J2000 Declination'),
                source_valid_samples=old['filterbank_metadata']['samples'],output_directory=str(directory/model),
                candidate_label='TGT%06d'%candidate_number,classification='targeted_inspection',rfi_flag='',target=item['target'],
                segments=item['segments'],peak_frequency_hz=item['peak_frequency_hz'],peak_loglike=item['peak_loglike'],
                coefficient_order_used=int(values.get('poly_order_used',['5'])[0]) if model=='polynomial' else 5,
                zero_padded_unfitted_orders=[k for k in range(6) if prefix+'_f'+str(k) not in values])
            row['zero_padded_unfitted_orders']=[k for k in range(6) if model=='polynomial' and k>order]
            duration=old['filterbank_metadata']['duration_s'];frequencies=[]
            if finite:
                for j in range(1001):
                    seconds=j*duration/1000-float(tref)
                    frequencies.append(sum(float(c)*seconds**k/math.factorial(k) for k,c in enumerate(coefficients)))
            row['frequency_range_hz']=[min(frequencies),max(frequencies)] if frequencies else None
            row['preflight_positive_frequency']=bool(frequencies and min(frequencies)>0)
            f0=float(coefficients[0]) if finite else item['peak_frequency_hz']
            row['estimated_bins']=min(64,max(1,math.ceil(1/(max(1e-9,f0)*float(old['filterbank_metadata']['header']['tsamp'])))))
            row['validation_reason']=('nonfinite_saved_coefficients' if not finite else ('nonpositive_saved_model_frequency' if min(frequencies)<=0 else 'preflight_passed'))
            Path(row['output_directory']).mkdir(exist_ok=True)
            if not row['preflight_positive_frequency']:
                save(Path(row['output_directory'])/'invalid_saved_model_v2.json',dict(row,reason='Saved model unavailable or non-positive; no Kepler coefficients invented, truncation or refitting.'))
            rows.append(row)
    csvsave(W/'thresholds.csv',thresholds)
    # Omit nested threshold metadata from the flat selected-run index.
    csvsave(W/'selected_runs.csv',[{k:v for k,v in r.items() if k!='threshold_metadata'} for r in selected])
    manifest=dict(rows=rows,plots=plots,selected_runs=selected,filterbank=old['filterbank'],filterbank_metadata=old['filterbank_metadata'],
        helper_sha256=sha(HELPER),threshold_code_sha256=sha(AGG),targets=sorted({r['target'] for r in plots}),
        selection='Strongest saved DM independently at each target and segment count; no threshold-based exclusions',
        input_models_refitted=False,all_saved_frequency_coefficients_retained=True,settings=old['settings'],historical_time_grid_warnings=warnings)
    save(W/'manifest.json',manifest);save(W/'original_fold_files.json',originals)
    # Render one curve to measure this exact plotting task before scheduling all.
    pilot=next(r for r in plots if r['target']=='C' and r['segments']==4)
    elapsed=make_plot(pilot);save(W/'plot_pilot.json',dict(output=pilot['output'],wall_s=elapsed,curve=pilot['curve']))
    (W/'plots_complete_pilot.txt').write_text(pilot['curve']+'\n')
    print(json.dumps(dict(targets=manifest['targets'],saved_curves=len(plots),selected_runs=len(selected),fold_models=len(rows),
        valid_models=sum(r['preflight_positive_frequency'] for r in rows),plot_pilot_s=elapsed,historical_time_grid_warnings=len(warnings)),indent=2),flush=True)

def plot(target):
    manifest=json.loads((W/'manifest.json').read_text());n=0
    for row in manifest['plots']:
        if row['target']!=target:continue
        if (Path(row['output'])/'likelihood_threshold.json').exists():continue
        make_plot(row);n+=1
    save(W/('plots_'+target+'_complete.json'),dict(target=target,new_plots=n))
    print('Rendered target',target,'plots',n,flush=True)

def native(index,validate=False):
    import numpy as np
    from presto import psr_utils
    def odd_safe_fft_rotate(arr,bins):
        arr=np.asarray(arr);f=np.arange(arr.size//2+1,dtype=float)
        return np.fft.irfft(np.exp(2j*np.pi*f*bins/arr.size)*np.fft.rfft(arr),n=arr.size)
    psr_utils.fft_rotate=odd_safe_fft_rotate
    manifest=json.loads((W/'manifest.json').read_text());assert sha(HELPER)==manifest['helper_sha256']
    h=loadmodule(HELPER,'target_native_helper');h.WORKFLOW=W
    row=manifest['rows'][index]
    if not row['preflight_positive_frequency']:print('Invalid saved model recorded',index);return
    h.run(index,validate_only=validate)

def report():
    manifest=json.loads((W/'manifest.json').read_text());rows=[]
    for r in manifest['rows']:
        d=Path(r['output_directory']);prefix='presto_'+r['model']+'_f0_f5'
        marker=d/(prefix+'_complete.json');record=json.loads(marker.read_text()) if marker.exists() else None
        failure=W/('failure_%04d_v7.json'%r['index'])
        status='complete' if record else ('invalid_saved_model' if not r['preflight_positive_frequency'] else ('validation_or_fold_failure' if failure.exists() else 'pending'))
        products=record['product_sha256'] if record else {}
        png=next((str(d/n) for n in products if n.endswith('.pfd.png')),'')
        pdf=next((str(d/n) for n in products if n.endswith('.pfd.pdf')),'')
        rows.append(dict(target=r['target'],segments=r['segments'],dm=r['dm'],model=r['model'],status=status,plot=png,pdf=pdf,
            likelihood_plot=str(d.parent/'likelihood_with_threshold.png'),directory=str(d)))
    version=int(time.time());csvsave(W/('fold_status_%d.csv'%version),rows)
    result=dict(status=dict(Counter(r['status'] for r in rows)),plots_complete=sum((Path(r['output'])/'likelihood_threshold.json').exists() for r in manifest['plots']),
        plots_expected=len(manifest['plots']),selected_runs=len(manifest['selected_runs']),targets=manifest['targets'],historical_time_grid_warnings=manifest['historical_time_grid_warnings'])
    save(W/('status_%d.json'%version),result)
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','plot','native','validate','report']);p.add_argument('--index',type=int);p.add_argument('--target');a=p.parse_args()
    if a.mode=='prepare':prepare()
    elif a.mode=='plot':plot(a.target)
    elif a.mode in ['native','validate']:native(a.index,a.mode=='validate')
    else:report()
