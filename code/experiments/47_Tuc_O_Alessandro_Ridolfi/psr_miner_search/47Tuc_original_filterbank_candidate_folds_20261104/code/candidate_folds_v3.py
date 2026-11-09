"""Re-fold saved HMM and PRESTO candidates on the original filterbank.

No search rerun, fit replacement, catalogue timing substitution, or old-file
overwrite. HMM polynomials retain F0--F5, including Kepler-derived Taylor fits.
PRESTO acceleration candidates retain the exact saved section and phase model.
"""
from pathlib import Path
from decimal import Decimal
from collections import Counter
from datetime import datetime,timezone
import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
import shutil
import stat
import subprocess
import tempfile
import time
import traceback
import native_helper as h

W=Path(__file__).resolve().parents[1]
E=W.parents[1]
P=E/'psr_miner_search/hmm_presto_v2_20261027T060742Z/production'
FIL=E/'data/47Tuc_22UL_1of2_L.fil'
sha=h.sha
def save(path,obj):
    path=Path(path);temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(obj,indent=2)+'\n');temp.replace(path)
def extrema(coeff,origin,start,stop):
    import numpy as np
    p=np.polynomial.Polynomial([float(v)/math.factorial(k) for k,v in enumerate(coeff)])
    roots=p.deriv().roots();low=start-origin;high=stop-origin
    times=[low,high]+[z.real for z in roots if abs(z.imag)<1e-7 and low<=z.real<=high]
    values=[float(p(x)) for x in times]
    return min(values),max(values)

def prepare():
    from presto import sigproc
    assert not (W/'manifest.json').exists()
    for name in ['rows','row_batches','records','logs','sources']:(W/name).mkdir(exist_ok=True)
    digests={}
    def cached_sha(path):
        key=str(path)
        if key not in digests:digests[key]=sha(path)
        return digests[key]
    header,hbytes=sigproc.read_header(str(FIL));assert header['nbits']==8 and header['nchans']==128 and header['nifs']==1
    samples=(FIL.stat().st_size-hbytes)//128;dt=float(header['tsamp']);duration=samples*dt
    sources=[('hmm_blind',P/'followup/restored_all_415/native_presto_hmm_v1/manifest.json'),
             ('presto',P/'followup/native_presto_results_v1/manifest.json'),
             ('hmm_targeted',P/'followup/targeted_native_v1/manifest.json'),
             ('hmm_supplementary_Nt16',P/'followup/supplementary_Nt16_v1/manifest.json')]
    rows=[];fullrows=[];source_records=[];snapshots={}
    for pipeline,path in sources:
        original=json.loads(path.read_text());archive=W/'sources'/(pipeline+'_manifest.json')
        if archive.exists():assert cached_sha(archive)==cached_sha(path)
        else:shutil.copyfile(path,archive)
        oldmeta=original.get('filterbank_metadata',{})
        oldheader=oldmeta.get('header',original.get('filterbank_header'))
        for key in ['tstart','tsamp','fch1','foff','nchans','nbits','nifs']:
            assert oldheader[key]==header[key],(pipeline,key,oldheader[key],header[key])
        source_records.append(dict(pipeline=pipeline,path=str(path),sha256=cached_sha(path),rows=len(original['rows'])))
        for r in original['rows']:
            r=dict(r);r['source_index']=r['index'];r['index']=len(rows);r['pipeline']=pipeline
            r['source_manifest']=str(path);r['original_input']=original['filterbank']
            output=Path(r['output_directory']);output.mkdir(parents=True,exist_ok=True)
            prefix='presto_original_filterbank_refined' if pipeline=='presto' else 'presto_original_filterbank_'+r['model']+'_f0_f5'
            assert not any(output.glob(prefix+'*')),output
            r['prefix']=prefix
            if str(output) not in snapshots:
                before={}
                with os.scandir(output) as entries:
                    for entry in entries:
                        s=entry.stat(follow_symlinks=False)
                        if stat.S_ISREG(s.st_mode):before[entry.name]=dict(size=s.st_size,mtime_ns=s.st_mtime_ns)
                snapshots[str(output)]=before
            if pipeline=='presto':
                # Absolute seconds from the common original epoch describe the
                # same searched interval. Fractions must change when the input
                # file duration changes. Never reuse cleaned-file fractions.
                start=float(r['requested_topocentric_start_s']);stop=float(r['requested_topocentric_stop_s'])
                assert 0<=start<stop<=duration+dt
                r['start_fraction']=start/duration;r['end_fraction']=min(1.,stop/duration)
                r['interval_start_s']=start;r['interval_stop_s']=stop
                assert cached_sha(r['candidate']['source_binary'])==r['candidate']['source_binary_sha256']
            else:
                r['start_fraction']=0.;r['end_fraction']=1.
                r['interval_start_s']=0.;r['interval_stop_s']=duration
                assert cached_sha(r['summary'])==r['summary_sha256'] and cached_sha(r['inf'])==r['inf_sha256']
            assert len(r['coefficient_strings'])==6
            finite=all(math.isfinite(float(c)) for c in r['coefficient_strings'])
            ref=float(r['reference_time_s'])
            # PRESTO section polynomials use an epoch at their own section
            # start. HMM epochs are saved seconds from the observation epoch.
            if not finite:freqrange=[None,None]
            elif pipeline=='presto':freqrange=extrema(r['coefficient_strings'],0,0,r['section_duration_s'])
            else:freqrange=extrema(r['coefficient_strings'],ref,0,duration)
            r['frequency_range_hz']=freqrange;r['preflight_positive_frequency']=finite and freqrange[0]>0
            r['estimated_bins']=min(64,max(1,math.ceil(1/(freqrange[0]*(1-2e-4)*dt)))) if r['preflight_positive_frequency'] else 64
            fullrows.append(r)
            rows.append(dict(index=r['index'],pipeline=pipeline,candidate_id=r['candidate_id'],model=r.get('model','acceleration'),
                output_directory=r['output_directory'],prefix=prefix,estimated_bins=r['estimated_bins'],
                duration_s=r['interval_stop_s']-r['interval_start_s'],valid=r['preflight_positive_frequency']))
    assert len({(r['output_directory'],r['prefix']) for r in rows})==len(rows)
    for offset in range(0,len(fullrows),500):save(W/'row_batches'/('%03d.json'%(offset//500)),fullrows[offset:offset+500])
    save(W/'old_files_snapshot.json',snapshots)
    save(W/'manifest.json',dict(rows=rows,sources=source_records,filterbank=str(FIL),filterbank_header=header,
        filterbank_size=FIL.stat().st_size,filterbank_mtime_ns=FIL.stat().st_mtime_ns,filterbank_samples=samples,
        filterbank_duration_s=duration,helper_sha256=sha(W/'code/native_helper.py'),old_snapshot_sha256=sha(W/'old_files_snapshot.json'),
        settings=dict(npart=64,nsub=32,ncpus=1,npfact=2,ndmfact=3,searchfdd=True),
        interval_policy='Full original observation for HMM; original searched intervals for PRESTO, with fractions recomputed for original input length.',
        search_rerun=False,coefficients_refitted=False,catalogue_timing_substituted=False))
    (W/'README.md').write_text('''# Original-filterbank HMM and PRESTO candidate folds

All native folds use the original `data/47Tuc_22UL_1of2_L.fil`.
Saved HMM F0--F5 and exact PRESTO acceleration-candidate coefficients are
retained. No search, candidate selection or track fitting is rerun.

HMM blind, targeted and supplementary Nt16 candidate sets are all included.
HMM models fold the full original observation. PRESTO models fold the same
searched intervals as previously; start/end fractions are recomputed against
the original input duration. Search inputs came from the cleaned search data,
which remains explicitly identified in every row's provenance.

All old plots and records remain unchanged. New products are added to each
existing model directory with `presto_original_filterbank_` filenames.
The full native PRESTO local search uses `-npfact 2 -ndmfact 3 -searchfdd`,
64 subintegrations and 32 subbands. Search parameters do not come from a
catalogue ephemeris. Invalid saved Taylor models remain recorded as invalid;
no coefficients are silently truncated or replaced.

`manifest.json`, immutable source snapshots, per-row records, old-file size/time
snapshots, preserved code, logs, final native PFDs and a consolidated result
index make the rerun retrievable. Profile significance is not a calibrated
blind-search significance. All candidates still require diagnostic review.
''')
    print(json.dumps(dict(rows=len(rows),counts=dict(Counter(r['pipeline'] for r in rows)),invalid=sum(not r['valid'] for r in rows)),indent=2),flush=True)

def run(index,pilot=False):
    import numpy as np
    from presto import psr_utils
    def rotate(arr,bins):
        arr=np.asarray(arr);freq=np.arange(arr.size//2+1,dtype=float)
        return np.fft.irfft(np.exp(2j*np.pi*freq*bins/arr.size)*np.fft.rfft(arr),n=arr.size)
    psr_utils.fft_rotate=rotate
    from presto.prepfold import pfd
    m=json.loads((W/'manifest.json').read_text());r=json.loads((W/'row_batches'/('%03d.json'%(index//500))).read_text())[index%500]
    record_path=W/'records'/('%05d%s.json'%(index,'_pilot' if pilot else ''))
    if record_path.exists():
        prior=json.loads(record_path.read_text());assert prior['status'] in ('complete','invalid_saved_model'),prior
        return
    target=(W/'pilots'/('%05d'%index)) if pilot else Path(r['output_directory'])
    target.mkdir(parents=True,exist_ok=True);prefix=r['prefix']
    assert not any(target.glob(prefix+'*')),'Preserve existing attempt: '+str(target)
    assert FIL.stat().st_size==m['filterbank_size'] and FIL.stat().st_mtime_ns==m['filterbank_mtime_ns']
    assert sha(W/'code/native_helper.py')==m['helper_sha256']
    if r['pipeline']=='presto':assert sha(r['candidate']['source_binary'])==r['candidate']['source_binary_sha256']
    else:assert sha(r['summary'])==r['summary_sha256'] and sha(r['inf'])==r['inf_sha256']
    rec=dict(index=index,row=r,pilot=pilot,filterbank=str(FIL),job_id=os.environ.get('SLURM_JOB_ID'),
        started_utc=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),status='running')
    if not r['preflight_positive_frequency']:
        rec.update(status='invalid_saved_model',reason='Saved six-coefficient model has non-finite coefficients or non-positive frequency within the observation; no truncation or refitting.')
        save(record_path,rec);save(target/(prefix+'_invalid_model.json'),rec);return
    started=time.monotonic();scratch=Path(tempfile.mkdtemp(prefix='original_%05d_'%index,dir=os.environ.get('SLURM_TMPDIR','/tmp')))
    save(record_path,rec)
    try:
        duration=r['section_duration_s'] if r['pipeline']=='presto' else m['filterbank_duration_s']
        bary_start=float(r['input_epoch_mjd'])
        bary=h.tempo_polycos(scratch/'bary_check',r,'@',0.,0.,bary_start,bary_start+duration/86400)
        rec['transfer_check']=h.verify_bary_transfer(bary,r,duration)
        header=m['filterbank_header'];freq=float(header['fch1'])+127*float(header['foff'])/2
        start=float(header['tstart'])+r['interval_start_s']/86400
        stop=float(header['tstart'])+r['interval_stop_s']/86400
        topo=h.tempo_polycos(scratch/'topocentric',r,'m',freq,r['dm'],start-2/86400,stop+2/86400)
        values=[h.evaluate(topo,Decimal(str(start))+Decimal(str(r['interval_stop_s']-r['interval_start_s']))*k/Decimal(86400000))[1] for k in range(1001)]
        assert all(math.isfinite(v) and v>0 for v in values)
        rec['topocentric_frequency_range_hz']=[min(values),max(values)]
        text=(scratch/'topocentric/polyco_for_observation.dat').read_text()
        text=re.sub(r'(?m)^(\s*)0000\+0000(?=\s)',lambda match:match.group(1)+r['candidate_label'],text)
        (scratch/'model.polycos').write_text(text);(scratch/'input.fil').symlink_to(FIL)
        command=['prepfold','-filterbank','-topo','-noxwin','-ncpus','1','-dm',r['dm'],'-nsub','32','-npart','64',
            '-psr',r['candidate_label'],'-polycos','model.polycos','-start',format(r['start_fraction'],'.17g'),
            '-end',format(r['end_fraction'],'.17g'),'-o',prefix,'input.fil']
        command+=['-nosearch'] if pilot else ['-npfact','2','-ndmfact','3','-searchfdd']
        rec['command']=command;save(record_path,rec)
        print('Starting',index,r['pipeline'],r['candidate_id'],r.get('model','acceleration'),'pilot',pilot,flush=True)
        native_start=time.monotonic();optimizer_start=None;logpath=scratch/'prepfold.log'
        with logpath.open('w') as log:
            process=subprocess.Popen(command,cwd=scratch,stdout=log,stderr=subprocess.STDOUT)
            while process.poll() is None:
                if optimizer_start is None and not pilot and 'Optimizing...' in logpath.read_text(errors='replace'):
                    optimizer_start=time.monotonic();rec['optimizer_start_after_s']=optimizer_start-native_start;save(record_path,rec)
                time.sleep(1)
        rec['returncode']=process.returncode;assert process.returncode==0
        native_stop=time.monotonic();rec['raw_fold_s']=(optimizer_start or native_stop)-native_start
        rec['optimization_and_render_s']=native_stop-(optimizer_start or native_stop)
        paths=list(scratch.glob(prefix+'*.pfd'));assert len(paths)==1
        obj=pfd(str(paths[0]));assert obj.npart==64 and obj.nsub==32 and np.isfinite(obj.profs).all()
        rec['profile_shape']=list(obj.profs.shape);rec['folded_samples']=int(obj.stats[:,:,0].sum()/obj.nsub)
        rec['folded_topocentric_epoch_mjd']=float(obj.tepoch)
        grid=re.search(r'Searching (\d+) DMs, (\d+) periods, (\d+) p-dots, and (\d+) p-dotdots',logpath.read_text(errors='replace'))
        if not pilot:
            assert grid and optimizer_start is not None
            rec['trial_dimensions']=[int(x) for x in grid.groups()];rec['trial_count']=math.prod(rec['trial_dimensions'])
        rec['bestprof_header']=[l for l in Path(str(paths[0])+'.bestprof').read_text().splitlines() if l.startswith('#')]
        rec['sigma']=None
        for line in rec['bestprof_header']:
            match=re.search(r'\(~?\s*([\d.eE+-]+)\s+sigma\)',line)
            if match:rec['sigma']=float(match[1])
        h.render_ps(Path(str(paths[0])+'.ps'))
        products={prefix+'.pfd':paths[0]}
        products.update({prefix+'.pfd'+suffix:Path(str(paths[0])+suffix) for suffix in ['.bestprof','.ps','.pdf','.png']})
        products.update({prefix+'.par':scratch/'topocentric/model.par',prefix+'.polycos':scratch/'model.polycos',
            prefix+'_prepfold.log':logpath,prefix+'_tempo.log':scratch/'topocentric/tempo.log',
            prefix+'_bary_check.polycos':scratch/'bary_check/polyco_for_observation.dat',prefix+'_tempo_raw.polycos':scratch/'topocentric/polyco.dat'})
        assert not any((target/name).exists() for name in products)
        for name,path in products.items():
            with (target/name).open('xb') as dst,path.open('rb') as src:shutil.copyfileobj(src,dst)
        rec['products']={name:sha(target/name) for name in products}
        rec['pdf']=str(target/(prefix+'.pfd.pdf'));rec['png']=str(target/(prefix+'.pfd.png'))
        # Verify the pre-existing files in this candidate folder were untouched.
        snapshot=json.loads((W/'old_files_snapshot.json').read_text())[r['output_directory']]
        for name,before in snapshot.items():
            stat=(Path(r['output_directory'])/name).stat()
            assert stat.st_size==before['size'] and stat.st_mtime_ns==before['mtime_ns'],name
        assert FIL.stat().st_size==m['filterbank_size'] and FIL.stat().st_mtime_ns==m['filterbank_mtime_ns']
        rec.update(status='complete',wall_s=time.monotonic()-started,completed_utc=datetime.now(timezone.utc).isoformat())
        save(target/(prefix+'_complete.json'),rec);save(record_path,rec)
        shutil.rmtree(scratch);print('Complete',index,'sigma',rec['sigma'],'seconds',rec['wall_s'],flush=True)
    except Exception:
        rec.update(status='failed',wall_s=time.monotonic()-started,error=traceback.format_exc())
        failure=W/'failures'/('%05d%s'%(index,'_pilot' if pilot else ''));failure.parent.mkdir(exist_ok=True)
        shutil.copytree(scratch,failure);rec['failure_details']=str(failure)
        save(record_path,rec);raise

def report():
    m=json.loads((W/'manifest.json').read_text());rows=[]
    for r in m['rows']:
        p=W/'records'/('%05d.json'%r['index']);rec=json.loads(p.read_text()) if p.exists() else {}
        rows.append(dict(index=r['index'],pipeline=r['pipeline'],candidate_id=r['candidate_id'],model=r['model'],
            status=rec.get('status','pending'),nominal_sigma=rec.get('sigma'),wall_s=rec.get('wall_s'),
            pdf=rec.get('pdf'),png=rec.get('png'),record=str(p),error=rec.get('error','')[-1000:]))
    save(W/'results.json',dict(updated_utc=datetime.now(timezone.utc).isoformat(),filterbank=str(FIL),
        expected=len(rows),status=dict(Counter(r['status'] for r in rows)),rows=rows))
    stream=io.StringIO();writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (W/'original_filterbank_fold_index.csv').write_text(stream.getvalue())
    print(json.dumps(dict(expected=len(rows),status=dict(Counter(r['status'] for r in rows))),indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','fold','report']);p.add_argument('--index',type=int);p.add_argument('--pilot',action='store_true');a=p.parse_args()
    prepare() if a.mode=='prepare' else report() if a.mode=='report' else run(a.index,a.pilot)
