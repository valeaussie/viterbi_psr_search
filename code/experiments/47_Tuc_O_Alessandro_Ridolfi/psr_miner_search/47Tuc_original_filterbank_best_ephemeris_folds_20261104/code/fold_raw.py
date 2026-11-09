"""Matched original-filterbank folds retaining every independent timing model.

Full PRESTO local refinement: period, first and second period derivatives, DM.
This measures the best local fold around a supplied ephemeris, not a newly
phase-connected timing solution. Timing models are never fitted to HMM tracks.
"""
from pathlib import Path
from decimal import Decimal
from datetime import datetime,timezone
import argparse
import json
import math
import os
import re
import shutil
import subprocess
import time
import traceback
import fold_engine as f

W=Path(__file__).resolve().parents[1]
E=W.parents[1]
FIL=E/'data/47Tuc_22UL_1of2_L.fil'

def prepare():
    from presto import sigproc
    assert not (W/'manifest.json').exists()
    source=json.loads((W/'source_manifest.json').read_text())
    header,hbytes=sigproc.read_header(str(FIL))
    assert header['nbits']==8 and header['nchans']==128 and header['nifs']==1
    assert (FIL.stat().st_size-hbytes)%128==0
    samples=(FIL.stat().st_size-hbytes)//128
    duration=samples*float(header['tsamp']);start=float(header['tstart']);stop=start+duration/86400
    rows=[]
    for source_row in source['models']:
        r=dict(source_row);fields=r['fields'];original=Path(r['path'])
        assert f.sha(original)==r['sha256']
        key=next(k for k in ('PSR','PSRJ','PSRB') if k in fields)
        name=fields[key][0].removeprefix('J')
        assert re.match(r'^\d{4}',name),name
        out=W/'folds'/('47Tuc_'+r['target'])/r['label'];out.mkdir(parents=True)
        adapted=re.sub(r'(?m)^(?:PSR|PSRJ|PSRB)\s+\S+',lambda m:'PSR '+name,original.read_text())
        model=out/'model_tempo.par';model.write_text(adapted)
        r.update(output=str(out),psr_name=name,par_original=str(original),par_sha256=r['sha256'],
            model_tempo=str(model),model_sha256=f.sha(model),binary=fields.get('BINARY',['isolated'])[0],
            dm=f.number(fields['DM'][0]),f0=f.number(fields['F0'][0]),
            fitted_start_mjd=f.number(fields['START'][0]) if 'START' in fields else None,
            fitted_finish_mjd=f.number(fields['FINISH'][0]) if 'FINISH' in fields else None,
            warning='Uses the published second-harmonic phase convention, not a verified fundamental-period model.' if r['target']=='aa' else '')
        # Upper bound on the native default bin count; the actual count is
        # read from the PFD and its exact trial dimensions are recorded.
        r['estimated_bins']=min(64,math.ceil(1/(r['f0']*float(header['tsamp']))))
        rows.append(r)
    f.save(W/'manifest.json',dict(rows=rows,filterbank=str(FIL),filterbank_size=FIL.stat().st_size,
        filterbank_mtime_ns=FIL.stat().st_mtime_ns,filterbank_header=header,filterbank_samples=samples,
        filterbank_header_bytes=hbytes,filterbank_duration_s=duration,observation_start_topocentric_mjd=start,
        observation_stop_topocentric_mjd=stop,source_manifest_sha256=f.sha(W/'source_manifest.json'),
        helper_sha256=f.sha(W/'code/verified_native_helper.py'),unavailable=source['unavailable'],
        full_refinement=dict(npfact=2,ndmfact=3,searchfdd=True,npart=64,nsub=32,ncpus=1,
        profile_bins='Native PRESTO default; actual value and search grid saved for every fold.'),
        comparison='Original uncleaned channelized filterbank. Full independently supplied ephemeris plus identical native PRESTO local refinement. Report all alternatives and the strongest fold for review.'))
    (W/'logs').mkdir(exist_ok=True)
    (W/'README.md').write_text('''# Original-filterbank ephemeris reference folds

All folds here use `data/47Tuc_22UL_1of2_L.fil`, the original channelized
filterbank. Earlier RFI-cleaned fixed-model results remain in their own folder.

`source_manifest.json` preserves latest and older Freire-site files, pre-existing
reference files, duplicates, SHA256 hashes and provenance. The paper O
observation-epoch BT model is an explicit alternative. No HMM-derived model is
used as an independent reference. No new catalogue-portal query was made.

Full folds use 64 subintegrations, 32 subbands, the native profile-bin default,
and `-npfact 2 -ndmfact 3 -searchfdd`. All physical timing-model coefficients
enter through phase-validated TEMPO predictors. Only the internal PSR label is
adapted for TEMPO compatibility. PRESTO refines residual spin parameters and DM;
it does not derive a new phase-connected orbital ephemeris.

`pilots/` contains fixed-model original-data I/O measurements, not final refined
folds. `folds/47Tuc_<name>/<model>/` contains the final native PFD, best profile,
PS, PDF and PNG, model, predictor, logs and a complete record. Existing products
are never overwritten. Every failure remains saved in its original directory.

The report ranks the nominal profile significances among completed models. A
maximum across models and refinement trials is not a calibrated detection
significance. Inspect phase stability, pulse profile, frequency consistency and
DM before accepting a best reference or claiming recovery. Missing public
ephemerides remain explicit, rather than being counted as negative detections.
''')
    print(json.dumps(dict(models=len(rows),pulsars=len({r['target'] for r in rows}),filterbank=str(FIL),duration_s=duration),indent=2))

def run(index,pilot=False):
    import numpy as np
    from presto import psr_utils
    def rotate(arr,bins):
        arr=np.asarray(arr);freq=np.arange(arr.size//2+1,dtype=float)
        return np.fft.irfft(np.exp(2j*np.pi*freq*bins/arr.size)*np.fft.rfft(arr),n=arr.size)
    psr_utils.fft_rotate=rotate
    from presto.prepfold import pfd
    m=json.loads((W/'manifest.json').read_text());r=m['rows'][index]
    out=W/'pilots'/r['model_id'] if pilot else Path(r['output'])
    out.mkdir(parents=True,exist_ok=True)
    marker=out/'complete.json'
    if marker.exists():
        prior=json.loads(marker.read_text());assert all(f.sha(out/k)==v for k,v in prior['products'].items());return
    assert not (out/'running.json').exists(),'Preserve prior attempt and use a new attempt directory.'
    assert FIL.stat().st_size==m['filterbank_size'] and FIL.stat().st_mtime_ns==m['filterbank_mtime_ns']
    assert f.sha(r['par_original'])==r['par_sha256'] and f.sha(r['model_tempo'])==r['model_sha256']
    rec=dict(row=r,pilot=pilot,filterbank=str(FIL),slurm_job_id=os.environ.get('SLURM_JOB_ID'),
        refinement=not pilot,started_utc=datetime.now(timezone.utc).isoformat(),worker_sha256=f.sha(__file__),
        engine_sha256=f.sha(W/'code/fold_engine.py'),helper_sha256=f.sha(W/'code/verified_native_helper.py'))
    started=time.monotonic()
    try:
        h=m['filterbank_header'];start=m['observation_start_topocentric_mjd'];stop=m['observation_stop_topocentric_mjd']
        freq=float(h['fch1'])+(int(h['nchans'])-1)*float(h['foff'])/2
        a=f.generate(r,out/'tempo_30min_12coeff',30,start,stop,freq,12)
        b=f.generate(r,out/'tempo_30min_15coeff_validation',30,start,stop,freq,15)
        phase=f.phase_extrema_check(a,b,start,stop)
        assert phase['maximum_phase_difference_cycles']<1e-4,phase
        frequencies=[]
        for k in range(1001):
            t=Decimal(str(start))+Decimal(str(m['filterbank_duration_s']))*k/Decimal(86400000)
            frequencies.append(f.evaluate_decimal(a,t)[1])
        assert all(math.isfinite(v) and v>0 for v in frequencies)
        rec['phase_check']=phase;rec['topocentric_frequency_range_hz']=[min(frequencies),max(frequencies)]
        rec['default_bin_upper_bound']=min(64,math.ceil(1/(min(frequencies)*float(h['tsamp']))))
        f.save(out/'phase_check.json',dict(phase_check=phase,frequency_range_hz=rec['topocentric_frequency_range_hz']))
        (out/'input.fil').symlink_to(FIL)
        command=['prepfold','-filterbank','-topo','-noxwin','-ncpus','1','-npart','64','-nsub','32',
            '-psr',r['psr_name'][:10],'-polycos','tempo_30min_12coeff/observation.polycos','-o','reference','input.fil']
        if pilot:command+=['-nosearch']
        else:command+=['-npfact','2','-ndmfact','3','-searchfdd']
        rec['command']=command;f.save(out/'running.json',rec)
        print('Starting',index,r['model_id'],'pilot',pilot,flush=True)
        native_start=time.monotonic();optimization_start=None;logfile=out/'prepfold.log'
        with logfile.open('w') as log:
            process=subprocess.Popen(command,cwd=out,stdout=log,stderr=subprocess.STDOUT)
            while process.poll() is None:
                if optimization_start is None and not pilot and 'Optimizing...' in logfile.read_text(errors='replace'):
                    optimization_start=time.monotonic()
                time.sleep(1)
        rec['returncode']=process.returncode;assert process.returncode==0
        rec['native_wall_s']=time.monotonic()-native_start
        rec['raw_fold_s']=(optimization_start or time.monotonic())-native_start
        rec['optimization_and_render_s']=time.monotonic()-(optimization_start or time.monotonic())
        paths=list(out.glob('reference*.pfd'));assert len(paths)==1
        folded=pfd(str(paths[0]));assert folded.npart==64 and folded.nsub==32
        assert np.isfinite(folded.profs).all()
        rec['profile_shape']=list(folded.profs.shape)
        rec['folded_samples']=int(folded.stats[:,:,0].sum()/folded.nsub)
        rec['folded_topocentric_epoch_mjd']=float(folded.tepoch)
        text=logfile.read_text(errors='replace')
        grid=re.search(r'Searching (\d+) DMs, (\d+) periods, (\d+) p-dots, and (\d+) p-dotdots',text)
        if not pilot:
            assert grid and optimization_start is not None,'Verify full native refinement grid'
            rec['trial_dimensions']=[int(v) for v in grid.groups()];rec['trial_count']=math.prod(rec['trial_dimensions'])
        headers=[l for l in Path(str(paths[0])+'.bestprof').read_text().splitlines() if l.startswith('#')]
        rec['bestprof_header']=headers;rec['sigma']=None;rec['reduced_chi2']=None
        for line in headers:
            match=re.search(r'\(~?\s*([\d.eE+-]+)\s+sigma\)',line)
            if match:rec['sigma']=float(match[1])
            if line.startswith('# Reduced chi-sqr'):rec['reduced_chi2']=float(line.split('=')[-1])
        pdf,png=f.helper().render_ps(Path(str(paths[0])+'.ps'))
        rec['pdf']=str(pdf);rec['png']=str(png)
        for label,path in [('fold.pdf',pdf),('fold.png',png),('fold.pfd',paths[0])]:
            (out/label).symlink_to(Path(path).name)
        assert f.sha(r['par_original'])==r['par_sha256'] and f.sha(r['model_tempo'])==r['model_sha256']
        assert FIL.stat().st_size==m['filterbank_size'] and FIL.stat().st_mtime_ns==m['filterbank_mtime_ns']
        rec['wall_s']=time.monotonic()-started;rec['completed_utc']=datetime.now(timezone.utc).isoformat()
        rec['products']={str(p.relative_to(out)):f.sha(p) for p in out.rglob('*') if p.is_file() and not p.is_symlink() and p.name not in ('running.json','complete.json')}
        f.save(marker,rec)
        print('Complete',r['model_id'],'sigma',rec['sigma'],'seconds',rec['wall_s'],flush=True)
    except Exception:
        rec['wall_s']=time.monotonic()-started;rec['traceback']=traceback.format_exc()
        f.save(out/'failure.json',rec);raise

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','fold']);p.add_argument('--index',type=int);p.add_argument('--pilot',action='store_true');a=p.parse_args()
    prepare() if a.mode=='prepare' else run(a.index,a.pilot)
