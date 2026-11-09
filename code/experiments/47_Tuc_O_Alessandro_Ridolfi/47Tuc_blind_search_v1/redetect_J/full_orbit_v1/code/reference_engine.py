"""Fold the original published Freire-site timing models without parameter searches.

Only the internal PSR name is stripped of an initial J for TEMPO TZ compatibility.
Every numeric timing parameter, clock, ephemeris and binary model is retained.
Original downloads are immutable. A finer independent polyco grid checks the
numerical interpolation, not the physical validity of extrapolated ephemerides.
"""
from pathlib import Path
from decimal import Decimal, getcontext
from datetime import datetime, timezone
import argparse
import hashlib
import importlib.util
import json
import math
import os
import re
import shutil
import subprocess
import time
import traceback

getcontext().prec=40
W=Path(__file__).resolve().parents[1]
E=W.parents[1]
P=E/'psr_miner_search/hmm_presto_v2_20261027T060742Z/production'
FIL=E/'47Tuc_blind_search_v1/stage1_clean/47Tuc_22UL_1of2_L_RFIcleaned_01.fil'
HELPER=P/'followup/restored_all_415/code/native_hmm_prepfold_v7.py'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def save(path,data):
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(path)

def helper():
    path=W/'code/verified_native_helper.py'
    s=importlib.util.spec_from_file_location('published_verified_native_helper',path)
    h=importlib.util.module_from_spec(s);s.loader.exec_module(h);return h

def number(value):return float(value.replace('D','E'))

def prepare():
    from presto import sigproc
    assert not (W/'manifest.json').exists()
    archive=W/'code/verified_native_helper.py'
    assert not archive.exists();shutil.copyfile(HELPER,archive)
    header,header_bytes=sigproc.read_header(str(FIL))
    assert header['nifs']==1 and header['nbits']==8 and header['nchans']==128
    samples=(FIL.stat().st_size-header_bytes)//128
    duration=samples*float(header['tsamp'])
    start=float(header['tstart']);stop=start+duration/86400
    freire=json.loads((W/'source_manifest.json').read_text())
    rows=[]
    for i,r in enumerate(freire['available']):
        assert sha(r['path'])==r['sha256']
        fields=r['fields'];name=fields['PSR'][0].removeprefix('J')
        # TEMPO TZ mode requires a numeric HHMM at the start of its PSR label.
        assert re.match(r'^\d{4}',name),name
        text=Path(r['path']).read_text()
        adapted=re.sub(r'(?m)^(PSR\s+)J',r'\g<1>',text)
        out=W/'folds'/('47Tuc_'+r['target'])
        out.mkdir(parents=True,exist_ok=True)
        model=out/'model_tempo.par';assert not model.exists();model.write_text(adapted)
        span_start=number(fields['START'][0]) if 'START' in fields else None
        span_stop=number(fields['FINISH'][0]) if 'FINISH' in fields else None
        rows.append(dict(index=i,target=r['target'],psr_name=name,par_original=r['path'],par_sha256=r['sha256'],
            source_url=r['url'],model_tempo=str(model),model_sha256=sha(model),output=str(out),
            changes=['Removed leading J from internal PSR label for TEMPO TZ compatibility.'] if text!=adapted else [],
            binary=fields.get('BINARY',['isolated'])[0],dm=number(fields['DM'][0]),f0=number(fields['F0'][0]),
            ephem=fields.get('EPHEM',[''])[0],clock=fields.get('CLK',[''])[0],
            fitted_start_mjd=span_start,fitted_finish_mjd=span_stop,
            extrapolated=(span_start is not None and start<span_start) or (span_stop is not None and stop>span_stop),
            days_after_fitted_finish=max(0,start-span_stop) if span_stop else None,
            source_warning=('The linked aa file uses 541.89365 Hz although the website says the true spin period is twice the original value. This fold retains the linked file as supplied, hence its second-harmonic phase model.' if r['target']=='aa' else '')))
    save(W/'manifest.json',dict(rows=rows,source_manifest_sha256=sha(W/'source_manifest.json'),
        filterbank=str(FIL),filterbank_size=FIL.stat().st_size,filterbank_mtime_ns=FIL.stat().st_mtime_ns,
        filterbank_header=header,filterbank_duration_s=duration,filterbank_samples=samples,
        observation_start_topocentric_mjd=start,observation_stop_topocentric_mjd=stop,
        helper_sha256=sha(archive),folding=dict(refinement=False,npart=64,nsub=32,profile_bins='Native PRESTO default',
        data='Same four-hour RFI-cleaned channelized filterbank used for the HMM/PRESTO comparison.',
        native_sigma='Native PRESTO profile chi-square significance at fixed published timing model and DM; not a blind-search significance.'),
        unavailable=freire['unavailable']))
    print(json.dumps(dict(targets=len(rows),duration_s=duration,all_extrapolated=all(r['extrapolated'] for r in rows)),indent=2),flush=True)

def generate(row,path,span,start,stop,frequency,ncoeff=12):
    path.mkdir()
    shutil.copyfile(row['model_tempo'],path/'model.par')
    (path/'tz.in').write_text('m 12 %d %d 430\n\n\n%s %d %d 12 %.9f\n'%(span,ncoeff,row['psr_name'],span,ncoeff,frequency))
    with (path/'tempo.log').open('w') as log:
        result=subprocess.run(['tempo','-z','-f','model.par'],cwd=path,
            input='%d %d\n'%(math.floor(start)-1,math.ceil(stop)),text=True,stdout=log,stderr=subprocess.STDOUT,timeout=120)
    if result.returncode or not (path/'polyco.dat').exists():
        raise RuntimeError('TEMPO failed: '+str(path/'tempo.log')+'\n'+(path/'tempo.log').read_text()[-4000:])
    blocks=helper().read_polycos(path/'polyco.dat',start,stop,path/'observation.polycos')
    # TEMPO's polyco output field is 10 characters, although its par-file PSR
    # parser accepts the longer aa/ab names. No model parameter is changed.
    assert blocks and all(b['name']==row['psr_name'][:10] for b in blocks)
    return blocks

def evaluate_decimal(blocks,epoch):
    possible=[b for b in blocks if abs(epoch-b['midpoint'])*Decimal(1440)<Decimal(b['span_minutes'])/2]
    if not possible:raise ValueError('Polyco coverage gap at '+str(epoch))
    b=min(possible,key=lambda b:abs(epoch-b['midpoint']))
    t=(epoch-b['midpoint'])*Decimal(1440)
    co=[Decimal(str(c)) for c in b['coefficients']]
    phase=b['phase']+t*60*Decimal(str(b['frequency']))+sum(c*t**k for k,c in enumerate(co))
    freq=Decimal(str(b['frequency']))+sum(k*c*t**(k-1)/60 for k,c in enumerate(co) if k)
    return phase,float(freq)

def phase_extrema_check(a,b,start,stop):
    """Check the difference polynomial at every stationary point and endpoint.

    On equal-span TEMPO grids the midpoint epochs match. The difference of two
    phase interpolants is itself a polynomial, so its extrema bound the phase
    discrepancy over each complete interval, including between sample epochs.
    A derivative discrepancy alone need not accumulate: every polyco is anchored
    to the full physical timing model. Folding requires accurate phase.
    """
    import numpy as np
    other={block['midpoint']:block for block in b};maximum=0.;intervals=[]
    for first in a:
        second=other[first['midpoint']]
        assert first['span_minutes']==second['span_minutes']==30
        midpoint=first['midpoint'];half=Decimal(first['span_minutes'])/2
        low=max(-half,(Decimal(str(start))-midpoint)*1440)
        high=min(half,(Decimal(str(stop))-midpoint)*1440)
        if low>high:continue
        coefficients=[Decimal(str(x))-Decimal(str(y)) for x,y in
            zip(first['coefficients']+[0.]*(len(second['coefficients'])-len(first['coefficients'])),second['coefficients'])]
        offset=first['phase']-second['phase'];offset-=offset.to_integral_value()
        coefficients[0]+=offset
        coefficients[1]+=60*(Decimal(str(first['frequency']))-Decimal(str(second['frequency'])))
        polynomial=np.polynomial.Polynomial([float(c) for c in coefficients])
        roots=polynomial.deriv().roots()
        points=[float(low),float(high)]+[float(z.real) for z in roots if abs(z.imag)<1e-7 and float(low)<=z.real<=float(high)]
        # A dense grid adds a numerical cross-check on the polynomial root solver.
        points+=np.linspace(float(low),float(high),101).tolist()
        peak=max(abs(float(polynomial(t))) for t in points)
        assert math.isfinite(peak)
        maximum=max(maximum,peak)
        intervals.append(dict(midpoint_mjd=str(midpoint),maximum_phase_difference_cycles=peak,stationary_points_checked=len(roots)))
    assert intervals
    return dict(maximum_phase_difference_cycles=maximum,intervals=intervals,
        method='Difference polynomial extrema at all real derivative roots and interval endpoints, with an additional dense grid check.')

def run(target):
    import numpy as np
    from presto import psr_utils
    # Preserve odd profile length when reading PFDs with this installed Python API.
    def odd_safe_fft_rotate(arr,bins):
        arr=np.asarray(arr);f=np.arange(arr.size//2+1,dtype=float)
        return np.fft.irfft(np.exp(2j*np.pi*f*bins/arr.size)*np.fft.rfft(arr),n=arr.size)
    psr_utils.fft_rotate=odd_safe_fft_rotate
    from presto.prepfold import pfd
    m=json.loads((W/'manifest.json').read_text());row=next(r for r in m['rows'] if r['target']==target)
    original_out=Path(row['output'])
    out=original_out;attempt=1
    while (out/'failure.json').exists():
        attempt+=1;out=original_out/('attempt_%02d'%attempt)
    out.mkdir(exist_ok=True)
    marker=out/'complete.json'
    if marker.exists():
        rec=json.loads(marker.read_text());assert all(sha(out/n)==h for n,h in rec['products'].items());print('Already complete',target);return
    assert not (out/'failure.json').exists(),'Failed prior attempt must remain preserved; use a new attempt directory.'
    assert FIL.stat().st_size==m['filterbank_size'] and FIL.stat().st_mtime_ns==m['filterbank_mtime_ns']
    assert sha(row['par_original'])==row['par_sha256'] and sha(row['model_tempo'])==row['model_sha256']
    record=dict(target=target,row=row,command=None,started_utc=datetime.now(timezone.utc).isoformat(),
        slurm_job_id=os.environ.get('SLURM_JOB_ID'),worker_sha256=sha(__file__),helper_sha256=sha(W/'code/verified_native_helper.py'),
        prepfold_executable=shutil.which('prepfold'),tempo_executable=shutil.which('tempo'),refinement=False)
    started=time.perf_counter()
    try:
        header=m['filterbank_header'];start=m['observation_start_topocentric_mjd'];stop=m['observation_stop_topocentric_mjd']
        freqs=[float(header['fch1'])+k*float(header['foff']) for k in range(int(header['nchans']))]
        fobs=(min(freqs)+max(freqs))/2
        a=generate(row,out/'tempo_30min_12coeff',30,start,stop,fobs,12)
        b=generate(row,out/'tempo_30min_15coeff_validation',30,start,stop,fobs,15)
        phase_error=[];freq_error=[];predicted=[]
        for k in range(1001):
            t=Decimal(str(start))+Decimal(str(m['filterbank_duration_s']))*k/Decimal(1000*86400)
            pa,fa=evaluate_decimal(a,t);pb,fb=evaluate_decimal(b,t)
            diff=pa-pb;diff-=diff.to_integral_value()
            phase_error.append(abs(float(diff)));freq_error.append(abs(fa-fb));predicted.append(fa)
        interval_check=phase_extrema_check(a,b,start,stop)
        validation=dict(samples=1001,max_phase_difference_cycles=max(phase_error),
            max_frequency_difference_hz=max(freq_error),phase_tolerance_cycles=1e-4,
            full_interval_phase_check=interval_check,
            frequency_difference_is_diagnostic=True,
            frequency_range_hz=[min(predicted),max(predicted)],method='Compare independently generated 12-coefficient and 15-coefficient TEMPO interpolants on 30-minute spans; retains every original timing parameter.',
            physical_ephemeris_validity_verified=False)
        assert validation['max_phase_difference_cycles']<1e-4,validation
        assert interval_check['maximum_phase_difference_cycles']<1e-4,validation
        assert all(math.isfinite(f) and f>0 for f in predicted)
        record['polyco_validation']=validation;save(out/'polyco_validation.json',validation)
        source=out/'input.fil';source.symlink_to(FIL)
        command=['prepfold','-filterbank','-topo','-noxwin','-ncpus','1','-npart','64','-nsub','32',
            '-nosearch','-psr',row['psr_name'][:10],'-polycos','tempo_30min_12coeff/observation.polycos',
            '-o','published_ephemeris','input.fil']
        record['command']=command;save(out/'running.json',record)
        print('Starting',target,json.dumps(validation),flush=True)
        with (out/'prepfold.log').open('w') as log:
            result=subprocess.run(command,cwd=out,stdout=log,stderr=subprocess.STDOUT)
        record['returncode']=result.returncode
        paths=list(out.glob('published_ephemeris*.pfd'));ps=list(out.glob('published_ephemeris*.pfd.ps'))
        assert result.returncode==0 and len(paths)==1 and len(ps)==1,(result.returncode,paths,ps)
        folded=pfd(str(paths[0]));assert folded.npart==64 and folded.nsub==32
        assert math.isfinite(float(folded.profs.sum()))
        record['profile_shape']=list(folded.profs.shape)
        best=Path(str(paths[0])+'.bestprof');assert best.exists()
        record['bestprof_header']=[s for s in best.read_text().splitlines() if s.startswith('#')]
        pdf,png=helper().render_ps(ps[0]);record['pdf']=str(pdf);record['png']=str(png)
        record['sigma']=None;record['reduced_chi2']=None
        for line in record['bestprof_header']:
            match=re.search(r'\(~?\s*([\d.eE+-]+)\s+sigma\)',line)
            if match:record['sigma']=float(match[1])
            if line.startswith('# Reduced chi-sqr'):record['reduced_chi2']=float(line.split('=')[-1])
        assert sha(row['par_original'])==row['par_sha256']
        assert FIL.stat().st_size==m['filterbank_size'] and FIL.stat().st_mtime_ns==m['filterbank_mtime_ns']
        record['wall_s']=time.perf_counter()-started;record['completed_utc']=datetime.now(timezone.utc).isoformat()
        record['products']={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file() and not p.is_symlink() and p.name not in ('running.json','complete.json')}
        save(marker,record);print('Complete',target,'wall_s',record['wall_s'],'sigma',record['sigma'],flush=True)
    except Exception:
        record['traceback']=traceback.format_exc();record['wall_s']=time.perf_counter()-started
        save(out/'failure.json',record);raise

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','fold']);p.add_argument('--target');args=p.parse_args()
    prepare() if args.mode=='prepare' else run(args.target)
