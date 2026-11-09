"""Fold the saved HMM circular frequency curve without a global Taylor expansion.

TEMPO is used only for the clock, dispersion and barycentric time mapping.
The saved sinusoidal frequency curve is integrated analytically, then encoded
in short polyco intervals, checked against that analytic phase and an
independent shorter interpolation. No published O orbital model is used.
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
import subprocess
import time
import traceback

getcontext().prec=40
W=Path(__file__).resolve().parents[1]
J=W.parent


def save(path,value):
    with Path(path).open('x') as out:json.dump(value,out,indent=2,allow_nan=False);out.write('\n')


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()


def loadmodule(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def helper():return loadmodule(W/'code/native_helper.py','J_full_helper')


def prepare():
    h=helper();old=json.loads((J/'search_manifest.json').read_text())
    selected=json.loads((J/'selected_candidate.json').read_text()); search=next(r for r in old['rows'] if r['index']==selected['search_index'])
    directory=Path(search['directory']);fit=Path(search['prefix']+'_fit_summary.dat');v=h.summary(fit.read_text())
    orbit={k:float(v[name][0]) for k,name in [('spin_hz','kepler_f_spin_hz'),('amplitude_hz','kepler_K_hz'),('period_s','kepler_Pb_s'),('zero_crossing_s','kepler_T0_s'),('reference_s','kepler_tref_s')]}
    orbit['spin_formal_error_hz']=float(v['kepler_f_spin_hz'][1])
    assert orbit['spin_hz']>0 and orbit['amplitude_hz']>=0 and orbit['period_s']>0
    source=Path(old['filterbank']);meta=old['filterbank_metadata']
    assert source.stat().st_size==meta['size_bytes'] and source.stat().st_mtime_ns==meta['mtime_ns']
    row=dict(candidate_id=search['candidate_id'],dm=search['dm'],ra=search['ra'],dec=search['dec'],
        input_epoch_mjd=search['input_epoch_mjd'],reference_epoch_mjd=search['input_epoch_mjd'],
        coefficient_strings=['1024','0','0','0','0','0'])
    save(W/'manifest.json',dict(created_utc=datetime.now(timezone.utc).isoformat(),orbit=orbit,time_mapping_row=row,
        fit=str(fit),fit_sha256=sha(fit),track=search['prefix']+'_track.dat',track_sha256=sha(search['prefix']+'_track.dat'),
        filterbank=str(source),filterbank_metadata=meta,search_row=search,mapping_clock_hz=1024.,
        source='Saved targeted HMM track circular frequency fit; no catalogue orbital parameters or refitting.',
        purpose='Replace the single full-observation frequency Taylor expansion with the analytic periodic phase.',
        code_sha256={p.name:sha(p) for p in (W/'code').glob('*.py')},
        variants=['fixed','standard_refinement','frequency_refinement']))


def evaluate(blocks,epoch):
    candidates=[b for b in blocks if abs(epoch-b['midpoint'])*1440<Decimal(b['span_minutes'])/2]
    if not candidates:raise ValueError('Missing polyco coverage: '+str(epoch))
    b=min(candidates,key=lambda b:abs(epoch-b['midpoint']));x=(epoch-b['midpoint'])*1440
    c=[Decimal(str(v)) for v in b['coefficients']]
    correction=c[-1];derivative=Decimal(0)
    for value in c[-2::-1]:
        derivative=derivative*x+correction
        correction=correction*x+value
    phase=b['phase']+x*60*Decimal(str(b['frequency']))+correction
    freq=Decimal(str(b['frequency']))+derivative/60
    return phase,float(freq),b


def make_mapping(m):
    h=helper();row=m['time_mapping_row'];header=m['filterbank_metadata']['header']
    start=float(header['tstart']);duration=m['filterbank_metadata']['duration_s'];stop=start+duration/86400
    epoch=float(row['input_epoch_mjd']);margin=600/86400
    bary=h.tempo_polycos(W/'mapping_bary',row,'@',0.,0.,epoch-margin,epoch+duration/86400+margin)
    zero=evaluate(bary,Decimal(row['input_epoch_mjd']))[0]
    err=0.
    for i in range(1001):
        t=Decimal(row['input_epoch_mjd'])+Decimal(str(duration))*i/Decimal(1000*86400)
        p,f,_=evaluate(bary,t);expected=(t-Decimal(row['input_epoch_mjd']))*86400
        err=max(err,abs(float((p-zero)/Decimal(str(m['mapping_clock_hz']))-expected)));assert abs(f/m['mapping_clock_hz']-1)<1e-10
    assert err<1e-6,err
    frequencies=[float(header['fch1'])+k*float(header['foff']) for k in range(int(header['nchans']))];fobs=(min(frequencies)+max(frequencies))/2
    topo=h.tempo_polycos(W/'mapping_topo',row,'m',fobs,row['dm'],start-margin,stop+margin)
    # Independent TEMPO interpolation of the same unit-frequency time mapping.
    checkdir=W/'mapping_topo_15coeff';checkdir.mkdir()
    (checkdir/'model.par').write_text(h.par_text(row,row['dm']))
    (checkdir/'tz.in').write_text('m 12 30 15 430\n\n\n0000+0000 30 15 12 %.9f\n'%fobs)
    with (checkdir/'tempo.log').open('w') as out:
        subprocess.run(['tempo','-z','-f','model.par'],cwd=checkdir,input='%d %d\n'%(math.floor(start-margin)-1,math.ceil(stop+margin)),text=True,stdout=out,stderr=subprocess.STDOUT,check=True)
    independent=h.read_polycos(checkdir/'polyco.dat',start-margin,stop+margin,checkdir/'observation.polycos')
    engine=loadmodule(W/'code/reference_engine.py','J_full_numerical_check')
    interval=engine.phase_extrema_check(topo,independent,start,stop)
    assert interval['maximum_phase_difference_cycles']<1e-6,interval
    save(W/'mapping_validation.json',dict(barycentric_linear_time_error_s=err,synthetic_clock_hz=m['mapping_clock_hz'],
        topo_12_vs_15coefficient_interval_check=interval,phase_zero_at_barycentric_input_epoch=str(zero),
        inferred_barycentric_seconds_at_filterbank_start=float(evaluate(topo,Decimal(str(start)))[0]-zero)/m['mapping_clock_hz'],
        observing_frequency_mhz=fobs))
    return topo,zero,fobs


def analytic(m,mapping,zero,epoch):
    p,rate,_=evaluate(mapping,epoch);seconds=float(p-zero)/m['mapping_clock_hz'];rate/=m['mapping_clock_hz']
    o=m['orbit'];w=2*math.pi/o['period_s'];x=seconds-o['reference_s']
    phase=o['spin_hz']*x+o['amplitude_hz']/w*(math.cos(w*(seconds-o['zero_crossing_s']))-math.cos(w*(o['reference_s']-o['zero_crossing_s'])))
    frequency=(o['spin_hz']-o['amplitude_hz']*math.sin(w*(seconds-o['zero_crossing_s'])))*rate
    return phase,frequency


def interpolate(m,mapping,zero,fobs,span,ncoeff,path):
    import numpy as np
    from numpy.polynomial.chebyshev import chebfit,cheb2poly
    start=Decimal(str(m['filterbank_metadata']['header']['tstart']));duration=m['filterbank_metadata']['duration_s']
    step=span*60*.8;half=span/2;lines=[];offsets={};max_fit_residual=0.
    for k in range(math.ceil(duration/step)+1):
        mid=start+Decimal(str(step*k))/86400
        phase,frequency=analytic(m,mapping,zero,mid)
        nodes=np.cos(np.pi*(np.arange(2*ncoeff)+.5)/(2*ncoeff));values=[]
        for node in nodes:
            dt=float(node)*half*60;p,_=analytic(m,mapping,zero,mid+Decimal(str(dt))/86400)
            values.append(p-phase-frequency*dt)
        power=cheb2poly(chebfit(nodes,np.asarray(values),ncoeff-1))
        coefficients=[float(c)/half**i for i,c in enumerate(power)]
        cycle=math.floor(phase);fraction=phase-cycle;midstring=str(mid.quantize(Decimal('0.000000000000001')))
        offsets[midstring]=cycle
        fraction_day=float(mid%1);seconds_day=fraction_day*86400;hh=int(seconds_day//3600);mm=int((seconds_day-hh*3600)//60);ss=seconds_day-hh*3600-mm*60
        utc=hh*10000+mm*100+ss
        line1='OFULL0001 05-NOV-26 %.2f %s %.6f 0.000 0.000'%(utc,midstring,float(m['time_mapping_row']['dm']))
        line2=f'{fraction:20.12f}{frequency:18.12f}{0:5d}{span:5d}{ncoeff:5d}{fobs:10.3f}'+' '*16
        assert len(line1)<=88 and len(line2)==79,(len(line1),len(line2))
        lines.extend([line1+'\n',line2+'\n'])
        for i in range(0,ncoeff,3):lines.append(''.join('%25.17e'%c for c in coefficients[i:i+3])+'\n')
    path.write_text(''.join(lines));blocks=helper().read_polycos(path)
    return blocks,offsets


def build():
    m=json.loads((W/'manifest.json').read_text());assert sha(m['fit'])==m['fit_sha256']
    started=time.perf_counter();mapping,zero,fobs=make_mapping(m)
    first,offsets=interpolate(m,mapping,zero,fobs,10,12,W/'full_orbit.polycos')
    second,other_offsets=interpolate(m,mapping,zero,fobs,5,15,W/'independent_full_orbit.polycos')
    duration=m['filterbank_metadata']['duration_s'];start=Decimal(str(m['filterbank_metadata']['header']['tstart']))
    phase_error=frequency_error=independent_error=0.;freqs=[]
    for i in range(10001):
        epoch=start+Decimal(str(duration))*i/Decimal(10000*86400)
        expected,expected_f=analytic(m,mapping,zero,epoch)
        phase,f,b=evaluate(first,epoch);phase+=offsets[str(b['midpoint'])]
        other,_,c=evaluate(second,epoch);other+=other_offsets[str(c['midpoint'])]
        phase_error=max(phase_error,abs(float(phase)-expected));frequency_error=max(frequency_error,abs(f-expected_f));independent_error=max(independent_error,abs(float(phase-other)))
        freqs.append(f)
    result=dict(samples=10001,max_phase_error_cycles=phase_error,max_frequency_error_hz=frequency_error,
        independent_interpolation_phase_error_cycles=independent_error,phase_tolerance_cycles=1e-4,
        frequency_tolerance_hz=1e-7,positive_frequency_range_hz=[min(freqs),max(freqs)],
        production_span_minutes=10,production_coefficients=12,validation_span_minutes=5,validation_coefficients=15,
        method='Analytic integral of the saved sinusoidal HMM frequency model evaluated through a validated 1024-Hz synthetic-clock TEMPO time mapping; no global Taylor approximation.',
        wall_s=time.perf_counter()-started,model_parameters_changed=bool(m.get('refinement')),source_fit_sha256=sha(m['fit']))
    result['passed']=phase_error<1e-4 and frequency_error<1e-7 and independent_error<1e-4 and min(freqs)>0
    save(W/'phase_validation.json',result);assert result['passed'],result
    save(W/'polyco_cycle_offsets.json',dict(production=offsets,independent=other_offsets))
    print(json.dumps(result),flush=True)


def fold(mode):
    import numpy as np
    from presto import psr_utils
    def rotate(arr,bins):
        arr=np.asarray(arr);f=np.arange(arr.size//2+1,dtype=float)
        return np.fft.irfft(np.exp(2j*np.pi*f*bins/arr.size)*np.fft.rfft(arr),n=arr.size)
    psr_utils.fft_rotate=rotate
    m=json.loads((W/'manifest.json').read_text());validation=json.loads((W/'phase_validation.json').read_text());assert validation['passed']
    assert sha(m['fit'])==m['fit_sha256'] and sha(m['track'])==m['track_sha256']
    source=Path(m['filterbank']);meta=m['filterbank_metadata'];assert source.stat().st_size==meta['size_bytes'] and source.stat().st_mtime_ns==meta['mtime_ns']
    out=W/'folds'/mode;out.mkdir()
    (out/'input.fil').symlink_to(source);(out/'model.polycos').symlink_to(W/'full_orbit.polycos')
    command=['prepfold','-filterbank','-topo','-noxwin','-ncpus','1','-dm',m.get('fold_dm',m['time_mapping_row']['dm']),'-npart','64','-nsub','32']
    if mode=='fixed':command+=['-nosearch']
    elif mode=='standard_refinement':command+=['-npfact','2','-ndmfact','3','-searchfdd']
    elif mode=='frequency_refinement':command+=['-npfact','2','-ndmfact','3','-nopdsearch']
    else:raise ValueError(mode)
    command+=['-psr','OFULL0001','-polycos','model.polycos','-o','full_orbit','input.fil']
    record=dict(mode=mode,command=command,started_utc=datetime.now(timezone.utc).isoformat(),job_id=os.environ.get('SLURM_JOB_ID'),
        folding_orbit=m['orbit'],refinement=m.get('refinement'),input_fit_sha256=m['fit_sha256'],polycos_sha256=sha(W/'full_orbit.polycos'),phase_validation=validation,
        worker_sha256=sha(__file__),spin_only_wide_search=(mode=='frequency_refinement'),published_orbit_used=False)
    started=time.perf_counter();save(out/'running.json',record)
    try:
        with (out/'prepfold.log').open('w') as log:subprocess.run(command,cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True)
        from presto.prepfold import pfd
        path=next(out.glob('full_orbit*.pfd'));folded=pfd(str(path));assert folded.profs.shape[:2]==(64,32) and math.isfinite(float(folded.profs.sum()))
        record['profile_shape']=list(folded.profs.shape);best=Path(str(path)+'.bestprof')
        record['bestprof_header']=[s for s in best.read_text().splitlines() if s.startswith('#')]
        record['sigma']=next(float(q[1]) for line in record['bestprof_header'] if (q:=re.search(r'\(~?\s*([\d.eE+-]+)\s+sigma\)',line)))
        ps=next(out.glob('full_orbit*.pfd.ps'));helper().render_ps(ps)
        record['wall_s']=time.perf_counter()-started;record['completed_utc']=datetime.now(timezone.utc).isoformat();record['returncode']=0
        record['products']={p.name:sha(p) for p in out.iterdir() if p.is_file() and not p.is_symlink() and p.name!='running.json'}
        assert sha(m['fit'])==m['fit_sha256'] and sha(m['track'])==m['track_sha256']
        save(out/'complete.json',record);print(json.dumps(dict(mode=mode,sigma=record['sigma'],wall_s=record['wall_s'])),flush=True)
    except Exception:
        record['error']=traceback.format_exc();record['wall_s']=time.perf_counter()-started;save(out/'failure.json',record);raise


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','build','fixed','standard_refinement','frequency_refinement']);a=p.parse_args()
    if a.mode=='prepare':prepare()
    elif a.mode=='build':build()
    else:fold(a.mode)
