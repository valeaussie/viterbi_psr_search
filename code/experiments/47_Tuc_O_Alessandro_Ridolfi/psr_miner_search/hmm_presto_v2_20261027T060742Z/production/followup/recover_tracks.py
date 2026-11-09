#!/usr/bin/env python3
"""Recover selected original HMM paths without changing the blind-search grid."""
import argparse,ast,csv,hashlib,html,json,os,subprocess,sys,time
from pathlib import Path
import numpy as np

def backtrace_delta(delta,end):
    """The saved recurrence uniquely supplies the original predecessor choices."""
    nf,nt=delta.shape
    if not 0<=end<nf:raise IndexError(end)
    indices=np.empty(nt,dtype=np.int64);indices[-1]=end
    for t in range(nt-1,0,-1):
        i=int(indices[t]);first=max(0,i-1);last=min(nf,i+2)
        indices[t-1]=first+np.argmax(delta[first:last,t-1])
    return indices

def emissions_from_delta(delta):
    spec=np.empty_like(delta);spec[:,0]=delta[:,0]
    for t in range(1,delta.shape[1]):
        prev=delta[:,t-1]
        maximum=prev.copy()
        maximum[1:]=np.maximum(maximum[1:],prev[:-1])
        maximum[:-1]=np.maximum(maximum[:-1],prev[1:])
        spec[:,t]=delta[:,t]-maximum
    return spec

def load_params(path):
    result={}
    for line in path.read_text().splitlines():
        if '=' not in line or line.lstrip().startswith('#'):continue
        k,v=line.split('=',1);result[k.strip()]=ast.literal_eval(v.strip())
    return result

def plot_candidate(path,frequencies,spec,tsft,row,directory):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    df=float(frequencies[1]-frequencies[0]);margin=max(.03,10*df)
    lo=max(frequencies[0],min(path)-margin);hi=min(frequencies[-1],max(path)+margin)
    mask=(frequencies>=lo)&(frequencies<=hi)
    view=spec[mask];freq=frequencies[mask]
    times=(np.arange(len(path))+.5)*tsft
    fig,(ax,curveax)=plt.subplots(2,1,figsize=(10,7),gridspec_kw={'height_ratios':[3,1.2]})
    ax.imshow(view,origin='lower',aspect='auto',interpolation='nearest',
              extent=[0,len(path)*tsft/3600,freq[0]-df/2,freq[-1]+df/2],
              cmap='gray_r',vmin=np.percentile(view,20),vmax=np.percentile(view,99.5))
    ax.plot(times/3600,path,color='tab:red',lw=1,drawstyle='steps-mid',label='Selected Viterbi path')
    ax.set(xlabel='Time from observation start (hours)',ylabel='Barycentric frequency (Hz)')
    ax.ticklabel_format(axis='y',style='plain',useOffset=False)
    ax.legend(loc='best',fontsize=8)
    cid=row['candidate_id'];classification=row.get('classification','unclassified')
    flag=row.get('rfi_flag','') or 'none'
    fig.suptitle(f"{cid}\n{classification} {row.get('parent_pulsar','')} | RFI frequency flag: {flag} (retained)",fontsize=9)
    curve=np.loadtxt(directory/'selected_source_curve.dat',ndmin=2)
    cmask=(curve[:,0]>=lo)&(curve[:,0]<=hi)
    curveax.plot(curve[cmask,0],curve[cmask,1],lw=.8,color='black')
    curveax.axhline(float(row['threshold']),color='tab:blue',lw=.8,ls='--',label='Selection threshold')
    curveax.plot(float(row['peak_freq_hz']),float(row['peak_loglike']),'o',color='tab:red',ms=4)
    curveax.set(xlabel='Final-state frequency (Hz)',ylabel='Final path score')
    curveax.ticklabel_format(axis='x',style='plain',useOffset=False)
    curveax.legend(fontsize=8)
    fig.tight_layout(rect=[0,0,1,.94])
    fig.savefig(directory/(cid+'_spectrogram.png'),dpi=150)
    fig.savefig(directory/(cid+'_spectrogram.pdf'),dpi=300)
    plt.close(fig)

def recover(args):
    manifest=json.loads(args.manifest.read_text());task=manifest['sources'][args.index]
    cache=Path(manifest['followup_root'])/'source_cache'/f'source_{args.index:03d}'
    cache.mkdir(parents=True,exist_ok=True)
    done=cache/'complete.json'
    if done.exists():raise FileExistsError(done)
    original=Path(manifest['blind_dir'])/task['source_curve']
    source_params=Path(manifest['blind_dir'])/task['source_params']
    values=load_params(source_params);values.update(out_prefix=str(cache/'source'),save_delta=True)
    params=cache/'source.params'
    params.write_text('\n'.join(f'{k} = {v!r}' for k,v in values.items())+'\n')
    frozen=Path(manifest['frozen_code']);timer=frozen/'experiments/47_Tuc_O_Alessandro_Ridolfi/jobs/measure_search.py'
    command=[sys.executable,str(timer),'--out',str(cache/'source_timing.json'),'--',sys.executable,str(frozen/'viterbi_pipeline.py'),'--params',str(params),'--lean']
    with (cache/'recovery.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    curve=np.loadtxt(cache/'source_loglike_curve.dat',ndmin=2);old=np.loadtxt(original,ndmin=2)
    if curve.shape!=old.shape or not np.allclose(curve[:,0],old[:,0],rtol=0,atol=1e-10):
        raise ValueError('Recovered frequency grid differs from frozen blind output')
    error=float(np.max(np.abs(curve[:,1]-old[:,1])))
    if error>1e-6:raise ValueError(f'Recovered path scores differ from original: {error}')
    delta=np.loadtxt(cache/'source_delta.dat',ndmin=2)
    if not np.allclose(delta[:,-1],old[:,1],rtol=0,atol=1e-6):raise ValueError('Delta/curve mismatch')
    sampling=json.loads((cache/'source_sampling.json').read_text());tsft=sampling['segment_duration_s'];nt=sampling['actual_segments']
    if delta.shape!=(len(curve),nt):raise ValueError('Incorrect dynamic-programming dimensions')
    spec=emissions_from_delta(delta)
    if np.min(spec)<-1e-7:raise ValueError('Negative reconstructed emission beyond rounding tolerance')
    np.savez_compressed(cache/'original_search_state.npz',delta=delta,spectrogram=spec,frequency_hz=curve[:,0],segment_duration_s=tsft)
    produced=[]
    for row in task['candidates']:
        cid=row['candidate_id'];destination=Path(manifest['candidate_dir'])/cid;destination.mkdir(exist_ok=True)
        end=int(row['source_bin'])
        if abs(curve[end,0]-float(row['peak_freq_hz']))>1e-9 or abs(delta[end,-1]-float(row['peak_loglike']))>1e-6:
            raise ValueError('Candidate endpoint is not its original selected peak: '+cid)
        indices=backtrace_delta(delta,end);frequencies=curve[indices,0];tmid=(np.arange(nt)+.5)*tsft
        np.savetxt(destination/(cid+'_track.dat'),np.c_[tmid,frequencies],header='time_s viterbi_frequency_hz (barycentric; exact original selected final state)')
        np.save(destination/'original_state_indices.npy',indices)
        np.savetxt(destination/'selected_source_curve.dat',old,header='frequency_hz final_path_score (frozen original blind output)')
        metadata=dict(row,exact_original_path=True,source_curve=str(original),source_params=str(source_params),source_bin=end,
            source_recovery_index=args.index,source_curve_max_score_difference=error,sampling=sampling,
            source_cache=str(cache),phase_reference='seconds from barycentric input epoch',
            frequency_frame='barycentric',selected_candidate_retained=True)
        (destination/'candidate.json').write_text(json.dumps(metadata,indent=2)+'\n')
        plot_candidate(frequencies,curve[:,0],spec,tsft,row,destination)
        produced.append(cid)
    done.write_text(json.dumps(dict(status='complete',candidate_ids=produced,source_curve=str(original),max_score_difference=error),indent=2)+'\n')
    print(f'Recovered {len(produced)} exact candidate paths from source {args.index}',flush=True)

def gallery(args):
    m=json.loads(args.manifest.read_text());root=Path(m['candidate_dir']);rows=[r for g in m['sources'] for r in g['candidates']]
    rows.sort(key=lambda r:int(r['cand_index']));errors=[];cards=[]
    for r in rows:
        cid=r['candidate_id'];directory=root/cid
        for suffix in ['_track.dat','_spectrogram.png','_spectrogram.pdf']:
            if not (directory/(cid+suffix)).is_file():errors.append(cid+suffix)
        if not (directory/'candidate.json').is_file():errors.append(cid+'/candidate.json')
        else:
            meta=json.loads((directory/'candidate.json').read_text())
            if not meta['exact_original_path']:errors.append(cid+' non-original path')
        text=f"{cid} {r.get('classification','')} {r.get('parent_pulsar','')} {r.get('rfi_flag','')}"
        esc=html.escape
        cards.append(f'<article data-search="{esc(text.lower(),quote=True)}"><h3>{esc(cid)}</h3><p>{esc(r.get("classification",""))} {esc(r.get("parent_pulsar",""))}; RFI flag: {esc(r.get("rfi_flag","") or "none")}, retained</p><a href="{cid}/{cid}_spectrogram.pdf">PDF</a> · <a href="{cid}/{cid}_track.dat">Track data</a><br><img loading="lazy" src="{cid}/{cid}_spectrogram.png" alt="{esc(cid)}"></article>')
    if errors:raise RuntimeError(f'Incomplete candidate diagnostics ({len(errors)}): {errors[:10]}')
    page='''<!doctype html><meta charset="utf-8"><title>47 Tuc blind HMM candidate tracks</title>
<style>body{font:15px system-ui;margin:24px;max-width:1200px;background:#f6f7f9;color:#18202c}input{position:sticky;top:0;padding:12px;width:95%;font-size:16px}article{background:white;border:1px solid #ddd;border-radius:8px;margin:18px 0;padding:12px}img{max-width:100%;height:auto}h3{overflow-wrap:anywhere;font-size:14px}</style>
<h1>47 Tuc: all 415 selected blind HMM candidate paths</h1><p>Exact paths for the originally selected final states. Every selected candidate is retained, including RFI flags. Catalogue labels are associations, not confirmations.</p>
<input id="q" placeholder="Filter by candidate ID, DM, Nt, pulsar name or RFI flag" oninput="document.querySelectorAll('article').forEach(a=>a.hidden=!a.dataset.search.includes(this.value.toLowerCase()))">
'''+''.join(cards)
    (root/'index.html').write_text(page)
    (root/'track_recovery_audit.json').write_text(json.dumps(dict(status='complete',expected_candidates=len(rows),source_runs=len(m['sources']),all_original_selected_paths=True,excluded_candidates=0,source_manifest=str(args.manifest)),indent=2)+'\n')
    print('All selected candidate tracks and plots are present:',len(rows))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='mode',required=True)
    for name in ['recover','gallery']:
        t=sub.add_parser(name);t.add_argument('--manifest',required=True,type=Path)
        if name=='recover':t.add_argument('--index',required=True,type=int)
    a=p.parse_args();recover(a) if a.mode=='recover' else gallery(a)
