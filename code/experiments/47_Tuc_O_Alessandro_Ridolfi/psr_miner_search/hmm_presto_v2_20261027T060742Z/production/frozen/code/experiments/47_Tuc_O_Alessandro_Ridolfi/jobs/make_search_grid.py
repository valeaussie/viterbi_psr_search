#!/usr/bin/env python3
"""Freeze the DM/Nt/subband grid or verify every measured search completed."""
import argparse, hashlib, json
from pathlib import Path

def values(path, cast):
    return [cast(s.strip()) for s in path.read_text().splitlines() if s.strip() and not s.lstrip().startswith('#')]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--jobs-dir',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--validate-results',type=Path)
    a=p.parse_args()
    if a.validate_results:
        manifest=json.loads(a.out.read_text()); failures=[]; total_cpu=total_wall=0.; hosts={}
        for task in manifest['tasks']:
            base=a.validate_results/task['curve'].replace('_loglike_curve.dat','')
            try:
                if not Path(str(base)+'_loglike_curve.dat').is_file(): raise ValueError('missing curve')
                m=json.loads(Path(str(base)+'_timing.json').read_text())
                s=json.loads(Path(str(base)+'_sampling.json').read_text())
                if m['returncode']!=0 or s['actual_segments']!=task['Nt'] or s['frequency_axis_version']!=2 or s['mode']!='exact_segments':
                    raise ValueError('failed or incorrect search metadata')
                total_cpu+=m['cpu_s']; total_wall+=m['wall_s']
                key=m['cpu_model']; h=hosts.setdefault(key,dict(tasks=0,cpu_s=0.,wall_s=0.))
                h['tasks']+=1; h['cpu_s']+=m['cpu_s']; h['wall_s']+=m['wall_s']
            except Exception as exc: failures.append(dict(index=task['index'],error=str(exc)))
        report=dict(expected=len(manifest['tasks']),failed=len(failures),failures=failures,
            total_cpu_s=total_cpu,sum_task_wall_s=total_wall,by_cpu_model=hosts)
        (a.validate_results/'grid_validation.json').write_text(json.dumps(report,indent=2)+'\n')
        if failures: raise ValueError(f"{len(failures)} incomplete or invalid grid tasks")
        print(json.dumps(report,indent=2)); return
    if a.out.exists(): raise FileExistsError(a.out)
    dms=values(a.jobs_dir/'dm_list.txt',float); nts=values(a.jobs_dir/'nt_list.txt',int)
    if len(set(dms))!=len(dms) or len(set(nts))!=len(nts): raise ValueError('Duplicate grid entries')
    tasks=[]
    for dm in dms:
        for nt in nts:
            for sub in range(50,1050,9):
                tasks.append(dict(index=len(tasks),dm=dm,Nt=nt,subband_f0=sub,
                    curve=f'DM{dm:05.2f}/Nt{nt}/f0_{sub:.3f}/blind_Nt{nt}_f0_{sub:.3f}_loglike_curve.dat'))
    manifest=dict(dm_values=dms,Nt_values=nts,band_start_hz=50,band_start_stop_exclusive_hz=1050,
        band_width_hz=10,band_step_hz=9,effective_frequency_hz=[50,1059],frame='bary',
        grid_size=len(tasks),tasks=tasks)
    a.out.write_text(json.dumps(manifest,indent=2)+'\n')
    print(f'Frozen {len(tasks)} grid configurations in {a.out}')
if __name__=='__main__': main()
