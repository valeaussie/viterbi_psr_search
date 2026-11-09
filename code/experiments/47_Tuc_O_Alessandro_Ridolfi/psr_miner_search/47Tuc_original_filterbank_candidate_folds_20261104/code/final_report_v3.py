"""Final status after every candidate worker ends, including killed jobs.

Read independent small records concurrently to avoid serial Lustre latency.
This limits metadata-reading threads only, never Slurm job concurrency.
"""
from pathlib import Path
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import csv
import io
import json
import subprocess
W=Path(__file__).resolve().parents[1]
def main():
    m=json.loads((W/'manifest.json').read_text());state=json.loads((W/'submitted_production.json').read_text())
    replacement=W/'production_replacements.json'
    if replacement.exists():state['jobs']=state['jobs']+json.loads(replacement.read_text())['jobs']
    jobids=[j['job_id'] for j in state['jobs']]
    accounting=subprocess.check_output(['sacct','-X','--array','-nP','-j',','.join(jobids),'-o','JobID%60,State,ExitCode,ElapsedRaw,TimelimitRaw'],text=True)
    (W/'final_slurm_accounting.txt').write_text(accounting)
    slurm={}
    for line in accounting.splitlines():
        parts=line.split('|')
        if len(parts)>=5:slurm[parts[0]]=dict(job_state=parts[1],exit_code=parts[2],elapsed_s=parts[3],requested_minutes=parts[4])
    assignments={}
    for j in state['jobs']:
        for pos,index in enumerate(j['indices']):assignments[index]=str(j['job_id'])+'_'+str(pos)
    def read(r):
        path=W/'records'/('%05d.json'%r['index'])
        rec=json.loads(path.read_text()) if path.exists() else {}
        status=rec.get('status','not_started');job=assignments[r['index']];acct=slurm.get(job,{})
        if status in ('running','not_started'):
            status='incomplete_'+acct.get('job_state','worker_finished_without_complete_record').lower().replace(' ','_')
        if status=='complete' and not (Path(rec['pdf']).is_file() and Path(rec['png']).is_file()):status='missing_native_plot'
        return dict(index=r['index'],pipeline=r['pipeline'],candidate_id=r['candidate_id'],model=r['model'],
            status=status,nominal_sigma=rec.get('sigma'),wall_s=rec.get('wall_s'),pdf=rec.get('pdf'),png=rec.get('png'),
            record=str(path),slurm_job=job,slurm_state=acct.get('job_state'),error=rec.get('error','')[-1500:])
    with ThreadPoolExecutor(max_workers=32) as pool:rows=list(pool.map(read,m['rows']))
    result=dict(updated_utc=datetime.now(timezone.utc).isoformat(),filterbank=m['filterbank'],expected=len(rows),
        status=dict(Counter(r['status'] for r in rows)),counts_by_pipeline=dict(Counter(r['pipeline'] for r in rows)),
        all_workers_finished=True,needs_followup=any(r['status'] not in ('complete','invalid_saved_model') for r in rows),rows=rows)
    temp=W/'results.json.tmp';temp.write_text(json.dumps(result,indent=2)+'\n');temp.replace(W/'results.json')
    stream=io.StringIO();writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (W/'original_filterbank_fold_index.csv').write_text(stream.getvalue())
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
if __name__=='__main__':main()
