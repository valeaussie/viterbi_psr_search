"""Dispatch the remaining J folds only after every search has completed."""
from pathlib import Path
import importlib.util
import json
import subprocess
import hashlib

W=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('J_workflow_dispatch',W/'code/workflow.py')
w=importlib.util.module_from_spec(spec);spec.loader.exec_module(w)
m=json.loads((W/'search_manifest.json').read_text())
assert all((Path(r['directory'])/'search_complete.json').exists() for r in m['rows'])
assert not (W/'submitted_folds.json').exists()
w.fold_prepare()
pilot={'jobs': []}
already={i for j in pilot['jobs'] for i in j['indices']}
indices=[i for i in range(2*len(m['rows'])) if i not in already]
job=w.submit('fold',indices,28,2048,label='production')
dependencies=':'.join([job['job_id']]+[j['job_id'] for j in pilot['jobs']])
command=['sbatch','--parsable','--account=oz022','--qos=normal','--job-name=J-report',
    '--ntasks=1','--cpus-per-task=1','--mem=512M','--time=00:02:00',
    '--dependency=afterany:'+dependencies,
    '--output='+str(W/'logs/final-report-%j.out'),'--error='+str(W/'logs/final-report-%j.err'),
    '--wrap=python '+str(W/'code/workflow.py')+' status']
report=subprocess.check_output(command,text=True).strip().split(';')[0]
w.save(W/'submitted_folds.json',dict(created_utc=w.now(),jobs=[job],pilot_jobs=pilot['jobs'],
    report_job_id=report,report_command=command,maximum_concurrent_tasks=None,
    dispatcher_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    wall_basis='Measured O original-filterbank native fold: 1317 s, same 35 pulse bins and full PRESTO refinement; 28 minutes includes 25 percent plus startup margin.'))
w.status()
print(json.dumps(dict(fold_job_id=job['job_id'],remaining_fold_tasks=len(indices),report_job_id=report)),flush=True)
