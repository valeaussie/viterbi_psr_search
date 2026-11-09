from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json

W=Path(__file__).resolve().parents[1]
rows=[];records=[]
for mode in ['fixed','standard_refinement','frequency_refinement']:
    p=W/'folds'/mode/'complete.json'
    if not p.exists():rows.append(dict(mode=mode,status='not_complete'));continue
    d=json.loads(p.read_text());checks={}
    for name,value in d['products'].items():checks[name]=hashlib.sha256((p.parent/name).read_bytes()).hexdigest()==value
    assert all(checks.values())
    row=dict(mode=mode,status='complete',sigma=d['sigma'],wall_s=d['wall_s'],all_hashes_passed=True)
    for line in d['bestprof_header']:
        if line.startswith('# Best DM'):row['best_dm']=float(line.split('=')[1])
    for ext in ['png','pdf']:
        row[ext]=str(next(p.parent.glob('*.pfd.'+ext)))
    rows.append(row);records.append(dict(completion_path=str(p),completion=d,product_checks=checks))
result=dict(checked_utc=datetime.now(timezone.utc).isoformat(),rows=rows,records=records,
    phase_validation=json.loads((W/'phase_validation.json').read_text()),
    full_model_source='Saved targeted HMM NT40 DM24.50 circular frequency fit',
    selection='Printed statistics require time, frequency and dispersion diagnostic review.',
    table_changed=False)
(W/'results.json').write_text(json.dumps(result,indent=2)+'\n')
with (W/'fold_summary.csv').open('w',newline='') as out:
    fields=list(dict.fromkeys(k for r in rows for k in r));wr=csv.DictWriter(out,fieldnames=fields);wr.writeheader();wr.writerows(rows)
print(json.dumps(rows,indent=2),flush=True)
