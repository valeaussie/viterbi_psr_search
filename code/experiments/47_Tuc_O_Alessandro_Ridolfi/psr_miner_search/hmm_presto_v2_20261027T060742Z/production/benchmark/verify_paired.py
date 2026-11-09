"""Fail the paired benchmark if any expected run or candidate record is missing."""
import argparse, csv, json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--expected-dms',type=int,required=True);p.add_argument('--expected-repeats',type=int,required=True);a=p.parse_args()
h=json.loads((a.root/'hmm/summary.json').read_text()); s=json.loads((a.root/'presto/summary.json').read_text())
assert h['status']==s['status']=='complete'
assert h['runs']==7*a.expected_repeats and h['corrected_segments_verified']
assert s['sections']==25*a.expected_dms and s['dm_trials']==a.expected_dms
with (a.root/'presto/candidates_all.csv').open(newline='') as f:rows=list(csv.DictReader(f))
assert len(rows)==s['retained_candidates'] and len({r['candidate_id'] for r in rows})==len(rows)
assert all(r['retained']=='True' and r['rfi_policy']=='flag_only' for r in rows)
# Every valid sample belongs to exactly one section at each partition level.
by={}
for path in (a.root/'presto').glob('DM*/parts*/section*/section.section.json'):
 r=json.loads(path.read_text());by.setdefault((r['dm'],r['parts']),[]).append(r)
for (dm,parts),sections in by.items():
 sections.sort(key=lambda r:r['index']); assert len(sections)==parts
 assert sections[0]['start_sample']==0
 assert all(x['end_sample_exclusive']==y['start_sample'] for x,y in zip(sections,sections[1:]))
 assert len({x['source'] for x in sections})==1
assert len(by)==3*a.expected_dms
record=dict(status='validated',hmm_runs=h['runs'],presto_sections=s['sections'],retained_presto_candidates=len(rows),sample_coverage='exact continuous valid input at each partition level',dedispersion_reused=True)
(a.root/'validation.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
