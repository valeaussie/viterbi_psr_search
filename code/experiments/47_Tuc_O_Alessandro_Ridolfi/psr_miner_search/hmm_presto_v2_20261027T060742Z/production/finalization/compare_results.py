"""Describe candidate associations and measured costs at the declared operating points."""
import collections,csv,json,subprocess,sys
from pathlib import Path
D=Path(__file__).resolve().parent.parent
submission=json.loads((D/'submitted_production.json').read_text())
C=Path(submission['candidate_output']);H=Path(submission['hmm_output'])
P=D/('paired_'+submission['paired_job_id'])/'presto'
J=D/'frozen/code/experiments/47_Tuc_O_Alessandro_Ridolfi/jobs'
K=D/'frozen/code/experiments/47_Tuc_O_Alessandro_Ridolfi/47Tuc_blind_search_v1/config/known_47tuc_pulsars.yaml'
O=D/'comparison';O.mkdir(exist_ok=True)
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def write(p,rows,fields):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
def call(name,*args):subprocess.check_call([sys.executable,str(J/name),*map(str,args)])
# The final PRESTO tables are additional parsed products; original candidate products remain untouched.
call('classify_known.py','--cands',P/'candidates_all.csv','--known-yaml',K,'--out',P/'classified_all.csv')
call('flag_rfi.py','--classified-csv',P/'classified_all.csv','--known-yaml',K,'--cand-dir',P,'--out-csv',P/'rfi_flags.csv')
manifest=json.loads((C/'candidate_manifest.json').read_text()); selected=read(C/'classified_timescale_alternatives.csv')
assert [r['candidate_id'] for r in selected]==manifest['candidate_ids']
raw=read(C/'candidates_raw.csv');members=read(C/'candidate_members.csv')
assert {r['peak_id'] for r in raw}=={r['peak_id'] for r in members}
assert len(raw)==len(members)==manifest['raw_peak_count']
assert len({r['peak_id'] for r in members})==len(raw)
group_by_peak={r['peak_id']:r['candidate_id'] for r in members}
raw_for_match=[dict(r,dm_best=r['dm'],candidate_id=r['peak_id'],group_candidate_id=group_by_peak[r['peak_id']]) for r in raw]
rawtable=O/'hmm_raw_peak_catalogue_input.csv'
write(rawtable,raw_for_match,list(dict.fromkeys([k for r in raw_for_match for k in r])) or ['candidate_id','dm_best','peak_freq_hz'])
call('classify_known.py','--cands',rawtable,'--known-yaml',K,'--out',O/'hmm_raw_peak_classified.csv')
raw_matched=read(O/'hmm_raw_peak_classified.csv'); presto=read(P/'classified_all.csv')
expected_presto=json.loads((P/'summary.json').read_text())['retained_candidates'];assert len(presto)==expected_presto
# Keep representative and raw-member associations separate. A member's best DM can differ from its group's representative.
reps={r['candidate_id']:r for r in selected}; extras=[]
for r in raw_matched:
 rep=reps[r['group_candidate_id']]
 if r['classification']=='known_main' and (rep['classification']!='known_main' or rep['parent_pulsar']!=r['parent_pulsar']):extras.append(r)
write(O/'raw_member_associations_not_in_representative.csv',extras,list(raw_matched[0]) if raw_matched else ['candidate_id'])
import yaml
catalogue=yaml.safe_load(K.read_text())['pulsars']
rows=[]
for entry in catalogue:
 name=entry.get('name','')
 def count(table,cls,**filters):
  return sum(r['parent_pulsar']==name and r['classification']==cls and all(r.get(k)==str(v) for k,v in filters.items()) for r in table)
 r=dict(pulsar=name,hmm_raw_fundamental_associations=count(raw_matched,'known_main'),hmm_group_fundamental_associations=count(selected,'known_main'),hmm_raw_harmonic_associations=count(raw_matched,'known_harmonic'))
 for parts in [1,8,16]:r[f'presto_parts{parts}_fundamental_associations']=count(presto,'known_main',parts=parts)
 for parts in [1,8,16]:r[f'presto_parts{parts}_harmonic_associations']=count(presto,'known_harmonic',parts=parts)
 rows.append(r)
write(O/'catalogue_candidate_associations.csv',rows,list(rows[0]))
grid=json.loads((H/'grid_validation.json').read_text());assert grid['failed']==0 and grid['expected']==16464
pt=read(P/'presto_timings.csv');assert len(pt)==525
costs=dict(hmm_grid=grid,presto_by_parts={str(parts):dict(sections=sum(r['parts']==str(parts) for r in pt),
 total_cpu_s=sum(float(r['total_cpu_s']) for r in pt if r['parts']==str(parts)),
 sum_section_wall_s=sum(float(r['total_wall_s']) for r in pt if r['parts']==str(parts)),
 preparation_cpu_s=sum(float(r['preparation_cpu_s']) for r in pt if r['parts']==str(parts)),
 search_cpu_s=sum(float(r['search_cpu_s']) for r in pt if r['parts']==str(parts))) for parts in [1,8,16]},
 presto_machine=json.loads((P/'machine.json').read_text()),
 scopes='Whole search cost including preparation; excludes dedispersion and follow-up. CPU totals and sums of command elapsed times are distinct from whole-job elapsed time.')
(O/'measured_search_costs.json').write_text(json.dumps(costs,indent=2)+'\n')
flags_h=read(C/'rfi_flags.csv');flags_p=read(P/'rfi_flags.csv')
assert len(flags_h)==len(selected) and len(flags_p)==len(presto)
assert all(r['retained']=='True' for r in flags_h+flags_p)
counts=dict(status='complete',hmm_raw_selected_peaks=len(raw),hmm_selected_timescale_groups=len(selected),
 presto_final_candidates=len(presto),hmm_group_classes=dict(collections.Counter(r['classification'] for r in selected)),
 hmm_raw_classes=dict(collections.Counter(r['classification'] for r in raw_matched)),
 presto_classes=dict(collections.Counter(r['classification'] for r in presto)),
 rfi_policy='flag_only',hmm_frequency_flags=sum(bool(r['rfi_flag']) for r in flags_h),presto_frequency_flags=sum(bool(r['rfi_flag']) for r in flags_p),
 association_warning='Catalogue associations at declared thresholds; folding and inspection remain necessary to confirm recoveries. Thresholds and signal models are not matched.',
 raw_member_associations_not_represented=len(extras),hmm_grid_complete=True)
(O/'summary.json').write_text(json.dumps(counts,indent=2)+'\n')
ids=[submission[k] for k in ['paired_job_id','hmm_job_id','aggregate_job_id']]
(O/'slurm_accounting.psv').write_text(subprocess.check_output(['sacct','-j',','.join(ids),'-P','--format=JobID,JobName,State,Elapsed,TotalCPU,AllocCPUS,CPUTime,MaxRSS,NodeList,ExitCode'],text=True))
print(json.dumps(counts,indent=2))
