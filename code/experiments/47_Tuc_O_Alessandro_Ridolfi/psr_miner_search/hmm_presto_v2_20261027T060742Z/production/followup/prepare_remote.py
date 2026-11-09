from pathlib import Path
import csv,hashlib,json,os,subprocess,collections,ast,datetime
F=Path(__file__).resolve().parent;D=F.parent;s=json.loads((D/'submitted_production.json').read_text());C=Path(s['candidate_output']);H=Path(s['hmm_output'])
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
rows=read(C/'candidates_timescale_alternatives.csv');cl={r['candidate_id']:r for r in read(C/'classified_timescale_alternatives.csv')};flags={r['candidate_id']:r for r in read(C/'rfi_flags.csv')}
assert len(rows)==len(cl)==len(flags)==415
assert hashlib.sha256((C/'candidates_timescale_alternatives.csv').read_bytes()).hexdigest()==json.loads((C/'candidate_manifest.json').read_text())['table_sha256']
groups=collections.defaultdict(list)
for r in rows:
 cid=r['candidate_id'];r=dict(r,**{k:v for k,v in cl[cid].items() if k not in r});r['rfi_flag']=flags[cid]['rfi_flag']
 p=C/cid/(cid+'_track.dat')
 if p.exists():raise FileExistsError('Existing candidate track; refuse overwrite: '+str(p))
 groups[r['source_curve']].append(r)
sources=[dict(source_curve=k,source_params=v[0]['source_params'],candidates=v) for k,v in sorted(groups.items())]
assert len(sources)==166
manifest=dict(followup_root=str(F),blind_dir=str(H),candidate_dir=str(C),frozen_code=str(D/'frozen/code'),sources=sources,expected_candidates=415,source_table_sha256=hashlib.sha256((C/'candidates_timescale_alternatives.csv').read_bytes()).hexdigest())
(F/'track_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
for p in F.iterdir():
 if p.suffix=='.py':ast.parse(p.read_text())
 if p.suffix=='.sbatch':subprocess.check_call(['bash','-n',str(p)])
env=os.environ.copy();env['FOLLOWUP_ROOT']=str(F)
cmd=['sbatch','--parsable','--account=oz022','--array=0-165','--export=ALL','--output='+str(F/'recover-%A_%a.out'),'--error='+str(F/'recover-%A_%a.err'),str(F/'recover.sbatch')]
job=subprocess.check_output(cmd,env=env,text=True).strip().split(';')[0]
(F/'submitted_tracks.json').write_text(json.dumps(dict(recovery_job_id=job),indent=2)+'\n')
cmd2=['sbatch','--parsable','--account=oz022','--kill-on-invalid-dep=yes','--dependency=afterok:'+job,'--export=ALL','--output='+str(F/'gallery-%j.out'),'--error='+str(F/'gallery-%j.err'),str(F/'gallery.sbatch')]
gallery=subprocess.check_output(cmd2,env=env,text=True).strip().split(';')[0]
record=dict(recovery_job_id=job,gallery_job_id=gallery,sources=166,candidates=415,followup_root=str(F),candidate_root=str(C),gallery=str(C/'index.html'),commands=[cmd,cmd2],scope='Exact original blind candidate paths, no targeted redetection, no changed search bands',utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
(F/'submitted_tracks.json').write_text(json.dumps(record,indent=2)+'\n')
(F/'code_manifest.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in F.iterdir() if p.is_file() and p.suffix in ('.py','.sbatch')},indent=2)+'\n')
print(json.dumps(record,indent=2))
