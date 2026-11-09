"""Prepare immutable candidate membership and small folding validation arrays."""
import collections
import datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess

F = Path(__file__).resolve().parent
D = F.parent
tracks = json.loads((F/'track_manifest.json').read_text())
presto = json.loads((F/'presto_parameters_fullprecision.json').read_text())['candidates']
C = Path(tracks['candidate_dir'])
blind = C.parent.parent
hmm = [dict(candidate_id=r['candidate_id'], candidate_directory=str(C/r['candidate_id']))
       for r in sorted([r for group in tracks['sources'] for r in group['candidates']], key=lambda r:r['candidate_id'])]
groups = collections.defaultdict(list)
for row in presto:
    groups[row['source_table']].append(row)
# At most four candidates per worker keeps short resource requests practical.
# Jobs are distributed freely by Slurm, without an extra concurrency cap.
presto_jobs = [dict(source_table=k, candidates=rows[i:i+4])
               for k, rows in sorted(groups.items()) for i in range(0, len(rows), 4)]
manifest = dict(output_root=str(blind/'stage4_fold/comparison_v1'),
                dedispersion_root=str(blind/'stage2_dedisp'), frozen_code=tracks['frozen_code'],
                hmm=hmm, presto=presto_jobs, expected_hmm_candidates=415,
                expected_presto_candidates=5538, expected_models_per_hmm=2,
                rfi_policy='flag_only, preserve every selected candidate',
                folding_policy='Barycentric dedispersed time series, fixed original model plus identical F0-only refinement; 64 phase bins, 64 subintegrations, 1025 refinement trials spanning +/-8 cycles per native interval. PRESTO uses original valid section; HMM uses full valid observation. No catalogue parameters supplied.',
                separation_from_benchmark='All fits, track recovery, folding and refinement are follow-up costs, excluded from blind-search timing.')
assert len(hmm)==415 and len({r['candidate_id'] for r in hmm})==415
assert sum(len(j['candidates']) for j in presto_jobs)==5538
(F/'fold_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
rows = [json.loads((Path(r['candidate_directory'])/'candidate.json').read_text()) for r in hmm]
best_j = min([v for v in rows if v['classification']=='known_main' and v['parent_pulsar']=='J0024-7204J' and v['Nt_best']=='64'], key=lambda v:float(v['freq_residual_hz']))
test_hmm = [next(r for r in hmm if r['candidate_id']==best_j['candidate_id']),
            next(r for r in hmm if r['candidate_id']==next(v['candidate_id'] for v in rows
             if v['classification']=='unclassified' and v['Nt_best']=='16'))]
short = sorted([r for r in presto if r['classification']=='known_main' and r['parts']=='16'
                and r['parent_pulsar']=='J0024-7204J'], key=lambda r:-float(r['sigma']))[0]
full = next(r for r in presto if r['parts']=='1' and r['classification']=='unclassified')
pilot = dict(manifest, output_root=str(F/'fold_pilot'), hmm=test_hmm,
             presto=[dict(candidates=[short]), dict(candidates=[full])])
(F/'fold_pilot_manifest.json').write_text(json.dumps(pilot, indent=2)+'\n')
subprocess.check_call(['bash', '-n', str(F/'fold_worker.sbatch')])
env = os.environ.copy()
env['FOLLOWUP_ROOT'] = str(F)
env['FOLD_MANIFEST'] = str(F/'fold_pilot_manifest.json')
jobs = {}
for method in ['hmm', 'presto']:
    env['FOLD_METHOD'] = method
    command = ['sbatch', '--parsable', '--account=oz022', '--array=0-1', '--export=ALL',
               '--output='+str(F/('fold-pilot-'+method+'-%A_%a.out')),
               '--error='+str(F/('fold-pilot-'+method+'-%A_%a.err')), str(F/'fold_worker.sbatch')]
    jobs[method] = subprocess.check_output(command, env=env, text=True).strip().split(';')[0]
record = dict(jobs=jobs, utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
              production_workers=dict(hmm=len(hmm), presto=len(presto_jobs)),
              pilot_hmm=[r['candidate_id'] for r in test_hmm],
              pilot_presto=[short['candidate_id'], full['candidate_id']],
              manifest_sha256=hashlib.sha256((F/'fold_manifest.json').read_bytes()).hexdigest())
(F/'submitted_fold_pilot.json').write_text(json.dumps(record, indent=2)+'\n')
print(json.dumps(record, indent=2))
