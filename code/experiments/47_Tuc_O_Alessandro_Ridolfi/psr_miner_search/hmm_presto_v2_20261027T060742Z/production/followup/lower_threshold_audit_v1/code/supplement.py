"""Inspect practical Nt16 cuts and locate target peaks in existing audit groups."""
from pathlib import Path
import bisect,csv,json,os,resource,time
import numpy as np

W=Path('/fred/oz022/vdimarco/software/install/viterbi_psr_search/code/experiments/47_Tuc_O_Alessandro_Ridolfi/psr_miner_search/hmm_presto_v2_20261027T060742Z/production/followup/lower_threshold_audit_v1')

def rows(path):return list(csv.DictReader(path.open()))
def save(name,data):
    with (W/name).open('x',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=list(dict.fromkeys(k for r in data for k in r)))
        writer.writeheader();writer.writerows(data)

start=time.perf_counter()
cfg=json.loads((W/'inputs.json').read_text())
assert (W/'complete.json').exists()
families=rows(W/'Nt16/families_above_median.csv')
scores=np.load(W/'Nt16/peaks_above_median.npz')['score']
cuts=[float(cfg['thresholds']['16']),80.,75.,70.,65.,60.,59.,58.5,58.,55.,50.,48.,47.,46.,45.,44.]
counts=[]
for cut in cuts:
    n=sum(float(r['representative_score'])>cut for r in families)
    counts.append(dict(Nt=16,score_cut=cut,raw_peaks_Nt16=int(np.count_nonzero(scores>cut)),
                       groups_Nt16=n,additional_groups_Nt16=n-38,
                       total_groups_if_other_timescales_unchanged=415-38+n))
save('practical_Nt16_cuts.csv',counts)
del scores
selected=rows(W/'target_best_by_timescale.csv')
mapped=[]
for nt in [16,32,64,128,256,512,1024]:
    ff=rows(W/f'Nt{nt}/families_above_median.csv')
    ff.sort(key=lambda r:float(r['representative_hz']))
    freqs=[float(r['representative_hz']) for r in ff]
    for row in selected:
        if int(row['Nt'])!=nt:continue
        freq=float(row['window_peak_hz']);score=float(row['window_peak_score'])
        lo,hi=bisect.bisect_left(freqs,freq-.035),bisect.bisect_right(freqs,freq+.035)
        candidates=ff[lo:hi]
        # Minimum creation ID is the earliest/highest-ranked matching anchor.
        rep=min(candidates,key=lambda r:int(r['family_id'].split('family')[1])) if candidates else None
        out=dict(row)
        if rep is not None and score>float(row['median_at_Nt']):
            out.update(family_id=rep['family_id'],family_hz=rep['representative_hz'],
                       family_score=rep['representative_score'],family_dm=rep['representative_dm'],
                       family_source_curve=rep['source_curve'],family_source_bin=rep['source_bin'],
                       family_originally_selected=rep['originally_selected'],
                       is_group_representative=(row['source_curve']==rep['source_curve'] and int(row['window_peak_bin'])==int(rep['source_bin'])))
        else:out.update(family_id='',is_group_representative=False)
        mapped.append(out)
save('target_group_membership.csv',mapped)
summary=[]
for target in ['C','H','L','M','E']:
    for h in [1,2,3]:
        candidates=[r for r in mapped if r['target']==target and int(r['harmonic'])==h and int(r['Nt'])==16 and r['selection_scope']=='within_0.3_of_catalogue_DM']
        if candidates:summary.append(candidates[0])
save('Nt16_target_inspection.csv',summary)
lines=['# Practical lower-threshold inspection', '',
       'All 16,464 original curves were audited. The original 4,397 peaks and 415 groups were reproduced exactly. No searches, track recoveries or folds were rerun.', '',
       'The following round cuts apply only to the 16-segment grid, leaving the six other timescales at their original cuts. Grouping is the historical greedy 0.035-Hz grouping. These are workload counts, not calibrated detection significances or counts of distinct sources.', '',
       '| Nt16 score cut | Raw peaks at Nt16 | Groups at Nt16 | Additional groups | Total with other timescales unchanged |',
       '|---:|---:|---:|---:|---:|']
for r in counts:
    lines.append(f"| {r['score_cut']:.3f} | {r['raw_peaks_Nt16']:,} | {r['groups_Nt16']:,} | {r['additional_groups_Nt16']:,} | {r['total_groups_if_other_timescales_unchanged']:,} |")
lines += ['', '## Target probes at 16 segments', '',
          'DMs below are restricted to within 0.3 of the reference catalogue DM. Each value is the highest interior local maximum in the declared frequency window. Harmonic number 1 means the fundamental. A nearby peak is not a confirmed pulsar track, particularly for the wider binary windows of H and E.', '',
          '| Target | Harmonic | DM | Peak Hz | Score | Group representative? | Representative Hz |',
          '|---|---:|---:|---:|---:|---|---:|']
for r in summary:
    lines.append(f"| {r['target']} | {r['harmonic']} | {float(r['dm']):.2f} | {float(r['window_peak_hz']):.9f} | {float(r['window_peak_score']):.3f} | {r['is_group_representative']} | {float(r.get('family_hz') or 0):.9f} |")
lines += ['', 'A cut of 58 at Nt16 admits the recorded peaks near C at its fundamental and second harmonic, E in its fundamental Doppler window, and L at its third harmonic. The source curves and exact endpoint bins are in Nt16_target_inspection.csv. This is a post-hoc exploratory operating point applied uniformly across the band, not a calibrated new blind-detection threshold.', '',
          'The peaks near the fundamentals of L and M are much weaker. H has a broader Doppler window and the highest score in that window need not be its actual frequency. They need track/dispersion inspection before any recovery claim. A useful candidate review should keep member alternatives when a stronger nearby endpoint becomes the group representative.', '',
          'All other harmonics, DMs, segment counts and individual source curves remain in the core audit tables. Original inputs and the original candidate selection were preserved.', '']
with (W/'INSPECTION_SUMMARY.md').open('x') as out:out.write('\n'.join(lines))
result=dict(wall_s=time.perf_counter()-start,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            job_id=os.environ.get('SLURM_JOB_ID'),cuts=counts)
with (W/'supplement_complete.json').open('x') as out:json.dump(result,out,indent=2);out.write('\n')
print(json.dumps(result),flush=True)
