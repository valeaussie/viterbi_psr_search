"""Export the exact phase convention used by the installed PRESTO prepfold."""
import argparse,csv,hashlib,json,sys
from pathlib import Path
import presto.presto as p

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--production',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args()
 root=a.production/'paired_19348703/presto'
 with (root/'classified_all.csv').open(newline='') as f:rows=list(csv.DictReader(f))
 with (root/'rfi_flags.csv').open(newline='') as f:flags={r['candidate_id']:r for r in csv.DictReader(f)}
 output=[]
 for r in rows:
  source=Path(r['source_table']);binary=Path(str(source)+'.cand');c=p.fourierprops()
  p.get_rzw_cand(str(binary),int(r['presto_candidate_number']),c)
  info=json.loads((source.parent/'section.section.json').read_text());T=info['fft_duration_s']
  z0=float(c.z)-.5*float(c.w);r0=float(c.r)-.5*z0-float(c.w)/6
  if abs(float(c.r)/T-float(r['peak_freq_hz']))>.002:raise ValueError('Binary/text candidate mismatch '+r['candidate_id'])
  output.append(dict(r,rfi_flag=flags[r['candidate_id']]['rfi_flag'],binary_r=float(c.r),binary_z=float(c.z),binary_w=float(c.w),
    f_start_hz=r0/T,fdot_start_hz_s=z0/(T*T),fddot_hz_s2=float(c.w)/(T*T*T),
    phase_reference='start of the searched section in barycentric samples',section_metadata=info,
    source_binary=str(binary),source_binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest()))
 if len(output)!=5538 or len({r['candidate_id'] for r in output})!=5538:raise ValueError('Candidate count/identity mismatch')
 a.out.write_text(json.dumps(dict(presto_root=str(root),candidates=output,conversion='z0=z-w/2; r0=r-z0/2-w/6; f=r0/T; fdot=z0/T^2; fddot=w/T^3, as installed prepfold.c',source_path='/software/presto/src/prepfold.c',retained_candidates=len(output)),indent=2)+'\n')
 print('Exported exact binary phase parameters for all candidates:',len(output))
if __name__=='__main__':main()
