"""Apply J's measured frequency-correction step, retaining all O orbital terms."""
from pathlib import Path
from decimal import Decimal
import copy,json,math,re
from presto.prepfold import pfd
import full_orbit as f
W=Path(__file__).resolve().parents[1];S=W/'spin_refined'
m=json.loads((W/'manifest.json').read_text());record_path=W/'folds/frequency_refinement/complete.json'
r=json.loads(record_path.read_text());source=next((W/'folds/frequency_refinement').glob('*.pfd'))
obj=pfd(str(source));offsets=[float(v) for v in obj.freq_offsets()]
assert abs(offsets[1])<1e-18 and abs(offsets[2])<1e-23,offsets
mapping=f.helper().read_polycos(W/'mapping_topo/polyco_for_observation.dat')
epoch=Decimal(str(m['filterbank_metadata']['header']['tstart']))
_,rate,_=f.evaluate(mapping,epoch);rate/=m['mapping_clock_hz'];assert .999<rate<1.001
df=offsets[0]/rate
new=copy.deepcopy(m);new['orbit']['spin_hz']+=df
for name in ['amplitude_hz','period_s','zero_crossing_s','reference_s']:assert new['orbit'][name]==m['orbit'][name]
bestdm=float(obj.bestdm);new['fold_dm']='%.17g'%bestdm
new['refinement']=dict(source_complete=str(record_path),source_complete_sha256=f.sha(record_path),
    source_pfd=str(source),source_pfd_sha256=f.sha(source),initial_orbit=m['orbit'],
    presto_initial_topo_frequency_hz=float(obj.fold_p1),presto_refined_topo_period_ms=float(obj.topo_p1)*1000,
    topocentric_frequency_adjustment_hz=offsets[0],barycentric_time_rate=rate,barycentric_spin_adjustment_hz=df,
    fold_dm=bestdm,orbital_parameters_unchanged=True,published_ephemeris_used=False,
    formal_spin_error_units=abs(df)/m['orbit']['spin_formal_error_hz'] if m['orbit']['spin_formal_error_hz']>0 else None)
for name in ['folds','logs']:(S/name).mkdir(exist_ok=True)
f.save(S/'manifest.json',new)
print(json.dumps(new['refinement'],indent=2),flush=True)
