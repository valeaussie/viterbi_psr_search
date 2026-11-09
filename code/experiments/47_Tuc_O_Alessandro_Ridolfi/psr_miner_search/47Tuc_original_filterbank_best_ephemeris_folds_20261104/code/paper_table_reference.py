"""Build fixed-ephemeris fold summaries and a standalone extension of supplied Table 3."""
from pathlib import Path
from collections import Counter
import csv
import io
import json
from datetime import datetime, timezone

W=Path(__file__).resolve().parents[1]
# Supplied screenshot: preserve its row order and values, adding only sigma_eph.
TABLE=[
 ('J','0.54(6)','15','yes','737.5','','blind recovery'),
 ('C','0.36(4)','13','no','','22.7','targeted, detected'),
 ('D','0.22(3)','10','no','1099.1','','blind recovery'),
 ('E','0.21(3)','25','yes','','12.4','targeted, detected'),
 ('F','0.15(2)','14','no','192.9','','blind recovery'),
 ('O','0.10(1)','12','yes','82.9','','blind recovery'),
 ('H','0.09(2)','13','yes','','5.2','targeted, ambiguous'),
 ('I','0.09(1)','14','yes','6.8','','blind recovery'),
 ('M','0.07(2)','25','no','','6.1','targeted, ambiguous'),
 ('U','0.06(1)','8','yes','','','not recovered'),
 ('G','0.05(2)','13','no','','','not recovered'),
 ('Q','0.05(2)','11','yes','','','not recovered'),
 ('L','0.04(1)','15','no','','','targeted, probable'),
 ('N','0.03(1)','11','no','','','not recovered')]

def main():
    m=json.loads((W/'manifest.json').read_text());rows=[]
    for r in m['rows']:
        original_out=Path(r['output']);out=original_out;attempt=1
        while (out/'failure.json').exists():
            attempt+=1;candidate=original_out/('attempt_%02d'%attempt)
            if not candidate.exists():break
            out=candidate
        marker=out/'complete.json';failed=out/'failure.json'
        rec=json.loads(marker.read_text()) if marker.exists() else {}
        rows.append(dict(pulsar='47 Tuc '+r['target'],target=r['target'],
            status='complete' if marker.exists() else 'failed' if failed.exists() else 'pending_or_running',
            native_sigma=rec.get('sigma'),reduced_chi2=rec.get('reduced_chi2'),dm=r['dm'],
            binary_model=r['binary'],ephemeris_url=r['source_url'],par_sha256=r['par_sha256'],
            fitted_finish_mjd=r['fitted_finish_mjd'],extrapolated=r['extrapolated'],
            days_after_fitted_finish=r['days_after_fitted_finish'],source_warning=r['source_warning'],
            wall_s=rec.get('wall_s'),pdf=rec.get('pdf',''),png=rec.get('png',''),
            failure_record=str(failed) if failed.exists() else '',
            failure_message=json.loads(failed.read_text()).get('traceback','')[-1200:] if failed.exists() else ''))
    for r in m['unavailable']:
        rows.append(dict(pulsar='47 Tuc '+r['target'],target=r['target'],status='no_public_par_file',source_warning=r['reason']))
    fields=list(dict.fromkeys(k for row in rows for k in row))
    s=io.StringIO();writer=csv.DictWriter(s,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    (W/'published_ephemeris_folds.csv').write_text(s.getvalue())
    summary=dict(updated_utc=datetime.now(timezone.utc).isoformat(),counts=dict(Counter(r['status'] for r in rows)),
        available_ephemerides=len(m['rows']),known_pulsars_recorded=len(rows),rows=rows)
    (W/'results.json').write_text(json.dumps(summary,indent=2)+'\n')
    by={r['target']:r for r in rows}
    tex=['% Standalone extension of the user-supplied Table 3. Existing values transcribed unchanged.',
        '% Add this tabular inside the existing table environment. No citation key is invented here.',
        r'\begin{tabular}{lrrrrrrl}',r'\hline',
        r'Pulsar & $S_{1400}$ & $w_{50}$ & Binary & $L-L_{\rm th}$ & $Z$ & $\sigma_{\rm eph}$ & Outcome \\',
        r' & (mJy) & (per cent) & & (blind) & (targeted) & (published ephemeris) & \\',r'\hline']
    column=[]
    for target,flux,duty,binary,like,z,outcome in TABLE:
        row=by[target];sigma=row.get('native_sigma')
        value=('--' if sigma is None else '%.1f'%sigma)
        tex.append(' & '.join(['47 Tuc '+target,flux,duty,binary,like,z,value,outcome])+r' \\')
        column.append(dict(pulsar='47 Tuc '+target,native_sigma=sigma,status=row['status'],pdf=row.get('pdf',''),png=row.get('png','')))
    tex += [r'\hline',r'\end{tabular}']
    (W/'table3_with_published_ephemeris_column.tex').write_text('\n'.join(tex)+'\n')
    (W/'table3_ephemeris_column.json').write_text(json.dumps(column,indent=2)+'\n')
    # Each prose sentence is kept on one source line for Overleaf copying.
    note=(r'The column $\sigma_{\rm eph}$ gives the nominal PRESTO profile significance obtained by folding the same cleaned observation with the latest timing file linked on Paulo Freire\textquotesingle s 47 Tuc website at the time of retrieval (2026 November 4), with all timing parameters and the dispersion measure held fixed.'+'\n'
          r'All downloaded timing solutions end before the observation epoch, so these folds extrapolate the published models; a weak fold does not establish that the pulsar was absent.'+'\n'
          r'These profile significances are distinct from the blind-search likelihood excess and the targeted statistic $Z$.'+'\n')
    (W/'table3_column_note.tex').write_text(note)
    text=['# 47 Tuc published ephemerides: fixed-model PRESTO folds','',
        'Source: https://www3.mpifr-bonn.mpg.de/staff/pfreire/47Tuc/index.html',
        'Retrieved 2026 November 4. The latest downloadable parameter-file column contains 23 timing files; 19 further pulsars have no linked timing file. No ATNF file is used. Original bytes and HTTP headers are archived.',
        '',
        'The input is the same four-hour RFI-cleaned, 128-channel filterbank used in the search comparison. Every published spin, astrometric and binary parameter is retained, with DE421 and UTC(NIST) as supplied. Only a leading J in the internal PSR name is removed for TEMPO TZ compatibility. No time, radio channel or candidate is selected on the basis of the fold result.',
        '',
        'TEMPO produces MeerKAT topocentric phase predictors from the full models. Independently generated 12-coefficient and 15-coefficient interpolants on 30-minute spans are compared across 1001 epochs before folding. This checks numerical interpolation and full coverage, not the physical accuracy of extrapolating a model beyond its fitted data span. The initial pilots stopped because TEMPO rejects the requested shorter spans; a second attempt exposed a PRESTO limit on input filename length. All attempts and logs remain preserved; the corrected pilot folds use attempt_03.',
        '',
        'PRESTO uses 64 time parts, 32 frequency subbands and its native default profile length. The published DM and timing model are held fixed with -nosearch. Native PFD, best-profile text, PostScript, PDF, PNG, polycos, logs and checksums are retained in each pulsar folder. Native sigma values are nominal profile statistics, not calibrated probabilities of discovery.',
        '',
        'All 23 timing spans end before this observation (MJD %.9f). Orbital extrapolation can be particularly unreliable for variable-period binaries. The aa file also retains the older 541.89365-Hz convention while the website says the physical spin period is twice the original value. Its result must be described as folding the supplied harmonic model.'%m['observation_start_topocentric_mjd'],
        '',
        'New column: sigma_eph, fixed published-ephemeris PRESTO profile significance. The supplied Table 3 row order and existing numbers are transcribed without changing their interpretation. The standalone LaTeX table, short note and all-42-pulsar CSV are regenerated from completion records.',
        '',
        'Current status: '+json.dumps(summary['counts']),
        '',
        '| Pulsar | Status | Native sigma | Published DM | Model | Native fold PDF |',
        '|---|---|---:|---:|---|---|']
    for r in rows:
        sig='--' if r.get('native_sigma') is None else '%.1f'%r['native_sigma']
        link='[PDF]('+r['pdf']+')' if r.get('pdf') else '--'
        text.append('| '+r['pulsar']+' | '+r['status']+' | '+sig+' | '+str(r.get('dm','--'))+' | '+r.get('binary_model','--')+' | '+link+' |')
    (W/'README.md').write_text('\n'.join(text)+'\n')
    print(json.dumps(dict(counts=summary['counts'],updated_utc=summary['updated_utc'],table_column=column),indent=2))

if __name__=='__main__':main()
