"""Complete reference-column output, keeping all original Table 3 entries."""
from pathlib import Path
import json
import report
from paper_table_reference import TABLE

W=Path(__file__).resolve().parents[1]
def main():
    report.main()
    result=json.loads((W/'results.json').read_text())
    by={r['target']:r for r in result['ranked_reference_folds']}
    tex=['% Original Table 3 entries retained. New reference folds use the original filterbank.',
        '% Nominal sigma maxima are provisional until the native fold diagnostics are reviewed.',
        r'\begin{tabular}{lrrrrrrl}',r'\hline',
        r'Pulsar & $S_{1400}$ & $w_{50}$ & Binary & $L-L_{\rm th}$ & $Z$ & $\sigma_{\rm ref}$ & Outcome \\',
        r' & (mJy) & (per cent) & & (blind) & (targeted) & (reference fold) & \\',r'\hline']
    rows=[]
    for target,flux,duty,binary,like,z,outcome in TABLE:
        item=by[target];r=item['strongest_completed'];sigma=None if r is None else r['nominal_sigma']
        tex.append(' & '.join(['47 Tuc '+target,flux,duty,binary,like,z,'--' if sigma is None else '%.1f'%sigma,outcome])+r' \\')
        rows.append(dict(target=target,nominal_sigma=sigma,ranked_model=r['model'] if r else None,
            pdf=r['pdf'] if r else None,model_file=r['ephemeris'] if r else None,
            all_models_finished=item['all_models_finished'],review_status='Pending visual review; no revised detection claim.'))
    tex.extend([r'\hline',r'\end{tabular}'])
    (W/'table3_with_original_data_reference_column.tex').write_text('\n'.join(tex)+'\n')
    (W/'table3_reference_column.json').write_text(json.dumps(rows,indent=2)+'\n')
    note=(r'The column $\sigma_{\rm ref}$ gives the nominal PRESTO profile significance of the strongest credible reference fold among the independently supplied timing models tested on the original filterbank.'+'\n'
        r'The tested models include the timing files linked on Paulo Freire\textquotesingle s 47 Tuc website and previously saved reference ephemerides, including the observation-epoch model used for 47 Tuc O in Figure~3.'+'\n'
        r'Each model is folded using all supplied timing parameters, followed by the same local PRESTO search over dispersion measure, period and its first two derivatives.'+'\n'
        r'The reference-fold significances are distinct from the HMM likelihood excess and targeted statistic $Z$; maximisation over models and refinement trials does not give a calibrated detection significance.'+'\n')
    (W/'table3_reference_column_note.tex').write_text(note)
    (W/'TABLE_REVIEW_REQUIRED.md').write_text('The table contains the largest nominal sigma among completed alternatives. Before copying the column note or final values into the paper, review the time, frequency, DM and pulse-profile panels and confirm that the selected fold is credible. All candidate alternatives, including any failed models, remain in results.json. Existing Outcome entries have not been changed.\n')
if __name__=='__main__':main()
