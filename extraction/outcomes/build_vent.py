from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from .vent_engine import worker
from extraction.common import DATA_DIR, OUTPUT_DIR
import os

O=DATA_DIR/'outcomes'
V=O/'raw'
QUALITY=OUTPUT_DIR/'outcomes'
O.mkdir(parents=True,exist_ok=True)
QUALITY.mkdir(parents=True,exist_ok=True)

def main():
    raw=pd.read_csv(DATA_DIR/'analysis_input.csv',low_memory=False)
    raw['t0']=raw.echo_time
    raw['t0_risk_eligible']=True
    inputs={}
    for key,name in [('raw','airway_chartevents_raw.csv'),('settings','ventilator_setting_raw.csv'),('derived','ventilation_intervals_raw.csv'),('procs','airway_procedureevents_raw.csv'),('stays','all_icu_stays_30d_raw.csv'),('adms','all_admissions_30d_raw.csv')]:
        d=pd.read_csv(V/name,low_memory=False,keep_default_na=False,na_values=[''])
        inputs[key]={int(k):g.to_dict('records') for k,g in d.groupby('cohort_subject_id')}
    jobs=[(r,*[inputs[k].get(int(r['subject_id']),[]) for k in ['raw','settings','derived','procs','stays','adms']]) for r in raw.to_dict('records')]
    with ProcessPoolExecutor(max_workers=int(os.environ.get('MIMIC_WORKERS','12'))) as pool:
        results=list(pool.map(worker,jobs,chunksize=8))
    d=pd.DataFrame(r[0] for r in results)
    assert len(d)==len(raw) and d.subject_id.is_unique
    d.to_csv(O/'vent_adjudication.csv',index=False)
    for days in [7,28]:
        assert d[f'vfd{days}'].dropna().between(0,days).all()
        d[f'vfd{days}_class'].value_counts().rename_axis('reason').reset_index(name='n').to_csv(QUALITY/f'vfd{days}_counts.csv',index=False)
    pd.DataFrame(z for r in results for z in r[1]).to_csv(O/'vent_state_intervals.csv',index=False)
    print(d[['vfd7','vfd28','death7']].notna().sum().to_string())

if __name__=='__main__': main()
