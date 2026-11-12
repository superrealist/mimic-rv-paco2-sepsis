"""Verify arterial samples and select the first eligible recorded echo pair.

The internal t0 is the later paired measurement for eligibility at this stage.
finalize_cohort.py preserves it as pair_end and sets final T0 to echo_time.
"""
from .common import ROOT
import pandas as pd,numpy as np,json,hashlib

NORMAL={'Nl RV function','Low normal function','Hyperdynamic'}
DYS={'Mild global RV hypo','Moderate global RV hypo','Severe global hypo','RV function depressed'}
ITEMS={50818:'pco2',50820:'ph',50821:'po2',50816:'fio2_sample',50813:'index_lactate',50803:'bicarbonate'}

def read(name,ts=()):
    d=pd.read_csv(ROOT/'data'/name)
    for c in ts:d[c]=pd.to_datetime(d[c],format='mixed')
    return d

def main():
    raw=read('raw_gas_labevents.csv',['charttime','storetime'])
    raw.groupby(['itemid','valueuom'],dropna=False).size().rename('rows').to_csv(ROOT/'tables'/'raw_gas_units.csv')
    samples=[];amb=[]
    for key,g in raw.groupby(['subject_id','hadm_id','specimen_id'],dropna=False):
        spec=set(g.loc[g.itemid.eq(52033),'value'].dropna().astype(str).str.upper().str.replace('.','',regex=False).str.strip())
        if spec!={'ART'}:continue
        chart=g.charttime.dropna().unique()
        if len(chart)!=1:
            amb.append(dict(subject_id=key[0],hadm_id=key[1],specimen_id=key[2],reason='specimen_charttime_conflict'));continue
        d=dict(subject_id=key[0],hadm_id=key[1],specimen_id=key[2],bg_time=pd.Timestamp(chart[0]),bg_storetime=g.loc[g.itemid.eq(50818),'storetime'].max())
        conflict=False
        for item,col in ITEMS.items():
            z=g.loc[g.itemid.eq(item),'valuenum'].dropna().unique()
            if len(z)>1:
                if col=='pco2':conflict=True
                d[col]=np.nan
            else:d[col]=z[0] if len(z) else np.nan
        pcunits=g.loc[g.itemid.eq(50818),'valueuom'].dropna().astype(str).str.lower().str.replace(' ','',regex=False)
        if len(pcunits) and not pcunits.isin(['mmhg']).all():conflict=True
        if conflict:amb.append(dict(subject_id=key[0],hadm_id=key[1],specimen_id=key[2],reason='within_specimen_core_or_unit_conflict'));continue
        if pd.isna(d['pco2']) or d['pco2']<=0:continue
        if d['pco2']>300:
            amb.append(dict(subject_id=key[0],hadm_id=key[1],specimen_id=key[2],reason='PaCO2_above300_requires_source_review'));continue
        samples.append(d)
    gas=pd.DataFrame(samples)
    # PaCO2 conflicts prevent an exposure definition. Accessory conflicts only
    # invalidate that accessory variable; they do not remove reliable PaCO2.
    keys=['subject_id','hadm_id','bg_time']
    nc=gas.groupby(keys)[['pco2','ph','po2']].nunique()
    (nc.reset_index()).to_csv(ROOT/'logs'/'raw_arterial_timestamp_value_counts.csv',index=False)
    conflict_keys=nc[nc.pco2.gt(1)].reset_index()[keys]
    conflict_keys.to_csv(ROOT/'logs'/'raw_arterial_timestamp_conflicts.csv',index=False)
    gas=gas.merge(conflict_keys.assign(ambiguous=True),on=keys,how='left')
    gas=gas[gas.ambiguous.isna()].copy()
    for col in ITEMS.values():
        counts=gas.groupby(keys)[col].transform('nunique')
        gas.loc[counts.gt(1),col]=np.nan
    gas=gas.sort_values('specimen_id').drop_duplicates(keys)
    gas.to_csv(ROOT/'data'/'verified_arterial_samples.csv',index=False)
    pd.DataFrame(amb).to_csv(ROOT/'logs'/'raw_arterial_specimen_exceptions.csv',index=False)
    a=read('echo_anchors.csv',['intime','outtime','echo_time','vent_start','vent_end'])
    rv=read('rv_raw.csv',['measurement_datetime'])
    r=rv[rv.measurement.eq('rv_function')].copy();r['result']=r.result.fillna('').str.strip()
    rr=r.groupby('measurement_id').result.agg(lambda z:'|'.join(sorted(set(z))))
    a=a.merge(rr.rename('rv_function_raw'),on='measurement_id',how='left')
    a['rv_binary']=a.rv_function_raw.map({**{s:0 for s in NORMAL},**{s:1 for s in DYS}})
    a['echo_midnight']=a.echo_time.dt.time.eq(pd.Timestamp('2000-01-01').time())
    a['echo_time_precision']='recorded_diagnostic_datetime'
    a.loc[a.echo_midnight,'echo_time_precision']='midnight_placeholder_status_undocumented'
    tm=read('time_sepsis_precision_by_stay.csv',['culture_lower','culture_upper','sepsis_complete_upper'])
    a=a.merge(tm[['stay_id','culture_lower','culture_upper','precision','unresolved_link','sepsis_complete_upper']],on='stay_id',how='left',validate='many_to_one')
    gby={k:v for k,v in gas.groupby('hadm_id')}
    pairs=[]
    for x in a[a.rv_binary.notna()&a.icu_matches.eq(1)&a.sepsis_complete_upper.notna()&~a.unresolved_link.fillna(True)].itertuples(index=False):
        g=gby.get(x.hadm_id)
        if g is None:continue
        g=g[(g.subject_id==x.subject_id)&g.bg_time.ge(max(x.intime,x.vent_start))&g.bg_time.lt(min(x.outtime,x.vent_end))]
        if g.empty:continue
        g=g.assign(delta_hours=(g.bg_time-x.echo_time).dt.total_seconds()/3600)
        g=g[g.delta_hours.abs().le(12)]
        if g.empty:continue
        # Qualify the pair before ranking: sepsis may become ascertainable within the window.
        g=g[g.bg_time.map(lambda t:max(t,x.echo_time)>=x.sepsis_complete_upper)]
        if g.empty:continue
        g=g.assign(absdelta=g.delta_hours.abs(),after=g.delta_hours.gt(0)).sort_values(['absdelta','after','bg_time','specimen_id'])
        for h in [6,12]:
            gh=g[g.absdelta.le(h)]
            if gh.empty:continue
            b=gh.iloc[0];d=x._asdict();d.update({c:b[c] for c in ['specimen_id','bg_time','bg_storetime','pco2','ph','po2','fio2_sample','index_lactate','bicarbonate','delta_hours']})
            d.update(window_hours=h,t0=max(x.echo_time,b.bg_time),prebase=min(x.echo_time,b.bg_time))
            pairs.append(d)
    pairs=pd.DataFrame(pairs)
    pairs.to_csv(ROOT/'data'/'all_qualified_pairs_12h.csv',index=False)
    audit=[]
    for h in [6,12]:
        p=pairs[pairs.window_hours.eq(h)].sort_values(['echo_time','measurement_id','bg_time'])
        for precise in [False,True]:
            q=p[~p.echo_midnight] if precise else p
            q=q.drop_duplicates('subject_id')
            label=f'cohort_{h}h_'+('time_precise' if precise else 'recorded_time_with_midnight')
            q.to_csv(ROOT/'data'/f'{label}.csv',index=False)
            audit.append(dict(cohort=label,n=len(q),rv_not_reduced=int(q.rv_binary.eq(0).sum()),rv_reduced=int(q.rv_binary.eq(1).sum()),tte_before_bg=int(q.delta_hours.gt(0).sum()),midnight=int(q.echo_midnight.sum()),date_only_culture=int(q.precision.astype(str).str.contains('date',case=False).sum()),median_abs_gap_h=q.delta_hours.abs().median()))
    q=pd.read_csv(ROOT/'data'/'cohort_12h_recorded_time_with_midnight.csv')
    q['age_90plus_proxy']=q.age_derived.ge(90).astype(int)
    q['anchor_age_censored']=q.anchor_age.gt(89).astype(int)
    q['age']=q.age_derived.clip(upper=90)
    q['temporal_B']=q.delta_hours.gt(0)
    q.to_csv(ROOT/'data'/'cohort_primary_pre_outcome_validation.csv',index=False)
    pd.DataFrame(audit).to_csv(ROOT/'tables'/'cohort_time_window_preflight.csv',index=False)
    freeze={'stage':'preflight_pairing_no_outcome_associations','primary_hours':12,'sensitivity_hours':6,'midnight_main':'retain_source_recorded_datetime','sepsis':'source_consistent_upper_bound_at_or_before_later_paired_measurement','age90':'min(anchor_age + admission_year - anchor_year, 90), deidentified proxy, censoring descriptors not extra main predictors','gas_source':'explicit arterial raw specimen, contradictory PaCO2 excluded; accessory conflicts only accessory NA; PCO2-specific storetime','main_n_before_mortality_validation':len(q),'sha256':hashlib.sha256((ROOT/'data'/'cohort_primary_pre_outcome_validation.csv').read_bytes()).hexdigest()}
    (ROOT/'logs'/'cohort_preflight_freeze.json').write_text(json.dumps(freeze,indent=2),encoding='utf-8')
    (ROOT/'logs'/'cohort_builder_complete.json').write_text(json.dumps(freeze,indent=2),encoding='utf-8')
    print(pd.DataFrame(audit).to_string(index=False))

if __name__=='__main__':main()
