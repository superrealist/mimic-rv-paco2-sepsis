"""Assemble V2 source-aware preflight sets; no outcome association analyses."""
from .common import ROOT
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,datetime

T=lambda x:pd.to_datetime(x,format='mixed',errors='coerce')

def load(name,times=()):
    d=pd.read_csv(ROOT/'data'/name)
    for k in times:
        if k in d:d[k]=T(d[k])
    return d

def grouped(d,key='cohort_subject_id'):
    return {int(k):g.copy() for k,g in d.groupby(key)}

def last(g,col,time='charttime',valid=None):
    if g is None or g.empty or col not in g:return np.nan,pd.NaT,False
    z=g[g[col].notna()].copy()
    if valid is not None:z=z[valid(z[col])]
    if z.empty:return np.nan,pd.NaT,False
    when=z[time].max();z=z[z[time].eq(when)]
    conflict=z[col].nunique()>1
    return (np.nan if conflict else z[col].iloc[0]),when,bool(conflict)

def add_measure(out,g,col,new,valid,time='charttime',zero_valid=False):
    raw,rt,rc=last(g,col,time)
    value,vt,vc=last(g,col,time,valid)
    out[new]=value;out[new+'_time']=vt
    out[new+'_raw_latest']=raw;out[new+'_raw_latest_time']=rt
    out[new+'_raw_latest_tie_conflict']=rc;out[new+'_valid_tie_conflict']=vc
    if rc:status='same_time_value_conflict'
    elif pd.isna(raw):status='no_nonmissing_record'
    elif raw==0:status='observed_zero_valid' if zero_valid else 'observed_zero_not_valid_for_analysis'
    elif bool(valid(pd.Series([raw])).iloc[0]):status='valid_observed'
    else:status='observed_out_of_range'
    out[new+'_latest_raw_status']=status
    out[new+'_invalid_or_zero_after_selected']=bool(pd.notna(vt) and pd.notna(rt) and rt>vt and status not in ['valid_observed','observed_zero_valid'])

PRESSORS={221906:'norepinephrine',221289:'epinephrine',221749:'phenylephrine',221662:'dopamine',222315:'vasopressin',229709:'angiotensin_ii',229764:'angiotensin_ii'}
INOTROPES={221653:'dobutamine',221986:'milrinone'}
MASS_RATE={'mcg/kg/min','mcg/min','mcg/kg/hour','mcg/hour','mg/kg/min','mg/min','mg/kg/hour','mg/hour','ng/kg/min','ng/min','ng/kg/hour','ng/hour'}

def fio2_normalize(z):
    z=pd.to_numeric(z,errors='coerce')
    return z.where(~z.between(.21,1.0),z*100).where(lambda v:v.between(21,100))

def main():
    inp=ROOT/'data/cohort_primary_pre_outcome_validation.csv'
    marker=ROOT/'logs/baseline_preflight_extraction_complete.json'
    if not marker.exists():raise RuntimeError('Do not read incomplete source exports.')
    provenance=json.loads(marker.read_text(encoding='utf-8'))
    if provenance['input_sha256']!=hashlib.sha256(inp.read_bytes()).hexdigest():raise RuntimeError('Core input changed after extraction.')
    final_marker=ROOT/'logs/baseline_preflight_complete.json'
    if final_marker.exists():final_marker.unlink()
    q=load(inp.name,['intime','outtime','echo_time','bg_time','bg_storetime','vent_start','vent_end','prebase','t0'])
    assert q.subject_id.is_unique and q.stay_id.is_unique
    sources={}
    for key,ts in [('sofa',['starttime','endtime']),('vitals',['charttime']),('vent',['charttime']),('nee',['starttime','endtime']),('lactate',['charttime']),('chemistry',['charttime']),('height',['charttime']),('weight',['charttime','storetime'])]:
        sources[key]=grouped(load('preflight_'+key+'_raw.csv',ts))
    smarker=ROOT/'logs/sofa_component_availability_complete.json'
    if not smarker.exists():raise RuntimeError('Run python -m extraction.extract_baseline --sofa-availability-only before source-aware SOFA flags.')
    if json.loads(smarker.read_text(encoding='utf-8'))['input_sha256']!=provenance['input_sha256']:raise RuntimeError('SOFA source audit input mismatch.')
    sofa_sources=grouped(load('preflight_sofa_component_sources.csv',['sofa_end','charttime']),'subject_id')
    press=load('preflight_pressors_raw.csv',['starttime','endtime','storetime'])
    press=press.merge(q[['subject_id','prebase']],left_on='cohort_subject_id',right_on='subject_id',suffixes=('','_cohort'),validate='many_to_one')
    press['drug']=press.itemid.map({**PRESSORS,**INOTROPES})
    press['unit_norm']=press.rateuom.fillna('').str.strip().str.lower()
    press['status_norm']=press.statusdescription.fillna('').str.strip().str.lower()
    press['known_rate_unit']=press.unit_norm.isin(MASS_RATE)
    press.loc[press.drug.eq('vasopressin'),'known_rate_unit']=press.loc[press.drug.eq('vasopressin'),'unit_norm'].isin(['units/min','units/hour'])
    press['usable_status']=~press.status_norm.isin(['rewritten','cancelled','canceled','bolus'])
    press['continuous_interval']=press.endtime.gt(press.starttime)&~press.ordercategorydescription.fillna('').str.contains('bolus|drugpush',case=False,regex=True)
    press['positive_valid_rate']=press.rate.gt(0)&np.isfinite(press.rate)&press.known_rate_unit
    press['active_left']=press.starttime.lt(press.prebase)&press.endtime.ge(press.prebase)&press.usable_status&press.continuous_interval
    press['active_halfopen']=press.starttime.le(press.prebase)&press.endtime.gt(press.prebase)&press.usable_status&press.continuous_interval
    press['support_left']=press.active_left&press.positive_valid_rate
    press['support_halfopen']=press.active_halfopen&press.positive_valid_rate
    press.to_csv(ROOT/'data/preflight_pressors_audited.csv',index=False)
    pg=grouped(press)
    rows=[]
    for x in q.itertuples(index=False):
        sid=int(x.subject_id);o={'subject_id':sid}
        g=sources['sofa'].get(sid)
        for k in ['sofa','noncardiopulmonary_sofa','ncp_sofa_observed_all4','ncp_sofa_source_aware','sofa_zero_assumed_n_organs','sofa_window_rows']:o[k]=np.nan
        o['sofa_time_used']=pd.NaT
        for organ in ['coagulation','liver','cns','renal']:
            o[organ+'_sofa_24h']=np.nan;o[organ+'_source_observed']=False
        if g is not None and not g.empty:
            g=g[g.endtime.le(x.prebase)&g.endtime.ge(x.prebase-pd.Timedelta(hours=24))]
            if len(g):
                z=g.sort_values(['endtime','hr']).iloc[-1]
                window=g[g.hr.between(z.hr-23,z.hr)&g.endtime.le(z.endtime)]
                vals=[];obs=[]
                for organ in ['coagulation','liver','cns','renal']:
                    sg=sofa_sources.get(sid)
                    value=z[organ+'_24hours'];seen=bool(sg is not None and sg.component.eq(organ).any())
                    o[organ+'_sofa_24h']=value;o[organ+'_source_observed']=seen
                    vals.append(value);obs.append(seen)
                total=sum(vals) if all(pd.notna(v) for v in vals) else np.nan
                o.update(sofa=z.sofa_24hours,noncardiopulmonary_sofa=total,sofa_time_used=z.endtime,ncp_sofa_observed_all4=all(obs),ncp_sofa_source_aware=total if all(obs) else np.nan,sofa_zero_assumed_n_organs=sum(not seen and v==0 for seen,v in zip(obs,vals)),sofa_window_rows=len(window),sofa_selected_hr=z.hr,sofa_latest_gap_hours=(x.prebase-z.endtime).total_seconds()/3600,sofa_source_window_start=z.endtime-pd.Timedelta(hours=24),sofa_source_window_end=z.endtime)
        g=sources['vitals'].get(sid)
        if g is not None:
            g=g[g.charttime.lt(x.prebase)&g.charttime.ge(x.prebase-pd.Timedelta(hours=2))].copy()
            g['map_valid_source']=g.mbp.where(g.mbp.between(20,200),g.mbp_ni.where(g.mbp_ni.between(20,200)))
            g['map_raw']=g.mbp.fillna(g.mbp_ni)
        add_measure(o,g,'map_valid_source','map',lambda z:z.between(20,200))
        if g is not None:
            o['map_raw_latest'],o['map_raw_latest_time'],o['map_raw_latest_tie_conflict']=last(g,'map_raw')
            rv=o['map_raw_latest'];rt=o['map_raw_latest_time']
            o['map_latest_raw_status']='same_time_value_conflict' if o['map_raw_latest_tie_conflict'] else ('no_nonmissing_record' if pd.isna(rv) else ('valid_observed' if 20<=rv<=200 else 'observed_zero_not_valid_for_analysis' if rv==0 else 'observed_out_of_range'))
            o['map_invalid_or_zero_after_selected']=bool(pd.notna(o['map_time']) and pd.notna(rt) and rt>o['map_time'] and o['map_latest_raw_status']!='valid_observed')
        for src,col,new,cond in [('lactate','lactate','lactate',lambda z:z.gt(0)),('chemistry','creatinine','creatinine',lambda z:z.gt(0))]:
            add_measure(o,sources[src].get(sid),col,new,cond)
        vg=sources['vent'].get(sid)
        pre=None if vg is None else vg[vg.charttime.lt(x.prebase)&vg.charttime.ge(max(x.vent_start,x.prebase-pd.Timedelta(hours=6)))].copy()
        if pre is not None:
            pre['rr']=pre.respiratory_rate_total.fillna(pre.respiratory_rate_set)
            pre['driving_pressure']=(pre.plateau_pressure-pre.peep).where(pre.plateau_pressure.between(5,80)&pre.peep.between(0,40)&pre.plateau_pressure.gt(pre.peep))
        for col,new,cond in [('peep','peep',lambda z:z.between(0,40)),('minute_volume','minute_volume',lambda z:z.gt(0)&z.le(50)),('respiratory_rate_total','rr_total',lambda z:z.between(2,80)),('respiratory_rate_set','rr_set',lambda z:z.between(2,80)),('rr','rr',lambda z:z.between(2,80)),('tidal_volume_observed','tidal_volume',lambda z:z.between(50,2000)),('plateau_pressure','plateau',lambda z:z.between(5,80)),('driving_pressure','driving_pressure',lambda z:z.between(1,70))]:
            add_measure(o,pre,col,new,cond,zero_valid=new=='peep')
        o['vent_hours']=(x.prebase-x.vent_start).total_seconds()/3600
        # Only preceding/current same-episode ventilator settings can supply FiO2.
        bf=None if vg is None else vg[vg.charttime.le(x.bg_time)&vg.charttime.ge(max(x.vent_start,x.bg_time-pd.Timedelta(hours=2)))&vg.charttime.le(x.vent_end)].copy()
        if bf is not None:bf['fio2_percent']=fio2_normalize(bf.fio2)
        bf30=None if bf is None else bf[bf.charttime.ge(x.bg_time-pd.Timedelta(minutes=30))]
        v30,t30,cf30=last(bf30,'fio2_percent',valid=lambda z:z.between(21,100))
        v2,t2,cf2=last(bf,'fio2_percent',valid=lambda z:z.between(21,100))
        sample=fio2_normalize(pd.Series([x.fio2_sample])).iloc[0]
        # Pragmatic main: recorded FiO2 up to2h before ART, same episode.
        # The common-pressure physical approximation is a flag, not a new
        # observed-value correction or automatic primary missingness rule.
        used=sample if pd.notna(sample) else v2
        source='same_arterial_sample' if pd.notna(sample) else ('same_episode_vent_le2h' if pd.notna(v2) else 'unknown')
        po2=x.po2
        valid_po2=bool(pd.notna(po2) and 10<=po2<=800)
        incompat=bool(pd.notna(used) and pd.notna(po2) and po2>used/100*(760-47))
        pf=po2/(used/100) if pd.notna(used) and valid_po2 else np.nan
        used30=sample if pd.notna(sample) else v30
        source30='same_arterial_sample' if pd.notna(sample) else ('same_episode_vent_le30min' if pd.notna(v30) else 'unknown')
        flag30=bool(pd.notna(used30) and pd.notna(po2) and po2>used30/100*(760-47))
        pf30=po2/(used30/100) if pd.notna(used30) and valid_po2 else np.nan
        gap=0 if source=='same_arterial_sample' else ((x.bg_time-t2).total_seconds()/60 if pd.notna(t2) else np.nan)
        gap30=0 if source30=='same_arterial_sample' else ((x.bg_time-t30).total_seconds()/60 if pd.notna(t30) else np.nan)
        o.update(fio2_sample_raw=x.fio2_sample,fio2_sample_percent=sample,fio2_vent30=v30,fio2_vent30_time=t30,fio2_vent30_tie_conflict=cf30,fio2_vent2h=v2,fio2_vent2h_time=t2,fio2_vent2h_tie_conflict=cf2,fio2_used=used,fio2_source=source,fio2_gap_min=gap,pf=pf,pf_pair_incompatible=incompat,pf_missing_reason='paO2_unavailable_or_invalid' if not valid_po2 else ('fio2_unavailable' if pd.isna(used) else ('known_flagged_physical_approximation' if incompat else 'known')))
        o.update(pf_2h_raw=pf,pf_2h_descriptor=pf,pf_2h_qc_sensitivity=pf if not incompat else np.nan,fio2_30m_used=used30,fio2_30m_source=source30,fio2_30m_gap_min=gap30,pf_30m_raw=pf30,pf_30m_sensitivity=pf30,pf_30m_pair_incompatible=flag30,pf_30m_qc_sensitivity=pf30 if not flag30 else np.nan)
        for mins,frame in [(30,bf30),(120,bf)]:
            suffix='30m' if mins==30 else '2h'
            o['fio2_valid_values_n'+suffix]=0 if frame is None else int(frame.fio2_percent.nunique())
            o['fio2_recorded_change_'+suffix]=False if frame is None else bool(frame.fio2_percent.nunique()>1)
            for setting in ['peep','minute_volume','respiratory_rate_total','tidal_volume_observed']:
                o[setting+'_recorded_change_'+suffix]=False if frame is None else bool(frame[setting].dropna().nunique()>1)
        ng=sources['nee'].get(sid)
        o['nee']=np.nan;o['nee_start']=pd.NaT;o['nee_end']=pd.NaT;o['nee_active_intervals']=0;o['nee_conflict']=False;o['nee_imputable']=False
        if ng is not None:
            ng=ng[ng.starttime.lt(x.prebase)&ng.endtime.ge(x.prebase)]
            valid=ng[ng.norepinephrine_equivalent_dose.ge(0)&np.isfinite(ng.norepinephrine_equivalent_dose)]
            o['nee_active_intervals']=len(ng);o['nee_conflict']=valid.norepinephrine_equivalent_dose.nunique()>1
            if len(valid) and not o['nee_conflict']:
                o['nee']=valid.norepinephrine_equivalent_dose.iloc[0];o['nee_start']=valid.starttime.min();o['nee_end']=valid.endtime.max()
        ag=pg.get(sid,press.iloc[0:0]);active=ag[ag.active_left]
        pos=ag[ag.support_left&ag.itemid.isin(PRESSORS)]
        officialpos=pos[pos.itemid.isin([221906,221289,221749,221662,222315])]
        angi=pos[pos.drug.eq('angiotensin_ii')]
        o.update(s_doc=int(len(pos)>0),s_doc_meaning='documented_positive_support' if len(pos) else 'support_not_verified',s_doc_positive_drugs=';'.join(sorted(set(pos.drug))),s_doc_positive_rows=len(pos),s_doc_raw_active_rows=int(active.itemid.isin(PRESSORS).sum()),s_doc_unknown_rate_unit_rows=int((active.itemid.isin(PRESSORS)&(active.rate.isna()|~active.known_rate_unit)).sum()),s_doc_zero_rate_rows=int((active.itemid.isin(PRESSORS)&active.rate.eq(0)&active.known_rate_unit).sum()),s_doc_official5=int(len(officialpos)>0),s_doc_angiotensin_ii=int(len(angi)>0),s_doc_halfopen=int((ag.support_halfopen&ag.itemid.isin(PRESSORS)).any()),documented_inotrope=int((ag.support_left&ag.itemid.isin(INOTROPES)).any()),all_pressor_off_verified=False)
        h=sources['height'].get(sid);w=sources['weight'].get(sid)
        if w is not None:
            w=w.copy();w['weight']=w.valuenum.where(w.valueuom.fillna('kg').str.lower().isin(['kg','kilogram','kilograms']))
        add_measure(o,h,'height','height',lambda z:z.between(100,250))
        add_measure(o,w,'weight','weight',lambda z:z.between(20,400))
        bmi=o['weight']/(o['height']/100)**2 if pd.notna(o['height']) and pd.notna(o['weight']) else np.nan
        o['bmi_source_value']=bmi;o['bmi']=bmi if pd.notna(bmi) and bmi>0 else np.nan
        o['bmi_extreme_review_flag']=bool(pd.notna(bmi) and not 10<=bmi<=80)
        o['bmi_invalid_pair']=False  # Bounds alone cannot identify a measured extreme as an error.
        o['height_age_days']=(x.prebase-o['height_time']).total_seconds()/86400 if pd.notna(o['height_time']) else np.nan
        o['weight_age_hours']=(x.prebase-o['weight_time']).total_seconds()/3600 if pd.notna(o['weight_time']) else np.nan
        pbw=(50 if x.gender=='M' else 45.5)+.91*(o['height']-152.4) if x.gender in ['M','F'] and pd.notna(o['height']) else np.nan
        o['predicted_body_weight']=pbw if pd.notna(pbw) and pbw>0 else np.nan
        o['vt_per_predicted_kg']=o['tidal_volume']/o['predicted_body_weight'] if pd.notna(o['tidal_volume']) and pd.notna(o['predicted_body_weight']) else np.nan
        rows.append(o)
    b=pd.DataFrame(rows)
    d=q.merge(b,on='subject_id',validate='one_to_one')
    d['sex']=d.gender
    d['cohort_C']=d.pf.between(100,300)
    cov=load('preflight_prior_history_coverage.csv',['prior_last_dischtime'])
    d=d.merge(cov,on='subject_id',validate='one_to_one')
    d['prior_history_observed']=d.prior_coded_admissions.gt(0).astype(int)
    diag=load('preflight_prior_diagnoses_raw.csv',['prior_dischtime'])
    diag.icd_code=diag.icd_code.astype(str).str.replace('.','',regex=False).str.upper()
    dg=grouped(diag,'subject_id')
    prefix={
        'copd':(['491','492','496'],['J41','J42','J43','J44']),
        'chf':(['428','40201','40211','40291','40401','40403','40411','40413','40491','40493'],['I50','I110','I130','I132']),
        'pulmonary_hypertension':(['4160','4168','4169'],['I270','I272']),
        'cad':(['410','411','412','413','414'],['I20','I21','I22','I23','I24','I25']),
        'ckd':(['585'],['N18'])}
    for name,(i9,i10) in prefix.items():
        flags=[]
        for sid in d.subject_id:
            g=dg.get(int(sid));present=False
            if g is not None:present=bool(g.loc[g.icd_version.eq(9),'icd_code'].str.startswith(tuple(i9)).any() or g.loc[g.icd_version.eq(10),'icd_code'].str.startswith(tuple(i10)).any())
            flags.append(int(present))
        d['prior_documented_'+name]=flags;d[name]=flags
    d['true_disease_status_verified']=False
    outcomes=load('preflight_outcomes_raw.csv',['dod','admittime','dischtime','deathtime','last_dischtime','all_admissions_deathtime','last_recorded_deathtime'])
    d=d.merge(outcomes,on='subject_id',validate='one_to_one')
    records=[]
    for x in d.itertuples(index=False):
        end=x.t0+pd.Timedelta(days=28);exact=x.all_admissions_deathtime
        y=np.nan;lo=np.nan;hi=np.nan;invalid=False;reason='followup_unknown'
        conflict=bool(pd.notna(exact) and pd.notna(x.dod) and exact.date()!=x.dod.date())
        if pd.notna(exact):
            if exact<=x.t0:invalid=True;reason='precise_death_not_after_T0'
            else:y=float(exact<=end);lo=y;hi=y;reason='precise_all_admission_death_timestamp'
        elif pd.notna(x.dod):
            if x.dod.date()<x.t0.date():invalid=True;reason='registry_death_date_before_T0'
            elif x.dod.date()==x.t0.date():invalid=True;reason='same_day_date_only_landmark_order_unknown'
            elif x.dod.date()==end.date():lo=0.;hi=1.;reason='day28_registry_date_boundary_unknown'
            else:y=float(x.dod.date()<end.date());lo=y;hi=y;reason='registry_date_unambiguous'
        elif pd.notna(x.last_dischtime) and x.last_dischtime+pd.Timedelta(days=365)>=end:
            y=0.;lo=0.;hi=0.;reason='no_registry_death_within_documented_followup'
        records.append({'subject_id':x.subject_id,'death28':y,'death28_lower':lo,'death28_upper':hi,'invalid_landmark':invalid,'mortality_status_reason':reason,'death_source_date_conflict':conflict})
    d=d.merge(pd.DataFrame(records),on='subject_id',validate='one_to_one')
    d['hospital_death']=d.hospital_expire_flag
    # Outcome validity is applied once to the frozen pair. No subsequent TTE is selected.
    d[d.invalid_landmark].to_csv(ROOT/'data/preflight_invalid_landmarks.csv',index=False)
    allset=d[~d.invalid_landmark].copy();known=allset[allset.death28.notna()].copy()
    allset.to_csv(ROOT/'data/preflight_n_all.csv',index=False)
    known.to_csv(ROOT/'data/preflight_yknown.csv',index=False)
    # No outcome-stratified descriptions or association/effect calculations.
    variables=['age','sex','bmi','copd','chf','pulmonary_hypertension','cad','ckd','prior_history_observed','noncardiopulmonary_sofa','ncp_sofa_source_aware','map','lactate','peep','minute_volume','vent_hours','pf','s_doc','nee','rr','tidal_volume','plateau','driving_pressure','creatinine']
    missing=pd.DataFrame([{'variable':k,'n':len(allset),'observed':int(allset[k].notna().sum()),'missing':int(allset[k].isna().sum()),'missing_percent':float(allset[k].isna().mean()*100)} for k in variables])
    missing.to_csv(ROOT/'tables/preflight_covariate_missingness.csv',index=False)
    d.groupby(['invalid_landmark','mortality_status_reason'],dropna=False).size().rename('n').reset_index().to_csv(ROOT/'tables/preflight_outcome_eligibility_reasons.csv',index=False)
    allset.groupby(['fio2_source','pf_missing_reason'],dropna=False).size().rename('n').reset_index().to_csv(ROOT/'tables/preflight_pf_source_counts.csv',index=False)
    press.groupby(['drug','rateuom','statusdescription','known_rate_unit'],dropna=False).agg(rows=('itemid','size'),positive_active_rows=('support_left','sum')).reset_index().to_csv(ROOT/'tables/preflight_pressor_units_status.csv',index=False)
    checks={}
    for col in ['map_time','lactate_time','creatinine_time','peep_time','minute_volume_time','plateau_time','rr_time','tidal_volume_time','driving_pressure_time','height_time','weight_time']:
        checks[col+'_not_before_L']=int((T(d[col])>=d.prebase).sum())
    checks['sofa_after_L']=int((T(d.sofa_time_used)>d.prebase).sum())
    checks['fio2_30m_fallback_after_BG']=int((T(d.fio2_vent30_time)>d.bg_time).sum())
    checks['fio2_2h_fallback_after_BG']=int((T(d.fio2_vent2h_time)>d.bg_time).sum())
    checks['fio2_primary_gap_gt120m']=int((d.fio2_source.eq('same_episode_vent_le2h')&d.fio2_gap_min.gt(120)).sum())
    checks['fio2_30m_sensitivity_gap_gt30m']=int((d.fio2_30m_source.eq('same_episode_vent_le30min')&d.fio2_30m_gap_min.gt(30)).sum())
    checks['negative_vent_hours']=int(d.vent_hours.lt(0).sum())
    checks['prior_discharge_not_before_L']=int((T(d.prior_last_dischtime)>=d.prebase).sum())
    assert not any(checks.values()),checks
    summary={'input_N':len(q),'N_all':len(allset),'N_Yknown':len(known),'invalid_landmark_exclusions':int(d.invalid_landmark.sum()),'eligible_unknown_Y':int(allset.death28.isna().sum()),'death_source_date_conflicts':int(d.death_source_date_conflict.sum()),'recorded_support_positive':int(allset.s_doc.sum()),'recorded_support_unverified':int(allset.s_doc.eq(0).sum()),'nee_observed':int(allset.nee.notna().sum()),'nee_overlap_conflicts':int(allset.nee_conflict.sum()),'angiotensin_ii_support':int(allset.s_doc_angiotensin_ii.sum()),'physically_incompatible_pf':int(allset.pf_pair_incompatible.sum()),'no_outcome_association_models':True,'no_outcome_stratified_rates':True,'checks':checks,'input_sha256':provenance['input_sha256']}
    summary.update(primary_fio2_window_minutes=120,midnight_primary_retained=True,literal_midnight_Nall=int(allset.echo_midnight.sum()),pf_observed=int(allset.pf.notna().sum()),pf_30m_raw_observed=int(allset.pf_30m_raw.notna().sum()),pf_2h_qc_observed=int(allset.pf_2h_qc_sensitivity.notna().sum()),pf_30m_qc_observed=int(allset.pf_30m_qc_sensitivity.notna().sum()))
    summary['output_sha256']={name:hashlib.sha256((ROOT/'data'/name).read_bytes()).hexdigest() for name in ['preflight_n_all.csv','preflight_yknown.csv']}
    (ROOT/'logs/baseline_preflight_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    component_qc(allset)
    (ROOT/'logs/baseline_preflight_complete.json').write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),**summary},indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2),flush=True)
    print(missing.to_string(index=False),flush=True)

def component_qc(d):
    h=load('preflight_height_raw.csv');w=load('preflight_weight_raw.csv');v=load('preflight_vent_raw.csv')
    height_bad=h.height.notna()&~h.height.between(100,250)
    weight_bad=w.valuenum.notna()&~w.valuenum.between(20,400)
    rows=[
        {'metric':'selected_fio2_exactly20_percent','n':int(d.fio2_used.eq(20).sum())},
        {'metric':'same_sample_fio2_exactly20_percent','n':int(d.fio2_sample_percent.eq(20).sum())},
        {'metric':'vent30_fio2_exactly20_percent','n':int(d.fio2_vent30.eq(20).sum())},
        {'metric':'any_extracted_vent_fio2_exactly20_percent_rows','n':int(v.fio2.eq(20).sum())},
        {'metric':'height_source_outside100_250_rows','n':int(height_bad.sum())},
        {'metric':'weight_source_outside20_400_rows','n':int(weight_bad.sum())},
        {'metric':'weight_source_outside20_400_persons','n':int(w.loc[weight_bad,'cohort_subject_id'].nunique())},
        {'metric':'Nall_latest_raw_weight_outside20_400','n':int((d.weight_raw_latest.notna()&~d.weight_raw_latest.between(20,400)).sum())},
        {'metric':'Nall_component_BMI_outside10_80','n':int((d.bmi_source_value.notna()&~d.bmi_source_value.between(10,80)).sum())}]
    pd.DataFrame(rows).to_csv(ROOT/'tables/preflight_component_qc_counts.csv',index=False)
    w[weight_bad].groupby('valuenum',dropna=False).agg(source_rows=('valuenum','size'),persons=('cohort_subject_id','nunique')).reset_index().to_csv(ROOT/'tables/preflight_weight_outlier_values.csv',index=False)
    ids=d.loc[(d.bmi_source_value.notna()&~d.bmi_source_value.between(10,80))|d.weight_raw_latest.gt(400),'subject_id']
    columns=['subject_id','height','height_time','weight','weight_time','weight_raw_latest','weight_raw_latest_time','bmi_source_value','bmi']
    review=d[d.subject_id.isin(ids)][columns].copy()
    review['classification']='source_review_required_not_automatically_error'
    review.to_csv(ROOT/'data/preflight_component_qc_review.csv',index=False)
    w[w.cohort_subject_id.isin(ids)].to_csv(ROOT/'data/preflight_component_review_weight_sources.csv',index=False)

if __name__=='__main__':main()
