"""Batched read-only V2 baseline and outcome-ascertainment exports, no models."""
from .common import ROOT,copy,run_sql
import pandas as pd,json,datetime,hashlib,sys

def main():
    if '--sofa-availability-only' in sys.argv:
        return audit_sofa_availability()
    path=ROOT/'data/cohort_primary_pre_outcome_validation.csv'
    d=pd.read_csv(path)
    assert d.subject_id.is_unique and d.stay_id.is_unique and len(d)>0
    cols=['subject_id','hadm_id','stay_id','measurement_id','prebase','t0','bg_time','vent_start','vent_end']
    assert d[cols].notna().all().all()
    values=','.join("(%d,%d,%d,%d,timestamp '%s',timestamp '%s',timestamp '%s',timestamp '%s',timestamp '%s')" % tuple(r) for r in d[cols].itertuples(index=False,name=None))
    c=f'c(subject_id,hadm_id,stay_id,measurement_id,prebase,t0,bg_time,vent_start,vent_end) AS (VALUES {values})'
    marker=ROOT/'logs/baseline_preflight_extraction_complete.json'
    if marker.exists():marker.unlink()
    sql=''
    for table,condition,target in [
        ('sofa',"z.endtime BETWEEN c.prebase-interval '24 hour' AND c.prebase",'sofa'),
        ('vitalsign',"z.charttime>=c.prebase-interval '2 hour' AND z.charttime<c.prebase",'vitals'),
        ('ventilator_setting',"z.charttime>=greatest(c.vent_start,c.prebase-interval '6 hour') AND z.charttime<=c.bg_time AND z.charttime<=c.vent_end",'vent'),
        ('norepinephrine_equivalent_dose',"z.starttime<c.prebase AND z.endtime>=c.prebase",'nee')]:
        sql+=copy(f'WITH {c} SELECT c.subject_id cohort_subject_id,c.stay_id cohort_stay_id,z.* FROM c JOIN mimiciv_derived.{table} z ON z.stay_id=c.stay_id AND {condition}',f'data/preflight_{target}_raw.csv')
    sql+=copy(f"WITH {c} SELECT c.subject_id cohort_subject_id,z.* FROM c JOIN mimiciv_derived.bg z ON z.subject_id=c.subject_id AND z.hadm_id=c.hadm_id AND z.charttime<c.prebase AND z.charttime>=c.prebase-interval '12 hour' WHERE z.lactate IS NOT NULL",'data/preflight_lactate_raw.csv')
    sql+=copy(f"WITH {c} SELECT c.subject_id cohort_subject_id,z.* FROM c JOIN mimiciv_derived.chemistry z ON z.hadm_id=c.hadm_id AND z.charttime<c.prebase AND z.charttime>=c.prebase-interval '24 hour' WHERE z.creatinine IS NOT NULL",'data/preflight_chemistry_raw.csv')
    sql+=copy(f"WITH {c} SELECT c.subject_id cohort_subject_id,z.subject_id,z.hadm_id,z.stay_id,z.itemid,z.starttime,z.endtime,z.storetime,z.rate,z.rateuom,z.amount,z.amountuom,z.orderid,z.linkorderid,z.statusdescription,z.ordercategoryname,z.ordercategorydescription,z.patientweight,z.originalrate FROM c JOIN mimiciv_icu.inputevents z ON z.stay_id=c.stay_id AND z.starttime<=c.prebase AND z.endtime>=c.prebase WHERE z.itemid IN (221906,221289,221749,221662,222315,229709,229764,221653,221986)",'data/preflight_pressors_raw.csv')
    sql+=copy(f"WITH {c} SELECT c.subject_id,count(DISTINCT a.hadm_id) prior_admissions,count(DISTINCT a.hadm_id) FILTER(WHERE d.icd_code IS NOT NULL) prior_coded_admissions,count(d.icd_code) prior_icd_rows,max(a.dischtime) prior_last_dischtime FROM c LEFT JOIN mimiciv_hosp.admissions a ON a.subject_id=c.subject_id AND a.hadm_id<>c.hadm_id AND a.dischtime<c.prebase LEFT JOIN mimiciv_hosp.diagnoses_icd d ON d.hadm_id=a.hadm_id GROUP BY c.subject_id",'data/preflight_prior_history_coverage.csv')
    sql+=copy(f"WITH {c} SELECT c.subject_id,a.hadm_id prior_hadm_id,a.dischtime prior_dischtime,d.icd_code,d.icd_version FROM c JOIN mimiciv_hosp.admissions a ON a.subject_id=c.subject_id AND a.hadm_id<>c.hadm_id AND a.dischtime<c.prebase JOIN mimiciv_hosp.diagnoses_icd d ON d.hadm_id=a.hadm_id",'data/preflight_prior_diagnoses_raw.csv')
    sql+=copy(f"WITH {c} SELECT c.subject_id cohort_subject_id,h.* FROM c JOIN mimiciv_derived.height h ON h.subject_id=c.subject_id AND h.charttime<c.prebase",'data/preflight_height_raw.csv')
    sql+=copy(f"WITH {c} SELECT c.subject_id cohort_subject_id,z.subject_id,z.stay_id,z.charttime,z.storetime,z.itemid,z.value,z.valuenum,z.valueuom FROM c JOIN mimiciv_icu.chartevents z ON z.stay_id=c.stay_id AND z.charttime<c.prebase WHERE z.itemid IN(226512,224639)",'data/preflight_weight_raw.csv')
    sql+=copy(f"WITH {c},all_a AS (SELECT a.subject_id,max(a.dischtime) last_dischtime,min(a.deathtime) all_admissions_deathtime,max(a.deathtime) last_recorded_deathtime,count(*) FILTER(WHERE a.deathtime IS NOT NULL) exact_death_rows FROM mimiciv_hosp.admissions a JOIN c ON c.subject_id=a.subject_id GROUP BY a.subject_id) SELECT c.subject_id,p.dod,a.admittime,a.dischtime,a.deathtime,a.hospital_expire_flag,all_a.last_dischtime,all_a.all_admissions_deathtime,all_a.last_recorded_deathtime,all_a.exact_death_rows FROM c JOIN mimiciv_hosp.patients p USING(subject_id) JOIN mimiciv_hosp.admissions a ON a.hadm_id=c.hadm_id LEFT JOIN all_a ON all_a.subject_id=c.subject_id",'data/preflight_outcomes_raw.csv')
    sql+=copy("SELECT itemid,label,unitname,linksto,category FROM mimiciv_icu.d_items WHERE itemid IN (221906,221289,221749,221662,222315,229709,229764,221653,221986,226512,224639) ORDER BY itemid",'tables/preflight_medication_weight_dictionary.csv')
    if '--resume-tail' in sys.argv:
        for stem in ['sofa','vitals','vent','nee','lactate','chemistry','pressors','height','weight']:
            p=ROOT/f'data/preflight_{stem}_raw.csv'
            assert p.exists() and p.stat().st_size>0,str(p)
            pd.read_csv(p,nrows=1)
        for name in ['preflight_prior_history_coverage.csv','preflight_prior_diagnoses_raw.csv']:
            assert (ROOT/'data'/name).stat().st_size>0
            pd.read_csv(ROOT/'data'/name,nrows=1)
        sql='\n'.join(line for line in sql.splitlines() if 'preflight_outcomes_raw.csv' in line or 'preflight_medication_weight_dictionary.csv' in line)+'\n'
        run_sql('01b_baseline_preflight_tail',sql)
    else:run_sql('01_baseline_preflight_extraction',sql)
    marker.write_text(json.dumps({'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input_rows':len(d),'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'no_association_models':True},indent=2),encoding='utf-8')
    print('All baseline exports complete; input rows',len(d),flush=True)

def audit_sofa_availability():
    q=pd.read_csv(ROOT/'data/cohort_primary_pre_outcome_validation.csv')
    s=pd.read_csv(ROOT/'data/preflight_sofa_raw.csv',parse_dates=['endtime'])
    last=s.sort_values(['endtime','hr']).drop_duplicates('cohort_subject_id',keep='last')
    last=q[['subject_id','hadm_id','stay_id']].merge(last[['cohort_subject_id','endtime']],left_on='subject_id',right_on='cohort_subject_id',validate='one_to_one')
    values=','.join("(%d,%d,%d,timestamp '%s')" % tuple(r) for r in last[['subject_id','hadm_id','stay_id','endtime']].itertuples(index=False,name=None))
    c=f'c(subject_id,hadm_id,stay_id,sofa_end) AS (VALUES {values})'
    parts=[]
    for table,component,col,key in [('complete_blood_count','coagulation','platelet','hadm_id'),('enzyme','liver','bilirubin_total','hadm_id'),('gcs','cns','gcs','stay_id'),('chemistry','renal','creatinine','hadm_id')]:
        parts.append(f"SELECT c.subject_id,c.stay_id,c.sofa_end,'{component}'::text component,'{table}'::text source_table,z.charttime,z.{col}::double precision source_value,NULL::double precision uo_tm_24hr FROM c JOIN mimiciv_derived.{table} z ON z.{key}=c.{key} AND z.charttime>c.sofa_end-interval '24 hour' AND z.charttime<=c.sofa_end WHERE z.{col} IS NOT NULL")
    parts.append("SELECT c.subject_id,c.stay_id,c.sofa_end,'renal'::text component,'urine_output_rate'::text source_table,z.charttime,(z.urineoutput_24hr/z.uo_tm_24hr*24)::double precision source_value,z.uo_tm_24hr::double precision FROM c JOIN mimiciv_derived.urine_output_rate z ON z.stay_id=c.stay_id AND z.charttime>c.sofa_end-interval '24 hour' AND z.charttime<=c.sofa_end WHERE z.urineoutput_24hr IS NOT NULL AND z.uo_tm_24hr BETWEEN 22 AND 30")
    run_sql('01c_sofa_component_availability',copy(f"WITH {c} "+' UNION ALL '.join(parts),'data/preflight_sofa_component_sources.csv'))
    (ROOT/'logs/sofa_component_availability_complete.json').write_text(json.dumps({'input_sha256':hashlib.sha256((ROOT/'data/cohort_primary_pre_outcome_validation.csv').read_bytes()).hexdigest(),'sofa_rows':len(last),'scope':'full selected SOFA preceding24h including pre-ICU contributors','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2),encoding='utf-8')

if __name__=='__main__':main()
