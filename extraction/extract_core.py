from .common import ROOT,copy,run_sql

BASE="""SELECT i.subject_id,i.hadm_id,i.stay_id,i.intime,i.outtime,
p.gender,p.anchor_age,p.anchor_year,p.anchor_year_group,
p.anchor_age+extract(year from i.intime)::int-p.anchor_year age_derived,
s.suspected_infection_time,s.antibiotic_time,s.culture_time,s.sofa_time
FROM mimiciv_icu.icustays i JOIN mimiciv_hosp.patients p USING(subject_id)
JOIN mimiciv_derived.sepsis3 s ON s.stay_id=i.stay_id AND s.sepsis3
WHERE p.anchor_age+extract(year from i.intime)::int-p.anchor_year>=18"""
ANCHORS=f"""WITH b AS({BASE}) SELECT b.*,e.measurement_id,e.measurement_datetime echo_time,
v.starttime vent_start,v.endtime vent_end,
(SELECT count(*) FROM mimiciv_icu.icustays ii WHERE ii.subject_id=b.subject_id AND e.measurement_datetime>=ii.intime AND e.measurement_datetime<ii.outtime) icu_matches
FROM b JOIN mimiciv_derived.ventilation v ON v.stay_id=b.stay_id AND v.ventilation_status='InvasiveVent'
JOIN mimiciv_echo.measurement_events e ON e.subject_id=b.subject_id AND e.test_type='tte'
AND e.measurement_datetime>=b.intime AND e.measurement_datetime<b.outtime
AND e.measurement_datetime>=v.starttime AND e.measurement_datetime<v.endtime"""

def main():
    sql=""
    sql+=copy(BASE,'data/adult_sepsis_stays.csv')
    sql+=copy(ANCHORS,'data/echo_anchors.csv')
    sql+=copy(f"WITH a AS({ANCHORS}) SELECT sm.* FROM mimiciv_echo.structured_measurement sm WHERE measurement_id IN(SELECT measurement_id FROM a) AND measurement IN('rv_function','tapse','rv_chamber_size','rv_diam','tr_velocity','tr_mmhg')",'data/rv_raw.csv')
    span=f"SELECT hadm_id,min(echo_time)-interval '12 hour' lo,max(echo_time)+interval '12 hour' hi FROM ({ANCHORS}) a GROUP BY hadm_id"
    sql+=copy(f"WITH s AS({span}) SELECT b.* FROM mimiciv_derived.bg b JOIN s USING(hadm_id) WHERE b.charttime BETWEEN s.lo AND s.hi",'data/blood_gases.csv')
    sql+=copy("SELECT table_schema,table_name,column_name,data_type FROM information_schema.columns WHERE (table_schema='mimiciv_hosp' AND table_name IN('patients','admissions','labevents','microbiologyevents','diagnoses_icd')) OR (table_schema='mimiciv_derived' AND table_name IN('sepsis3','suspicion_of_infection','sofa','bg','ventilation','ventilator_setting','norepinephrine_equivalent_dose','chemistry','height','weight_durations')) OR (table_schema='mimiciv_echo' AND table_name IN('structured_measurement','measurement_events')) ORDER BY 1,2,ordinal_position",'tables/source_columns.csv')
    run_sql('00_core_preflight',sql)

if __name__=='__main__':main()
