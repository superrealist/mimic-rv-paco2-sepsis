"""Read-only reconstruction of Sepsis culture-time precision; no outcomes/models."""
from pathlib import Path
import subprocess,json,hashlib,datetime
import pandas as pd,numpy as np

from .common import ROOT,copy,run_sql

BASE="""SELECT i.subject_id,i.hadm_id,i.stay_id,i.intime,i.outtime,
p.anchor_age+extract(year from i.intime)::int-p.anchor_year age,
s.antibiotic_time,s.culture_time,s.suspected_infection_time,s.sofa_time,s.sepsis3
FROM mimiciv_icu.icustays i JOIN mimiciv_hosp.patients p USING(subject_id)
JOIN mimiciv_derived.sepsis3 s ON s.stay_id=i.stay_id AND s.sepsis3
WHERE p.anchor_age+extract(year from i.intime)::int-p.anchor_year>=18"""

RECON=f"""WITH b AS MATERIALIZED ({BASE}),
me AS MATERIALIZED (
 SELECT m.micro_specimen_id,max(m.subject_id) subject_id,max(m.hadm_id) micro_hadm_id,
 max(m.chartdate)::date chartdate,max(m.charttime) charttime,
 max(m.spec_type_desc) specimen,
 max(CASE WHEN m.org_name IS NOT NULL AND m.org_itemid<>90856 AND m.org_name<>'' THEN 1 ELSE 0 END) positiveculture,
 count(*) raw_rows,count(DISTINCT m.subject_id) n_subjects,count(DISTINCT m.hadm_id) n_hadm,
 count(DISTINCT m.chartdate) n_chartdate,count(DISTINCT m.charttime) n_charttime,
 min(m.charttime) min_charttime,max(m.charttime) max_charttime
 FROM mimiciv_hosp.microbiologyevents m
 WHERE m.subject_id IN (SELECT subject_id FROM b)
 GROUP BY m.micro_specimen_id
),
before_rank AS (
 SELECT b.stay_id,me.*,'before_or_date_window'::text matching_phase,
 rank() OVER(PARTITION BY b.stay_id ORDER BY me.chartdate,me.charttime NULLS LAST) source_rank
 FROM b JOIN me ON me.subject_id=b.subject_id AND (
 (me.charttime IS NOT NULL AND b.antibiotic_time>me.charttime AND b.antibiotic_time<=me.charttime+interval '72 hour')
 OR (me.charttime IS NULL AND b.antibiotic_time::date>=me.chartdate AND b.antibiotic_time::date<=me.chartdate+3)
 )
), before_first AS MATERIALIZED(SELECT * FROM before_rank WHERE source_rank=1),
after_rank AS (
 SELECT b.stay_id,me.*,'after_or_date_window'::text matching_phase,
 rank() OVER(PARTITION BY b.stay_id ORDER BY me.chartdate,me.charttime NULLS LAST) source_rank
 FROM b JOIN me ON me.subject_id=b.subject_id AND (
 (me.charttime IS NOT NULL AND b.antibiotic_time>=me.charttime-interval '24 hour' AND b.antibiotic_time<me.charttime)
 OR (me.charttime IS NULL AND b.antibiotic_time::date>=me.chartdate-1 AND b.antibiotic_time::date<=me.chartdate)
 )
), chosen AS (
 SELECT * FROM before_first
 UNION ALL SELECT * FROM after_rank a WHERE a.source_rank=1 AND NOT EXISTS(SELECT 1 FROM before_first f WHERE f.stay_id=a.stay_id)
)
SELECT b.*,c.micro_specimen_id,c.micro_hadm_id,c.chartdate source_chartdate,c.charttime source_charttime,
c.specimen source_specimen,c.positiveculture,c.raw_rows,c.n_subjects,c.n_hadm,c.n_chartdate,c.n_charttime,
c.min_charttime,c.max_charttime,c.matching_phase,
coalesce(c.charttime,c.chartdate::timestamp)=b.culture_time culture_time_matches,
CASE WHEN c.charttime IS NOT NULL THEN c.charttime ELSE c.chartdate::timestamp END candidate_lower,
CASE WHEN c.charttime IS NOT NULL THEN c.charttime ELSE c.chartdate::timestamp+interval '1 day' END candidate_upper
FROM b LEFT JOIN chosen c USING(stay_id)
ORDER BY b.stay_id,c.micro_specimen_id"""

def extract():
    sql=copy(BASE,'data/time_adult_sepsis_stays.csv')
    sql+=copy(RECON,'data/time_culture_reconstruction.csv')
    sql+=copy(f"WITH b AS ({BASE}) SELECT DISTINCT s.subject_id,s.stay_id,s.hadm_id,s.ab_id,s.antibiotic,s.antibiotic_time,s.culture_time,s.suspected_infection_time,s.specimen,s.positive_culture FROM mimiciv_derived.suspicion_of_infection s JOIN b ON s.subject_id=b.subject_id AND s.stay_id=b.stay_id AND s.antibiotic_time=b.antibiotic_time AND s.culture_time IS NOT DISTINCT FROM b.culture_time AND s.suspected_infection_time IS NOT DISTINCT FROM b.suspected_infection_time",'data/time_matching_suspicion_rows.csv')
    sql+=copy("SELECT test_type,count(*) studies,count(*) FILTER(WHERE measurement_datetime::time=time '00:00:00') midnight_studies,count(*) FILTER(WHERE measurement_datetime IS NULL) missing_time FROM mimiciv_echo.measurement_events GROUP BY test_type ORDER BY test_type",'tables/time_echo_time_distribution.csv')
    run_sql('source_extract',sql)

def summarize():
    c=pd.read_csv(ROOT/'data'/'time_culture_reconstruction.csv')
    for k in ['intime','outtime','antibiotic_time','culture_time','suspected_infection_time','sofa_time','source_chartdate','source_charttime','candidate_lower','candidate_upper']:c[k]=pd.to_datetime(c[k])
    soi=pd.read_csv(ROOT/'data'/'time_matching_suspicion_rows.csv')
    soi['specimen']=soi.specimen.fillna('')
    specs={k:set(g.specimen) for k,g in soi.groupby('stay_id')}
    c['specimen_matches']=c.apply(lambda r:r.source_specimen in specs.get(r.stay_id,set()),axis=1)
    c['time_match']=c.culture_time_matches.fillna(False).astype(bool)
    c['source_metadata_conflict']=c[['n_subjects','n_hadm','n_chartdate','n_charttime']].fillna(1).gt(1).any(axis=1)
    c['cross_hadm']=c.micro_hadm_id.notna() & c.micro_hadm_id.ne(c.hadm_id)
    rows=[]
    for sid,g in c.groupby('stay_id',sort=False):
        b=g.iloc[0]
        # Preserve all tied source identities. Different specimen labels at the same
        # canonical time do not imply ambiguity in that time. Never select a
        # favorable exact-time candidate from conflicting point/date candidates.
        selected=g[g.time_match & ~g.source_metadata_conflict & g.micro_specimen_id.notna()]
        alltopvalid=(len(selected)==len(g)) and selected.micro_specimen_id.notna().all()
        if selected.empty:
            status='source_unresolved';lo=pd.NaT;hi=pd.NaT
        elif not alltopvalid:
            status='ambiguous_top_source';lo=selected.candidate_lower.min();hi=selected.candidate_upper.max()
        elif selected.source_charttime.notna().all():
            status='exact' if len(selected)==1 else 'exact_time_tied_specimens'
            lo=selected.candidate_lower.min();hi=selected.candidate_upper.max()
        elif selected.source_charttime.isna().all():
            status='date_only' if len(selected)==1 else 'date_only_tied_specimens'
            lo=selected.candidate_lower.min();hi=selected.candidate_upper.max()
        else:
            status='mixed_precision_tied_sources';lo=selected.candidate_lower.min();hi=selected.candidate_upper.max()
        trusted=alltopvalid and status not in ['source_unresolved','ambiguous_top_source']
        # Source ties may be identity-ambiguous but their combined time interval remains auditable.
        onsetlo=max([b.antibiotic_time,b.suspected_infection_time,b.sofa_time,lo]) if pd.notna(lo) else pd.NaT
        onsethi=max([b.antibiotic_time,b.suspected_infection_time,b.sofa_time,hi]) if pd.notna(hi) else pd.NaT
        rows.append(dict(subject_id=int(b.subject_id),hadm_id=int(b.hadm_id),stay_id=int(sid),intime=b.intime,outtime=b.outtime,age=b.age,antibiotic_time=b.antibiotic_time,culture_time=b.culture_time,suspected_infection_time=b.suspected_infection_time,sofa_time=b.sofa_time,culture_time_lower=lo,culture_time_upper=hi,sepsis_components_lower=onsetlo,sepsis_components_upper=onsethi,culture_precision=status,culture_source_trusted=trusted,top_source_candidates=len(g),matched_valid_candidates=len(selected),canonical_specimen_matches=int(g.specimen_matches.sum()),source_specimen_ids=';'.join(str(int(s)) for s in selected.micro_specimen_id.dropna()),cross_hadm_source=bool(selected.cross_hadm.any()),source_metadata_conflict=bool(g.source_metadata_conflict.any()),culture_lower=lo,culture_upper=hi,precision=status,unresolved_link=not trusted,sepsis_complete_upper=onsethi,culture_upper_inclusive=bool(selected.source_charttime.notna().all()) if trusted else False))
    d=pd.DataFrame(rows)
    d.to_csv(ROOT/'data'/'time_sepsis_precision_by_stay.csv',index=False)
    c.to_csv(ROOT/'data'/'time_culture_candidates_audited.csv',index=False)
    d.groupby(['culture_precision','culture_source_trusted']).agg(stays=('stay_id','size'),subjects=('subject_id','nunique'),cross_hadm=('cross_hadm_source','sum')).reset_index().to_csv(ROOT/'tables'/'time_sepsis_precision_counts.csv',index=False)
    d[d.culture_time.dt.time.eq(datetime.time(0,0))].groupby('culture_precision').agg(stays=('stay_id','size'),subjects=('subject_id','nunique')).reset_index().to_csv(ROOT/'tables'/'time_midnight_culture_precision.csv',index=False)
    summary={'stays':len(d),'subjects':int(d.subject_id.nunique()),'trusted_interval_stays':int(d.culture_source_trusted.sum()),'unresolved_stays':int((~d.culture_source_trusted).sum()),'date_only_stays':int(d.culture_precision.str.startswith('date_only').sum()),'exact_time_stays':int(d.culture_precision.str.startswith('exact').sum()),'cross_hadm_stays':int(d.cross_hadm_source.sum()),'metadata_conflict_stays':int(d.source_metadata_conflict.sum()),'no_outcomes_read':True,'no_association_models_run':True}
    (ROOT/'logs'/'time_precision_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))
    print(d.groupby('culture_precision').size().to_string())

if __name__=='__main__':
    import sys
    if '--summarize-only' not in sys.argv:extract()
    summarize()
