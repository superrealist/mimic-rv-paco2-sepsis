"""Export source records used to adjudicate organ-support-free days.

Run from the repository root after extraction has produced analysis_input.csv.
The clinical time origin is echocardiography. The ventilation and vasopressor
source-query windows retain the original later-pair extraction anchor; CRRT
source-query windows use echocardiography.
"""
from pathlib import Path
import argparse
import pandas as pd

from extraction.common import DATA_DIR, OUTPUT_DIR, copy_query


def cohort_cte(frame, anchor, columns):
    rows = []
    for row in frame.itertuples():
        stamp = pd.Timestamp(getattr(row, anchor))
        if pd.isna(stamp):
            raise ValueError(f'Missing {anchor}')
        values = [str(int(getattr(row, c))) for c in columns]
        values.append("timestamp '" + stamp.isoformat(sep=' ') + "'")
        rows.append('(' + ','.join(values) + ')')
    names = ','.join(columns + ['t0'])
    return f'WITH c({names}) AS (VALUES ' + ','.join(rows) + ') '


def build_queries(frame):
    frame = frame.copy()
    if 't0_prior_later_pair' not in frame:
        frame['t0_prior_later_pair'] = frame[['echo_time', 'bg_time']].apply(pd.to_datetime).max(axis=1)
    assert frame.subject_id.is_unique
    v = frame.rename(columns={'subject_id':'cohort_subject_id','hadm_id':'index_hadm_id','stay_id':'index_stay_id'})
    prefix = cohort_cte(v, 't0_prior_later_pair', ['cohort_subject_id','index_hadm_id','index_stay_id'])
    ad = "a.subject_id=c.cohort_subject_id AND (a.hadm_id=c.index_hadm_id OR (a.admittime<c.t0+interval '30 day' AND a.dischtime>=c.t0-interval '24 hour'))"
    st = "s.subject_id=c.cohort_subject_id AND s.intime<c.t0+interval '30 day' AND s.outtime>=c.t0-interval '24 hour'"
    queries = {
        'ventilation_intervals_raw': prefix + f"SELECT c.cohort_subject_id,c.index_hadm_id,c.index_stay_id,s.subject_id,s.hadm_id,v.* FROM c JOIN mimiciv_icu.icustays s ON {st} JOIN mimiciv_derived.ventilation v ON v.stay_id=s.stay_id WHERE v.endtime>=c.t0-interval '24 hour' AND v.starttime<=c.t0+interval '30 day'",
        'ventilator_setting_raw': prefix + f"SELECT c.cohort_subject_id,c.index_hadm_id,c.index_stay_id,s.hadm_id,v.* FROM c JOIN mimiciv_icu.icustays s ON {st} JOIN mimiciv_derived.ventilator_setting v ON v.stay_id=s.stay_id WHERE v.charttime>=c.t0-interval '24 hour' AND v.charttime<=c.t0+interval '30 day'",
        'airway_chartevents_raw': prefix + f"SELECT c.cohort_subject_id,c.index_hadm_id,c.index_stay_id,v.* FROM c JOIN mimiciv_icu.icustays s ON {st} JOIN mimiciv_icu.chartevents v ON v.stay_id=s.stay_id WHERE v.charttime>=c.t0-interval '24 hour' AND v.charttime<=c.t0+interval '30 day' AND v.itemid IN (223836,223848,223849,229314,226260,226732,224829,223059,227567,227568,227577,227578,227579,227580,227583,226814,226815)",
        'airway_procedureevents_raw': prefix + f"SELECT c.cohort_subject_id,c.index_hadm_id,c.index_stay_id,v.* FROM c JOIN mimiciv_hosp.admissions a ON {ad} JOIN mimiciv_icu.procedureevents v ON v.hadm_id=a.hadm_id JOIN mimiciv_icu.d_items di ON di.itemid=v.itemid WHERE v.starttime<=c.t0+interval '30 day' AND v.endtime>=c.t0-interval '24 hour' AND (di.category IN ('1-Intubation/Extubation','2-Ventilation') OR lower(di.label) ~ 'extubat|intubat|tracheo|tracheot|ventilat')",
        'all_icu_stays_30d_raw': prefix + f"SELECT c.cohort_subject_id,c.index_hadm_id,c.index_stay_id,s.* FROM c JOIN mimiciv_icu.icustays s ON {st}",
        'all_admissions_30d_raw': prefix + f"SELECT c.cohort_subject_id,c.index_hadm_id,c.index_stay_id,a.* FROM c JOIN mimiciv_hosp.admissions a ON {ad}",
    }
    pressor_prefix = cohort_cte(frame, 't0_prior_later_pair', ['subject_id'])
    queries['pressors_28'] = pressor_prefix + """SELECT c.subject_id,c.t0,z.stay_id,z.starttime,z.endtime,s.intime,s.outtime,z.itemid,z.rate,z.rateuom,z.statusdescription,z.ordercategorydescription,z.continueinnextdept
FROM c JOIN mimiciv_icu.icustays s ON s.subject_id=c.subject_id AND s.outtime>=c.t0 AND s.intime<=c.t0+interval '28 day'
JOIN mimiciv_icu.inputevents z ON z.stay_id=s.stay_id AND z.subject_id=c.subject_id AND z.hadm_id=s.hadm_id WHERE z.itemid IN (221906,221289,229617,221749,229630,229631,229632,221662,222315,229709,229764) AND z.starttime<=c.t0+interval '28 day' AND z.endtime>=s.intime"""
    emar = """SELECT DISTINCT c.subject_id,e.charttime,e.medication,e.event_txt,d.infusion_rate,d.infusion_rate_unit,d.route,d.new_iv_bag_hung,d.infusion_complete FROM c JOIN mimiciv_hosp.emar e ON e.subject_id=c.subject_id AND e.charttime>=c.t0 AND e.charttime<=c.t0+interval '28 days' JOIN mimiciv_hosp.emar_detail d ON d.subject_id=e.subject_id AND d.emar_id=e.emar_id AND d.emar_seq=e.emar_seq WHERE lower(e.medication)~'(norepinephrine|epinephrine|vasopressin|phenylephrine|dopamine|angiotensin ii|levophed|neosynephrine|vasostrict|giapreza)' AND (nullif(d.infusion_rate,'') IS NOT NULL OR d.new_iv_bag_hung='Yes')"""
    queries['emar_28'] = pressor_prefix + emar
    queries['emar_28_bag'] = pressor_prefix + emar.replace("(nullif(d.infusion_rate,'') IS NOT NULL OR d.new_iv_bag_hung='Yes')", "nullif(d.infusion_rate,'') IS NULL AND d.new_iv_bag_hung='Y'")
    renal_prefix = cohort_cte(frame.rename(columns={'subject_id':'sid'}), 'echo_time', ['sid'])
    renal_st = "s.subject_id=c.sid AND s.outtime>=c.t0-interval '24 hour' AND s.intime<=c.t0+interval '28 day'"
    queries['crrt'] = renal_prefix + f"SELECT c.sid as subject_id,c.t0,s.intime,s.outtime,z.* FROM c JOIN mimiciv_icu.icustays s ON {renal_st} JOIN mimiciv_derived.crrt z ON z.stay_id=s.stay_id WHERE z.charttime BETWEEN c.t0-interval '24 hour' AND c.t0+interval '28 day'"
    queries['rrt'] = renal_prefix + f"SELECT c.sid as subject_id,c.t0,z.* FROM c JOIN mimiciv_icu.icustays s ON {renal_st} JOIN mimiciv_derived.rrt z ON z.stay_id=s.stay_id WHERE z.charttime BETWEEN c.t0-interval '24 hour' AND c.t0+interval '28 day'"
    queries['crrt_proc'] = renal_prefix + f"SELECT c.sid as subject_id,c.t0,s.intime,s.outtime,z.stay_id,z.itemid,z.starttime,z.endtime,z.statusdescription,z.continueinnextdept FROM c JOIN mimiciv_icu.icustays s ON {renal_st} JOIN mimiciv_icu.procedureevents z ON z.stay_id=s.stay_id WHERE z.itemid IN(225802,225803,225809,225955,225441) AND z.endtime>=c.t0-interval '24 hour' AND z.starttime<=c.t0+interval '28 day'"
    queries['crrt_system'] = renal_prefix + f"SELECT c.sid as subject_id,z.stay_id,z.charttime,z.itemid,z.value,z.valuenum FROM c JOIN mimiciv_icu.icustays s ON {renal_st} JOIN mimiciv_icu.chartevents z ON z.stay_id=s.stay_id JOIN mimiciv_icu.d_items di ON z.itemid=di.itemid WHERE z.charttime BETWEEN c.t0-interval '24 hour' AND c.t0+interval '28 day' AND z.itemid=224146"
    return queries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--only', nargs='*', help='Optional source names to export')
    args = parser.parse_args()
    frame = pd.read_csv(DATA_DIR/'analysis_input.csv', low_memory=False)
    queries = build_queries(frame)
    names = args.only or list(queries)
    unknown = set(names) - set(queries)
    if unknown:
        raise ValueError(f'Unknown source names: {sorted(unknown)}')
    target_dir = DATA_DIR/'outcomes'/'raw'
    target_dir.mkdir(parents=True, exist_ok=True)
    for name in names:
        copy_query(queries[name], target_dir/f'{name}.csv', name=name)


if __name__ == '__main__':
    main()
