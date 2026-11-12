from .common import copy,run_sql
from .extract_core import ANCHORS

def main():
    span=f"SELECT hadm_id,min(echo_time)-interval '12 hour' lo,max(echo_time)+interval '12 hour' hi FROM({ANCHORS})a GROUP BY hadm_id"
    sql=copy(f"WITH s AS({span}) SELECT l.labevent_id,l.subject_id,l.hadm_id,l.specimen_id,l.itemid,l.charttime,l.storetime,l.value,l.valuenum,l.valueuom FROM mimiciv_hosp.labevents l JOIN s USING(hadm_id) WHERE l.charttime BETWEEN s.lo AND s.hi AND l.itemid IN(52033,50818,50820,50821,50816,50813,50803)",'data/raw_gas_labevents.csv')
    sql+=copy("SELECT itemid,label,fluid,category FROM mimiciv_hosp.d_labitems WHERE itemid IN(52033,50818,50820,50821,50816,50813,50803)",'tables/blood_gas_item_dictionary.csv')
    run_sql('01_raw_gas_and_lv_context',sql)

if __name__=='__main__':main()
