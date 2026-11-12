"""Horizon-specific terminal support-free days under documented cessation rules."""
from pathlib import Path
import numpy as np
import pandas as pd

from extraction.common import DATA_DIR, OUTPUT_DIR
O=DATA_DIR/'outcomes'
V=O/'raw'
OLD=V
SUPPORT=DATA_DIR/'support'
QUALITY=OUTPUT_DIR/'outcomes'
SUPPORT.mkdir(parents=True,exist_ok=True)
QUALITY.mkdir(parents=True,exist_ok=True)
main=pd.read_csv(DATA_DIR/'analysis_input.csv',low_memory=False)
SAFE={'HOME','HOME HEALTH CARE','REHAB','SKILLED NURSING FACILITY','PSYCH FACILITY'}

def dates(d,cols):
    for col in cols:d[col]=pd.to_datetime(d[col],format='mixed')
    return d

ad=dates(pd.read_csv(V/'all_admissions_30d_raw.csv'),['admittime','dischtime'])
adm={k:g for k,g in ad.groupby('subject_id')}
st=dates(pd.read_csv(V/'all_icu_stays_30d_raw.csv'),['intime','outtime'])
stays={k:g for k,g in st.groupby('subject_id')}

def observed(pid,start,end,strict_reentry=True):
    a=adm[pid]
    if strict_reentry:
        if ((stays[pid].intime>start)&(stays[pid].intime<=end)).any():return 'ICU_reentry_without_new_stop'
        if ((a.admittime>start)&(a.admittime<=end)).any():return 'readmission_without_new_stop'
    else:
        discharges=a[(a.dischtime>=start)&(a.dischtime<end)]
        if (~discharges.discharge_location.isin(SAFE)).any():return 'intervening_disposition_not_observable'
    present=a[(a.admittime<=end)&(a.dischtime>=start)]
    if present.empty:return 'disposition_not_observable'
    q=present.sort_values('dischtime').iloc[-1]
    if q.dischtime<end and q.discharge_location not in SAFE:return 'disposition_not_observable'
    return None

drug={221906:'ne',221289:'epi',229617:'epi',221749:'phen',229630:'phen',229631:'phen',229632:'phen',221662:'dop',222315:'vaso',229709:'ang',229764:'ang'}
units={'ne':{'mcg/kg/min'},'epi':{'mcg/kg/min'},'phen':{'mcg/kg/min'},'dop':{'mcg/kg/min'},'vaso':{'units/hour'},'ang':{'ng/kg/min','mcg/kg/min'}}
press=dates(pd.read_csv(OLD/'pressors_28.csv',low_memory=False),['starttime','endtime','intime','outtime'])
press=press[press.ordercategorydescription.eq('Continuous Med') & ~press.statusdescription.isin(['Rewritten','Bolus'])].copy()
press['drug']=press.itemid.map(drug)
press['unit_valid']=[u in units[g] for g,u in zip(press.drug,press.rateuom)]
press['pos']=press.rate.gt(0)&press.unit_valid
press['zero']=press.rate.eq(0)&press.unit_valid
press['stop']=press.statusdescription.isin(['Stopped','Paused'])&press.continueinnextdept.eq(0)&press.endtime.le(press.outtime)
press['s']=press[['starttime','intime']].max(axis=1)
press['e']=press[['endtime','outtime']].min(axis=1)
press=press[press.e>press.s]
pdct={k:g for k,g in press.groupby('subject_id')}
emar=dates(pd.concat([pd.read_csv(OLD/f,low_memory=False) for f in ['emar_28.csv','emar_28_bag.csv']]).drop_duplicates(),['charttime'])
route=emar.route.fillna('').str.lower(); name=emar.medication.fillna('').str.lower()
iv=route.isin(['iv','iv drip','intravenous','iv infusion','iv (intravenous)','iv central line','iv peripheral'])|route.str.contains('intravenous')
liquid=name.str.contains('infusion|premix|mcg/ml|mg/250|units/100')
topical=name.str.contains('lidocaine|bupivacaine|ophth|nasal|ophthalm|topical|racemic|auto-inject|pen injector|desmopressin')
neg=emar.event_txt.fillna('').str.lower().str.contains('not given|held|refused|cancel|stop|complete|due')
rate=pd.to_numeric(emar.infusion_rate,errors='coerce');bag=emar.new_iv_bag_hung.isin(['Y','Yes'])
exact=name.isin(['norepinephrine','epinephrine','vasopressin','phenylephrine','dopamine','angiotensin ii'])
drip=emar.infusion_rate_unit.fillna('').str.lower().isin(['unit/hr','mcg/kg/min','nanograms/kg/minute'])
emar=emar[(iv|liquid|(exact&(drip|bag)))&~topical&~neg&((rate>0)|bag)]
edct={k:g.charttime for k,g in emar.groupby('subject_id')}

crrt=dates(pd.read_csv(V/'crrt.csv',low_memory=False),['charttime','intime','outtime'])
crrt['pos']=crrt.system_active.eq(1)|crrt.blood_flow.gt(0)
cdct={k:g for k,g in crrt.groupby('subject_id')}
system=dates(pd.read_csv(V/'crrt_system.csv'),['charttime'])
sdct={k:g for k,g in system.groupby('subject_id')}
rrt=dates(pd.read_csv(V/'rrt.csv'),['charttime'])
potential=rrt[rrt.dialysis_active.eq(1)&rrt.dialysis_type.isin(['CRRT','CVVH','CVVHD','CVVHDF','SCUF'])][['subject_id','charttime']]
flows=crrt.crrt_mode.notna()|crrt.dialysate_rate.gt(0)|crrt.prefilter_replacement_rate.gt(0)|crrt.postfilter_replacement_rate.gt(0)|crrt.replacement_rate.gt(0)
potential=pd.concat([potential,crrt.loc[flows,['subject_id','charttime']]],ignore_index=True).drop_duplicates()
potential_dict={k:g.charttime for k,g in potential.groupby('subject_id')}
proc=dates(pd.read_csv(V/'crrt_proc.csv'),['starttime','endtime','intime','outtime'])
proc=proc[proc.itemid.isin([225802,225803,225809,225955])&~proc.statusdescription.eq('Rewritten')&proc.endtime.gt(proc.starttime)].copy()
proc['s']=proc[['starttime','intime']].max(axis=1);proc['e']=proc[['endtime','outtime']].min(axis=1)
proc['stop']=proc.statusdescription.isin(['Stopped','Paused','FinishedRunning'])&proc.continueinnextdept.eq(0)&proc.endtime.lt(proc.outtime)
cpdct={k:g for k,g in proc.groupby('subject_id')}

def pressor(r,days):
    t=pd.Timestamp(r.echo_time);h=t+pd.Timedelta(days=days)
    d=getattr(r,f'death{days}_echo')
    if pd.isna(d):return np.nan,'unknown_death_at_horizon',np.nan
    if d:return 0.,'death_by_horizon',np.nan
    ev=pdct.get(r.subject_id,press.iloc[0:0])
    ev=ev[(ev.s<h)&((ev.stay_id==r.stay_id)|(ev.e>t))]
    on=ev[ev.pos]
    if on.empty:return np.nan,'no_confirmed_baseline_or_later_infusion',np.nan
    if on.e.max()>h-pd.Timedelta(hours=48):return 0.,'infusion_in_last_48h',np.nan
    stops=[]
    for _,g in ev.groupby('drug'):
        act=g[g.pos]
        if act.empty:
            if (~g.zero).any():return np.nan,'uncertified_terminal_drug_state',np.nan
            continue
        end=act.e.max();tie=act[act.e.eq(end)]
        zero=g[g.zero&g.s.ge(end)]
        if tie.stop.all():cert=end
        elif len(zero):cert=zero.s.min()
        else:return np.nan,'uncertified_terminal_drug_state',np.nan
        if ((~g.pos)&(~g.zero)&g.e.gt(cert)).any():return np.nan,'uncertified_terminal_drug_state',np.nan
        stops.append(cert)
    final=max(max(stops),t)
    ee=edct.get(r.subject_id,pd.Series(dtype='datetime64[ns]'))
    if ((ee>final)&(ee<=h)).any():return np.nan,'later_emar_infusion_unresolved',np.nan
    reason=observed(r.subject_id,final,h)
    if reason:return np.nan,reason,np.nan
    n=(h-final).total_seconds()/86400
    return n if n>=2 else 0.,'documented_durable_stop' if n>=2 else 'off_less_than_48h',(final-t).total_seconds()/3600

def renal(r,days):
    t=pd.Timestamp(r.echo_time);h=t+pd.Timedelta(days=days)
    d=getattr(r,f'death{days}_echo')
    if pd.isna(d):return np.nan,'unknown_death_at_horizon',np.nan
    if d:return 0.,'death_by_horizon',np.nan
    c=cdct.get(r.subject_id,crrt.iloc[0:0]);c=c[c.charttime<h]
    p=cpdct.get(r.subject_id,proc.iloc[0:0]);p=p[p.s<h]
    positive=c[c.pos];events=p[p.e>t]
    active=positive[positive.charttime>=t]
    sys=sdct.get(r.subject_id,system.iloc[0:0]);sys=sys[sys.charttime<h]
    off=sys[sys.value.eq('Discontinued')&~sys.charttime.isin(positive.charttime)]
    possible=potential_dict.get(r.subject_id,pd.Series(dtype='datetime64[ns]'))
    possible=possible[(possible>=t)&(possible<h)]
    if events.empty and active.empty:
        if len(possible):return np.nan,'unresolved_continuous_RRT_evidence',np.nan
        prior=positive[positive.charttime<t]
        # The baseline pre-T0 buffer must not hide an unresolved ongoing course.
        if len(prior):
            previous=prior.charttime.max()
            prior_stop=((off.charttime>=previous)&(off.charttime<=t)).any() or ((p.stop)&(p.e>=previous)&(p.e<=t)).any()
            if not prior_stop:
                later=off[off.charttime>t]
                if later.empty:return np.nan,'baseline_CRRT_state_uncertain',np.nan
                final=later.charttime.min()
                reason=observed(r.subject_id,final,h)
                if reason:return np.nan,reason,np.nan
                n=(h-final).total_seconds()/86400
                return n if n>=2 else 0.,'documented_durable_stop' if n>=2 else 'off_less_than_48h',(final-t).total_seconds()/3600
        reason=observed(r.subject_id,t,h,strict_reentry=False)
        return (np.nan,reason,np.nan) if reason else (float(days),'no_CRRT_recorded_with_observable_followup',0.)
    last=max([t]+active.charttime.tolist()+events.e.tolist())
    if last>h-pd.Timedelta(hours=48):return 0.,'CRRT_in_last_48h',np.nan
    # Broad derived RRT activity can inherit a catheter or auxiliary infusion.
    # Use it only to block an unsupported OFF claim, never as definite ON.
    last=max([last]+possible.tolist())
    off=off[off.charttime.ge(last)]
    certified=p[p.stop&p.e.ge(last)]
    candidates=off.charttime.tolist()+certified.e.tolist()
    if not candidates:return np.nan,'terminal_CRRT_stop_uncertified',np.nan
    final=max(t,min(candidates))
    reason=observed(r.subject_id,final,h)
    if reason:return np.nan,reason,np.nan
    n=(h-final).total_seconds()/86400
    return n if n>=2 else 0.,'documented_durable_stop' if n>=2 else 'off_less_than_48h',(final-t).total_seconds()/3600

rows=[]
for r in main.itertuples():
    rec={'subject_id':r.subject_id}
    for stem,fn in [('vpfd',pressor),('crrtfd',renal)]:
        for days in [7,28]:
            value,why,stop=fn(r,days)
            rec[f'{stem}{days}']=value;rec[f'{stem}{days}_reason']=why;rec[f'{stem}{days}_stop_h']=stop
    rows.append(rec)
d=pd.DataFrame(rows)
vent=pd.read_csv(O/'vent_adjudication.csv')
d=d.merge(vent[['subject_id','vfd7','vfd28','vfd7_class','vfd28_class']],validate='1:1')
d.to_csv(O/'support_adjudication.csv',index=False)
summ=[]
for stem in ['vpfd','vfd','crrtfd']:
    for days in [7,28]:
        k=f'{stem}{days}';v=d[k];reason=k+('_class' if stem=='vfd' else '_reason')
        assert v.dropna().between(0,days).all()
        assert v[main[f'death{days}_echo'].eq(1)].eq(0).all()
        summ.append({'outcome':k,'n':int(v.notna().sum()),'missing':int(v.isna().sum()),'zero':int(v.eq(0).sum()),'positive':int(v.gt(0).sum())})
        d[reason].value_counts().rename_axis('reason').reset_index(name='n').to_csv(QUALITY/f'{k}_counts.csv',index=False)
        x=main.merge(d[['subject_id',k]].dropna(),validate='1:1')
        x['analysis_y']=x[k];x['recorded_L_drug_count']=x.s_doc_positive_drugs.str.split(';').map(len)
        x.to_csv(SUPPORT/f'{k}_input.csv',index=False)
pd.DataFrame(summ).to_csv(QUALITY/'endpoint_counts.csv',index=False)
print(pd.DataFrame(summ).to_string(index=False))
