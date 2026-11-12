"""Reconstruct documented invasive-ventilation states and terminal free days."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import json, datetime
import numpy as np
import pandas as pd

O=Path(__file__).resolve().parents[1]
SAFE={'HOME','HOME HEALTH CARE','REHAB','SKILLED NURSING FACILITY','PSYCH FACILITY'}
OFFDEV={'None','Nasal cannula','Trach mask','Aerosol-cool','Face tent','High flow nasal cannula','High flow neb','T-piece','Bipap mask','Oxymizer','Non-rebreather','CPAP mask','Venti mask','Medium conc mask','Vapomist','Ultrasonic neb'}
NIVMODES={'NIV','NIV-ST','nCPAP-PS'}
NONACTIVE={'Standby','Ambient'}
EXTUB={227194,225468,225477}
def stamp(v):return pd.Timestamp(v) if pd.notna(v) else None
def finite(v):return v is not None and pd.notna(v) and np.isfinite(v)
def worker(args):
    row,raw,settings,derived,procs,stays,adms=args
    pid=int(row['subject_id']);t0=stamp(row['t0']);ixstay=int(row['stay_id'])
    hh=lambda x:(stamp(x)-t0).total_seconds()/3600 if pd.notna(x) else np.nan
    icu={int(z['stay_id']):(hh(z['intime']),hh(z['outtime']),int(z['hadm_id'])) for z in stays}
    events={};positive_intervals=[];proc_off=[];conflicts=[]
    def add(t,s,k,origin):
        if not finite(t) or not -100000<t<=744:return
        e=events.setdefault((float(t),int(s)),dict(dev=set(),mode=set(),extub=False,et=False,newtrach=False,derived=False,origin=set()))
        e['origin'].add(origin)
        if k=='extub':e['extub']=True
        elif k=='et':e['et']=True
        elif k=='newtrach':e['newtrach']=True
        elif k=='derived':e['derived']=True
        else:e[k[0]].add(k[1])
    for z in raw:
        t=hh(z['charttime']);s=int(z['stay_id']);i=int(z['itemid']);v=str(z.get('value','')).strip()
        if i==226732:add(t,s,('dev',v),'raw_device')
        elif i in [223849,229314] and v and v!='nan':add(t,s,('mode',v),'raw_mode')
    for z in settings:
        for k in ['ventilator_mode','ventilator_mode_hamilton']:
            v=str(z.get(k,'')).strip()
            if v and v!='nan':add(hh(z['charttime']),int(z['stay_id']),('mode',v),'derived_settings_mode')
    for z in procs:
        s=int(z['stay_id']);a=hh(z['starttime']);b=hh(z['endtime']);i=int(z['itemid'])
        if z.get('statusdescription')=='Rewritten':continue
        if i in EXTUB:add(a,s,'extub','extubation_procedure');proc_off.append(a)
        elif i==224385:add(a,s,'et','intubation_procedure')
        elif i in [225448,226237]:add(a,s,'newtrach','tracheostomy_procedure')
        elif i==225792 and finite(b) and b>a and s in icu:
            a=max(a,icu[s][0]);b=min(b,icu[s][1]);
            if b>a:positive_intervals.append((a,b,s,'procedure_IMV'));add(a,s,'et','procedure_IMV_start')
    for z in derived:
        if z['ventilation_status']!='InvasiveVent':continue
        s=int(z['stay_id']);a=hh(z['starttime']);b=hh(z['endtime'])
        if s in icu and finite(b):
            a=max(a,icu[s][0]);b=min(b,icu[s][1])
            if b>a:positive_intervals.append((a,b,s,'derived_IMV'));add(a,s,'derived','derived_IMV_start')
    typed=[];routes={};trach_seen=set()
    for (t,s),e in sorted(events.items()):
        if s not in icu or not icu[s][0]-1e-8<=t<=icu[s][1]+1e-8:continue
        d=e['dev'];m=e['mode'];active=m-NIVMODES-NONACTIVE;noninv=bool(d&OFFDEV);tube='Endotracheal tube' in d;trach='Tracheostomy tube' in d
        # NIV modes or a simultaneous NIV interface identify a non-invasive
        # route; an invasive tube/IMV procedure at the same instant conflicts.
        previous,rt=routes.get(s,('unknown',-np.inf));recent=previous if t-rt<=14 else 'unknown'
        first_trach=trach and s not in trach_seen
        interface=bool(d&{'Bipap mask','CPAP mask'}) or (not d and recent=='mask')
        if e['et'] or tube:route='tube'
        elif e['newtrach'] or trach:route='trach'
        elif interface:route='mask'
        elif d&{'Trach mask','T-piece'}:route='trach_off'
        elif noninv or e['extub']:route='off'
        else:route=recent
        if trach or d&{'Trach mask','T-piece'} or e['newtrach']:trach_seen.add(s)
        ppv=bool(active or m&NIVMODES)
        via_tube=route in ['tube','trach'] and ppv and not interface
        inferred_mode=bool(active) and route=='unknown'
        weak_derived=e['derived'] and route=='unknown'
        on=tube or e['et'] or via_tube or inferred_mode or weak_derived
        off=e['extub'] or noninv or (interface and ppv) or (bool(m&NIVMODES) and route not in ['tube','trach'])
        conflict=(on and off) or (bool(active) and noninv and not interface)
        if d or e['extub'] or e['et'] or e['newtrach']:routes[s]=(route,t)
        # A new tracheostomy does not by itself terminate an ongoing,
        # source-supported IMV episode. Do not refresh the PPV timestamp.
        if (e['newtrach'] or first_trach) and not on and not off and typed:
            prev=typed[-1]
            prioroff=max([z['h'] for z in typed if z['state']=='OFF'],default=-np.inf)
            ongoing=any(a<=t<b and ss==s and a>prioroff for a,b,ss,_ in positive_intervals)
            if prev['stay']==s and prev['state']=='ON' and (t-prev['h']<=14 or ongoing):continue
        if conflict:
            state='CONFLICT';conflicts.append(dict(subject_id=pid,hour=t,stay_id=s,devices=';'.join(sorted(d)),modes=';'.join(sorted(m)),extub=e['extub'],invasive=e['et']))
        elif on:state='ON'
        elif off:state='OFF'
        elif e['newtrach'] or first_trach:state='UNKNOWN_TRACH'
        elif ppv:state='CONFLICT'
        else:continue
        typed.append(dict(h=t,stay=s,state=state,devices=';'.join(sorted(d)),modes=';'.join(sorted(m)),extub=e['extub']))
    typed.sort(key=lambda z:(z['h'],z['stay']));times=np.array([z['h'] for z in typed])
    admissions=[]
    for z in adms:
        admissions.append(dict(a=hh(z['admittime']),b=hh(z['dischtime']),dest=str(z.get('discharge_location','')).strip(),hadm=int(z['hadm_id'])))
    admissions.sort(key=lambda z:z['a'])
    exact=[]
    for key in ['all_admissions_deathtime','last_recorded_deathtime','deathtime']:
        if pd.notna(row.get(key)):exact.append(hh(row[key]))
    for z in adms:
        if pd.notna(z.get('deathtime')):exact.append(hh(z['deathtime']))
    ex=min(exact) if exact else np.nan
    dlo=hh(row.get('dod'));dhi=dlo+24 if finite(dlo) else np.nan
    def death(h):
        if finite(ex):return int(ex<=h+1e-10)
        if finite(dlo):
            if dhi<=h:return 1
            if dlo>h:return 0
            return None
        return 0
    def state(h,strict=False):
        inside=[s for s,(a,b,hadm) in icu.items() if a<=h<b]
        p=int(np.searchsorted(times,h,side='right'))-1
        last=typed[p] if p>=0 else None
        # New ICU/hospital entry cannot inherit an old ON beyond its source.
        if last:
            tag=last['state'];lt=last['h'];s=last['stay']
            if tag in ['ON','CONFLICT']:
                lastoff=max([z['h'] for z in typed[:p+1] if z['state']=='OFF'],default=-np.inf)
                coverage=any(a<=h<b and ss in inside and a>lastoff for a,b,ss,_ in positive_intervals)
                if s in inside and (h-lt<=14 or coverage):return ('ON' if tag=='ON' else 'UNKNOWN'),lt
                return 'UNKNOWN',lt
            if tag=='UNKNOWN_TRACH':return 'UNKNOWN',lt
            # OFF persists only with positive source evidence, not table silence.
            if strict and not inside:return 'UNKNOWN',lt
            # A later ICU admission permanently ends this old OFF certificate;
            # leaving that later ICU must not resurrect the former OFF state.
            if any(a>lt+1e-10 and a<=h for a,b,hadm in icu.values()):return 'UNKNOWN',lt
            current=[z for z in admissions if z['a']<=h<z['b']]
            after=[z for z in admissions if z['a']>lt and z['a']<=h]
            if after:return 'UNKNOWN',lt
            if current:return 'OFF',lt
            discharged=[z for z in admissions if z['b']<=h and z['b']>=lt]
            if discharged:
                q=max(discharged,key=lambda z:z['b'])
                if q['dest'] in SAFE:return 'OFF',lt
            return 'UNKNOWN',lt
        active=[x for x in positive_intervals if x[0]<=h<x[1] and x[2] in inside]
        return ('ON',max(x[0] for x in active)) if active else ('UNKNOWN',np.nan)
    edges={0.,24.,120.,144.,168.,192.,216.,624.,648.,672.,696.,720.}
    for z in typed:
        if 0<z['h']<720:edges.update([z['h'],min(720.,z['h']+14)])
    for a,b,s,origin in positive_intervals:
        if 0<a<720:edges.add(a)
        if 0<b<720:edges.add(b)
    for a,b,h in icu.values():
        if 0<a<720:edges.add(a)
        if 0<b<720:edges.add(b)
    for z in admissions:
        if 0<z['a']<720:edges.add(z['a'])
        if 0<z['b']<720:edges.add(z['b'])
    edges=sorted(edges);segments=[]
    for a,b in zip(edges,edges[1:]):
        ss,_=state((a+b)/2);strict,_=state((a+b)/2,True)
        segments.append(dict(subject_id=pid,start_h=a,end_h=b,state=ss,strict_state=strict))
    baseline,_=state(0);eligible=death(0)==0
    def trailing(h,key='state'):
        ss,_=state(h,key=='strict_state')
        if ss!='OFF':return np.nan
        zz=[z for z in segments if z['start_h']<h and z['end_h']>0]
        start=h
        for z in reversed(zz):
            if z[key]=='OFF':start=max(0,z['start_h'])
            else:break
        return start
    def y7(required=48,strict=False,grace=False):
        if not eligible:return np.nan,'not_both_T0_supported'
        dh=death(168);ss,_=state(168,strict);u=trailing(168,'strict_state' if strict else 'state')
        if dh==1:return 0,'death_by_day7'
        if ss=='ON':return 0,'IMV_at_day7'
        if ss=='OFF' and finite(u):
            if not grace:
                if 168-u<required-1e-8:return 0,'OFF_less_than_required_at_day7'
                if dh==0:return 1,'alive_day7_durable_OFF'
            else:
                end=max(168,u+required);dd=death(end)
                z=[z for z in segments if z['start_h']<end and z['end_h']>u]
                if dd==1 or any(x['state']=='ON' for x in z):return 0,'failed_48h_postcessation_validation'
                if dd==0 and all(x['state']=='OFF' for x in z) and end<=216:return 1,'cessation_by7_validated_48h'
        return np.nan,'unknown_survival_or_final_OFF'
    def vfd(required=48,strict=False,horizon=28):
        if not eligible:return np.nan,'not_both_T0_supported',0.,float(horizon)
        if death(horizon*24)==1:return 0.,'death_by_horizon',0.,0.
        if death(horizon*24) is None:return np.nan,'unknown_death_at_horizon',0.,float(horizon)
        ss,_=state((horizon*24),strict);u=trailing((horizon*24),'strict_state' if strict else 'state')
        if ss=='ON':return 0.,'IMV_at_day28',0.,0.
        if ss=='OFF' and finite(u):
            if (horizon*24)-u>=required-1e-8:
                value=((horizon*24)-u)/24;return value,'alive28_documented_durable_OFF',value,value
            return 0.,'late_OFF_unconfirmed_duration',0.,0.
        laston=max([min(z['end_h'],(horizon*24)) for z in segments if z['start_h']<(horizon*24) and z['state']=='ON'],default=0)
        upper=((horizon*24)-laston)/24
        if upper<required/24-1e-8:return 0.,'late_IMV_excludes_durable_OFF',0.,0.
        return np.nan,'unknown_final_IMV_or_observation',0.,upper
    y,cl=y7();ya,cla=y7(24);yb,clb=y7(grace=True);ys,cls=y7(strict=True)
    v7,c7,_,_=vfd(horizon=7)
    vf,vc,lo,hi=vfd();vf24,vc24,_,_=vfd(24);vfs,vcs,_,_=vfd(strict=True)
    result=dict(vfd7=v7,vfd7_class=c7,subject_id=pid,measurement_id=row['measurement_id'],rv_binary=row['rv_binary'],pco2=row['pco2'],pressor_T0=bool(row['t0_risk_eligible']),IMV_T0_state=baseline,both_T0_eligible=eligible,death7=death(168),death28=row['death28'],day7_state=state(168)[0],day28_state=state(672)[0],liberation7=y,liberation7_class=cl,liberation7_24h=ya,liberation7_grace=yb,liberation7_grace_class=clb,liberation7_strict=ys,liberation7_strict_class=cls,vfd28=vf,vfd28_class=vc,vfd28_lower=lo,vfd28_upper=hi,vfd28_24h=vf24,vfd28_24h_class=vc24,vfd28_strict=vfs,vfd28_strict_class=vcs,final_OFF7_start_h=trailing(168),final_OFF28_start_h=trailing(672),source_conflicts=len(conflicts),typed_points=len(typed))
    for horizon in [168,672]:result[f'unknown_hours{horizon}']=sum(min(z['end_h'],horizon)-z['start_h'] for z in segments if z['start_h']<horizon and z['state']=='UNKNOWN')
    return result,segments,conflicts,[(dict(subject_id=pid,**z)) for z in typed]


