"""Valid continuous pressor records and 60-minute merged-course membership."""
import pandas as pd
import numpy as np

PRESSORS={221906:'norepinephrine',221289:'epinephrine',229617:'epinephrine',221749:'phenylephrine',229630:'phenylephrine',229631:'phenylephrine',229632:'phenylephrine',221662:'dopamine',222315:'vasopressin',229709:'angiotensin_ii',229764:'angiotensin_ii'}
OLD_ITEMS={221906,221289,221749,221662,222315,229709,229764}
MASS={'mcg/kg/min','mcg/min','mcg/kg/hour','mcg/hour','mg/kg/min','mg/min','mg/kg/hour','mg/hour','ng/kg/min','ng/min','ng/kg/hour','ng/hour'}
def episodes(g,gap,icu_end):
 rows=[]
 for a,b in g[['start_clipped','end_clipped']].sort_values(['start_clipped','end_clipped']).itertuples(index=False,name=None):
  if not rows or a>rows[-1]['recorded_end']+pd.Timedelta(minutes=gap):rows.append({'start':a,'recorded_end':b,'records':1,'bridges':0})
  else:
   if a>rows[-1]['recorded_end']:rows[-1]['bridges']+=1
   rows[-1]['recorded_end']=max(rows[-1]['recorded_end'],b);rows[-1]['records']+=1
 for z in rows:z['operational_end']=min(z['recorded_end']+pd.Timedelta(minutes=gap),icu_end)
 return rows
def assess(row,g,gap):
 L=row.prebase;E=row.echo_time
 eps=episodes(g,gap,row.outtime)
 anchors=[z for z in eps if z['start']<L<=z['recorded_end']]
 assert len(anchors)==1,(row.subject_id,gap,len(anchors))
 ep=anchors[0];eligible=ep['start']<=E<ep['operational_end']
 # Future resumption must not be required to certify eligibility at echo.
 prefix=g.loc[g.start_clipped.le(E)]
 pp=[z for z in episodes(prefix,gap,row.outtime) if z['start']<L<=z['recorded_end']]
 assert len(pp)==1 and eligible==(pp[0]['start']<=E<pp[0]['operational_end'])
 raw=bool((g.start_clipped.le(E)&g.end_clipped.gt(E)).any())
 earlier=g.loc[g.start_clipped.le(E)];end_before=earlier.end_clipped.max()
 nexts=g.loc[g.start_clipped.gt(E),'start_clipped']
 resumed=bool(len(nexts) and nexts.min()<=end_before+pd.Timedelta(minutes=gap))
 return {'subject_id':int(row.subject_id),'grace_minutes':gap,'eligible':eligible,'raw_positive_at_echo':raw,
  'grace_at_echo':bool(eligible and not raw),'gap_with_documented_resumption':bool(eligible and not raw and resumed),
  'tail_grace_without_resumption':bool(eligible and not raw and not resumed),
  'episode_start':ep['start'],'recorded_episode_end':ep['recorded_end'],'operational_episode_end':ep['operational_end'],
  'episode_records':ep['records'],'bridged_gaps':ep['bridges'],'same_eligibility_without_future_records':True}
