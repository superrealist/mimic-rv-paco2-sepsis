"""All-cause death ascertainment from the recorded echo index."""
import pandas as pd
import numpy as np

def endpoint(frame,days):
    rows=[]
    for r in frame.itertuples():
        end = r.echo_time + pd.Timedelta(days=days)
        y=np.nan
        invalid=False
        reason='followup_unknown'
        if pd.notna(r.all_admissions_deathtime):
            if r.all_admissions_deathtime<=r.echo_time:
                invalid=True
                reason='precise_death_not_after_echo'
            else:
                y=float(r.all_admissions_deathtime<=end)
                reason='precise_all_admission_death_timestamp'
        elif pd.notna(r.dod):
            if r.dod.date()<=r.echo_time.date():
                invalid=True
                reason='death_date_before_or_same_as_echo'
            elif r.dod.date()==end.date():
                reason=f'day{days}_registry_date_boundary_unknown'
            else:
                y=float(r.dod.date()<end.date())
                reason='registry_date_unambiguous'
        elif pd.notna(r.last_dischtime) and r.last_dischtime+pd.Timedelta(days=365)>=end:
            y=0.
            reason='no_registry_death_within_documented_followup'
        rows.append(dict(subject_id=r.subject_id,**{f'death{days}_echo':y,
                f'invalid_echo{days}':invalid,f'echo{days}_status_reason':reason}))
    return pd.DataFrame(rows)
