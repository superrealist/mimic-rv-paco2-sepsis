"""Apply quality-screened P/F, pressor-course eligibility, and echo time origin."""
import hashlib
import json
import numpy as np
import pandas as pd
from .common import ROOT, DATA_DIR, copy, run_sql
from .pressor_course import PRESSORS, OLD_ITEMS, MASS, assess
from .mortality import endpoint


def prepare_input():
    d = pd.read_csv(ROOT / "data/preflight_yknown.csv", low_memory=False)
    assert d.subject_id.is_unique and d.death28.isin([0, 1]).all()
    d["pf_recorded_raw"] = d.pf
    d["pf"] = d.pf_2h_qc_sensitivity
    assert d.loc[d.pf_pair_incompatible, "pf"].isna().all()
    d.to_csv(ROOT / "data/paired_analysis_input.csv", index=False)
    d.loc[d.s_doc.eq(1)].to_csv(ROOT / "data/pressor_baseline_population.csv", index=False)


def extract_pressor_intervals():
    q = pd.read_csv(ROOT / "data/pressor_baseline_population.csv", low_memory=False)
    assert len(q) > 0 and q.subject_id.is_unique
    values = ",".join(
        f"({int(r.subject_id)},{int(r.hadm_id)},{int(r.stay_id)},"
        f"timestamp '{pd.Timestamp(r.prebase)}',timestamp '{pd.Timestamp(r.t0)}')"
        for r in q.itertuples(index=False)
    )
    query = f"""WITH c(subject_id,hadm_id,index_stay_id,l,t0) AS (VALUES {values})
    SELECT c.subject_id cohort_subject_id,c.index_stay_id,z.*,
           s.intime event_icu_intime,s.outtime event_icu_outtime
    FROM c JOIN mimiciv_icu.inputevents z
      ON z.subject_id=c.subject_id AND z.hadm_id=c.hadm_id
    JOIN mimiciv_icu.icustays s ON s.stay_id=z.stay_id
      AND s.subject_id=c.subject_id AND s.hadm_id=c.hadm_id
    WHERE z.endtime>=c.l-interval '6 hour'
      AND z.starttime<=c.t0+interval '54 hour'
      AND z.itemid IN(221289,229617,221749,229630,229631,229632,
                     221906,221662,222315,229709,229764,221653,221986)"""
    run_sql("pressor_course_sources", copy(query, "data/pressor_intervals_raw.csv"))


def build_final():
    a = pd.read_csv(ROOT / "data/pressor_baseline_population.csv", low_memory=False)
    for c in ("t0", "prebase", "echo_time", "bg_time", "intime", "outtime"):
        a[c] = pd.to_datetime(a[c], format="mixed")
    assert a.subject_id.is_unique
    assert a.t0.eq(a[["echo_time", "bg_time"]].max(axis=1)).all()
    assert a.prebase.eq(a[["echo_time", "bg_time"]].min(axis=1)).all()
    p = pd.read_csv(ROOT / "data/pressor_intervals_raw.csv", low_memory=False)
    for c in ("starttime", "endtime", "event_icu_intime", "event_icu_outtime"):
        p[c] = pd.to_datetime(p[c], format="mixed")
    p = p.merge(a[["subject_id", "stay_id", "prebase", "echo_time", "outtime"]].rename(
        columns={"subject_id": "fixed_subject", "stay_id": "fixed_stay"}),
        left_on="cohort_subject_id", right_on="fixed_subject", validate="many_to_one")
    p["drug"] = p.itemid.map(PRESSORS)
    p["unit"] = p.rateuom.fillna("").str.strip().str.lower()
    p["status"] = p.statusdescription.fillna("").str.strip().str.lower()
    known = p.unit.isin(MASS)
    known.loc[p.drug.eq("vasopressin")] = p.loc[p.drug.eq("vasopressin"), "unit"].isin(["units/min", "units/hour"])
    valid = (p.drug.notna() & p.rate.gt(0) & np.isfinite(p.rate) & known
        & ~p.status.isin(["rewritten", "cancelled", "canceled", "bolus"])
        & p.endtime.gt(p.starttime)
        & ~p.ordercategorydescription.fillna("").str.contains("bolus|drugpush", case=False, regex=True)
        & p.stay_id.eq(p.fixed_stay))
    witness = valid & p.itemid.isin(OLD_ITEMS) & p.starttime.lt(p.prebase) & p.endtime.ge(p.prebase)
    assert set(p.loc[witness, "cohort_subject_id"]) == set(a.subject_id)
    p = p.loc[valid].copy()
    p["start_clipped"] = p[["starttime", "event_icu_intime"]].max(axis=1)
    p["end_clipped"] = p[["endtime", "event_icu_outtime"]].min(axis=1)
    p = p.loc[p.end_clipped.gt(p.start_clipped)]
    p.to_csv(ROOT / "data/valid_pressor_source_intervals.csv", index=False)
    groups = {int(k): g for k, g in p.groupby("cohort_subject_id")}
    audit = pd.DataFrame([assess(r, groups[int(r.subject_id)], 60) for r in a.itertuples()])
    audit.to_csv(ROOT / "data/echo_course_audit.csv", index=False)
    eligible = set(audit.loc[audit.eligible, "subject_id"])
    d = a.loc[a.subject_id.isin(eligible)].copy().sort_values("subject_id")
    d["abs_pair_gap_h"] = (d.bg_time - d.echo_time).abs().dt.total_seconds() / 3600
    assert d.abs_pair_gap_h.le(12).all()
    for c in ("sepsis_complete_upper", "dod", "last_dischtime", "all_admissions_deathtime"):
        d[c] = pd.to_datetime(d[c], format="mixed")
    y28 = endpoint(d, 28)
    y7 = endpoint(d, 7)
    # This is a cohort-preserving re-anchoring check used by the final analysis.
    if y28.invalid_echo28.any() or not np.array_equal(y28.death28_echo.to_numpy(), d.death28.to_numpy()):
        raise ValueError("Echo-origin death ascertainment differs from the fixed paired cohort; review source version before analysis")
    d["t0_prior_later_pair"] = d.t0
    d["pair_end"] = d.t0
    d["t0"] = d.echo_time
    d["sepsis_complete_after_echo"] = d.sepsis_complete_upper > d.echo_time
    d = d.merge(y28, on="subject_id", validate="1:1").merge(y7, on="subject_id", validate="1:1")
    for c in ("map_time", "lactate_time", "peep_time", "minute_volume_time", "height_time", "weight_time", "sofa_time_used"):
        assert not (pd.to_datetime(d[c], format="mixed") > d.echo_time).any()
    target = DATA_DIR / "analysis_input.csv"
    d.to_csv(target, index=False)
    d.loc[d.abs_pair_gap_h.le(6)].to_csv(DATA_DIR / "analysis_input_6h.csv", index=False)
    summary = {"main_n": len(d), "sensitivity_6h_n": int(d.abs_pair_gap_h.le(6).sum()),
               "rv_reduced_n": int(d.rv_binary.sum()), "death28_n": int(d.death28.sum()),
               "analysis_input_sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
               "T0": "echo_time", "L": "minimum of echo_time and bg_time"}
    (ROOT / "logs/final_cohort.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    prepare_input()
    extract_pressor_intervals()
    build_final()
