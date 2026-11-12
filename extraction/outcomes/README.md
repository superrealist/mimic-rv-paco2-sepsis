# Organ-support outcomes

This module extracts records and reconstructs vasopressor-free, invasive-ventilation-free, and continuous renal replacement therapy (CRRT)-free days at 7 and 28 days after the index echocardiogram.

## Run order

Run from the repository root after the main cohort and baseline extraction has created `MIMIC_DATA_DIR/analysis_input.csv`:

```bash
python -m extraction.outcomes.extract_sources
python -m extraction.outcomes.build_vent
python -m extraction.outcomes.build_support
```

`MIMIC_WORKERS` controls the number of ventilation-adjudication workers (default 12). `extract_sources --only SOURCE_NAME` exports selected sources. Database access uses the shared read-only extraction helper and standard PostgreSQL environment settings.

## Inputs and source windows

The cohort input retains the original identifiers and extraction fields, including `subject_id`, `hadm_id`, `stay_id`, `measurement_id`, `echo_time`, `bg_time`, `t0_prior_later_pair`, `rv_binary`, `pco2`, `death7_echo`, `death28_echo`, `death28`, `dod`, and `s_doc_positive_drugs`. Available precise death timestamps from the cohort and subsequent admissions are used by the ventilation engine. Baseline covariates are carried through to the endpoint-specific model files.

The outcome time origin is `echo_time`. The source-query windows preserve the original extraction anchors:

- Ventilation records, ICU stays, and admissions: from 24 hours before the later paired measurement through 30 days after it; the index hospital admission is retained in full.
- Vasopressor and electronic medication administration records: the original later-pair 28-day extraction window, including overlapping ICU infusion records.
- CRRT records: from 24 hours before echocardiography through 28 days after it.

The extractor uses `mimiciv_icu` tables `chartevents`, `procedureevents`, `inputevents`, `icustays`, and `d_items`; `mimiciv_hosp` tables `admissions`, `emar`, and `emar_detail`; and the established `mimiciv_derived` tables `ventilation`, `ventilator_setting`, `crrt`, and `rrt`.

## Outcome construction

Deaths within each outcome horizon receive zero days. For survivors, terminal free days accrue from the final documented discontinuation through the horizon, with at least 48 hours of sustained liberation. The adjudicators retain uncertain support or follow-up status as missing. Ventilation states combine documented devices, ventilator modes, intubation and extubation procedures, and derived invasive-ventilation intervals. Vasopressor adjudication combines infusion intervals, documented stopping or zero-rate evidence, and subsequent electronic medication administration. CRRT adjudication combines active-system and flow evidence, procedure intervals, and documented discontinuation.

Subsequent ICU admission, hospital readmission, and discharge destination enter the ascertainment rules. Each endpoint-specific model file contains only patients with an ascertainable endpoint and retains the baseline model columns.

## Outputs

Raw source exports and patient-level adjudication files are written below `MIMIC_DATA_DIR/outcomes`. Endpoint inputs are written to:

```text
MIMIC_DATA_DIR/support/vpfd7_input.csv
MIMIC_DATA_DIR/support/vfd7_input.csv
MIMIC_DATA_DIR/support/crrtfd7_input.csv
MIMIC_DATA_DIR/support/vpfd28_input.csv
MIMIC_DATA_DIR/support/vfd28_input.csv
MIMIC_DATA_DIR/support/crrtfd28_input.csv
```

Each file includes `analysis_y`. Aggregate endpoint counts and ascertainment-reason counts are written to `MIMIC_OUTPUT_DIR/outcomes`. Runtime source exports and patient-level files belong in the configured local data directory and are excluded from the repository.
