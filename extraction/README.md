# MIMIC extraction

This directory contains the SQL query builders and record-level preparation used
to construct the MIMIC analysis input. It performs no model fitting. The database
is queried in read-only mode; all exports stay on the authorized user's machine.

## Required source tables

Use an authorized PostgreSQL installation with these schemas:

- `mimiciv_hosp`: `patients`, `admissions`, `labevents`, `d_labitems`,
  `microbiologyevents`, `diagnoses_icd`.
- `mimiciv_icu`: `icustays`, `inputevents`, `chartevents`, `d_items`.
- `mimiciv_derived`: `sepsis3`, `suspicion_of_infection`, `ventilation`, `bg`,
  `sofa`, `vitalsign`, `ventilator_setting`, `norepinephrine_equivalent_dose`,
  `chemistry`, `height`, `complete_blood_count`, `enzyme`, `gcs`,
  `urine_output_rate`.
- `mimiciv_echo`: the two study-level standardized tables below.

The derived tables follow the PostgreSQL concepts in the
[MIMIC code repository](https://github.com/MIT-LCP/mimic-code/tree/main/mimic-iv/concepts).
They must be built before this pipeline. Record the MIMIC release and the concept
code revision locally. The repository does not distribute source records or
rebuild the official concepts.

### Echo input contract

The source environment used a normalized echo import, not the raw download layout.
Its original import implementation is not included. Preparing the following
tables from authorized echo files is an explicit upstream dependency. In
particular, study identifiers and timestamp precision must come from that import;
this pipeline does not infer them from patient identifiers or assign a new clock
time to date-only records.

| Table | Required fields | Contract |
|---|---|---|
| `measurement_events` | `measurement_id`, `subject_id`, `test_type`, `measurement_datetime` | One row per study; stable integer study ID; `test_type='tte'`; source-recorded diagnostic timestamp. |
| `structured_measurement` | `measurement_id`, `subject_id`, `test_type`, `measurement_datetime`, `measurement`, `result` | Each measurement links to the same study ID and timestamp; values retain their source labels. |

`measurement_id` is used for deterministic tie resolution and joins. Its mapping
must remain stable. `validate_echo_inputs.py` can validate CSV exports of these
two prepared tables before database loading; it does not perform raw-file parsing
or create study identifiers.

For `measurement='rv_function'`, the recognized labels are:

- No recorded reduction: `Nl RV function`, `Low normal function`, `Hyperdynamic`.
- Reduced: `Mild global RV hypo`, `Moderate global RV hypo`,
  `Severe global hypo`, `RV function depressed`.

Conflicting or other labels remain unclassified. Other extracted echo variables
are `tapse`, `rv_chamber_size`, `rv_diam`, `tr_velocity`, and `tr_mmhg`.

## Configuration and execution

Install Python dependencies from the repository requirements and install the
PostgreSQL `psql` client. Use standard libpq settings (`PGHOST`, `PGPORT`,
`PGDATABASE`, `PGUSER`) and a local password file or `PGPASSWORD` environment
variable. Never put credentials in a command, source file, or committed config.
`MIMIC_PSQL` optionally specifies the `psql` executable.

`MIMIC_DATA_DIR` defaults to `data/` at the repository root and
`MIMIC_OUTPUT_DIR` defaults to `outputs/`. Keep either directory private even when
using a location outside the repository.

From the repository root:

```text
python -m extraction.run_primary
```

The stages can also be invoked individually with `--stage core`,
`--stage baseline`, and `--stage course`, in that order. Each depends on the
preceding stage's files. Generated SQL, patient-level intermediate CSVs, and logs
are written beneath `MIMIC_DATA_DIR/stages/primary/`; they are not public artifacts.
The pipeline exports `MIMIC_DATA_DIR/analysis_input.csv` and
`MIMIC_DATA_DIR/analysis_input_6h.csv`. See `outcomes/README.md` for organ-support
outcomes and the six additional analysis inputs.

## Definitions implemented in code

1. Adult Sepsis-3 ICU stays are linked to transthoracic echo examinations during
   an invasive ventilation interval. Microbiology source records establish exact
   culture times or date-only intervals; the upper bound is used for sepsis
   ascertainment at the later paired measurement.
2. Arterial gases require an unambiguous `ART` specimen and consistent PaCO2
   measurements in mmHg. The nearest eligible gas within 12 hours on either side
   of echo is selected in the same ICU and ventilation interval. An equidistant
   tie selects the earlier gas. The first eligible echo per patient is retained.
   Source-recorded midnight echo timestamps are retained and flagged.
3. Baseline reference `L=min(echo_time,bg_time)`. MAP uses the most recent valid
   observation before L within 2 hours. Lactate uses 12 hours; PEEP and minute
   ventilation use 6 hours within the same ventilation interval. SOFA uses the
   latest available rolling 24-hour assessment ending in the 24 hours before L.
   The noncardiopulmonary score sums coagulation, liver, CNS, and renal components.
   Full SOFA and component availability are retained separately.
4. P/F uses the arterial sample FiO2 when available, otherwise the most recent
   setting in the preceding 2 hours in the same ventilation interval. Values
   flagged by `PaO2 > FiO2_fraction*(760-47)` are missing in the analysis P/F field;
   the recorded ratio is retained as `pf_recorded_raw`. Height and weight precede
   L. The five history indicators use coded diagnoses from admissions discharged
   before L. Age is capped at 90 years.
5. A valid positive continuous pressor record must start before L and reach L.
   Valid records in the index ICU are combined across pressor drugs. Gaps up to
   60 minutes are merged and the terminal interval receives 60 minutes of grace,
   bounded by ICU discharge. Echo must lie in that same treatment course. A
   built-in check verifies that course membership at echo does not require
   future-starting records.
6. Final `t0=echo_time`. During the initial fixed-pair eligibility and baseline
   extraction stages, the internal `t0` field denotes the later paired measurement.
   `finalize_cohort.py` preserves it as `pair_end` and `t0_prior_later_pair`, then
   sets the analysis time origin to echo. It recomputes 7- and 28-day mortality
   and verifies unchanged 28-day classifications for the retained cohort.
7. The 6-hour cohort restricts the original selected pair to an absolute gap of
   at most 6 hours; it does not select another echo or gas.

The saved study cohort contains 2,051 patients, including 779 with reduced RV
systolic function and 998 deaths by day 28. The 6-hour subset contains 1,812.
These are reference counts for the documented source version, not substitutes
for inspecting the extraction logs when using another database build.

## Scope of reproducibility

The published code exposes the cohort, pairing, covariate, pressor-course, death,
and support-outcome algorithms. Running it requires authorized source data,
official derived concepts, and the standardized echo import above. No database
extraction was rerun while preparing this public package. Local validation of
the packaged transformations is recorded outside the public repository.
