# Analysis inputs

All input files are local and excluded from version control. See `extraction/` for preparation steps and limits of the extraction utilities.

## Main input: `analysis_input.csv`

One row per patient, in a fixed order. Required fields include:

| Fields | Meaning |
|---|---|
| `subject_id` | Unique patient key, retained locally to check alignment |
| `death28`, `death28_echo` | Known binary 28-day death label from the echo origin; these must agree |
| `death7_echo` | Echo-origin 7-day death label, needed for support outcomes |
| `rv_binary` | 0: no recorded global RV systolic reduction; 1: reduced |
| `pco2`, `po2` | PaCO₂ and PaO₂ of the selected arterial sample, mmHg |
| `echo_time`, `bg_time`, `t0` | Consistent timestamps; `t0` must equal `echo_time` |
| `age`, `sex` | Age in years; sex encoded as `M` or `F` |
| `height`, `weight`, `bmi` | cm, kg, kg/m²; missing values remain missing before imputation |
| `map`, `pf`, `lactate` | Mean arterial pressure (mmHg), PaO₂/FiO₂ (mmHg), lactate (mmol/L) |
| `peep`, `minute_volume` | PEEP (cmH₂O), minute ventilation (L/min) |
| `noncardiopulmonary_sofa` | Sum of coagulation, liver, CNS and renal SOFA components (0–16) |
| `vent_hours` | Index ventilation duration before L, in hours |
| `copd`, `chf`, `pulmonary_hypertension`, `cad`, `ckd` | Binary pre-index coded comorbidities |
| `prior_history_observed`, `s_doc` | Binary history-observability auxiliary and qualifying baseline pressor-record indicator |

L is the earlier of echo and blood-gas time. The exposure, RV group, death label, age, sex, prior ventilation duration and coded comorbidity indicators are observed, not imputed. The categorical input helper generates `co2lo=I(pco2<35)` and `co2hi=I(pco2>45)` deterministically; 35 and 45 belong to the middle group.

The published cohort is fixed at 2,051 patients with 998 deaths; scripts check these values. These assertions verify the published study cohort.

## Support outcome inputs

Files in `data/support/`:

```text
vpfd7_input.csv       vfd7_input.csv       crrtfd7_input.csv
vpfd28_input.csv      vfd28_input.csv      crrtfd28_input.csv
```

Each contains only patients whose endpoint is ascertainable, the main-input variables, `analysis_y` (observed free days), and `recorded_L_drug_count`. `death7_echo` or `death28_echo` must be known for the selected horizon. Unresolved support cessation/follow-up is not converted to zero. No endpoint is imputed.

## Output contracts

| Folder | Main outputs |
|---|---|
| `outputs/primary/` | 100 checkpoints, fitted primary models, coefficients and D1 tests |
| `outputs/category/` | Categorical imputations and within-RV category OR table |
| `outputs/contrasts/tables/` | `mimic_risks.csv`, `mimic_relative40.csv`, `mimic_sliding.csv`, fixed contrasts and tests |
| `outputs/sensitivity_6h/tables/` | Corresponding `mimic6h_*.csv` subset outputs |
| `outputs/support/tables/` | Six endpoint contrasts and `forest.csv` |
| `outputs/figures/` | PDF/SVG/PNG rendered from the pooled estimates |

The risk and contrast files use probability units; multiply by 100 for percentages or percentage points. Variances require multiplication by 10,000. Model objects, imputation checkpoints and cross-grid inference objects are not public data exports.
