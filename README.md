# MIMIC-IV: RV function and PaCO₂-associated risk

Code supporting the MIMIC-IV component of:

**Heterogeneity in PaCO₂ Associations by Right Ventricular Function in Mechanically Ventilated Patients with Vasopressor-Treated Sepsis**

Repository: <https://github.com/superrealist/mimic-rv-paco2-sepsis>

The study evaluates the association between paired arterial PaCO₂ and 28-day mortality across qualitative right ventricular systolic-function groups in adults with sepsis receiving invasive ventilation and vasopressor support. The index time is echocardiography. The main MIMIC cohort contains 2,051 patients; the fixed-pair ±6-hour sensitivity subset contains 1,812.

This repository contains MIMIC-specific extraction, analysis, inference, and rendering code. External-hospital data and models are not included. No patient-level data, imputed datasets, fitted objects, or database credentials are distributed.

## Layout

```text
extraction/             Credentialed extraction and preparation utilities
R/config.R              Data/output paths, fixed knots and model specifications
R/functions/            Imputation, pooling and full-covariance risk contrasts
R/01–09_*.R             Analysis entry points
plots/                  Figures from saved pooled estimates
docs/                   Input contracts and methodological details
data/                   Local inputs; ignored by version control
outputs/                Local results, fits and diagnostics; ignored by version control
```

## Requirements and paths

The analysis was developed with R 4.5.2 on Windows. Required R packages are checked by `Rscript R/check_dependencies.R`. The numerical imputation wrappers inspect the installed `smcfcs` function body and include structural assertions; compatibility with a different package version must be checked before analysis. Python dependencies are listed in `requirements.txt`. This is not a complete environment lock.

Run commands from the repository root. Environment variables configure local paths:

| Variable | Default | Purpose |
|---|---|---|
| `MIMIC_CODE_ROOT` | Current repository directory | R source location |
| `MIMIC_DATA_DIR` | `data/` | Locally prepared inputs and extraction intermediates |
| `MIMIC_OUTPUT_DIR` | `outputs/` | Estimates, fitted objects and figures |
| `MIMIC_ANALYSIS_INPUT` | `data/analysis_input.csv` under the selected data directory | Final echo-anchored analysis input |
| `MIMIC_WORKERS` | `4` | Parallel R workers |

Use operating-system environment settings to provide these values. Keep inputs and outputs in ignored directories or outside this repository. Database settings and extraction order are described under `extraction/`.

With the documented database and linked echocardiography prerequisites in place, the extraction entry points are:

```sh
python -m extraction.run_primary
python -m extraction.outcomes.extract_sources
python -m extraction.outcomes.build_vent
python -m extraction.outcomes.build_support
```

The primary extraction writes `analysis_input.csv` and `analysis_input_6h.csv`; the support steps write the six endpoint inputs. The linked echocardiography preprocessing contract is documented in `extraction/README.md`.

## Analysis order

Prepare inputs that satisfy the cohort and variable definitions in [docs/input_contract.md](docs/input_contract.md) before running patient analyses.

```sh
Rscript R/check_dependencies.R
Rscript R/01_primary_mi.R
Rscript R/03_categorical_mi.R
Rscript R/04_categorical_models.R
Rscript R/05_risk_curves.R
Rscript R/05b_curve_minimum.R
Rscript R/06_simultaneous_bands.R
Rscript R/07_pairing_sensitivity.R
Rscript R/08_support_mi.R
Rscript R/09_support_models.R
python plots/render_figures.py all
```

`01_primary_mi.R` creates 100 independently seeded completed datasets and then invokes `02_primary_models.R`. The latter can also be run separately against complete, matching checkpoints.

Categorical mortality inference uses separately generated categorical-model-compatible imputations. If exactly one categorical chain is missing after a numerical convergence failure, `03b_category_retry.R` reproduces the study's deterministic retry and assembly procedure using the same seed. This conditional recovery script is not an unconditional step after a successful 100-chain run.

The six support outcomes require their separate input files. `08_support_mi.R` processes all six by default; an endpoint name can select an individual outcome. Each uses 40 completed datasets. Outcome ascertainment is performed before covariate imputation.

## Methods retained

- Fixed three-knot PaCO₂ RCS at 30, 40 and 55 mmHg, with RV-by-spline interactions.
- Adjustment for age, sex, BMI, five recorded comorbidities, mean arterial pressure, PaO₂/FiO₂, lactate, noncardiopulmonary SOFA, prior index-episode ventilation duration, PEEP and minute ventilation. Fixed functional forms are defined in `R/config.R` and the input helpers.
- Bounded substantive-model-compatible imputation for mortality analyses. The prior-history observability indicator is an auxiliary in the continuous primary model. Categorical imputations contain that indicator in their substantive model; final categorical regressions omit it.
- Rubin pooling, Barnard–Rubin scalar degrees of freedom, and D1 joint tests. Group-specific nonlinearity and global interaction are evaluated on the log-odds scale.
- Cohort-standardized risks, differences relative to 40 mmHg, and changes over a continuous +10-mmHg interval. Delta-method variances retain the covariance among predictions and shared reference values.
- Pointwise intervals and an approximate simultaneous band over 301 interaction contrasts on a fixed PaCO₂ grid, obtained from 200,000 correlated Gaussian simulations using the pooled covariance.
- The ±6-hour analysis restricts the original pair and refits the model in the corresponding subsets of the 100 primary completed datasets. It does not select a new blood gas or create new imputations.
- Endpoint-specific outcome-enriched PMM/FCS for support outcomes and linear models with HC3 covariance. Death within the horizon gives zero free days; free days follow documented final discontinuation sustained for at least 48 hours.

See [docs/methods.md](docs/methods.md) for estimands and units. Figure layouts are adapted for a single MIMIC cohort; they do not reproduce the two-cohort manuscript layout.

## Reproducibility scope

This public version reorganizes source scripts, removes external-cohort branches, and replaces private workspace paths with configurable paths. Syntax and synthetic mathematical checks are recorded in [docs/validation.md](docs/validation.md). The public entry points have not been executed end to end on patient data. Authorized MIMIC access, the linked echocardiography source, the documented derived tables, and complete input preparation are required to reproduce the patient analyses.

Checkpoint and model files contain patient-derived information and remain local. Do not publish them. Figures and aggregate tables should be reviewed for the applicable data-use requirements before sharing.

No license file is included. Third-party dependencies retain their respective licenses.

## 中文说明

本仓库公开 MIMIC 部分的提取、分析及绘图代码，不含本地医院队列、患者数据或数据库凭据。路径通过环境变量设置。连续模型、分类模型和器官支持结局采用各自的插补流程；6小时敏感性分析在主模型完成数据的固定配对子集中重拟合。公开版已进行语法及合成数学检查，未重新执行患者数据的完整分析。
