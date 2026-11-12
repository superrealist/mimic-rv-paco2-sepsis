# Validation scope

The public code is a reorganization of the study source. Validation completed during preparation included:

- R syntax parsing and Python syntax checks.
- Synthetic comparisons of the public and source implementations of scalar pooling, D1 joint inference, marginal-risk gradients, and risk-contrast matrices. The checked calculations agreed exactly.
- Synthetic checks of shared-reference and full-covariance propagation. A reusable patient-free check is available as `Rscript tests/synthetic_inference.R`.
- Rendering of all four single-cohort plotting interfaces from synthetic aggregate estimates.
- Reconstruction of the prepared primary cohort from existing local source extracts: all 2,051 patient identities and 276 common fields matched the frozen study input. The fixed-pair six-hour subset contained 1,812 patients. This was a local preparation comparison, not a new database extraction or refitting of patient models.
- Inspection of the publication tree for patient-data artifacts, private absolute paths, contact details, embedded credentials, and external-hospital code branches.

The public extraction and modeling entry points have not been run end to end from a newly installed MIMIC database. Database source contracts, linked echocardiography preprocessing, package versions, and local input preparation remain required.

## R packages available during the mathematical checks

| Package | Version |
|---|---|
| smcfcs | 2.0.2 |
| mice | 3.19.0 |
| jsonlite | 2.0.0 |
| digest | 0.6.39 |
| abind | 1.4.8 |
| posterior | 1.7.0 |
| sandwich | 3.1.1 |
| mitml | 0.4.5 |
| detectseparation | 0.4.0 |
| matrixStats | 1.5.0 |

These values document the check environment; they are not an environment lock or a claim that every imputation chain has been rerun in the public layout.
