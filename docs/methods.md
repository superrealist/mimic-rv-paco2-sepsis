# Estimands and computation

For RV group g, `R_g(x)` is the predicted 28-day mortality probability averaged over every patient's completed covariates after assigning RV group g and PaCO₂ x to all patients. Both RV curves use the same cohort-wide reference population.

- Relative-to-40 contrast: `R_g(x) - R_g(40)`.
- Continuous 10-mmHg change: `D_g(x) = R_g(x + 10) - R_g(x)`.
- Additive interaction contrast: `I(x) = D_1(x) - D_0(x)`.
- Fixed contrasts: 30→40, 40→50, 50→60 and 40→60 mmHg.

These are cross-sectional exposure contrasts, not observed within-patient trajectories or intervention effects. The fitted curve minimum is a descriptive model estimate within 25–65 mmHg, not an ROC cutoff or treatment target.

For each imputation, a combined contrast gradient is multiplied by the full coefficient covariance matrix. Rubin pooling incorporates within- and between-imputation variance. The simultaneous band uses the full cross-grid total covariance for the 301 starting values from 25 to 55 mmHg. Its coverage is approximate, asymptotic and specific to that grid and cohort. It is not a joint confidence statement across all endpoints or a guarantee between grid points.

The scalar intervals condition on each completed empirical covariate distribution. That distribution is not bootstrapped. Global and nonlinear interaction P values test logistic coefficients and therefore are not the same null hypothesis as an additive risk-change contrast.

## Statistical references

- Bartlett JW et al. Substantive-model-compatible fully conditional specification. *Stat Methods Med Res*. 2015;24:462–487. [doi:10.1177/0962280214521348](https://doi.org/10.1177/0962280214521348).
- Barnard J, Rubin DB. Small-sample degrees of freedom with multiple imputation. *Biometrika*. 1999;86:948–955. [doi:10.1093/biomet/86.4.948](https://doi.org/10.1093/biomet/86.4.948).
- Reiter JP. Multi-component significance tests with multiply imputed data. *Biometrika*. 2007;94:502–508. [doi:10.1093/biomet/asm028](https://doi.org/10.1093/biomet/asm028).
- Muller CJ, MacLehose RF. Predicted probabilities and target populations. *Int J Epidemiol*. 2014;43:962–970. [doi:10.1093/ije/dyu029](https://doi.org/10.1093/ije/dyu029).
