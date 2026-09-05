# jointBC prototype

This package fits a cross-sectional joint model for one Gaussian continuous
outcome and one binary outcome. The two outcomes may have different formulas.
Dependence is induced by a shared normal latent effect, and the subject-level
likelihood is integrated by Gauss-Hermite quadrature with log-sum-exp
stabilization.

## Model used in this prototype

For subject \(i\),

\[
Y_{ci}\mid u_i \sim N(x_{ci}'\beta_c + \sigma_c u_i,\;\sigma_c^2),
\qquad
P(Y_{bi}=1\mid u_i)=G(x_{bi}'\beta_b+u_i),
\qquad
u_i\sim N(0,\sigma_u^2).
\]

`G` is logistic by default and can be changed to the standard normal CDF.
Both positive scale parameters are optimized on the log scale. Standard errors
come from the observed Hessian and are transformed to the natural scale with
the delta method.

**Parameterization status:** the scaling above agrees with Appendix 2 of the
[open-access Teixeira-Pinto and Normand paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC2818753/):
the binary predictor contains `+ u`, while the continuous likelihood contains
`u * sigma2` and has residual SD `sigma2`. The SAS statement
`random u ~ normal(0, sigmab)` uses `sigmab` as a variance. This prototype exposes
the equivalent `sigma_u = sqrt(sigmab)` as a standard deviation, so the probit
latent correlation is `sigma_u^2 / (1 + sigma_u^2)`.

**Still to verify before freezing the public API:** decide whether the exported
association parameter should remain the more explicit SD `sigma_u` or instead
use the paper's variance-scale `sigmab`. The logit model is an intentional
generalization, not the paper's original probit model. Its reported correlation
uses logistic latent residual variance `pi^2 / 3` and is not an observed
binary-continuous Pearson correlation. The fixed positive loading also directly
represents positive association; negative association requires reversing one
outcome's coding/sign, as in the paper's discussion.

## Quick start

```r
devtools::load_all("path/to/jointBC")

dat <- simulate_joint_bc(1000, link = "logit", seed = 1)
fit <- joint_bc(
  y_cont ~ treatment + x,
  y_bin ~ treatment,
  data = dat,
  link = "logit",
  n_quad = 25
)

summary(fit)
joint_test(fit, "treatment", method = "both")
```

The omnibus test jointly tests every coefficient produced by each requested
formula term in both submodels. `method="wald"` uses the full joint covariance
matrix. `method="lrt"` removes the terms from both formulas and refits.

## Run checks

```sh
R CMD check jointBC
Rscript jointBC/inst/simulation/run_simulation_tests.R jointBC simulation-results.csv
```
