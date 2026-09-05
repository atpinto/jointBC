#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
package_dir <- if (length(args)) normalizePath(args[[1]]) else normalizePath(".")
result_file <- if (length(args) >= 2L) args[[2L]] else NULL
devtools::load_all(package_dir, quiet = TRUE)

truth <- c(
  "continuous.(Intercept)" = 0.5,
  "continuous.treatment" = 0.7,
  "continuous.x" = -0.4,
  "binary.(Intercept)" = -0.3,
  "binary.treatment" = 0.8,
  "binary.x" = 0.35,
  sigma_c = 0.8,
  sigma_u = 0.65
)

run_link <- function(link, seeds = 1:8, n = 650, n_quad = 17) {
  estimates <- matrix(NA_real_, length(seeds), length(truth),
                      dimnames = list(NULL, names(truth)))
  convergence <- integer(length(seeds))
  finite_se <- logical(length(seeds))
  for (i in seq_along(seeds)) {
    dat <- simulate_joint_bc(n, link = link, seed = 10000 + seeds[[i]])
    fit <- joint_bc(y_cont ~ treatment + x, y_bin ~ treatment + x,
                    dat, link = link, n_quad = n_quad)
    estimates[i, ] <- coef(fit)
    convergence[i] <- fit$convergence
    finite_se[i] <- all(is.finite(sqrt(diag(vcov(fit)))))
  }
  data.frame(
    link = link,
    parameter = names(truth),
    truth = unname(truth),
    mean_estimate = colMeans(estimates),
    bias = colMeans(estimates) - truth,
    rmse = sqrt(colMeans((sweep(estimates, 2, truth, "-"))^2)),
    convergence_rate = mean(convergence == 0L),
    finite_se_rate = mean(finite_se),
    row.names = NULL
  )
}

logit_results <- run_link("logit")
probit_results <- run_link("probit")
results <- rbind(logit_results, probit_results)
print(results, digits = 4, row.names = FALSE)
if (!is.null(result_file)) utils::write.csv(results, result_file, row.names = FALSE)

# Deliberately large predictors stress the Bernoulli tail calculations.
stress <- simulate_joint_bc(400, seed = 991)
stress$x <- stress$x * 25
stress_fit <- joint_bc(y_cont ~ treatment + x, y_bin ~ treatment + x,
                       stress, n_quad = 35, control = list(maxit = 1500))
cat("\nStress fit: convergence=", stress_fit$convergence,
    ", finite_logLik=", is.finite(as.numeric(logLik(stress_fit))),
    ", finite_SE=", all(is.finite(sqrt(diag(vcov(stress_fit))))), "\n", sep = "")

if (any(results$convergence_rate < 0.875)) stop("Convergence rate below test threshold.")
if (any(results$finite_se_rate < 0.875)) stop("Finite-SE rate below test threshold.")
if (max(abs(results$bias[results$parameter != "sigma_u"])) > 0.25) stop("Recovery bias exceeded threshold.")
if (max(abs(results$bias[results$parameter == "sigma_u"])) > 0.30) stop("sigma_u recovery bias exceeded threshold.")
if (stress_fit$convergence != 0L || !is.finite(as.numeric(logLik(stress_fit)))) stop("Stress fit failed.")
