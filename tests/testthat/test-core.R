test_that("Gauss-Hermite weights integrate a standard normal density", {
  gh <- jointBC:::.gauss_hermite(25)
  expect_equal(sum(gh$weights) / sqrt(pi), 1, tolerance = 1e-12)
  expect_equal(sum(gh$weights * (sqrt(2) * gh$nodes)^2) / sqrt(pi), 1, tolerance = 1e-12)
})

test_that("logit fit recovers a simulated signal and returns finite inference", {
  dat <- simulate_joint_bc(1800, link = "logit", seed = 20260905)
  fit <- joint_bc(y_cont ~ treatment + x, y_bin ~ treatment + x,
                  dat, link = "logit", n_quad = 25)
  truth <- c(0.5, 0.7, -0.4, -0.3, 0.8, 0.35, 0.8, 0.65)
  expect_equal(unname(coef(fit)), truth, tolerance = 0.25)
  expect_equal(fit$convergence, 0)
  expect_true(all(is.finite(vcov(fit))))
  expect_gt(joint_test(fit, "treatment", "wald")$statistic, 10)
})

test_that("probit option and likelihood-ratio omnibus test run", {
  dat <- simulate_joint_bc(1400, link = "probit", seed = 9052026)
  fit <- joint_bc(y_cont ~ treatment + x, y_bin ~ treatment + x,
                  dat, link = "probit", n_quad = 25)
  truth <- c(0.5, 0.7, -0.4, -0.3, 0.8, 0.35, 0.8, 0.65)
  expect_equal(unname(coef(fit)), truth, tolerance = 0.25)
  ans <- joint_test(fit, "treatment", "both")
  expect_equal(ans$df, c(2, 2))
  expect_true(all(ans$p_value < 0.001))
})

test_that("extreme predictors retain a finite log likelihood", {
  dat <- simulate_joint_bc(300, seed = 42)
  dat$x <- dat$x * 20
  fit <- joint_bc(y_cont ~ treatment + x, y_bin ~ treatment + x,
                  dat, n_quad = 35, control = list(maxit = 1500))
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_true(all(is.finite(coef(fit))))
})

test_that("quadrature refinement leaves the fitted result stable", {
  dat <- simulate_joint_bc(500, sigma_u = 1.1, seed = 731)
  fit_15 <- joint_bc(y_cont ~ treatment + x, y_bin ~ treatment + x,
                     dat, n_quad = 15)
  fit_35 <- joint_bc(y_cont ~ treatment + x, y_bin ~ treatment + x,
                     dat, n_quad = 35)
  expect_lt(max(abs(coef(fit_15) - coef(fit_35))), 1e-3)
  expect_lt(abs(as.numeric(logLik(fit_15) - logLik(fit_35))), 1e-3)
})
