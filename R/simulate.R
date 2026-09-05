#' Simulate from the joint binary-continuous model
#'
#' Generates an intercept, a binary treatment, and one standard-normal
#' covariate. The coefficient vectors must correspond to those three columns.
#'
#' @param n Number of subjects.
#' @param beta_continuous Length-three continuous-model coefficient vector.
#' @param beta_binary Length-three binary-model coefficient vector.
#' @param sigma_c Positive continuous residual scale.
#' @param sigma_u Positive latent-effect standard deviation.
#' @param link Binary link.
#' @param seed Optional random seed.
#' @return A data frame with observed outcomes, covariates, and the simulated
#'   latent effect in column `u`.
#' @export
simulate_joint_bc <- function(n, beta_continuous = c(0.5, 0.7, -0.4),
                              beta_binary = c(-0.3, 0.8, 0.35),
                              sigma_c = 0.8, sigma_u = 0.65,
                              link = c("logit", "probit"), seed = NULL) {
  link <- match.arg(link)
  if (!is.null(seed)) set.seed(seed)
  if (length(n) != 1L || n < 1 || n != as.integer(n)) stop("`n` must be a positive integer.", call. = FALSE)
  if (length(beta_continuous) != 3L || length(beta_binary) != 3L) {
    stop("Both coefficient vectors must have length three: intercept, treatment, x.", call. = FALSE)
  }
  if (!is.finite(sigma_c) || !is.finite(sigma_u) || sigma_c <= 0 || sigma_u <= 0) {
    stop("`sigma_c` and `sigma_u` must be positive.", call. = FALSE)
  }
  treatment <- stats::rbinom(n, 1, 0.5)
  x <- stats::rnorm(n)
  X <- cbind(1, treatment, x)
  u <- stats::rnorm(n, sd = sigma_u)
  y_cont <- drop(X %*% beta_continuous) + sigma_c * u + stats::rnorm(n, sd = sigma_c)
  eta_bin <- drop(X %*% beta_binary) + u
  p <- if (link == "logit") stats::plogis(eta_bin) else stats::pnorm(eta_bin)
  y_bin <- stats::rbinom(n, 1, p)
  data.frame(y_cont = y_cont, y_bin = y_bin, treatment = treatment, x = x, u = u)
}

