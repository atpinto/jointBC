#' Fit a joint binary-continuous model
#'
#' Fits one Gaussian continuous outcome and one binary outcome with separate
#' formulas and a shared normal latent effect. The implemented conditional
#' model is
#' \deqn{Y_{ci}|u_i \sim N(x_{ci}'\beta_c + \sigma_c u_i, \sigma_c^2),}
#' \deqn{P(Y_{bi}=1|u_i)=G(x_{bi}'\beta_b+u_i),\quad
#' u_i\sim N(0,\sigma_u^2),}
#' where \eqn{G} is logistic by default or standard normal for `link="probit"`.
#' The subject likelihood is integrated with non-adaptive Gauss-Hermite
#' quadrature and accumulated on the log scale.
#'
#' @param formula_continuous Formula for the Gaussian outcome.
#' @param formula_binary Formula for the binary outcome.
#' @param data Data frame containing both outcomes and all covariates.
#' @param link Binary link, either `"logit"` (default) or `"probit"`.
#' @param n_quad Number of Gauss-Hermite nodes.
#' @param start Optional natural-scale starting vector in the same order as
#'   `coef()`; it may be named.
#' @param control Control list passed to [stats::optim()].
#' @param keep_data Keep the model data so likelihood-ratio joint tests can
#'   refit a reduced model.
#' @return An object of class `joint_bc`.
#' @export
joint_bc <- function(formula_continuous, formula_binary, data,
                     link = c("logit", "probit"), n_quad = 25L,
                     start = NULL, control = list(), keep_data = TRUE) {
  link <- match.arg(link)
  call <- match.call()
  if (!is.data.frame(data)) data <- as.data.frame(data)

  mf_c <- stats::model.frame(
    formula_continuous, data = data, na.action = stats::na.fail,
    drop.unused.levels = TRUE
  )
  mf_b <- stats::model.frame(
    formula_binary, data = data, na.action = stats::na.fail,
    drop.unused.levels = TRUE
  )
  if (nrow(mf_c) != nrow(mf_b) || !identical(row.names(mf_c), row.names(mf_b))) {
    stop("The two formulas must evaluate on the same subjects in the same order.", call. = FALSE)
  }

  terms_c <- stats::terms(mf_c)
  terms_b <- stats::terms(mf_b)
  y_c <- stats::model.response(mf_c)
  y_b_raw <- stats::model.response(mf_b)
  if (!is.numeric(y_c)) stop("The continuous response must be numeric.", call. = FALSE)
  if (is.factor(y_b_raw)) {
    if (nlevels(y_b_raw) != 2L) stop("The binary response factor must have two levels.", call. = FALSE)
    y_b <- as.integer(y_b_raw) - 1L
  } else if (is.logical(y_b_raw)) {
    y_b <- as.integer(y_b_raw)
  } else {
    y_b <- as.numeric(y_b_raw)
  }
  if (any(!is.finite(y_c)) || any(!y_b %in% c(0, 1))) {
    stop("Responses must be finite, and the binary response must be coded 0/1 or be a two-level factor.", call. = FALSE)
  }

  X_c <- stats::model.matrix(terms_c, mf_c)
  X_b <- stats::model.matrix(terms_b, mf_b)
  p_c <- ncol(X_c)
  p_b <- ncol(X_b)
  gh <- .gauss_hermite(n_quad)
  log_weights <- log(gh$weights) - 0.5 * log(pi)

  natural_names <- c(
    paste0("continuous.", colnames(X_c)),
    paste0("binary.", colnames(X_b)),
    "sigma_c", "sigma_u"
  )

  if (is.null(start)) {
    lm_fit <- stats::lm.fit(X_c, y_c)
    beta_c0 <- lm_fit$coefficients
    beta_c0[!is.finite(beta_c0)] <- 0
    sigma_c0 <- sqrt(sum(lm_fit$residuals^2) / max(1, nrow(X_c) - p_c))
    if (!is.finite(sigma_c0) || sigma_c0 <= 0) sigma_c0 <- stats::sd(y_c)
    if (!is.finite(sigma_c0) || sigma_c0 <= 0) sigma_c0 <- 1
    glm_fit <- suppressWarnings(stats::glm.fit(X_b, y_b, family = stats::binomial(link = link)))
    beta_b0 <- glm_fit$coefficients
    beta_b0[!is.finite(beta_b0)] <- 0
    start_natural <- c(beta_c0, beta_b0, sigma_c0, 0.35)
  } else {
    if (!is.numeric(start) || length(start) != length(natural_names)) {
      stop("`start` must be a numeric vector with one value per model parameter.", call. = FALSE)
    }
    if (!is.null(names(start))) {
      missing_names <- setdiff(natural_names, names(start))
      if (length(missing_names)) stop("Named `start` is missing: ", paste(missing_names, collapse = ", "), call. = FALSE)
      start <- start[natural_names]
    }
    start_natural <- unname(start)
  }
  if (any(start_natural[c(p_c + p_b + 1L, p_c + p_b + 2L)] <= 0)) {
    stop("Starting values for `sigma_c` and `sigma_u` must be positive.", call. = FALSE)
  }

  to_working <- function(theta) {
    c(theta[seq_len(p_c + p_b)], log(theta[p_c + p_b + 1L]), log(theta[p_c + p_b + 2L]))
  }
  from_working <- function(par) {
    c(par[seq_len(p_c + p_b)], exp(par[p_c + p_b + 1L]), exp(par[p_c + p_b + 2L]))
  }

  objective <- function(par) {
    beta_c <- par[seq_len(p_c)]
    beta_b <- par[p_c + seq_len(p_b)]
    sigma_c <- exp(par[p_c + p_b + 1L])
    sigma_u <- exp(par[p_c + p_b + 2L])
    u_nodes <- sqrt(2) * sigma_u * gh$nodes

    eta_c <- outer(drop(X_c %*% beta_c), sigma_c * u_nodes, "+")
    eta_b <- outer(drop(X_b %*% beta_b), u_nodes, "+")
    ll_c <- matrix(
      stats::dnorm(rep(y_c, times = length(u_nodes)),
                   mean = as.vector(eta_c), sd = sigma_c, log = TRUE),
      nrow = length(y_c)
    )
    ll_b <- .binary_loglik(y_b, eta_b, link)
    log_integrand <- sweep(ll_c + ll_b, 2L, log_weights, "+")
    value <- -.row_log_sum_exp(log_integrand)
    total <- sum(value)
    if (!is.finite(total)) .Machine$double.xmax^(1 / 4) else total
  }

  defaults <- list(maxit = 1000L, reltol = 1e-9)
  control <- utils::modifyList(defaults, control)
  opt <- stats::optim(
    par = to_working(start_natural), fn = objective, method = "BFGS",
    control = control, hessian = FALSE
  )
  hessian <- stats::optimHess(opt$par, objective)
  vcov_working <- .invert_hessian(hessian)
  estimate <- from_working(opt$par)
  names(estimate) <- natural_names

  jacobian <- diag(length(estimate))
  jacobian[p_c + p_b + 1L, p_c + p_b + 1L] <- estimate[p_c + p_b + 1L]
  jacobian[p_c + p_b + 2L, p_c + p_b + 2L] <- estimate[p_c + p_b + 2L]
  vcov_natural <- jacobian %*% vcov_working %*% jacobian
  dimnames(vcov_natural) <- list(natural_names, natural_names)
  dimnames(vcov_working) <- list(natural_names, natural_names)

  structure(list(
    call = call,
    formula_continuous = formula_continuous,
    formula_binary = formula_binary,
    terms_continuous = terms_c,
    terms_binary = terms_b,
    link = link,
    n_quad = as.integer(n_quad),
    nobs = length(y_c),
    coefficients = estimate,
    coefficients_working = stats::setNames(opt$par, natural_names),
    vcov = vcov_natural,
    vcov_working = vcov_working,
    hessian = hessian,
    logLik = -opt$value,
    convergence = opt$convergence,
    message = opt$message,
    counts = opt$counts,
    X_continuous = X_c,
    X_binary = X_b,
    y_continuous = y_c,
    y_binary = y_b,
    data = if (isTRUE(keep_data)) data else NULL,
    keep_data = isTRUE(keep_data),
    control = control
  ), class = "joint_bc")
}

#' Latent-scale correlation implied by a fitted model
#'
#' For probit this is the exact correlation between the latent continuous and
#' binary variables under the implemented parameterization. For logit it uses
#' the conventional logistic latent residual variance, \eqn{\pi^2/3}, and is
#' therefore a latent-variable interpretation rather than an observed Pearson
#' correlation.
#'
#' @param object A fitted `joint_bc` model, or a positive `sigma_u` value.
#' @param link Needed only when `object` is numeric.
#' @return The implied latent-scale correlation.
#' @export
latent_correlation <- function(object, link = c("logit", "probit")) {
  if (inherits(object, "joint_bc")) {
    sigma_u <- unname(object$coefficients["sigma_u"])
    link <- object$link
  } else {
    link <- match.arg(link)
    sigma_u <- as.numeric(object)
    if (length(sigma_u) != 1L || !is.finite(sigma_u) || sigma_u < 0) {
      stop("`object` must be a fitted model or one non-negative `sigma_u`.", call. = FALSE)
    }
  }
  v <- sigma_u^2
  if (link == "probit") v / (1 + v) else v / sqrt((1 + v) * (pi^2 / 3 + v))
}

