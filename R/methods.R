#' @export
coef.joint_bc <- function(object, ...) object$coefficients

#' @export
vcov.joint_bc <- function(object, ...) object$vcov

#' @export
logLik.joint_bc <- function(object, ...) {
  out <- object$logLik
  attr(out, "df") <- length(object$coefficients)
  attr(out, "nobs") <- object$nobs
  class(out) <- "logLik"
  out
}

#' @export
summary.joint_bc <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  z <- object$coefficients / se
  table <- cbind(
    Estimate = object$coefficients,
    `Std. Error` = se,
    `z value` = z,
    `Pr(>|z|)` = 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  )
  structure(list(
    call = object$call,
    coefficients = table,
    logLik = object$logLik,
    convergence = object$convergence,
    message = object$message,
    link = object$link,
    nobs = object$nobs,
    n_quad = object$n_quad,
    latent_correlation = latent_correlation(object)
  ), class = "summary.joint_bc")
}

#' @export
print.summary.joint_bc <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Joint binary-continuous model\n")
  cat("Binary link:", x$link, " | subjects:", x$nobs,
      " | quadrature nodes:", x$n_quad, "\n")
  cat("Log likelihood:", formatC(x$logLik, digits = digits, format = "f"),
      " | convergence code:", x$convergence, "\n")
  cat("Implied latent correlation:", formatC(x$latent_correlation, digits = digits), "\n\n")
  stats::printCoefmat(x$coefficients, digits = digits)
  invisible(x)
}

#' @export
print.joint_bc <- function(x, ...) {
  cat("Joint binary-continuous model (", x$link, ")\n", sep = "")
  cat("Subjects:", x$nobs, " | log likelihood:", format(x$logLik),
      " | convergence code:", x$convergence, "\n")
  print(x$coefficients)
  invisible(x)
}
