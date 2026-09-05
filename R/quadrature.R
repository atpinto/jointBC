# Internal numerical helpers -------------------------------------------------

.gauss_hermite <- function(n) {
  if (length(n) != 1L || !is.finite(n) || n < 1 || n != as.integer(n)) {
    stop("`n_quad` must be a positive integer.", call. = FALSE)
  }
  n <- as.integer(n)
  if (n == 1L) return(list(nodes = 0, weights = sqrt(pi)))

  # Golub-Welsch rule for the weight exp(-x^2).
  off_diag <- sqrt(seq_len(n - 1L) / 2)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(seq_len(n - 1L), 2:n)] <- off_diag
  jacobi[cbind(2:n, seq_len(n - 1L))] <- off_diag
  eig <- eigen(jacobi, symmetric = TRUE)
  ord <- order(eig$values)
  list(
    nodes = eig$values[ord],
    weights = sqrt(pi) * eig$vectors[1L, ord]^2
  )
}

.row_log_sum_exp <- function(x) {
  maxima <- apply(x, 1L, max)
  maxima + log(rowSums(exp(x - maxima)))
}

.binary_loglik <- function(y, eta, link) {
  if (link == "logit") {
    out <- matrix(0, nrow(eta), ncol(eta))
    is_one <- y == 1
    if (any(is_one)) out[is_one, ] <- stats::plogis(eta[is_one, , drop = FALSE], log.p = TRUE)
    if (any(!is_one)) out[!is_one, ] <- stats::plogis(-eta[!is_one, , drop = FALSE], log.p = TRUE)
    return(out)
  }

  out <- matrix(0, nrow(eta), ncol(eta))
  is_one <- y == 1
  if (any(is_one)) out[is_one, ] <- stats::pnorm(eta[is_one, , drop = FALSE], log.p = TRUE)
  if (any(!is_one)) {
    out[!is_one, ] <- stats::pnorm(
      eta[!is_one, , drop = FALSE], lower.tail = FALSE, log.p = TRUE
    )
  }
  out
}

.invert_hessian <- function(hessian) {
  hessian <- (hessian + t(hessian)) / 2
  ev <- eigen(hessian, symmetric = TRUE, only.values = TRUE)$values
  tol <- max(abs(ev)) * sqrt(.Machine$double.eps)
  if (any(!is.finite(ev)) || min(ev) <= tol) {
    warning(
      "The observed Hessian is not positive definite; standard errors are unavailable.",
      call. = FALSE
    )
    return(matrix(NA_real_, nrow(hessian), ncol(hessian)))
  }
  solve(hessian)
}

