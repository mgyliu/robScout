gridge <- function(x, cov_method, rho = 0, v = NULL, thetas = NULL, u = NULL) {
  x <- x / sqrt(nrow(x) - 1)
  p <- 2 * rho

  if (is.null(v) || is.null(thetas) || is.null(u)) {
    covarswap <- x %*% t(x) # est_cov(x, method = cov_method)
    eigenswap <- eigen(covarswap)
    keep <- ((eigenswap$values) > 1e-11)
    v <- t(diag(1 / sqrt(eigenswap$values[keep])) %*% t(eigenswap$vectors[, keep]) %*% x)
    thetas <- eigenswap$values[keep]
    u <- eigenswap$vectors[, keep]
  }

  lambda <- sqrt(4 * p) / 2 # = 2 * sqrt(p) / 2 = sqrt(p)
  diagmat <- (-lambda + (-thetas + sqrt(thetas^2 + 4 * p)) / 2)
  dbar <- .5 * (thetas + sqrt(thetas^2 + 4 * p)) - sqrt(p)

  svdstuff <- list(u = u, v = v, thetas = thetas)
  wistuff <- list(
    v = v, firstdiag = (1 / sqrt(p)),
    diagsandwich = (-((1 / p) * (1 / (1 / sqrt(p) + 1 / dbar))))
  )
  wstuff <- list(
    v = v, firstdiag = lambda,
    diagsandwich = (thetas + diagmat)
  )

  cov_est <- diag(rep(wstuff$firstdiag, ncol(x))) + wstuff$v %*% diag(wstuff$diagsandwich) %*% t(wstuff$v)
  return(list(svdstuff = svdstuff, wistuff = wistuff, wstuff = wstuff, cov_est = cov_est))
}

#' "Graphical Ridge" - TODO there are issues with signs in the svd
#' @description
#' `gridge2` computes the closed form solution of \eqn{\Theta^{-1}} which solves
#' \eqn{\Theta^{-1} - 2\lambda_1 \Theta = X^T X}
#' @param x a (n x p) data matrix. it should be standardized
#' @param rho a numeric vector of penalty terms
#' @param svd_x a list with names "d", "u", and "v". By default, gridge2
#' will compute the SVD of x. If svd_x is passed, then gridge will use
#' those values and not compute the SVD of x.
#' @return a list of (svd_x, ws, rho) where:
#' - `svd_x`` is the singular value decomposition of x
#' - `w` is a list of length length(rho) containing the inverse precision matrix
#'    estimates
#' - `rho` is a vector of the penalty parameters that were used
#' @export
gridge2 <- function(x, rho = c(0), svd_x = list()) {
  x <- x / sqrt(nrow(x) - 1)
  if (length(svd_x) == 0) {
    svd_x <- svd(x)
  }

  stopifnot(sort(names(svd_x)) == c("d", "u", "v"))

  # svd_x$d is a vector of length p = ncol(x)
  d_sqr <- svd_x$d^2
  # make a matrix of d_tilde_sqr where each column contains the result for each rho
  ws <- lapply(rho, function(r) {
    d_tilde_sqr <- sqrt(svd_x$d^4 + 8 * r) - sqrt(8 * r)
    w <- svd_x$v %*% diag(d_sqr + d_tilde_sqr) %*% t(svd_x$v)
    w
  })

  return(list(svd_x = svd_x, ws = ws, rho = rho))
}
