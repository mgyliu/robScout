#' @description computes an estimate of cov(X,Y) with adjusted
#' multivariate Winsorization, as described in Lafit et al. 2022.
#' @param X a nxp matrix
#' @param Y a nx1 vector or NULL
#' @return a pxp matrix if Y is NULL; otherwise, a length-p vector.
cov_winsor <- function(X, Y = NULL) {
    p <- ncol(X)
    dispersion_x <- unname(apply(X, MARGIN = 2, FUN = stats::mad))
    dispersion_y <- stats::mad(Y)
    # If Y is NULL, compute cov(X)
    if (is.null(Y)) {
        cormat <- matrix(NA, nrow = p, ncol = p)
        # j and k index the columns of X
        for (j in 1:p) {
            for (k in 1:p) {
                cormat[j, k] <- robustHD::corHuber(X[, j], X[, k], type = "bivariate", standardized = FALSE)
            }
        }
        winsor_cov <- diag(dispersion_x) %*% cormat %*% diag(dispersion_x)
        return(as.matrix(Matrix::nearPD(winsor_cov)$mat))
    }

    # If Y is not NULL, compute cov(X,Y)
    else {
        cormat <- rep(NA, p)
        for (j in 1:p) {
            cormat[j] <- robustHD::corHuber(X[, j], Y, type = "bivariate", standardized = FALSE)
        }
        return(cormat * dispersion_x * dispersion_y)
    }
}

#' @description Computes a covariance estimate based on the type
#' specificed by the user. This function assumes that the data
#' has already been standardized.
#' @param X (n x p) predictor matrix
#' @param Y (n x 1) response vector or NULL
#' @return a (p x p) covariance matrix or length-p vector
est_cov <- function(X, Y = NULL, method = c("default", "winsor")) {
    if (method == "default") {
        cov(X, Y)
    } else if (method == "winsor") {
        cov_winsor(X, Y)
    } else {
        warning(glue::glue("Cov method {method} is not implemented. Using default."))
        cov(X, Y)
    }

    # TODO implement other robust covariance estimates
    # if (type == "qn") {}
    # if (type == "tau") {}
    # if (type == "gauss") {}
    # if (type == "spearman") {}
    # if (type == "quadrant") {}
}
