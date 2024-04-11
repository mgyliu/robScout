#' @title cov_winsor
#' @description computes an estimate of cov(X,Y) with adjusted
#' multivariate Winsorization, as described in Lafit et al. 2022.
#' This function will standardize the data using median and mad.
#' Hence, if the data has already been standardized with those
#' robust center and scale functions, standardization here
#' will have no effect.
#' @param X a nxp matrix
#' @param Y a nx1 vector or NULL
#' @return a pxp matrix if Y is NULL; otherwise, a length-p vector.
#' @export
cov_winsor <- function(X, Y = NULL) {
    p <- ncol(X)
    dispersion_x <- unname(apply(X, MARGIN = 2, FUN = stats::mad))
    dispersion_y <- stats::mad(Y)
    X <- robustHD::robStandardize(X, median, mad)
    if (!is.null(Y)) Y <- robustHD::robStandardize(Y, median, mad)

    if (is.null(Y)) { # If Y is NULL, compute cov(X)
        cormat <- matrix(NA, nrow = p, ncol = p)
        # j and k index the columns of X
        for (j in 1:p) {
            for (k in 1:p) {
                cormat[j, k] <- robustHD::corHuber(X[, j], X[, k], type = "bivariate", standardized = TRUE)
            }
        }
        winsor_cov <- diag(dispersion_x) %*% cormat %*% diag(dispersion_x)
        return(as.matrix(Matrix::nearPD(winsor_cov)$mat))
    } else { # If Y is not NULL, compute cov(X,Y)
        cormat <- matrix(NA, nrow = p, ncol = 1)
        for (j in 1:p) {
            cormat[j, ] <- robustHD::corHuber(X[, j], Y, type = "bivariate", standardized = TRUE)
        }
        return(cormat * dispersion_x * dispersion_y)
    }
}

cov_wrap <- function(X, Y = NULL) {
    locScaleX <- cellWise::estLocScale(X)
    Xw <- cellWise::wrap(X, locScaleX$loc, locScaleX$scale)$Xw
    if (is.null(Y)) {
        return(cov(Xw))
    } else {
        locScaleY <- cellWise::estLocScale(Y)
        Yw <- cellWise::wrap(Y, locScaleY$loc, locScaleY$scale)$Xw
        return(cov(Xw, Yw))
    }
}

#' @title est_cov
#' @description Computes a covariance estimate based on the type
#' specificed by the user. This function assumes that the data
#' has already been standardized.
#' @param X (n x p) predictor matrix
#' @param Y (n x 1) response vector or NULL
#' @param method string indicating which covariance estimator
#' to compute. One of "default", "winsor", or "wrap".
#' @return a (p x p) covariance matrix or length-p vector
#' @export
est_cov <- function(X, Y = NULL, method = c("default", "winsor", "wrap")) {
    if (method == "default") {
        cov(X, Y)
    } else if (method == "winsor") {
        cov_winsor(X, Y)
    } else if (method == "wrap") {
        cov_wrap(X, Y)
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
