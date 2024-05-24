# This code was copied from the implemented of scout::scout and modified.

#' Covariance-regularized regression, aka the Scout.
#' @description
#' The main function of the "scout" package. Performs
#' covariance-regularized regression. Required inputs are an x matrix of
#' features (the columns are the features) and a y vector of
#' observations. By default, Scout(2,1) is performed; however, $p_1$ and
#' $p_2$ can be specified (in which case Scout($p_1$, $p_2$) is
#' performed). Also, by default Scout is performed over a grid of lambda1
#' and lambda2 values, but a different grid of values (or individual
#' values, rather than an entire grid) can be specified.
#' @param x A matrix of predictors, where the rows are the samples and
#' the columns are the predictors
#' @param y A matrix of observations, where length(y) should equal nrow(x)
#' @param newx An *optional* argument, consisting of a matrix with
#' ncol(x) columns, at which one wishes to make predictions for each
#' (lam1,lam2) pair.
#' @param p1 The $L_p$ penalty for the covariance regularization. Must be
#' one of 1, 2, or NULL. NULL corresponds to no covariance
#' regularization.  WARNING: When p1=1, and ncol(x)>500, Scout can be
#' SLOW. We recommend that for very large data sets, you use Scout with
#' p1=2. Also, when ncol(x)>nrow(x) and p1=1, then very small values of
#' lambda1 (lambda1 < 1e-4) will cause problems with graphical lasso,
#' and so those values will be automatically increased to 1e-4.
#' @param p2 The $L_p$ penalty for the estimation of the regression
#' coefficients based on the regularized covariance matrix. Must be one
#' of 1 (for $L_1$ regularization) or NULL (for no regularization).
#' @param lam1s The (vector of) tuning parameters for regularization of the
#' covariance matrix. Can be NULL if p1=NULL, since then no covariance
#' regularization is taking place. If p1=1 and nrow(x)<ncol(x), then the no value in lam1s
#' should be smaller than 1e-3, because this will cause graphical lasso
#' to take too long. Also, if ncol(x)>500 then we really do not
#' recommend using p1=1, as graphical lasso can be uncomfortably slow.
#' @param lam2s The (vector of) tuning parameters for the $L_1$ regularization of
#' the regression coefficients, using the regularized covariance
#' matrix. Can be NULL if p2=NULL. (If p2=NULL, then non-zero lam2s
#' have no effect). A value of 0 will result in no
#' regularization.
#' @param rescale Should coefficients beta obtained by
#' covariance-regularized regression be re-scaled by a constant, given
#' by regressing $y$ onto $x beta$? This is done in Witten and
#' Tibshirani (2008) and is important for good performance. Default is
#' TRUE.
#' @param standardize Should the columns of x be scaled to have standard deviation
#' 1, and should y be scaled to have standard deviation 1, before
#' covariance-regularized regression is performed? This affects the
#' meaning of the penalties that are applied. In general,
#' standardization should be performed. Default is TRUE.
#' @return
#' * intercepts - Returns a matrix of intercepts, of dimension
#'                length(lam1s) x length(lam2s)
#' * coefficients - Returns an array of coefficients, of dimension
#'                  length(lam1s) x length(lam2s) x ncol(x)
#' * p1 - p1 value used
#' * p2 - p2 value used
#' * lam1s - lam1s used
#' * lam2s - lam2s used
#' @export
scout <- function(x, y, newx = NULL, p1 = 2, p2 = 1,
                  lam1s = seq(.001, .2, len = 10),
                  lam2s = seq(.001, .2, len = 10),
                  rescale = TRUE, standardize = TRUE,
                  centerFun = mean, scaleFun = sd,
                  cov_method = c("default", "winsor")) {
    call <- match.call()
    if (!is.null(p1) && p1 != 1 && p1 != 2) stop("p1 must be 1, 2, or NULL.")
    if (!is.null(p2) && p2 != 1) stop("p1 must be 1 or NULL.")
    if ((sum(is.na(x)) + sum(is.na(y))) > 0) stop("Please fix the NAs in your data set first. Missing values can be imputed using library 'impute'.")
    x <- as.matrix(x)
    if (min(apply(x, 2, sd)) == 0) stop("Please do not enter an x matrix with variables that are constant.")

    # Need to center and scale x,y
    meany <- centerFun(y)
    meanx <- apply(x, 2, centerFun)
    if (standardize) {
        sdy <- scaleFun(y)
        sdx <- apply(x, 2, scaleFun)
        x <- robustHD::standardize(x, centerFun, scaleFun)
    } else {
        # x <- scale(x, T, F)
        sdx <- rep(1, ncol(x))
        sdy <- 1
    }
    y <- (y - meany) / sdy

    # Want to re-order lam2s in increasing order
    if (length(lam2s) > 0) {
        lam2s.orig <- lam2s
        s2.out <- sort(lam2s.orig, index.return = TRUE)
        lam2s <- lam2s.orig[s2.out$ix]
    }
    # Done re-ordering lam2s
    # Want to re-order lam1s in increasing order
    if (length(lam1s) > 0) {
        lam1s.orig <- lam1s
        s1.out <- sort(lam1s.orig, index.return = TRUE)
        lam1s <- lam1s.orig[s1.out$ix]
    }
    # Done re-ordering lam1s
    if (is.null(lam1s) || is.null(p1) || sum(lam1s) == 0) {
        p1 <- 0
        lam1s <- c(0)
        lam1s.orig <- c(0)
    }
    if (is.null(lam2s) || is.null(p2) || sum(lam2s) == 0) {
        p2 <- 0
        lam2s <- c(0)
        lam2s.orig <- c(0)
    }
    if (p1 != 0 && p1 != 1 && p1 != 2) stop("p1 must be 1, 2, or 0")
    if (p2 != 0 && p2 != 1) stop("p2 must be 1 or 0")
    if (!is.null(lam1s) && min(lam1s) < 0) stop("lam1s cannot be negative")
    if (!is.null(lam2s) && min(lam2s) < 0) stop("lam2s cannot be negative")
    if (ncol(x) >= nrow(x) && p1 == 0 && p2 == 0) {
        stop("p1 and p2 cannot both be zero when p>n.")
    }
    c.betamat <- NA
    if (p1 == 0) {
        betamat <- array(NA, dim = c(1, length(lam2s), ncol(x)))
        for (j in 1:length(lam2s)) {
            if (lam2s[j] == 0) {
                if (ncol(x) >= nrow(x)) stop("Cannot do Least Squares when ncol(x)>=nrow(x)")
                beta <- lsfit(x, y, intercept = FALSE)$coef
                betamat[1, j, ] <- beta
            } else {
                if (j == 1) beta <- scout::crossProdLasso(est_cov(x, method = cov_method), est_cov(x, y, method = cov_method), rho = lam2s[j])$beta
                if (j != 1) {
                    if (sum(abs(beta)) != 0 || lam2s[j] < lam2s[j - 1]) {
                        beta <- scout::crossProdLasso(est_cov(x, method = cov_method), est_cov(x, y, method = cov_method), rho = lam2s[j], beta.init = beta)$beta
                        # if got zero for smaller value of lambda 2,
                        # then no need to keep computing!!!
                    }
                }
                if (rescale && sum(abs(beta)) != 0) beta <- beta * lsfit(x %*% beta, y, intercept = FALSE)$coef
                betamat[1, j, ] <- beta
            }
        }
    } else if (p1 == 1) {
        # betamat <- scout1something(x, y, p2, lam1s, lam2s, rescale, cov_method)
        betamat <- scout1something.huge(x, y, p2, lam1s, lam2s, rescale, cov_method)
    } else if (p1 == 2) {
        betamat <- scout2something(x, y, p2, lam1s, lam2s, rescale, cov_method)
    }

    interceptmat <- matrix(meany, nrow = length(lam1s), ncol = length(lam2s))
    for (i in 1:length(lam1s)) {
        for (j in 1:length(lam2s)) {
            interceptmat[i, j] <- interceptmat[i, j] - sum((sdy * meanx / sdx) * betamat[i, j, ])
        }
    }
    betamat <- sweep(betamat, 3, sdy / sdx, "*")
    betamat <- array(betamat[rank(lam1s.orig), rank(lam2s.orig), ], dim = c(length(lam1s.orig), length(lam2s.orig), ncol(x)))
    interceptmat <- matrix(interceptmat[rank(lam1s.orig), rank(lam2s.orig)], nrow = length(lam1s.orig), ncol = length(lam2s.orig))
    scout.obj <- (list(intercepts = interceptmat, coefficients = betamat, p1 = p1, p2 = p2, lam1s = lam1s.orig, lam2s = lam2s.orig, yhat = NULL, call = call))
    class(scout.obj) <- "scoutobject"
    if (!is.null(newx)) {
        yhat <- predict(scout.obj, newx)
        scout.obj$yhat <- yhat
    }
    return(scout.obj)
}
