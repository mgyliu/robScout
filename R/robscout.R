#' @title robscout
#' @description
#' @param X
#' @param Y
#' @param lam1
#' @param lam2
#' @param standardize
robscout <- function(X, Y, est_type = c("glasso", "gridge"), lam1, lam2, centerFun = stats::median, scaleFun = stats::mad) {

}


# This code is mostly copied from https://github.com/cran/scout/blob/master/R/scout.R#L237-326
# I've added formatting improvements and removed code related to plotting
# Calls to `cov` are replaced with `est_cov`, the robustified version.
cv.robscout <- function(X, Y, K = 10, lam1s, lam2s, p1 = 1, p2 = 1, rescale = TRUE, ...) {
    call <- match.call()
    if (K == 1) stop("You can't do 1-fold cross-validation! Please use K > 1.")
    if (K > length(y) / 2) stop("Please choose a value of K between 2 and length(y)/2.")
    if (p1 == 0 && p2 == 0) stop("Why would you want to cross-validate least squares?")
    if (is.null(p1)) lam1s <- 0
    if (is.null(p2)) lam2s <- 0
    lam1s <- c(lam1s)
    lam2s <- c(lam2s)
    if (length(lam1s) < 2 && length(lam2s) < 2) stop("Not a reasonable range of lambdas over which to be cross-validating")

    all.folds <- cv.folds(length(y), K)

    if (length(lam1s) > 1 && length(lam2s) > 1) {
        residmat <- array(0, dim = c(length(lam1s), length(lam2s), K))
        for (i in seq(K)) {
            if (trace) cat("\n CV Fold", i, "\t")
            omit <- all.folds[[i]]
            fit <- scout(x[-omit, ], y[-omit], newx = x[omit, ], p1 = p1, p2 = p2, lam1s = lam1s, lam2s = lam2s, rescale = rescale, trace = trace)
            residmat[, , i] <- apply((sweep(fit$yhat, 3, y[omit], "-"))^2, c(1, 2), mean)
        }
        cv <- apply(residmat, c(1, 2), mean)
        cv.error <- sqrt(apply(residmat, c(1, 2), var) / K)
        object <- list(p1 = p1, p2 = p2, lam1s = lam1s, lam2s = lam2s, cv = cv, cv.error = cv.error, call = call, bestlam1 = lam1s[which.min(apply(cv, 1, min))], bestlam2 = lam2s[which.min(apply(cv, 2, min))])
    } else if (length(lam1s) == 1) {
        lam1 <- lam1s[1]
        residmat <- matrix(0, nrow = length(lam2s), ncol = K)
        for (i in seq(K)) {
            if (trace) cat("\n CV Fold", i, "\t")
            omit <- all.folds[[i]]
            fit <- scout(x[-omit, ], y[-omit], newx = x[omit, ], p1 = p1, p2 = p2, lam1s = lam1, lam2s = lam2s)
            residmat[, i] <- apply(sweep(fit$yhat[1, , ], 2, y[omit], "-")^2, 1, mean)
        }
        cv <- apply(residmat, 1, mean)
        cv.error <- sqrt(apply(residmat, 1, var) / K)
        object <- list(p1 = p1, p2 = p2, lam1s = lam1, lam2s = lam2s, cv = cv, cv.error = cv.error, call = call, bestlam1 = lam1, bestlam2 = lam2s[which.min(cv)])
    } else if (length(lam2s) == 1) {
        lam2 <- lam2s[1]
        residmat <- matrix(0, nrow = length(lam1s), ncol = K)
        for (i in seq(K)) {
            if (trace) cat("\n CV Fold", i, "\t")
            omit <- all.folds[[i]]
            fit <- scout(x[-omit, ], y[-omit], newx = x[omit, ], p1 = p1, p2 = p2, lam1s = lam1s, lam2s = lam2)
            residmat[, i] <- apply(sweep(fit$yhat[, 1, ], 2, y[omit], "-")^2, 1, mean)
        }
        cv <- apply(residmat, 1, mean)
        cv.error <- sqrt(apply(residmat, 1, var) / K)
        object <- list(p1 = p1, p2 = p2, lam1s = lam1s, lam2s = lam2, cv = cv, cv.error = cv.error, call = call, bestlam1 = lam1s[which.min(cv)], bestlam2 = lam2)
    }
    if (trace) cat("\n")
    object$call <- call
    object$folds <- all.folds
    class(object) <- "cvobject"
    invisible(object)
}
