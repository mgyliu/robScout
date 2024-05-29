#' Scout(1,.) Stepwise
#'
#' @description
#' Computes the best inverse covariance estimate using the glasso and uses that
#' in the second step of the Scout.
#'
#' @param X n x p data matrix
#' @param Y n x 1 response vector
#' @param p2 regularization type in second step. one of NULL or 1
#' @param K number of folds to use in CV, only applicable for p2 = 1
#' @param nlambda1 number of regularization params to use in glasso step
#' @param lambda1.min.ratio smallest value of lambda1 as a fraction of lambda1_max
#' @param nlambda2 number of regularization params to use in lasso step
#' @param lambda2.min.ratio smallest value of lambda2 as a fraction of lambda2_max
#' @param cov_method a string indicating which covariance matrix to use. See below for details
#' @param glasso_crit a string indicating which criteria to minimize in glasso. One of
#' "bic" or "loglik"
#' @param ddc_first whether or not to detect and impute cellwise outliers using
#' DDC before computing estimates
#' @param ddc_with_response whether or not to concatenate Y before running DDC. Only
#' applicable if `ddc_first` is `TRUE`
#' @param standardize whether or not to standardize variables before
#' computing estiamtes
#' @param centerFun a function to compute an estiamte of the center of a variable.
#' This is ignored if standardize is `FALSE`
#' @param scaleFun a function to compute an estimate of the scale of a variable.
#' This is ignored if standardize is `FALSE`
#' @param rescale should coefficients be re-scaled by a constant?
#' @details `cov_method` can be one of the following:
#' * `default`: computes the default covariance with \code{\link{cov}}
#' * `winsor`: computes the covariance matrix based on adjusted multivariate Winsorization,
#' as seen in Lafit et al. 2022.
#' @return list of
#' * `g.res`: result of graphical lasso step
#' * `cv.res`: result of cross-validated lasso step, or NA if p2 = NULL
#' * `mod`: final scout model using cross-validated parameters
#' @export
scout_1something_stepwise <- function(X, Y, p2, K = 5,
                                      nlambda1 = 100, lambda1.min.ratio = 0.1,
                                      nlambda2 = 100, lambda2.min.ratio = 0.1,
                                      cov_method = "default",
                                      glasso_crit = "bic",
                                      ddc_first = FALSE,
                                      ddc_with_response = FALSE,
                                      standardize = TRUE,
                                      centerFun = mean, scaleFun = sd,
                                      rescale = TRUE) {
    if (ddc_first) {
        if (ddc_with_response) {
            ddc_res <- run_ddc(X, Y)
            Y <- ddc_res$Y
        } else {
            ddc_res <- run_ddc(X)
        }
        X <- ddc_res$X
    }

    # Use standardized X to find best lambda1
    X_std <- robustHD::robStandardize(X, centerFun = centerFun, scaleFun = scaleFun)
    g.res <- glasso_select(
        X_std, X_std, # X and X test are the same, assuming we're using a BIC or loglik
        cov_method = cov_method, crit = glasso_crit,
        nlambda = nlambda1, lambda.min.ratio = lambda1.min.ratio,
    )

    if (is.null(p2)) {
        cv.res <- NA
        mod <- scout(
            X, Y,
            p1 = 1, p2 = p2,
            lam1s = g.res$best_lambda, standardize = standardize,
            centerFun = centerFun, scaleFun = scaleFun,
            rescale = rescale,
            cov_method = cov_method
        )
    } else if (p2 == 1) {
        cv.res <- cv.scout(
            X, Y,
            K = K,
            p1 = 1, p2 = p2,
            lam1s = c(g.res$best_lambda),
            nlambda2 = nlambda2,
            lambda2.min.ratio = lambda2.min.ratio,
            standardize = standardize,
            rescale = rescale,
            centerFun = centerFun, scaleFun = scaleFun,
            cov_method = cov_method
        )
        mod <- scout(
            X, Y,
            p1 = 1, p2 = p2,
            lam1s = g.res$best_lambda,
            lam2s = cv.res$bestlam2,
            standardize = standardize,
            rescale = rescale,
            centerFun = centerFun, scaleFun = scaleFun,
            cov_method = cov_method
        )
    } else {
        stop(glue::glue("scout_lasso_stepwise not implemented for p2 = {deparse(p2)}"))
    }

    list(g.res = g.res, cv.res = cv.res, mod = mod)
}
