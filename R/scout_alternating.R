#' Computes covariance matrix estimate using either L1 or L2 penalty
#' @param X n x p data matrix. It is assumed to be standardized
#' @param p1 numeric value, either 1 or 2, corresponding to L1 and L2
#' regularization respectively. If p1 = 1, huge::huge.glasso is used to estimate
#' a sparse covariance matrix. If p1 = 2, scout::gridge2 is used to estimate
#' a covariance matrix.
#' @param lam1s numeric vector of penalty parameters for the regularization
#' @return a list of length length(lam1s) where each entry is a covariance
#' matrix estimate. List indices correspond to lam1s indices.
get_xtxs <- function(X, p1, lam1s, cov_method) {
  if (p1 == 1) {
    g.out <- huge::huge.glasso(est_cov(X, method = cov_method), lambda = lam1s, verbose = FALSE, cov.output = TRUE)
    return(g.out$cov)
  } else if (p1 == 2) {
    g.out.svd <- gridge(X, cov_method, rho = lam1s[1])$svdstuff
    ws <- lapply(lam1s, function(lam1) {
      gridge(X, cov_method, rho = lam1, v = g.out.svd$v, thetas = g.out.svd$thetas, u = g.out.svd$u)$cov_est
    })
    return(ws)
  } else {
    stop(glue::glue("get_xtxs is not implemented for p1 = {p1}"))
  }
}

#' betas_to_original_scale
#' @description rescales betas to original scale of training data using the
#' standard deviation estimates of the X and Y training data.
#' @param betamat a matrix where each column is a beta_hat estimate for some
#' regularization parameter
#' @param sdx a numeric vector. column SDs of original training X data
#' @param sdy a numeric value. SD of original training Y data
#' @return betamat where each column is scaled by sdy/sdx
betas_to_original_scale <- function(betamat, sdx, sdy) {
  sweep(betamat, 1, sdy / sdx, "*")
}

#' compute_intercepts
#' @description computes the intercept term for each beta_hat estimate.
#' The beta matrix should be such that each column contains a different beta_hat
#' @param betamat a matrix where each column is a beta_hat estimate for some
#' regularization parameter
#' @param meanx a numeric vector. column means of original training X data
#' @param sdx a numeric vector. column SDs of original training X data
#' @param meany a numeric value. mean of original training Y data
#' @param sdy a numeric value. SD of original training Y data
compute_intercepts <- function(betamat, meanx, sdx, meany, sdy) {
  apply(betamat, MARGIN = 2, FUN = function(beta_hat) {
    meany + sum((sdy * meanx / sdx) * beta_hat)
  })
}

#' compute_rmspe_betamat
#' @description computes the rmspe between the provided Y_test and the predicted
#' X_test %*% beta_hat for each beta_hat in the provided beta matrix.
#' The beta matrix should be such that each column contains a different beta_hat
#' @param X_test n_test x p data matrix (test data, unstandardized)
#' @param Y_test n_test x 1 response vector (test data, unstandardized)
#' @param intercepts numeric vector of length nlambda
#' @param betamat numeric matrix of dimensions p x nlambda
compute_rmspe_betamat <- function(X_test, Y_test, intercepts, betamat) {
  errors <- apply(rbind(intercepts, betamat), MARGIN = 2, FUN = function(beta_hat) {
    yhat <- beta_hat[1] + X_test %*% beta_hat[-1]
    perry::rmspe(yhat, Y_test)
  })
}

#' Rescales beta estimates as in Step 4 of the Scout algo (Witten & Tibshirani)
#' @param X standardized training data
#' @param Y standardized training response
#' @param betamat a matrix where each column is a beta_hat estimate for some
#' regularization parameter
rescale_betas <- function(X, Y, betamat) {
  rescaled <- apply(betamat, MARGIN = 2, FUN = function(beta_hat) {
    if (all(beta_hat == 0)) {
      rep(0, length(beta_hat))
    } else {
      beta_hat * lsfit(X %*% beta_hat, Y, intercept = FALSE)$coef
    }
  })
  rescaled
}

# assume: p2 = 1
# get_lam2_sequence(cov_x_est)
get_best_lam2_lasso <- function(X, Y, X_test, Y_test,
                                meanx, meany, sdx, sdy,
                                cov_x_est, lam2s, rescale, cov_method) {
  # Initialize beta_hat estimates. The beta_hat for a given lambda2 goes into
  # that lambda2's column index.
  betamat <- matrix(NA, nrow = ncol(X), ncol = length(lam2s))
  beta <- NA
  for (i in 1:length(lam2s)) {
    if (i == 1) {
      beta <- scout::crossProdLasso(cov_x_est, est_cov(X, Y, method = cov_method), rho = lam2s[i])$beta
    } else {
      beta <- scout::crossProdLasso(cov_x_est, est_cov(X, Y, method = cov_method), rho = lam2s[i], beta.init = beta)$beta
    }
    betamat[, i] <- beta
  }

  # Rescale betas if needed (Scout procedure Step 4)
  if (rescale) {
    betamat <- rescale_betas(X, Y, betamat)
  }

  # Estimate intercepts
  intercepts <- compute_intercepts(betamat, meanx, sdx, meany, sdy)

  # Rescale betas to the original scale of the training data
  betamat <- betas_to_original_scale(betamat, sdx, sdy)

  # Compute yhats based on estimated intercepts, beta matrix, and test data
  errors <- compute_rmspe_betamat(X_test, Y_test, intercepts, betamat)

  # Get best lam2 and corresponding coefficients
  list(
    best_lam2_idx = which.min(errors),
    best_lam2 = lam2s[which.min(errors)],
    intercept = intercepts[which.min(errors)],
    beta_hat = betamat[, which.min(errors)],
    error = min(errors)
  )
}

get_best_lam1 <- function(X, Y, X_test, Y_test,
                          meanx, meany, sdx, sdy,
                          p1, lam1s, xtxs, lam2, rescale, cov_method) {
  # xtxs <- get_xtxs(X, p1, lam1s)
  # Initialize beta_hat estimates. The beta_hat for a given lambda2 goes into
  # that lambda2's column index.
  betamat <- matrix(NA, nrow = ncol(X), ncol = length(lam1s))
  beta <- NA
  for (i in 1:length(lam1s)) {
    if (i == 1) {
      beta <- scout::crossProdLasso(xtxs[[i]], est_cov(X, Y, cov_method), rho = lam2)$beta
    } else {
      beta <- scout::crossProdLasso(xtxs[[i]], est_cov(X, Y, cov_method), rho = lam2, beta.init = beta)$beta
    }
    betamat[, i] <- beta
  }

  # Rescale betas if needed (Scout procedure Step 4)
  if (rescale) {
    betamat <- rescale_betas(X, Y, betamat)
  }

  # Estimate intercepts
  intercepts <- compute_intercepts(betamat, meanx, sdx, meany, sdy)

  # Rescale betas to the original scale of the training data
  betamat <- betas_to_original_scale(betamat, sdx, sdy)

  # Compute yhats based on estimated intercepts, beta matrix, and test data
  errors <- compute_rmspe_betamat(X_test, Y_test, intercepts, betamat)

  # Get best lam1 and corresponding coefficients
  list(
    best_lam1_idx = which.min(errors),
    best_lam1 = lam1s[which.min(errors)],
    best_xtx = xtxs[[which.min(errors)]],
    intercept = intercepts[which.min(errors)],
    beta_hat = betamat[, which.min(errors)],
    error = min(errors)
  )
}

#' @title scout_alternating_lasso
#' @description implements an alternating search algorithm to efficiently
#' search over a grid of lambda1 x lambda2. L_p penalty for coefficient
#' regularization is assumed to be 1 (i.e., p2 = 1)
#' @param X A matrix of predictors. Rows are observations and columns
#' are variables
#' @param Y A vector of outcomes
#' @param X_test A matrix of predictors (test set)
#' @param Y_test A vector of outcomes (test set)
#' @param p1 L_p penalty for covariance regularization. must be 1 or 2
#' @param nlambda1 number of regularization terms to use in covariance
#' regularization step
#' @param nlambda2 number of regularization terms to use in coefficient
#' regularization step
#' @param lambda1.min.ratio smallest fraction of lambda1_max to include in
#' lambda1 sequence
#' @param lambda2.min.ratio smallest fraction of lambda2_max to include in
#' lambda2 sequence
#' @param tol convergence tolerance for difference between previous estimated
#' rmspe and current estimated rmspe. if the difference becomes smaller than
#' this value, the algorithm stops and returns the current lam1 and lam2
#' combination
#' @param max_iter maximum number of alternating iterations to compute before
#' stopping. the loop will stop when either the number of iterations exceeds
#' max_iter, or when the errors converge to have a diff below the tolerance
#' @param rescale Should coefficients beta obtained by
#' covariance-regularized regression be re-scaled by a constant, given
#' by regressing $y$ onto $x beta$? This is done in Witten and
#' Tibshirani (2008) and is important for good performance. Default is
#' TRUE.
#' @param standardize whether or not to scale the training X and Y data
#' before performing estimation
#' @param lam1_init one of "random" or "max"
#' @return list with items:
#' * errors: vector of numeric RMSPE values
#' * betas: list of vectors of estimated beta_hats
#' * intercepts: vector of numeric intercept values
#' * lambda_pairs: list of vectors of length 2, each vector is (lam1, lam2)
#' * lambda2_paths: list of vectors of lambda2 paths. user usually doesn't need
#' this but it's provided as a convenience
#' Incides of each list or vector (other than lambda2_paths) correspond to
#' the lambda pair at that index in lambda_pairs
#' @export
scout_alternating_lasso <- function(X, Y, X_test, Y_test, p1,
                                    nlambda1 = 100, nlambda2 = 100,
                                    lambda1.min.ratio = 0.1, lambda2.min.ratio = 0.001,
                                    tol = 1e-4, max_iter = 10,
                                    rescale = TRUE, standardize = TRUE,
                                    centerFun = mean, scaleFun = sd,
                                    lam1_init = "random",
                                    cov_method = "default") {
  # Can only pass in one lam1_init value
  stopifnot(lam1_init %in% c("random", "max", "min"))

  # Standardize training data if needed
  meany <- centerFun(Y)
  meanx <- apply(X, 2, centerFun)
  if (standardize) {
    sdy <- scaleFun(Y)
    sdx <- apply(X, 2, scaleFun)
    X <- robustHD::standardize(X, centerFun = centerFun, scaleFun = scaleFun)
  } else {
    # center with centerFun
    X <- apply(X, 2, function(xi) xi - centerFun(xi))
    sdx <- rep(1, ncol(X))
    sdy <- 1
  }
  Y <- (Y - meany) / sdy

  # diff keeps track of the difference between current RMSE and previous RMSE
  # niter keeps track of how many loops we did so far
  diff <- tol + 1
  niter <- 0

  # Get lambda 1 sequence
  lam1s <- get_lambda1_path(est_cov(X, method = cov_method), p1, nlambda1, lambda1.min.ratio)
  xtxs <- get_xtxs(X, p1, lam1s, cov_method)

  # Compute inital lambda1
  if (lam1_init == "random") {
    lam1 <- lam1s[sample(1:nlambda1, 1)]
  } else if (lam1_init == "max") {
    lam1 <- lam1s[1]
  } else {
    lam1 <- tail(lam1s, 1)
  }

  # Compute initial covariance estimate
  cov_x_est <- xtxs[[1]]

  # More to keep track of: errors, betas, intercepts, lambda pairs
  errors <- c()
  betas <- list()
  intercepts <- c()
  lambda_pairs <- list()
  lam2_paths <- list()

  # Iterate until either:
  # (a) difference in curr and prev errors is <= tol
  # (b) we ran for more than max_iter iterations
  # In each iteration:
  # - get best lam2 for this lam1
  # - get best lam1 for the lam2 in prev step
  # - update covariance estimate for next iteration
  while (diff > tol | niter <= max_iter) {
    niter <- niter + 1

    # Compute best lambda2 for the current cov_x_est
    cov_xy_est <- est_cov(X, Y, method = cov_method)
    lam2s <- get_lambda2_path(X, Y, cov_x_est, cov_xy_est, 1, nlambda2, lambda2.min.ratio)
    best_lam2_res <- get_best_lam2_lasso(
      X, Y, X_test, Y_test,
      meanx, meany, sdx, sdy,
      cov_x_est, lam2s, rescale, cov_method
    )
    lambda_pairs <- c(lambda_pairs, list(c(lam1, best_lam2_res$best_lam2)))
    errors <- c(errors, best_lam2_res$error)
    intercepts <- c(intercepts, best_lam2_res$intercept)
    betas <- c(betas, list(best_lam2_res$beta_hat))
    lam2_paths <- c(lam2_paths, list(lam2s))

    # Update rmspe diff
    if (niter > 1) {
      diff <- abs(diff(tail(errors, 2)))
    }

    if (diff <= tol | niter > max_iter) {
      break
    }

    niter <- niter + 1

    # Update covariance est for next iteration
    best_lam1_res <- get_best_lam1(
      X, Y, X_test, Y_test,
      meanx, meany, sdx, sdy,
      p1, lam1s, xtxs, best_lam2_res$best_lam2, rescale, cov_method
    )
    lambda_pairs <- c(lambda_pairs, list(c(best_lam1_res$best_lam1, best_lam2_res$best_lam2)))
    errors <- c(errors, best_lam1_res$error)
    intercepts <- c(intercepts, best_lam1_res$intercept)
    betas <- c(betas, list(best_lam1_res$beta_hat))
    lam2_paths <- c(lam2_paths, list(lam2s))

    # Update rmspe diff
    if (niter > 1) {
      diff <- abs(diff(tail(errors, 2)))
    }

    cov_x_est <- best_lam1_res$best_xtx
    lam1 <- best_lam1_res$best_lam1
  }

  return(list(
    errors = errors,
    betas = betas,
    intercepts = intercepts,
    lambda_pairs = lambda_pairs,
    lam2_paths = lam2_paths,
    lam1s = lam1s
  ))
}

#' cv.scout_alternating_lasso
#' @description implements an alternating search algorithm to efficiently
#' search over a grid of lambda1 x lambda2. L_p penalty for coefficient
#' regularization is assumed to be 1 (i.e., p2 = 1). splits training data into
#' K folds and picks the result corresponding to fold with lowest RMSPE on
#' validation data.
#' @param X A matrix of predictors. Rows are observations and columns
#' are variables
#' @param Y A vector of outcomes
#' @param p1 L_p penalty for covariance regularization. must be 1 or 2
#' @param nlambda1 number of regularization terms to use in covariance
#' regularization step
#' @param nlambda2 number of regularization terms to use in coefficient
#' regularization step
#' @param lambda1.min.ratio smallest fraction of lambda1_max to include in
#' lambda1 sequence
#' @param lambda2.min.ratio smallest fraction of lambda2_max to include in
#' lambda2 sequence
#' @param K number of folds to use in cross-validation
#' @param tol convergence tolerance for difference between previous estimated
#' rmspe and current estimated rmspe. if the difference becomes smaller than
#' this value, the algorithm stops and returns the current lam1 and lam2
#' combination
#' @param max_iter maximum number of alternating iterations to compute before
#' stopping. the loop will stop when either the number of iterations exceeds
#' max_iter, or when the errors converge to have a diff below the tolerance
#' @param rescale Should coefficients beta obtained by
#' covariance-regularized regression be re-scaled by a constant, given
#' by regressing $y$ onto $x beta$? This is done in Witten and
#' Tibshirani (2008) and is important for good performance. Default is
#' TRUE.
#' @param ddc_first whether or not to detect and impute cellwise outliers using
#' DDC before computing estimates
#' @param ddc_with_response whether or not to concatenate Y before running DDC. Only
#' applicable if `ddc_first` is `TRUE`
#' @param standardize whether or not to scale the training X and Y data
#' before performing estimation
#' @param lam1_init one of "random" or "max"
#' @return a list with two objects:
#' * mod: the scout object from fitting the full X and Y_test with
#'   scout(p1, 1) using the best lambdas from the cross-validation
#' * best_cv_res: a list with the results of the cross-validation fold with
#'   the lowest cross-validation error
#' @export
cv.scout_alternating_lasso <- function(X, Y, p1,
                                       nlambda1 = 100, nlambda2 = 100,
                                       lambda1.min.ratio = 0.1, lambda2.min.ratio = 0.001,
                                       K = 5, tol = 1e-4, max_iter = 10,
                                       rescale = TRUE, standardize = TRUE,
                                       centerFun = mean, scaleFun = sd,
                                       ddc_first = FALSE, ddc_with_response = FALSE,
                                       lam1_init = "random", cov_method = "default") {
  if (ddc_first) {
    ddc_input <- if (ddc_with_response) cbind(X, Y) else X
    ddc_res <- cellWise::DDC(ddc_input, DDCpars = list(fastDDC = TRUE, silent = TRUE))
    ddc_imp <- ddc_res$Ximp

    p <- ncol(X)
    names_X <- colnames(X) # save old column names
    X <- ddc_imp[, 1:p]
    colnames(X) <- names_X
    if (ddc_with_response) {
      names_Y <- colnames(Y) # save old column names
      Y <- as.matrix(ddc_imp[, (p + 1)], ncol = 1)
      colnames(Y) <- names(Y)
    }
  }

  folds <- cv.folds(nrow(X), K)
  # cv.folds returns list of length K where each item is a vector containing
  # row indices of X for the validation set

  apply_res <- lapply(folds, function(validation_idx) {
    X_cv <- X[-validation_idx, ]
    Y_cv <- Y[-validation_idx, ]
    X_val <- X[validation_idx, ]
    Y_val <- Y[validation_idx, ]
    scout_alternating_lasso(
      X_cv, Y_cv, X_val, Y_val, p1,
      nlambda1, nlambda2, lambda1.min.ratio, lambda2.min.ratio,
      tol, max_iter, rescale, standardize, centerFun, scaleFun,
      lam1_init, cov_method
    )
  })

  cv_min_errors <- sapply(apply_res, function(scout_alt_res) tail(scout_alt_res$errors, 1))
  best_res_idx <- which.min(cv_min_errors)
  best_res <- apply_res[[best_res_idx]]
  best_lambda_pair <- unlist(tail(best_res$lambda_pairs, 1))

  list(
    mod = scout(
      X, Y,
      p1 = p1, p2 = 1,
      lam1s = best_lambda_pair[1], lam2s = best_lambda_pair[2],
      centerFun = mean, scaleFun = sd,
      cov_method = cov_method
    ),
    best_cv_res = best_res
  )
}
