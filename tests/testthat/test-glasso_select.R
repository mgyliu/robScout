# Helper function to generate data for tests
gen_cov <- function(p, rho = 0.5) {
  outer(1:p, 1:p, FUN = function(r, c) {
    rho^(abs(r - c))
  })
}

gen_data <- function(n = 5, p = 5) {
  beta <- 1:p
  Sigma <- gen_cov(p)
  X <- MASS::mvrnorm(n, rep(0, p), Sigma)
  Y <- X %*% beta + rnorm(n)
  list(X = X, Y = Y)
}

test_that("huge_glasso_lambda_seq produces a lambda path with the correct length and values", {
  p <- 5
  cov_X <- matrix(rnorm(p^2), nrow = 5)
  lambda.min.ratio <- 0.01
  nlambda <- 10

  test_seq <- huge_glasso_lambda_seq(cov_X, nlambda, lambda.min.ratio)
  # Length should be correct
  expect_equal(length(test_seq), nlambda)
  # Should be decreasing, and ratio should be correct
  expect_equal(test_seq[1] * lambda.min.ratio, tail(test_seq, 1))
})

test_that("icov_eval works", {
  p <- 5
  cov <- gen_cov(p)
  # Generate inverse of covariance with small amount of error
  icov <- solve(cov) + rnorm(p^2, 0, 0.001)
  expect_equal(length(icov_eval(icov, cov, 50, "bic")), 1)
})

test_that("glasso_select returns a list with the right items", {
  X <- gen_data()$X
  X.test <- gen_data(n = 10)$X
  gs_res <- glasso_select(X, X.test,
    cov_method = "default", crit = "ebic",
    nlambda = 5, lambda.min.ratio = 0.1
  )

  # it returns the named items we expect
  expect_true(setequal(
    names(gs_res),
    c("icovx", "covx", "best_lambda", "lambda", "errors")
  ))
  # it returns the correct best_lambda
  expect_equal(
    gs_res$errors[which(gs_res$best_lambda == gs_res$lambda)],
    min(gs_res$errors)
  )
})

test_that("glasso_cv works", {
  X <- gen_data(n = 10)$X
  expect_no_error(glasso_cv(X, 2, "default", "bic"))
})

test_that("glasso_cv should use provided lambda seq if it's not NULL", {
  X <- gen_data(n = 15)$X
  my_lambdas <- seq(0.1, 0.9, 0.1)
  cv_K <- 3
  res <- glasso_cv(X, cv_K, "default", "loglik", lambdas = my_lambdas)
  expect_identical(res$lambda_seq, my_lambdas)
  expect_true(all(dim(res$crit) == c(length(my_lambdas), cv_K)))
})

test_that("glasso_cv should use the provided cv folds if it's not NULL", {
  X <- gen_data(n = 10)$X
  my_folds <- list(c(1:5), c(6:10))
  expect_no_error(res <- glasso_cv(X, 2, "default", "loglik", folds = my_folds))
  expect_identical(res$cv.folds, my_folds)

  # without providing folds
  set.seed(2024)
  expected_folds <- cv.folds(10, 2)
  set.seed(2024)
  expect_no_error(res2 <- glasso_cv(X, 2, "default", "loglik"))
  expect_identical(res2$cv.folds, expected_folds)
})

test_that("glasso_cv returns the correct covariance estimate that was used for huge", {
  X <- gen_data(n = 10)$X
  res <- glasso_cv(X, 2, "winsor", "loglik")
  expect_identical(res$input.cov, est_cov(X, method = "winsor"))
})
