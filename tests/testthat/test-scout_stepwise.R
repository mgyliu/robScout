make_test_xy_data <- function(n = 10, p = 10) {
  x <- MASS::mvrnorm(n, rep(0, p), diag(p) * 0.5)
  beta <- rep(1, p)
  y <- x %*% beta + rnorm(n, 0, 0.1)
  list(x = x, y = y)
}

test_that("scout_1something_stepwise works for p2 = NULL", {
  data <- make_test_xy_data()

  res <- scout_1something_stepwise(data$x, data$y, p2 = NULL, nlambda1 = 5, nlambda2 = 5, cov_method = "default", glasso_crit = "bic")
  expect_named(res, c("g.res", "cv.res", "mod"))
})

test_that("scout_1something_stepwise works for p2 = 1", {
  data <- make_test_xy_data()

  res <- scout_1something_stepwise(data$x, data$y, p2 = 1, K = 3, nlambda1 = 5, nlambda2 = 5, cov_method = "default", glasso_crit = "bic")
  expect_named(res, c("g.res", "cv.res", "mod"))
})

test_that("scout_1something_stepwise errors if p2 is not 1 or NULL", {
  data <- make_test_xy_data()
  expect_error(
    scout_1something_stepwise(data$x, data$y, p2 = 2, cov_method = "default", glasso_crit = "bic"),
    "not implemented"
  )
})
