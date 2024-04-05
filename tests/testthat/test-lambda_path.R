test_that("get_lambda1_path runs and returns correct number of values for different p1 values", {
  test_X <- matrix(rnorm(50), nrow = 10)
  cov_x <- est_cov(test_X, method = "default")

  # p1 == NULL or 0
  expect_equal(get_lambda1_path(cov_x, NULL), 0)
  expect_equal(get_lambda1_path(cov_x, 0), 0)

  # p2 == 1
  expect_equal(length(get_lambda1_path(cov_x, 1, nlambda = 10)), 10)

  # p2 == 2
  expect_equal(length(get_lambda1_path(cov_x, 2, nlambda = 10)), 10)

  # any other p2
  expect_warning({
    lambdas <- get_lambda1_path(cov_x, 3, nlambda = 10)
  })
  expect_equal(length(lambdas), 10)
})

test_that("get_lambda2_path runs and returns correct number of values for different p2 values", {
  test_X <- matrix(rnorm(50), nrow = 10)
  test_Y <- rnorm(10)

  cov_x <- est_cov(test_X, method = "default")
  cov_xy <- est_cov(test_X, test_Y, method = "default")

  # p2 == NULL
  expect_equal(get_lambda2_path(test_X, test_Y, cov_x, cov_xy, NULL), 0)

  # p2 == 1
  expect_equal(length(get_lambda2_path(test_X, test_Y, cov_x, cov_xy, 1, nlambda = 10)), 10)

  # p1 == 2
  expect_warning({
    lambdas <- get_lambda2_path(test_X, test_Y, cov_x, cov_xy, 2, nlambda = 10)
  })
  expect_equal(length(lambdas), 10)
})

test_that("is_off_diag does the right thing", {
  expect_true(is_off_diag_zero(diag(10)))
  expect_false(is_off_diag_zero(matrix(c(rep(3, 5), rep(0, 20)), nrow = 5)))
})

test_that("get_lambda1_max_glasso runs and returns a numeric", {
  test_X <- matrix(rnorm(50), nrow = 10)
  expect_equal(
    class(get_lambda1_max_glasso(test_X)), "numeric"
  )
})

test_that("get_lambda2_max_lasso runs and returns a numeric", {
  test_X <- matrix(rnorm(50), nrow = 10)
  test_Y <- rnorm(10)
  expect_equal(
    class(get_lambda2_max_lasso(test_X, test_Y, est_cov(test_X, method = "default"), est_cov(test_X, test_Y, method = "default"))),
    "numeric"
  )
})

test_that("get_lambda_path runs and returns correct number of values", {
  expect_equal(length(get_lambda_path(5, 5, 0.01)), 5)
})
