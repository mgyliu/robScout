test_that("cov_winsor works for cov(X) or cov(X,Y)", {
  p <- 5
  X <- MASS::mvrnorm(10, rep(0, p), diag(p))
  Y <- rnorm(10)

  expect_true("matrix" %in% class(cov_winsor(X)))
  expect_equal(dim(cov_winsor(X)), c(p, p))
  expect_true("numeric" %in% class(cov_winsor(X, Y)))
  expect_equal(length(cov_winsor(X, Y)), p)
})

test_that("est_cov works", {
  p <- 5
  X <- MASS::mvrnorm(10, rep(0, p), diag(p))
  Y <- rnorm(10)

  cov_default <- est_cov(X, method = "default")
  expect_true("matrix" %in% class(cov_default))
  expect_equal(dim(cov_default), c(p, p))
  cov_winsor <- est_cov(X, method = "winsor")
  expect_identical(cov_winsor, cov_winsor(X))
  expect_warning(est_cov(X, method = "not implemented"), "not implemented")
})
