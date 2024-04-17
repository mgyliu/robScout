test_that("cov_winsor works for cov(X) or cov(X,Y)", {
  p <- 5
  X <- MASS::mvrnorm(10, rep(0, p), diag(p))
  Y <- rnorm(10)

  expect_true("matrix" %in% class(cov_winsor(X)))
  expect_equal(dim(cov_winsor(X)), c(p, p))
  expect_true("matrix" %in% class(cov_winsor(X, Y)))
  expect_equal(length(cov_winsor(X, Y)), p)
})

test_that("est_cov works for all available options", {
  p <- 5
  X <- MASS::mvrnorm(10, rep(0, p), diag(p))
  Y <- rnorm(10)

  cov_default <- est_cov(X, method = "default")
  expect_true("matrix" %in% class(cov_default))
  expect_equal(dim(cov_default), c(p, p))

  cov_winsor <- est_cov(X, method = "winsor")
  expect_identical(cov_winsor, cov_winsor(X))
  expect_equal(dim(cov_winsor), c(p, p))

  cov_wrap <- est_cov(X, method = "wrap")
  expect_identical(cov_wrap, cov_wrap(X))
  expect_equal(dim(cov_wrap), c(p, p))

  cov_ddc <- est_cov(X, method = "ddc")
  expect_identical(cov_ddc, cov_ddc(X))
  expect_equal(dim(cov_ddc), c(p, p))

  expect_warning(est_cov(X, method = "not implemented"), "not implemented")
})

test_that("cov_winsor works when X is univariate", {
  X <- matrix(rnorm(10), ncol = 1)
  expect_equal(dim(est_cov(X, method = "default")), c(1, 1))
  expect_equal(dim(est_cov(X, method = "winsor")), c(1, 1))
})

test_that("cov_ddc warns if input data has too few columns", {
  X <- matrix(rnorm(10), ncol = 1)
  expect_warning(est <- est_cov(X, method = "ddc"))
  expect_equal(dim(est), c(1, 1))
})

test_that("est_cov works for all methods when computing cov(X,Y)", {
  p <- 5
  X <- MASS::mvrnorm(10, rep(0, p), diag(p))
  Y <- rnorm(10)

  cov_default <- est_cov(X, Y, method = "default")
  expect_true("matrix" %in% class(cov_default))
  expect_equal(dim(cov_default), c(p, 1))

  cov_winsor <- est_cov(X, Y, method = "winsor")
  expect_identical(cov_winsor, cov_winsor(X, Y))
  expect_equal(dim(cov_winsor), c(p, 1))

  cov_ddc <- est_cov(X, Y, method = "ddc")
  expect_equal(dim(cov_ddc), c(p, 1))

  cov_wrap <- est_cov(X, Y, method = "wrap")
  expect_identical(cov_wrap, cov_wrap(X, Y))
  expect_equal(dim(cov_wrap), c(p, 1))
})
