test_that("run_ddc works and returns imputed data", {
  X <- MASS::mvrnorm(Sigma = diag(5), mu = rep(0, 5), n = 20)
  # X only
  res <- run_ddc(X)
  expect_false(all(X == res$X))
  expect_identical(dim(X), dim(res$X))
  expect_identical(colnames(X), colnames(res$X))
  expect_true(is.na(res$Y))

  # X with Y
  Y <- X %*% rep(1, 5) + rnorm(n = 20)
  res_with_Y <- run_ddc(X, Y)
  expect_false(all(X == res_with_Y$X))
  expect_identical(dim(X), dim(res_with_Y$X))
  expect_identical(colnames(X), colnames(res_with_Y$X))

  expect_false(all(is.na(res_with_Y$Y)))
  expect_identical(dim(Y), dim(res_with_Y$Y))
  expect_identical(colnames(Y), colnames(res_with_Y$Y))
})
