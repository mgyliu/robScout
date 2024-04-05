# A helper function to generate mock data for testing
make_test_data <- function() {
  n <- 10
  p <- 5
  rho <- 0.7
  mu <- rep(0, p)
  true_beta <- c(rep(3, 3), rep(0, p - 3))
  Sigma <- outer(1:p, 1:p,
    FUN = function(r, c) rho^(abs(r - c))
  )

  X_train <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  Y_train <- X_train %*% true_beta + rnorm(n, 0, 1)

  X_test <- MASS::mvrnorm(n = 20, mu = mu, Sigma = Sigma)
  Y_test <- X_test %*% true_beta + rnorm(20, 0, 1)

  meany <- mean(Y_train)
  meanx <- apply(X_train, 2, mean)
  sdy <- sd(Y_train)
  sdx <- apply(X_train, 2, sd)
  X_train_std <- scale(X_train, TRUE, TRUE)
  Y_train_std <- (Y_train - meany) / sdy

  list(
    meany = meany, meanx = meanx, sdy = sdy, sdx = sdx,
    X_train = X_train_std, Y_train = Y_train_std,
    X_test = X_test, Y_test = Y_test,
    true_beta = true_beta
  )
}

test_that("get_best_lam2 returns what you expect", {
  test_data <- make_test_data()
  hg <- huge::huge.glasso(cov(test_data$X_train), nlambda = 7, cov.output = TRUE, verbose = FALSE)
  lam1s <- hg$lambda
  lam2s <- lam1s
  cov_x_est <- hg$cov[[1]]

  res <- get_best_lam2_lasso(
    test_data$X_train, test_data$Y_train,
    test_data$X_test, test_data$Y_test,
    test_data$meanx, test_data$meany,
    test_data$sdx, test_data$sdy,
    cov_x_est,
    lam2s = lam2s, rescale = TRUE,
    cov_method = "default"
  )

  expect_equal(length(res$beta_hat), ncol(test_data$X_train))
  expected_names <- c("best_lam2", "best_lam2_idx", "beta_hat", "error", "intercept")
  expect_equal(sort(names(res)), expected_names)
})

test_that("get_best_lam1 returns what you expect", {
  test_data <- make_test_data()
  hg <- huge::huge.glasso(cov(test_data$X_train), nlambda = 7, cov.output = TRUE, verbose = FALSE)
  lam1s <- hg$lambda
  xtxs <- get_xtxs(test_data$X_train, lam1s, p1 = 1, cov_method = "default")
  lam2 <- lam1s[2]

  res <- get_best_lam1(
    test_data$X_train, test_data$Y_train,
    test_data$X_test, test_data$Y_test,
    test_data$meanx, test_data$meany,
    test_data$sdx, test_data$sdy,
    1, lam1s, xtxs, lam2,
    rescale = TRUE, cov_method = "default"
  )

  expect_equal(length(res$beta_hat), ncol(test_data$X_train))
  expected_names <- c("best_lam1", "best_lam1_idx", "best_xtx", "beta_hat", "error", "intercept")
  expect_equal(sort(names(res)), expected_names)
})

test_that("scout_alternating_lasso works for p1 = 1", {
  test_data <- make_test_data()
  X_train <- test_data$X_train
  Y_train <- test_data$Y_train
  X_test <- test_data$X_test
  Y_test <- test_data$Y_test

  out.random <- scout_alternating_lasso(X_train, Y_train, X_test, Y_test, p1 = 1, lam1_init = "random", cov_method = "default")
  out.max <- scout_alternating_lasso(X_train, Y_train, X_test, Y_test, p1 = 1, lam1_init = "max", cov_method = "default")
  expected_names <- c("errors", "betas", "intercepts", "lambda_pairs", "lam2_paths", "lam1s")
  expect_equal(sort(names(out.random)), sort(expected_names))
  expect_equal(sort(names(out.max)), sort(expected_names))
})

test_that("scout_alternating_lasso works for p1 = 2", {
  test_data <- make_test_data()
  X_train <- test_data$X_train
  Y_train <- test_data$Y_train
  X_test <- test_data$X_test
  Y_test <- test_data$Y_test

  out.random <- scout_alternating_lasso(X_train, Y_train, X_test, Y_test, p1 = 2, lam1_init = "random", cov_method = "default")
  out.max <- scout_alternating_lasso(X_train, Y_train, X_test, Y_test, p1 = 2, lam1_init = "max", cov_method = "default")
  expected_names <- c("errors", "betas", "intercepts", "lambda_pairs", "lam2_paths", "lam1s")
  expect_equal(sort(names(out.random)), sort(expected_names))
  expect_equal(sort(names(out.max)), sort(expected_names))
})


test_that("cv.scout_alternating_lasso works", {
  test_data <- make_test_data()
  X_train <- test_data$X_train
  Y_train <- test_data$Y_train
  # X_test <- test_data$X_test
  # Y_test <- test_data$Y_test

  cv_alt_res <- cv.scout_alternating_lasso(X_train, Y_train, p1 = 1, K = 2, lam1_init = "random", cov_method = "default")

  expected_names <- c("mod", "best_cv_res")
  expect_equal(sort(names(cv_alt_res)), sort(expected_names))
  expect_equal(class(cv_alt_res$mod), "scoutobject")
})
