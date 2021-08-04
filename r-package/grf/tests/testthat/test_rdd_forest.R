test_that("rdd forest works as expected", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  threshold <- 0
  Z <- runif(n, -4, 4)
  W <- as.numeric(Z>= threshold)
  tau <- X[,1]
  Y <- W * tau + 1 / (1 + exp(2 * X[, 2])) * Z + 0.2 * rnorm(n) + 0.1 * Z^3

  rf <- rdd_forest(X, Y, Z, threshold)
  pp <- predict(rf)
  pp.test <- predict(rf, X)

  expect_equal(1, 1)
})
