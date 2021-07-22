# install With:
# devtools::install_github("erikcs/grf", subdir = "r-package/grf", ref = "RATE"))

test_that("rank_average_treatment_effect works as expected", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- tau * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  cf <- causal_forest(X, Y, W, num.trees = 250)

  prio <- get_scores(cf)
  rate <- rank_average_treatment_effect(cf, prio)
  q.length <- nrow(rate$TOC)
  expect_equal(rate$TOC[q.length, "estimate"], 0, tolerance = 1e-10) # Last TOC curve entry = zero.
  expect_equal(rate$TOC[q.length, "std.err"], 0, tolerance = 1e-10) # Last TOC curve entry = zero.
  expect_true(all(rate$TOC$estimate[1:(q.length - 1)] >= 0)) # for this prio all points on curve are > 0

  rate.lo <- rank_average_treatment_effect(cf, prio, subset = prio < 0)
  expect_gt(rate[["estimate"]], rate.lo[["estimate"]])

  rate.all.eq <- rank_average_treatment_effect(cf, rep(1, n))
  expect_equal(rate.all.eq[["estimate"]], 0, tolerance = 1e-10)

  rand.prio <- sample(1:100, n, TRUE)
  autoc.rand <- rank_average_treatment_effect(cf, rand.prio, R = 150)
  expect_equal(autoc.rand[["estimate"]], 0, tolerance = 3 * autoc.rand[["std.err"]])

  qini.rand <- rank_average_treatment_effect(cf, rand.prio, method = "QINI", R = 150)
  expect_equal(qini.rand[["estimate"]], 0, tolerance = 3 * qini.rand[["std.err"]])

  print(rate)
  plot(rate)
})

test_that("rank_average_treatment_effect agrees with plain brute-force calculation", {
  # ...
  # make a handful test cases with/wo ties where you do the straight forward
  # brute force calculation of all estimates and check they agree with
  # the `rank_average_treatment_effect` output


  expect_equal(1, 1)
})

test_that("cluster robust rank_average_treatment_effect is consistent", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- tau * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)

  Xc <- rbind(X, X, X, X, X)
  Wc <- c(W, W, W, W, W)
  Yc <- c(Y, Y, Y, Y, Y)
  clust <- rep(1:n, 5)

  cf <- causal_forest(X, Y, W, num.trees = 250)
  cf.clust <- causal_forest(Xc, Yc, Wc, clusters = clust, num.trees = 250)
  prio <- runif(n)
  prio.clust <- rep(prio, 5)

  autoc <- rank_average_treatment_effect(cf, prio, method = "AUTOC", R = 150)
  qini <- rank_average_treatment_effect(cf, prio, method = "QINI", R = 150)

  autoc.clust <- rank_average_treatment_effect(cf.clust, prio.clust, method = "AUTOC", R = 150)
  qini.clust <- rank_average_treatment_effect(cf.clust, prio.clust, method = "QINI", R = 150)

  expect_equal(autoc[["estimate"]], autoc.clust[["estimate"]], tolerance = 0.05)
  expect_equal(autoc[["std.err"]], autoc.clust[["std.err"]], tolerance = 0.02)

  expect_equal(qini[["estimate"]], qini.clust[["estimate"]], tolerance = 0.05)
  expect_equal(qini[["std.err"]], qini.clust[["std.err"]], tolerance = 0.02)
})
