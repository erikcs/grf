test_that("rank_average_treatment_effect works as expected", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- tau * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  cf <- causal_forest(X, Y, W, W.hat = 0.5, num.trees = 250)

  prio <- get_scores(cf)
  rate <- rank_average_treatment_effect(cf, prio)
  capture.output(print(rate))
  plot(rate)

  q.length <- nrow(rate$TOC)
  expect_equal(rate$TOC[q.length, "estimate"], 0, tolerance = 1e-10) # Last TOC curve entry (q=1) = zero.
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

  cf.survival <- causal_survival_forest(X, Y, W, rep(1, n), W.hat = 0.5, num.trees = 250)
  autoc.cfs <- rank_average_treatment_effect(cf.survival, rand.prio, R = 150)
  expect_equal(autoc.cfs[["estimate"]], 0, tolerance = 3 * autoc.cfs[["std.err"]])
})

test_that("TOC grid works as expected", {
  n <- 500
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- tau * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  cf <- causal_forest(X, Y, W, W.hat = 0.5, num.trees = 250)
  prio <- tau

  # Computing TOC on grid 1/n <= q <= 1 agrees exactly with AUTOC.
  # (Disable half-sample bootstrap (R=0) to be able to pass in such a granular grid)
  q.full <- seq(1/n, 1, by = 1/n)
  rate.full <-rank_average_treatment_effect(cf, prio, q = q.full, R = 0)
  autoc <- rate.full$estimate
  expect_equal(mean(rate.full$TOC$estimate), autoc, tolerance = 1e-10)

  # Can ask for reasonable coverage of the entire curve.
  rand.prio <- sample(1:100, n, TRUE)
  rate <- rank_average_treatment_effect(cf, rand.prio, R = 250)
  TOC <- rate$TOC$estimate
  TOC.se <- rate$TOC$std.err
  expect_equal(TOC, rep(0, length(TOC)), tolerance = 3 * TOC.se)

  q5 <- seq(0.05, 1, by = 0.05)
  rate.q5 <- rank_average_treatment_effect(cf, rand.prio, q = q5, R = 250)
  TOC.q5 <- rate.q5$TOC$estimate
  TOC.q5.se <- rate.q5$TOC$std.err
  expect_equal(TOC.q5, rep(0, length(TOC.q5)), tolerance = 3 * TOC.q5.se)
})

test_that("rank_average_treatment_effect agrees with plain brute-force calculation", {
  n <- 50
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  W <- rbinom(n, 1, 0.5)
  tau <- pmax(X[, 1], 0)
  Y <- tau * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
  cf <- causal_forest(X, Y, W, W.hat = 0.5, num.trees = 250)
  DR.scores <- get_scores(cf)

  # 1. Unique priorities
  prio <- runif(n)
  autoc <- rank_average_treatment_effect(cf, prio, R = 0)
  qini <- rank_average_treatment_effect(cf, prio, method = "QINI", R = 0)

  sort.idx <- order(prio, decreasing = TRUE)
  TOC <- rep(NA, n)
  for (i in 1:n) {
    TOC[i] <- mean(DR.scores[sort.idx[1:i]]) - mean(DR.scores)
  }
  AUTOC <- mean(TOC)
  QINI <- weighted.mean(TOC, 1:n)

  expect_equal(autoc[["estimate"]], AUTOC, tolerance = 1e-10)
  expect_equal(qini[["estimate"]], QINI, tolerance = 1e-10)

  # 2. Duplicate priorities
  prio.dup <- sample(1:20, n, replace = TRUE)
  autoc.dup <- rank_average_treatment_effect(cf, prio.dup, R = 0)
  qini.dup <- rank_average_treatment_effect(cf, prio.dup, method = "QINI", R = 0)

  # average the doubly robust scores within tied groups
  scores.by.prio <- split(DR.scores, prio.dup) # orders by prio in increasing order
  ties.count.by.prio <- lapply(scores.by.prio, length)
  mean.scores.by.prio <- lapply(scores.by.prio, mean)
  # scores in decreasing priority order
  scores.order <- rev(unlist(rep(mean.scores.by.prio, ties.count.by.prio)))

  TOC.dup <- rep(NA, n)
  for (i in 1:n) {
    TOC.dup[i] <- mean(scores.order[1:i]) - mean(DR.scores)
  }
  AUTOC.dup <- mean(TOC.dup)
  QINI.dup <- weighted.mean(TOC.dup, 1:n)

  expect_equal(autoc.dup[["estimate"]], AUTOC.dup, tolerance = 1e-10)
  expect_equal(qini.dup[["estimate"]], QINI.dup, tolerance = 1e-10)
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

  cf <- causal_forest(X, Y, W, W.hat = 0.5, num.trees = 250)
  cf.clust <- causal_forest(Xc, Yc, Wc, W.hat = 0.5, clusters = clust, num.trees = 250)
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

test_that("internal bootstrap function `boot` works as expected", {
  n <- 50
  mu <- 1 + rnorm(n)
  clust <- 1:n

  statistic <- function(data, indices) {
    mean(data[indices, 1])
  }

  lm <- lm(mu ~ 1)
  lmt <- lmtest::coeftest(lm, vcov = sandwich::vcovCL, type = "HC", cluster = clust)
  boot.output <- boot(data.frame(mu), statistic, R = 250, clusters = clust, half.sample = FALSE)
  expect_equal(sd(boot.output$t), lmt[1, "Std. Error"], tolerance = 0.02)

  mu.cl <- c(mu, mu, mu, mu)
  clust.cl <- rep(1:n, 4)

  lm.cl <- lm(mu.cl ~ 1)
  lmt.cl <- lmtest::coeftest(lm.cl, vcov = sandwich::vcovCL, type = "HC", cluster = clust.cl)
  boot.output.cl <- boot(data.frame(mu.cl), statistic, R = 250, clusters = clust.cl, half.sample = FALSE)
  expect_equal(sd(boot.output.cl$t), lmt.cl[1, "Std. Error"], tolerance = 0.02)
})
