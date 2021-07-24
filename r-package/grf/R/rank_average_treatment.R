#' Estimate a Rank-Weighted Average Treatment Effect (RATE).
#'
#' Consider a rule S(Xi) assigning scores to units in decreasing order of treatment prioritization.
#' In the case of a forest with binary treatment, we provide estimates of the following, where
#' 1/n < q <= 1 represents the fraction of treated units:
#' \itemize{
#'   \item The Rank-Weighted Average Treatment Effect (RATE):
#'    \eqn{\int_{0}^{1} alpha(q) TOC(q; S) dq}, where alpha is a weighting method
#'    corresponding to either `AUTOC` (identity-weighting) or `QINI` (linear weighting).
#'   \item The Targeting Operating Characteristic (TOC):
#'     \eqn{E[Y(1) - Y(0) | F(S(Xi)) >= 1 - q] - E[Y(1) - Y(0)]}, where F(.) is the distribution function of S(Xi).
#' }
#' The Targeting Operating Characteristic (TOC) is a curve comparing the benefit of treating only a certain
#' fraction q of units (as prioritized by S(Xi)), to the overall average treatment effect.
#' The Rank-Weighted Average Treatment Effect (RATE) is a weighted sum of this curve,
#' and is a measure designed to identify prioritization rules that effectively targets treatment.
#'
#' TODO This paragraph can contain details. Etc etc. Breifly say ties are fine. Briefly
#' say when you want AUTOC and when you want QINI.
#'
#' @param forest The trained forest.
#' @param priorities A vector of treatment prioritorization scores S(Xi) for the units used to train the forest.
#'  WARNING: for valid statistical performance, these scores should be constructed independently from the forest
#'  training data.
#' @param method The type of RATE estimate, options are `AUTOC` or `QINI`, corresponding to
#'  identity or linear weighting. Default is `AUTOC`.
#' @param q The grid q to compute the TOC curve on. Default is
#'  (10\%, 20\%, ..., 100\%).
#' @param R Number of bootstrap replicates for SEs. Default is 150.
#' @param subset Specifies subset of the training examples over which we
#'               estimate the ATE. WARNING: For valid statistical performance,
#'               the subset should be defined only using features Xi, not using
#'               the treatment Wi or the outcome Yi.
#' @param debiasing.weights A vector of length n (or the subset length) of debiasing weights.
#'               If NULL (default) these are obtained via the appropriate doubly robust score
#'               construction, e.g., in the case of causal_forests with a binary treatment, they
#'               are obtained via inverse-propensity weighting.
#' @param compliance.score Only used with instrumental forests. An estimate of the causal
#'               effect of Z on W, i.e., Delta(X) = E[W | X, Z = 1] - E[W | X, Z = 0],
#'               which can then be used to produce debiasing.weights. If not provided,
#'               this is estimated via an auxiliary causal forest.
#' @param num.trees.for.weights In some cases (e.g., with causal forests with a continuous
#'               treatment), we need to train auxiliary forests to learn debiasing weights.
#'               This is the number of trees used for this task. Note: this argument is only
#'               used when debiasing.weights = NULL.
#'
#' @references TODO a ref to your paper Athey, Susan, and Stefan Wager. "Policy Learning With Observational Data."
#'             Econometrica 89.1 (2021): 133-161.
#'
#' @examples
#' \donttest{
#' # Train a causal forest to estimate a priority ranking.
#' n <- 1500
#' p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' Y <- pmax(X[, 1], 0) * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
#' train <- sample(1:n, n / 2)
#' cf.priority <- causal_forest(X[train, ], Y[train], W[train])
#'
#' # Compute a prioritorization based on estimated effect quantiles.
#' tau.hats <- predict(cf.priority, X[-train, ])$predictions
#' priority <- cut(tau.hats, breaks = quantile(tau.hats), include.lowest = TRUE)
#'
#' # Estimate RATE on held out data.
#' cf <- causal_forest(X[-train, ], Y[-train], W[-train])
#' rate <- rank_average_treatment_effect(cf, priority)
#' rate
#'
#' # Plot the Targeting Operator Characteristic curve.
#' plot(rate)
#' }
#'
#' @return A list of class `rank_average_treatment_effect` with elements \itemize{
#'  \item estimate: the RATE estimate.
#'  \item std.err: bootstrapped standard error of RATE.
#'  \item TOC: a data.frame with the Targeting Operator Characteristic curve
#'    estimated on grid q, along with bootstrapped SEs.
#'  \item method: the type of RATE.
#' }
#' @export
rank_average_treatment_effect <- function(forest,
                                          priorities,
                                          method = c("AUTOC", "QINI"),
                                          q = seq(0.1, 1, by = 0.1),
                                          R = 150,
                                          subset = NULL,
                                          debiasing.weights = NULL,
                                          compliance.score = NULL,
                                          num.trees.for.weights = 500) {
  method <- match.arg(method)
  cluster.se <- length(forest$clusters) > 0
  clusters <- if (cluster.se) {
    forest$clusters
  } else {
    1:NROW(forest$Y.orig)
  }
  observation.weight <- observation_weights(forest)

  subset <- validate_subset(forest, subset)
  subset.clusters <- clusters[subset]
  subset.weights <- observation.weight[subset]

  if (length(unique(subset.clusters)) <= 1) {
    stop("The specified subset must contain units from more than one cluster.")
  }

  if (!is.null(debiasing.weights)) {
    if (length(debiasing.weights) == NROW(forest$Y.orig)) {
      debiasing.weights <- debiasing.weights[subset]
    } else if (length(debiasing.weights) != length(subset)) {
      stop("If specified, debiasing.weights must be a vector of length n or the subset length.")
    }
  }
  if (anyNA(priorities)) {
    stop("`priorities` contains missing values.")
  }
  if (length(priorities) == NROW(forest$Y.orig)) {
    # store as factor, more efficient `tabulate()` for large |S| with many ties.
    priorities <- as.factor(priorities[subset])
  } else if (length(priorities) != length(subset)) {
    stop("`priorities` must be a vector of length n or the subset length.")
  }

  if (!any(c("causal_forest", "instrumental_forest", "causal_survival_forest") %in% class(forest))) {
    stop("`rank_average_treatment_effect` is not implemented for this forest type.")
  }
  if (!all(forest$W.orig %in% c(0, 1))) {
    stop("Rank-weighted average treatment effect estimation only implemented for binary treatment.")
  }
  if (is.unsorted(q) || anyDuplicated(q) || min(q) <= 0 || max(q) != 1) {
    stop("`q` should correspond to a grid of fractions on the interval (0, 1].")
  }
  samples.by.cluster <- split(seq_along(subset.clusters), subset.clusters)
  n.half <- floor(length(samples.by.cluster) / 2)
  # lower bound on number of units in half-sample bootstrap. Equal to n.half when each unit its own cluster.
  smallest.bs.n <- sum(lengths(samples.by.cluster)[order(lengths(samples.by.cluster))][1:n.half])
  smallest.grid.bucket.size <- floor(smallest.bs.n * q)
  if (min(smallest.grid.bucket.size) == 0 || anyDuplicated(smallest.grid.bucket.size)) {
    stop(paste0("Provided `q` grid is too dense to uniquely assign each unit to a bucket ",
                "(or some clusters contains too few units)."))
  }

  # For all supported forest types, DR.scores is a n-length vector
  # (future support for multi_arm_causal_forest may be hairy...)
  DR.scores <- get_scores(forest, subset = subset, debiasing.weights = debiasing.weights,
                          compliance.score = compliance.score, num.trees.for.weights = num.trees.for.weights)

  # *** Compute the TOC and RATE ***

  # TODO: bake in sample weight ATE <- weighted.mean(DR.scores, subset.weights) sample.weights

  if (method == "AUTOC") {
    wtd.mean <- function(x) mean(x)
  } else if (method == "QINI") {
    wtd.mean <- function(x) weighted.mean(x, seq.int(1, length(x)))
  }

  # Compute estimates, a function to be passed on to boostrap routine.
  # @data: a data.frame with the original data. column 1: DR.scores, column 2: priority scores (factor)
  # @indices: a vector of indices which define the bootstrap sample.
  # @returns: an estimate of RATE, together with the TOC curve.
  estimate <- function(data, indices) {
    # let q be a fraction in (0, 1].
    # we have 1) TOC(q; Sj) = 1/[qn] sum_{i=1}^{[qn]} Gamma_{i(j)} - ATE
    # and 2) RATE = 1/n sum_{i=1}^{n} TOC(i/n; Sj)
    # For boostrapping the TOC curve, we fix q on some grid q'.
    # For estimating the RATE we set 1/n <= q <= 1.
    # Observations:
    # a) the entire TOC curve can be computed as a cumulative sum of sorted DR.scores.
    # b) taking ties into account amounts to using the average DR.scores within each
    # tied group instead of the individual DR.scores.
    # So the steps are:
    # Compute average DR.scores by priority group in increasing order: DR.scores',
    # repeat the possibly tied entries with the number of duplicates.
    # Take the cumsum of this divided by ni to get the TOC curve. Take a (weighted) average
    # of this to get RATE. This is what is being done below, using R's fastest primitives for
    # quick aggregation and grouping (using for example "rowsum"),
    nq <- floor(length(indices) * q)
    grid.id <- rep.int(nq, c(nq[1], diff(nq)))

    prio <- data[indices, 2]
    group.length <- tabulate(prio, nlevels(prio))
    group.length <- group.length[group.length != 0] # ignore potential levels not present in BS sample
    grp.means <- rowsum(data[indices, 1], as.integer(prio)) / group.length

    DR.scores.sorted <- rev(rep.int(grp.means, group.length))
    DR.scores.sorted.grid <- rowsum(DR.scores.sorted, grid.id)

    TOC <- cumsum(DR.scores.sorted) / seq_along(DR.scores.sorted) - mean(DR.scores.sorted)
    TOC.grid <- cumsum(DR.scores.sorted.grid) / nq - mean(DR.scores.sorted)

    RATE <- wtd.mean(TOC)
    c(RATE, TOC.grid)
  }
  boot.output <- boot(data.frame(DR.scores, priorities), estimate, R, subset.clusters, half.sample = TRUE)
  point.estimate <- boot.output[["t0"]]
  std.errors <- apply(boot.output[["t"]], 2, sd)
  # ensure invariance: always >= 0
  if (R < 2) {
    std.errors[] <- 0
  }

  output <- list()
  class(output) <- "rank_average_treatment_effect"
  output[["estimate"]] <- point.estimate[1]
  output[["std.err"]] <- std.errors[1]
  output[["TOC"]] <- data.frame(estimate = point.estimate[-1], std.err = std.errors[-1], q = q)
  output[["method"]] <- method

  output
}

#' Plot the Targeting Operator Characteristic curve.
#' @param x The output of rank_average_treatment_effect.
#' @param ... Additional arguments passed to plot.default.
#'
#' @method plot rank_average_treatment_effect
#' @export
plot.rank_average_treatment_effect <- function(x, ...) {
  n <- NROW(x[["TOC"]])
  TOC <- x[["TOC"]]
  q <- TOC[, "q"]
  ub <- TOC[, "estimate"] + 1.96 * TOC[, "std.err"]
  lb <- TOC[, "estimate"] - 1.96 * TOC[, "std.err"]
  plot(q, TOC[, "estimate"], type = "l", ylim = c(min(lb), max(ub)),
       ylab = "",
       main = "Targeting Operator Characteristic",
       sub = "(95 % confidence bars in dashed lines)",
       ...)
  abline(h = 0, lty = 3)
  lines(q, ub, col = "black", lty = 2)
  lines(q, lb, col = "black", lty = 2)
}

#' Print the Rank-Weighted Average Treatment Effect (RATE).
#' @param x The output of rank_average_treatment_effect.
#' @param ... Additional arguments (currently ignored).
#'
#' @method print rank_average_treatment_effect
#' @export
print.rank_average_treatment_effect <- function(x, ...) {
  print(c(estimate = x[["estimate"]], std.err = x[["std.err"]]))
}

#' Simple clustered bootstrap.
#'
#' Adopted from the `boot` function in the boostrap package with clusters added.
#' A future TODO could be to add parallel (not necessarily worth it)
#' https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf
#'
#' @param data A data frame with the original data.
#' @param statistic A function computing estimate(s) with signature (data, indices, ...) where
#' data is the original data, and indices a vector which defines the bootstrap sample.
#' @param R The number of bootstrap replications.
#' @param clusters Integer vector of cluster assignment, setting to 1:N corresponds to an ordinary
#'  unclustered bootstrap.
#' @param half.sample Whether to do half sample boostrap (half the clusters are drawn). Default is TRUE.
#' @param ... Additional arguments passed on to statistic.
#' @return A list with the original estimate t0, and bootstrap estimates t.
#'
#' @references Angelo Canty and Brian Ripley (2021). boot: Bootstrap R (S-Plus) Functions.
#'
#' @keywords internal
boot <- function(data, statistic, R, clusters, half.sample = TRUE, ...) {
  samples.by.cluster <- split(seq_along(clusters), clusters)
  n <- length(samples.by.cluster) # number of clusters
  if (n <= 1 || (half.sample && floor(n / 2) <= 1)) {
    stop("Cannot bootstrap sample with only one effective unit.")
  }
  if (half.sample) {
    n.bs <- floor(n / 2)
    index.list <- replicate(R, unlist(samples.by.cluster[sample.int(n, n.bs, replace = FALSE)], use.names = FALSE), simplify = FALSE)
  } else {
    index.list <- replicate(R, unlist(samples.by.cluster[sample.int(n, replace = TRUE)], use.names = FALSE), simplify = FALSE)
  }

  t0 <- statistic(data, seq_len(NROW(data)), ...)

  res <- lapply(seq_len(R), function(i) statistic(data, index.list[[i]], ...))
  t <- matrix(, R, length(t0))
  for (r in seq_len(R)) {
    t[r, ] <- res[[r]]
  }

  list(t0 = t0, t = t)
}
