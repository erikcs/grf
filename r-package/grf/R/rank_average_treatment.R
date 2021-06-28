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
#'  WARNING: for valid statistical performance, these scores should be obtained independently from the training
#'  data.
#' @param method The type of RATE estimate, options are `AUTOC` or `QINI`, corresponding to
#'  identity or linear weighting. Default is `AUTOC`.
#' @param R Optional number of bootstrap replicates for SEs (default is 0, no SEs. Requires)
#'  the optional `boot` library. todo
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
#' # Train a causal forest.
#' n <- 500
#' p <- 5
#' X <- matrix(rnorm(n * p), n, p)
#' W <- rbinom(n, 1, 0.5)
#' tau <- pmax(X[, 1], 0)
#' Y <- tau * W + X[, 2] + pmin(X[, 3], 0) + rnorm(n)
#' c.forest <- causal_forest(X, Y, W)
#' S <- tau
#'
#' # Estimate the RATE with AUTOC weighting.
#' rate <- rank_average_treatment_effect(c.forest, S)
#' # Print the RATE estimate.
#' rate
#' # Plot the Targeting Operator Characteristic.
#' plot(rate)
#' }
#'
#' @return A list with elements `estimate`: the RATE estimate, `std.err`: boostrapped standard error
#'  of RATE, `TOC`: the Targeting Operator Characteristic curve, along with bootstrapped standard errors.
#' @export
rank_average_treatment_effect <- function(forest,
                                          priorities,
                                          method = c("AUTOC", "QINI"),
                                          R = 5,
                                          subset = NULL,
                                          debiasing.weights = NULL,
                                          compliance.score = NULL,
                                          num.trees.for.weights = 500) {
  # TODO: ALPHA in docstring is an unfortunate name collision with GRF forest weights...
  # TODO: allow a vector of scores instead of forest as well?
  # TODO: could estimate both AUTOC and QINI at the same time at practically zero cost...?
  method <- match.arg(method)
  if (R > 0 && !requireNamespace('boot', quietly = TRUE)) {
    stop("Bootstrapped standard error estimates requires the `boot` library.")
  }
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
  if (length(priorities) == NROW(forest$Y.orig)) {
    priorities <- as.factor(priorities[subset]) # store as factor, more efficient for large |S| with many ties.
  } else if (length(priorities) != length(subset)) {
    stop("`priorities` must be a vector of length n or the subset length.")
  }

  if (!any(c("causal_forest", "instrumental_forest", "causal_survival_forest") %in% class(forest))) {
    stop("`rank_average_treatment_effect` is not implemented for this forest type.")
  }
  if (!all(forest$W.orig %in% c(0, 1))) {
    stop("Rank-weighted average treatment effect estimation only implemented for binary treatment.")
  }

  # For all supported forest types, DR.scores is a n-length vector
  # (future support for multi_arm_causal_forest may be hairy...)
  DR.scores <- get_scores(forest, subset = subset, debiasing.weights = debiasing.weights,
                          compliance.score = compliance.score, num.trees.for.weights = num.trees.for.weights)

  ATE <- weighted.mean(DR.scores, subset.weights)
  n <- length(priorities)

  # *** Compute the TOC and RATE ***

  # TODO: bake in sample weight
  # TODO: complete rest of this
  # TODO: is the bootsrapped TOC correct?

  if (method == "AUTOC") {
    alpha <- rep(1, n)
  } else {
    alpha <- seq.int(1, n)
  }

  # Compute estimates, a function to be passed on to boostrap routine.
  # @data: a data.frame with the original data, DR.scores and priority scores.
  # @indices: a vector of indices which define the bootstrap sample.
  # @returns: an estimate of RATE, together with the TOC curve.
  estimate <- function(data, indices) {
    prio <- data[indices, 2]
    group.length <- tabulate(prio, nlevels(prio))
    group.length <- group.length[group.length != 0] # ignore potential levels not present in BS sample
    grp.means <- rowsum(data[indices, 1], as.integer(prio)) / group.length
    DR.scores.sorted <- rev(grp.means)
    TOC <- cumsum(DR.scores.sorted) / seq_along(DR.scores.sorted) - ATE
    TOC <- rep.int(TOC, rev(group.length))

    RATE <- weighted.mean(TOC, alpha)
    c(RATE, TOC)
  }
  # TODO: write custom bootstrap function, `boot` doesnt do clustering.
  # should use parallell: https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf
  # boot.output <- boot::boot(data.frame(DR.scores, priorities), estimate, R)
  boot.output <- boot(data.frame(DR.scores, priorities), estimate, R, half.sample = F)
  point.estimate <- boot.output[["t0"]]
  std.errors <- apply(boot.output[["t"]], 2, sd) # ensure invariance: should always be >= 0.

  output <- list()
  class(output) <- "rank_average_treatment_effect"
  output[["estimate"]] <- point.estimate[1]
  output[["std.err"]] <- std.errors[1]
  output[["TOC"]] <- data.frame(estimate = point.estimate[-1], std.err = std.errors[-1])
  output[["method"]] <- method

  output
}


# TODO: move to plot.R when done

#' Plot the Targeting Operator Characteristic curve.
#' @param x The output of rank_average_treatment_effect.
#' @param ... Additional arguments passed to plot.default.
#'
#' @method plot rank_average_treatment_effect
#' @export
plot.rank_average_treatment_effect <- function(x, ...) {
  n <- NROW(x[["TOC"]])
  q <- seq(1 / n, 1, by = 1 / n)
  TOC <- x[["TOC"]]
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

# TODO: move to print.R when done

#' Print the Rank-Weighted Average Treatment Effect (RATE).
#' @param x The output of rank_average_treatment_effect.
#' @param ... Additional arguments (currently ignored).
#'
#' @method print rank_average_treatment_effect
#' @export
print.rank_average_treatment_effect <- function(x, ...) {
  print(c(estimate = x[["estimate"]], std.err = x[["std.err"]]))
}

#' Simple clustered bootstrap
#' Adopted from Bryan Ripley's `boot` boostrap package
#' @param x The output of rank_average_treatment_effect.
#' @param x The output of rank_average_treatment_effect.
#' @param x The output of rank_average_treatment_effect.
#' @param ... Additional arguments (currently ignored).
#'
#' @references Angelo Canty and Brian Ripley (2021). boot: Bootstrap R (S-Plus) Functions. R package version 1.3-28.
#' @keywords internal
boot <- function(data, statistic, R, clusters, half.sample = TRUE) {
  n <- NROW(data)
  if (n <= 1) {
    stop("Cannot bootstrap sample of dim 1.")
  }
  if (half.sample) {
    n.bs <- floor(n / 2)
  } else {
    n.bs <- n
  }
  index.array <- sample.int(n, n.bs * R, replace = TRUE)
  # check: fortran order speed?
  dim(index.array) <- c(R, n.bs)
  # dim(index.array) <- c(n.bs, R)

  t0 <- statistic(data, seq_len(n))

  res <- lapply(seq_len(R), function(i) statistic(data, index.array[i, ]))
  # res <- lapply(seq_len(R), function(i) statistic(data, index.array[, i]))
  t <- matrix(, R, length(t0))
  for (r in seq_len(R)) {
    t[r, ] <- res[[r]]
  }

  list(t0 = t0, t = t)
}
