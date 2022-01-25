#' Simulate survival data
#'
#' Survival benchmarks for estimating the conditional expected survival time E[T | X].
#' The DGPs are scenarios 1 to 4 from Zhu, R., & Kosorok, M. R. (2012). Recursively imputed survival trees.
#' Journal of the American Statistical Association, 107(497), 331-340.
#'
#' @param n The number of samples.
#' @param Y.max The maximum follow-up time (optional).
#' @param n.mc The number of monte carlo draws to estimate the mean with. Default is 10000.
#' @param dgp The type of DGP.
#'
#' @return A list with entries:
#'  `X`: the covariates, `Y`: the event times, `D`: the censoring indicator,
#'  `ET`: the expected survival time E[T | X] estimated by monte carlo,
#'  `dgp`: the dgp name, `Y.max`: the maximum follow-up time.
#'
#' @examples
#' \donttest{
#' # Generate data
#' n <- 1000
#' data <- generate_survival_data(n, dgp = "ZK1")
#' }
#'
#' @export
generate_survival_data <- function(n, Y.max = NULL, n.mc = 10000,
                                   dgp = c("ZK1", "ZK2", "ZK3", "ZK4")) {
  dgp <- match.arg(dgp)
  if (dgp != "ZK2") {
    if (!("MASS" %in% utils::installed.packages())) {
      stop("This DGP requires the MASS library.")
    }
  }
  if (dgp == "ZK1") {
    p <- 25
    if (is.null(Y.max)) {
      Y.max <- 4
    }
    X <- pnorm(MASS::mvrnorm(n, rep(0, p), stats::toeplitz(0.9^seq(0, p - 1))))
    muT <- 0.1 * rowSums(X[, 11:20])
    muC <- mean(muT) / 2
    ft <- rexp(n, muT)
    failure.time <- pmin(ft, Y.max)
    censor.time <- rexp(n, muC)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    ET <- rep(NA, n)
    for (i in 1:n) {
      ft <- rpois(n.mc, muT[i])
      ET[i] <- mean(pmin(ft, Y.max))
    }
  } else if (dgp == "ZK2") {
    p <- 10
    if (is.null(Y.max)) {
      Y.max <- 6
    }
    X <- matrix(runif(n * p), n, p)
    muT <- sin(0.5 * X[, 1]) + 2*abs(X[, 2] - 0.5) + X[, 3]^3
    ft <- rexp(n, muT)
    failure.time <- pmin(ft, Y.max)
    censor.time <- runif(n, 0, Y.max)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    ET <- rep(NA, n)
    for (i in 1:n) {
      ft <- rpois(n.mc, muT[i])
      ET[i] <- mean(pmin(ft, Y.max))
    }
  } else if (dgp == "ZK3") {
    p <- 25
    if (is.null(Y.max)) {
      Y.max <- 10
    }
    X <- pnorm(MASS::mvrnorm(n, rep(0, p), stats::toeplitz(0.75^seq(0, p - 1))))
    muT <- 0.5 + 0.3 * abs(rowSums(X[, 11:15]))
    ft <- rgamma(n, shape = muT, scale = 2)
    failure.time <- pmin(ft, Y.max)
    censor.time <- runif(n, 0, 1.5 * Y.max)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    ET <- rep(NA, n)
    for (i in 1:n) {
      ft <- rgamma(n.mc, shape = muT[i], scale = 2)
      ET[i] <- mean(pmin(ft, Y.max))
    }
  } else if (dgp == "ZK4") {
    p <- 25
    if (is.null(Y.max)) {
      Y.max <- 4
    }
    X <- pnorm(MASS::mvrnorm(n, rep(0, p), stats::toeplitz(0.75^seq(0, p - 1))))
    muT <- 0.1 * abs(rowSums(X[, 1:5])) + 0.1 * abs(rowSums(X[, 21:25]))
    ft <- rlnorm(n, muT)
    failure.time <- pmin(ft, Y.max)
    censor.time <- rlnorm(n, muT + 0.5)
    Y <- pmin(failure.time, censor.time)
    D <- as.integer(failure.time <= censor.time)
    ET <- rep(NA, n)
    for (i in 1:n) {
      ft <- rlnorm(n.mc, muT[i])
      ET[i] <- mean(pmin(ft, Y.max))
    }
  }

  list(X = X, Y = Y, D = D, ET = ET, dgp = dgp, Y.max = Y.max)
}
