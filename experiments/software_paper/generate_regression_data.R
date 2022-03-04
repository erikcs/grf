#' Simulate regression data
#'
#' Regression benchmarks based on simulations DGPs in grf::generate_causal_data as well as
#'
#' 1)
#' 2)
#'
#' @param n The number of samples.
#' @param p The number of covariates.
#' @param dgp The type of DGP.
#'
#' @return A list with entries:
#'  `X`: the covariates, `Y`: the event times, `D`: the censoring indicator,
#'  `ET`: the expected survival time E[T | X] estimated by monte carlo,
#'  `dgp`: the dgp name, `Y.max`: the maximum follow-up time.
generate_regression_data <- function(n, p,
                                     dgp = c("kunzel", "ai1", "ai2", "NSLM", "AngristEvans")) {
  dgp <- match.arg(dgp)
  if (dgp %in% c("kunzel", "ai1", "ai2")) {
    return (grf::generate_causal_data(n, p, dgp = dgp))
  }
  1

  list(X = X, Y = Y, D = D, ET = ET, dgp = dgp, Y.max = Y.max)
}
