% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi_regression_forest.R
\name{predict.multi_regression_forest}
\alias{predict.multi_regression_forest}
\title{Predict with a multi regression forest}
\usage{
\method{predict}{multi_regression_forest}(object, newdata = NULL, num.threads = NULL, drop = FALSE, ...)
}
\arguments{
\item{object}{The trained forest.}

\item{newdata}{Points at which predictions should be made. If NULL, makes out-of-bag
predictions on the training set instead (i.e., provides predictions at
Xi using only trees that did not use the i-th training example). Note
that this matrix should have the number of columns as the training
matrix, and that the columns must appear in the same order.}

\item{num.threads}{Number of threads used in prediction. If set to NULL, the software
automatically selects an appropriate amount.}

\item{drop}{If TRUE, coerce the prediction result to the lowest possible dimension. Default is FALSE.}

\item{...}{Additional arguments (currently ignored).}
}
\value{
A list containing `predictions`: a matrix of predictions for each outcome.
}
\description{
Gets estimates of E[Y_i | X = x] using a trained multi regression forest.
}
\examples{
\donttest{
# Train a standard regression forest.
n <- 500
p <- 5
X <- matrix(rnorm(n * p), n, p)
Y <-  X[, 1, drop = FALSE] \%*\% cbind(1, 2) + rnorm(n)
mr.forest <- multi_regression_forest(X, Y)

# Predict using the forest.
X.test <- matrix(0, 101, p)
X.test[, 1] <- seq(-2, 2, length.out = 101)
mr.pred <- predict(mr.forest, X.test)

# Predict on out-of-bag training samples.
mr.pred <- predict(mr.forest)
}

}
