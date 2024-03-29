% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{perm.toscca}
\alias{perm.toscca}
\title{Permutation testing for toscca}
\usage{
perm.toscca(
  A,
  B,
  nonzero_a,
  nonzero_b,
  K,
  alpha_init = c("eigen", "random", "uniform"),
  folds = 1,
  toPlot = FALSE,
  draws = 20,
  cancor,
  bootCCA = NULL,
  silent = TRUE,
  parallel_logic = TRUE,
  nuisanceVar = 0,
  testStatType = "CC"
)
}
\arguments{
\item{A, B}{Data matrices.}

\item{nonzero_a, nonzero_b}{Numeric. Scalar or vector over the number of nonzeroes allowed for a correlation estimate.}

\item{K}{Numeric. Number of components to be computed.}

\item{folds}{Numeric. Number of folds for the cross-validation process.}

\item{toPlot}{Logical. If TRUE, plot will be generated automatically showing the estimated canonical weights. Default is TRUE.}

\item{draws}{Numeric. Number of permutations for each component.}

\item{cancor}{Numeric. Scalar or vector: anonical correlation estimate(s).}

\item{silent}{Logical. If FALSE, a progress bar will appear on the console. Default is FALSE.}

\item{parallel_logic}{Logical. If TRUE, cross-validation is done in parallel.Default is FALSE.}

\item{nuisanceVar}{Data with nuisance variables. For statistic type.}

\item{testStatType}{Character. Choice of statistic. Options are CC (default), Wilks and Roy.}

\item{combination}{Logical. If TRUE, the algorithm will search for the best combination of sparsity choice nonzero_a and nonzero_b for each component. This should be used for exploratory analysis. Default is FALSE.}
}
\value{
Matrix with permutation estimates.
}
\description{
This function performs permutation testing on CC estimates.
}
\details{
For a exploratory analysis nonzero_a and nonzero_b can be vectors. The algorithm will then search for the best combination of sparsity choice nonzero_a and nonzero_b for each component.
}
