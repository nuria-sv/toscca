% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{getCanSubspace}
\alias{getCanSubspace}
\title{Performs matrix residualisation over estimated canonical vectors. There are
three types: basic (subtracts scaled estimated latent variable from data),
null (uses the null space of the estimated canonical vector to construct a
new matrix) and LV (uses SVD to residualise).}
\usage{
getCanSubspace(mat, vec)
}
\arguments{
\item{mat}{An nxp matrix.}

\item{vec}{A vector of dimensions nxk.}
}
\value{
An nxk matrix.
}
\description{
Performs matrix residualisation over estimated canonical vectors. There are
three types: basic (subtracts scaled estimated latent variable from data),
null (uses the null space of the estimated canonical vector to construct a
new matrix) and LV (uses SVD to residualise).
}
\details{
For nxp matrix \deqn{\mathbf{A}} and pxk vector \deqn{\mathbf{\alpha}}, the canonical
is compute as \eqn{\mathbf{A}_{sub} = \mathbf{A}\mathbf{\alpha}(\mathbf{\alpha}^T\mathbf{\alpha})\mathbf{\alpha}^T}.
}
