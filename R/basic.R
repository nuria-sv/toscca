
#' Progress bar
#'
#'Shows progress of a process.
#'
#' @param end maximum number of times a process will run.
#' @param round current round
#' @return Display in consol of current status.
progressBar <- function(end, round) {
  width = round(50*0.75)
  speed = width/end
  step  = round*speed

  cat("|", rep(".", step), rep(" ", (width - step)), "|", round((round/end)*100), "%         \r")

}

modes = function(d) {
  i <- which(diff(sign(diff(d$y))) < 0 ) + 1
  data.frame(x = d$x[i], y = d$y[i])
}

numToBi = function(data) {
  r = nrow(data)
  c = ncol(data)

  logicalMatrix = matrix(data != 0, nrow = r, ncol = c)
  binaryMatrix  = matrix(as.numeric(logicalMatrix), nrow = r, ncol = c)

  binaryMatrix
}


getWhich = function(data, fun) {
  fun = match.fun(fun)

  position = (which(data == fun(data)))

  position
}


triangleArea = function(x, y, z = 0, type = c("rectangle")) { #assume a vertice is origin

  h = x - z
  b = y - z
  tArea = (b*h)/2

  tArea


}
#' Stardardise a matrix
#'
#' This function stardardises a matrix or a vector and gives the option to
#' centre or normalise (only vectors).
#'
#' @param X Matrix or vector to be standardise.
#' @param centre Logical, if true, cetre to mean zero.
#' @param normalise Logical, if true, performs vector normalisation.
#' @return A matrix or vector with the preferred standardarisation
#' @export

standardVar = function(mat, centre = TRUE, normalise = FALSE) {
  if(isFALSE(normalise)) {
    if(length(dim(mat)) == 2) {
      XMat = matrix(mat, nrow = nrow(mat), ncol = ncol(mat))

    }
  } else {
    XMat = matrix(mat, nrow = length(mat), ncol = 1)

  }

  if(centre) {
    center <- colMeans(XMat, na.rm=TRUE)
    XMat <- sweep(XMat, 2L, center)

  }

  dist.f = function(A) {
    a = A[!is.na(A)]
    a = sqrt(sum(a^2 / max(1, length(a) - 1L)))

  }


  dist = apply(XMat, 2L, dist.f)
  Xstan = sweep(XMat, 2L, dist, "/")


  if(isTRUE(normalise)) {
    vec <- matrix(XMat, ncol = 1)
    a <- sqrt(sum(vec^2))
    if(a==0) a <- .05
    Xstan <-  matrix(vec/a, ncol=1)

  }


  Xstan
}


#' Plot heatmap
#'
#' This function generated a heatmap. Withing the package it is used to
#' provide with a visualisation of the exploration of optimal sparsity levels.
#'
#' @param mat Matrix containing values along a grid.
#' @param palette Choice of palette. Default is Teal.
#' @param coln Number of columns for the grid distribution. Default is 12.
#' @param xlab Lable for X axis.
#' @param ylab Lable for Y axis.
#' @param axes Logical. Have axes between 0 and 1. Default is FALSE.
#' @return Grid plot.
#' @export
myHeatmap = function (mat, palette = "Teal", coln = 12, xlab = "", ylab = "", axes = FALSE)
{
  par(fig = c(0, 7/10, 0, 1))
  image(mat, col = hcl.colors(length(mat), palette, rev = TRUE), axes = axes, xlab = xlab, ylab = ylab)
  if(isFALSE(axes)) {
    axis(1, seq(0, 1, length = nrow(mat)), rownames(mat))
    axis(2, seq(0, 1, length = ncol(mat)), colnames(mat), las = 1)
  }

  # if (as == "i") {
  #   axis(1, seq(0, 1, length = nrow(mat)), c(1:nrow(mat)),
  #        tck = FALSE)
  #   axis(2, seq(0, 1, length = ncol(mat)), c(1:ncol(mat)),
  #        tck = FALSE)
  # }

  par(fig = c(7/10, 1, 0, 1), new = TRUE)
  colvec = matrix(seq(min(mat), max(mat), length = coln))
  image(t(colvec), col = hcl.colors(12, palette, rev = TRUE), axes = FALSE)
  axis(2, seq(0.0625, 1, length = 10/2), format(colvec[seq(2,
                                                           10, 2)], digits = 3, nsmall = 2), las = 1)
  par(fig = c(0, 1, 0, 1))
}



#' Performs power method.
#'
#'
#' @param mat A square matrix nxn.
#' @param vec A vector of dimensions nx1.
#' @param tol Convergence criterion. Default is 10^(-6).
#' @param matIter Maximum iterations. Default is 500.
#' @param silent Logical. If TRUE, convergence performance will be printed.
#' @return List: vec: eigen vector; lambda: eigen value; t: total iterations.
#' @export
powerMethod = function(mat, vec, tol = 10^(-6), maxIter = 500, silent = TRUE) {
  vec = matrix((vec/norm(v, "2")), ncol = 1)
  e = 10
  t = 1

  while (e > tol | t >= maxIter) {

    vNew = mat%*%vec
    vNew = vNew/norm(vNew, "2")

    lambda = as.numeric(t(vNew)%*%mat%*%vNew)

    # update
    vec = vNew
    t = t+1

    e = norm(mat%*%vec - lambda*vec, "2")

    if(isFALSE(silent)) cat("Iterations ", t - 1, "with error: ", e, "      \r")
  }

  return(list(v = vec,
              lambda = lambda,
              iterations = t))
}
