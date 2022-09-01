#' Performs eigen decomposition of a matrix in PS space.
#'
#'
#' @param A A square matrix nxn.
#' @return Matrix. Positive definite matrix.
eigenDecompostion = function(A) {
  if(nrow(A) != ncol(A)) stop("Matrix should be square.")

  e = eigen(A)
  E = e$vector%*%diag(sqrt(pmax(0,e$values)))%*%t(e$vectors) # make positive definite matrix
  # transpose of eigenvectors is equat to inversy => symmetry

  return(E)
}

#' Initialised the canonical vector for the iterative process based on
#' positive eigen values. Then, SVD is performed on that PS matrix.
#'
#'
#' @param A An nxp matrix.
#' @param B An nxq matrix.
#' @return An pzp vector.
#' @export
initialiseCanVar = function(A, B) {
  bbt = B%*%t(B) # positive definite matrix
  # A%*%t(A) is symmetric
  E   = eigenDecompostion(bbt)
  H   = t(A)%*%E

  H_sdv  = svd(H)
  alpha0 = H_sdv$u

  return(alpha0)

}

#' Performs scalling for matrix residualisation based on calculated coefficients.
#'
#'
#' @param A An nxp matrix.
#' @return scaled matrix.
scaledResidualMat = function(A) {
  invA = solve(t(A)%*%A)
  resA = A%*%invA%*%t(A)

  return(resA)
}

#' Performs matrix residualisation over estimated canonical vectors. There are
#' three types: basic (subtracts scaled estimated latent variable from data),
#' null (uses the null space of the estimated canonical vector to construct a
#' new matrix) and LV (uses SVD to residualise).
#'
#'
#' @param mat An nxp matrix.
#' @param vec A vector of dimensions nxk.
#' @param spaceMat Only for "null" type residualisation. Default is NULL.
#' @param type Character. It can be LV, null or basic depending on which type of residualisation will be performed.
#' @param na.allow Logical. If TRUE, NAs will be allowed.
#' @return Matrix.
#' @export
residualisation = function(mat, vec, spaceMat = NULL, type = c("LV", "null", "basic"), na.allow = TRUE) {

  if(type == "basic") {
    lv     = mat%*%vec
    lv_res = scaledResidualMat(lv)

    mat_res = mat - lv_res%*%mat

  }
  if(type == "null") {
    require(mcompanion)
    vec_null  = mcompanion::null_complement(vec, universe = spaceMat, na.allow = na.allow)
    vec_com   = cbind(vec, vec_null)
    mat_res = mat%*%vec_com

  }

  if(type == "LV") {
    lv     = mat%*%vec
    lv_svd = svd(lv)
    s      = lv_svd$d
    s_inv  = s; s_inv[s_inv!=0] = 1/s_inv[s_inv!=0] # using solve can bring issues with zeroes
    s_inv  = matrix(s_inv*diag(length(s_inv)), ncol = length(s_inv))
    lv_inv = t(lv_svd$u%*%s_inv%*%t(lv_svd$v))

    R = diag(length(lv)) - lv%*%lv_inv

    mat_res = R%*%mat

  }

  return(mat_res)

}

#' Calculated cummulative percentage of explained variance.
#'
#' @param mat An nxp matrix.
#' @param maxK An nxk matrix. Each column corresponds to a latent variable.
#' @return Scalar.
#' @export
cpev.fun = function(mat, matK) {
  var_mat = t(mat)%*%mat
  traceX  = sum(rowSums(var_mat))

  var_matK = t(matK)%*%matK
  traceK   = sum(rowSums(var_matK))

  cpev = traceK/traceX

  return(cpev)

}
#' Performs matrix residualisation over estimated canonical vectors. There are
#' three types: basic (subtracts scaled estimated latent variable from data),
#' null (uses the null space of the estimated canonical vector to construct a
#' new matrix) and LV (uses SVD to residualise).
#'
#' @param mat An nxp matrix.
#' @param vec A vector of dimensions nxk.
#' @details
#' For nxp matrix \deqn{\mathbf{A}} and pxk vector \deqn{\mathbf{\alpha}}, the canonical
#'  is compute as \eqn{\mathbf{A}_{sub} = \mathbf{A}\mathbf{\alpha}(\mathbf{\alpha}^T\mathbf{\alpha})\mathbf{\alpha}^T}.
#' @return An nxk matrix.
#' @export
getCanSubspace = function(mat, vec) {
  mat_sub = mat%*%vec%*%solve(t(vec)%*%vec)%*%t(vec)

  return(mat_sub)
}

#' Get the estatistic for the permutations.
#'
#' @param cancor Numeric. Canonical Correlation estimate.
#' @param A An nxp matrix.
#' @param B An nxq matrix.
#' @param C An nxs matrix. Confounding variables.
#' @param type Character. Choice of statistic: Canonical correlation, Wilks'statistic or Roy's statistic.
#' @return Statistic
CCAtStat = function(cancor, A, B, C = 0, type = c("CC", "Wilks", "Roy")) {
  N = nrow(A)
  p = ncol(A)
  q = ncol(Y)
  K = length(cancor)


  # out-out sample cannonical correlation
  if(type == "CC") {
    tStat = cancor

    return(list(tStatistic = tStat))
  }

  # Wilks'statisitc
  if(type == "Wilks") {
    secondTerm = sapply(1:K, function(k) prod((1-cancor[k]^2)))
    tStat = sapply(1:K, function(k) - (N - C - (p + q + 3)/2)*secondTerm[k])
    df    = sapply(1:K, FUN =  function(k) (p - k + 1)*(q - k + 1))

    return(list(tStatistic = tStat, df = df))
  }

  # Roy's statistic
  if(type == "Roy") {
    tStat = cancor^2

    s = sapply(1:K, FUN =  function(k) min(p, q) - k + 1)
    m = sapply(1:K, FUN =  function(k) (abs(p - q) - 1)/2)
    n = sapply(1:K, FUN =  function(k) (N - C - p - q - 2)/2)

    param = list(s = s, m = m, n = n)

    return(list(tStatistic = tStat, distribParameters = param))

  }

}

#' Sparse Canonical Correlation Analysis. Computation of CC via NIPALS with soft thresholding.
#'
#' @param alphaInit Character. Type initialisation for \deqn{\mathbf{\alpha}}.
#' @param A,B Data matrices.
#' @param nonzero_a,nonzero_b Numeric. Scalar or vector over the number of nonzeroes allowed for a correlation estimate.
#' @param iter Numeric. Maximum number of iterations. Default is 20.
#' @param tol Numeric. Tolerance threshold. Default is 10^6.
#' @param silent Logical. If FALSE, a progress bar will appear on the console. Default is FALSE.
#' @return a list with the following elements:
#' \itemize{
#' \item{alpha}{Canonical vector for matrix \deqn{\mathbf{A}}, for each combination of sparsity value specified.}
#' \item{beta}{Canonical vector for matrix \deqn{\mathbf{B}}, for each combination of sparsity value specified.}
#' \item{cancor}{Max. canonical correlation estimate.}
#' \item{cancor_all}{Call canonical correlations calculated for each sparsity levels.}
SCCA = function(alphaInit, A, B, nonzero_a, nonzero_b, iter = 20, tol = 10^(-6), silent = FALSE)
{
  if(ncol(Y) <= max(nonzeroY)) {
    message("At least one of the nonzero options for Y is not sparse. Changing to meet criteria")
    nonzeroY = nonzeroY[nonzeroY < ncol(Y)]
  }

  if(ncol(X) <= max(nonzeroX)) {
    message("At least one of the nonzero options for X is not sparse. Changing to meet criteria")
    nonzeroX = nonzeroX[nonzeroX < ncol(X)]
  }

  #Create the matrix A
  alpha = sapply(nonzeroX, function(x) c(alphaInit))

  varTol1 = matrix(0, nrow = nrow(X), ncol = length(nonzeroX))
  varTol2 = matrix(0, nrow = nrow(Y), ncol = length(nonzeroX))
  i = 0
  e = 10
  while (e > tol & i <= iter) {
    i = i +1

    # refresh
    if(i > 1) varTol1 = gamma
    if(i > 1) varTol2 = zeta


    gamma =  X %*% alpha
    dist = sqrt(colSums(gamma^2))
    gamma = sweep(gamma, 2, dist, "/")

    beta = t(Y) %*% gamma

    beta = apply(rbind(beta,nonzeroY), 2, function(x)
    {
      nonzero1 = x[length(x)]
      y = x[-length(x)]
      thres = abs(y)[order(abs(y), decreasing=TRUE)[nonzero1+1]]
      tmp = (abs(y) - thres)
      tmp[tmp<=0] = 0
      sign(y) * tmp
    })

    zeta = Y %*% beta
    dist = sqrt(colSums(zeta^2))
    zeta = sweep(zeta, 2, dist, "/")

    alpha = t(X) %*% zeta

    alpha = apply(rbind(alpha,nonzeroX), 2, function(x)
    {
      nonzero1 = x[length(x)]
      y = x[-length(x)]
      thres = abs(y)[order(abs(y), decreasing=TRUE)[nonzero1+1]]
      tmp = (abs(y) - thres)
      tmp[tmp<=0] = 0
      sign(y) * tmp
    })

    if(length(nonzeroX) == 1) e = mean(abs(gamma - varTol1)) + mean(abs(zeta - varTol2))
    if(length(nonzeroX) > 1) e = mean(colMeans(abs(gamma - varTol1))) + mean(colMeans(abs(zeta - varTol2)))

    textSCCA = paste0(" Common convergence error: ", round(e, 5), " & Iterations: ", i)
    if(isFALSE(silent) & (e<= tol || i > iter)) cat(textSCCA, "\r")

  }

  alpha = sapply(1:ncol(alpha), function(p) standardVar(alpha[,p], centre = FALSE, normalise = TRUE))
  beta  = sapply(1:ncol(beta), function(q) standardVar(beta[,q], centre = FALSE, normalise = TRUE))

  #Return the CCA
  cancor_out = abs(sapply(1:ncol(gamma), function(j) cor(gamma[,j], zeta[,j])))


  out = list( alpha  = alpha,
              beta   = beta ,
              cancor = max(cancor_out),
              cancor_all = cancor_out)
}

#' Sparse Canonical Correlation Analysis. Computation of CC via NIPALS with soft thresholding.
#'
#' @param alphaInit Character. Type initialisation for \deqn{\mathbf{\alpha}}. Default is "eigen".
#' @param A,B Data matrices.
#' @param nonzero_a,nonzero_b Numeric. Scalar or vector over the number of nonzeroes allowed for a correlation estimate.
#' @param iter Numeric. Maximum number of iterations. Default is 20.
#' @param tol Numeric. Tolerance threshold. Default is 10^6.
#' @param silent Logical. If FALSE, a progress bar will appear on the console. Default is FALSE.
#' @return a list with the following elements:
#' \itemize{
#' \item{alpha}{Canonical vector for matrix \deqn{\mathbf{A}}, for each combination of sparsity value specified.}
#' \item{beta}{Canonical vector for matrix \deqn{\mathbf{B}}, for each combination of sparsity value specified.}
#' \item{cancor}{Max. canonical correlation estimate.}
#' \item{nonzero_a,nonzero_b}{Optimal nonzero values for each canonical vector.}
KFoldSCCA = function(A, B, nonzero_a, nonzero_b, alphaStart = "eigen", folds = 10, parallel_logic = FALSE, silent = FALSE, toPlot = TRUE, XTest_res = NULL, YTest_res = NULL) {
  N = nrow(A) # observations
  p = ncol(A) # predictor variables (not really since CCA is symmetric)
  q = ncol(B) # response variables (not really since CCA is symmetric)
  s = rep(1:folds, times=ceiling(N/folds))
  s = s[1:N]
  s = s[sample(1:length(s), length(s))]
  nonzeroGrid = expand.grid(nonzero_a, nonzero_b)
  h = nrow(nonzeroGrid)
  canCor = matrix(NA, folds, h)


  alphaMat <- list()
  betaMat  <- list()
  if(isTRUE(parallel_logic)) {
    # create env
    myCluster <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK")
    doParallel::registerDoParallel(cl = myCluster)
    if(isFALSE(foreach::getDoParRegistered())) stop("DoPar not registered. Check cores")

    suppressWarnings(
      canCor <- foreach(f = 1:folds, .combine = "rbind", .export = ls(globalenv())) %dopar% {
        ATrain = A[s!=f, ]
        BTrain = B[s!=f, ]
        if(!is.null(ATest_res)) ATest = ATest_res[s==f, ]
        if(is.null(ATest_res)) ATest  = A[s==f, ]
        if(!is.null(BTest_res)) BTest = BTest_res[s==f, ]
        if(is.null(BTest_res)) BTest  = B[s==f, ]


        # check between selected vector vs. one with higher cancor or follow up
        if(alphaStart == "eigen") alphaInit = initialiseCanVar(A = ATrain, B = BTrain)[,1]
        if(alphaStart == "random") alphaInit = standardVar(replicate(1, rnorm(p)), normalise = TRUE)
        if(alphaStart == "uniform") alphaInit = standardVar(matrix(rep(1, p), nrow = p, ncol = 1))


        if(isFALSE(silent)) progressBar(folds, f)

        resultKFold = SCCA(alphaInit = alphaInit, A = ATrain, B = BTrain, nonzero_a = nonzeroGrid[,1], nonzero_b = nonzeroGrid[,2], silent = silent)

        alphaMat[[f]] <- resultKFold$alpha
        betaMat[[f]]  <- resultKFold$beta

        gamma = standardVar(ATest%*%alphaMat[[f]], centre = FALSE)
        zeta  = standardVar(BTest%*%betaMat[[f]], centre = FALSE)

        if(ncol(gamma) + ncol(zeta) > 2) canCor[f,]  = abs(sapply(1:ncol(gamma), function(j) cor(gamma[,j], zeta[,j])))
        if(ncol(gamma) + ncol(zeta) == 2) canCor[f,] = abs(cor(gamma, zeta))

        # cat(canCor[1:k,], "\r")
        # cat(table(s), "\r")
        if(any(is.na(canCor[f,]))) stop("oneis NA")

        if(f == folds) resultKFold <<- resultKFold

        canCor[f,]

      } )
  } else {
    for (f in 1:folds) {
      ATrain = A[s!=f, ]
      BTrain = B[s!=f, ]
      if(!is.null(ATest_res)) ATest = ATest_res[s==f, ]
      if(is.null(ATest_res)) ATest  = A[s==f, ]
      if(!is.null(BTest_res)) BTest = BTest_res[s==f, ]
      if(is.null(BTest_res)) BTest  = B[s==f, ]


      # check between selected vector vs. one with higher cancor or follow up
      if(alphaStart == "eigen") alphaInit = initialiseCanVar(A = ATrain, B = BTrain)[,1]
      if(alphaStart == "random") alphaInit = standardVar(replicate(1, rnorm(p)), normalise = TRUE)
      if(alphaStart == "uniform") alphaInit = standardVar(matrix(rep(1, p), nrow = p, ncol = 1))


      if(isFALSE(silent)) progressBar(folds, f)

      resultKFold = SCCA(alphaInit = alphaInit, A = ATrain, B = BTrain, nonzero_a = nonzeroGrid[,1], nonzero_b = nonzeroGrid[,2], silent = silent)

      alphaMat[[f]] <- resultKFold$alpha
      betaMat[[f]]  <- resultKFold$beta

      gamma = standardVar(ATest%*%alphaMat[[f]], centre = FALSE)
      zeta  = standardVar(BTest%*%betaMat[[f]], centre = FALSE)

      if(ncol(gamma) + ncol(zeta) > 2) canCor[f,]  = abs(sapply(1:ncol(gamma), function(j) cor(gamma[,j], zeta[,j])))
      if(ncol(gamma) + ncol(zeta) == 2) canCor[f,] = abs(cor(gamma, zeta))

      # cat(canCor[1:k,], "\r")
      # cat(table(s), "\r")
      if(any(is.na(canCor[f,]))) stop("oneis NA")

    }

  }

  if(isTRUE(parallel_logic)) parallel::stopCluster(cl = myCluster)


  canCorKMeans = colSums(abs(canCor))/folds
  select       = getWhich(abs(canCorKMeans), max)
  canCorPrint  = canCorKMeans[select]

  names(canCorPrint) <- ("k-fold cv max. cancor")
  if(isFALSE(silent)) cat("\n")
  if(isFALSE(silent)) print(canCorPrint)
  if(isFALSE(silent)) {
    cat("\n ........................................ \n",
        paste0("# nonzero A: ", nonzeroGrid[select, 1],   "\n", "
               # nonzero B: ", nonzeroGrid[select, 2],
               "\n ........................................ \n"))
  }



  if(toPlot & isFALSE(parallel_logic)) {
    par(mfrow = c(2, 2))

    matplot(matrix(resultKFold$alpha, nrow = p), type = "l", ylab = "alpha", xlab = "p")
    title("Full set")

    matplot(resultKFold$alpha[,select], type = "l", ylab = "alpha", xlab = "p")
    title("K-fold CV best")

    matplot(matrix(resultKFold$beta, nrow = q), type = "l", ylab = "beta", xlab = "q")
    matplot(resultKFold$beta[,select], type = "l", ylab = "beta", xlab = "q")


  }

  if(isTRUE(parallel_logic))  alphaInit = initialiseCanVar(B, A)[,1]
  if(isFALSE(parallel_logic)) alphaInit = resultKFold$alpha[, select]

  result     = SCCA(alphaInit = alphaInit, A = A, B = B, nonzero_a = nonzeroGrid[select, 1], nonzero_b = nonzeroGrid[select, 2], silent = silent)


  resultSCCA = list(cancor = canCorPrint,
                    alpha  = result$alpha,
                    beta   = result$beta,
                    # alphaMat       = alphaMat,
                    # betaMat        = betaMat,
                    # position       = select,
                    nonzero_a = nonzeroGrid[select, 1],
                    nonzero_b = nonzeroGrid[select, 2])

  return(resultSCCA)
}

#' Sparse Canonical Correlation Analysis. Computation of CC via NIPALS with soft thresholding.
#'
#' @description This function performs CCA on matrices \deqn{\mathbf{A}} and \deqn{\mathbf{B}} via Non-Iterative PArtial Least Squares (NIPALS) algorithm
#' imposing sparsity over a fixed number of variables especified.
#' @param A,B Data matrices.
#' @param nonzero_a,nonzero_b Numeric. Scalar or vector over the number of nonzeroes allowed for a correlation estimate.
#' @param K Numeric. Number of components to be computed.
#' @param alphaStart Character. Type initialisation for \deqn{\mathbf{\alpha}}. Default is "eigen".
#' @param folds Numeric. Number of folds for the cross-validation process.
#' @param silent Logical. If FALSE, a progress bar will appear on the console. Default is FALSE.
#' @param toPlot Logical. If TRUE, plot will be generated automatically showing the estimated canonical weights. Default is TRUE.
#' @param typeResid Character. Choice of residualisation technique. Options are basic (default), null and LV.
#' @param combination Logical. If TRUE, the algorithm will search for the best combination of sparsity choice nonzero_a and nonzero_b for each component. This should be used for exploratory analysis. Default is FALSE.
#' @param parallel_logic Logical. If TRUE, cross-validation is done in parallel.Default is FALSE.
#' @details For a exploratory analysis nonzero_a and nonzero_b can be vectors. The algorithm will then search for the best combination of sparsity choice nonzero_a and nonzero_b for each component.
#' @return a list with the following elements:
#' \itemize{
#' \item{alpha}{Canonical vector for matrix \deqn{\mathbf{A}}, for each combination of sparsity value specified.}
#' \item{beta}{Canonical vector for matrix \deqn{\mathbf{B}}, for each combination of sparsity value specified.}
#' \item{cancor}{Max. canonical correlation estimate.}
#' @export
MSCCA = function(A, B, nonzero_a, nonzero_b, K = 1, alphaStart = "eigen", folds = 10, silent = FALSE, toPlot = TRUE, typeResid = "basic", combination = FALSE, parallel_logic = FALSE) {

  cancorComponents = matrix(NA, nrow = K, ncol = 1)
  alphaComponents  = matrix(NA, nrow = ncol(A), K)
  betaComponents   = matrix(NA, nrow = ncol(B), K)
  InitComponents   = matrix(NA, nrow = ncol(A), K)

  for(k in 1: K) {
    if(isFALSE(silent)) cat("\n__________________________________________ \n For component K = ", k, ": \n")
    if(k==1) {
      Ea = A
      Eb = B

    } else if(k > 1) {


      Ea = residualisation(vec = matrix(alphaComponents[, k - 1], nrow = ncol(Ea)),  mat = Ea, type = typeResid)
      Eb = residualisation(vec = matrix(betaComponents[, k - 1], nrow = ncol(Eb)),  mat = Eb, type = typeResid)

      Ea = standardVar(Ea)
      Eb = standardVar(Eb)
    }
    if(combination) {
      result = KFoldSCCA(A = Ea, B = Eb, nonzero_a, nonzero_b, alphaStart, folds, silent = silent, toPlot = toPlot, XTest_res = X, YTest_res = Y, parallel_logic = parallel_logic)

    } else {
      nonzero_aK = nonzero_a[k]
      nonzero_bK = nonzero_b[k]
      result = KFoldSCCA(A = Ea, B = Eb, nonzero_aK, nonzero_bK, alphaStart, folds, silent = silent, toPlot = toPlot, XTest_res = X, YTest_res = Y, parallel_logic = parallel_logic)

    }


    cancorComponents[k]  = result$cancor
    alphaComponents[, k] = result$alpha
    betaComponents[, k]  = result$beta


  }


  resultSCCA = list(cancor = cancorComponents,
                    alpha  = alphaComponents,
                    beta   = betaComponents
  )

  return(resultSCCA)
}

#' Permutation testing for MSCCA
#'
#' @description This function performs permutation testing on CC estimates.
#' @param A,B Data matrices.
#' @param nonzero_a,nonzero_b Numeric. Scalar or vector over the number of nonzeroes allowed for a correlation estimate.
#' @param K Numeric. Number of components to be computed.
#' @param draws Numeric. Number of permutations for each component.
#' @param folds Numeric. Number of folds for the cross-validation process.
#' @param silent Logical. If FALSE, a progress bar will appear on the console. Default is FALSE.
#' @param toPlot Logical. If TRUE, plot will be generated automatically showing the estimated canonical weights. Default is TRUE.
#' @param cancor Numeric. Scalar or vector: anonical correlation estimate(s).
#' @param parallel_logic Logical. If TRUE, cross-validation is done in parallel.Default is FALSE.
#' @param nuisanceVar Data with nuisance variables. For statistic type.
#' @param testStatType Character. Choice of statistic. Options are CC (default), Wilks and Roy.
#' @param combination Logical. If TRUE, the algorithm will search for the best combination of sparsity choice nonzero_a and nonzero_b for each component. This should be used for exploratory analysis. Default is FALSE.
#' @details For a exploratory analysis nonzero_a and nonzero_b can be vectors. The algorithm will then search for the best combination of sparsity choice nonzero_a and nonzero_b for each component.
#' @return Matrix with permutation estimates.
#' @export
permcvscca = function(A, B, nonzero_a, nonzero_b, K, folds = 10, toPlot = FALSE, draws = 20, cancor, bootCCA = NULL, silent = TRUE, parallel_logic = TRUE, nuisanceVar = 0, testStatType = "CC", combination = TRUE) {
  perm = matrix(NA, nrow = draws, ncol = K)

  if(isTRUE(parallel_logic)) {
    # create env
    myCluster <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK")
    doParallel::registerDoParallel(cl = myCluster)
    if(isFALSE(foreach::getDoParRegistered())) stop("DoPar not registered. Check cores")

    perm <- foreach(d = 1:draws, .combine = "rbind", .export = ls(globalenv())) %dopar% {
      if(isFALSE(silent)) progressBar(draws, d)

      ASample = A[sample(1:nrow(A), nrow(A)),]
      t(CCAtStat(MSCCA(A=ASample, B=B, K = K, combination = combination, nonzero_a=nonzero_a, nonzero_b=nonzero_b, toPlot = FALSE, silent = TRUE, parallel_logic = FALSE)$cancor, ASample, B, C = nuisanceVar, type = testStatType)[["tStatistic"]]) #off-sample cancor

    }

  } else {
    for (d in 1:draws) {
      # cat("|", rep(".", d), rep(" ", (draws-d)), "|", (d/draws)*100, "%\r")
      if(isFALSE(silent)) progressBar(draws, d)
      ASample = A[sample(1:nrow(A), nrow(A)),]
      perm[d,] = CCAtStat(MSCCA(A=ASample, B=B, K = K, combination = combination, nonzero_a=nonzero_a, nonzero_b=nonzero_b, toPlot = FALSE, silent = TRUE)$cancor, ASample, B, C = nuisanceVar, type = testStatType)[["tStatistic"]] #off-sample cancor


    }
  }


  if(isTRUE(parallel_logic)) parallel::stopCluster(cl = myCluster)


  margin = 0.05
  testStatistic = CCAtStat(cancor, A, B, C = nuisanceVar, type = testStatType)[["tStatistic"]]
  h =  hist(perm[, getWhich(testStatistic, max)], breaks = draws/2)

  if(is.null(bootCCA)) {
    permDensity = density(perm[, getWhich(testStatistic, max)])

    xlim = c(min(min(perm) - margin, testStatistic - margin), max(max(perm) + margin, testStatistic + margin))
    ylim = c(min((modes(permDensity)$y), min(perm[,getWhich(testStatistic, max)])), max((modes(permDensity)$y) + margin, max(h$counts) + margin))

    par(mfrow=c(1,1))
    hist(perm[, getWhich(testStatistic, max)], breaks = draws/2, xlab = "Canonical Correlation",
         xlim = xlim, ylim =  ylim, main = paste0("Distribution under de Null - ", testStatType, " Statistic") , col = "#93D9D9")
    lines(permDensity, col="black", lwd = 2)
    lines(seq(xlim[1], xlim[2], 0.01), dnorm(seq(xlim[1], xlim[2], 0.01), mean(perm[, getWhich(testStatistic, max)]), sd(perm[, getWhich(testStatistic, max)])), col="steelblue", lwd = 2, lty = 4)
    abline(v=testStatistic, col = "red", lwd = 2)
    legend("topleft", c("normal", "density", "model canCor"), col = c("steelblue", "black", "red"), lty=c(2, 1, 1), cex=0.8)
    text(x = as.character(testStatistic), y = 0.9*par('usr')[4], labels = as.character(1:K), cex = 0.9)

  } else {
    li = bootCCA$ccaInterval[1]
    ui = bootCCA$ccaInterval[2]

    xlim = c(min(min(perm) - margin, modelCanCor - margin, li - margin), max(max(perm) + margin, modelCanCor + margin, ui + margin))
    li = bootCCA$ccaInterval[1]
    ui = bootCCA$ccaInterval[2]

    xlim = c(min(min(perm) - margin, modelCanCor - margin), max(max(perm) + margin, modelCanCor + margin))
    permDensity = density(perm[, getWhich(testStatistic, max)])

    par(mfrow=c(1,1))
    hist(perm[, getWhich(testStatistic, max)], breaks = draws/2, xlab = "Canonical Correlation", xlim = xlim, ylim =  sort(modes(permDensity)$y), main = paste0("Distribution under de Null - ", testStatType, " Statistic") , col = "#93D9D9")
    lines(permDensity, col="black", lwd = 2)
    lines(seq(xlim[1], xlim[2], 0.01), dnorm(seq(xlim[1], xlim[2], 0.01), mean(perm[, getWhich(modelCanCor, max)]), sd(perm[, getWhich(modelCanCor, max)])), col="steelblue", lwd = 2, lty = 4)
    lines(x = c(li, li), y = c(max(histNullCCA$density)*0.5/2, max(histNullCCA$density)*0.6/2))
    lines(x = c(ui, ui), y = c(max(histNullCCA$density)*0.5/2, max(histNullCCA$density)*0.6/2))
    lines(x = c(li, ui), y = rep((max(histNullCCA$density)*0.5/2 + max(histNullCCA$density)*0.6/2)/2, 2))
    abline(v = modelCanCor, col = "purple", lwd = 2)
    abline(v=mean(bootCCA$boostCanCor), col = "red", lwd = 2)
    legend("topleft", c("normal", "density", "model canCor", "bootstrap avg canCor"), col = c("steelblue", "black", "purple", "red"), lty=c(2, 2, 1, 1), cex=0.8)
    text(x = as.character(testStatistic), y = 0.9*par('usr')[4], labels = as.character(1:K), cex = 0.9)
  }

  pValues = sapply(1:K, function (k) round(mean(testStatistic[k]<=(perm[, getWhich(testStatistic, max)])), 6))

  print(paste0("Empirical p-values: \n", pValues))

  return(permCanCor)

}


latentVariable <- function(data, alpha, beta) {
  X = data[[1]]
  Y = data[[2]]
  gamma  = matrix(X%*%alpha, nrow = nrow(X))
  zeta   = matrix(Y%*%beta, nrow = nrow(X))
  cancor = cor(gamma, zeta)

  return(cancor)
}

#' Bootstrap to get empirical intervals of canonical correlation.
#'
#' @param A,B Data matrices.
#' @param nonzero_a,nonzero_b Numeric. Scalar or vector over the number of nonzeroes allowed for a correlation estimate.
#' @param K Numeric. Number of components to be computed.
#' @param draws Numeric. Number of permutations for each component.
#' @param folds Numeric. Number of folds for the cross-validation process.
#' @param n Numeric. Times of bootstrapping.
#' @param ci_quant Numeric. Between 0 and 1. quantile of intervale.
#' @param silent Logical. If FALSE, a progress bar will appear on the console. Default is FALSE.
#' @param toPlot Logical. If TRUE, plot will be generated automatically showing the estimated canonical weights. Default is TRUE.
#' @param cancor Numeric. Scalar or vector: anonical correlation estimate(s).
#' @param parallel_logic Logical. If TRUE, cross-validation is done in parallel.Default is FALSE.
#' @param nuisanceVar Data with nuisance variables. For statistic type.
#' @param testStatType Character. Choice of statistic. Options are CC (default), Wilks and Roy.
#' @details For a exploratory analysis nonzero_a and nonzero_b can be vectors. The algorithm will then search for the best combination of sparsity choice nonzero_a and nonzero_b for each component.
#' @return Matrix with permutation estimates.
#' @export
boostrapCCA = function(A, B, nonzero_a, nonzero_b, cancor, folds = 10, n = 100, ci_quant = 0.01, silent = TRUE, toPlot = FALSE, parallel_logic = TRUE, nuisanceVar = 0, testStatType = "CC") {
  N = nrow(A)
  p = ncol(A)
  q = ncol(B)

  t = rep(NA, n)

  if(isTRUE(parallel_logic)) {
    # create env
    myCluster <- parallel::makeCluster(parallel::detectCores() - 1, type = "PSOCK")
    doParallel::registerDoParallel(cl = myCluster)
    if(isFALSE(foreach::getDoParRegistered())) stop("DoPar not registered. Check cores")

    t <- foreach(i = 1:n, .combine = "c", .export = ls(globalenv())) %dopar% {
      if(isFALSE(silent)) progressBar(n, i)

      s       = sample(1:N, N, replace = TRUE)
      ASample = A[s, ]
      BSample = B[s, ]
      data    = list(ASample, BSample)
      t[i] = CCAtStat(KFoldSCCA(A, B, alphaStart = "eigen", nonzero_a = nonzero_a, nonzero_b = nonzero_b, folds = folds, silent = TRUE, toPlot = FALSE)$cancor, ASample, B, C = nuisanceVar, type = testStatType)

    }

  } else {
    for(i in 1:n) {
      if(isFALSE(silent)) progressBar(n, i)
      s       = sample(1:N, N, replace = TRUE)
      ASample = A[s, ]
      BSample = B[s, ]
      data    = list(ASample, ASample)
      t[i] = CCAtStat(KFoldSCCA(A, B, alphaStart = "eigen", nonzero_a = nonzero_a, nonzero_b = nonzero_b, folds = folds, silent = TRUE, toPlot = FALSE)$cancor, XSample, Y, C = nuisanceVar, type = testStatType)


    }
  }

  parallel::stopCluster(cl = myCluster)


  d = t - cancor
  d = sort(d)
  d_quant = quantile(d, c(ci_quant, 1 - ci_quant))
  ci = cancor + d_quant

  if(toPlot) {
    par(mfrow =c(1, 1))
    plot(t, ylim = c(min(0.5, min(t)), 1), col = "grey")
    abline(h=cancor, col = "royalblue", lwd=2)
    abline(h=ci[1], col = "red", lwd=2, lty = 2)
    abline(h=ci[2], col = "red", lwd=2, lty = 2)
    abline(h=mean(t), col = "purple", lwd=2, lty = 3)
    legend("bottomright",  c("CI", "bootstrap canCor", "model canCor"), col = c("red", "purple", "royalblue"), lty=c(2, 3, 1), cex=0.8)

  }

  ciList = list(boostCanCor = t,
                ccaQuant = d_quant,
                ccaInterval = ci)
}

