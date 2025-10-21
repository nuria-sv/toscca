
# toscca ------------------------------------------------------------------
#' Performs eigen decomposition of a matrix in PS space.
#'
#'
#' @param A A square matrix nxn.
#' @return Matrix. Positive definite matrix.
fastEigen = function(A) {
  if(nrow(A) != ncol(A)) stop("Matrix should be square.")

  e = eigen(A)
  E = e$vector%*%diag(sqrt(pmax(0,e$values)))%*%t(e$vectors)   # make positive definite matrix
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
#' @examples
#' #sample size etc
#' N = 10
#' p = 25
#' q = 5
#' # noise
#' X0 = sapply(1:p, function(x) rnorm(N))
#' Y0 = sapply(1:q, function(x) rnorm(N))
#'
#' colnames(X0) = paste0("x", 1:p)
#' colnames(Y0) = paste0("y", 1:q)
#'
#' # signal
#' Z1 = rnorm(N,0,1)
#'
#'
#' #Some associations with the true signal
#' alpha = (6:10) / 10
#' beta  = -(2:3) / 10
#'
#' loc_alpha = 1:length(alpha)
#' loc_beta  = 1:length(beta)
#'
#' for(j in 1:length(alpha))
#'   X0[, loc_alpha[j]] =  alpha[j] * Z1 + rnorm(N,0,0.3)
#'
#' for(j in 1:length(beta))
#'   Y0[, loc_beta[j]] =  beta[j] * Z1 + rnorm(N,0,0.3)
#' cv = initialiseCanVar(X0, Y0)
#' @export
initialiseCanVar = function(A, B) {
  bbt = B%*%t(B) # positive definite matrix
  # A%*%t(A) is symmetric
  E   = fastEigen(bbt)
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
    # require(mcompanion)
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


#' Performs matrix residualisation over estimated canonical vectors by
#' using the null space of the estimated canonical vector to construct a
#' new matrix.
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
#' @export
toscca.tStat = function(cancor, A, B, C = 0, type = c("CC", "Wilks", "Roy")) {
  N = nrow(A)
  p = ncol(A)
  q = ncol(B)
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

toscca.core = function(alphaInit, A, B, nonzero_a, nonzero_b, iter = 20, tol = 10^(-6), silent = FALSE)
{

  # checks
  if(ncol(B) <= max(nonzero_b)) {
    message("At least one of the nonzero options for B is not sparse. Changing to meet criteria")
    nonzero_b = nonzero_b[nonzero_b < ncol(B)]
  }

  if(ncol(A) <= max(nonzero_a)) {
    message("At least one of the nonzero options for A is not sparse. Changing to meet criteria")
    nonzero_a = nonzero_a[nonzero_a < ncol(A)]
  }


  #Create the matrix A
  alpha = sapply(nonzero_a, function(x) c(alphaInit))

  varTol1 = matrix(0, nrow = nrow(A), ncol = length(nonzero_a))
  varTol2 = matrix(0, nrow = nrow(A), ncol = length(nonzero_b))
  i = 0
  e = 10
  while (e > tol & i <= iter) {
    i = i +1

    # refresh
    if(i > 1) varTol1 = gamma
    if(i > 1) varTol2 = zeta


    gamma =  A %*% alpha
    dist  = sqrt(colSums(gamma^2))
    gamma = sweep(gamma, 2, dist, "/")


    beta = t(B) %*% gamma

    beta = apply(rbind(beta,nonzero_b), 2, function(x)
    {
      nonzero1 = x[length(x)]
      y = x[-length(x)]
      thres = abs(y)[order(abs(y), decreasing=TRUE)[nonzero1+1]]
      tmp = (abs(y) - thres)
      tmp[tmp<=0] = 0
      sign(y) * tmp
    })

    zeta = B %*% beta
    dist = sqrt(colSums(zeta^2))
    zeta = sweep(zeta, 2, dist, "/")


    alpha = t(A) %*% zeta

    alpha = apply(rbind(alpha,nonzero_a), 2, function(x)
    {
      nonzero1 = x[length(x)]
      y = x[-length(x)]
      thres = abs(y)[order(abs(y), decreasing=TRUE)[nonzero1+1]]
      tmp = (abs(y) - thres)
      tmp[tmp<=0] = 0
      sign(y) * tmp
    })


    if(length(nonzero_a) == 1) e = mean(abs(gamma - varTol1)) + mean(abs(zeta - varTol2))
    if(length(nonzero_a) > 1) e  = mean(colMeans(abs(gamma - varTol1))) + mean(colMeans(abs(zeta - varTol2)))

    textSCCA = paste0(" Common convergence error: ", round(e, 5), " & Iterations: ", i)
    if(isFALSE(silent) & (e<= tol || i > iter)) cat(textSCCA, "\r")

  }


  #Return the CCA
  cancor_out = abs(sapply(1:ncol(gamma), function(j) cor(gamma[,j], zeta[,j])))


  out = list( alpha  = alpha,
              beta   = beta ,
              cancor = max(cancor_out),
              cancor_all = cancor_out)
}

#' Sparse Canonical Correlation Analysis. Computation of CC via NIPALS with soft thresholding.
#'
#' @param alpha_init Character. Type initialisation for \deqn{\mathbf{\alpha}}. Default is "eigen".
#' @param A,B Data matrices.
#' @param nonzero_a,nonzero_b Numeric. Scalar or vector over the number of nonzeroes allowed for a correlation estimate.
#' @param folds Integer. Indicates number of folds to perform.
#' @param parallel_logic Logical. TRUE to parallelise folds. Default is FALSE.
#' @param silent Logical. TRUE to keep silent output messages. Default is FALSE.
#' @param toPlot Logical. TRUE to plot results.
#' @param ATest_res NULL. Keep NULL.
#' @param BTest_res NULL. Keep NULL.
#' @importFrom foreach foreach
#' @importFrom stats rnorm
#' @importFrom graphics matplot
#' @importFrom graphics title
#' @return a list with the following elements:

toscca.folds = function(A, B, nonzero_a, nonzero_b, alpha_init, folds = 1, parallel_logic = FALSE, silent = FALSE, toPlot = TRUE, ATest_res = NULL, BTest_res = NULL) {
  N = nrow(A) # observations
  p = ncol(A) # predictor variables (not really since CCA is symmetric)
  q = ncol(B) # response variables (not really since CCA is symmetric)
  s = rep(1:folds, times=ceiling(N/folds))
  s = s[1:N]
  s = s[sample(1:length(s), length(s))]
  if(folds == 1) s[sample(1:N, 0.25*N)] = 2
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
        if(alpha_init == "eigen") alphaInit = initialiseCanVar(A = ATrain, B = BTrain)[,1]
        if(alpha_init == "random") alphaInit = standardVar(replicate(1, rnorm(p)), normalise = TRUE)
        if(alpha_init == "uniform") alphaInit = standardVar(replicate(1, runif(p)), normalise = TRUE)


        if(isFALSE(silent)) progressBar(folds, f)

        resultKFold = toscca.core(alphaInit = alphaInit, A = ATrain, B = BTrain, nonzero_a = nonzeroGrid[,1], nonzero_b = nonzeroGrid[,2], silent = silent)

        alphaMat[[f]] <- resultKFold$alpha
        betaMat[[f]]  <- resultKFold$beta

        gamma = standardVar(ATest%*%alphaMat[[f]], centre = FALSE)
        zeta  = standardVar(BTest%*%betaMat[[f]], centre = FALSE)

        if(ncol(gamma) + ncol(zeta) > 2) canCor[f,]  = abs(sapply(1:ncol(gamma), function(j) cor(gamma[,j], zeta[,j])))
        if(ncol(gamma) + ncol(zeta) == 2) canCor[f,] = abs(cor(gamma, zeta))


        if(any(is.na(canCor[f,]))) stop("oneis NA")

        if(f == folds) resultKFold <- resultKFold

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
      if(alpha_init == "eigen") alphaInit = initialiseCanVar(A = ATrain, B = BTrain)[,1]
      if(alpha_init == "random") alphaInit = standardVar(replicate(1, rnorm(p)), normalise = TRUE)
      if(alpha_init == "uniform") alphaInit = standardVar(replicate(1, runif(p)), normalise = TRUE)


      if(isFALSE(silent)) progressBar(folds, f)

      resultKFold = toscca.core(alphaInit = alphaInit, A = ATrain, B = BTrain, nonzero_a = nonzeroGrid[,1], nonzero_b = nonzeroGrid[,2], silent = silent)

      alphaMat[[f]] <- resultKFold$alpha
      betaMat[[f]]  <- resultKFold$beta

      gamma = standardVar(ATest%*%alphaMat[[f]], centre = FALSE)
      zeta  = standardVar(BTest%*%betaMat[[f]], centre = FALSE)

      if(ncol(gamma) + ncol(zeta) > 2) canCor[f,]  = abs(sapply(1:ncol(gamma), function(j) cor(gamma[,j], zeta[,j])))
      if(ncol(gamma) + ncol(zeta) == 2) canCor[f,] = abs(cor(gamma, zeta))


      if(any(is.na(canCor[f,]))) stop("oneis NA")

    }

  }

  if(isTRUE(parallel_logic)) parallel::stopCluster(cl = myCluster)


  canCorKMeans = colSums(abs(canCor))/folds
  select       = getWhich(abs(canCorKMeans), max)
  canCorPrint  = canCorKMeans[select]

  names(canCorPrint) <- ("k-fold cv max. cancor")
  if(isFALSE(silent)) cat("\n")
  if(isFALSE(silent)) message(canCorPrint)
  if(isFALSE(silent)) {
    cat("\n ........................................ \n",
        paste0("# nonzero A: ", nonzeroGrid[select, 1],   "\n",
               " # nonzero B: ", nonzeroGrid[select, 2],
               "\n ........................................ \n"))
  }
  mat = matrix(canCor, nrow = length(nonzero_a), ncol = length(nonzero_b))

  if(toPlot & nrow(nonzeroGrid) > 1) {
    mat = matrix(canCor, nrow = length(nonzero_a), ncol = length(nonzero_b))
    rownames(mat) = nonzero_a
    colnames(mat) = nonzero_b
    myHeatmap(mat)

  }

  if(toPlot & isFALSE(parallel_logic)) {
    oldpar <- par(no.readonly = TRUE) # code line i
    on.exit(par(oldpar))

    par(mfrow = c(2, 2))

    matplot(matrix(resultKFold$alpha, nrow = p), type = "l", ylab = "alpha", xlab = "p")
    title("Full set")

    matplot(resultKFold$alpha[,select], type = "l", ylab = "alpha", xlab = "p")
    title("K-fold CV best")

    matplot(matrix(resultKFold$beta, nrow = q), type = "l", ylab = "beta", xlab = "q")
    matplot(resultKFold$beta[,select], type = "l", ylab = "beta", xlab = "q")


  }

  if(isTRUE(parallel_logic))  alphaInit = initialiseCanVar(A, B)[,1]
  if(isFALSE(parallel_logic)) alphaInit = resultKFold$alpha[, select]

  result     = toscca.core(alphaInit = alphaInit, A = A, B = B, nonzero_a = nonzeroGrid[select, 1], nonzero_b = nonzeroGrid[select, 2], silent = TRUE)


  if(nrow(nonzeroGrid) > 1) {
    resultSCCA = list(cancor = canCorPrint,
                    alpha  = result$alpha,
                    beta   = result$beta,
                    # alphaMat       = alphaMat,
                    # betaMat        = betaMat,
                    # position       = select,
                    nonzero_a = nonzeroGrid[select, 1],
                    nonzero_b = nonzeroGrid[select, 2],
                    canCor_grid = mat)
  } else {
    resultSCCA = list(cancor = canCorPrint,
                      alpha  = result$alpha,
                      beta   = result$beta,
                      # alphaMat       = alphaMat,
                      # betaMat        = betaMat,
                      # position       = select,
                      nonzero_a = nonzeroGrid[select, 1],
                      nonzero_b = nonzeroGrid[select, 2])
    }

  return(resultSCCA)
}

#' Sparse Canonical Correlation Analysis. Computation of CC via NIPALS with soft thresholding.
#'
#' @param A,B Data matrices.
#' @param nonzero_a,nonzero_b Numeric. Scalar or vector over the number of nonzeroes allowed for a correlation estimate.
#' @param K Numeric. Number of components to be computed.
#' @param alpha_init Character. Type initialisation for \deqn{\mathbf{\alpha}}. Default is "eigen".
#' @param folds Numeric. Number of folds for the cross-validation process.
#' @param silent Logical. If FALSE, a progress bar will appear on the console. Default is FALSE.
#' @param toPlot Logical. If TRUE, plot will be generated automatically showing the estimated canonical weights. Default is TRUE.
#' @param typeResid Character. Choice of residualisation technique. Options are basic (default), null and LV.
#' @param combination Logical. If TRUE, the algorithm will search for the best combination of sparsity choice nonzero_a and nonzero_b for each component. This should be used for exploratory analysis. Default is FALSE.
#' @param parallel_logic Logical. If TRUE, cross-validation is done in parallel.Default is FALSE.
#' @return a list with the following elements:
#' @examples
#' #sample size etc
#' N = 10
#' p = 25
#' q = 5
#' # noise
#' X0 = sapply(1:p, function(x) rnorm(N))
#' Y0 = sapply(1:q, function(x) rnorm(N))
#'
#' colnames(X0) = paste0("x", 1:p)
#' colnames(Y0) = paste0("y", 1:q)
#'
#' # signal
#' Z1 = rnorm(N,0,1)
#'
#'
#' #Some associations with the true signal
#' alpha = (6:10) / 10
#' beta  = -(2:3) / 10
#'
#' loc_alpha = 1:length(alpha)
#' loc_beta  = 1:length(beta)
#'
#' for(j in 1:length(alpha))
#'   X0[, loc_alpha[j]] =  alpha[j] * Z1 + rnorm(N,0,0.3)
#'
#' for(j in 1:length(beta))
#'   Y0[, loc_beta[j]] =  beta[j] * Z1 + rnorm(N,0,0.3)
#'
# performa toscca
#' X = standardVar(X0)
#' Y = standardVar(Y0)
#' K = 2                                       # number of components to be estimated
#' nonz_x = c(2,5, 10, 20)                     # number of nonzero variables for X
#' nonz_y = c(1, 2, 3, 4)                      # number of nonzero variables for Y
#' init   = "uniform"                          # type of initialisation
#' cca_toscca  = toscca(X, Y, nonz_x, nonz_y, K, alpha_init = init, silent = TRUE, toPlot = FALSE)
#'
#' @return List with estimated toscca parameters.
#' @export
toscca = function(A, B, nonzero_a, nonzero_b, K = 1, alpha_init = c("eigen", "random", "uniform"), folds = 1, silent = FALSE, toPlot = TRUE, typeResid = "basic", combination = TRUE, parallel_logic = FALSE) {

  # checks


  if(K != length(nonzero_a) & isFALSE(combination)) {
    if(K > 1 & length(nonzero_a) == 1) {
      nonzero_a = rep(nonzero_a, K)
    } else {
      message("nonzero_a must have length 1 or K.")
    }
  }

  if(K != length(nonzero_b)& isFALSE(combination)) {
    if(K > 1 & length(nonzero_b) == 1) {
      nonzero_b = rep(nonzero_b, K)
    } else {
      message("nonzero_b must have length 1 or K.")
    }
  }


  # set up
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
      result = toscca.folds(A = Ea, B = Eb, nonzero_a, nonzero_b, alpha_init, folds, silent = silent, toPlot = toPlot, ATest_res = A, BTest_res = B, parallel_logic = parallel_logic)

    } else {
      nonzero_aK = nonzero_a[k]
      nonzero_bK = nonzero_b[k]
      result = toscca.folds(A = Ea, B = Eb, nonzero_aK, nonzero_bK, alpha_init, folds, silent = silent, toPlot = toPlot, ATest_res = A, BTest_res = B, parallel_logic = parallel_logic)

    }


    cancorComponents[k]  = result$cancor
    alphaComponents[, k] = result$alpha
    betaComponents[, k]  = result$beta


  }


  alpha = sapply(1:K, function(p) standardVar(alphaComponents[,p], centre = FALSE, normalise = TRUE))
  beta  = sapply(1:K, function(q) standardVar(betaComponents[,q], centre = FALSE, normalise = TRUE))


  resultSCCA = list(cancor = cancorComponents,
                    alpha  = alpha,
                    beta   = beta
  )

  return(resultSCCA)
}

#' Permutation testing for toscca
#'
#' @description This function performs permutation testing on CC estimates.
#' @param A,B Data matrices.
#' @param nonzero_a,nonzero_b Numeric. Scalar or vector over the number of nonzeroes allowed for a correlation estimate.
#' @param alpha_init Character. Type initialisation for \deqn{\mathbf{\alpha}}. Default is "eigen".
#' @param K Numeric. Number of components to be computed.
#' @param draws Numeric. Number of permutations for each component.
#' @param folds Numeric. Number of folds for the cross-validation process.
#' @param silent Logical. If FALSE, a progress bar will appear on the console. Default is FALSE.
#' @param toPlot Logical. If TRUE, plot will be generated automatically showing the estimated canonical weights. Default is TRUE.
#' @param cancor Numeric. Scalar or vector: anonical correlation estimate(s).
#' @param parallel_logic Logical. If TRUE, cross-validation is done in parallel.Default is FALSE.
#' @param ncores numeric. Number of cores to use in parallelisation. Default is detectCores() -1.
#' @param nuisanceVar Data with nuisance variables. For statistic type.
#' @param testStatType Character. Choice of statistic. Options are CC (default), Wilks and Roy.
#' @importFrom foreach foreach
#' @importFrom EnvStats  epdfPlot
#' @details For a exploratory analysis nonzero_a and nonzero_b can be vectors. The algorithm will then search for the best combination of sparsity choice nonzero_a and nonzero_b for each component.
#' @return Matrix with permutation estimates.
#' @examples
#' #sample size etc
#' N = 10
#' p = 25
#' q = 5
#' # noise
#' X0 = sapply(1:p, function(x) rnorm(N))
#' Y0 = sapply(1:q, function(x) rnorm(N))
#'
#' colnames(X0) = paste0("x", 1:p)
#' colnames(Y0) = paste0("y", 1:q)
#'
#' # signal
#' Z1 = rnorm(N,0,1)
#'
#'
#' #Some associations with the true signal
#' alpha = (6:10) / 10
#' beta  = -(2:3) / 10
#'
#' loc_alpha = 1:length(alpha)
#' loc_beta  = 1:length(beta)
#'
#' for(j in 1:length(alpha))
#'   X0[, loc_alpha[j]] =  alpha[j] * Z1 + rnorm(N,0,0.3)
#'
#' for(j in 1:length(beta))
#'   Y0[, loc_beta[j]] =  beta[j] * Z1 + rnorm(N,0,0.3)
#'
# performa toscca
#' X = standardVar(X0)
#' Y = standardVar(Y0)
#' K = 2                                       # number of components to be estimated
#' nonz_x = c(2,5, 10, 20)                     # number of nonzero variables for X
#' nonz_y = c(1, 2, 3, 4)                      # number of nonzero variables for Y
#' init   = "uniform"                          # type of initialisation
#' cca_toscca  = toscca(X, Y, nonz_x, nonz_y, K, alpha_init = init, silent = TRUE, toPlot = FALSE)
#' \donttest{
#' #dont run due to parallelisation.
#' cc = cca_toscca$cancor
#' perm_toscca = toscca.perm(X, Y, nonz_x, nonz_y, K = K, init, draws = 10, cancor = cc, ncores = 2)
#' }
#' @return List with permuted correlations and p-values.
#' @export
toscca.perm = function(A, B, nonzero_a, nonzero_b, K, alpha_init = c("eigen", "random", "uniform"), folds = 1, toPlot = FALSE, draws = 20, cancor, silent = TRUE, parallel_logic = TRUE, nuisanceVar = 0, testStatType = "CC", ncores = NULL) {

  perm = matrix(NA, nrow = draws, ncol = K)

  if(isTRUE(parallel_logic)) {
    # create env
    if(is.null(ncores)) ncores = parallel::detectCores() - 1
    myCluster <- parallel::makeCluster(ncores, type = "PSOCK")
    doParallel::registerDoParallel(cl = myCluster)
    if(isFALSE(foreach::getDoParRegistered())) stop("DoPar not registered. Check cores")

    perm <- foreach(d = 1:draws, .combine = "rbind", .export = ls(globalenv())) %dopar% {
      if(isFALSE(silent)) progressBar(draws, d)

      ASample = A[sample(1:nrow(A), nrow(A)),]
      t(toscca.tStat(toscca(A = ASample, B = B, K = K, alpha_init = alpha_init, combination = FALSE, nonzero_a=nonzero_a, nonzero_b=nonzero_b, toPlot = FALSE, silent = TRUE, parallel_logic = FALSE)$cancor, ASample, B, C = nuisanceVar, type = testStatType)[["tStatistic"]]) #off-sample cancor

    }

  } else {
    for (d in 1:draws) {
      # cat("|", rep(".", d), rep(" ", (draws-d)), "|", (d/draws)*100, "%\r")
      if(isFALSE(silent)) progressBar(draws, d)
      ASample = A[sample(1:nrow(A), nrow(A)),]
      perm[d,] = toscca.tStat(toscca(A = ASample, B = B, K = K, alpha_init = alpha_init, combination = FALSE, nonzero_a=nonzero_a, nonzero_b=nonzero_b, toPlot = FALSE, silent = TRUE)$cancor, ASample, B, C = nuisanceVar, type = testStatType)[["tStatistic"]] #off-sample cancor


    }
  }


  if(isTRUE(parallel_logic)) parallel::stopCluster(cl = myCluster)


  margin = 0.05
  testStatistic = toscca.tStat(cancor, A, B, C = nuisanceVar, type = testStatType)[["tStatistic"]]
  h =  hist(perm[, getWhich(testStatistic, max)], breaks = draws/2, plot = FALSE)

    permDensity = density(perm[, getWhich(testStatistic, max)])

    xlim = c(min(min(perm) - margin, testStatistic - margin), max(max(perm) + margin, testStatistic + margin))
    ylim = c(min((modes(permDensity)$y), min(perm[,getWhich(testStatistic, max)])), max((modes(permDensity)$y) + margin, max(h$counts) + margin))


    oldpar <- par(no.readonly = TRUE) # code line i
    on.exit(par(oldpar))

    par(mfrow=c(1,1))
    hist(perm[, getWhich(testStatistic, max)], breaks = draws/2, xlab = "Canonical Correlation",
         xlim = xlim, ylim =  ylim, main = paste0("Distribution under de Null - ", testStatType, " Statistic") , col =  scales::alpha("#01768c", 0.2))
    lines(permDensity, col="black", lwd = 2)
    epdfPlot(perm[,getWhich(testStatistic, max)], discrete = FALSE, density.arg.list = NULL, plot.it = TRUE,
             add = TRUE, epdf.col = "steelblue", epdf.lwd = 3 * par("cex"), epdf.lty = 1,
             curve.fill = FALSE, curve.fill.col = "steelblue", main = NULL, xlab = NULL, ylab = NULL)
    abline(v=testStatistic, col = "#e86d07", lwd = 2)
    legend("topleft", c("Empirical pdf", "density", "model canCor"), col = c("steelblue", "black", "red"), lty=c(2, 1, 1), cex=0.8)
    text(x = as.character(testStatistic), y = 0.9*par('usr')[4], labels = as.character(1:K), cex = 0.9)


  pValues = sapply(1:K, function (k) round(mean(testStatistic[k]<=(perm[, getWhich(testStatistic, max)])), 6))

  if(isFALSE(silent)) message(cat("Empirical p-values:", pValues, sep = "\n"))

  return(list(perm_estimate = perm, p_values = pValues))

}


toscca.lv <- function(data, alpha, beta) {
  X = data[[1]]
  Y = data[[2]]
  gamma  = matrix(X%*%alpha, nrow = nrow(X))
  zeta   = matrix(Y%*%beta, nrow = nrow(X))
  cancor = cor(gamma, zeta)

  return(cancor)
}




# tosccamm ----------------------------------------------------------------

# toscca-mm core
tosccamm.core = function(alphaInit, A, B, nonzero_a, nonzero_b, iter = 20, tol = 10^(-6), silent = FALSE, model = c("arima", "lme"), arformula = c(1,0,0), lmeformula = " ~ -1 + time + (1|id)")
{



  # checks
  if(ncol(B) <= max(nonzero_b)) {
    message("At least one of the nonzero options for B is not sparse. Changing to meet criteria")
    nonzeroB = nonzero_b[nonzero_b < (ncol(B) - 2)]
  }

  if(ncol(A) <= max(nonzero_a)) {
    message("At least one of the nonzero options for A is not sparse. Changing to meet criteria")
    nonzero_a = nonzero_a[nonzero_a < (ncol(A) - 2)]
  }


  #Create the matrix A
  alpha = sapply(nonzero_a, function(x) c(alphaInit))

  varTol1 = matrix(0, nrow = nrow(A), ncol = length(nonzero_a))
  varTol2 = matrix(0, nrow = nrow(B), ncol = length(nonzero_b))
  i = 0
  e = 10


  # format data
  id_a = A[,"id"]
  id_b = B[,"id"]
  time_a = A[,"time"]
  time_b = B[,"time"]

  A = as.matrix(A[, 3:ncol(A)])
  B = as.matrix(B[, 3:ncol(B)])

  while (e > tol & i <= iter) {
    i = i +1

    # refresh
    if(i > 1) varTol1 = gamma
    if(i > 1) varTol2 = zeta


    dist = sqrt(colSums(alpha^2))
    alpha = sweep(alpha, 2, dist, "/")

    gamma =  A %*% alpha
    # dist  = sqrt(colSums(gamma^2))
    # gamma = sweep(gamma, 2, dist, "/")
    gamma = as.matrix(scale_rm(data.frame(id = id_a, time = time_a, gamma = gamma))[,-c(1,2)])

    if(length(model) != 1 | !(model %in% c("arima", "lme"))) {
      model = "lme"
      lmeformula = " ~ -1 + time + (1|id)"

      warning("Model not correctly specify. \n Default is lme with formula ~ -1 + time + (1|id)")
    }

    if(model == "arima") {
      me = list()
      pred_me = matrix(NA, nrow(B), 1)

      if(is.null(arformula)){
        for(n in unique(id_a)){
          me[[n]] = forecast::auto.arima(gamma[n == id_a],max.p = 5,max.q = 5,max.P = 5,max.Q = 5,max.d = 3,seasonal = FALSE,ic = 'aicc')
          pred_me[which(id_b ==n)] = as.numeric(forecast::forecast(me[[n]], h = length(time_b[which(id_b == n)]))$fitted)

        }
      } else {
        for(n in unique(id_a)){
          me[[n]] = stats::arima(gamma[n == id_a], order = arformula, method = "ML")
          pred_me[which(id_b ==n)] = as.numeric(stats::predict(me[[n]], n.ahead = length(time_b[which(id_b == n)]))$pred)

        }
      }

      # for(n in unique(id_a)){
      #   pred_me[which(id_b ==n)] = as.numeric(predict(me[[n]], n.ahead = length(time_b[which(id_b == n)]))$pred)
      # }
    }

    if(model == "lme") {
      me = sapply(1:ncol(alpha), function(j) lme4::lmer(as.formula(paste("gamma", lmeformula)), data = data.frame(gamma = gamma[,j], time = time_a, id = id_a), REML = TRUE))
      pred_me = sapply(1:ncol(alpha), function(j) stats::predict(me[[j]], newdata = data.frame(time = time_b, id = id_b), allow.new.levels = TRUE, re.form = NULL))

    }

    beta = t(B) %*% pred_me #NAS IN PREDME

    rm(me, pred_me)

    beta = apply(rbind(beta,nonzero_b), 2, function(x)
    {
      nonzero1 = x[length(x)]
      y = x[-length(x)]
      thres = abs(y)[order(abs(y), decreasing=TRUE)[nonzero1+1]]
      tmp = (abs(y) - thres)
      tmp[tmp<=0] = 0
      sign(y) * tmp
    })

    dist = sqrt(colSums(beta^2))
    beta = sweep(beta, 2, dist, "/")

    zeta = B %*% beta
    # dist = sqrt(colSums(zeta^2))
    # zeta = sweep(zeta, 2, dist, "/")
    zeta = as.matrix(scale_rm(data.frame(id = id_b, time = time_b, zeta = zeta))[,-c(1,2)])


    if(model == "arima") {
      me = list()
      pred_me = matrix(NA, length(gamma), 1)

      if(is.null(arformula)){
        for(n in unique(id_b)){
          me[[n]] = forecast::auto.arima(zeta[n == id_b],max.p = 5,max.q = 5,max.P = 5,max.Q = 5,max.d = 3,seasonal = FALSE,ic = 'aicc')
          pred_me[which(id_a ==n)] = as.numeric(forecast::forecast(me[[n]], h = length(time_a[which(id_a == n)]))$fitted)

        }
      } else {
        for(n in unique(id_b)){
          me[[n]] = stats::arima(zeta[n == id_b], order = arformula, method = "ML")
          pred_me[which(id_a ==n)] = as.numeric(stats::predict(me[[n]], n.ahead = length(time_a[which(id_a == n)]))$pred)

        }
      }

      # for(n in unique(id_b)){
      #   pred_me[which(id_a ==n)] = as.numeric(predict(me[[n]], n.ahead = length(time_a[which(id_a == n)]))$pred)
      # }
    }

    if(model == "lme") {
      me = sapply(1:ncol(beta), function(j) lme4::lmer(as.formula(paste("zeta", lmeformula)), data = data.frame(zeta = zeta[,j], time = time_b, id = id_b), REML = TRUE))
      pred_me = sapply(1:ncol(beta), function(j) stats::predict(me[[j]], newdata = data.frame(time = time_a, id = id_a), allow.new.levels = TRUE, re.form = NULL))

    }
    alpha = t(A) %*% pred_me

    rm(me, pred_me)

    alpha = apply(rbind(alpha,nonzero_a), 2, function(x)
    {
      nonzero1 = x[length(x)]
      y = x[-length(x)]
      thres = abs(y)[order(abs(y), decreasing=TRUE)[nonzero1+1]]
      tmp = (abs(y) - thres)
      tmp[tmp<=0] = 0
      sign(y) * tmp
    })


    if(length(nonzero_a) == 1) e = mean(abs(gamma - varTol1)) + mean(abs(zeta - varTol2))
    if(length(nonzero_a) > 1) e  = mean(colMeans(abs(gamma - varTol1))) + mean(colMeans(abs(zeta - varTol2)))

    textSCCA = paste0(" Common convergence error: ", round(e, 5), " & Iterations: ", i)
    if(isFALSE(silent) & (e<= tol || i > iter)) cat(textSCCA, "\r")

  }

  if(model == "lme") {
    me_x = sapply(1:ncol(alpha), function(j) lme4::lmer(as.formula(paste("gamma", lmeformula)), data = data.frame(gamma = gamma[,j], time = time_a, id = id_a), REML = TRUE))
    # pred_x = sapply(1:ncol(alpha), function(j) predict(me_x[[j]], newdata = data.frame(time = time_b, id = id_b), allow.new.levels = TRUE, re.form = NULL))
    # lmer(as.formula(paste("gamma", lmeformula)), data = data.frame(gamma = gamma, time = time_a, id = id_a), REML = TRUE)
    me_y = sapply(1:ncol(beta), function(j) lme4::lmer(as.formula(paste("zeta", lmeformula)), data = data.frame(zeta = zeta[,j], time = time_b, id = id_b), REML = TRUE))
    # pred_y = sapply(1:ncol(beta), function(j) predict(me_y[[j]], newdata = data.frame(time = time_a, id = id_a), allow.new.levels = TRUE, re.form = NULL))
    # me_y = lmer(as.formula(paste("zeta", lmeformula)), data = data.frame(zeta = zeta, time = time_b, id = id_b), REML = TRUE)
  } else {
    me_x = NULL
    me_y = NULL
  }

  # alpha[alpha!=0] = scale(alpha[alpha!=0])
  # beta[beta!=0]   = scale(beta[beta!=0])
  return(list(a = alpha, b = beta, conv = e, iter = i, me_x = me_x, me_y = me_y))
}

# toscca-mm folds
#' Computes TOSCCA-MM
#'
#' This function estimates sparse canonical vectors for matrices with multiple measurements
#' and the trajectories of the latent variables.
#'
#' @param A A data.frame with id and time as first two columns.
#' @param B A data.frame with id and time as first two columns.
#' @param nonzero_a Integer. Threshold parameter for A.
#' @param nonzero_b Integer. Threshold parameter for B.
#' @param folds Integer. Indicates number of folds to perform.
#' @param parallel_logic Logical. TRUE to parallelise folds. Default is FALSE.
#' @param silent Logical. TRUE to keep silent output messages. Default is FALSE.
#' @param ATest_res NULL. Keep NULL.
#' @param BTest_res NULL. Keep NULL.
#' @param model Character. c("lme", "ar"). Model to fit longitudinal latent space.
#' @param lmeformula Character. LME formula. Default is " ~ -1 + time + (1|id)".
#' @param arformula Numeric vector. Choice of ARIMA. Default is c(1,0,0).
#' @importFrom MASS mvrnorm
#' @return Canonical vectors for k components.
#' @examples
#' \donttest{
#' # example code
#' #sample size etc
#' N = 10
#' p = 25
#' q = 5
#' X0 = list()
#' Y0 = list()
#'
#' #Some associations with the true signal
#' cwa = (6:10) / 10
#' cwb  = -(2:3) / 10
#'
#' alpha = rep(0, p)
#' beta = rep(0, q)
#'
#' loc_alpha = 1:length(alpha)
#' loc_beta  = 1:length(beta)
#'
#' alpha[loc_alpha] = cwa
#' beta[loc_beta] = cwb
#'
#' sg = matrix(c(1, 0.6, 0.3, rep(0, 2),
#'               0.6, 1, 0.6, 0.3, rep(0, 1),
#'               0.3, 0.6, 1, 0.6, 0.3,
#'               rep(0,1), 0.3, 0.6, 1, 0.6,
#'               rep(0,2), 0.3, 0.6, 1), ncol = 5)
#' for(i in 1:N)
#' {
#'   times = 1:5
#'   Zi1 = (sin(100*times))^times +   times * 0.65 +rnorm(1,0,0.95)
#'   Zi = cbind(Zi1)
#'   #Simulate data and add some noise
#'   X0i = sapply(1:p, function(a) MASS::mvrnorm(1, (Zi %*% t(alpha))[,a], Sigma = sg))
#'   Y0i = sapply(1:q, function(a) MASS::mvrnorm(1, (Zi %*% t(beta))[,a], Sigma = sg))
#'
#'   colnames(X0i) = paste0("X", 1:ncol(X0i))
#'   colnames(Y0i) = paste0("Y", 1:ncol(Y0i))
#'   #Check the simulated cross correlation
#'   #image(cor(X0i, Y0i))
#'
#'   #Remove some observations
#'   # p_observed = 1
#'   X0i = cbind(id=i, time=times, X0i)#[rbinom(length(times),1,p_observed)==1,]
#'   Y0i = cbind(id=i, time=times, Y0i)#[rbinom(length(times),1,p_observed)==1,]
#'
#'   X0[[i]] = X0i
#'   Y0[[i]] = Y0i
#' }
#'
#' X0 = do.call("rbind", X0)
#' Y0 = do.call("rbind", Y0)
#'
#' X = data.frame(X0); Y = data.frame(Y0)
#' nonz_a = c(2, 5, 10, 20)
#' nonz_b =  c(2, 3, 4)
#'
#' mod <- tosccamm(X, Y, folds = 2, nonzero_a = nonz_a, nonzero_b = nonz_b, silent = TRUE)
#' }
#' @return List with estimated tosccamm parameters.
#' @export
tosccamm = function(A, B, nonzero_a, nonzero_b, folds = 1, parallel_logic = FALSE, silent = FALSE, ATest_res = NULL, BTest_res = NULL, model = "lme", lmeformula = " ~ -1 + time + (1|id)", arformula = c(1,0,0)) {
  # N = min(nrow(A), nrow(B)) # observations
  # p = ncol(A) - 2 # predictor variables (not really since CCA is symmetric)
  # q = ncol(B) - 2# response variables (not really since CCA is symmetric)
  # s = rep(1:folds, length.out=length(unique(A$id)))
  # # s = s[1:N]
  # # s = s[sample(1:length(s), length(s))]
  # if(folds == 1) s[sample(1:N, 0.25*N)] = 2
  # nonzeroGrid = expand.grid(nonzero_a, nonzero_b)
  # h = nrow(nonzeroGrid)
  # canCor_a = matrix(NA, folds, h)
  # canCor_b = matrix(NA, folds, h)
  # canCor = matrix(NA, folds, h)
  #
  #
  # alphaMat <- list()
  # betaMat  <- list()
  #
  # for (f in 1:folds) {
  #   fold_ids <- unique(A$id)[s == f]
  #   ATrain = A[A$id %in% fold_ids, ]
  #   BTrain = B[B$id %in% fold_ids, ]
  #   if(!is.null(ATest_res)) ATest = ATest_res[!(A$id %in% fold_ids), ]
  #   if(is.null(ATest_res)) ATest  = A[!(A$id %in% fold_ids), ]
  #   if(!is.null(BTest_res)) BTest = BTest_res[!(B$id %in% fold_ids), ]
  #   if(is.null(BTest_res)) BTest  = B[!(B$id %in% fold_ids), ]
  #
  #   #
  #   #       # check between selected vector vs. one with higher cancor or follow up
  #   #       if(alpha_init == "eigen") alphaInit = initialiseCanVar(A = ATrain, B = BTrain)[,1]
  #   #       if(alpha_init == "random") alphaInit = standardVar(replicate(1, rnorm(p)), normalise = TRUE)
  #   #       if(alpha_init == "uniform") alphaInit = standardVar(replicate(1, runif(p)), normalise = TRUE)
  #
  #
  #   # if(isFALSE(silent)) progressBar(folds, f)
  #
  #   resultKFold = tosccamm.core(alphaInit = runif(ncol(A)-2), A = ATrain, B = BTrain, nonzero_a = nonzeroGrid[,1], nonzero_b = nonzeroGrid[,2], model = model, lmeformula = lmeformula, silent = silent)
  #
  #   alphaMat[[f]] <- resultKFold$a
  #   betaMat[[f]]  <- resultKFold$b
  #
  #   # gamma_a = sapply(1:ncol(alphaMat[[f]]), function(j) predict(resultKFold$me_x[[j]], newdata = data.frame(time = ATest[,"time"], id = ATest[,"id"]), allow.new.levels = TRUE, re.form = NULL))
  #   # gamma_a = data.frame(id = ATest$id, time = ATest$time, x = gamma_a)
  #   # gamma_a = data.frame(time = unique(sort(round(gamma_a$time,2))), sapply(1:nrow(nonzeroGrid), function(j) aggregate(gamma_a[,j+2], by = list(round(gamma_a$time,2)), FUN =mean)$x)) #scale_rm(gamma_a, centre = FALSE)
  #   # zeta_a = sapply(1:ncol(betaMat[[f]]), function(j) predict(resultKFold$me_y[[j]], newdata = data.frame(time = ATest[,"time"], id = ATest[,"id"]), allow.new.levels = TRUE, re.form = NULL))
  #   # zeta_a = data.frame(id = ATest$id, time = ATest$time, x = zeta_a)
  #   # zeta_a  = data.frame(time = unique(sort(round(zeta_a$time,2))), sapply(1:nrow(nonzeroGrid), function(j) aggregate(zeta_a[,j+2], by = list(round(zeta_a$time,2)), FUN =mean)$x)) #scale_rm(zeta_a, centre = TRUE)
  #   #
  #   # gamma_b = sapply(1:ncol(alphaMat[[f]]), function(j) predict(resultKFold$me_x[[j]], newdata = data.frame(time = BTest[,"time"], id = BTest[,"id"]), allow.new.levels = TRUE, re.form = NULL))
  #   # gamma_b = data.frame(id = BTest$id, time = BTest$time, x = gamma_b)
  #   # gamma_b = data.frame(time = unique(sort(round(gamma_b$time,2))), sapply(1:nrow(nonzeroGrid), function(j) aggregate(gamma_b[,j+2], by = list(round(gamma_b$time,2)), FUN =mean)$x)) #scale_rm(gamma_b, centre = FALSE)
  #   # zeta_b = sapply(1:ncol(betaMat[[f]]), function(j) predict(resultKFold$me_y[[j]], newdata = data.frame(time = BTest[,"time"], id = BTest[,"id"]), allow.new.levels = TRUE, re.form = NULL))
  #   # zeta_b = data.frame(id = BTest$id, time = BTest$time, x = zeta_b)
  #   # zeta_b  = data.frame(time = unique(sort(round(zeta_b$time,2))), sapply(1:nrow(nonzeroGrid), function(j) aggregate(zeta_b[,j+2], by = list(round(zeta_b$time,2)), FUN =mean)$x)) #scale_rm(zeta_b, centre = TRUE)
  #
  #   getlv = function(mat, long_mod, new_data, grid) {
  #     lv = sapply(1:ncol(mat), function(j) predict(long_mod[[j]], newdata = new_data, allow.new.levels = TRUE, re.form = NULL))
  #     lv =  data.frame(id = new_data$id, time = new_data$time, x = lv)
  #     lv =  data.frame(time = unique(sort(round(lv$time,2))), sapply(1:nrow(grid), function(j) aggregate(lv[,j+2], by = list(round(lv$time,2)), FUN =mean)$x))
  #
  #     return(lv)
  #
  #   }
  #
  #   gamma_a = getlv(alphaMat[[f]], resultKFold$me_x, data.frame(time = ATest[,"time"], id = ATest[,"id"]), nonzeroGrid)
  #   zeta_a  = getlv(betaMat[[f]], resultKFold$me_y, data.frame(time = ATest[,"time"], id = ATest[,"id"]), nonzeroGrid)
  #   gamma_b = getlv(alphaMat[[f]], resultKFold$me_x, data.frame(time = BTest[,"time"], id = BTest[,"id"]), nonzeroGrid)
  #   zeta_b  = getlv(betaMat[[f]], resultKFold$me_y, data.frame(time = BTest[,"time"], id = BTest[,"id"]), nonzeroGrid)
  #   # gammaT <- aggregate(. ~ id + time, data = gamma, FUN = mean)
  #   # zettaT <- aggregate(. ~ id + time, data = zeta, FUN = mean)
  #
  #
  #   # cc_time = matrix(NA, ncol = ncol(gamma_a)-1, nrow =length(unique(A$time)))
  #   #   for (m in unique(A$time)) {
  #   #     w = unique(intersect(A[A$time==m,]$id, B[B$time==m,]$id))
  #   #     g = as.matrix(A[w, -c(1,2)])%*%resultKFold$a; e = as.matrix(B[w, -c(1,2)])%*%resultKFold$b;
  #   #     cc_time[m, ] = sapply(1:(ncol(gamma_a) -1), function(x) cor(g[,x], e[,x]))
  #   #   }
  #   cc = 1# colMeans(cc_time)
  #
  #   if(ncol(gamma_a) + ncol(zeta_a) > 2) canCor_a[f,]  = abs(sapply(2:ncol(gamma_a), function(j) cor(gamma_a[,j], zeta_a[,j])*cc[j-1]))
  #   if(ncol(gamma_a) + ncol(zeta_a) == 2) canCor_a[f,] = abs(cor(gamma_a, zeta_a_b)*cc)
  #
  #   if(ncol(gamma_b) + ncol(zeta_b) > 2) canCor_b[f,]  = abs(sapply(2:ncol(gamma_b), function(j) cor(gamma_b[,j], zeta_b[,j])*cc[j-1]))
  #   if(ncol(gamma_b) + ncol(zeta_b) == 2) canCor_b[f,] = abs(cor(gamma_b, zeta_b)*cc)
  #
  #   canCor[f,] = sapply(1:length(canCor_a[f,]), function(x) mean(canCor_a[f,x], canCor_b[f,x]))
  #
  #   if(any(is.na(canCor[f,]))) stop("oneis NA")
  #
  # }
  #
  # canCorKMeans = colSums(abs(canCor))/folds
  # select       = getWhich(abs(canCorKMeans), max)
  # if(length(select)>1) select <- select[1]
  # canCorPrint  = canCorKMeans[select]
  #
  # cat("📊 K-Fold CV Max. Canonical Correlation:  ", canCorPrint, "                      ")
  #
  # if (isFALSE(silent)) cat("\n")
  # # if (isFALSE(silent)) ca(canCorPrint)
  #
  # if (isFALSE(silent)) {
  #   a <- nonzeroGrid[select, 1]
  #   b <- nonzeroGrid[select, 2]
  #
  #   cat("\n",
  #       "   ──────────────────────────\n",
  #       sprintf("      Nonzero A: %8d \n", a),
  #       sprintf("      Nonzero B: %8d \n", b),
  #       "   ──────────────────────────\n",
  #       sep = "")
  # }
  #
  # if(nrow(nonzeroGrid) > 1) {
  #   mat = matrix(canCorKMeans, nrow = length(nonzero_a), byrow = F)
  #   # mat = matrix(canCorKMeans, nrow = length(nonzero_a), ncol = length(nonzero_b))
  #   rownames(mat) = nonzero_a
  #   colnames(mat) = nonzero_b
  #   # myHeatmap(mat)
  #
  # }
  #
  # result     = tosccamm.core(alphaInit =  runif(ncol(A)-2), A = A, B = B, nonzero_a = nonzeroGrid[select, 1], nonzero_b = nonzeroGrid[select, 2], silent = TRUE,  model = model, arformula = arformula, lmeformula = lmeformula)
  #
  # cat("done")
  # if(nrow(nonzeroGrid) > 1) {
  #   resultSCCA = list(cancor = canCorPrint,
  #                     alpha  = result$a,
  #                     beta   = result$b,
  #                     me_x   = result$me_x,
  #                     me_y   = result$me_y,
  #                     cv_logical = nrow(nonzeroGrid) > 1,
  #                     mat_cc = mat,
  #                     # alphaMat       = alphaMat,
  #                     # betaMat        = betaMat,
  #                     # position       = select,
  #                     nonzero_a = nonzeroGrid[select, 1],
  #                     nonzero_b = nonzeroGrid[select, 2],
  #                     canCor_grid = mat)
  # } else {
  #   resultSCCA = list(cancor = canCorPrint,
  #                     alpha  = result$a,
  #                     beta   = result$b,
  #                     me_x   = result$me_x,
  #                     me_y   = result$me_y,
  #                     cv_logical = nrow(nonzeroGrid) > 1,
  #                     # alphaMat       = alphaMat,
  #                     # betaMat        = betaMat,
  #                     # position       = select,
  #                     nonzero_a = nonzeroGrid[select, 1],
  #                     nonzero_b = nonzeroGrid[select, 2])
  # }

  N = min(nrow(A), nrow(B)) # observations
  p = ncol(A) - 2 # predictor variables (not really since CCA is symmetric)
  q = ncol(B) - 2# response variables (not really since CCA is symmetric)
  s = rep(1:folds, length.out=length(unique(A$id)))
  # s = s[1:N]
  # s = s[sample(1:length(s), length(s))]
  if(folds == 1) s[sample(1:N, 0.25*N)] = 2
  nonzeroGrid = expand.grid(nonzero_a, nonzero_b)
  h = nrow(nonzeroGrid)
  canCor_a = matrix(NA, folds, h)
  canCor_b = matrix(NA, folds, h)
  canCor = matrix(NA, folds, h)


  alphaMat <- list()
  betaMat  <- list()

  for (f in 1:folds) {
    fold_ids <- unique(A$id)[s == f]
    ATrain = A[A$id %in% fold_ids, ]
    BTrain = B[B$id %in% fold_ids, ]
    if(!is.null(ATest_res)) ATest = ATest_res[!(A$id %in% fold_ids), ]
    if(is.null(ATest_res)) ATest  = A[!(A$id %in% fold_ids), ]
    if(!is.null(BTest_res)) BTest = BTest_res[!(B$id %in% fold_ids), ]
    if(is.null(BTest_res)) BTest  = B[!(B$id %in% fold_ids), ]

    #
    #       # check between selected vector vs. one with higher cancor or follow up
    #       if(alpha_init == "eigen") alphaInit = initialiseCanVar(A = ATrain, B = BTrain)[,1]
    #       if(alpha_init == "random") alphaInit = standardVar(replicate(1, rnorm(p)), normalise = TRUE)
    #       if(alpha_init == "uniform") alphaInit = standardVar(replicate(1, runif(p)), normalise = TRUE)


    # if(isFALSE(silent)) progressBar(folds, f)

    resultKFold = tosccamm.core(alphaInit = runif(ncol(A)-2), A = ATrain, B = BTrain, nonzero_a = nonzeroGrid[,1], nonzero_b = nonzeroGrid[,2], model = model, lmeformula = lmeformula, silent = silent)

    alphaMat[[f]] <- resultKFold$a
    betaMat[[f]]  <- resultKFold$b

    gamma_a = sapply(1:ncol(alphaMat[[f]]), function(j) predict(resultKFold$me_x[[j]], newdata = data.frame(time = ATest[,"time"], id = ATest[,"id"]), allow.new.levels = TRUE, re.form = NULL))
    gamma_a = data.frame(id = ATest$id, time = ATest$time, x = gamma_a)
    gamma_a = data.frame(time = unique(sort(round(gamma_a$time,2))), sapply(1:nrow(nonzeroGrid), function(j) aggregate(gamma_a[,j+2], by = list(round(gamma_a$time,2)), FUN =mean)$x)) #scale_rm(gamma_a, centre = FALSE)
    zeta_a = sapply(1:ncol(betaMat[[f]]), function(j) predict(resultKFold$me_y[[j]], newdata = data.frame(time = ATest[,"time"], id = ATest[,"id"]), allow.new.levels = TRUE, re.form = NULL))
    zeta_a = data.frame(id = ATest$id, time = ATest$time, x = zeta_a)
    zeta_a  = data.frame(time = unique(sort(round(zeta_a$time,2))), sapply(1:nrow(nonzeroGrid), function(j) aggregate(zeta_a[,j+2], by = list(round(zeta_a$time,2)), FUN =mean)$x)) #scale_rm(zeta_a, centre = TRUE)

    gamma_b = sapply(1:ncol(alphaMat[[f]]), function(j) predict(resultKFold$me_x[[j]], newdata = data.frame(time = BTest[,"time"], id = BTest[,"id"]), allow.new.levels = TRUE, re.form = NULL))
    gamma_b = data.frame(id = BTest$id, time = BTest$time, x = gamma_b)
    gamma_b = data.frame(time = unique(sort(round(gamma_b$time,2))), sapply(1:nrow(nonzeroGrid), function(j) aggregate(gamma_b[,j+2], by = list(round(gamma_b$time,2)), FUN =mean)$x)) #scale_rm(gamma_b, centre = FALSE)
    zeta_b = sapply(1:ncol(betaMat[[f]]), function(j) predict(resultKFold$me_y[[j]], newdata = data.frame(time = BTest[,"time"], id = BTest[,"id"]), allow.new.levels = TRUE, re.form = NULL))
    zeta_b = data.frame(id = BTest$id, time = BTest$time, x = zeta_b)
    zeta_b  = data.frame(time = unique(sort(round(zeta_b$time,2))), sapply(1:nrow(nonzeroGrid), function(j) aggregate(zeta_b[,j+2], by = list(round(zeta_b$time,2)), FUN =mean)$x)) #scale_rm(zeta_b, centre = TRUE)

    # gammaT <- aggregate(. ~ id + time, data = gamma, FUN = mean)
    # zettaT <- aggregate(. ~ id + time, data = zeta, FUN = mean)

    tv = sort(union(unique(ATest$time), unique(BTest$time)))
    l = sapply(1:length(tv), function(v) length(intersect(ATest[ATest$time==v, ]$id, BTest[BTest$time==v, ]$id)))
    v = which(l==max(l))[1]
    iv = intersect(ATest[ATest$time==v, ]$id, BTest[BTest$time==v, ]$id)
    g = data.frame(id = ATest$id, time = ATest$time, as.matrix(ATest[,-c(1,2)])%*%alphaMat[[f]]); e = data.frame(id = BTest$id, time = BTest$time, as.matrix(BTest[,-c(1,2)])%*%betaMat[[f]])
    m = sapply(1:ncol(alphaMat[[f]]), function(x) abs(cor(g[g$id %in% iv & g$time==v, x+2], e[e$id %in% iv & e$time==v, x+2])))

    if(ncol(gamma_a) + ncol(zeta_a) > 2) canCor_a[f,]  = abs(sapply(2:ncol(gamma_a), function(j) cor(gamma_a[,j], zeta_a[,j])))
    if(ncol(gamma_a) + ncol(zeta_a) == 2) canCor_a[f,] = abs(cor(gamma_a, zeta_a))

    if(ncol(gamma_b) + ncol(zeta_b) > 2) canCor_b[f,]  = abs(sapply(2:ncol(gamma_b), function(j) cor(gamma_b[,j], zeta_b[,j])))
    if(ncol(gamma_b) + ncol(zeta_b) == 2) canCor_b[f,] = abs(cor(gamma_b, zeta_b))

    canCor[f,] = sapply(1:length(canCor_a[f,]), function(x) mean(canCor_a[f,x], canCor_b[f,x])*m[x]) #

    if(any(is.na(canCor[f,]))) stop("oneis NA")

  }

  canCorKMeans = colSums(abs(canCor))/folds
  select       = getWhich(abs(canCorKMeans), max)
  canCorPrint  = canCorKMeans[select]

  names(canCorPrint) <- ("k-fold cv max. cancor")
  if(isFALSE(silent)) cat("\n")
  if(isFALSE(silent)) message(canCorPrint)
  if(isFALSE(silent)) {
    cat("\n ........................................ \n",
        paste0("# nonzero A: ", nonzeroGrid[select, 1],   "\n",
               " # nonzero B: ", nonzeroGrid[select, 2],
               "\n ........................................ \n"))
  }

  if(nrow(nonzeroGrid) > 1) {
    mat = matrix(canCorKMeans, nrow = length(nonzero_a), byrow = FALSE)
    # mat = matrix(canCorKMeans, nrow = length(nonzero_a), ncol = length(nonzero_b))
    rownames(mat) = nonzero_a
    colnames(mat) = nonzero_b
    myHeatmap(mat)

  }

  result     = tosccamm.core(alphaInit =  runif(ncol(A)-2), A = A, B = B, nonzero_a = nonzeroGrid[select, 1], nonzero_b = nonzeroGrid[select, 2], silent = TRUE,  model = model, arformula = arformula, lmeformula = lmeformula)


  if(nrow(nonzeroGrid) > 1) {
    resultSCCA = list(cancor = canCorPrint,
                      alpha  = result$a,
                      beta   = result$b,
                      me_x   = result$me_x,
                      me_y   = result$me_y,
                      # alphaMat       = alphaMat,
                      # betaMat        = betaMat,
                      # position       = select,
                      nonzero_a = nonzeroGrid[select, 1],
                      nonzero_b = nonzeroGrid[select, 2],
                      canCor_grid = mat)
  } else {
    resultSCCA = list(cancor = canCorPrint,
                      alpha  = result$a,
                      beta   = result$b,
                      me_x   = result$me_x,
                      me_y   = result$me_y,
                      # alphaMat       = alphaMat,
                      # betaMat        = betaMat,
                      # position       = select,
                      nonzero_a = nonzeroGrid[select, 1],
                      nonzero_b = nonzeroGrid[select, 2])
  }

  return(resultSCCA)
}



# permutation testing tosccamm
#' Computes permutatied cc fot TOSCCA-MM
#'
#' This function estimates sparse canonical vectors for permuted matrices with multiple measurements.
#' @param A A matrix.
#' @param B A matrix.
#' @param nonzero_a Integer. Threshold parameter for A.
#' @param nonzero_b Integer. Threshold parameter for B.
#' @param K Interger. Numner of components.
#' @param toPlot Logical. Indicates if plots will be produced. Default is False.
#' @param draws Integer. Number of draws in permutation.
#' @param cancor Numeric vector of length K with estimated canonical correlations.
#' @param folds Integer. Indicates number of folds to perform.
#' @param parallel_logic Logical. TRUE to parallelise folds. Default is FALSE.
#' @param ncores numeric. Number of cores to use in parallelisation. Default is detectCores() -1.
#' @param nuisanceVar Numeric. Number of nuisance variables.
#' @param bootCCA deprecated.
#' @param testStatType Character. Choice of test-statistic c("CC", "Wilks", "Roy"),
#' @param silent Logical. TRUE to keep silent output messages. Default is FALSE.
#' @param model Character. c("lme", "ar"). Model to fit longitudinal latent space.
#' @param lmeformula Character. LME formula. Default is " ~ -1 + time + (1|id)".
#' @param arformula Numeric vector. Choice of ARIMA. Default is c(1,0,0).
#' @return Permuted canonical correlation for ell K and p-values.
#'
#' @importFrom graphics abline hist legend lines par text
#' @importFrom utils install.packages installed.packages
#' @importFrom stats aggregate arima as.formula coefficients cor density predict runif
#' @importFrom foreach %dopar%
#' @examples
#' # example code
#' \donttest{
#' # dont run due to parallel processes
#' #sample size etc
#' N = 10
#' p = 25
#' q = 5
#' X0 = list()
#' Y0 = list()
#'
#' #Some associations with the true signal
#' cwa = (6:10) / 10
#' cwb  = -(2:3) / 10
#'
#' alpha = rep(0, p)
#' beta = rep(0, q)
#'
#' loc_alpha = 1:length(alpha)
#' loc_beta  = 1:length(beta)
#'
#' alpha[loc_alpha] = cwa
#' beta[loc_beta] = cwb
#'
#' sg = matrix(c(1, 0.6, 0.3, rep(0, 2),
#'               0.6, 1, 0.6, 0.3, rep(0, 1),
#'               0.3, 0.6, 1, 0.6, 0.3,
#'               rep(0,1), 0.3, 0.6, 1, 0.6,
#'               rep(0,2), 0.3, 0.6, 1), ncol = 5)
#' for(i in 1:N)
#' {
#'   times = 1:5
#'   Zi1 = (sin(100*times))^times +   times * 0.65 +rnorm(1,0,0.95)
#'   Zi = cbind(Zi1)
#'   #Simulate data and add some noise
#'   X0i = sapply(1:p, function(a) MASS::mvrnorm(1, (Zi %*% t(alpha))[,a], Sigma = sg))
#'   Y0i = sapply(1:q, function(a) MASS::mvrnorm(1, (Zi %*% t(beta))[,a], Sigma = sg))
#'
#'   colnames(X0i) = paste0("X", 1:ncol(X0i))
#'   colnames(Y0i) = paste0("Y", 1:ncol(Y0i))
#'   #Check the simulated cross correlation
#'   #image(cor(X0i, Y0i))
#'
#'   #Remove some observations
#'   # p_observed = 1
#'   X0i = cbind(id=i, time=times, X0i)#[rbinom(length(times),1,p_observed)==1,]
#'   Y0i = cbind(id=i, time=times, Y0i)#[rbinom(length(times),1,p_observed)==1,]
#'
#'   X0[[i]] = X0i
#'   Y0[[i]] = Y0i
#' }
#'
#' X0 = do.call("rbind", X0)
#' Y0 = do.call("rbind", Y0)
#'
#' X = data.frame(X0); Y = data.frame(Y0)
#' nonz_a = c(2, 5, 10, 20)
#' nonz_b =  c(2, 3, 4)
#'
#' mod <- tosccamm(X, Y, folds = 2, nonzero_a = nonz_a, nonzero_b = nonz_b, silent = TRUE)
#' nza <- mod$nonzero_a
#' nzb <- mod$nonzero_b
#' cc  <- mod$cancor
#' perm_cc <- toscamm.perm(X,Y, nonzero_a=nza, nonzero_b=nzb,cancor=cc, ncores=2, draws = 10)
#' }
#' @return List with permuted correlations and p-values.
#' @export

toscamm.perm = function (A, B, nonzero_a, nonzero_b, K=1, folds = 1, toPlot = FALSE, draws = 1000,
                         cancor, bootCCA = NULL, silent = TRUE, parallel_logic = TRUE,
                         nuisanceVar = 0, testStatType = "CC", model = "lme", lmeformula = " ~ 0 + poly(time,3) + (1|id)", arformula = NULL, ncores = NULL)
{
  histNullCCA <- modelCanCor <- x <- y <- `..count..` <- NULL

  perm = matrix(NA, nrow = draws, ncol = K)
  if (isTRUE(parallel_logic)) {
    if(is.null(ncores)) ncores = parallel::detectCores() - 1
    myCluster <- parallel::makeCluster(ncores, type = "PSOCK")
    doParallel::registerDoParallel(cl = myCluster)
    if (isFALSE(foreach::getDoParRegistered()))
      stop("DoPar not registered. Check cores")
    perm <- foreach::foreach(d = 1:draws, .combine = "rbind", .export = ls(globalenv())) %dopar%
      {
        res_perm = list()
        ASample = data.frame(A[,1:2], A[sample(1:nrow(A), nrow(A)), -c(1,2)])
        X.temp = ASample
        Y.temp = B
        cc_time = matrix(NA, length(unique(Y.temp$time)), K)
        for (k in 1:K) {
          if(k > 1) {
            # residualise for subsequent components
            X.temp = data.frame(X.temp[,c(1,2)],toscca::residualisation(as.matrix(X.temp[,-c(1,2)]), res_perm[[k-1]]$alpha, type = "basic") )
            Y.temp = data.frame(Y.temp[,c(1,2)],toscca::residualisation(as.matrix(Y.temp[,-c(1,2)]), res_perm[[k-1]]$beta, type = "basic") )

            nz_a_gen = as.numeric(table(res_perm[[k-1]]$alpha != 0)[2])
            nz_b_gen = as.numeric(table(res_perm[[k-1]]$beta != 0)[2])
          }

          res_perm[[k]] <- tosccamm(X.temp, Y.temp, folds = 2,
                                    nonzero_a[k], nonzero_b[k],
                                    model = model, lmeformula = lmeformula)

          #  for (s in unique(X.temp$time)) {
          #    w = unique(intersect(X.temp[X.temp$time==s,]$id, Y.temp[Y.temp$time==s,]$id))
          #    g = as.matrix(X.temp[w, -c(1,2)])%*%res_perm[[k]]$alpha; e = as.matrix(Y.temp[w, -c(1,2)])%*%res_perm[[k]]$beta;
          #    cc_time[s, k] = abs(cor(g, e))
          #  }
          #
          # m= colMeans(cc_time)

        }
        sapply(1:K, function(k) res_perm[[k]]$cancor)

      }
  }
  else {
    for (d in 1:draws) {
      # if (isFALSE(silent))
      # progressBar(draws, d)

      res_perm = list()
      ASample = data.frame(A[,1:2], A[sample(1:nrow(A), nrow(A)), -c(1,2)])
      X.temp = ASample
      Y.temp = B
      cc_time = matrix(NA, length(unique(Y.temp$time)), K)
      for (k in 1:K) {
        if(k > 1) {
          # residualise for subsequent components
          X.temp = data.frame(X.temp[,c(1,2)],toscca::residualisation(as.matrix(X.temp[,-c(1,2)]), res_perm[[k-1]]$alpha, type = "basic") )
          Y.temp = data.frame(Y.temp[,c(1,2)],toscca::residualisation(as.matrix(Y.temp[,-c(1,2)]), res_perm[[k-1]]$beta, type = "basic") )

          nz_a_gen = as.numeric(table(res_perm[[k-1]]$alpha != 0)[2])
          nz_b_gen = as.numeric(table(res_perm[[k-1]]$beta != 0)[2])
        }

        res_perm[[k]] <- tosccamm(X.temp, Y.temp, folds = 2,
                                  nonzero_a[k], nonzero_b[k],
                                  model = model, lmeformula = lmeformula, silent = silent)

        #  for (s in unique(X.temp$time)) {
        #    w = unique(intersect(X.temp[X.temp$time==s,]$id, Y.temp[Y.temp$time==s,]$id))
        #    g = as.matrix(X.temp[w, -c(1,2)])%*%res_perm[[k]]$alpha; e = as.matrix(Y.temp[w, -c(1,2)])%*%res_perm[[k]]$beta;
        #    cc_time[s, k] = abs(cor(g, e))
        #  }
        #
        # m= colMeans(cc_time)

      }
      perm[d, ] = sapply(1:K, function(k) res_perm[[k]]$cancor)
    }
  }
  if (isTRUE(parallel_logic))
    parallel::stopCluster(cl = myCluster)
  margin = 0.05
  # cc_time = matrix(NA, length(unique(XX2$time)), K)
  # for (k in 1:K) {
  #   for (s in unique(A$time)) {
  #   w = unique(intersect(XX2[XX2$time==s,]$id, YY2[YY2$time==s,]$id))
  #   g = as.matrix(XX2[w, -c(1,2)])%*%res_perm[[k]]$alpha; e = as.matrix(YY2[w, -c(1,2)])%*%res_perm[[k]]$beta;
  #   cc_time[s, k] = abs(cor(g, e))
  #   }}
  m = (cancor) #sapply(1:K, function(k)mean(cc_time[,k])) #
  testStatistic = m #sapply(1:K, function(k) res_perm[[k]]$cancor) # CCAtStat(cancor, A, B, C = nuisanceVar, type = testStatType)[["tStatistic"]]
  h = hist(perm[, getWhich(testStatistic, max)], breaks = draws/2,
           plot = FALSE)
  if (is.null(bootCCA)) {
    permDensity = density(perm[, getWhich(testStatistic,
                                          max)])
    xlim = c(min(min(perm) - margin, testStatistic - margin),
             max(max(perm) + margin, testStatistic + margin))
    ylim = c(min((modes(permDensity)$y), min(perm[, getWhich(testStatistic,
                                                             max)])), max((modes(permDensity)$y) + margin, max(h$counts) +
                                                                            margin))

    oldpar <- par(no.readonly = TRUE) # code line i
    on.exit(par(oldpar))

    par(mfrow = c(1, 1))
    #   hist(perm[, getWhich(testStatistic, max)], breaks = draws/2,
    #        xlab = "Canonical Correlation", xlim = xlim, ylim = ylim,
    #        main = paste0("Distribution under de Null - ", testStatType,
    #                      " Statistic"), col = scales::alpha("#01768c",
    #                                                         0.2))
    #   lines(permDensity, col = "black", lwd = 2)
    #   EnvStats::epdfPlot(perm[, getWhich(testStatistic, max)],
    #                      discrete = FALSE, density.arg.list = NULL, plot.it = TRUE,
    #                      add = TRUE, epdf.col = "steelblue", epdf.lwd = 3 *
    #                        par("cex"), epdf.lty = 1, curve.fill = FALSE,
    #                      curve.fill.col = "steelblue", main = NULL, xlab = NULL,
    #                      ylab = NULL)
    #   abline(v = testStatistic, col = "#e86d07", lwd = 2)
    #   legend("topleft", c("Empirical pdf", "density", "model canCor"),
    #          col = c("steelblue", "black", "red"), lty = c(2,
    #                                                        1, 1), cex = 0.8)
    #   text(x = as.character(testStatistic), y = 0.9 * par("usr")[4],
    #        labels = as.character(1:K), cex = 0.9)

    col = mpalette[3:7] # viridis::magma(5)
    # library(ggplot2)
    # library(viridis)
    # library(data.table)
    # library(scales)

    # Assuming perm is a data.table or matrix, and testStatistic, margin, h, draws, testStatType, K are defined

    # Extract the relevant values
    perm_values <- perm[, getWhich(testStatistic, max)]
    testStat_val <- testStatistic

    # Compute density
    perm_density <- density(perm_values)

    # Create a data frame for plotting
    df <- data.frame(
      x = perm_values
    )

    density_df <- data.frame(
      x = perm_density$x,
      y = perm_density$y
    )

    # Compute histogram manually for better control in ggplot2
    hist_data <- hist(perm_values, breaks = draws / 2, plot = FALSE)
    hist_df <- data.frame(
      x = hist_data$mids,
      y = hist_data$counts
    )

    # Define x and y limits
    xlim <- c(min(min(perm) - margin, testStat_val - margin),
              max(max(perm) + margin, testStat_val + margin))

    ylim <- c(min(c(min(perm_values), min(perm_density$y))),
              max(c(max(hist_data$counts), max(perm_density$y))) + margin)

    # Plot

    plt= ggplot2::ggplot() +
      ggplot2::geom_histogram(ggplot2::aes(x = x, y = ..count..),
                              data = df,
                              bins = draws / 2,
                              fill = mpalette[5], # viridis::viridis(1, option = "magma", alpha = 0.2),
                              color = NA) +
      # geom_line(data = density_df, aes(x = x, y = y), color = col[5], size = 1) +
      ggplot2::geom_line(data = density_df, ggplot2::aes(x = x, y = y), color = col[1], size = 1.2, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = testStat_val, color = col[5], size = 1) +
      ggplot2::labs(x = "Canonical Correlation",
                    y = "Density / Count",
                    title = paste0("")) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
      ggplot2::theme_minimal(base_size = 14) +
      # ggplot2::scale_fill_viridis_d(option = "magma") +
      ggplot2::scale_fill_manual(values = mpalette) +
      ggplot2::annotate("text", x = testStat_val, y = 0.9 * ylim[2], label = as.character(1:K), size = 4) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
      ggplot2::theme(legend.position = "none")
    print(plt)

  }
  else {
    li = bootCCA$ccaInterval[1]
    ui = bootCCA$ccaInterval[2]
    xlim = c(min(min(perm) - margin, modelCanCor - margin,
                 li - margin), max(max(perm) + margin, modelCanCor +
                                     margin, ui + margin))
    li = bootCCA$ccaInterval[1]
    ui = bootCCA$ccaInterval[2]
    xlim = c(min(min(perm) - margin, modelCanCor - margin),
             max(max(perm) + margin, modelCanCor + margin))
    permDensity = density(perm[, getWhich(testStatistic,
                                          max)])

    oldpar <- par(no.readonly = TRUE) # code line i
    on.exit(par(oldpar))

    par(mfrow = c(1, 1))
    hist(perm[, getWhich(testStatistic, max)], breaks = draws/2,
         xlab = "Canonical Correlation", xlim = xlim, ylim = sort(modes(permDensity)$y),
         main = paste0("Distribution under de Null - ", testStatType,
                       " Statistic"), col = "#93D9D9")
    lines(permDensity, col = "black", lwd = 2)
    lines(x = c(li, li), y = c(max(histNullCCA$density) *
                                 0.5/2, max(histNullCCA$density) * 0.6/2))
    lines(x = c(ui, ui), y = c(max(histNullCCA$density) *
                                 0.5/2, max(histNullCCA$density) * 0.6/2))
    lines(x = c(li, ui), y = rep((max(histNullCCA$density) *
                                    0.5/2 + max(histNullCCA$density) * 0.6/2)/2, 2))
    abline(v = modelCanCor, col = "purple", lwd = 2)
    abline(v = mean(bootCCA$boostCanCor), col = "red", lwd = 2)
    legend("topleft", c("normal", "density", "model canCor",
                        "bootstrap avg canCor"), col = c("steelblue", "black",
                                                         "purple", "red"), lty = c(2, 2, 1, 1), cex = 0.8)
    text(x = as.character(testStatistic), y = 0.9 * par("usr")[4],
         labels = as.character(1:K), cex = 0.9)
  }
  pValues = sapply(1:K, function(k) round(mean(testStatistic[k] <=
                                                 (perm[, getWhich(testStatistic, max)])), 6))
  if(isFALSE(silent)) message(cat("Empirical p-values:", pValues, sep = "\n"))
  return(list(perm_estimate = perm, p_values = pValues))
}
