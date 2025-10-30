#' Covariance kernel with squared-exponential for \code{gpSphere()}
#'
#' @description
#' A \pkg{nimble} function that constructs a covariance matrix on
#' \eqn{t = X'\theta} using a squared–exponential Gaussian kernel with amplitude
#' \eqn{\eta} and length-scale parameter \eqn{l}.
#' Each entry is defined as \eqn{K_{ij} = \eta\exp\{-(t_i - t_j)^2 / l\}, \quad i, j = 1, \ldots, n},
#' symmetrized explicitly and stabilized with a small diagonal jitter term.
#'
#' @details
#' For the squared–exponential kernel construction, the covariance matrix is symmetrized
#' using \eqn{(K + K') / 2} and a small jitter term (\eqn{10^{-4}}) is added
#' to the diagonal to ensure positive-definiteness and numerical stability.
#' The parameters \code{amp} and \code{l} jointly control the amplitude
#' (vertical scale) and smoothness (horizontal scale) of the process.
#'
#' @param vec Numeric vector of input locations. \eqn{t = X'\theta} is the main input value for the single-index model.
#' @param l Positive numeric scalar controlling the length-scale of the kernel.
#'   Larger \code{l} yields slower decay of correlations.
#' @param amp Non-negative numeric scalar specifying the amplitude (variance scale)
#'   of the kernel.
#'
#' @return
#' A numeric \eqn{n \times n} covariance matrix with entries
#' \eqn{K_{ij} = \eta\,\exp\{-(t_i - t_j)^2 / l\}}
#' for \eqn{i, j = 1, \cdots, n}, symmetrized and stabilized with a
#' diagonal jitter term \eqn{10^{-4}}.
#'
#' @seealso
#' \code{\link{gpSphere}}, \code{\link{predict.bsimGp}}
#'
#' @export
expcov_gpSphere <- nimbleFunction(
  run = function(vec = double(1), l = double(0), amp = double(0)){
    returnType(double(2))
    n = length(vec)
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n){
      for(j in 1:n){
        result[i, j] <- amp * exp(-(vec[i] - vec[j])^2/l)
      }
    }
    result = (result + t(result))/2 + diag(rep(1e-4, n))

    return(result)

  }
)

# test expcov (n1 * n2)
expcovTest_gpSphere <- nimbleFunction(
  run = function(vec1 = double(1), vec2 = double(1), l = double(0), amp = double(0)){
    returnType(double(2))
    n1 <- length(vec1)
    n2 <- length(vec2)
    result <- matrix(nrow = n1, ncol = n2, init = FALSE)
    for(i in 1:n1){
      for(j in 1:n2){
        result[i, j] <- amp * exp(-(vec1[i] - vec2[j])^2/l)
      }
    }

    return(result)

  }
)

# beta
conBeta <- nimbleFunction(
  run = function(beta = double(1)){
    returnType(double(1))
    beta <- beta/sqrt(sum(beta^2))
    if (beta[1] < 0) beta <- beta*(-1)
    return (beta)
  }
)


# MAP
obj_btt <- nimbleFunction(
  setup = function(model, p, shape, rate,
                   a_amp, b_amp, a_sig, b_sig) {
    X <- model$X
    Y <- model$Y
  },
  run = function(par = double(1)) {
    returnType(double(0))
    # nimCat("params: ",par)
    # nimPrint("")
    index <- par[1:p]
    lengthscale  <- exp(par[(p+1)])
    amp    <- exp(par[(p+2)])
    # if (sigmaTrue) sigma2  <- exp(par[p+3])
    index <- index/sqrt(sum(index^2))
    # nimCat("index: ",index)

    if (abs(lengthscale) <= 1e-10){
      lengthscale <- lengthscale + 1e-6
    }

    n = nimDim(X)[1]
    Xlin <- numeric(n)
    for (i in 1:n){
      Xlin[i] <- sum(X[i,] * index)
    }

    # nimCat("sigma2:", model$sigma2)
    An <- expcov_gpSphere(Xlin, lengthscale, amp) +
      diag(rep(model$sigma2 + 1e-5, n))


    # if (!sigmaTrue){
    #
    # } else{
    #   An <- expcov_gpSphere(Xlin, lengthscale, amp) +
    #     diag(rep(sigma2 + 1e-10, n))
    # }
    y <- matrix(Y, ncol = 1)
    obj1 <- - (1/2) * t(y) %*% solve(An, y)
    # nimCat("obj1:",obj1)
    # nimPrint("")
    R <- chol(An)
    logdet_An <- 2 * sum(log(diag(R)))
    obj2 <- - logdet_An/2
    # nimCat("obj2:",obj2)
    # nimPrint("")
    objs <- obj1[1,1] + obj2
    # nimCat("lik:",objs)
    # nimPrint("")
    lp <- dunitSphere(index, dim = p, log = TRUE) +
      dgamma(lengthscale,shape, rate, log = TRUE) +
      dlnorm(amp, a_amp ,b_amp, log = TRUE)

    # if (sigmaTrue){
    #   lp <- lp + dinvgamma(sigma2, a_sig, b_sig, log = TRUE)
    # }

    # nimCat("prior:",lp)
    # nimPrint("")
    # nimPrint(-objs-lp)

    return(-objs-lp)
  }
)

# MAP
obj_btt_EB <- nimbleFunction(
  setup = function(X, Y, p, shape, rate,
                   a_amp, b_amp, a_sig, b_sig) {
  },
  run = function(par = double(1)) {
    returnType(double(0))
    # nimCat("params: ",par)
    # nimPrint("")
    index <- par[1:p]
    lengthscale  <- exp(par[(p+1)])
    amp    <- exp(par[(p+2)])
    sigma2  <- exp(par[p+3])
    # if (sigmaTrue) sigma2  <- exp(par[p+3])
    index <- index/sqrt(sum(index^2))
    # nimPrint(index)

    n = nimDim(X)[1]
    Xlin <- numeric(n)
    for (i in 1:n){
      Xlin[i] <- sum(X[i,] * index)
    }

    An <- expcov_gpSphere(Xlin, lengthscale, amp) +
      diag(rep(sigma2 + 1e-3, n))

    # if (!sigmaTrue){
    #
    # } else{
    #   An <- expcov_gpSphere(Xlin, lengthscale, amp) +
    #     diag(rep(sigma2 + 1e-10, n))
    # }
    y <- matrix(Y, ncol = 1)
    obj1 <- - (1/2) * t(y) %*% solve(An, y)
    # nimCat("obj1:",obj1)
    # nimPrint("")
    R <- chol(An)
    logdet_An <- 2 * sum(log(diag(R)))
    obj2 <- - logdet_An/2
    # nimCat("obj2:",obj2)
    # nimPrint("")
    objs <- obj1[1,1] + obj2
    # nimCat("lik:",objs)
    # nimPrint("")
    lp <- dunitSphere(index, dim = p, log = TRUE) +
      dgamma(lengthscale,shape, rate, log = TRUE) +
      dlnorm(amp, a_amp ,b_amp, log = TRUE)

    # if (sigmaTrue){
    #   lp <- lp + dinvgamma(sigma2, a_sig, b_sig, log = TRUE)
    # }

    # nimCat("prior:",lp)
    # nimPrint("")
    # nimPrint(-objs-lp)

    return(-objs-lp)
  }
)

## initfunction
initFunction_gpSphere <- function(N, p, X, Y,
                                index, lengthscale,
                                amp, sigma2, method,
                                setSeed){
  # init$index, init$lengthscale, init$amp, init$sigma2
  # N, p

  if (setSeed != FALSE){
    set.seed(setSeed)
  }

  # index
  init_index <- numeric(p)
  if (is.null(index)) {
    init_index <- rnorm(p)
  } else {
    stopifnot(length(index) == p)
    init_index <- index
  }
  init_index <- conBeta(init_index)
  init_Xlin <- as.vector(X %*% matrix(init_index, nrow = p))

  # lengthscale
  init_lengthscale <- lengthscale
  if (length(init_lengthscale) >= 2 || !is.numeric(init_lengthscale)){
    stop("Error: Initial value of lengthscale should be scalar.")
  }
  if (init_lengthscale < 0){
    stop("Error: Initial value of lengthscale should be positive.")
  }

  # amp (eta)
  init_amp <- amp
  if (length(init_amp) >= 2 || !is.numeric(init_amp)){
    stop("Error: Initial value of amp should be scalar.")
  }
  if (init_amp < 0){
    stop("Error: Initial value of amp should be positive.")
  }
  init_log_amp <- log(amp)

  init_cov <- expcov_gpSphere(init_Xlin, init_lengthscale, init_amp)
  init_linkFunction <- mvtnorm::rmvnorm(1, rep(0, N), sigma = init_cov)[1,]

  # sigma2
  init_sigma2 <- sigma2
  if (length(init_sigma2) >= 2 || !is.numeric(init_sigma2)){
    stop("Error: Initial value of sigma2 should be scalar.")
  }
  if (init_sigma2 < 0){
    stop("Error: Initial value of sigma2 should be positive.")
  }
  init_Simga <- init_sigma2 * diag(1, N)

  # EB : MAP


  if (method == "FB"){
    return(list(index0 = init_index, index = init_index,
                lengthscale = init_lengthscale, amp = init_amp,
                Xlin = init_Xlin,cov = init_cov,
                sigma2 = init_sigma2,
                linkFunction = init_linkFunction, Sigma = init_Simga))
  } else{
    return(list(index = init_index,
                lengthscale = init_lengthscale, amp = init_amp,
                Xlin = init_Xlin,cov = init_cov,
                sigma2 = init_sigma2,
                linkFunction = init_linkFunction, Sigma = init_Simga))
  }



}

## predict test function
pred_gpSphere <- nimbleFunction(
  run = function(newdata = double(2), nsamp = integer(0), y = double(1),
                 indexSample = double(2),
                 XlinSample = double(2), sigma2_samples = double(1),
                 lengthSample = double(1), ampSample = double(1)){
    returnType(double(2))

    new_ncol <- nimDim(newdata)[1]
    orig_ncol <- length(y)

    testPred <- matrix(0, nrow = nsamp, ncol = new_ncol)
    for (i in 1:nsamp){
      currLength <- lengthSample[i]
      currAmp <- ampSample[i]
      currSigma <- sigma2_samples[i]
      # 1. compute index
      sampleZ <- newdata %*% matrix(indexSample[i, ], ncol = 1)

      # 2. compute covariance
      ## 1) C(ori, ori)
      cov_ori_ori <- expcov_gpSphere(XlinSample[i,], currLength, currAmp)
      ## 2) C(new, ori)
      cov_new_ori <- expcovTest_gpSphere(sampleZ[,1], XlinSample[i,], currLength, currAmp)
      ## 3) c(new, new)
      cov_new_new <- expcov_gpSphere(sampleZ[,1], currLength, currAmp)

      # 3. compute mu, Sigma
      midMatrix <- inverse(cov_ori_ori + diag(rep(currSigma, orig_ncol)))
      mu <- cov_new_ori %*% midMatrix %*% matrix(y, ncol = 1)
      Sigcov <- cov_new_new - cov_new_ori %*% midMatrix %*%  t(cov_new_ori)

      # 4. Sampling
      predcov <- Sigcov + diag(rep(currSigma, new_ncol))
      cholpredcov <- chol(predcov)
      testPred[i,] <- t(rmnorm_chol(1, mean = mu[,1], cholesky = cholpredcov,
                                  prec_param  = FALSE))
    }
    return(testPred)

  }
)

