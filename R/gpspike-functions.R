
#' Covariance kernel with squared-exponential and unit nugget for \code{gpSpike()}
#'
#' @description
#' A \pkg{nimble} function that constructs a covariance matrix on the \eqn{t = X'\theta}
#'  using a squared–exponential Gaussian kernel scaled by
#' \eqn{\lambda^{-1}}, with an added unit nugget on the diagonal.
#' Covariance matrix is \eqn{I + K_{ij}} where \eqn{K_{ij} = \lambda^{-1}\exp\{-(\mathrm{t}_i - \mathrm{t}_j)^2\}} for \eqn{i, j = 1, \cdots, n}.
#'
#' @details
#' The off-diagonal structure follows the squared–exponential kernel, producing
#' rapidly decaying correlations as the squared distance grows. The matrix is
#' filled in a symmetric manner and then a unit nugget \eqn{I} is added to
#' the diagonal. The parameter \code{invlambda} controls the overall signal
#' scale of the kernel component. If a different nugget is desired, adjust externally.
#'
#' @param vec Numeric vector of input locations. \eqn{t = X'\theta} is the main input value for the single-index model.
#' @param invlambda Non-negative numeric scalar scaling the kernel amplitude.
#'   Larger values increase both diagonal and off-diagonal
#'   entries proportionally.
#'
#' @return
#' A numeric \eqn{n \times n} covariance matrix with entries \eqn{K_{ij} = \lambda^{-1}\,\exp\{-(\mathrm{t}_i - \mathrm{t}_j)^2\}}
#' for \eqn{i \ne j}, and \eqn{K_{ii} = \lambda^{-1} + 1} where \eqn{i, j = 1, \cdots, n}.
#'
#' @seealso
#' \code{\link{gpSpike}}, \code{\link{predict.bsimGp}}
#'
#' @export
# expcov: I + lambda^(-1)*K
expcov_gpSpike <- nimbleFunction(
  run = function(vec = double(1), invlambda= double(0)){
    returnType(double(2))
    n = length(vec)
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        result[i, j] <- exp(-(vec[i] - vec[j])^2) * invlambda
        result[j, i] <- result[i, j]
      }
      result[i, i] <- invlambda
    }
    result[n, n] <- invlambda
    result = result + diag(rep(1, n))

    return(result)

  }
)


## 2. pred cov - inverse lambda kernel
expcovnn_gpSpike <- nimbleFunction(
  run = function(vec1 = double(1), vec2 = double(1), invlambda= double(0)){
    returnType(double(2))
    n <- length(vec1); k <- length(vec2)
    result <- matrix(nrow = n, ncol = k, init = FALSE)
    for(i in 1:n){
      for(j in 1:k){
        result[i, j] <- exp(-(vec1[i] - vec2[j])^2) * invlambda
      }
    }
    # result = (result + t(result))/2

    return(result)

  })



computeA1 <- nimbleFunction(
  run = function(y = double(2), expcov = double(2)){
    returnType(double(0))
    computeA1 <- t(y) %*% inverse(expcov) %*% y

    return(computeA1[1,1])
  }
)


computeA2 <- nimbleFunction(
  run = function(y = double(2), z = double(2), expcov = double(2)){
    returnType(double(0))
    invCov <- inverse(expcov)
    part1 <- t(y) %*% invCov %*% z
    part2 <- inverse(t(z) %*% invCov %*% z)
    part3 <- t(z) %*% invCov %*% y

    computeA2 <- part1 %*% part2 %*% part3
    result <- computeA2[1,1]
    return(result)
  }
)
## 2. lambda posterior log - likelihood
llFunLambda <- nimbleFunction(
  setup = function(model){
    a_sig <- model$getConstants()$a_sig
    b_sig <- model$getConstants()$b_sig
    a_lam <- model$getConstants()$a_lam
    b_lam <- model$getConstants()$b_lam
    N <- model$getConstants()$N
    # pz <- model$getConstants()$z
  },
  run = function(){
    returnType(double(0))
    invlambda <- model$invlambda # target
    # z <- model$Z
    y <- model$Y
    expcovtemp <- model$Ki
    sqrt_expcovtemp <- log(diag(chol(expcovtemp)))

    # matPart3 <- t(z) %*% inverse(expcovtemp) %*% z
    # sqrt_matPart3 <- log(diag(chol(matPart3)))


    part1 <- log(invlambda) + dgamma(invlambda, shape = a_lam, rate = b_lam, log = TRUE)
    part2 <- 0
    # part3 <- 0
    for (i in 1:N){
      part2 <- part2 + sqrt_expcovtemp[i]
    }
    # for (j in 1:pz){
    #   part3 <- part3 + sqrt_matPart3[j]
    # }
    part4 <- -(a_sig + (N)/2)*log(b_sig + 0.5*(computeA1(y, expcovtemp)))

    total <- part1 - part2 + part4

    # nimPrint(total)

    return(total[1])
  }
)

# 3. theta and delta posterior log - likelihood
llFunThetaV1 <- nimbleFunction(
  setup = function(model){
    p <- model$getConstants()$p
    a0 <- model$getConstants()$a0
    b0 <- model$getConstants()$b0
    a_sig <- model$getConstants()$a_sig
    b_sig <- model$getConstants()$b_sig
    a_lam <- model$getConstants()$a_lam
    b_lam <- model$getConstants()$b_lam
    theta_sig <- model$getConstants()$sigma_theta
    N <- model$getConstants()$N
    # pz <- model$getConstants()$z
  },
  run = function(idx = integer(0)){
    returnType(double(0))
    # z <- model$Z
    y <- model$Y
    vi <- model$nu; sumvi <- 0
    theta_raw <- model$index_raw
    for (i in 1:length(vi)){
      sumvi <- sumvi + vi[i]
    }
    expcovtemp <- model$Ki
    sqrt_expcovtemp <- log(diag(chol(expcovtemp)))
    # matPart3 <- t(z) %*% inverse(expcovtemp) %*% z
    # sqrt_matPart3 <- log(diag(chol(matPart3)))

    part1 <- lgamma(sumvi + a0)+lgamma(p-sumvi+b0)
    # print(part1)
    part2 <- dnorm(theta_raw[idx], 0, theta_sig, log = TRUE) * vi[idx]
    # print(part2)
    part3 <- 0;
    # part4 <- 0
    for (i in 1:N){
      part3 <- part3 + sqrt_expcovtemp[i]
    }
    # print(part3)
    # for (j in 1:pz){
    #   part4 <- part4 + sqrt_matPart3[j]
    # }
    # print(part4)
    part5 <- -(a_sig + (N)*0.5)*log(b_sig + 0.5*(computeA1(y, expcovtemp)))
    # print(part5)
    total <- part1 + part2 - part3 + part5
    # total <- part2 - part3 - part4 + part5

    # nimPrint(total)

    return(total)
  }
)

llFunThetaV2 <- nimbleFunction(
  setup = function(model){
    p <- model$getConstants()$p
    a0 <- model$getConstants()$a0
    b0 <- model$getConstants()$b0
    a_sig <- model$getConstants()$a_sig
    b_sig <- model$getConstants()$b_sig
    a_lam <- model$getConstants()$a_lam
    b_lam <- model$getConstants()$b_lam
    theta_sig <- model$getConstants()$sigma_theta
    N <- model$getConstants()$N
    # pz <- model$getConstants()$z
  },
  run = function(idx = integer(0)){
    returnType(double(0))
    # z <- model$Z
    y <- model$Y
    vi <- model$nu; sumvi <- 0
    theta_raw <- model$index_raw
    for (i in 1:length(vi)){
      sumvi <- sumvi + vi[i]
    }
    expcovtemp <- model$Ki
    sqrt_expcovtemp <- log(diag(chol(expcovtemp)))
    # matPart3 <- t(z) %*% inverse(expcovtemp) %*% z
    # sqrt_matPart3 <- log(diag(chol(matPart3)))

    # part1 <- lgamma(sumvi + a0)+lgamma(p-sumvi+b0)
    part2 <- dnorm(theta_raw[idx], 0, theta_sig, log = TRUE)
    # print(part2)
    part3 <- 0;
    # part4 <- 0
    for (i in 1:N){
      part3 <- part3 + sqrt_expcovtemp[i]
    }
    # print(part3)
    # for (j in 1:pz){
    #   part4 <- part4 + sqrt_matPart3[j]
    # }
    # print(part4)
    part5 <- -(a_sig + (N)*0.5)*log(b_sig + 0.5*(computeA1(y, expcovtemp)))
    # print(part5)
    # total <- part1 + part2 - part3 - part4 + part5
    total <- part2 - part3 + part5

    # nimPrint(total)

    return(total)
  }
)


# 4. transition kernel

transitionTheta <- nimbleFunction(
  run = function(theta = double(0), theta_sig = double(0), delta = integer(0), sumv = double(0)){
    returnType(double(0))
    part1 <- dnorm(theta, 0, theta_sig, log = TRUE) * delta
    part2 <- 0
    if (sumv == 0){
      part2 <- log(1)
    } else{
      part2 <- log(0.5)
    }

    total <- part1 + part2
    return(total)
  }
)


# 5. multiply scalar in matrix
multiplyMatrixByConstant <- nimbleFunction(
  run = function(mat = double(2), c = double(0)) {
    returnType(double(2))

    nrow <- nimDim(mat)[1]
    ncol <- nimDim(mat)[2]
    result <- matrix(0, nrow, ncol)

    for(i in 1:nrow){
      for(j in 1:ncol){
        result[i, j] <- mat[i, j] * c
      }
    }
    return(result)
  }
)

# Initial value functions
initfunction_gpSpike <- function(x, y, p,
                                 pi, nu, index, inv_lambda, sigma2,
                                 setSeed){
  if (setSeed != FALSE){
    set.seed(setSeed)
  }
 # Initialize
  ## pi, nu, theta
  init_pi <- pi
  if (length(pi) >= 2 || !is.numeric(pi)){
    stop("Error: Initial value of pi should be scalar.")
  }
  if (pi < 0){
    stop("Error: Initial value of pi should be positive.")
  }

  ## nu, theta
  if (!is.null(nu) & (length(nu) != p)){
    stop("Error: Incorrect dimension on initial value of nu.")
  }

  if (!is.null(nu) & (!all(nu %in% c(0, 1)))){
    stop("Error: Delta should consist of 0 or 1.")
  }

  if (is.null(nu) & is.null(index)){
    init_nu <- rbinom(p, 1 ,prob = init_pi)
  } else if (is.null(nu) & !is.null(index)){
    init_nu <- as.numeric(index!=0)
  } else{
    init_nu <- nu
  }


  if (!is.null(index) & (length(index) != p)){
    stop("Error: Incorrect dimension on initial value of index")
  }

  if (is.null(index)){
    nonzero <- sum(init_nu)
    nonzeroidx <- which(init_nu == 1)
    temp <- rtruncnorm(n = sum(nonzeroidx), 0, 1, lower = 0, upper = Inf)
    init_index <- rep(0, p)
    init_index[nonzeroidx] <-temp[nonzeroidx]
  } else{
    bad_zero  <- which(init_nu == 0 & index != 0)
    bad_nonzero <- which(init_nu == 1 & index == 0)

    if (length(bad_zero) > 0) {
      warning(sprintf(
        "Theta is nonzero at positions where nu = 0: %s",
        paste(bad_zero, collapse = ", ")
      ))
    }

    if (length(bad_nonzero) > 0) {
      warning(sprintf(
        "Theta is zero at positions where nu = 1: %s",
        paste(bad_nonzero, collapse = ", ")
      ))
    }

    init_index <- index * init_nu
  }

  if (sum(init_index) < 0){
    init_index <- (-1) * init_index
  }

  init_Xlin <- as.vector(x %*% matrix(init_index, nrow = p))

  ## inv_lambda
  init_invlambda <- ifelse(is.null(inv_lambda), exp(rnorm(1)*sqrt(1)/10), inv_lambda)

  if (length(init_invlambda) >= 2 || !is.numeric(init_invlambda)){
    stop("Error: Initial value of inverse lambda should be scalar.")
  }
  if ((init_invlambda < 0)){
    stop("Error: Initial value of inverse lambda should be positive.")
  }


  ## sigma2
  init_sigma2 <- ifelse(is.null(sigma2), 1, sigma2)
  if (length(init_sigma2) >= 2 || !is.numeric(init_sigma2)){
    stop("Error: Initial value of sigma2 should be scalar.")
  }
  if (init_sigma2 < 0){
    stop("Error: Initial value of sigma2 should be positive.")
  }

  init_Ki <- expcov_gpSpike(init_Xlin, init_invlambda)
  init_cov <- init_Ki * init_sigma2

  return(list(pi = init_pi, nu = init_nu,
              index_raw = init_index,
              indexstar = init_index, Xlin = init_Xlin,
              Ki = init_Ki, invlambda = init_invlambda,
              sigma2 = init_sigma2, cov = init_cov))
}


# Prediction
pred_gpSpike <- nimbleFunction(
  run = function(newdata = double(2), nsamp = integer(0), y = double(1),
                 indexSample = double(2),
                 XlinSample = double(2), sigma2_samples = double(1),
                 invlambdaSample = double(1)){
    returnType(double(2))

    new_ncol <- nimDim(newdata)[1]
    orig_ncol <- length(y)

    testPred <- matrix(0, nrow = nsamp, ncol = new_ncol)
    for (i in 1:nsamp){
      currInvlambda <- invlambdaSample[i]
      currSigma <- sigma2_samples[i]
      # 1. compute index
      sampleZ <- newdata %*% matrix(indexSample[i, ], ncol = 1)

      # 2. compute covariance
      ## 1) C(ori, ori) - I + invlambda*Koo
      cov_ori_ori <- expcov_gpSpike(XlinSample[i,], currInvlambda)
      ## 2) C(new, ori) - Kno * invlambda
      cov_new_ori <- expcovnn_gpSpike(sampleZ[,1], XlinSample[i,], currInvlambda)
      ## 3) c(new, new) - knn * invlambda
      cov_new_new <- expcovnn_gpSpike(sampleZ[,1], sampleZ[,1], currInvlambda)

      # 3. compute mu, Sigma
      midMatrix <- inverse(cov_ori_ori)
      mu <- cov_new_ori %*% midMatrix %*% matrix(y, ncol = 1)
      Sigcov <- cov_new_new - cov_new_ori %*% midMatrix %*%  t(cov_new_ori)

      # 4. Sampling
      predcov <- Sigcov * currSigma + diag(rep(currSigma, new_ncol))
      cholpredcov <- chol(predcov)
      testPred[i,] <- t(rmnorm_chol(1, mean = mu[,1], cholesky = cholpredcov,
                                    prec_param  = FALSE))
    }
    return(testPred)
  }
)
