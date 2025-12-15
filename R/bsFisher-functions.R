nimNorm <- nimble::nimbleFunction(
  run = function(vec = double(1)){
    returnType(double(1))
    sumALL <- sqrt(sum(vec^2))
    result <- vec/sumALL

    return(result)
  }
)


rW <- nimble::nimbleFunction(
  run = function(size = integer(0), lambda = double(0), d = integer(0)){
    returnType(double(1))
    resultW <- numeric(size)
    b <- d / (sqrt(4 * lambda * lambda + d * d) + 2 * lambda)
    x0 <- (1 - b) / (1 + b)
    c <- lambda * x0 + d * log(1 - x0^2)

    for(i in 1:size){
      continue_sampling <- TRUE
      while (continue_sampling) {
        Z = rbeta(1, d/2, d/2)
        w = (1 - (1 + b) * Z) / (1 - (1 - b) * Z)
        U = runif(1, 0, 1)

        if (lambda * w + d * log(1 - x0 * w) - c < log(U)) {
          resultW[i] <- w
          continue_sampling <- FALSE
        }
      }
    }

    return(resultW)
  }
)


# modified Bessel function of the first kind
besselI_nimble <- nimble::nimbleFunction(
  run = function(x = double(0), nu = double(0), tol = double(0, default = 1e-10), max_iter = integer(0, default = 1e4)) {
    returnType(double(0))

    sum <- 0.0
    m <- 0
    term <- 1.0

    while (term > tol & m < max_iter) {
      logpow_val <- (2 * m + nu) * log(x / 2)
      logdenom <- lgamma(m + 1.0) + lgamma(m + nu + 1.0)
      logterm <- logpow_val - logdenom
      term <- exp(logterm)
      sum <- sum + term
      m <- m + 1
    }

    return(sum)
  }
)


dvMFnim <- nimble::nimbleFunction(
  run = function(x = double(1), theta = double(1), log = integer(0, default = 0)) {
    returnType(double(0))
    d <- nimDim(x)[1]
    nu_1 <- (d / 2)
    kappa <- sqrt(sum(theta^2))

    ## log density
    part1 <- inprod(theta, x)
    part2 <- log(besselI_nimble(x = kappa, nu = nu_1))
    part3 <- lgamma(nu_1)
    part4 <- (nu_1 - 1) * log(kappa / 2)

    logden <- part1 - (part2 + part3 - part4)

    if (log)
      return(logden)
    else
      return(exp(logden))
  }
)


rvMFnim <- nimble::nimbleFunction(
  run = function(n = integer(0), theta = double(1)) {
    returnType(double(1))
    p <- nimDim(theta)[1]
    kappa <- sqrt(sum(theta^2))

    if (kappa == 0) {
      X <- nimMatrix(rnorm(n * p), n, p)
      for (j in 1:n)
        X[j, 1:p] <- nimNorm(X[j, 1:p])
    } else {
      d <- p - 1
      W <- rW(n, kappa, d)

      mu <- nimMatrix(nimNorm(theta), ncol = 1)
      Z <- nimMatrix(rnorm(n * p), nrow = n, ncol = p)
      Vtemp <- Z - Z %*% mu %*% t(mu)
      X <- matrix(0, nrow = n, ncol = p)
      for (i in 1:n) {
        Vtemp[i, 1:p] <- nimNorm(Vtemp[i, 1:p])
        sqrt_term <- sqrt(1 - W[i]^2)
        for (j in 1:p) {
          X[i, j] <- W[i] * mu[j, 1] +
            sqrt_term * Vtemp[i, j]
        }
      }
    }
    return(X[1, ])
  }
)


Stheta <- nimble::nimbleFunction(
  run = function(Xmat = double(2), Y = double(2), rho = double(0)){
    returnType(double(0))
    p <- nimDim(Xmat)[2]
    Y_mat <- matrix(Y, ncol = 1)

    temp_xmat <- t(Xmat) %*% Xmat
    temp_rho <- diag(nimRep(rho, p))
    Sigma_0 <- inverse(temp_xmat)
    temp_sum <- temp_xmat + temp_rho
    Sigma_rho <- inverse(temp_sum)

    part1 <- t(Y_mat) %*% Y_mat
    part2 <- t(Y_mat) %*% Xmat
    part3 <- Sigma_0 + 0.5 * (Sigma_rho %*% (diag(p) - inverse(Sigma_0) %*% Sigma_rho))

    first <- part1[1,1]
    second <- part2 %*% part3 %*% t(part2)
    result <- first - second[1,1]
    return(result)

  }
)

# marginalized posterior(D)
postll_bspline_fisher <- nimble::nimbleFunction(
  setup = function(model){
    A <- model$getConstants()$a_sig
    B <- model$getConstants()$b_sig
    n <- model$getConstants()$N
    Y <- model$Y
    index_prior <- model$getConstants()$index_prior
    consprior <- sqrt(sum(index_prior^2))
  },
  run = function(rho = double(0)){
    returnType(double(0))

    index <- model$index
    Xmat <- model$Xmat

    Stheta0 <- Stheta(Xmat, Y, rho)
    part1 <- (-A + n/2)*log(Stheta0 +(2/B))

    angle <- inprod(index, index_prior)
    part2 <- consprior * angle

    logresult <- part1 + part2
    return(logresult)

  }
)

# estBeta
estBeta_fisher <- nimble::nimbleFunction(
  run = function(Xmat = double(2), Y = double(2), rho = double(0)){
    returnType(double(1))
    p <- nimDim(Xmat)[2]
    part1 <- inverse(t(Xmat) %*% Xmat + diag(nimRep(rho, p)))
    part2 <- t(Xmat) %*% Y
    result <- part1 %*% part2
    return (result[,1])
  }
)



gvcCV <- nimble::nimbleFunction(
  run = function(Xmat = double(2), Y = double(2)){
    returnType(double(0))
    grid <- nimSeq(from = 0.1, to = 4, length.out = 20)
    n <- nimDim(Xmat)[1]
    p <- nimDim(Xmat)[2]

    for (i in 1:20){
      rho <- grid[i]
      # print(paste("time: ", i))
      beta <- estBeta_fisher(Xmat, Y, rho) ## beta
      norm <- (Y - Xmat %*% beta)^2
      RSS <- mean(norm[,1])

      Sigma_rho <- inverse(t(Xmat) %*% Xmat + diag(nimRep(rho, p)))
      Ls <- Xmat %*% Sigma_rho %*% t(Xmat) # N*N
      part1 <- (diag(n) - Ls) %*% Xmat # n*p
      part2 <- inverse(t(Xmat) %*% (diag(n) - Ls) %*% Xmat)
      part3 <- t(Xmat) %*% (diag(n) - Ls)
      L_rho <- Ls + part1 %*% part2 %*% part3 # N*N
      df <- (1-(sum(diag(L_rho))/n))^2
      currgcv <- RSS/df
      # print(currgcv)

      if (i == 1){
        mingcv <- currgcv
        optRho <- rho
      } else{
        if (currgcv < mingcv){
          # Find min
          mingcv <- currgcv
          optRho <- rho
        }
      }
    }

    return(optRho)

  }
)



# Design Matrix of bsFisher(), bsPolar(), bsSpike()
transX_fisher <- nimble::nimbleFunction(
  run = function(Xlin = double(1), df = integer(0), degree = integer(0), delta = double(0)){
    returnType(double(2))
    n_internal_knots <- df - degree
    idxmax <- n_internal_knots + 1
    # prob_vec <- nimSeq(0, 1, length.out = n_internal_knots+2)[2:idxmax]

    seq_full <- nimSeq(0, 1, length.out = n_internal_knots + 2)
    prob_vec <- seq_full[2:idxmax]
    internal_knots <- quantile_nimble(Xlin, prob_vec)

    X <- bsBasis(Xlin, df = df, degree = degree,
                 knots = internal_knots,
                 boundary_knots = c(min(Xlin)-delta, max(Xlin)+delta))
    #
    # N <- nimDim(X)[1]
    # result <- matrix(X, nrow = N, ncol = df)
    return(X)

  }
)

## initial value
initFunction_bF <- function(X, Y,
                            index, sigma2, beta,
                            df, degree, delta,
                            index_direction, setSeed){
  if (setSeed != FALSE){
    set.seed(setSeed)
  }

  # index
  if (is.null(index)){
    init_index <- rvMFnim(1, index_direction)
  } else{
    init_index <- index
  }
  if (sum((init_index)^2) != 1){
    init_index <- init_index/sqrt(sum(init_index^2))
  }

  # sigma2
  init_sigma2 <- sigma2
  init_Xlin <- as.vector(X %*% matrix(init_index, ncol = 1))
  init_Xmat <- transX_fisher(init_Xlin, df = df, degree = degree, delta = delta)
  init_rho <- gvcCV(init_Xmat, Y)

  # beta
  if (is.null(beta)){
    init_beta <- estBeta_fisher(init_Xmat, Y, init_rho)
  } else{
    init_beta <- beta
  }
  init_linkFunction <- init_Xmat %*%  init_beta

  return(list(index = init_index, Xlin = init_Xlin,
              Xmat = init_Xmat,
              sigma2 = init_sigma2, beta = init_beta,
              linkFunction = init_linkFunction))

}


pred_bsplineFisher <- nimbleFunction(
  run = function(newdata = double(2), indexSample = double(2),
                 knotsSample = double(2), XlinSample = double(2),
                 betaSample = double(2), sigma2_samples = double(1), nsamp = double(0),
                 df = double(0), degree = double(0), delta = double(0),
                 prediction = integer(0)){
    # prediction = 1: latent, prediction = 2: response
    returnType(double(2))
    new_ncol <- nimDim(newdata)[1]
    testPred <- nimMatrix(0, nrow = nsamp, ncol = new_ncol)

    for (i in 1:nsamp){
      sampleZ <- newdata %*% matrix(indexSample[i, ], ncol = 1) # linear pred
      Xmat <- bsBasis(sampleZ[,1], df = df,
                      degree = degree, knots = knotsSample[i,],
                      boundary_knots = c(min(XlinSample[i,])-delta,
                                         max(XlinSample[i,])+delta))
      linkPred <- Xmat %*% betaSample[i,]

      if (prediction == 1){ # latent
        testPred[i,] <- linkPred[,1]

      } else if (prediction == 2){
        predsigma <- chol(diag(rep(sigma2_samples[i], new_ncol)))
        pred <- t(rmnorm_chol(1, mean = linkPred[,1], cholesky = predsigma, prec_param  = FALSE))

        # computation for data over boundaries: substitute to last values
        idxmin <- (sampleZ[,1] < (min(XlinSample[i,]) - delta))
        idxmax <- (sampleZ[,1] > (max(XlinSample[i,]) + delta))

        if (sum(idxmin) != 0){
          minVal_idx <- which(XlinSample[i,] == min(XlinSample[i,]))
          predmin <- min(pred[1, minVal_idx])
          pred[1, idxmin] <- nimRep(predmin, sum(idxmin))
        }

        if (sum(idxmax) != 0){
          maxVal_idx <- which(XlinSample[i,] == max(XlinSample[i,]))
          predmax <- max(pred[1, maxVal_idx])
          pred[1, idxmax] <- nimRep(predmax, sum(idxmax))
        }

        testPred[i,] <- pred[1,]
      }

    }
    return(testPred)
  }
)
