
dKnotsSimple <- nimble::nimbleFunction(
  run = function(x = double(1), a = double(0), b = double(0), k = double(0),
                 alpha = double(1), log = integer(0, default = 0)) {
    returnType(double(0))
    K <- k
    real_alpha <- alpha[1:K]

    gaps <- numeric(K + 1)
    gaps[1] <- x[1] - a
    if (K > 1){
      for (i in 2:K) {
        gaps[i] <- x[i] - x[i-1]
      }
      gaps[K+1] <- b - x[K]
    }
    gaps <- gaps / (b - a)

    alpha_full <- c(real_alpha, 1.0)
    logProb <- ddirch(gaps, alpha_full, log = TRUE)
    if (log) return(logProb)
    else return(exp(logProb))
  }
)


rKnotsSimple <- nimble::nimbleFunction(
  run = function(n = integer(0), a = double(0), b = double(0), k = double(0),
                 alpha = double(1)) {
    returnType(double(1))
    K <- k
    real_alpha <- alpha[1:K]
    alpha_full <- c(alpha, 1.0)
    gaps <- rdirch(1, alpha_full)
    out <- numeric(K)
    cs <- 0.0
    for (i in 1:K) {
      cs <- cs + gaps[i]
      out[i] <- a + cs * (b - a)
    }
    return(out)
  }
)



dunitSphere <- nimble::nimbleFunction(
  run = function(x = double(1), dim =double(), log = integer(0, default = 0)){
    returnType(double(0))

    norm_x <- 0
    for (i in 1:dim) {
      norm_x <- norm_x + x[i] * x[i]
    }
    norm_x <- sqrt(norm_x)

    ## Uniform Density on the Unit Sphere
    density <- (2*pi^(length(x)/2))/gamma(length(x)/2)

    if(log) return(log(density))
    return(density)
  }
)

runitSphere <- nimble::nimbleFunction(
  run = function(n = integer(0),dim = double()) {
    returnType(double(1))
    if (dim <= 0) return(numeric(0))

    z <- rnorm(dim, 0, 1)
    z <- z / sqrt(sum(z * z))
    return(z)
  }
)

transX_sp <- nimble::nimbleFunction(
 run = function(Xlin = double(1), degree = integer(0), knots = double(1),
                 k = integer(0), maxBasis = integer(0), a_alpha = double(0),
                 b_alpha = double(0)){
    returnType(double(2))
    realKnots <- knots[1:k]
    boundary_knots <- nimC(a_alpha, b_alpha)
    X <- bsBasis(Xlin, knots = realKnots, degree = degree, intercept = TRUE,
                 boundary_knots = boundary_knots)
    N <- nimDim(X)[1]; ncol <- nimDim(X)[2]
    result <- nimMatrix(0, nrow = N, ncol = maxBasis)
    result[1:N, 1:ncol] <- X
    return(result)

  }
)

estBetaInit_sp <-  nimble::nimbleFunction(
  run = function(Xmat = double(2), Y = double(2), numBasis = integer(0)){
    returnType(double(2))
    maxBasis <- nimDim(Xmat)[2]
    Xalpha <- Xmat[,1:numBasis]
    lambda <- 1e-8
    part1 <- inverse(t(Xalpha) %*% Xalpha + diag(rep(lambda, numBasis)))
    part2 <- t(Xalpha) %*% Y
    result <- part1 %*% part2 # p*1

    total <- matrix(0, nrow = maxBasis, ncol = 1)
    total[1:numBasis, 1] <- result

    return (total)
  }
)


# (numBasis＊numBasis)
computeSig <- nimble::nimbleFunction(
  run = function(tau = double(0), Sig0 = double(2), xmat = double(2), numBasis = integer(0)){
    returnType(double(2))
    part1 <- inverse(Sig0) / tau
    part2 <- t(xmat) %*% xmat
    result <- inverse(part1 + part2)
    result2 <- result[1:numBasis, 1:numBasis]
    return(result2)
  }
)


ComputeS <- nimble::nimbleFunction(
  run = function(y = double(2), xmat = double(2), Sig = double(2), numBasis = integer(0)){
    returnType(double(0))
    corXmat <- xmat[, 1:numBasis]
    part1 <- t(y) %*% y
    part2 <- t(y) %*% corXmat
    part3 <- part2 %*% Sig %*% t(part2)
    result <- part1 - part3
    return(result[1,1])
  }
)

# About sampling - ll
logdet_nim <- nimble::nimbleFunction(
  run = function(M = double(2)) {
    returnType(double(0))
    U <- chol(M)
    out <- 2 * sum(log(diag(U)))
    return(out)
  }
)


postll_bspline_sphere <- nimble::nimbleFunction(
  setup = function(model){
    tau2 <- model$getConstants()$tau
    Sigma0 <- model$getConstants()$Sigma0
    y <- model$Y
    N <- model$getConstants()$N
  },
  run = function(){
    returnType(double(0))
    k <- model$k[1]
    numBasis <- model$numBasis
    xmat <- model$Xmat
    b_alpha <- model$b_alpha[1]
    a_alpha <- model$a_alpha[1]

    numBasisI <- as.integer(numBasis[1])
    Sigma <- computeSig(tau2, Sigma0, xmat, numBasisI)
    part1 <- logdet_nim(Sigma)/2
    part2 <- logdet_nim(Sigma0[1:numBasisI, 1:numBasisI])/2
    part3 <- log(ComputeS(y, xmat, Sigma, numBasisI)) * (-(N/2))
    part4 <- -k * log(b_alpha - a_alpha)

    result <- part1 - part2 + part3 + part4


    return(result)
  }
)



postll_knots <- nimble::nimbleFunction(
  setup = function(model){
    tau2 <- model$getConstants()$tau
    Sigma0 <- model$getConstants()$Sigma0
    y <- model$Y
    N <- model$getConstants()$N
    k <- model$k
  },
  run = function(){
    returnType(double(0))
    numBasis <- model$numBasis
    xmat <- model$Xmat
    b_alpha <- model$b_alpha
    a_alpha <- model$a_alpha

    numBasisI <- as.integer(numBasis[1])
    Sigma <- computeSig(tau2, Sigma0, xmat, numBasisI)
    part1 <- logdet_nim(Sigma)/2
    part2 <- logdet_nim(Sigma0[1:numBasisI, 1:numBasisI])/2
    part3 <- log(ComputeS(y, xmat, Sigma, numBasisI)) * (-(N/2))

    result <- part1 - part2 + part3

    return(result)
    # return(0.0)

  }
)

# Compute A1
betaFunction <- nimble::nimbleFunction(
  run = function(a = double(0), b = double(0), log = logical(0, default = TRUE)){
    returnType(double(0))
    result <- lgamma(a) + lgamma(b) - lgamma(a+b)
    if (log){
      return (result)
    } else{
      return(result)
    }
  }
)


computeA1_nu1 <- nimble::nimbleFunction(
  run = function(r1 = double(0), r2 = double(0), nnu = double(0),
                 p = double(0), r3 = double(0), r4 = double(0),
                 alpha = double(0)){
    returnType(double(0))

    part1 <- log(r2 + p - nnu) - log(nnu + r1 - 1)
    part2 <- dbeta(alpha^2, r3, r4, log = TRUE)
    part3 <- betaFunction(0.5, (nnu-1)/2)
    part4 <- -0.5*log(alpha^2)+((nnu-1)/2)*log(1-alpha^2)
    result <- part1 + part2 + part3 - part4
    return(result)
  }
)


computeA1_nu0 <- nimble::nimbleFunction(
  run = function(r1 = double(0), r2 = double(0), nnu = double(0),
                 p = double(0), r3 = double(0), r4 = double(0),
                 eta = double(0)){
    returnType(double(0))

    part1 <- log(nnu + r1) - log(r2 + p - nnu -1)
    part2 <- -dbeta(eta, r3, r4, log = TRUE)
    part3 <- -0.5*log(eta)+(nnu/2)*log(1-eta)
    part4 <- betaFunction(0.5, (nnu)/2)

    result <- part1 + part2 + part3 - part4
    return(result)
  }
)

# compute A2
computeA2_add <- nimble::nimbleFunction(
  run = function(k = double(0), minKnot = double(0), maxKnot = double(0),
                 minBound = double(0), maxBound = double(0)){
    returnType(double(0))
    part1 <- log(k+1)
    part2 <- log(maxKnot - minKnot)
    part3 <- log(maxBound - minBound)
    result <- part1 + part2 - part3
    return(result)
  }
)



computeA2_delete <- nimble::nimbleFunction(
  run = function(k = double(0), minKnot = double(0), maxKnot = double(0),
                 minBound = double(0), maxBound = double(0)){
    returnType(double(0))
    part1 <- -log(k)
    part2 <- log(maxKnot - minKnot)
    part3 <- log(maxBound - minBound)
    result <- part1 + part3 - part2
    return(result)
  }
)


# sampling truncated normal
rnormTrun <- nimble::nimbleFunction(
  run = function(mean = double(0), sd = double(0)){
    returnType(double(0))
    samp <- rnorm(1, mean, sd)

    while(samp <= -3.141593 | samp >= 3.141593){
      samp <- rnorm(1, mean, sd)
    }

    return(samp)

  }
)

# calculate knots
newKnots <- nimble::nimbleFunction(
  run = function(oldknots = double(1), olda = double(0), oldb = double(0),
                 newa = double(0), newb = double(0), k = double(0)){
    returnType(double(1))
    newknots <- numeric(nimDim(oldknots)[1])
    for (i in 1:k){
      old <- (oldknots[i] - olda) / (oldb - olda)
      newknots[i] <- newa + old * (newb - newa)
    }
    return(newknots)
  }
)

# nimSort - in bSpline2
# bk, dk, mk
prop_add <- nimble::nimbleFunction(
  run = function(lambda = double(0), k = double(0), tau = double(0)){
    returnType(double(0))
    if (k == 0){
      return(1)
    } else{
      calc <- lambda/((k+1)*sqrt(tau))
      result <- 0.4 * min(1, calc)
      return(result)
    }
  }
)


prop_delete <- nimble::nimbleFunction(
  run = function(lambda = double(0), k = double(0), tau = double(0)){
    returnType(double(0))
    if (k == 0){
      return(0)
    } else{
      calc <- (k*sqrt(tau))/lambda
      result <- 0.4 * min(1, calc)
      return(result)
    }
  }
)

# truncated normal
rtruncnorm <- nimble::nimbleFunction(
  run = function(n = integer(0),
                 mean = double(0), sigma = double(0),
                 lower = double(0), upper = double(0)) {
    returnType(double(1))

    # if (sigma <= 0) stop("sigma must be > 0")
    # if (lower >= upper) stop("need lower < upper")

    out <- numeric(n)

    zL <- (lower - mean) / sigma
    zU <- (upper - mean) / sigma

    pL <- 0.0
    if (lower <= -1.0e300) {
      pL <- 0.0
    } else {
      pL <- pnorm(zL, 0.0, 1.0)
    }

    pU <- 1.0
    if (upper >=  1.0e300) {
      pU <- 1.0
    } else {
      pU <- pnorm(zU, 0.0, 1.0)
    }

    for (i in 1:n) {
      u <- runif(1, pL, pU)
      out[i] <- mean + sigma * qnorm(u, 0.0, 1.0)
    }
    return(out)
  }
)

# initFunction_bS
initFunction_bS <- function(nu, index, k, beta, sigma2, knots,
                            p, X, Y, lambda_k, maxknots, degree, maxBasis,
                            tau, Sigma0, setSeed){
  init <- 0
  if (setSeed != FALSE){
    set.seed(setSeed)
  }

  # init_nu
  if (is.null(nu) & is.null(index)){
    init_nu <- c(1, rbinom(p-1, 1, 0.5))
  } else if (is.null(nu) & !is.null(index)){
    init_nu <- as.numeric(index != 0)
  } else{
    init_nu <- init$nu
  }

  # index
  if (is.null(index)){
    nonzero <- sum(init_nu)
    nonzeroidx <- which(init_nu == 1)
    temp <- runitSphere(1, nonzero)
    init_index <- rep(0, p)
    init_index[nonzeroidx] <-temp
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

  if (sum((init_index)^2) != 1){
    init_index <- init_index/sqrt(sum(init_index^2))
  }

  if(init_index[1] < 0){
    init_index[1] <- -init_index[1]
  }

  init_Xlin <- as.vector(X %*% matrix(init_index, ncol = 1))
  init_a_alpha <- min(init_Xlin); init_b_alpha <- max(init_Xlin)

  ## knots
  ## k (# of knots)
  if (is.null(k)){
    init_k <- rpois(1, lambda_k) + 1
  } else{
    init_k <- k
  }

  if ((!is.null(k)) & (length(init_k) >= 2 || !is.numeric(init_k))){
    stop("Error: Initial value of the number of knots(k) should be scalar.")
  }
  if ((!is.null(k)) & (init_k < 0)){
    stop("Error: Initial value of the number of knots(k) should be positive.")
  }
  if (!is.null(init_k) & !is.null(knots)){
    if(length(knots) != init_k){
      stop("Error: Number of knots and k is not the same.")
    }
  }

  ## location of knots (weights)
  if (!is.null(knots) & (knots <= 0 || knots >= 1)){
    stop("Error: Initial knots should be bigger than 0, less than 1.")
  }

  if (is.null(knots)){
    init_weights <- rdirch(1, rep(1, init_k))
  } else{
    init_weights <- knots
  }
  cumsum_w <- cumsum(init_weights); epsilon <- 0.05
  init_knots <- init_a_alpha +
    (init_b_alpha - init_a_alpha) * (epsilon + (1 - 2 * epsilon) * cumsum_w)
  init_knots <- c(init_knots, rep(0, maxknots - init_k))
  init_numBasis <- init_k + degree + 1
  init_Xmat <- transX_sp(init_Xlin, degree = degree, knots = init_knots, init_k, maxBasis, init_a_alpha, init_b_alpha)

  if (is.null(beta)){
    init_beta <- estBetaInit_sp(init_Xmat, Y, init_numBasis)
  } else{
    init_beta <- c(beta, rep(0, maxBasis - init_numBasis))
  }

  if ((!is.null(beta) & (length(init_beta) > maxBasis)) ||
      (!is.null(beta) & (length(init_beta) != init_numBasis))){
    stop("Error: Incorrect dimension on initial value of beta.")
  }

  ## sigma2
  init_sigma2 <- sigma2
  if (length(sigma2) >= 2 || !is.numeric(sigma2)){
    stop("Error: Initial value of sigma2 should be scalar.")
  }
  if (sigma2 < 0){
    stop("Error: Initial value of sigma2 should be positive.")
  }
  init_linkFunction <- init_Xmat %*%  init_beta
  init_covbeta <- init_sigma2 * tau * Sigma0
  N <- nrow(X)
  init_covY <- diag(rep(init_sigma2, N))

  return(list(nu = init_nu, index = init_index,
              Xlin = init_Xlin,a_alpha = init_a_alpha,
              b_alpha = init_b_alpha, k = init_k,
              knots = init_knots, numBasis = init_numBasis,
              Xmat = init_Xmat, sigma2 = init_sigma2,
              beta = init_beta, covBeta = init_covbeta,
              linkFunction = init_linkFunction))

}

# prediction
pred_bsplineSphere <- nimbleFunction(
  run = function(newdata = double(2), indexSample = double(2),
                 knotsSample = double(2), XlinSample = double(2),
                 betaSample = double(2), sigma2_samples = double(1),
                 k = double(1), dfSample = double(1),
                 mina = double(1), maxb = double(1), nsamp = double(0),
                 degree = double(0), prediction = integer(0)){
    returnType(double(2))

    # prediction = 1: latent, prediction = 2: response
    new_ncol <- nimDim(newdata)[1]
    testPred <- nimMatrix(0, nrow = nsamp, ncol = new_ncol)

    for (i in 1:nsamp){
      sampleZ <- newdata %*% matrix(indexSample[i, ], ncol = 1)
      numKnots <- k[i]
      Xmat <- bsBasis(sampleZ[,1],
                      degree = degree, knots = knotsSample[i,1:numKnots],
                      intercept = TRUE,
                      boundary_knots = c(mina[i], maxb[i]))
      linkPred <- Xmat %*% betaSample[i,1:dfSample[i]]
      if (prediction == 1){ # latent
        testPred[i,] <- linkPred[ ,1]
      } else if (prediction == 2){ # response
        predsigma <- chol(diag(rep(sigma2_samples[i], new_ncol)))
        pred <- t(rmnorm_chol(1, mean = linkPred[,1], cholesky = predsigma, prec_param  = FALSE))
        idxmin <- (sampleZ[,1] < (mina[i]))
        idxmax <- (sampleZ[,1] > (maxb[i]))

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


