## Define samplers
## Gibbs sigma
gibbsSampler_sigma2 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    N <- model$getConstants()$N
    a <- model$getConstants()$a
    b <- model$getConstants()$b
  },
  run = function() {
    orderf <- nimOrder(model$Xlin)# index
    y <- model$y; y <- y[orderf]
    eta <- model$linkFunction; eta <- eta[orderf]
    err <- sum((y - eta)^2)/2
    newSigma <- 1/rgamma(1, shape = a + (N/2), rate = (b + err))

    model[[target]] <<- newSigma
    model$calculate(calcNodes)

    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {})
)

## Gibbs kappa
gibbsSampler_kappa <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    kappas <- nimSeq(model$getConstants()$kappa_a,
                     model$getConstants()$kappa_b,
                     control$grid.width)
    N <- model$getConstants()$N
  },
  run = function() {
    probs <- nimRep(0.0, length(kappas))
    idx <- nimOrder(model$Xlin, decreasing = FALSE)
    eta <- nimMatrix(model$linkFunction, ncol = 1); eta <- eta[idx,] # n * 1

    for (i in 1:length(kappas)) {
      C_inv <- invcov(model$Xlin, kappas[i])
      lp1 <- matrix(0.5 * log(det(C_inv)), nrow = 1, ncol = 1)
      lp2 <- matrix(- 0.5 * (t(eta) %*% C_inv %*% eta), nrow = 1, ncol = 1)
      lp <- lp1 + lp2
      if (any(is.nan(lp))) {
        nimPrint("NaN detected in log posterior at kappa index ", i)
        lp <- matrix(1e8, nrow = 1, ncol = 1)
      }
      probs[i] <- lp[1, 1]
    }

    maxlp <- max(probs)
    probs <- exp(probs-maxlp)

    chosen <- rcat(1, probs / sum(probs))
    storeKappa <- kappas[chosen:chosen]
    model$kappa <<- storeKappa

    model$calculate(calcNodes)
    nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list(reset = function() {})
)

## Update theta/eta
MH_thetaeta <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    p <- model$getConstants()$p
    N <- model$getConstants()$N
    paramalpha <- model$getConstants()$c
    parambeta <- model$d
    obj_fn     <- obj_btt_theta(model, p-1)
    thetas <- seq(0.1, 1-0.1, 0.05)
    grid.df  <- do.call(expand.grid,
                         replicate(p-1, thetas, simplify = FALSE))
    grid.df <- as.matrix(grid.df)
  },
  run = function() {
    # 1) Old eta, theta, theta, alpha(beta)
    oldtheta <- model$psi
    oldBeta <- model$index
    orderf <- nimOrder(model$Xlin) # index
    oldeta <- model$linkFunction; oldeta <- oldeta[orderf] # sort f
    Stheta <- diag(nimRep(0.1, N))

    # surrogate data
    g <- rmnorm_chol(1, mean = oldeta, cholesky = chol(Stheta),
                     prec_param= FALSE)

    z <- nimMatrix(model$Xlin, ncol = 1)
    z <- z[orderf, ] # (n*1)
    R1 <- nimMatrix(rep(1, N), ncol = 1) %*% t(z) # (n*1) * (1*n)
    R <- (abs(R1-t(R1))) # n*n matrix
    kappa <- nimMatrix(rep(model$kappa, N*N), nrow = N, ncol = N)
    C <- exp(-kappa * R)
    Sigmatheta <- C

    Rtheta <- Stheta - Stheta %*% inverse(Stheta+Sigmatheta) %*% Stheta
    Rtheta <- (Rtheta + t(Rtheta))/2
    LRtheta <- t(chol(Rtheta))
    mtheta_g <- Rtheta %*% inverse(Stheta) %*% g

    # latent variable
    etaLatent <- inverse(LRtheta) %*% (oldeta-mtheta_g) # (100, 1)

    # 2) New eta, theta
    current <- oldtheta
    num.calc <- nimDim(grid.df)[1]
    PP1 <- numeric(num.calc)
    for (i in 1:num.calc){
      input <- grid.df[i, ]
      PP1[i] <- obj_fn$run(input)
    }

    idxTEMP <- (1:num.calc)[PP1 == max(PP1)]
    t <- grid.df[idxTEMP, ] # theta
    newparamebeta <- numeric(p-1) # d

    for(i in 1:(p-1)){
      newparamebeta[i] <- (paramalpha[i]-1-paramalpha[i]*t[1,i]+2*t[1,i])/t[1,i]
    }

    newtheta <- nimMatrix(rep(0, p-1), ncol = 1) #New possible theta
    for (i in 1:(p-1)){
      newtheta[i,1] <- rbeta(1,paramalpha[i],newparamebeta[i])
    }

    newbeta <- alphaTheta(newtheta[,1]*pi)
    z1 <- nimMatrix(0, nrow = N, ncol = 1)
    xData <- model$x
    for (i in 1:N) {
      for (j in 1:length(newbeta)) {
        z1[i, 1] <- z1[i, 1]  + xData[i, j] * newbeta[j]
      }
    }
    z1 <- z1[orderf,]

    R1 <- nimMatrix(rep(1, N), ncol = 1) %*% t(z1)
    R <- (abs(R1-t(R1)))
    C <- exp(-kappa * R)
    Sigmathetaprime <- C

    Rthetaprime <- Stheta-Stheta %*% inverse(Stheta+Sigmathetaprime) %*% Stheta
    Rthetaprime <- (Rthetaprime + t(Rthetaprime))/2
    LRthetaprime <- t(chol(Rthetaprime))
    mtheta_gprime <- Rtheta %*% inverse(Stheta) %*% g
    newetaMat <- LRthetaprime %*% etaLatent + mtheta_gprime
    neweta   <- newetaMat[ ,1]

    ## Beta transition kernel
    righttrans <- 0.0 # (t) -> (t+1)
    reversetrans <- 0.0 # (t+1) -> (t)
    for (i in 1:(p-1)){
      reversetrans <- reversetrans +
        dbeta(model$psi[i], paramalpha[i], newparamebeta[i], log = TRUE)
      righttrans <- righttrans +
        dbeta(newtheta[i,1], paramalpha[i], parambeta[i], log = TRUE)
    }

    # update
    model$d <<- newparamebeta
    model$psi <<- newtheta[,1]
    storeEta <- neweta[nimOrder(orderf*1.0)]
    model$linkFunction <<- storeEta
    model$calculate(calcNodes)

    # acceptance rate
    L <- (dmnorm_chol(neweta, mean=rep(0,N),
                      cholesky = chol(Sigmathetaprime), prec_param = FALSE, log=TRUE) -
            dmnorm_chol(oldeta, mean=rep(0,N),
                        cholesky = chol(Sigmatheta), prec_param = FALSE, log=TRUE) +
            dmnorm_chol(g, mean=rep(0,N),
                        cholesky = chol(Sigmathetaprime+Rthetaprime),
                        prec_param = FALSE, log=TRUE) -
            dmnorm_chol(g, mean=rep(0,N),
                        cholesky = chol(Sigmatheta+Rtheta), prec_param = FALSE, log=TRUE)) +
      reversetrans - righttrans + thetaPrior(newtheta[,1]*pi) - thetaPrior(model$psi*pi)

    ratio <- exp(L)
    u <- runif(1)

    # accept/reject
    if(u < ratio) { # accept
      nimCopy(from = model, to = mvSaved, row = 1, nodes  = calcNodes, logProb = TRUE)
    }
    if(u > ratio){
      nimCopy(from = mvSaved, to = model, row = 1, nodes  = calcNodes, logProb = TRUE)
    }
    ###########################################################################
    # Gibbs for linkFunction
    orderf <- nimOrder(model$Xlin) # integer index
    y <- nimMatrix(model$y, ncol = 1); y <- y[orderf,] # sort y
    sigma2     <- model$sigma2
    sigma2Inv  <- 1.0 / sigma2
    kappa_inv <- model$kappa[1]
    Cinv       <- invcov(model$Xlin, kappa_inv)
    A <- Cinv + diag(rep(sigma2Inv, N))
    A <- (A + t(A))/2
    Sigma <- inverse(A) # matrix
    mu    <- (Sigma %*% y) # matrix (n*n) %*% (n*1)
    mu0 <- mu[,1] # vector
    mu0 <- mu0 * rep(sigma2Inv, N)

    L1 <- chol(Sigma)
    if (any(is.nan(L1))) {
      nimPrint("Cholesky decomposition failed in linkFunction sampler.")
    }


    f <- rmnorm_chol(1, mean = mu0, cholesky = L1, prec_param = FALSE)
    reorder <- orderf * rep(1.0, N)
    idx <- nimOrder(reorder)
    f <- f[idx]

    ## update
    model$linkFunction  <<- f
    # nimPrint(f)
    model$calculate(c("linkFunction", "y"))


  },
  methods = list( reset = function () {} )

)
