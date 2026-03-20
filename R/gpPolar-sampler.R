## Define samplers
## Gibbs sigma
gibbsSampler_sigma2 <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes <- model$getDependencies(target)
    N <- model$getConstants()$N
    a <- model$getConstants()$a
    b <- model$getConstants()$b
    y <- model$y
  },
  run = function() {
    # nimPrint("============Sigma=================")
    eta <- model$linkFunction
    # nimPrint(eta)
    err <- sum((y - eta)^2)/2
    # nimCat("err:", err, "\n")
    newSigma <- 1/rgamma(1, shape = a + (N/2), rate = (b + err))
    # nimCat("sigma2: ", newSigma, "\n")

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
    # nimPrint("============kappa=================")
    # nimPrint("Xlinear")
    # nimPrint(model$Xlin)
    # nimPrint("linkFunction")
    # nimPrint(model$linkFunction)

    probs <- nimNumeric(length(kappas), init = TRUE)
    idx <- nimOrder(model$Xlin)
    sortXlin <- model$Xlin[idx]
    eta <- model$linkFunction ; etaMat <- nimMatrix(eta[idx], ncol = 1) # n * 1
    # eta <- model$linkFunction; eta <- eta[idx]

    # nimPrint("sorted Xlinear")
    # nimPrint(sortXlin[1:10])
    # nimPrint("sorted linkFunction")
    # nimPrint(etaMat[1:10, 1])

    for (i in 1:length(kappas)) {
      C_inv <- invcov(sortXlin, kappas[i])
      # cat("kappa: ", kappas[i], "\n")
      # print(C_inv[1:5, 1:5])
      lp1 <- matrix(0.5 * log(det(C_inv)), nrow = 1, ncol = 1)
      lp2 <- -0.5 * (t(etaMat) %*% C_inv %*% etaMat)
      lp <- lp1 + lp2
      if (any(is.nan(lp))) {
        nimPrint("NaN detected in log posterior at kappa index ", i)
        lp <- matrix(-1e8, nrow = 1, ncol = 1)
      }
      if (lp[1,1] > 1e10){
        lp <- matrix(1e10, nrow = 1, ncol = 1)
      }
      # if (any(is.infinite(lp))){ ## 여기
      #   lp <- matrix(1e8, nrow = 1, ncol = 1)
      # }
      probs[i] <- lp[1, 1]
      # probs[i] <- lp
    }
    # nimCat("chosenIdx: ", probs, "\n")
    # nimCat("probability", probs[1:20], "\n")
    maxlp <- max(probs)
    probs <- exp(probs-maxlp)
    # probs <- exp(probs)

    chosen <- rcat(1, probs)

    storeKappa <- kappas[chosen:chosen]
    model$kappa <<- storeKappa
    # nimCat("kappa", storeKappa, "\n")

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

    obj_fn     <- obj_btt_theta(model)
    tempList <- nimSeq(0.04, 0.09, 0.01)
    num <- rcat(1, nimRep(1, 6))
    idxx <-tempList[num] + 0.01
    # idxx <- 0.02
    thetas <- seq(idxx, pi-0.01, 0.1)
    # thetas <- seq(0.1, 1-0.1, 0.05)*pi # orig
    # thetas <- seq(0.01, pi-0.01, 0.1)
    grid.df  <- do.call(expand.grid,
                         replicate(p-1, thetas, simplify = FALSE))
    grid.df <- as.matrix(grid.df)
    num.calc <- nimDim(grid.df)[1]

  },
  run = function() {
    # nimPrint("============rest=================")
    # nimPrint("Xliner")
    # nimPrint(model$Xlin)
    # nimPrint("linkFunction")
    # nimPrint(model$linkFunction)

    # 1) Old eta, theta, theta, alpha(beta)
    oldtheta <- model$psi
    parambeta <- model$d
    oldBeta <- model$index
    orderf <- nimOrder(model$Xlin) # index
    oldeta <- model$linkFunction; oldeta <- oldeta[orderf] # sort f
    Stheta <- diag(nimRep(0.1, N))

    # current: oldtheta
    PP1 <- nimNumeric(num.calc, init = TRUE)
    for (i in 1:num.calc){
      input <- grid.df[i, ]
      PP1[i] <- obj_fn$run(input, model$kappa[1], model$linkFunction)
    }

    # nimCat("theta probs:", PP1[1:10], "\n")

    idxTEMP <- nimC(1:num.calc)[PP1 == max(PP1)]
    # idxTEMP <- which.max(PP1)
    t <- grid.df[idxTEMP, ]/pi # theta
    # nimCat("theta: ", t*pi, "\n")
    # nimCat("Chosen theta:", grid.df[idxTEMP, ], "\n")

    # surrogate data
    g <- rmnorm_chol(1, mean = oldeta, cholesky = chol(Stheta),
                     prec_param= FALSE)
    # nimCat("Surrogate eta:", g[1:10], "\n")

    z <- nimMatrix(model$Xlin, ncol = 1)
    z <- z[orderf, ] # (n*1)
    R1 <- nimMatrix(rep(1, N), ncol = 1) %*% t(z) # (n*1) * (1*n)
    R <- (abs(R1-t(R1))) # n*n matrix
    # kappa <- model$kappa
    kappa <- nimMatrix(rep(model$kappa, N*N), nrow = N, ncol = N)
    C <- exp(-kappa * R)
    Sigmatheta <- C

    Rtheta <- Stheta - Stheta %*% inverse(Stheta+Sigmatheta) %*% Stheta
    Rtheta <- (Rtheta + t(Rtheta))/2
    LRtheta <- t(chol(Rtheta))
    mtheta_g <- Rtheta %*% inverse(Stheta) %*% g

    # latent variable
    etaLatent <- inverse(LRtheta) %*% (oldeta-mtheta_g) # (100, 1)

    # nimCat("Surrogate eta:", etaLatent[1:10,], "\n")

    # 2) New eta, theta
    newparamebeta <- numeric(p-1) # d
    for(i in 1:(p-1)){
      newparamebeta[i] <- (paramalpha[i]-1-paramalpha[i]*t[1, i]+2*t[1, i])/t[1, i]
    }
    # nimCat("newparamebeta: ", newparamebeta, "\n")

    newtheta <- nimMatrix(rep(0, p-1), ncol = 1) # New possible theta
    for (i in 1:(p-1)){
      newtheta[i,1] <- rbeta(1, paramalpha[i], newparamebeta[i])
    }
    # nimCat("newtheta: ", newtheta[,1]*pi, "\n")

    newbeta <- alphaTheta(newtheta[,1] * pi)
    z1 <- nimMatrix(0, nrow = N, ncol = 1)
    xData <- model$x
    z1 <- xData %*% nimMatrix(newbeta, ncol = 1)
    model$Xlin <<- z1[,1]
    # for (i in 1:N) {
    #   for (j in 1:length(newbeta)) {
    #     z1[i, 1] <- z1[i, 1]  + xData[i, j] * newbeta[j]
    #   }
    # }
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
    # nimCat("New eta:", neweta[1:10], "\n")

    ## Beta transition kernel
    righttrans <- 0.0 # (t) -> (t+1)
    reversetrans <- 0.0 # (t+1) -> (t)
    for (i in 1:(p-1)){
      reversetrans <- reversetrans +
        dbeta(model$psi[i], paramalpha[i], newparamebeta[i], log = TRUE)
      righttrans <- righttrans +
        dbeta(newtheta[i,1], paramalpha[i], parambeta[i], log = TRUE)
    }

    # nimCat("New reversetrans:", reversetrans, "\n")
    # nimCat("New righttrans:", righttrans, "\n")

    # update
    model$d <<- newparamebeta
    model$psi <<- newtheta[,1]
    storeEta <- neweta[nimOrder(orderf*1.0)]
    model$linkFunction <<- storeEta
    model$index <<- newbeta

    # model$calculate(calcNodes)
    # model$calculate(c("index", "Xlin", "cov"))
    model$calculate(c("cov"))
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
      reversetrans - righttrans + thetaPrior(newtheta[,1]*pi) - thetaPrior((model$psi)*pi)

    ratio <- exp(L)
    # nimCat("ratio: ", ratio, "\n")
    u <- runif(1)
    # accept/reject
    if(u < ratio) { # accept
      # nimPrint("accept")
      nimCopy(from = model, to = mvSaved, row = 1, nodes  = calcNodes, logProb = TRUE)
    } else {
      # nimPrint("reject")
      nimCopy(from = mvSaved, to = model, row = 1, nodes  = calcNodes, logProb = TRUE)
    }

    ###########################################################################
    # Gibbs for linkFunction
    orderf <- nimOrder(model$Xlin) # integer index
    y <- nimMatrix(model$y, ncol = 1); y <- y[orderf,] # sort y
    sigma2     <- model$sigma2
    sigma2Inv  <- 1.0 / sigma2
    kappa_inv <- model$kappa[1]
    sortXlin <- model$Xlin[orderf]
    Cinv       <- invcov(sortXlin, kappa_inv)
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

    # nimPrint("mu0")
    # nimPrint(mu0[1:10])
    # nimPrint("Sigma")
    # nimPrint(Sigma[1:5, 1:5])


    f <- rmnorm_chol(1, mean = mu0, cholesky = L1, prec_param = FALSE)
    reorder <- orderf * rep(1.0, N)
    idx <- nimOrder(reorder)
    f <- f[idx]
    # nimPrint("f - after")
    # nimPrint(f[1:10])

    ## update
    model$linkFunction  <<- f
    # nimPrint(f)
    # model$calculate(c("linkFunction"))
    nimCopy(from = model, to = mvSaved, row = 1,
            nodes  = "linkFunction", logProb = TRUE)

  },
  methods = list( reset = function () {} )

)
