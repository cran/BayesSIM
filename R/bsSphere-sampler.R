# Sampler
## 1. nuSampler
nuSampler_bspline_sphere <- nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control){
    calcNodes <- model$getDependencies(c(target, "index", "knots"))
    p <- model$getConstants()$p
    r1 <- model$getConstants()$r1
    r2 <- model$getConstants()$r2
    r3 <- 0.5
    Rpostll <- postll_bspline_sphere(model)
  },
  run = function(){
    k <- 0; olda <- 0.0; oldb <- 0.0
    k <- model$k[1]
    currnu <- model$nu # vector[1:p]
    currTheta <- model$index # vector[1:p]
    oldknots <- model$knots # vector[1:maxknots]
    olda <- model$a_alpha[1]; oldb <- model$b_alpha[1]
    # nimPrint("################# Start iteration #################")
    beforell <- Rpostll$run() # compute
    # nimCat("beforell:", beforell, "\n")

    # proposalidx <- integer(1) # initialize
    proposalidx <- rcat(1, prob = rep(1,p)/p)
    # nimCat("Propsed idx: ", proposalidx, "\n")
    currDeltIdx <- 0; curralpha <- 0; nnu <- 0
    currDeltIdx <- currnu[proposalidx] ## current nu value
    curralpha <- currTheta[proposalidx]
    nnu <- sum(currnu)
    propnu <- currnu
    propnu[proposalidx] <- 1-currnu[proposalidx]

    if (proposalidx == 1){
      # nimPrint("Reject assumption")
    } else{
      model$nu <<- propnu

      if (currDeltIdx == 1){
        # nimPrint("change to 0")
        changedTheta <- currTheta
        changedTheta[proposalidx] <- 0
        changedTheta <- changedTheta/sqrt(sum(changedTheta^2))
        model$index <<- changedTheta
        r4 <- 0.0; r4 <- (nnu-1)/2

        # computeA1
        logA1 <- computeA1_nu1(r1, r2, nnu, p, r3, r4, curralpha)

      } else{
        # nimPrint("change to 1")
        r4 <- 0.0; eta <- 0.0; signTheta <- 0.0
        r4 <- nnu/2
        eta <- rbeta(1, r3, r4)
        changedTheta <- currTheta * sqrt(1-eta)
        signTheta <- 2*rbinom(1, 1, 1/2) -1
        changedTheta[proposalidx] <- signTheta * sqrt(eta)
        model$index <<- changedTheta

        # computeA1
        logA1 <- computeA1_nu0(r1, r2, nnu, p, r3, r4, eta)
      }

      model$calculate(c("Xlin", "a_alpha", "b_alpha")) ##
      newa <- 0.0; newb <- 0.0
      newa <- model$a_alpha[1]; newb <- model$b_alpha[1]
      changesKnots <- newKnots(oldknots, olda, oldb, newa, newb, k)
      model$knots <<- changesKnots
      model$calculate(calcNodes) ## Change all

      # nimCat("logA1:", logA1, "\n")
      afterll <- Rpostll$run()
      # nimCat("afterll:", afterll, "\n")

      L <- (afterll - beforell) + logA1
      # nimCat("logratio:", L, "\n")
      # nimPrint(L)
      ratio <- min(c(1, exp(L)))
      # nimCat("ratio:", ratio, "\n")
      u <- runif(1)

      # accept/reject
      if(u < ratio) { # accept
        # nimPrint("accept")
        # calculate knots
        tempTheta <- 0.0; tempTheta <- changedTheta[1]
        if (tempTheta < 0 ){
          # nimPrint("change sign")
          newa <- 0.0; newb <- 0.0
          newa <- model$a_alpha[1]; newb <- model$b_alpha[1]
          model$index <<- (model$index) * (-1)
          temp <- (model$knots) * (-1)
          tempSort <- nimSort(temp)
          model$nu <<- tempSort
          model$a_alpha <<- nimC((newa) * (-1))
          model$b_alpha <<- nimC((newb) * (-1))
          model$calculate(calcNodes)
        }

        nimCopy(from = model, to = mvSaved, row = 1, nodes  = calcNodes, logProb = TRUE)
      }
      if(u > ratio | changedTheta[1] == 0){
        # nimPrint("reject")
        nimCopy(from = mvSaved, to = model, row = 1, nodes  = calcNodes, logProb = TRUE)
      }
    }

  },
  methods = list(reset = function(){})
)

## 2. indexSampler
indexSampler_bspline_sphere <-  nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control){
    calcNodes <- model$getDependencies(c(target, "knots"))
    p <- model$getConstants()$p
    Rpostll <- postll_bspline_sphere(model)
    sig1_tune <- 0.1
  },
  run = function(){
    # nimPrint("################# Sampling index #################")
    k <- 0; olda <- 0.0; oldb <- 0.0
    oldknots <- model$knots
    k <- model$k[1]
    olda <- model$a_alpha[1]
    oldb <- model$b_alpha[1]
    currTheta <- model$index
    # nimCat("currTheta:", currTheta, "\n")

    beforell <- Rpostll$run() # compute
    # nimCat("beforell:", beforell, "\n")

    totalP <- which(currTheta != 0)
    numtotalP <- 0; numtotalP <- length(totalP)
    # nimCat("totalP:", totalP, "\n")
    if (numtotalP < 2){
      # nimCat("no variable selected \n")
    } else{
      idxs <- integer(2)
      idxs[1] <- rcat(1, nimRep(1 / numtotalP, numtotalP))
      idxs[2] <- rcat(1, nimRep(1 / numtotalP, numtotalP))

      while (idxs[1] == idxs[2]) {
        idxs[1] <- rcat(1, nimRep(1 / numtotalP, numtotalP))
        idxs[2] <- rcat(1, nimRep(1 / numtotalP, numtotalP))
      }

      eta <- 0.0
      idx_i <- totalP[idxs[1]]
      idx_j <- totalP[idxs[2]]
      eta <- rnormTrun(0, sig1_tune)
      # nimCat("idxs: ", idxs, "\n")
      # nimCat("idx_i: ", idx_i, "\n")
      # nimCat("idx_j: ", idx_j, "\n")
      # nimCat("eta: ", eta, "\n")

      proposedTheta <- currTheta
      proposedTheta[idx_i] <- currTheta[idx_i] * cos(eta) + currTheta[idx_j] * sin(eta)
      proposedTheta[idx_j] <- -currTheta[idx_i] * sin(eta) + currTheta[idx_j] * cos(eta)
      # nimCat("proposedTheta:", proposedTheta, "\n")

      model$index <<- proposedTheta
      model$calculate(c("Xlin", "a_alpha", "b_alpha")) ##
      newa <- 0.0; newb <- 0.0
      newa <- model$a_alpha[1]; newb <- model$b_alpha[1]
      changesKnots <- newKnots(oldknots, olda, oldb, newa, newb, k)
      model$knots <<- changesKnots

      model$calculate(calcNodes) ##
      afterll <- Rpostll$run()
      # nimCat("afterll:", afterll, "\n")

      L <- (afterll - beforell)
      # nimCat("logratio:", L, "\n")
      # nimPrint(L)
      ratio <- min(c(1, exp(L)))
      # nimPrint(ratio)
      u <- runif(1)

      # accept/reject
      if(u < ratio) { # accept
        # nimPrint("accept")
        # calculate knots
        tempTheta <- 0.0; tempTheta <- proposedTheta[1]
        if (tempTheta < 0 ){
          # nimPrint("change sign")
          newa <- 0.0; newb <- 0.0
          newa <- model$a_alpha[1]; newb <- model$b_alpha[1]
          model$index <<- (model$index) * (-1)
          temp <- (model$knots) * (-1)
          tempSort <- nimSort(temp)
          model$nu <<- tempSort
          model$a_alpha <<- nimC((newa) * (-1))
          model$b_alpha <<- nimC((newb) * (-1))
          model$calculate(calcNodes)
        }
        ## change knots
        nimCopy(from = model, to = mvSaved, row = 1, nodes  = calcNodes, logProb = TRUE)
      }
      if(u > ratio | proposedTheta[1] == 0){
        # nimPrint("reject")
        nimCopy(from = mvSaved, to = model, row = 1, nodes  = calcNodes, logProb = TRUE)
      }
    }



  },
  methods = list(reset = function(){})
)

## 3. sampling knots
knotsSampler_bspline_sphere <-  nimble::nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control){
    calcNodes <- model$getDependencies(c(target, "k"))
    p <- model$getConstants()$p[1]
    Rpostll <- postll_knots(model)
    lambda <- model$getConstants()$lambda_k[1]
    tau <- model$getConstants()$tau[1]
    r5 <- 1
    r6 <- 1
    maxknots <- model$getConstants()$maxknots[1]
    sig2_tune <- 0.5
  },
  run = function(){
    # nimPrint("################# Sampling knots #################")
    # Calculate probabilities: bk, dk, mk
    # nimCopy(from = model, to = mvSaved, row = 1, nodes = saveNodes, logProb = TRUE)
    currk <- 0; a_alpha <- 0.0; b_alpha <- 0.0
    currk <- model$k[1]
    a_alpha <- model$a_alpha[1]
    b_alpha <- model$b_alpha[1]
    bk <- prop_add(lambda, currk, tau)
    dk <- prop_delete(lambda, currk, tau)
    mk <- 1 - bk - dk
    currknots <- model$knots
    # nimCat("curr knots:", currknots, "\n")
    beforell <- Rpostll$run() # compute
    # nimCat("beforell:", beforell, "\n")
    # nimCat("probabilities:, ",  nimC(bk, dk, mk), "\n")

    choseu <- rcat(1, prob = nimC(bk, dk, mk))
    newknots <- nimNumeric(maxknots)
    # nimCat("choseu:, ", choseu, "\n")
    if (choseu == 1){
      # nimPrint("Add knots")
      # add knots
      idx <- rcat(1, prob = nimRep(1/(currk+1), (currk+1))) # interval index
      # nimCat("Location:", idx, "\n")
      minIdx <- 0.0; maxIdx <- 0.0

      if (idx == 1){
        minIdx <- a_alpha; maxIdx <- currknots[idx]
      } else if (idx == (currk + 1)){
        minIdx <- currknots[(idx-1)]; maxIdx <- b_alpha
      } else {
        minIdx <- currknots[(idx-1)]; maxIdx <- currknots[idx]
      }
      # nimCat("minIdx:", minIdx, ", maxIdx:", maxIdx, "\n")

      augEta <- 0.0
      augEta <- rbeta(1, r5, r6) * (maxIdx - minIdx) + minIdx

      for (i in 1:(maxknots-1)){
        if (i < idx){
          newknots[i] <- currknots[i]
        } else if (i == idx){
          newknots[i] <- augEta
        } else{
          newknots[i] <- currknots[(i-1)] ## index
        }
      }
      # nimCat("newknots:", newknots, "\n")
      newk <- currk + 1
      logA2 <- computeA2_add(currk, minIdx, maxIdx, a_alpha, b_alpha)
      # nimCat("logA2", logA2, "\n")
    } else if (choseu == 2){
      # nimPrint("Delete knots")
      # delete
      idx <- rcat(1, prob = nimRep(1/(currk), currk))
      # nimCat("Location:", idx, "\n")
      minIdx <- 0.0; maxIdx <- 0.0

      if (currk == 1){
        minIdx <- a_alpha; maxIdx <- b_alpha
      } else{
        if (idx == 1){
          minIdx <- a_alpha; maxIdx <- currknots[(idx+1)]
        } else if (idx == (currk)){
          minIdx <- currknots[(idx-1)]; maxIdx <- b_alpha
        } else {
          minIdx <- currknots[(idx-1)]; maxIdx <- currknots[(idx+1)]
        }
      }

      # nimCat("minIdx:", minIdx, ", maxIdx:", maxIdx, "\n")

      for (i in 1:maxknots){
        if (i < idx){
          newknots[i] <- currknots[i]
        } else if (i >= idx & i < maxknots){
          newknots[i] <- currknots[i+1]
        } else if (i == maxknots){
          newknots[i] <- 0
        }
      }
      # nimCat("newknots:", newknots, "\n")
      newk <- currk - 1
      logA2 <- computeA2_delete(currk, minIdx, maxIdx, a_alpha, b_alpha)
      # nimCat("logA2", logA2, "\n")
    } else{
      # nimPrint("Relocate knots")
      # relocate
      idx <- rcat(1, prob = nimRep(1/(currk), currk))
      # nimCat("Location:", idx, "\n")
      minIdx <- 0.0; maxIdx <- 0.0

      if (currk == 1){
        minIdx <- a_alpha; maxIdx <- b_alpha
      } else{
        if (idx == 1){
          minIdx <- a_alpha; maxIdx <- currknots[(idx+1)]
        } else if (idx == (currk)){
          minIdx <- currknots[(idx-1)]; maxIdx <- b_alpha
        } else {
          minIdx <- currknots[(idx-1)]; maxIdx <- currknots[(idx+1)]
        }
      }

      # nimCat("minIdx:", minIdx, ", maxIdx:", maxIdx, "\n")

      newLoc <- 0.0
      newLoc <- rtruncnorm(1, mean = currknots[idx],
                           sigma = (maxIdx - minIdx)*sig2_tune,
                           lower = minIdx, upper = maxIdx)[1]
      # nimCat("newLoc:", newLoc, "\n")
      newknots <- currknots
      newknots[idx] <- newLoc
      # nimCat("newknots:", newknots, "\n")
      newk <- currk
      logA2 <- 0
    }

    if (newk == 0){
      # nimPrint("no knots!")
    } else{
      model$k <<- nimC(newk)
      model$knots <<- newknots
      model$calculate(calcNodes)

      afterll <- Rpostll$run()
      # nimCat("afterll:", afterll, "\n")

      L <- (afterll - beforell) + logA2
      # nimCat("logratio:", L, "\n")
      # nimPrint(L)
      ratio <- min(c(1, exp(L)))
      # nimCat("ratio:", ratio, "\n")
      u <- runif(1)

      # accept/reject
      if(u < ratio) { # accept
        # nimPrint("accept")
        ## change knots
        nimCopy(from = model, to = mvSaved, row = 1, nodes  = calcNodes, logProb = TRUE)
      }
      if(u > ratio){
        # nimPrint("reject")
        nimCopy(from = mvSaved, to = model, row = 1, nodes  = calcNodes, logProb = TRUE)
      }
    }


  },
  methods = list(reset = function(){})
)


## 4. beta Sampler
betaSampler_bspline_sphere <- nimble::nimbleFunction(
    contains = nimble::sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      calcNodes <- model$getDependencies(target)
      N <- model$getConstants()$N
      y <- model$Y
      tau2 <- model$getConstants()$tau
      Sigma0 <- model$getConstants()$Sigma0
      maxBasis <- model$getConstants()$maxBasis
    },
    run = function(){
      numBasis <- 0; sigma2 <- 0.0; alpha <- 0.0
      numBasis <- model$numBasis[1]
      xmat <- model$Xmat
      sigma2 <- model$sigma2[1]
      alpha <- N/2
      numBasisI <- as.integer(numBasis)
      Sigma <- computeSig(tau2, Sigma0, xmat, numBasisI)
      realxmat <- xmat[,1:numBasis]

      meanBeta <- Sigma %*% t(realxmat) %*% y
      covBeta <- sigma2 * Sigma
      cholBeta <- chol(covBeta)
      newBeta <- rmnorm_chol(1, mean = meanBeta[,1], cholesky  = cholBeta)
      betaMat <- nimMatrix(0, nrow = maxBasis, ncol = 1)
      for (i in 1:numBasis){
        betaMat[i, 1] <- newBeta[i]
      }

      model$beta <<- betaMat
      model$calculate(calcNodes)
      nimCopy(from = model, to = mvSaved, row = 1, nodes  = calcNodes, logProb = TRUE)

    },
    methods = list( reset = function () {} )

  )


## 5. sigma sampler
sigma2Sampler_bspline_sphere <-  nimble::nimbleFunction(
    contains = nimble::sampler_BASE,
    setup = function(model, mvSaved, target, control) {
      calcNodes <- model$getDependencies(target)
      N <- model$getConstants()$N
      y <- model$Y
      tau2 <- model$getConstants()$tau
      Sigma0 <- model$getConstants()$Sigma0
    },
    run = function(){
      numBasis <- model$numBasis
      xmat <- model$Xmat
      alpha <- N/2
      numBasisI <- as.integer(numBasis[1])
      Sigma <- computeSig(tau2, Sigma0, xmat, numBasisI)
      SS <- ComputeS(y, xmat, Sigma, numBasisI)
      beta <- SS/2
      newSigma <- rinvgamma(1, alpha, beta)
      # nimPrint(newSigma)

      model$sigma2 <<- nimC(newSigma)
      model$calculate(calcNodes)
      nimCopy(from = model, to = mvSaved, row = 1, nodes  = calcNodes, logProb = TRUE)

    },
    methods = list( reset = function () {} )

  )

