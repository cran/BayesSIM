SamplingLambda_gp_spike<- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    calcNodes <- model$getDependencies(target)
    p <- model$getConstants()$p
    N <- model$getConstants()$N
    RllFunLambda <- llFunLambda(model)
    oldinvlambda <- model$invlambda
  },
  run = function(){
    # nimPrint("################# Start iteration #################")
    # nimPrint("Sampling Lambda")
    beforell <- RllFunLambda$run() # compute
    # nimCat("beforell:", beforell, "\n")
    logproposal <- log(oldinvlambda) + rnorm(1)
    # nimCat("Log proposal:", logproposal, "\n")
    model$invlambda <<- nimC(exp(logproposal))
    model$calculate(calcNodes)
    afterll <- RllFunLambda$run()
    # nimCat("afterll:", afterll, "\n")

    L <- (afterll - beforell)
    # nimCat("logratio:", L, "\n")
    # nimPrint(L)
    ratio <- exp(L)
    # nimCat("ratio:", ratio, "\n")
    u <- runif(1)

    # accept/reject
    if(u < ratio) { # accept
      # nimPrint("accept")
      # nimPrint(nimMatrix(nimC(model$thetastar, mvSaved[["thetastar"]][[1]]),
      #                    nrow = 18))
      # nimCat("model:", model$thetastar, "\n")
      # nimPrint("mvSaved: ")
      # nimPrint(mvSaved[["thetastar"]][[1]])
      nimCopy(from = model, to = mvSaved, row = 1, nodes  = calcNodes, logProb = TRUE)
    }
    if(u > ratio){
      # nimPrint("reject")
      # nimPrint(nimMatrix(nimC(model$thetastar, mvSaved[["thetastar"]][[1]]),
      #                    nrow = 18))
      nimCopy(from = mvSaved, to = model, row = 1, nodes  = calcNodes, logProb = TRUE)
    }
  },
  methods = list(reset = function(){})
)

## 2. Draw Joint theta, v
SamplingThetaV_gp_spike <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    calcNodes <- model$getDependencies(target)
    p <- model$getConstants()$p
    N <- model$getConstants()$N
    v <- model$nu
    theta_raw <- model$index_raw
    theta_sig <- model$getConstants()$sigma_theta
    jump <- 0.35
    RllFunThetaV1 <- llFunThetaV1(model)
    RllFunThetaV2 <- llFunThetaV2(model)
  },
  run = function(){
    # nimPrint("--------------------------------------")
    # nimPrint("Sampling v and theta")
    # nimPrint(beforell)
    # initialize
    newV <- model$nu
    newTheta <- model$index_raw
    proposalTheta <- 0
    sumv <- 0
    for (i in 1:p){
      sumv = sumv + v[i]
    }
    # nimCat("old sumv: ", sumv)

    u1 <- runif(1)
    if (u1 <= 0.5 | sumv == 0){
      # nimPrint("move 1")
      # 1) MH
      idx <- rcat(1, prob = nimRep(1/p, p))
      beforell <- RllFunThetaV1$run(idx) # compute
      beforeTrans <- transitionTheta(theta_raw[idx], theta_sig, v[idx], sumv)
      newV[idx] <- 1 - newV[idx]

      if (newV[idx] == 1){
        proposalTheta <- rnorm(1, 0, theta_sig)
        sumv <- sumv + 1
      } else{
        proposalTheta <- 0
        sumv <- sumv - 1
      }
      newTheta[idx] <- proposalTheta

      if (sum(newTheta) < 0){
        newTheta <- (-1) * newTheta
      }

      model$nu <<- newV
      model$index_raw <<- newTheta
      model$indexstar <<- newTheta
      model$calculate(calcNodes)

      afterll <- RllFunThetaV1$run(idx)
      afterTrans <- transitionTheta(proposalTheta, theta_sig, newV[idx], sumv)
      # nimCat("idx:", idx, "\n")
      # nimCat("beforell:", beforell, "\n")
      # nimCat("beforeTrans:", beforeTrans, "\n")
      # nimCat("proposalTheta:", proposalTheta, "\n")
      # nimCat("new sumv:", sumv, "\n")
      # nimCat("afterll:", afterll, "\n")
      # nimCat("afterTrans:", afterTrans, "\n")

    } else{
      # 2) RW-MH
      # nimPrint("move 2")
      beforeTrans <- 1; afterTrans <- 1
      idx <- rcat(1, prob = v/sumv)
      # idx <- 1
      # count <- 0
      # while(idx <= p & count < idxtemp) {
      #   idx <- idx + 1
      #   # print(paste("idx:", idx))
      #   if (v[idx] == 1) {
      #     # print("in")
      #     count <- count + 1
      #     # print(paste("count:",count))
      #   }
      # }

      beforell <- RllFunThetaV2$run(idx) # compute
      proposalTheta <- rnorm(1, theta_raw[idx], jump)
      newTheta[idx] <- proposalTheta
      if (sum(newTheta) < 0){
        newTheta <- (-1) * newTheta
      }

      model$nu <<- newV
      model$index_raw <<- newTheta
      model$indexstar <<- newTheta
      model$calculate(calcNodes)
      afterll <- RllFunThetaV2$run(idx)
      # nimCat("idx:", idx, "\n")
      # nimCat("beforell:", beforell, "\n")
      # nimCat("beforeTrans:", beforeTrans, "\n")
      # nimCat("proposalTheta:", proposalTheta, "\n")
      # nimCat("new sumv:", sumv, "\n")
      # nimCat("afterll:", afterll, "\n")
      # nimCat("afterTrans:", afterTrans, "\n")

    }


    ## acceptance rate - log
    L <- (afterll - beforell) - (afterTrans - beforeTrans)
    # nimPrint(L)
    ratio <- exp(L)
    # nimPrint(ratio)
    u2 <- runif(1)

    # nimCat("log acceptance rate:", L, "\n")
    # nimCat("acceptance rate:", ratio, "\n")

    # accept/reject
    if(u2 < ratio) { # accept
      # nimPrint("accept")
      # nimPrint(nimMatrix(nimC(model$thetastar, mvSaved[["thetastar"]][[1]]),
      #                    nrow = 18))
      nimCopy(from = model, to = mvSaved, row = 1, nodes  = calcNodes, logProb = TRUE)
    }
    if(u2 > ratio){
      # nimPrint("reject")
      # nimPrint(nimMatrix(nimC(model$thetastar, mvSaved[["thetastar"]][[1]]),
      #                    nrow = 18))
      nimCopy(from = mvSaved, to = model, row = 1, nodes  = calcNodes, logProb = TRUE)
    }


  },
  methods = list(reset = function(){})
)

## 3. Draw sigma2
gibbsSigma2_gp_spike <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    calcNodes <- model$getDependencies(target)
    p <- model$getConstants()$p
    N <- model$getConstants()$N
    y <- model$Y
    # z <- model$Z
    # pz <- model$getConstants()$z
    a_sig <- model$getConstants()$a_sig
    b_sig <- model$getConstants()$b_sig
  },
  run = function(){
    # nimPrint("--------------------------------------")
    # nimPrint("Sampling sigma")
    v <- model$nu
    expcovtemp <- model$Ki
    alpha <- a_sig + (N)/2
    beta <- b_sig + 0.5*(computeA1(y, expcovtemp))
    newSig <- rgamma(1, alpha, rate = beta)
    propsal <- nimC(1/newSig)
    model$sigma2 <<- propsal
    model$calculate(calcNodes)
    # nimCat("\n")
    # nimPrint(nimMatrix(nimC(model$thetastar, mvSaved[["thetastar"]][[1]]),
    #                    nrow = 18))
    nimCopy(from = model, to = mvSaved, row = 1, nodes  = calcNodes, logProb = TRUE)


    # nimCat("alpha:", alpha, "\n")
    # nimCat("beta:", 1/beta, "\n")
    # nimCat("propsal:", propsal, "\n")
    # nimCat("stored sigma:", model$sigma2, "\n")
    },
  methods = list(reset = function(){})
)

## 4. Draw gam
gibbsGam_gp_spike <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    calcNodes <- model$getDependencies(target)
    p <- model$getConstants()$p
    N <- model$getConstants()$N
    y <- model$Y
    z <- model$Z
  },
  run = function(){
    # nimPrint("--------------------------------------")
    # nimPrint("Sampling gamma")
    expcovtemp <- model$Ki
    sigma2 <- model$sigma2
    invCov <- inverse(expcovtemp)
    part1 <- inverse(t(z) %*% invCov %*% z)
    part2 <- t(z) %*% invCov %*% y
    meanGam <- part1 %*% part2
    covGam <- multiplyMatrixByConstant(part1, sigma2[1])
    # covGam <- part1 * sigma2
    cholGam <- chol(covGam)


    newGamma <- rmnorm_chol(1, mean = meanGam[,1], cholesky = cholGam,
                            prec_param = FALSE)
    # nimPrint(newGamma)

    # nimCat("meanGam:", meanGam, "\n")
    # nimCat("covGam: \n")
    # nimCat(covGam, "\n")
    # nimCat("covGam: (no sigma) \n")
    # nimCat(part1, "\n")
    # nimCat("cholGam: \n")
    # nimCat(cholGam, "\n")


    model$gam <<- newGamma
    model$calculate(calcNodes)
    # nimCat("\n")
    # nimPrint(nimMatrix(nimC(model$thetastar, mvSaved[["thetastar"]][[1]]),
    #                    nrow = 18))
    nimCopy(from = model, to = mvSaved, row = 1, nodes  = calcNodes, logProb = TRUE)

    # nimPrint("################# End iteration #################")

  },
  methods = list(reset = function(){})
)


