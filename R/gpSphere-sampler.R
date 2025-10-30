# index - MH
indexSampler_gpSphere <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control){
    calcNodes <- model$getDependencies(target)
    p <- model$getConstants()$p
    N <- model$getConstants()$N
  },
  run = function(){
    oldIndex <- model$index
    beforell <- model$calculate(calcNodes)
    covMat <- chol(diag(nimRep(0.25, p)))
    proposal <- rmnorm_chol(n = 1, oldIndex,
                            cholesky = covMat,
                            prec_param = FALSE)
    proposal <- proposal/sqrt(sum(proposal^2))
    if (proposal[1] < 0){
      proposal <- proposal * (-1)
    }
    model$index <<- nimC(proposal)
    model$calculate(calcNodes)
    afterll <- model$calculate(calcNodes)

    L <- (afterll - beforell)
    ratio <- exp(L)
    u <- runif(1)

    # accept/reject
    if(u < ratio) { # accept
      nimCopy(from = model, to = mvSaved, row = 1, nodes  = calcNodes, logProb = TRUE)
    }
    if(u > ratio){
      nimCopy(from = mvSaved, to = model, row = 1, nodes  = calcNodes, logProb = TRUE)
    }
  },
  methods = list(reset = function(){})
)

# sigma2 - gibbs
sigma2Sampler_gpSphere <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes  <-  model$getDependencies(target)
    a_sig <- model$getConstants()$a_sig
    b_sig <- model$getConstants()$b_sig
    n <- model$getConstants()$N
    Y <- model$Y # not matrix
  },
  run = function() {
    eta <- model$linkFunction # not matrix
    # inverse gamma
    post_alpha <- a_sig  + (n/2)
    post_beta <- (1/b_sig) + 0.5 * sum((Y - eta)^2)
    proposal <- rinvgamma(n = 1, post_alpha, post_beta)
    model$sigma2 <<- nimC(proposal)
    model$calculate(calcNodes)
    nimCopy(from = model, to = mvSaved, row = 1, nodes  = target, logProb = TRUE)
  },
  methods = list( reset = function () {} )

)

# opt - MAP
optSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    calcNodes  <-  model$getDependencies(target)
    shape <- model$getConstants()$shape
    rate <- model$getConstants()$rate
    a_amp <- model$getConstants()$a_amp
    b_amp <- model$getConstants()$b_amp
    a_sig <- model$getConstants()$a_sig
    b_sig <- model$getConstants()$b_sig
    obj_fn <- obj_btt(model, length(model$index),
                      shape, rate, a_amp, b_amp,
                      a_sig, b_sig)
    lowerB <- control$lowerB
    upperB <- control$upperB
    # lowerB <- c(rep(-1, length(model$index)),-1e2, -1e2)
    # upperB <- c(rep(1, length(model$index)),1e2, 1e2)

  },
  run = function() {
    current <- c(model$index, log(model$lengthscale), log(model$amp))
    # if (sigmaTrue){
    #   current <- c(current, log(model$sigma2))
    # }
    # controlOpt <- optimDefaultControl()
    # controlOpt$factr <- 1e11
    # controlOpt$pgtol <- 1e-4
    # controlOpt$maxit <- 80


    optRes  <- nimOptim(par    = current,
                        fn     = obj_fn$run,
                        method = "L-BFGS-B",
                        lower  = lowerB,
                        upper  = upperB)

    ## update
    p1 <- length(model$index)
    proposedIndex <- optRes$par[1:p1]
    if (proposedIndex[1] < 0){
      proposedIndex <- (proposedIndex) * (-1)
    }
    model$index  <<- proposedIndex/sqrt(sum((proposedIndex)^2))
    lengthscale <- exp(optRes$par[(p1+1):(p1+1)])
    if (lengthscale[1] < 1e-10){
      lengthscale <- lengthscale + 1e-6
    }
    model$lengthscale <<- lengthscale
    model$amp   <<- exp(optRes$par[(p1+2):(p1+2)])
    # if (sigmaTrue) model$sigma2   <<- exp(optRes$par[(p1+3):(p1+3)])

    model$calculate(calcNodes)

    ## deterministic move → always accept
    nimCopy(from = model, to = mvSaved, row = 1, nodes  = target, logProb = TRUE)
  },
  methods = list( reset = function () {} )

)
